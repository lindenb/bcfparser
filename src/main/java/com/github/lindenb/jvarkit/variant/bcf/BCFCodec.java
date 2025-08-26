/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


*/
package com.github.lindenb.jvarkit.variant.bcf;

import java.io.ByteArrayInputStream;
import java.io.Closeable;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.seekablestream.ISeekableStreamFactory;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.Log;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.tribble.readers.SynchronousLineReader;
import htsjdk.variant.bcf2.BCFVersion;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class BCFCodec implements Closeable {
	private static final String HTS_IDX_DELIM = "##idx##";
    private  InputStream mInputStream;    
	private static final Log LOG=Log.getInstance(BCFCodec.class);
	public static final byte[] MAGIC_HEADER_START = "BCF".getBytes();
	private AbstractVersionCodec subCodec;
	
	

	private abstract class AbstractVersionCodec {
		protected VCFHeader header=null;
		protected final List<String> idx2word =new ArrayList<>();
		protected final List<GenotypeBuilder> genotypeBuilers = new ArrayList<>();
		public abstract VCFHeader readHeader()  throws IOException;
		public abstract VariantContext decode() throws IOException;
		public abstract void rewind() throws IOException;
		public Locatable decodeLoc() throws IOException {
			return decode();
			}
		}
	
	private class BCF2_2Codec extends AbstractVersionCodec {
		byte[] buffer=null;
		long firstVariantOffset = -1L;
		
		public  VCFHeader readHeader()  throws IOException {
			@SuppressWarnings("resource")
			final BinaryCodec bc1 = new BinaryCodec(BCFCodec.this.mInputStream);

		    final int headerSizeInBytes =  uint32ToInt(bc1.readUInt());
		    final byte[] headerBytes = new byte[headerSizeInBytes];
		    bc1.readBytes(headerBytes);
		   

		    try(final PositionalBufferedStream bps = new PositionalBufferedStream(
		    		new ByteArrayInputStream(headerBytes,0,headerBytes.length-1)//skip EOF
		    		)) {
		    	  final LineIterator lr = new LineIteratorImpl(new SynchronousLineReader(bps));
		    	super.header=(VCFHeader)new VCFCodec().readActualHeader(lr);
		    	}
		    for(int x=0;x< super.header.getNGenotypeSamples();++x)  {
		    	super.genotypeBuilers.add(new GenotypeBuilder(super.header.getGenotypeSamples().get(x)));
		    }
		    
		    final List<VCFContigHeaderLine> L=this.header.getContigLines();
		    if(L.isEmpty()) throw new IOException("no contig line");
		    
			final Set<String> seen=new HashSet<>();
		    idx2word.add(VCFConstants.PASSES_FILTERS_v4);
		    seen.add(VCFConstants.PASSES_FILTERS_v4);

		    for(VCFHeaderLine hl:this.header.getMetaDataInInputOrder()) {
		    	String s=null;
		    	if(hl instanceof VCFFilterHeaderLine) {
		    		s=VCFFilterHeaderLine.class.cast(hl).getID();
		    		}
		    	else if(hl instanceof VCFInfoHeaderLine) {
		    		s=VCFInfoHeaderLine.class.cast(hl).getID();
		    		}
		    	else if(hl instanceof VCFFormatHeaderLine) {
		    		s=VCFFormatHeaderLine.class.cast(hl).getID();
		    		}
		    	if(s==null || seen.contains(s)) continue;
		    	this.idx2word.add(s);
		    	seen.add(s);
		    	}
		   if(isSupportingRandomAccess()) {
			   	firstVariantOffset= getPosition();
		   		}
		    return header;
			}
		
		private byte[] fillBuffer(BinaryCodec bc,int n)
			{
			if(this.buffer==null || this.buffer.length<n) {
				this.buffer=new byte[n];
				}
			bc.readBytes(buffer, 0, n);
			return this.buffer;
			}
		
		@Override
		public void rewind() throws IOException {
			seek(this.firstVariantOffset);
			}
		
		private Interval decodeLoc(BinaryCodec bc) {
			final int tid = bc.readInt();
			final SAMSequenceRecord ssr = this.header.getSequenceDictionary().getSequence(tid);
			final int pos0= bc.readInt();
			int rlen= bc.readInt();
			return new Interval(ssr.getContig(),pos0+1,pos0+rlen);
			}
		
		private int readInfoLength(BinaryCodec bc) throws EOFException {
			final long to_info_length ;
			try {
				to_info_length= bc.readUInt();
				}
			catch(Throwable err) {
				throw new EOFException();
				}
			return uint32ToInt(to_info_length);
			}
		
		public Locatable decodeLoc( ) throws IOException {
			@SuppressWarnings("resource")
			final BinaryCodec binaryCodec1 = new BinaryCodec(mInputStream);
			final int to_info_length ;
			try {
				to_info_length=readInfoLength(binaryCodec1);
				}
			catch(EOFException err) {
				return null;
				}
			int format_length = uint32ToInt(binaryCodec1.readUInt());
			Interval interval = null;
			
			fillBuffer(binaryCodec1,to_info_length);
			try(ByteArrayInputStream is=new ByteArrayInputStream(this.buffer,0,to_info_length)) {
				final BinaryCodec bc2 = new BinaryCodec(is);
				interval = decodeLoc(bc2);
				}
			mInputStream.skipNBytes(format_length);
			return interval;
			}
		
		public VariantContext decode() throws IOException {
			@SuppressWarnings("resource")
			final BinaryCodec binaryCodec1 = new BinaryCodec(mInputStream);
			final int to_info_length ;
			
			final List<Allele> alleles;
			try {
				to_info_length=readInfoLength(binaryCodec1);
				}
			catch(EOFException err) {
				return null;
				}		
			LOG.debug(to_info_length);
			int format_length = uint32ToInt(binaryCodec1.readUInt());
			LOG.debug(format_length);

			
			final VariantContextBuilder vcb=new VariantContextBuilder();
		
			int n_fmt=0;
			
			fillBuffer(binaryCodec1,to_info_length);
			try(ByteArrayInputStream is=new ByteArrayInputStream(this.buffer,0,to_info_length)) {
				final BinaryCodec bc2 = new BinaryCodec(is);
				final Interval interval  = decodeLoc(bc2);
				vcb.chr(interval.getContig());
				vcb.start(interval.getStart());
				vcb.stop(interval.getEnd());
				
				LOG.debug("interval= "+interval);
				
				/** QUAL ********************************************/
				float qual= bc2.readFloat();
				if(!Float.isNaN(qual) && qual!=BCFTypedData.bcf_float_missing) {
					LOG.debug("qual ok="+qual);
					vcb.log10PError(qual/-10);
					}
				LOG.debug("qual="+qual);
				
				int n_info = bc2.readUShort();
				LOG.debug("n_info="+n_info);
				int n_allele = bc2.readUShort();
				LOG.debug("n_allele="+n_allele);
				
				// n samples is in 3 bytes but should be the same as header.n_samples
				bc2.readByte();bc2.readByte();bc2.readByte();
				int n_samples= this.header.getNGenotypeSamples();
				LOG.debug("n_samples"+n_samples);
				for(GenotypeBuilder gb:this.genotypeBuilers) {
					gb.reset(true);//keep sample name
					//gb.unfiltered();
					}
				
				n_fmt = bc2.readUByte();
				LOG.debug("n_fmt="+n_fmt);

				/** ID **/
				LOG.debug("read ID");
				final String id = BCFTypedData.readString(bc2);
				if(id!=null && !id.isEmpty()) {
					LOG.debug("ID "+id);
					vcb.id(id);
					}
				
				/** ALLELES **/
				LOG.debug("n-Alleles "+n_allele);
				alleles=new ArrayList<>(n_allele);
				for(int i=0;i< n_allele;i++) {
					final Allele a= Allele.create(BCFTypedData.readString(bc2), i==0);
					LOG.debug("allele["+i+"]="+a);
					alleles.add(a);
					}
				vcb.alleles(alleles);
				
				/** FILTERS **/
				LOG.debug("read filters...");
				int[] filter_idx=BCFTypedData.readIntArray(bc2);
				if(filter_idx!=null && filter_idx.length==1 && filter_idx[0]==0/* PASS is always Oth item */ ) {
					vcb.passFilters();
					}
				else if(filter_idx!=null && filter_idx.length>0) {
					vcb.filters(Arrays.stream(filter_idx).mapToObj(it->idx2word.get(it)).collect(Collectors.toSet()));
					}
				
				LOG.debug(Arrays.toString(filter_idx));
				LOG.debug();
				
				LOG.debug("read info...");
				/** INFO */
				for(int i=0;i< n_info;++i) {
					int tag_id = BCFTypedData.read(bc2).intValue();
					LOG.debug("tag_id="+tag_id);
					LOG.debug("tag="+idx2word.get(tag_id));
					BCFTypedData td=BCFTypedData.read(bc2);
					final String tag=idx2word.get(tag_id);
					final VCFInfoHeaderLine hinfo= this.header.getInfoHeaderLine(tag);
					if(hinfo==null) throw new IllegalStateException();
					LOG.debug("INFO/"+tag + "= "+td);
					//Object value= BCFTypedData.convertToVCFtype(hinfo.getType(), td.getValue());
					Object value= td.getValue();
					
					
					
					if(value instanceof String && hinfo.getType()==VCFHeaderLineType.String) {
						String[] ss= String.class.cast(value).split("[,]");
						if(ss.length==1) {
							value=ss[0];
							}
						else
							{
							value=Arrays.asList(ss);
							}
						}
					
					else if(value instanceof List && List.class.cast(value).size()==1) {
						value= List.class.cast(value).get(0);
						}
					
					else if(value==null) value=Boolean.TRUE;//FLAG
					
					
					
					vcb.attribute(tag, value);
					
					
					/*
					if(tag.equals("SVLEN")) System.err.println("SVLEN="+td+" "+interval);
					if(tag.equals("SVLEN") && interval.getStart()==interval.getEnd() && td.getCount()==1) {
						int svlen = Math.abs(Integer.parseInt(value.toString()));
						System.err.println("SVLEN="+td+" stop="+(interval.getEnd()+svlen-1));
						vcb.stop(interval.getEnd()+svlen-1);
						}*/
					

					}
				}
			LOG.debug("Attributes = #####################################################");
			fillBuffer(binaryCodec1,format_length);
			
			try(ByteArrayInputStream is=new ByteArrayInputStream(this.buffer,0,format_length)) {
				final BinaryCodec bc2 = new BinaryCodec(is);
				
				for(int i=0;i< n_fmt;i++) {
					int tag_id = BCFTypedData.read(bc2).intValue();
					final String tag = idx2word.get(tag_id);
					final VCFFormatHeaderLine hFmtLine= super.header.getFormatHeaderLine(tag);
					if(hFmtLine==null) throw new IllegalArgumentException("FORMAT/"+tag+" used but not defined int header");
					final byte b= bc2.readByte();
					LOG.debug("FORMAT/"+tag +" has fmt_type="+BCFTypedData.decodeType(b));
					int n_element= BCFTypedData.decodeCount(bc2,b);
					BCFTypedData.Type type = BCFTypedData.decodeType(b);
					for(int x=0;x< this.header.getNGenotypeSamples();++x) {
						final GenotypeBuilder gb=this.genotypeBuilers.get(x);
						boolean phased=false;
						final List<Allele> gt_alleles=new ArrayList<>(n_element);
						final List<Object> gt_values= new ArrayList<>(n_element);
						
						if(BCFTypedData.decodeType(b).equals(BCFTypedData.Type.CHAR)) {
							final String s=BCFTypedData.readString(bc2,n_element);
							LOG.debug("add string for "+tag +" length="+n_element+" value="+s);
							gt_values.add(s);
							}
						else {
							for(int j=0;j< n_element;++j) {
								final Object o = BCFTypedData.readAtomic(bc2,type);
								if(o.equals(type.getEndVector())) {
									LOG.debug("end vector, break v="+o+" "+type.getEndVector()+"/"+type.getMissing());
									break;
									}
								if(o.equals(type.getMissing()))
									{
									LOG.debug("[sample"+x+"]["+j+"]"+tag+"=MISING");
									if(tag.equals(VCFConstants.GENOTYPE_FILTER_KEY)) {
										gt_values.add(VCFConstants.UNFILTERED);
										}
									else
										{
										gt_values.add(null);
										}
									}
								if(o instanceof Integer) {
									final Integer v0 = Integer.class.cast(o);
									if(tag.equals(VCFConstants.GENOTYPE_KEY)) {
										final int v=v0.intValue();
										int allele_idx=((v>>1)-1);
										phased = (v & 0x01) == 1;
										LOG.debug("[sample"+x+"]["+j+"]="+allele_idx+" "+alleles);
										Allele a= allele_idx<0?Allele.NO_CALL:alleles.get(allele_idx);
										gt_alleles.add(a);
										}
									else
										{
										gt_values.add(v0);
										}
									//final boolean phased = ((encoded.length > 1 ? encoded[1] : encoded[0]) & 0x01) == 1;
									//LOG.debug("[sample"+x+"]["+j+"]="+(v>>1));
									}
								else
									{
									gt_values.add(castType(hFmtLine.getType(),o));
									}
								}
							}
						if(tag.equals(VCFConstants.GENOTYPE_KEY) && !gt_alleles.isEmpty()) {
							gb.alleles(gt_alleles);
							gb.phased(phased);
							}
						else if(tag.equals(VCFConstants.GENOTYPE_QUALITY_KEY)) {
							if(gt_values.isEmpty()) {
								//nothing
								}
							else if(gt_values.size()==1) {
								int gq = Integer.class.cast(gt_values.get(0));
								gb.GQ(gq);
								}
							else
								{
								throw new IllegalArgumentException();
								}
							}
						else if(tag.equals(VCFConstants.DEPTH_KEY)) {
							if(gt_values.isEmpty()) {
								//nothing
								}
							else if(gt_values.size()==1) {
								int dp = Integer.class.cast(gt_values.get(0));
								gb.DP(dp);
								}
							else
								{
								throw new IllegalArgumentException();
								}
							}
						else if(tag.equals(VCFConstants.GENOTYPE_ALLELE_DEPTHS)) {
							int[] ads = gt_values.stream().mapToInt(it->Integer.class.cast(it).intValue()).toArray();
							gb.AD(ads);
							}
						else if(tag.equals(VCFConstants.GENOTYPE_PL_KEY)) {
							int[] pls = gt_values.stream().mapToInt(it->Integer.class.cast(it).intValue()).toArray();
							gb.PL(pls);
							}
						else if(tag.equals(VCFConstants.GENOTYPE_FILTER_KEY)) {
							LOG.debug("X3 add string for "+tag +" as "+gt_values);
							if(gt_values.isEmpty()) {
								LOG.debug("X4 add string for "+tag +" as .");

								gb.unfiltered();
								}
							else if(gt_values.size()==1) {
								String filt =(String)gt_values.get(0);
								if(filt.equals(VCFConstants.UNFILTERED)) filt=null;
								LOG.debug(">>>43 gb.filter  for "+tag +" as "+filt);

								gb.unfiltered();
								if(filt!=null) gb.filter(filt);
								}
							else
								{
								throw new IllegalArgumentException("cannot set filter for "+gt_values);
								}
							}
						else
							{
							LOG.debug(">>>>DEFAULT "+tag+"="+gt_values);
							}
						LOG.debug(">>>>X1 "+tag+"="+gt_values);
						}
					}
				if(super.header.hasGenotypingData()) {
					vcb.genotypes(this.genotypeBuilers.stream().map(GB->GB.make()).collect(Collectors.toList()));
					}
				}
			LOG.debug("EOF #####################################################");
			
			return vcb.make();
			}

		
		
		
		}
	
	
private static Object castType(final VCFHeaderLineType type,Object o) {
	if(o==null) return null;
	if(o instanceof List) {
		List<?> L=new ArrayList<>();
		return L.stream().map(it->castType(type,it)).collect(Collectors.toList());
		}
	else if(o instanceof String) {
		if(type.equals(VCFHeaderLineType.Float)) return Float.valueOf(String.class.cast(o));
		if(type.equals(VCFHeaderLineType.Integer)) return Integer.valueOf(String.class.cast(o));
		return o;
	} else {
		return o;
	}
}
	

private BCFCodec(InputStream bci) {
	this.mInputStream = bci;
}
	
	
static BCFVersion readVersion(InputStream in)  throws IOException {
	@SuppressWarnings("resource")
	final BinaryCodec binaryCodec=new BinaryCodec(in);
	final byte[] magicBytes = new byte[MAGIC_HEADER_START.length];
    binaryCodec.readBytes(magicBytes);
    if (!Arrays.equals(magicBytes, MAGIC_HEADER_START) ) throw new IOException("Cannot read BCF MAGIC");
    final int majorByte =  binaryCodec.readUByte();
    final int minorByte =  binaryCodec.readUByte();
    return new BCFVersion(majorByte, minorByte);
	}

public VCFHeader readHeader()  throws IOException {
	@SuppressWarnings("resource")
	final BCFVersion version =readVersion(this.mInputStream);
	
    if(version.getMajorVersion()==2 && version.getMinorVersion()==2) {
    	this.subCodec=new BCF2_2Codec();
    	}
    else
    	{
    	throw new IOException("Bad BCF version. Not handled "+version);
    	}
    return this.subCodec.readHeader();
	}

private boolean isSupportingRandomAccess() {
	return this.mInputStream instanceof BlockCompressedInputStream;
}

private long getPosition() {
	if(this.mInputStream instanceof BlockCompressedInputStream) {
		return BlockCompressedInputStream.class.cast(this.mInputStream).getFilePointer();
	} else
	{
		throw new IllegalArgumentException("not a seekable stream");
	}
}

public void seek(long position) throws IOException {
	if(this.mInputStream instanceof BlockCompressedInputStream) {
		BlockCompressedInputStream.class.cast(this.mInputStream).seek(position);
	} else
	{
		throw new IllegalArgumentException("not a seekable stream");
	}
	}

public VariantContext decode() throws IOException {
	 return this.subCodec.decode();
	}


public Locatable decodeLoc() throws IOException {
	return this.subCodec.decodeLoc();
	}

private static int uint32ToInt(long v) {
	if(v<0L || v > Integer.MAX_VALUE) throw new IllegalArgumentException();
	return (int)v;
	}

public void rewind() throws IOException {
	this.subCodec.rewind();
}

@Override
public void close() {
	try {
		this.mInputStream.close();
		}
	catch(IOException err) {
		LOG.warn(err, "close");	
		}
	}

public static BCFCodec open(String path) throws IOException {
	path = splitPathAndIndex(path)[0];
	ISeekableStreamFactory ssf = SeekableStreamFactory.getInstance();
	SeekableStream ss = ssf.getStreamFor(path);
	ss= ssf.getBufferedStream(ss);
	BlockCompressedInputStream bci=new BlockCompressedInputStream(ss);
		
	return new BCFCodec(bci);
	}

public static BCFCodec open(InputStream in) throws IOException {
	return new BCFCodec(in);
	}

public static String[] splitPathAndIndex(String path) {
	int i = path.indexOf(HTS_IDX_DELIM);
	if(i==-1) return new String[] {path,null};
	return new String[] {
		path.substring(0, i),
		path.substring(i+HTS_IDX_DELIM.length())
		};
	}
}
