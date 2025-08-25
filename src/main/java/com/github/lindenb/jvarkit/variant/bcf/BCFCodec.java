package com.github.lindenb.jvarkit.variant.bcf;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.BinaryFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.SimpleFeature;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.tribble.readers.SynchronousLineReader;
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
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class BCFCodec extends BinaryFeatureCodec<VariantContext> {
	public static final byte[] MAGIC_HEADER_START = "BCF".getBytes();
	private AbstractVersionCodec subCodec;
	
	

	private abstract class AbstractVersionCodec {
		protected VCFHeader header=null;
		protected final List<String> idx2word =new ArrayList<>();
		protected final List<GenotypeBuilder> genotypeBuilers = new ArrayList<>();
		public abstract FeatureCodecHeader readHeader(PositionalBufferedStream in)  throws IOException;
		public abstract VariantContext decode(PositionalBufferedStream source) throws IOException;
		public Locatable decodeLoc(PositionalBufferedStream in) throws IOException {
			return decode(in);
			}
		}
	
	private class BCF2_2Codec extends AbstractVersionCodec {
		byte[] buffer=null;
		public  FeatureCodecHeader readHeader(PositionalBufferedStream in)  throws IOException {
			@SuppressWarnings("resource")
			final BinaryCodec bc1 = new BinaryCodec(in);

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
		    System.err.println("codec.POSITION="+in.getPosition());
		    return new FeatureCodecHeader(header, in.getPosition());
			}
		
		private byte[] fillBuffer(BinaryCodec bc,int n)
			{
			if(this.buffer==null || this.buffer.length<n) {
				this.buffer=new byte[n];
				}
			bc.readBytes(buffer, 0, n);
			return this.buffer;
			}
		
		private Interval decodeLoc(BinaryCodec bc) {
			final int tid = bc.readInt();
			final SAMSequenceRecord ssr = this.header.getSequenceDictionary().getSequence(tid);
			final int pos0= bc.readInt();
			int rlen= bc.readInt();
			return new Interval(ssr.getContig(),pos0+1,pos0+rlen);
			}
		
		public Locatable decodeLoc(PositionalBufferedStream in) throws IOException {
			@SuppressWarnings("resource")
			final BinaryCodec binaryCodec1 = new BinaryCodec(in);
			int to_info_length = uint32ToInt(binaryCodec1.readUInt());
			int format_length = uint32ToInt(binaryCodec1.readUInt());
			Interval interval = null;
			
			fillBuffer(binaryCodec1,to_info_length);
			try(ByteArrayInputStream is=new ByteArrayInputStream(this.buffer,0,to_info_length)) {
				final BinaryCodec bc2 = new BinaryCodec(is);
				interval = decodeLoc(bc2);
				}
			in.skipNBytes(format_length);
			return interval;
			}
		
		public VariantContext decode(PositionalBufferedStream in) throws IOException {
			@SuppressWarnings("resource")
			final BinaryCodec binaryCodec1 = new BinaryCodec(in);
			int to_info_length = uint32ToInt(binaryCodec1.readUInt());
			System.err.println(to_info_length);
			int format_length = uint32ToInt(binaryCodec1.readUInt());
			System.err.println(format_length);

			
			final VariantContextBuilder vcb=new VariantContextBuilder();
		
			int n_fmt=0;
			
			fillBuffer(binaryCodec1,to_info_length);
			try(ByteArrayInputStream is=new ByteArrayInputStream(this.buffer,0,to_info_length)) {
				final BinaryCodec bc2 = new BinaryCodec(is);
				final Interval interval  = decodeLoc(bc2);
				vcb.chr(interval.getContig());
				vcb.start(interval.getStart());
				vcb.stop(interval.getEnd());
				
				float qual= bc2.readFloat();
				if(qual==0x7F800001) qual=-1f;
				System.err.println("qual="+qual);
				
				int n_info = bc2.readUShort();
				System.err.println("n_info="+n_info);
				int n_allele = bc2.readUShort();
				System.err.println("n_allele="+n_allele);
				
				// n samples is in 3 bytes but should be the same as header.n_samples
				bc2.readByte();bc2.readByte();bc2.readByte();
				int n_samples= this.header.getNGenotypeSamples();
				System.err.println("n_samples"+n_samples);
				for(GenotypeBuilder gb:this.genotypeBuilers) {
					gb.reset(true);//keep sample name
					}
				
				n_fmt = bc2.readUByte();
				System.err.println("n_fmt="+n_fmt);

				/** ID **/
				final String id = BCFTypedData.readString(bc2);
				if(id!=null && !id.isEmpty()) vcb.id(id);
				
				/** ALLELES **/
				final List<Allele> alleles=new ArrayList<>(n_allele);
				for(int i=0;i< n_allele;i++) {
					final Allele a= Allele.create(BCFTypedData.readString(bc2), i==0);
					alleles.add(a);
					}
				vcb.alleles(alleles);
				
				/** FILTERS **/
				System.err.println("read filters...");
				int[] filter_idx=BCFTypedData.readIntArray(bc2);
				if(filter_idx!=null && filter_idx.length==1 && filter_idx[0]==0/* PASS is always Oth item */ ) {
					vcb.passFilters();
					}
				else if(filter_idx!=null && filter_idx.length>0) {
					vcb.filters(Arrays.stream(filter_idx).mapToObj(it->idx2word.get(it)).collect(Collectors.toSet()));
					}
				
				System.err.println(Arrays.toString(filter_idx));
				System.err.println();
				
				System.err.println("read info...");
				/** INFO */
				for(int i=0;i< n_info;++i) {
					int tag_id = BCFTypedData.read(bc2).intValue();
					System.err.println("tag_id="+tag_id);
					System.err.println("tag="+idx2word.get(tag_id));
					BCFTypedData td=BCFTypedData.read(bc2);
					final String tag=idx2word.get(tag_id);
					System.err.println("INFO/"+tag + "= "+td);
					Object value= td.getValue();
					if(value==null) value=Boolean.TRUE;//FLAG
					vcb.attribute(tag, value);
					}
				}
			System.err.println("GT = #####################################################");
			fillBuffer(binaryCodec1,format_length);
			
			try(ByteArrayInputStream is=new ByteArrayInputStream(this.buffer,0,format_length)) {
				final BinaryCodec bc2 = new BinaryCodec(is);
				
				for(int i=0;i< n_fmt;i++) {
					int tag_id = BCFTypedData.read(bc2).intValue();
					System.err.println("FMT-id="+tag_id);
					final String tag = idx2word.get(tag_id);
					byte b= bc2.readByte();
					System.err.println("FORMAT/"+tag +" has fmt_type="+BCFTypedData.decodeType(b));
					int n_element= BCFTypedData.decodeCount(bc2,b);
					BCFTypedData.Type type = BCFTypedData.decodeType(b);
					for(int x=0;x< this.header.getNGenotypeSamples();++x) {
						final GenotypeBuilder gb=this.genotypeBuilers.get(x);
						
						
						for(int j=0;j< n_element;++j) {
							Object o=BCFTypedData.readAtomic(bc2,type);
							if(o.equals(type.getEndVector())) {
								System.err.println("end vector, break v="+o+" "+type.getEndVector()+"/"+type.getMissing());
								break;
								}
							if(o.equals(type.getMissing()))
								{
								System.err.println("[sample"+x+"]["+j+"]=MISING");
								}
							if(o instanceof Integer) {
								int v = Integer.class.cast(o);
								//final boolean phased = ((encoded.length > 1 ? encoded[1] : encoded[0]) & 0x01) == 1;
								System.err.println("[sample"+x+"]["+j+"]="+(v>>1));
								}
							else
								{
								System.err.println("[sample"+x+"]["+j+"]="+(o));
								}
							}
						}
					}
				if(super.header.hasGenotypingData()) {
					vcb.genotypes(this.genotypeBuilers.stream().map(GB->GB.make()).collect(Collectors.toList()));
					}
				}
			System.err.println("EOF #####################################################");
			
			return vcb.make();
			}
		
		}
@Override
public Class<VariantContext> getFeatureType() {
    return VariantContext.class;
}

@Override
public FeatureCodecHeader readHeader(PositionalBufferedStream in)  throws IOException {
	@SuppressWarnings("resource")
	final BinaryCodec binaryCodec=new BinaryCodec(in);
	
	final byte[] magicBytes = new byte[MAGIC_HEADER_START.length];
    binaryCodec.readBytes(magicBytes);
    if (!Arrays.equals(magicBytes, MAGIC_HEADER_START) ) throw new IOException("Cannot read BCF MAGIC");
    final int majorByte =  binaryCodec.readUByte();
    final int minorByte =  binaryCodec.readUByte();
    if(majorByte==2 && minorByte==2) {
    	this.subCodec=new BCF2_2Codec();
    	}
    else
    	{
    	throw new IOException("Bad BCF version. Not handled "+majorByte+"."+minorByte);
    	}
    return this.subCodec.readHeader(in);
	}

@Override
public VariantContext decode(PositionalBufferedStream in) throws IOException {
	 return this.subCodec.decode(in);
	}

@Override
public Feature decodeLoc(PositionalBufferedStream in) throws IOException {
	final Locatable r= this.subCodec.decodeLoc(in);
	return new SimpleFeature(r.getContig(),r.getStart(),r.getEnd());
	}

private static int uint32ToInt(long v) {
	if(v<0L || v > Integer.MAX_VALUE) throw new IllegalArgumentException();
	return (int)v;
	}
@Override
public boolean canDecode(String path) {
	return path.endsWith(".bcf");
	}
}
