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

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;


import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**
 * Implementation of VariantContextWriter for BCF. 
 * Currently in dev, DO NOT USE
 */
class BCFWriter implements VariantContextWriter {
	private final BlockCompressedOutputStream bgzout;
	private VCFHeader header;
	private final Map<String,Integer> word2idx=new HashMap<>();

	BCFWriter(OutputStream os,final Path p) throws IOException {
		this.bgzout = (os instanceof BlockCompressedOutputStream? BlockCompressedOutputStream.class.cast(os): new BlockCompressedOutputStream(os,p));
		}

	
	BCFWriter(Path p) throws IOException {
		this(new BlockCompressedOutputStream(p, BlockCompressedOutputStream.getDefaultCompressionLevel(), BlockCompressedOutputStream.getDefaultDeflaterFactory()),p);
		}
	
	private void initHeader(final VCFHeader header) {
		if(this.header!=null) throw new IllegalStateException("header specified twice");
		this.header = header;
		this.word2idx.put(VCFConstants.PASSES_FILTERS_v4,0);
		
		
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
	    	if(s==null || this.word2idx.containsKey(s)) continue;
	    	this.word2idx.put(s, word2idx.size());
	    	}
		
		@SuppressWarnings("resource")
		final BinaryCodec bc=new BinaryCodec(this.bgzout);
		bc.writeBytes(BCFCodec.MAGIC_HEADER_START);
		bc.writeByte(2);
		bc.writeByte(2);
		try(ByteArrayOutputStream baos=new ByteArrayOutputStream()) {
			new VariantContextWriterBuilder().
				setOutputVCFStream(baos).
				build().
				writeHeader(header);
			baos.flush();
			byte[] array = baos.toByteArray();
			bc.writeInt(array.length);//headerSizeInBytes
			bc.writeBytes(array);
			}
		catch(IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	
	@Override
	public void writeHeader(final VCFHeader header) {
		initHeader(header);
		}

	@Override
	public void close() {
		try {
			this.bgzout.close();
		} catch (IOException e) {
			throw new RuntimeIOException(e);
		}

	}

	@Override
	public boolean checkError() {
		return false;
	}

	private byte[] toInfo(final VariantContext vc) {
		final SAMSequenceDictionary dict = this.header.getSequenceDictionary();
		try(ByteArrayOutputStream baos=new ByteArrayOutputStream()) {
			@SuppressWarnings("resource")
			final BinaryCodec bc1=new BinaryCodec(baos);
			
			final int tid = dict.getSequenceIndex(vc.getContig());
			if(tid<0) throw new IllegalStateException("undefined contig "+vc.getContig()+" in vcf dictionary");
			bc1.writeInt(tid);
			
			final int pos0= vc.getStart()-1;
			bc1.writeInt(pos0);
			final int rlen = vc.getLengthOnReference();
			bc1.writeInt(rlen);
			
			if(vc.hasLog10PError()) {
				bc1.writeFloat((float)vc.getPhredScaledQual());
				}
			else
				{
				bc1.writeFloat(BCFTypedData.bcf_float_missing);
				}
			final Map<String,Object> atts = vc.getAttributes();
			final int n_info = atts.size();
			bc1.writeUShort(n_info);
			final int n_alleles = vc.getNAlleles();
			bc1.writeUShort(n_alleles);
			final int n_samples = this.header.getNGenotypeSamples();
			//TODO write sample as 3 bytes
			
			//TODO write ID
			
			for(Allele a: vc.getAlleles()) {
				writeString(bc1,a.getDisplayString());
			}
			
			//write FILTERS
			
			for(final String key :atts.keySet()) {
				
				}
			
			
			
			return baos.toByteArray();
			}
		 catch (IOException e) {
			throw new RuntimeIOException(e);
			}
		}
	
	private byte[] toFormat(final VariantContext vc) {
		try(ByteArrayOutputStream baos=new ByteArrayOutputStream()) {
			@SuppressWarnings("resource")
			final BinaryCodec bc1=new BinaryCodec(baos);

			
			
			return baos.toByteArray();
			}
		 catch (IOException e) {
			throw new RuntimeIOException(e);
			}
		}
		
	@Override
	public void add(VariantContext vc) {
		if(this.header==null) throw new IllegalStateException("header was not specified");
		@SuppressWarnings("resource")
		final BinaryCodec bc1=new BinaryCodec(this.bgzout);
		final byte[] array1 = toInfo(vc);
		final byte[] array2 = toFormat(vc);
		bc1.writeUInt(array1.length);
		bc1.writeBytes(array1);
		bc1.writeUInt(array2.length);
		bc1.writeBytes(array2);
		}
	
	@Override
	public void setHeader(VCFHeader header) {
		initHeader(header);
	}

	
	private void writeString(BinaryCodec bc,final String s) {
		//TODO
		}
}
