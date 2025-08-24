package com.github.lindenb.jvarkit.variant.bcf;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class BCFCodecTest {
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][] {
			{"testdata/test01.bcf","testdata/test01.vcf"}
			};
		}
	static void compareLoc(Locatable f1,Locatable f2) {
		Assert.assertEquals(f1.getContig(), f2.getContig());
		Assert.assertEquals(f1.getStart(), f2.getStart());
		Assert.assertEquals(f1.getEnd(), f2.getEnd());
	}

	
	static void compareVariantContext(VariantContext f1,VariantContext f2) {
		compareLoc(f1,f2);
		Assert.assertEquals(f1.hasID(), f2.hasID());
		Assert.assertEquals(f1.getAlleles(), f2.getAlleles());
		Assert.assertEquals(f1.filtersWereApplied(), f2.filtersWereApplied());
		Assert.assertEquals(f1.getFilters(), f2.getFilters());
		if(f1.hasID()) {
			Assert.assertEquals(f1.getID(), f2.getID());
			}
		Assert.assertEquals(f1.getAttributes().size(), f2.getAttributes().size());
		Assert.assertEquals(f1.getAttributes().keySet(), f2.getAttributes().keySet());
		for(final String k:f1.getAttributes().keySet()) {
			Object o1 = f1.getAttribute(k);
			Assert.assertNotNull(o1,"key="+k);
			Object o2 = f2.getAttribute(k);
			Assert.assertNotNull(o2,"key="+k);
			Assert.assertEquals(o1.toString(), o2.toString(),"not the same["+k+"] ("+o1+":"+o1.getClass()+") vs ("+o2+":"+o2.getClass()+")");
			}
	}
	
	@Test(dataProvider="src1")
	public void read(String bcffname,String vcfname) throws IOException{
		List<VariantContext> L1=new ArrayList<>();

		BCFCodec codec=new BCFCodec();
		try(BlockCompressedInputStream in=new BlockCompressedInputStream(Paths.get(bcffname))) {
			try(PositionalBufferedStream pb=new PositionalBufferedStream(in)) {
				codec.readHeader(pb);
				while(!pb.isDone()) {
					L1.add(codec.decode(pb));
					}
				}
			}
		catch(Throwable err) {
			err.printStackTrace();
			throw new IOException(err);
		}
		
		List<VariantContext> L2=new ArrayList<>();
		VCFCodec codec2=new VCFCodec();
		try(PositionalBufferedStream pb=new PositionalBufferedStream(Files.newInputStream(Paths.get(vcfname)))) {
			LineIterator li= new LineIteratorImpl(AsciiLineReader.from(pb)); 
				codec2.readActualHeader(li);
				while(li.hasNext()) {
					L2.add(codec2.decode(li));
					}
			}
		catch(Throwable err) {
			err.printStackTrace();
			throw new IOException(err);
		}
		Assert.assertEquals(L1.size(), L2.size());
		for(int i=0;i< L1.size();i++) {
			VariantContext f1 = L1.get(i);
			VariantContext f2 = L2.get(i);
			compareVariantContext(f1, f2);
			}
		
	}
	
	@Test(dataProvider="src1")
	public void readLoc(String bcffname,String vcfname) throws IOException {
		
		List<Feature> L1=new ArrayList<>();
		BCFCodec codec=new BCFCodec();
		try(BlockCompressedInputStream in=new BlockCompressedInputStream(Paths.get(bcffname))) {
			try(PositionalBufferedStream pb=new PositionalBufferedStream(in)) {
				codec.readHeader(pb);
				while(!pb.isDone()) {
					L1.add(codec.decodeLoc(pb));
					}
				}
			}
		catch(Throwable err) {
			err.printStackTrace();
			throw new IOException(err);
		}
		
		List<Feature> L2=new ArrayList<>();
		VCFCodec codec2=new VCFCodec();
		try(PositionalBufferedStream pb=new PositionalBufferedStream(Files.newInputStream(Paths.get(vcfname)))) {
			LineIterator li= new LineIteratorImpl(AsciiLineReader.from(pb)); 
				codec2.readActualHeader(li);
				while(li.hasNext()) {
					L2.add(codec2.decodeLoc(li));
					}
			}
		catch(Throwable err) {
			err.printStackTrace();
			throw new IOException(err);
		}
		
		Assert.assertEquals(L1.size(), L2.size());
		for(int i=0;i< L1.size();i++) {
			Feature f1 = L1.get(i);
			Feature f2 = L2.get(i);
			compareLoc(f1, f2);
			}
	}
}
