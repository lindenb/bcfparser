package com.github.lindenb.jvarkit.variant.bcf;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;

public class BCFCodecTest {
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][] {
			{"testdata/test01.bcf","testdata/test01.vcf"}
			};
		}
	
	static List<VariantContext> readPlainVCF(final String vcfname) {
		List<VariantContext> L2=new ArrayList<>();
		try(VCFFileReader r=new VCFFileReader(Paths.get(vcfname),false)) {
			try(CloseableIterator<VariantContext> iter=r.iterator()) {
				while(iter.hasNext()) {
					final VariantContext ctx=iter.next();
					Assert.assertNotNull(ctx);
					L2.add(ctx);
					}
				}
			}
		return L2;
		}
	
	static void compareLoc(Locatable f1,Locatable f2) {
		Assert.assertEquals(f1.getContig(), f2.getContig());
		Assert.assertEquals(f1.getStart(), f2.getStart());
		Assert.assertEquals(f1.getEnd(), f2.getEnd());
	}

	static void compareLists(List<VariantContext> L1,List<VariantContext> L2) {
		Assert.assertEquals(L1.size(), L2.size());
		for(int i=0;i< L1.size();++i) {
			compareVariantContext(L1.get(i),L2.get(i));
		}
	}
	
	static void compareGenotypes(Genotype g1,Genotype g2) {
		Assert.assertEquals(g1.getSampleName(), g2.getSampleName());
		Assert.assertEquals(g1.getPloidy(), g2.getPloidy());
		Assert.assertEquals(g1.hasDP(), g2.hasDP());
		Assert.assertEquals(g1.hasGQ(), g2.hasGQ());
		Assert.assertEquals(g1.hasAD(), g2.hasAD());
		Assert.assertEquals(g1.hasPL(), g2.hasPL());
		Assert.assertEquals(g1.isPhased(), g2.isPhased(),""+g1.getGenotypeString()+" "+g2.getGenotypeString());
		Assert.assertEquals(g1.isFiltered(), g2.isFiltered(),
				"was:"+g1.getSampleName()+":("+g1.getFilters()+":"+g1.isFiltered()+") ("+g2.getFilters()+":"+g2.isFiltered()+")"
				);
		if(g1.hasDP()) {
			Assert.assertEquals(g1.getDP(), g2.getDP());
			}
		if(g1.hasGQ()) {
			Assert.assertEquals(g1.getGQ(), g2.getGQ());
			}
		if(g1.hasAD()) {
			Assert.assertTrue(Arrays.equals(g1.getAD(), g2.getAD()));
			}
		if(g1.hasPL()) {
			Assert.assertTrue(Arrays.equals(g1.getPL(), g2.getPL()));
			}
		Assert.assertTrue(g1.sameGenotype(g2, false));
		Map<String,Object> m1= g1.getExtendedAttributes();
		Map<String,Object> m2= g2.getExtendedAttributes();
		Assert.assertEquals(m1.keySet(), m2.keySet());
		for(final String k: m1.keySet()) {
			compareAttribute(k, m1.get(k), m2.get(k));
			}
		}
	
	static void compareVariantContext(VariantContext f1,VariantContext f2) {
		compareLoc(f1,f2);
		Assert.assertEquals(f1.hasID(), f2.hasID());
		Assert.assertEquals(f1.getNAlleles(), f2.getNAlleles());
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
			compareAttribute(k,o1,o2);
			}
		Assert.assertEquals(f1.getNSamples(), f2.getNSamples());
		for(int x=0;x< f1.getNSamples();++x) {
			Genotype g1 = f1.getGenotype(x);
			Genotype g2 = f2.getGenotype(x);

			compareGenotypes(g1,g2);
		}
	}
	
	static void compareAttribute(final String key,Object o1,Object o2) {
		if(o1!=null && o2!=null && !o1.getClass().equals(o2.getClass())) {
			if(o1 instanceof Float && o2 instanceof String) {
				compareAttribute(key, o1, Float.valueOf(String.class.cast(o2)));
				return;
				}
			}
		
		Assert.assertEquals(o1.toString(), o2.toString(),"not the same["+key+"] ("+o1+":"+o1.getClass()+") vs ("+o2+":"+o2.getClass()+")");
		}
	
	@Test(dataProvider="src1")
	public void read(String bcffname,String vcfname) throws IOException{
		List<VariantContext> L1=new ArrayList<>();

		try(BlockCompressedInputStream in=new BlockCompressedInputStream(Paths.get(bcffname))) {
			try(BCFCodec codec=BCFCodec.open(in)) {
				codec.readHeader();
				for(;;) {
					VariantContext ctx=codec.decode();
					if(ctx==null) break;
					L1.add(ctx);
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
		
		List<Locatable> L1=new ArrayList<>();
		try(BlockCompressedInputStream in=new BlockCompressedInputStream(Paths.get(bcffname))) {
			try(BCFCodec codec=BCFCodec.open(in)) {
				codec.readHeader();
				for(;;) {
					Locatable ctx=codec.decodeLoc();
					if(ctx==null) break;
					L1.add(ctx);
					}
				}
			}
		catch(Throwable err) {
			err.printStackTrace();
			throw new IOException(err);
		}
		
		List<Locatable> L2=new ArrayList<>();
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
			Locatable f1 = L1.get(i);
			Locatable f2 = L2.get(i);
			compareLoc(f1, f2);
			}
	}
}
