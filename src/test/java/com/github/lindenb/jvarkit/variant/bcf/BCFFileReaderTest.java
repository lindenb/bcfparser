package com.github.lindenb.jvarkit.variant.bcf;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class BCFFileReaderTest {
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][] {
			{"testdata/test01.bcf","testdata/test01.vcf"}
			};
		}
	
	@DataProvider(name = "src2")
	public Object[][] createData2() {
		return new Object[][] {
			{"chr1",1,100000000},
			{"chr5",1,1000000000}
			};
		}
	
	@DataProvider(name = "src3")
	public Object[][] createData3() {
		Object[][] o1=createData1();
		Object[][] o2=createData2();
		List<Object[]> L=new ArrayList<>();
		for(int x=0;x< o1.length;x++) {
			for(int y=0;y< o2.length;y++) {
				List<Object> L2=new ArrayList<>();
				for(Object i:o1[x]) L2.add(i);
				for(Object i:o2[y]) L2.add(i);
				L.add(L2.toArray(new Object[L2.size()]));
				}
			}
		Object[][] o3= new Object[L.size()][];
		for(int i=0;i< L.size();++i) {
			o3[i]=L.get(i);
			}
		return o3;
		}
	
	@Test(dataProvider="src1")
	public void testNoIndex(String bcffname,String vcfname) throws IOException{
		List<VariantContext> L1=new ArrayList<>();
		try(BCFFileReader reader=new BCFFileReader(Paths.get(bcffname), false)) {
			final VCFHeader h= reader.getHeader();
			Assert.assertNotNull(h);
			for(int i=0;i< 3;i++) {//rest rewind file
				L1.clear();
				try(CloseableIterator<VariantContext> iter=reader.iterator()) {
					while(iter.hasNext()) {
						final VariantContext ctx=iter.next();
						Assert.assertNotNull(ctx);
						L1.add(ctx);
						}
					}
				}
			}
		Assert.assertFalse(L1.isEmpty(),"NO variant was found");
		
		List<VariantContext> L2= BCFCodecTest.readPlainVCF(vcfname);
		 BCFCodecTest.compareLists(L1,L2);
		}
	
	@Test(dataProvider="src3")
	public void testWithIndex(String bcffname,String vcfname,String contig,int start,int end) throws IOException{
		List<VariantContext> L1=new ArrayList<>();
		Interval loc =new Interval(contig,start,end);
		try(BCFFileReader reader=new BCFFileReader(Paths.get(bcffname), true)) {
			final VCFHeader h= reader.getHeader();
			Assert.assertNotNull(h);
			try(CloseableIterator<VariantContext> iter=reader.query(loc)) {
				Assert.assertNotNull(iter);
				while(iter.hasNext()) {
					VariantContext vc=iter.next();
					Assert.assertNotNull(vc);
					Assert.assertTrue(vc.overlaps(loc),""+loc+" "+new Interval(vc));
					L1.add(vc);
				}
			}
		}
		
		List<VariantContext> L2= BCFCodecTest.readPlainVCF(vcfname).stream().filter(it->it.overlaps(loc)).collect(Collectors.toList());
		BCFCodecTest.compareLists(L1, L2);
	}
}
