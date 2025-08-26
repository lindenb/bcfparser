package com.github.lindenb.jvarkit.variant.bcf;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class BCFIteratorTest {
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][] {
			{"testdata/test01.bcf","testdata/test01.vcf"}
			};
		}
	@Test(dataProvider="src1")
	public void read(String bcffname,String vcfname) throws IOException{
		List<VariantContext> L1=new ArrayList<>();
		try(BCFIterator iter=BCFIterator.open(bcffname)) {
			while(iter.hasNext()) {
				final VariantContext ctx=iter.next();
				Assert.assertNotNull(ctx);
				L1.add(ctx);
				}
		}
		
		Assert.assertFalse(L1.isEmpty(),"NO variant was found");

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
		Assert.assertFalse(L2.isEmpty(),"NO variant was found");
		BCFCodecTest.compareLists(L1, L2);
	}
}
