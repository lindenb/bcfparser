package com.github.lindenb.jvarkit.variant.bcf;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class BCFFileReaderTest {
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][] {
			{"testdata/test01.bcf","testdata/test01.vcf"}
			};
		}
	@Test(dataProvider="src1")
	public void testNoIndex(String bcffname,String vcfname) throws IOException{
		List<VariantContext> L1=new ArrayList<>();
		try(BCFFileReader reader=new BCFFileReader(Paths.get(bcffname), false)) {
			final VCFHeader h= reader.getHeader();
			Assert.assertNotNull(h);
			try(CloseableIterator<VariantContext> iter=reader.iterator()) {
				while(iter.hasNext()) {
					final VariantContext ctx=iter.next();
					Assert.assertNotNull(ctx);
					L1.add(ctx);
					}
				}
			}
		Assert.assertFalse(L1.isEmpty());
		}
	
	@Test(dataProvider="src1")
	public void testWithIndex(String bcffname,String vcfname) throws IOException{
		try(BCFFileReader reader=new BCFFileReader(Paths.get(bcffname), true)) {
			final VCFHeader h= reader.getHeader();
			Assert.assertNotNull(h);
		}
	}
}
