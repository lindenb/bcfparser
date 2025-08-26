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

import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class BCFIteratorTest {
	@DataProvider(name = "src1")
	public Iterator<Object[]> createData1() {
		return BCFCodecTest.listVCFsFortest();
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
