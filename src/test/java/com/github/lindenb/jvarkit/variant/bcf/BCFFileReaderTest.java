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
import java.util.Arrays;
import java.util.Iterator;
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
	public Iterator<Object[]> createData1() {
		return BCFCodecTest.listVCFsFortest();
		}
	
	@DataProvider(name = "src2")
	public Iterator<Object[]>  createData2() {
		return Arrays.asList(
			new Object[]{"chr1",1,100000000},
			new Object[]{"chr5",1,1000000000}
			).iterator();
		}
	
	@DataProvider(name = "src3")
	public Object[][] createData3() {
		Iterator<Object[]> iter =createData1();
		List<Object[]> L1= new ArrayList<>();
		while(iter.hasNext()) L1.add(iter.next());
		
		iter =createData2();
		List<Object[]> L2= new ArrayList<>();
		while(iter.hasNext()) L2.add(iter.next());
		
		List<Object[]> L=new ArrayList<>();
		for(Object[] i1:L1 ) {
			for(Object[] i2:L2) {
				List<Object> L4=new ArrayList<>();
				for(Object i:i1) L4.add(i);
				for(Object i:i2) L4.add(i);
				L.add(L4.toArray(new Object[L4.size()]));
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
