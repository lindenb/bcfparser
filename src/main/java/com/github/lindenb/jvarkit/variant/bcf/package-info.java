/**

This package is a custom implemenation of a BCF 2.2 java reader

Usage:

<pre>
try(BCFFileReader r=new  BCFFileReader(path,true)) {
	VCFHeader h= r.getHeader();
	try(CloseableIterator&lt;VariantContext&gt; iter= r.query("chr1",100,200)) {
		while(iter.hasNext()) {
			VariantContext ctx = iter.next();
			}
		}
	}
</pre>


<pre>
try(VCFIterator iter = BCFIterator.open(pat)) {
	VCFHeader h=  iter.getHeader();
	while(iter.hasNext()) {
		VariantContext ctx = iter.next();
		}
	}
</pre>

@author Pierre Lindenbaum

*/
package com.github.lindenb.jvarkit.variant.bcf;
