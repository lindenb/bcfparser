A BCF 2.2 parser using HTSJDK.

# Motivation

The developpment for HTSJDK/vcf is stalled and I really need a BCF reader.
At the time of writing HTSJDK only support BCF2.1.
See (2016) https://github.com/samtools/htsjdk/issues/596

All the data is decoded, all the genotypes are decoded on the fly, there is no soft mecanism to decode the genotypes on demand. 

So far, the tests I ran with a few bcfs show no differences with the VCF+Htsjdk (The classes of the attributes might change though, e.g "1.5" as String vs "1.5" as Float)

# Usage

```
try(BCFFileReader r=new  BCFFileReader(path,true)) {
        VCFHeader h= r.getHeader();
        try(CloseableIterator<VariantContext> iter= r.query("chr1",100,200)) {
                while(iter.hasNext()) {
                        VariantContext ctx = iter.next();
                        }
                }
        }
```

```
try(VCFIterator iter = BCFIterator.open(pat)) {
        VCFHeader h=  iter.getHeader();
        while(iter.hasNext()) {
                VariantContext ctx = iter.next();
                }
        }

```

# Compilation

Compilation of the library `bcf.jar` requires a java compiler as well as the jar libraries for htsjdk , jcommander and testng
Yeah, I'm lazy, feel free to add a build.gradle.

```
make HTSJDK=/path/to/htsjdk.jar TESTNG=/path/to/jcommander.jar:/path/to/testng.jar

```

# Author

Pierre Lindenbaum PhD
Institut du Thorax
44000 Nantes
France
