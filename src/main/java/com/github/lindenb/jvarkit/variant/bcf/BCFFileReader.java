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

import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.BAMFileSpan;
import htsjdk.samtools.BAMIndexMetaData;
import htsjdk.samtools.CSIIndex;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;
/****************
 *
 * Implementation of a {@link VCFReader} for BCF
 * @author Pierre Lindenbaum
 */
public class BCFFileReader implements VCFReader {
private static final Log LOG=Log.getInstance(BCFFileReader.class);

private final VCFHeader header;
private final BCFCodec codec;
private CSIIndex index = null;
private BaseIterator currentIterator=null;

private BCFFileReader(final BCFCodec codec)  throws IOException  {
	this.codec= codec;
	this.header=(VCFHeader)this.codec.readHeader();
	}
/**
 * Open a new BCFFileReader
 * @param file the bcf file
 * @param requireIndex shall we load an index 
 * @throws IOException on I/O error
 */
public	BCFFileReader(final Path file,boolean requireIndex) throws IOException {
	this(BCFCodec.open(file.toString()));
	
	if(requireIndex) {
		final Path csiPath=findIndex(file);
		if(csiPath==null) throw new FileNotFoundException("cannot find csi index for "+file);
		this.index = new CSIIndex(csiPath, this.header.getSequenceDictionary());
		}
	else
		{
		this.index=null;
		}
	}

/** return the chromosomes carrying variants for this BCF file 
 * @throws IllegalStateException if there is no index
 * @return  the chromosomes carrying variants for this BCF file
 */
public List<String> getChromosomes() {
	if(this.index==null) throw new IllegalStateException("Cannot invoke getChromosomes() if there is no associated index");
	final SAMSequenceDictionary dict = this.header.getSequenceDictionary();
	final List<String> chroms = new ArrayList<>();
	for(int i=0;i< dict.size();i++) {
		final BAMIndexMetaData meta= this.index.getMetaData(i);
		if(meta==null || meta.getAlignedRecordCount()==0) continue;
		chroms.add(dict.getSequence(i).getContig());
		}
	return chroms;
	}

@Override
public CloseableIterator<VariantContext> iterator() {
	if(BCFFileReader.this.currentIterator!=null) throw new IllegalStateException("already iterating");
	BCFFileReader.this.currentIterator = new MyIterator();
	return BCFFileReader.this.currentIterator;
	}

@Override
public final CloseableIterator<VariantContext> query(String chrom, int start, int end) {
	return query(new Interval(chrom,start,end));
	}
@Override
public CloseableIterator<VariantContext> query(Locatable loc) {
	if(this.index==null) throw new IllegalStateException("no index is available");
	if(loc==null) throw new IllegalArgumentException("loc is null");
	if(BCFFileReader.this.currentIterator!=null) throw new IllegalStateException("already iterating");
    int tid = this.getHeader().getSequenceDictionary().getSequenceIndex(loc.getContig());
    LOG.debug("TID="+tid+" for "+loc);
    if(tid==-1) return new EmptyIterator();
	final BAMFileSpan span = this.index.getSpanOverlapping(tid,loc.getStart(),loc.getEnd());
    long[] array=span.toCoordinateArray();
    BCFFileReader.this.currentIterator =  new MyQueryIterator(loc,array);
    return BCFFileReader.this.currentIterator;
	}


private abstract class BaseIterator extends AbstractIterator<VariantContext>
implements CloseableIterator<VariantContext>
	{
	
	VariantContext getNextRecord() throws IOException {
		return BCFFileReader.this.codec.decode();
		}
	@Override
	public void close() {
		BCFFileReader.this.currentIterator=null;
		}	
	}

private class EmptyIterator extends BaseIterator {
	@Override
	protected VariantContext advance() {
		return null;
		}	
	}

private class MyIterator extends BaseIterator
	{
	MyIterator() {
		
		try {
			BCFFileReader.this.codec.rewind();
		} catch (IOException e) {
			throw new RuntimeIOException(e);
		}
	
		}
	
	@Override
	protected VariantContext advance() {
		try {
			return getNextRecord();
			}
		catch(IOException err) {
			throw new RuntimeIOException(err);
			}
		}	
	}

private class MyQueryIterator extends BaseIterator
	{
	final Locatable interval;
	final long[] mFilePointers;
    private int mFilePointerIndex = 0;
    private long mFilePointerLimit = -1;

	MyQueryIterator(Locatable interval,long[] mFilePointers) {
		this.interval=interval;
		this.mFilePointers=mFilePointers;
		}
	@Override
	protected VariantContext advance() {
		try {
		for(;;) {
			 // Advance to next file block if necessary
	        while (0 >= mFilePointerLimit) {
	            if (mFilePointers == null ||
	                    mFilePointerIndex >= mFilePointers.length) {
	                return null;
	            }
	            final long startOffset = mFilePointers[mFilePointerIndex++];
	            final long endOffset = mFilePointers[mFilePointerIndex++];
	            BCFFileReader.this.codec.seek(startOffset);
	            mFilePointerLimit = endOffset;
	        	}
	        
	      
	        
	        // Pull next record from stream
	        VariantContext ctx= getNextRecord();
	        if(ctx==null) return null;
	        if(!ctx.overlaps(this.interval)) continue;
	        return ctx;
	        }
			}
		catch(IOException err ) {
			throw new RuntimeIOException(err);
			}	
		}
	}

@Override
public VCFHeader getHeader() {
	return this.header;
	}
@Override
public boolean isQueryable() {
	return this.index!=null;
	}

@Override
public void close() throws IOException {
	if(currentIterator!=null) currentIterator.close();
	currentIterator=null;
	this.codec.close();
	}


private static Path findIndex(final Path vcfPath0)  {
	for(int side=0;side <2;++side) {
		final Path vcfPath;
		try {
			vcfPath = (side==0?vcfPath0:vcfPath0.toRealPath());
			}
		catch(IOException err) {
			continue;
			}
		final String fileName = vcfPath.getFileName().toString();
		LOG.debug(fileName);
		final String csi = fileName+ FileExtensions.CSI;
		LOG.debug(csi);
 		final Path indexPath = vcfPath.resolveSibling(csi);
		LOG.debug(indexPath);

        if (Files.isRegularFile(indexPath)) {
           return indexPath;
           }
		}
	return null;
	}
}
