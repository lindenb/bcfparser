package com.github.lindenb.jvarkit.variant.bcf;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import htsjdk.samtools.BAMFileSpan;
import htsjdk.samtools.CSIIndex;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;
/****************
 *
 * 
 *
 */
public class BCFFileReader implements VCFReader {
private final VCFHeader header;
private final BCFCodec codec;
private CSIIndex index = null;
private BaseIterator currentIterator=null;

private BCFFileReader(BCFCodec codec)  throws IOException  {
	this.codec= codec;;
	this.header=(VCFHeader)this.codec.readHeader();
	}

public	BCFFileReader(Path file,boolean requireIndex) throws IOException {
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



@Override
public CloseableIterator<VariantContext> iterator() {
	if(BCFFileReader.this.currentIterator!=null) throw new IllegalStateException("already iterating");
	BCFFileReader.this.currentIterator = new MyIterator();
	return BCFFileReader.this.currentIterator;
	}

public final CloseableIterator<VariantContext> query(String chrom, int start, int end) {
	return query(new Interval(chrom,start,end));
	}
@Override
	public CloseableIterator<VariantContext> query(Locatable loc) {
	if(BCFFileReader.this.currentIterator!=null) throw new IllegalStateException("already iterating");
    int tid = this.getHeader().getSequenceDictionary().getSequenceIndex(loc.getContig());
    if(tid==-1) return new EmptyIterator();
	final BAMFileSpan span = index.getSpanOverlapping(tid,loc.getStart(),loc.getEnd());
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
	        return getNextRecord();
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
		System.err.println(fileName);
		final String csi = fileName+ FileExtensions.CSI;
		System.err.println(csi);
 		final Path indexPath = vcfPath.resolveSibling(csi);
		System.err.println(indexPath);

        if (Files.isRegularFile(indexPath)) {
           return indexPath;
           }
		}
	return null;
	}
}
