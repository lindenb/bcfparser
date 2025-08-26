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
import java.io.InputStream;
import java.nio.file.Path;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;

/**
 * Implementation of {@link VCFIterator} for BCF
 * @author Pierre Lindenbaum
 *
 */
public class BCFIterator extends AbstractIterator<VariantContext> implements VCFIterator{
private static final Log LOG=Log.getInstance(BCFIterator.class);
private final BCFCodec codec;
private final VCFHeader header;
BCFIterator(final BCFCodec codec) throws IOException {
	this.codec = codec;
	this.header = codec.readHeader();
	}

/** 
 * open a BCFIterator from an InputStream
 * @param in the input stream
 * @return the new BCFIterator
 * @throws IOException on I/O error
 */
public static BCFIterator open(final InputStream in) throws IOException {
	return new BCFIterator(BCFCodec.open(in));
	}

/** open a VCFIterator from a PATH
 * 
 * @param fname the path 
 * @throws IOException on I/O error
 * @return the new BCFIterator
 */
public static BCFIterator open(Path fname) throws IOException {
	return new BCFIterator(BCFCodec.open(fname.toString()));
	}

/** 
 * open a BCFIterator from a file name
 * @param fname the input file
 * @throws IOException on I/O error
 * @return the new BCFIterator
 */
public static BCFIterator open(String fname) throws IOException {
	return new BCFIterator(BCFCodec.open(fname));
	}

@Override
protected VariantContext advance() {
	try {
		return this.codec.decode();
	} catch (IOException e) {
		LOG.error(e, "Cannot advance");
		throw new RuntimeIOException(e);
		}
	}
@Override
public VCFHeader getHeader() {
	return this.header;
	}

@Override
public void close() {
	this.codec.close();
	}
}
