package com.github.lindenb.jvarkit.variant.bcf;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Path;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;

public class BCFIterator extends AbstractIterator<VariantContext> implements VCFIterator{
private final BCFCodec codec;
private final VCFHeader header;
BCFIterator(final BCFCodec codec) throws IOException {
	this.codec = codec;
	this.header = codec.readHeader();
	}

public static BCFIterator open(InputStream in) throws IOException {
	return new BCFIterator(BCFCodec.open(in));
	}

public static BCFIterator open(Path fname) throws IOException {
	return new BCFIterator(BCFCodec.open(fname.toString()));
	}

public static BCFIterator open(String fname) throws IOException {
	return new BCFIterator(BCFCodec.open(fname));
	}

@Override
protected VariantContext advance() {
	try {
		return this.codec.decode();
	} catch (IOException e) {
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
