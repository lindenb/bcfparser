package com.github.lindenb.jvarkit.variant.bcf;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import htsjdk.samtools.util.BinaryCodec;

class BCFTypedData {
private static final int MAX_LENGTH=15;
private static final float bcf_float_missing    = 0x7F800001;
private static final float bcf_float_vector_end = 0x7F800002;
private static final int  bcf_int8_vector_end  = Byte.MIN_VALUE+1;
private static final int  bcf_int16_vector_end = Short.MIN_VALUE+1;
private static final int  bcf_int32_vector_end = Integer.MIN_VALUE+1;
private static final int  bcf_int8_missing  = Byte.MIN_VALUE;
private static final int  bcf_int16_missing = Short.MIN_VALUE;
private static final int  bcf_int32_missing = Integer.MIN_VALUE;
private static final byte  bcf_char_vector_end = (byte)0;



public static enum Type {
	MISSING(0,null,null),
	INT8(1,bcf_int8_missing,bcf_int8_vector_end),
	INT16(2,bcf_int16_missing,bcf_int16_vector_end),
	INT32(3,bcf_int32_missing,bcf_int32_vector_end),
	FLOAT(5,bcf_float_missing,bcf_float_vector_end),
	CHAR(7,null,null);
	
	final int opcode;
	final Object missing_value;
	final Object end_vector;
	Type(int opcode,Object missing_value,Object end_vector) {
		this.opcode=opcode;
		this.missing_value=missing_value;
		this.end_vector=end_vector;
		}
	public Object getMissing() { return this.missing_value;}
	public Object getEndVector() { return this.end_vector;}
	}
private final Type type;
private final int count;
private final Object value;

public BCFTypedData(Type type,int count,Object o) {
	this.type = type;
	this.count = count;
	this.value = o;
	}

public Type getType() {
	return type;
	}	
public int getCount() {
	return count;
	}

public Object getValue() {
	return this.value;
	}
public boolean isNull() {
	return this.value==null;
	}
public boolean isList() {
	return this.value!=null && this.value instanceof List;
	}

public int intValue() {
	if(this.value==null || !(value instanceof Integer)) {
		throw new IllegalArgumentException("not and integer "+this.toString());
		}
	return Integer.class.cast(this.value);
	}

	@Override
	public String toString() {
		String s= "{type="+type.name()+",count="+this.count+",value=";
		if(value==null)
			{
			s+="null";
			}
		else if(value instanceof List)
			{
			s+="[";
			List<?> L=(List<?>)this.value;
			s+= L.stream().map(O->O.toString()).collect(Collectors.joining(","));
			s+="]";
			}
		else
			{
			s+=String.valueOf(value);
			}
		s+="}";
		return s;
		}

public static String readString(BinaryCodec bc,int length) {
	int i;
	byte[] a=new byte[length];
	bc.readBytes(a);
	for(i=0;i< length;i++) {
		if(a[i]==bcf_char_vector_end) break;
		}
	return new String(a,0,i);
	}	
	
public static String readString(BinaryCodec bc) {
	final BCFTypedData td = read(bc);
	if(td.type!=Type.CHAR) throw new IllegalStateException("expected CHAR but got "+td.type.name());
	return String.class.cast(td.value);
	}

public static int[] readIntArray(BinaryCodec bc) {
	BCFTypedData td = read(bc);
	switch(td.type) {
		case MISSING: return new int[0];
		case INT8:case INT16: case INT32:break;//ok
		default:throw new IllegalStateException("expected INT TYPE but got "+td.type.name());
		}
	if(td.isNull()) {
		return new int[0];
		}
	else if(td.isList()) {
		List<?> L=(List<?>)td.value;
		int a[]=new int[L.size()];
		for(int i=0;i< L.size();i++) {
			a[i]=Integer.class.cast(L.get(i)).intValue();
			}
		return a;
		}
	else
		{
		return new int[] {Integer.class.cast(td.value).intValue()};
		}
	}

public static Object readAtomic(BinaryCodec bc,Type t) {
	switch(t) {
		case MISSING : return null;
		case INT8: return Integer.valueOf(bc.readByte());
		case INT16: return Integer.valueOf(bc.readShort());
		case INT32: return Integer.valueOf(bc.readInt());
		case FLOAT: return Float.valueOf(bc.readFloat());
		case CHAR: return String.valueOf((char)bc.readByte());
		default: throw new IllegalArgumentException("not int type:"+t);
	}
}

public static BCFTypedData read(BinaryCodec bc) {
	final byte b = bc.readByte();
	return read(bc,b);
	}

public static int decodeCount(BinaryCodec bc,byte b) {
	final int count0 = decodeSize(b);
	final int count;
	if(count0 >= MAX_LENGTH) {
		count = read(bc).intValue();
		}
	else
		{
		count = count0;
		}
	return count;
	}

public static BCFTypedData read(BinaryCodec bc,byte b) {
	System.err.println("byte ="+(int)b);
	final Type t = decodeType(b);
	System.err.println("type ="+t);
	final int count = decodeCount(bc,b);
	System.err.println("count ="+count);
	final Object o;
	switch(t) {
		case MISSING: {
			return new BCFTypedData(t,decodeSize(b),null);
			}
		case INT8:
			{
			if(count==0) {
				o = null;
				}
			else if(count==1) {
				o= Integer.valueOf(bc.readByte());
				}
			else
				{
				final List<Integer> L=new ArrayList<>(count);
				for(int i=0;i< count;i++) {
					Integer v= Integer.valueOf(bc.readByte());
					if(v.intValue()==bcf_int8_vector_end) break;
					if(v.intValue()==bcf_int8_missing) v=null;
					L.add(v);
					}	
				o = L;
				}
			return new BCFTypedData(t,count,o);
			}
		case INT16:
			{
			if(count==0) {
				o = null;
				}
			else if(count==1) {
				o= Integer.valueOf(bc.readShort());
				}
			else
				{
				final List<Integer> L=new ArrayList<>(count);
				for(int i=0;i< count;i++) {
					Integer v= Integer.valueOf(bc.readShort());
					if(v.intValue()==bcf_int16_vector_end) break;
					if(v.intValue()==bcf_int16_missing) v=null;
					L.add(v);
					}	
				o = L;
				}
			return new BCFTypedData(t,count,o);
			}
		case INT32:
			{
			if(count==0) {
				o = null;
				}
			else if(count==1) {
				o= Integer.valueOf(bc.readInt());
				}
			else
				{
				final List<Integer> L=new ArrayList<>(count);
				for(int i=0;i< count;i++) {
					Integer v= Integer.valueOf(bc.readInt());
					if(v.intValue()==bcf_int32_vector_end) break;
					if(v.intValue()==bcf_int32_missing) v=null;
					L.add(v);
					}	
				o = L;
				}
			return new BCFTypedData(t,count,o);
			}
		case FLOAT:
			{
				if(count==0) {
					o = null;
					}
				else if(count==1) {
					o= Float.valueOf(bc.readFloat());
					}
				else
					{
					final List<Float> L=new ArrayList<>(count);
					for(int i=0;i< count;i++) {
						Float v= Float.valueOf(bc.readFloat());
						if(v.floatValue()==bcf_float_vector_end) break;
						if(v.floatValue()==bcf_float_missing) v=null;
						L.add(v);
						}	
					o = L;
					}
				return new BCFTypedData(t,count,o);
				}
		case CHAR:
			{
			o = readString(bc,count);
			return new BCFTypedData(t,count,o);
			}
		default: throw new IllegalArgumentException("undefined type:"+t);
		}
	}

private static Float readFloat(BinaryCodec bc) {
	return bc.readFloat();
}

private static int decodeSize(final byte typeDescriptor) {
    return (0xF0 & typeDescriptor) >> 4;
}
static int decodeTypeID(final byte typeDescriptor) {
    return typeDescriptor & 0x0F;
}
static Type decodeType(final byte typeDescriptor) {
   final int t=decodeTypeID(typeDescriptor);
	   switch(t) {
	   case 0: return Type.MISSING;
	   case 1: return Type.INT8;
	   case 2 : return Type.INT16;
	   case 3: return Type.INT32;
	   case 5: return Type.FLOAT;
	   case 7: return Type.CHAR;
	   default: throw new IllegalArgumentException("undefined typeDescriptor:"+t);
	   }
	}
}
