package sketch;

import align2.Tools;
import structures.LongHeapSet;

/**
 * @author Brian Bushnell
 * @date July 7, 2016
 *
 */
public class Sketch implements Comparable<Sketch> {
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	public Sketch(long[] array_){
		this(array_, -1, null);
	}
	
	public Sketch(LongHeapSet heap){
		this(SketchTool.toSketchArray(heap.heap), heap.id, heap.name);
	}
	
	public Sketch(long[] array_, int taxID_, String name_){
		array=array_;
		taxID=taxID_;
		name=name_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	public void toBinary(final int bits){
		assert(binaryArray==null);
		binaryCardinality=0;
		final int len=(bits+63)/64;
		binaryArray=new long[len];
		for(long x : array){
			final int bitIndex=(int)(x%bits);
			final int index=bitIndex/64;
			final int shift=bitIndex-64*index;
			binaryArray[index]|=(1L<<shift);
		}
		for(long x : binaryArray){binaryCardinality+=Long.bitCount(x);}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Comparison          ----------------*/
	/*--------------------------------------------------------------*/
	
	public float identity(Sketch b){
		return identity(array, b.array);
	}
	
	public static float identity(long[] a, long[] b){
		int matches=countMatches(a, b);
		return matches/(float)(Tools.max(1, Tools.min(a.length, b.length)));
	}
	
	public float identityBinary(Sketch b){
		long matches=countMatchesBinary(binaryArray, b.binaryArray);
		return matches/(float)(Tools.max(1, Tools.min(binaryCardinality, b.binaryCardinality)));
	}
	
	public static int countMatches(long[] a, long[] b){
		int matches=0;
		for(int i=0, j=0; i<a.length && j<b.length; ){
			final long ka=a[i], kb=b[j];
			if(ka==kb){
				matches++;
				i++;
				j++;
			}else if(ka<kb){
				i++;
			}else{
				j++;
			}
		}
		return matches;
	}
	
	public static long countMatchesBinary(long[] a, long[] b){
		long matches=0;
		for(int i=0, j=0; i<a.length && j<b.length; ){
			final long ka=a[i], kb=b[j];
			matches+=Long.bitCount(ka&kb);
		}
		return matches;
	}
	
	@Override
	public int hashCode(){
		return taxID>=0 ? taxID : name!=null ? name.hashCode() : super.hashCode();
	}
	
	@Override
	public int compareTo(Sketch b){
		if(this==b){return 0;}
		if(taxID>-1 && b.taxID>-1){return taxID-b.taxID;}
		return name.compareTo(b.name);
	}
	
	@Override
	public boolean equals(Object b){
		if(this==b){return true;}
		if(b==null || this.getClass()!=b.getClass()){return false;}
		return equals((Sketch)b);
	}
	
	public boolean equals(Sketch b){
		return compareTo(b)==0;
	}
	
	public String toString(){
		long prev=0;
		StringBuilder sb=new StringBuilder();
		sb.append("#SIZE:"+array.length);
		if(taxID>=0){sb.append("\tTAXID:"+taxID);}
		if(name!=null){sb.append("\tNAME:"+name);}
		sb.append("\n");
		for(int i=0; i<array.length; i++){
			long key=array[i];
			sb.append(Long.toHexString(key-prev)).append('\n');
			//if(delta){prev=key;}
		}
		return sb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public long[] array;
	public int taxID;
	public String name;
	
	public long[] binaryArray;
	public long binaryCardinality;
	
	public static final boolean delta=true;
}
