package sketch;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.atomic.AtomicInteger;

import align2.Shared;
import align2.Tools;
import dna.Parser;
import dna.Timer;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import kmer.HashArray1D;
import kmer.HashForest;
import kmer.KmerNode;
import kmer.KmerTableSet;
import kmer.Primes;
import structures.LongHeap;
import structures.LongHeapSet;
import structures.LongList;

/**
 * @author Brian Bushnell
 * @date June 28, 2016
 *
 */
public final class SketchTool {
	
	/*--------------------------------------------------------------*/
	/*----------------         Main Method          ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			//printOptions();
			System.exit(0);
		}
		
		Timer t=new Timer();
		t.start();
		
		//Create a new CountKmersExact instance
		SketchTool mhs=new SketchTool(args);
		t.stop();
		System.err.println("Time: \t"+t);
	}
	
	public SketchTool(String[] args){
		System.err.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		/* Set global defaults */
		ReadWrite.ZIPLEVEL=2;
		ReadWrite.USE_UNPIGZ=true;
		
		/* Initialize local variables with defaults */
		Parser parser=new Parser();
		
		ArrayList<String> list=new ArrayList<String>();
		
		int k_=31;
		int size_=10000;
		int mincount_=1;
		int bitArrayBits_=0;
		float cutoff=0.02f;
		
		/* Parse arguments */
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if("null".equalsIgnoreCase(b)){b=null;}
			while(a.charAt(0)=='-' && (a.indexOf('.')<0 || i>1 || !new File(a).exists())){a=a.substring(1);}
			
			if(a.equals("in")){
				if(b!=null){
					for(String s : b.split(",")){
						list.add(s);
					}
				}
			}else if(a.equals("k")){
				k_=Integer.parseInt(b);
			}else if(a.equals("size") || a.equals("length") || a.equals("len")){
				size_=Integer.parseInt(b);
			}else if(a.equals("mincount")){
				mincount_=Integer.parseInt(b);
			}else if(a.equals("cutoff")){
				cutoff=Float.parseFloat(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(b==null){
				list.add(arg);
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		prime1=-1;
		prime2=-1;
		k=k_;
		size=size_;
		mincount=mincount_;
		if(bitArrayBits_<1){
			bitArrayBits_=(int)Primes.primeAtLeast(size*3);
		}
		bitArrayBits=bitArrayBits_;
		
		Timer t=new Timer();
		ArrayList<Sketch> sketches=loadSketches_MT(list);
		t.stop();
		System.err.println("Loaded "+sketches.size()+" sketches in \t"+t);
		t.start();
		Sketch sketch=sketches.get(0);
		for(int i=1; i<sketches.size(); i++){
			Sketch sketch2=sketches.get(i);
			float identity=sketch.identity(sketch2);
//			System.err.println(sketch+"\n"+sketch2);
			if(identity>=cutoff){
				System.out.println(String.format("%.2f%%", 100*identity)+" identity for "+sketch.name+" vs "+sketch2.name);
			}
		}
		t.stop();
		System.err.println("Compared "+(sketches.size()-1)+" sketches in \t"+t);
		
//		String fname=list.get(0);
//		Sketch sketch=loadSketches(fname).get(0);
//		for(int i=1; i<list.size(); i++){
//			String fname2=list.get(i);
//			Sketch sketch2=loadSketches(fname2).get(0);
//			float identity=sketch.identity(sketch2);
//			System.out.println("Identity for "+fname+" vs "+fname2+":\t"+String.format("%.2f%%", 100*identity));
//		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Normal Constructor      ----------------*/
	/*--------------------------------------------------------------*/
	
	public SketchTool(int size_, int k_, int mincount_){
		k=k_;
		size=size_;
		mincount=mincount_;
		prime1=Primes.primeAtLeast(1L<<k);
		prime2=Primes.primeAtMost((long)(prime1*0.9f));
		bitArrayBits=(int)Primes.primeAtLeast(size*3);

		assert(k>0 && k<32) : "Sketches require 0 < K < 32.";
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	public Sketch toSketch(KmerTableSet tables, boolean multithreaded){
		final int threads=(multithreaded ? Tools.mid(1, Shared.threads(), tables.ways()) : 1);
		return (threads<2 ? toSketch_ST(tables) : toSketch_MT(tables, threads));
	}
	
	private Sketch toSketch_ST(KmerTableSet tables){
		LongHeapSet heap=new LongHeapSet(size);

		KmerTableSet kts=(KmerTableSet)tables;
		for(int tnum=0; tnum<kts.ways; tnum++){
			HashArray1D table=kts.getTable(tnum);
			toHeap(table, heap);
		}
		
		return new Sketch(heap);
	}
	
	private Sketch toSketch_MT(KmerTableSet tables, final int threads){
		ArrayList<SketchThread> alst=new ArrayList<SketchThread>(threads);
		AtomicInteger ai=new AtomicInteger(0);
		for(int i=0; i<threads; i++){
			alst.add(new SketchThread(ai, tables));
		}

		//Start the threads
		for(SketchThread pt : alst){
			pt.start();
		}

		ArrayList<LongHeapSet> heaps=new ArrayList<LongHeapSet>(threads);
		for(SketchThread pt : alst){

			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			if(pt.heap.size()>0){
				heaps.add(pt.heap);
			}
		}
		alst.clear();
		return toSketch(heaps);
	}
	
//	public BitSet toBinarySketch(long[] sketch){
//		BitSet bs=new BitSet();
//		for(long value : sketch){
//			value&=MASK;
//			int idx=(int)(value%bitArrayBits);
//			bs.set(idx);
//		}
//		return bs;
//	}
	
	public static final long[] toSketchArray(LongHeap heap){
		long[] array=new long[heap.size()];
		for(int i=0; i<array.length; i++){
			array[i]=Long.MAX_VALUE-heap.poll();
		}
		Tools.reverseInPlace(array);
		assert(heap.size()==0);
		return array;
	}
	
	public LongHeapSet toHeap(HashArray1D table, LongHeapSet heap){
//		if(heap==null){heap=new LongHeap(size, true);}
		long[] kmers=table.array();
		int[] counts=table.values();
		for(int i=0; i<table.arrayLength(); i++){
			int count=counts[i];
			if(count>=mincount){
				long hash=hash(kmers[i]);
				heap.add(hash);
			}
		}
		HashForest forest=table.victims();
		if(forest!=null){
			for(KmerNode kn : forest.array()){
				if(kn!=null){addRecursive(heap, kn);}
			}
		}
		return heap;
	}
	
//	public long[] toSketchArray(ArrayList<LongHeap> heaps){
//		if(heaps.size()==1){return toSketchArray(heaps.get(0));}
//		LongList list=new LongList(size);
//		for(LongHeap heap : heaps){
//			while(heap.size()>0){list.add(Long.MAX_VALUE-heap.poll());}
//		}
//		list.sort();
//		list.shrinkToUnique();
//		list.size=Tools.min(size, list.size);
//		return list.toArray();
//	}
	
	public Sketch toSketch(ArrayList<LongHeapSet> heaps){
		LongHeapSet a=heaps.get(0);
		for(int i=1; i<heaps.size(); i++){
			LongHeapSet b=heaps.get(i);
			a.add(b);
		}
		return new Sketch(a);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/
	
	private void addRecursive(LongHeapSet heap, KmerNode kn){
		if(kn==null){return;}
		if(kn.count()>=mincount){
			long hash=hash(kn.pivot());
			heap.add(hash);
		}
		if(kn.left()!=null){addRecursive(heap, kn.left());}
		if(kn.right()!=null){addRecursive(heap, kn.right());}
	}
	
	public static long parseHex(byte[] line){
		if(line.length==0){return 0;}
		long x=0;
		for(byte b : line){
			x<<=4;
			x|=hexTable[b];
		}
		if(line[0]=='-'){x*=-1;}
		return x;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             I/O              ----------------*/
	/*--------------------------------------------------------------*/
	
	public static ArrayList<Sketch> loadSketches_ST(String...fnames){
		ArrayList<Sketch> sketches=null;
		for(String s : fnames){
			ArrayList<Sketch> temp;
			if(s.indexOf(',')<0 || s.startsWith("stdin") || new File(s).exists()){
				temp=loadSketches(s);
			}else{
				temp=loadSketches_ST(s.split(","));
			}
			if(sketches==null){sketches=temp;}
			else{sketches.addAll(temp);}
		}
		return sketches;
	}
	
	public static ArrayList<Sketch> loadSketches_MT(ArrayList<String> fnames){
		return loadSketches_MT(fnames.toArray(new String[0]));
	}
	
	public static ArrayList<Sketch> loadSketches_MT(String...fnames){
		ArrayList<String> decomposedFnames=new ArrayList<String>(fnames.length);
		for(String s : fnames){
			if(s.indexOf(',')<0 || s.startsWith("stdin") || new File(s).exists()){
				decomposedFnames.add(s);
			}else{
				for(String s2 : s.split(",")){
					decomposedFnames.add(s2);
				}
			}
		}

		if(decomposedFnames.size()==0){return null;}
		if(decomposedFnames.size()==1){return loadSketches(decomposedFnames.get(0));}
		
		
		//Determine how many threads may be used
		final int threads=Tools.min(Shared.threads(), decomposedFnames.size());

		//Fill a list with LoadThreads
		ArrayList<LoadThread> allt=new ArrayList<LoadThread>(threads);
		
		for(int i=0; i<threads; i++){
			allt.add(new LoadThread(decomposedFnames.get(i)));
		}
		
		ArrayList<Sketch> sketches=new ArrayList<Sketch>();
		
		//Start the threads
		for(LoadThread lt : allt){lt.start();}

		//Wait for completion of all threads
		boolean success=true;
		for(LoadThread lt : allt){

			//Wait until this thread has terminated
			while(lt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					lt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			sketches.addAll(lt.list);
			success&=lt.success;
			if(!success){System.err.println("Failure loading "+lt.fname);}
		}
		assert(success) : "Failure loading some files.";
		return sketches;
	}
	
	private static class LoadThread extends Thread{
		
		public LoadThread(String fname_){
			fname=fname_;
		}
		
		public void run(){
			success=false;
			list=loadSketches(fname);
			success=true;
		}
		
		final String fname;
		ArrayList<Sketch> list;
		boolean success=false;
		
	}
	
	public static ArrayList<Sketch> loadSketches(String fname){
		ArrayList<Sketch> sketches=new ArrayList<Sketch>();
		
		ByteFile bf=ByteFile.makeByteFile(fname, false, false);
		int size=10000;
		int taxID=-1;
		String name=null;
		LongList list=null;
		long sum=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length>0){
//				System.err.println("Processing line "+new String(line));
				if(line[0]=='#'){
					if(list!=null){
						assert(list.size==list.array.length);
						list.shrink();
						Sketch sketch=new Sketch(list.array, taxID, name);
						sketches.add(sketch);
//						System.err.println("Made sketch "+sketch);
					}
					name=null;
					list=null;
					sum=0;
					taxID=-1;

					if(line.length>1){
						String[] split=new String(line, 1, line.length-1).split("\t");
						for(String s : split){
							if(s.startsWith("SIZE")){
								size=Integer.parseInt(s.substring(5));
							}else if(s.startsWith("TAXID")){
								taxID=Integer.parseInt(s.substring(6));
							}else if(s.startsWith("ID")){
								taxID=Integer.parseInt(s.substring(3));
							}else if(s.startsWith("NAME")){
								name=s.substring(5);
							}
						}
					}
					if(size>0){list=new LongList(size);}
				}else{
					long x=parseHex(line);
//					System.err.println("sum="+sum+", x="+x+" -> "+(sum+x));
					sum+=x;
					assert(x>=0) : x+"\n"+new String(line);
					assert(sum>=0) : "The sketch was made with delta compression off.  Please regenerate it.";
					list.add(delta ? sum : x);
					//						System.err.println("List="+list);
				}
			}
		}
		if(list==null){list=new LongList(0);}
		assert(list.size==list.array.length);
		list.shrink();
		Sketch sketch=new Sketch(list.array, taxID, name);
		sketches.add(sketch);
		return sketches;
	}
	
	public void write(ArrayList<Sketch> sketches, FileFormat ff[]){
		final int len=ff.length;
		TextStreamWriter tsw[]=new TextStreamWriter[len];
		for(int i=0; i<len; i++){
			tsw[i]=new TextStreamWriter(ff[i]);
			tsw[i].start();
		}
		for(int i=0; i<sketches.size(); i++){
			write(sketches.get(i), tsw[i%len]);
		}
		for(int i=0; i<len; i++){
			tsw[i].poisonAndWait();
		}
	}
	
	public void write(ArrayList<Sketch> sketches, FileFormat ff){
		TextStreamWriter tsw=new TextStreamWriter(ff);
		tsw.start();
		for(Sketch sketch : sketches){
			write(sketch, tsw);
		}
		tsw.poisonAndWait();
	}
	
	public void write(Sketch sketch, FileFormat ff){
		TextStreamWriter tsw=new TextStreamWriter(ff);
		tsw.start();
		write(sketch, tsw);
		tsw.poisonAndWait();
	}
	
	public void write(Sketch sketch, TextStreamWriter tsw){
		long prev=0;
		long[] array=sketch.array;
		tsw.print("#SIZE:"+array.length);
		if(sketch.taxID>=0){tsw.print("\tTAXID:"+sketch.taxID);}
		if(sketch.name!=null){tsw.print("\tNAME:"+sketch.name);}
		tsw.print("\n");
		for(int i=0; i<array.length; i++){
			long key=array[i];
			tsw.println(Long.toHexString(key-prev));
			if(delta){prev=key;}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Converts KmerTableSets to Heaps */
	private class SketchThread extends Thread {

		SketchThread(AtomicInteger next_, KmerTableSet kts_){
			next=next_;
			kts=kts_;
		}

		public void run(){
			final int ways=kts.ways();
			int tnum=next.getAndIncrement();
			while(tnum<ways){
				if(heap==null){heap=new LongHeapSet(size);}
				HashArray1D table=kts.getTable(tnum);
				toHeap(table, heap);
				tnum=next.getAndIncrement();
			}
		}

		final AtomicInteger next;
		final KmerTableSet kts;
		LongHeapSet heap;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Hashing            ----------------*/
	/*--------------------------------------------------------------*/
	
	public long hashOld(long kmer){
		assert(prime1>0);
		long rot=Long.rotateRight(kmer^0x5555555555555555L, 9);
		long a=rot%prime1;
		long b=rot%prime2;
		return Long.rotateRight(kmer, 17)^Long.rotateLeft(a, 31)^Long.rotateRight(b, 31);
	}
	
	private static long[][] makeCodes(int symbols, int modes){
		Random randy=new Random(12345);
		long[][] r=new long[symbols][modes];
		for(int i=0; i<symbols; i++){
			for(int j=0; j<modes; j++){
				r[i][j]=randy.nextLong();
			}
		}
		return r;
	}
	
	public static final long hash(long kmer){
		long code=kmer;
		for(int i=0; i<8; i++){
			int x=(int)(kmer&0xFF);
			code^=codes[i][x];
		}
		return code;
	}
	
	private static final long[][] codes=makeCodes(8, 256);
		
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** For hashing */
	private final long prime1, prime2;
	/** For binary representation */
	private final int bitArrayBits;
//	private final int bitArrayLen;
	private final int k;
	private final int size;
	private final int mincount;
	
	public static final boolean delta=true;
	private static final long MASK=Long.MAX_VALUE;
	
	private static final byte[] hexTable=new byte[128];
	static {
		Arrays.fill(hexTable, (byte)-1);
		for(int i='0'; i<='9'; i++){
			hexTable[i]=(byte)(i-'0');
		}
		for(int i='A'; i<='F'; i++){
			hexTable[i]=hexTable[i+'a'-'A']=(byte)(i-'A'+10);
		}
		hexTable['x']=hexTable['X']=hexTable['-']=hexTable['+']=0;
	}
	
}
