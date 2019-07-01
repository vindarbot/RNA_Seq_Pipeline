package sketch;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

import align2.ReadStats;
import align2.Shared;
import align2.Tools;
import dna.AminoAcid;
import dna.Parser;
import dna.Timer;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import structures.LongHeapSet;
import tax.GiToNcbi;
import tax.TaxNode;
import tax.TaxTree;

/**
 * Creates MinHashSketches rapidly.
 * 
 * @author Brian Bushnell
 * @date July 6, 2016
 *
 */
public class SketchMaker {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		SketchMaker sm=new SketchMaker(args);
		
		//Run the object
		sm.process(t);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public SketchMaker(String[] args){
		
		//Process any config files
		args=Parser.parseConfig(args);
		
		//Detect whether the uses needs help
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		//Print the program name and arguments
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		//Set some shared static variables regarding PIGZ
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		Shared.READ_BUFFER_LENGTH=1;
		
		//Create a parser object
		Parser parser=new Parser();
		
		int size_=10000;
		int k_=31;
		int files_=1;
		boolean rcomp_=true;
		int mode_=ONE_SKETCH;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //Strip leading hyphens
			
			if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("files")){
				files_=Integer.parseInt(b);
			}else if(a.equals("size")){
				size_=(int)Tools.parseKMG(b);
			}else if(a.equals("k")){
				k_=Integer.parseInt(b);
			}else if(a.equals("rcomp")){
				rcomp_=Tools.parseBoolean(b);
			}else if(a.equals("name")){
				outName=b;
			}else if(a.equals("taxid")){
				outTaxid=Integer.parseInt(b);
			}else if(a.equals("mode")){
				if(b.equalsIgnoreCase("single")){
					mode_=ONE_SKETCH;
				}else if(b.equalsIgnoreCase("taxa")){
					mode_=PER_TAXA;
				}else if(b.equalsIgnoreCase("sequence")){
					mode_=PER_SEQUENCE;
				}
			}else if(b==null && a.equals("single")){
				mode_=ONE_SKETCH;
			}else if(b==null && a.equals("taxa")){
				mode_=PER_TAXA;
			}else if(b==null && a.equals("sequence")){
				mode_=PER_SEQUENCE;
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Tools.parseKMG(b);
				//Set a variable here
			}
			
			else if(a.equals("table") || a.equals("gi") || a.equals("gitable")){
				giTableFile=b;
				if("auto".equalsIgnoreCase(b)){giTableFile=TaxTree.DefaultTableFile;}
			}else if(a.equals("taxtree") || a.equals("tree")){
				taxTreeFile=b;
				if("auto".equalsIgnoreCase(b)){taxTreeFile=TaxTree.DefaultTreeFile;}
			}
			
			else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			} 
			
			else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		size=size_;
		k=k_;
		rcomp=rcomp_;
		mode=mode_;
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;
			in2=parser.in2;

			out1=parser.out1;
			
			extin=parser.extin;
		}
		files=(out1==null ? 0 : files_);
		
		assert(mode!=ONE_SKETCH || files<2) : "Multiple output files are not allowed in single-sketch mode.";
		
		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		
		//Adjust interleaved detection based on the number of input files
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		//Ensure there is an input file
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		ffout=makeFFArray(out1, files, overwrite, append);
		
//		//Ensure output files can be written
//		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
//			outstream.println((out1==null)+", "+out1);
//			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
//		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2, taxTreeFile, giTableFile)){
			throw new RuntimeException("\nCan't read to some input files.\n");
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, taxTreeFile, giTableFile)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
		tool=new SketchTool(size, k, 1);
		
		if(taxTreeFile!=null){
			taxtree=loadTaxTree();
		}else{
			taxtree=null;
		}
		
		if(giTableFile!=null){
			loadGiToNcbi();
		}
	}
	
	private static FileFormat[] makeFFArray(String fname0, int files, boolean overwrite, boolean append){
		if(files<1 || fname0==null){return null;}
		String[] fnames=new String[files];
		FileFormat[] ff=new FileFormat[files];
		for(int i=0; i<files; i++){
			String fname=fname0;
			if(files>1){
				assert(fname.indexOf('#')>-1) : "Output name requires # symbol for multiple files.";
				fname=fname.replaceFirst("#", ""+i);
			}
			fnames[i]=fname;
			ff[i]=FileFormat.testOutput(fname, FileFormat.TEXT, null, true, overwrite, append, false);
		}
		
		if(!Tools.testOutputFiles(overwrite, append, false, fnames)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+Arrays.toString(fnames)+"\n");
		}
		
		return ff;
	}
	
	private static TextStreamWriter[] makeTSWArray(FileFormat[] ff){
		if(ff==null || ff.length==0){return null;}
		TextStreamWriter[] tsw=new TextStreamWriter[ff.length];
		for(int i=0; i<ff.length; i++){
			tsw[i]=new TextStreamWriter(ff[i]);
			tsw[i].start();
		}
		return tsw;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the reads in separate threads
		spawnThreads(cris);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStream(cris);
		
		//TODO: Write sketch
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		{
			t.stop();
			
			//Calculate units per nanosecond
			double rpnano=readsProcessed/(double)(t.elapsed);
			double bpnano=basesProcessed/(double)(t.elapsed);
			
			//Add "k" and "m" for large numbers
			String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
			String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");
			
			//Format the strings so they have they are right-justified
			while(rpstring.length()<8){rpstring=" "+rpstring;}
			while(bpstring.length()<8){bpstring=" "+bpstring;}
			
			outstream.println("Time:                         \t"+t);
			outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
			outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		}
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		
		if(mode==PER_TAXA){
			map=new ConcurrentHashMap<Integer, LongHeapSet>();
		}else if(mode==PER_SEQUENCE){
			tsw=makeTSWArray(ffout);
		}
		
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, i));
		}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		
		//Wait for completion of all threads
		boolean success=true;
		LongHeapSet heap=null;
		for(ProcessThread pt : alpt){
			
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
			
			//Accumulate per-thread statistics
			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			kmersProcessed+=pt.kmersProcessedT;
//			System.err.println("pt.readsProcessedT="+pt.readsProcessedT);
			if(mode==ONE_SKETCH){
				LongHeapSet temp=pt.heap;
				
				if(temp==null){
					//do nothing
				}else if(heap==null){heap=pt.heap;}
				else{heap.add(pt.heap);}
				
				if(heap!=null){
					if(outTaxid>=0){heap.id=outTaxid;}
					if(outName!=null){heap.name=outName;}
				}
			}
			success&=pt.success;
		}
		
		if(heap!=null && heap.name==null){
			heap.name=ReadWrite.stripToCore(in1);
		}
		
		if(ffout!=null){
			if(mode==PER_TAXA){
				tsw=makeTSWArray(ffout);
				success&=writeMap(map);
			}else if(mode==ONE_SKETCH){
				tool.write(new Sketch(heap), ffout[0]);
			}
		}
		
		if(tsw!=null){
			for(int i=0; i<tsw.length; i++){tsw[i].poisonAndWait();}
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
		
		//Do anything necessary after processing
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          I/O Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	
	private boolean writeMap(ConcurrentHashMap<Integer, LongHeapSet> map){
		
		//Determine how many threads may be used
		final int threads=files;

		//Fill a list with WriteThreads
		ArrayList<WriteThread> alwt=new ArrayList<WriteThread>(threads);

		//Start the threads
		for(int i=0; i<threads; i++){
			WriteThread wt=new WriteThread(i);
			alwt.add(wt);
			wt.start();
		}
		
		//Wait for completion of all threads
		boolean success=true;
		for(WriteThread wt : alwt){

			//Wait until this thread has terminated
			while(wt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					wt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			success&=wt.success;
		}
		return success;
	}
	
	private class WriteThread extends Thread{
		
		WriteThread(int tnum_){
			tnum=tnum_;
		}
		
		public void run(){
			success=false;
			for(Entry<Integer, LongHeapSet> entry : map.entrySet()){
				if(entry.getKey().intValue()%files==tnum){
					LongHeapSet heap=entry.getValue();
					Sketch s=new Sketch(heap);
					tool.write(s, tsw[tnum]);
				}
			}
			success=true;
		}
		
		final int tnum;
		boolean success=false;
	}
	
//	private void writeOutput(ConcurrentHashMap<Integer, LongHeapSet> map){
//		TextStreamWriter tsw=new TextStreamWriter(ffout);
//		tsw.start();
//		KeySetView<Integer, LongHeapSet> y=map.keySet();
//		for(Integer x : map.keySet()){
//			LongHeapSet heap=map.get(x);
//			Sketch s=tool.toSketch(heap);
//			tool.write(s, tsw);
//		}
//		tsw.poisonAndWait();
//	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Tax Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	private void loadGiToNcbi(){
		Timer t=new Timer();
		outstream.println("Loading gi to taxa translation table.");
		GiToNcbi.initialize(giTableFile);
		t.stop();
		if(true){
			outstream.println("Time: \t"+t);
			Shared.printMemory();
			outstream.println();
		}
	}
	
	private TaxTree loadTaxTree(){
		Timer t=new Timer();
		outstream.print("\nLoading tax tree; ");
		final TaxTree tree=ReadWrite.read(TaxTree.class, taxTreeFile, true);
		t.stop();
		if(true){
			outstream.println("time: \t"+t);
			Shared.printMemory();
			outstream.println();
		}
		return tree;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This is called if the program runs with no parameters */
	private void printOptions(){
		throw new RuntimeException("TODO"); //TODO
	}
	
	private final long toValue(long kmer, long rkmer){
		long value=(rcomp ? Tools.max(kmer, rkmer) : kmer);
		return value;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final int tid_){
			cris=cris_;
			tid=tid_;
			
			shift=2*k;
			shift2=shift-2;
			mask=~((-1L)<<shift);
			
			heap=new LongHeapSet(size);
		}
		
		//Called by start()
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			processInner();
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Iterate through the reads */
		void processInner(){
			
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			//Check to ensure pairing is as expected
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired()); //Disabled due to non-static access
			}

			//As long as there is a nonempty read list...
			while(reads!=null && reads.size()>0){
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access

				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;

					//Track the initial length for statistics
					final int initialLength1=r1.length();
					final int initialLength2=r1.mateLength();

					//Increment counters
					readsProcessedT+=1+r1.mateCount();
					basesProcessedT+=initialLength1+initialLength2;

					processRead(r1);
					if(r2!=null){processRead(r2);}
					
					{
						int taxID=-1;
						if(taxtree!=null){
							TaxNode node=taxtree.getNode(r1.id);
							if(node!=null){taxID=node.id;}
//							System.err.println("Node: "+r1.id+"\n->\n"+node);
						}
						if(mode==ONE_SKETCH){
							if(heap.id<0){heap.id=taxID;}
							if(heap.name==null){heap.name=r1.id;}
						}else if(mode==PER_SEQUENCE){
							if(heap.size()>0){
								Sketch sketch=new Sketch(heap);
								if(tsw!=null){
									final int choice=(sketch.hashCode()&Integer.MAX_VALUE)%files;
									synchronized(tsw[choice]){
										tool.write(sketch, tsw[choice]);
									}
								}
							}
							heap.clear();
						}else if(mode==PER_TAXA){
							assert(mode==PER_TAXA);
							if(heap.size()>0){
								Integer key=taxID>-1 ? taxID : nextUnknown.getAndIncrement();
								LongHeapSet old=map.get(key);
								if(old==null){
									heap.id=taxID;
									heap.name=r1.id;
									map.put(key, heap);
									heap=new LongHeapSet(size);
								}else{
									synchronized(old){
										old.add(heap);
									}
									heap.clear();
								}
							}
						}
					}
				}

				//Notify the input stream that the list was used
				cris.returnList(ln.id, ln.list.isEmpty());
//				if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access

				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		void processRead(final Read r){
			final byte[] bases=r.bases;
			final byte[] quals=r.quality;
			long kmer=0;
			long rkmer=0;
			int len=0;
			
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				long x=AminoAcid.baseToNumber[b];
				long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=(rkmer>>>2)|(x2<<shift2);
				if(x<0){len=0;}else{len++;}
				if(len>=k){
					kmersProcessedT++;
					long z=toValue(kmer, rkmer);
					long hash=SketchTool.hash(z);
					if(hash>0){heap.add(hash);}
				}
			}
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		/** Number of kmers processed by this thread */
		protected long kmersProcessedT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Thread ID */
		final int tid;

		LongHeapSet heap;

		final int shift;
		final int shift2;
		final long mask;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Hashing            ----------------*/
	/*--------------------------------------------------------------*/
	
//	private static long[][] makeCodes(int symbols, int modes){
//		Random randy=new Random(12345);
//		long[][] r=new long[symbols][modes];
//		for(int i=0; i<symbols; i++){
//			for(int j=0; j<modes; j++){
//				r[i][j]=randy.nextLong();
//			}
//		}
//		return r;
//	}
//	
//	private static final long hash(long kmer){
//		long code=kmer;
//		for(int i=0; i<8; i++){
//			int x=(int)(kmer&0xFF);
//			code^=codes[i][x];
//		}
//		return code;
//	}
//	
//	private static final long[][] codes=makeCodes(8, 256);
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path */
	private String in2=null;

	/** Primary output file path */
	private String out1=null;
	
	/** Override input file extension */
	private String extin=null;
	
	private String giTableFile=null;
	private String taxTreeFile=null;
	
	private String outName=null;
	private int outTaxid=-1;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;
	/** Number of bases processed */
	protected long kmersProcessed=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	

	private ConcurrentHashMap<Integer, LongHeapSet> map;
	private TextStreamWriter tsw[];
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	/** Secondary input file */
	private final FileFormat ffin2;
	
	/** Primary output files */
	private final FileFormat ffout[];
	/** Number of output files */
	private final int files;
	
	private final boolean rcomp;
	private final int mode;
	
	private final SketchTool tool;
	private final int size;
	private final int k;
	
	private final TaxTree taxtree;
	
	private static final int ONE_SKETCH=1, PER_SEQUENCE=2, PER_TAXA=3;
	
	private final AtomicInteger nextUnknown=new AtomicInteger(2000000000);
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=false;
	/** Append to existing output files */
	private boolean append=false;
	
}
