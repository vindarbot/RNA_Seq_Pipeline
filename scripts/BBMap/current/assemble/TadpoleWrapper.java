package assemble;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import align2.Tools;

import jgi.AssemblyStats2;

/**
 * Assembles with multiple kmer lengths to find the best kmer length. 
 * @author Brian Bushnell
 * @date Oct 15, 2015
 *
 */
public class TadpoleWrapper {
	
	public static void main(String[] args){
		HashSet<Integer> set=new HashSet<Integer>();
		ArrayList<String> argList=new ArrayList<String>();
		String contigsName="contigs%.fa";
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if("null".equalsIgnoreCase(b)){b=null;}
			while(a.charAt(0)=='-' && (a.indexOf('.')<0 || i>1 || !new File(a).exists())){a=a.substring(1);}
			
			if(a.equals("k")){
				for(String s2 : b.split(",")){
					set.add(Integer.parseInt(s2));
				}
			}else if(a.equals("out")){
				contigsName=b;
				assert(b.contains("%")) : "Output name must contain % symbol.";
			}else{
				argList.add(arg);
			}
		}
		
		if(set.isEmpty()){
			kmers=new int[] {31};
		}else{
			kmers=new int[set.size()];
			int i=0;
			for(Integer x : set){
				kmers[i]=x;
				i++;
			}
			Arrays.sort(kmers);
		}

		long[] L50=new long[kmers.length];
		long[] contiglen=new long[kmers.length];
		long[] contigs=new long[kmers.length];
		long[] maxContig=new long[kmers.length];
		
		argList.add("");
		argList.add("");
		StringBuilder sb=new StringBuilder("in=");
		
		for(int i=0; i<kmers.length; i++){
			int k=kmers[i];
			argList.set(argList.size()-2, "k="+k);
			argList.set(argList.size()-1, "out="+contigsName.replace("%", ""+k));
			String[] args2=argList.toArray(new String[0]);
			System.gc();
			Tadpole.main(args2);
			
			L50[i]=AssemblyStats2.lastL50;
			contiglen[i]=AssemblyStats2.lastSize;
			contigs[i]=AssemblyStats2.lastContigs;
			maxContig[i]=AssemblyStats2.lastMaxContig;
		}
		
		int index=Tools.maxIndex(L50);
		
		System.err.println("Recommended K:\t"+kmers[index]);
		
	}
	
	private static int[] kmers;
	
}
