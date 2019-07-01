/***************************************************************

   The Subread software package is free software package: 
   you can redistribute it and/or modify it under the terms
   of the GNU General Public License as published by the 
   Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   Subread is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty
   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
   
   See the GNU General Public License for more details.

   Authors: Drs Yang Liao and Wei Shi

  ***************************************************************/
  
  
#include <ctype.h>
#include <errno.h>
#include <unistd.h>
#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include "subread.h"
#include "HelperFunctions.h"
#include "removeDupReads.h"
#include "input-files.h"
extern unsigned int BASE_BLOCK_LENGTH;

int is_read_selected(char *read_selection_list , unsigned int read_number)
{
	unsigned int byte_offset = read_number/8;
	int bit_offset = read_number & 7;
	char old_8values = read_selection_list[byte_offset];
	return ((old_8values >> bit_offset)&1) ? 1:0;

}

void unselect_read_in_list(char *read_selection_list , unsigned int read_number)
{
	unsigned int byte_offset = read_number/8;
	int bit_offset = read_number & 7;
	char old_8values = read_selection_list[byte_offset];
	old_8values &= ~(1 << bit_offset);
	read_selection_list[byte_offset] = old_8values;
}


int parse_base_blocks(char * temp_prefix, chromosome_t * chromosomes, char * read_selection_list, int threshold, int all_threads, int thread_no)
{
	int chromosome_block_no = 0;
	int chromosome_no;

	unsigned short * read_voting_table;

	read_voting_table = (unsigned short *)malloc(sizeof(short) * BASE_BLOCK_LENGTH);
	if(!read_voting_table)
	{
		fatal_memory_size();
		return 1;
	}

	for(chromosome_no=0 ; chromosomes[chromosome_no].chromosome_name[0]; chromosomes++ )
	{
		int block_in_chro = chromosomes[chromosome_no].known_length / BASE_BLOCK_LENGTH+1;
		int i;


		for(i=0;i<block_in_chro; i++)
		{
			char temp_file[300];
			unsigned int block_start_base = BASE_BLOCK_LENGTH * i;
			int rrtv;
			sprintf(temp_file , "%s%s-%04u.bin",temp_prefix, chromosomes[chromosome_no].chromosome_name, i);
			chromosome_block_no++;

			FILE * block_fp = f_subr_open(temp_file,"rb");
			if(!block_fp)
				continue;

			memset(read_voting_table, 0 , sizeof(short) * BASE_BLOCK_LENGTH);

			while(!feof(block_fp))
			{
				base_block_temp_read_t read_record;
				rrtv = fread(&read_record, sizeof(base_block_temp_read_t), 1, block_fp);
				if(rrtv <1) return -1;
				read_voting_table[read_record.pos-block_start_base] ++;
			}

			fseek(block_fp,0,SEEK_SET);

			while(!feof(block_fp))
			{
				base_block_temp_read_t read_record;
				rrtv = fread(&read_record, sizeof(base_block_temp_read_t), 1, block_fp);
				if(rrtv < 1) return -1;
				if(read_voting_table[read_record.pos-block_start_base]>= threshold)
					unselect_read_in_list(read_selection_list , read_record.read_number);
			}


			fclose(block_fp);
			unlink(temp_file);
		}
	}
	free(read_voting_table);
	return 0;
}


int parse_base_blocks_maybe_thread(char * temp_location, chromosome_t * chromosomes, char * read_selection_list, int threshold, int threads)
{
	return parse_base_blocks(temp_location,  chromosomes,   read_selection_list,threshold, threads, 0);
}


int report_remainder(char *in_SAM_file, char *out_SAM_file, char* read_selection_list)
{
	int read_no = 0;
	FILE * fp = f_subr_open(in_SAM_file,"r");
	FILE * out_fp = f_subr_open(out_SAM_file,"w");
	if(!out_fp)
	{
		SUBREADprintf("Unable to open the output file, '%s'.\n",out_SAM_file);
		return 1;
	}

	while(!feof(fp))
	{
		char line_buffer [3000];
		int linelen = read_line(2999, fp, line_buffer, 0);

		if(line_buffer[0]=='@')
		{
			fwrite(line_buffer, linelen, 1, out_fp);
			fputc('\n', out_fp);
		}
		else if((line_buffer[0] >='A' && line_buffer[0]<='Z') || (line_buffer[0] >='a' && line_buffer[0] <='z') || (line_buffer[0] >='0' && line_buffer[0] <='9') || line_buffer[0] =='_' || line_buffer[0] =='.')
		{
			/*
			if(is_read_selected(read_selection_list, read_no))
				fputs("SEL ", out_fp);
			else
				fputs("    ", out_fp);
			*/
			if(is_read_selected(read_selection_list, read_no)){
				fwrite(line_buffer, linelen, 1, out_fp);
				fputc('\n', out_fp);
			}
			read_no ++;
		}
	}
	fclose(fp);
	fclose(out_fp);
	return 0;
}


// Giving temp_location 'NULL' makes the function to use the current directory to store temporary files.
// Giving known_read_count '0' makas the function use its default value: known_read_count = 400,000,000, namely ~50MB memory is allocated to store the selection table (other parts of the program may use more memory)
// This function returns 0 if everything is OK and -1 if memory is overused.
int repeated_read_removal(char * in_SAM_file, int threshold, char * out_SAM_file, char * temp_location, unsigned int known_read_count, int threads)
{
	char temp_file_prefix[300];
	char * read_selection_list;
	unsigned int read_selection_list_size;
	unsigned int real_read_count = 0;
	chromosome_t * known_chromosomes;
	double start_time = miltime();
	unsigned short rand48_seed[3];
	
	// Checking parameters
	if(known_read_count<1) known_read_count = 400000000;
	read_selection_list_size = known_read_count / 8 +1;	// One bit for a read, one more byte for remainder.
	
	SUBREADprintf("Repeated Read Removal\nInput=%s\nOutput=%s\nTemporary Path=%s\nThreshold=%d\nMaximum Reads=%u\n\n", in_SAM_file, out_SAM_file, temp_location, threshold, known_read_count);


	// All parameters are validated.
	// Step 0: the read selection list and the chromosome list are initialized to be empty.
	read_selection_list = (char *) SUBREAD_malloc(sizeof(char) * read_selection_list_size);
	if(!read_selection_list)
	{
		fatal_memory_size();
		return -1;
	}
	memset(read_selection_list, 0xff, sizeof(char) * read_selection_list_size);

	known_chromosomes = (chromosome_t *) SUBREAD_malloc(sizeof(chromosome_t) * XOFFSET_TABLE_SIZE);
	if(!known_chromosomes)
	{
		fatal_memory_size();
		return -1;
	}
	// A chromosome_t is used if the first char of its name is not 0. 
	known_chromosomes[0].chromosome_name[0]=0;
	
	//  Step 1: scanning the SAM file, classify each read into a block. If a chro_name is unseen before, a new block is created for it and it is put into the chromosome list. Each block has a temporary file "temp-delrep-00000-XXXXXX-chrX-21" where chrX is the chromosome name, '21' is the block number, 00000 is the process id and XXXXXX is a randomly string.

	memcpy(rand48_seed, &start_time, 6);
	seed48(rand48_seed);

	char mac_rand[13];
	mac_or_rand_str(mac_rand);
	sprintf(temp_file_prefix, "%s/temp-delrep-%06u-%s-", temp_location==NULL?".":temp_location, getpid(), mac_rand);

	if(break_SAM_file(in_SAM_file, 0, temp_file_prefix, &real_read_count, NULL, known_chromosomes, 0 /* This 0 means that the sequence/quality/cigar fields are not needed in the temp files*/, 0, NULL, NULL, NULL, NULL, NULL)) return -1;

	// Step 2: initialize the read voting table, scanning each temporary file. Each read in the temporary file is a vote to a location in the voting table. Then, each read in the temporary file is scanned again against the voting table. If the mapping location of this read receives >=threshold votes in the table, the read is removed from the read selection list.

	if(parse_base_blocks_maybe_thread(temp_file_prefix, known_chromosomes, read_selection_list, threshold, threads)) return -1;


	// Step 3: rescanning the SAM file, copy a line to output file if its line number is marked 'selected' in the read selection list.

	if(report_remainder(in_SAM_file, out_SAM_file, read_selection_list)) return -1;
	free(read_selection_list);
	free(known_chromosomes);
	SUBREADprintf("Finished.\n");
	return 0;
}

void print_usage_rrr(char * myname)
{
	SUBREADprintf("\nremoveDup Version %s\n\n", SUBREAD_VERSION);
	SUBREADputs("  Remove duplicated reads.");
	SUBREADputs("");
	SUBREADputs("Usage:");
	SUBREADputs("  ./removeDup [options] -i <input_file> -o <output_file>\n");
	SUBREADputs("Required arguments:");
	SUBREADputs("");
	SUBREADputs("  -i <string> Name of input file in SAM format.");
	SUBREADputs("");
	SUBREADputs("  -o <string> Name of output SAM file including filtered reads.");
	SUBREADputs("");
	SUBREADputs("Aptional arguments:");
	SUBREADputs("");
	SUBREADputs("  -r <int>    Specify the duplication cutoff. All the reads mapped to a location");
	SUBREADputs("              are removed from the output if the number of reads mapped to this");
	SUBREADputs("              location is equal or higher than the cutoff. 10 by default.");
	SUBREADputs("");
	SUBREADputs("  -t <string> A directory storing temporary files generated by the program.");
	SUBREADputs("");
	SUBREADputs("  -c <int>    The maximum number of reads the input file can have. 40 million by");
	SUBREADputs("              default.");
	SUBREADputs("");
}

struct option rem_long_options[]={{0,0,0,0}};

#ifdef MAKE_STANDALONE
int main(int argc,char ** argv)
#else
// This function is actually useless because R uses repeated_read_removal() as the entry. However no main() function should be remained while building .so for R. 
int main_repeated_test(int argc,char ** argv)
#endif
{
	char c;
	char input_SAM_file[300];
	char output_SAM_file[300];
	char temp_path[300];
	int threshold, optindex = 0;
	int threads;
	unsigned int read_count;

	input_SAM_file [0] = 0;
	output_SAM_file[0] = 0;
	temp_path[0]=0;

	threshold = 10;
	read_count = 0;
	threads = 0;

	optind=0;
	opterr=1;
	optopt=63;

	if(argc<2)
	{
		print_usage_rrr(argv[0]);
		return 0;
	}
	while ((c = getopt_long (argc, argv, "i:o:r:t:c:?", rem_long_options, &optindex)) != -1)
	{
		switch (c)
		{
			case 'i':
				strncpy(input_SAM_file, optarg,299);
				break;

			case 'o':
				strncpy(output_SAM_file, optarg,299);
				break;


			case 't':
				strncpy(temp_path,  optarg,299);
				break;


			case 'T':
				threads = atoi(optarg);
				break;

			case 'r':
				threshold = atoi(optarg);
				break;

			case 'c':
				read_count = atoi(optarg);
				break;

			case '?':
			default:
				print_usage_rrr(argv[0]);
		}
	}

	
	return repeated_read_removal(input_SAM_file, threshold, output_SAM_file, temp_path[0]?temp_path:NULL, read_count, threads);
}
