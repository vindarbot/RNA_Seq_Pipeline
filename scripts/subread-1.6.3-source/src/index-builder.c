/***************************************************************

   The Subread and Rsubread software packages are free
   software packages:
 
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
  
  
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <getopt.h>
#include <signal.h>
#include <unistd.h>
#include <sys/types.h>
#include "hashtable.h"
#include "gene-value-index.h"
#include "HelperFunctions.h"
#include "gene-algorithms.h"
#include "sorted-hashtable.h"
#include "input-files.h"

#define NO_GENE_DEBUG_
#define _GENE_DEBUG_SIZE_ 40000
#define MIN_READ_SPLICING 2000000


#define MAX_KEY_MATCH GENE_VOTE_SPACE 
int GENE_SLIDING_STEP = 3;
int IS_COLOR_SPACE = 0;
int VALUE_ARRAY_INDEX = 1;
int MARK_NONINFORMATIVE_SUBREADS = 0;
int IS_FORCED_ONE_BLOCK = 0;

#define NEXT_READ 1
#define NEXT_FILE 2
#define FULL_PARTITION 4


void print_build_log(double finished_rate, double read_per_second, double expected_seconds, unsigned long long int total_reads)
{
	print_in_box( 81,0,0,"%4d%%%%, %3d mins elapsed, rate=%.1fk bps/s, total=%llum\r", (int)(finished_rate*100), (int)(miltime()-begin_ftime)/60, read_per_second/1000 ,total_reads/1000000);
}

void copy_non_informative_subread(gehash_t * index_table, gehash_t * noninf_table)
{
	int i,j;
	for (i=0; i< noninf_table -> buckets_number; i++)
	{
		struct gehash_bucket * current_bucket = &(noninf_table->buckets[i]);

		if(current_bucket -> current_items>=1)
		{
			for (j=0; j<current_bucket -> current_items; j++)
			{
				unsigned int noninf_subread = current_bucket -> item_keys[j];
				// the non-informative subreads in the index are marked as they are at 0xffffffff.
				gehash_insert(index_table, noninf_subread, 0xffffffffu);
			}
		}
	}

}


#define MAX_BASES_IN_INDEX 4294900000.0

int build_gene_index(const char index_prefix [], char ** chro_files, int chro_file_number, unsigned int memory_megabytes, int threshold, gehash_t * huge_table, unsigned int * chro_lens)
{
	int file_number, table_no;
	int status = NEXT_FILE;
	unsigned int offset, read_no;
	unsigned int segment_size = (unsigned int)(memory_megabytes * 1024.0 / 8.) * 1024 ;
	long long int all_bases = guess_gene_bases(chro_files,chro_file_number);
	double local_begin_ftime = 0.;

	int chro_table_maxsize=100, dump_res = 0;
	unsigned int * read_offsets = malloc(sizeof(unsigned int) * chro_table_maxsize);
	char * read_names = malloc(MAX_READ_NAME_LEN * chro_table_maxsize);
	gehash_t table;
//	gehash_t huge_table;
	gene_value_index_t value_array_index;

	gene_input_t ginp;

//	SUBREADprintf ("Index items per partition = %u\n\n", segment_size);

	if (chro_file_number > 199)
	{
		SUBREADprintf("ERROR: There are too many FASTA files. You may merge them into one FASTA file.\n");
		return -1;
	}

	if (strlen (index_prefix) > 290)
	{
		SUBREADprintf("ERROR: The path is too long. It should not be longer than 290 chars.\n");
		return -1;
	}

	if(all_bases<1)
	{
		SUBREADprintf("ERROR: File '%s' is inaccessible.\n", chro_files[-all_bases-1]);
		return -1;
	}

	int padding_around_contigs = MAX_READ_LENGTH;
	if(gehash_create_ex(& table, segment_size, 0, SUBINDEX_VER2, GENE_SLIDING_STEP, padding_around_contigs)) return 1;

	if(MARK_NONINFORMATIVE_SUBREADS)
		copy_non_informative_subread(&table, huge_table);

	if(VALUE_ARRAY_INDEX)
		if(gvindex_init(&value_array_index, 0)) return 1;

	file_number = 0;
	offset = table.padding;
	table_no = 0;
	read_no = 0;

	char * fn = malloc(3100);
	bzero(read_offsets, chro_table_maxsize*sizeof(int));
	sprintf(fn, "%s.files", index_prefix);
	unlink(fn);

	status = NEXT_FILE;

	print_in_box(80,0,0,"Build the index...");

	{
		char window [16], last_color_base=-1, last_last_color_base=-1;
		int i, read_len = 0;
		unsigned int int_key = 0, array_int_key = 0;
		int skips=0, all_skips = 0;

		//Pre-fill
		while(1)
		{
			//Subread Cycle
			char next_char;

			if (status == NEXT_FILE)
			{
				if(file_number == chro_file_number)
				{
					FILE * fp;

					geinput_close(&ginp);

					//SUBREADprintf ("Processing chromosome files ...\n");

					sprintf (fn, "%s.%02d.%c.tab", index_prefix, table_no, IS_COLOR_SPACE?'c':'b');
					SUBREADfflush(stdout);

					if(!dump_res)dump_res |= gehash_dump(&table, fn);

					if(VALUE_ARRAY_INDEX)
					{
						sprintf (fn, "%s.%02d.%c.array", index_prefix, table_no, IS_COLOR_SPACE?'c':'b');
						if(!dump_res)dump_res |= gvindex_dump(&value_array_index, fn);
						gvindex_destory(&value_array_index) ;
					}

					gehash_destory(&table);

					read_offsets[read_no-1] = offset + table.padding;

					for(i=table_no+1; i<100; i++)
					{
						sprintf(fn, "%s.%02d.%c.tab", index_prefix, i, IS_COLOR_SPACE?'c':'b');
						unlink(fn);
						sprintf(fn, "%s.%02d.%c.array", index_prefix, i, IS_COLOR_SPACE?'c':'b');
						unlink(fn);
					}

					sprintf (fn, "%s.reads", index_prefix);
					fp = f_subr_open(fn, "w");
					for (i=0; i<read_no; i++)
						fprintf(fp, "%u\t%s\n", read_offsets[i], read_names+i*MAX_READ_NAME_LEN);

					fclose(fp);

					break;
				}
				else
				{
					if (file_number)
						geinput_close(&ginp);
					geinput_open(chro_files[file_number++], &ginp);
					status = NEXT_READ;
				}
			}
			if (status == NEXT_READ)
			{

				geinput_readline(&ginp, fn, 0);

				if(offset==table.padding)
					local_begin_ftime = miltime();

				if(read_no>0){
					read_offsets[read_no-1] = offset + table.padding;
					offset += 2*table.padding;
				}

				//printf("TTTXT FN=%s\n",fn);

				for(i=0;(fn[i+1] != ' ' && fn[i+1] != '\0' && fn[i+1] != '\t' && i<MAX_CHROMOSOME_NAME_LEN - 1); i++)
					*(read_names + MAX_READ_NAME_LEN*read_no + i) = fn[i+1];

				*(read_names + MAX_READ_NAME_LEN*read_no + i) = 0;

				sprintf(fn, "%s.files", index_prefix);
				FILE * fname_fp = f_subr_open(fn, "a");
				fprintf(fname_fp, "%s\t%s\t%ld\n", read_names+read_no*MAX_READ_NAME_LEN, ginp.filename, ftell(ginp.input_fp));
				fclose(fname_fp);
				
				for (i=0; i<16; i++)
				{
					char nch = geinput_next_char(&ginp);
					if (nch == 'N') skips = 16;
					else if (skips>0) skips--;
					window[i] = nch;
				}

				read_len = 16;
				read_no ++;

				if(read_no >= chro_table_maxsize)
				{
					read_offsets = realloc(read_offsets, 2* chro_table_maxsize * sizeof(unsigned int));
					read_names = realloc(read_names, 2* chro_table_maxsize * MAX_READ_NAME_LEN);
					chro_table_maxsize *= 2;
				}

				if(IS_COLOR_SPACE)
				{
					int_key = genekey2color('A',window);
					if(VALUE_ARRAY_INDEX)
						array_int_key = genekey2int(window,GENE_SPACE_BASE);
					last_last_color_base = -1;
					last_color_base = window[15];
				}
				else
					array_int_key = int_key = genekey2int(window,GENE_SPACE_BASE);
	
			}
			if (status == FULL_PARTITION) 
			{
				int seek_back_reads ;

				if (read_len < MIN_READ_SPLICING)
				{
					seek_back_reads = read_len;
					offset -= seek_back_reads;
					offset += 16;
				}
				else
				{
					seek_back_reads = MIN_READ_SPLICING - 10;
					seek_back_reads -= seek_back_reads%3;
					seek_back_reads += 1;
					offset -= seek_back_reads;
					offset += 16;
				}

				sprintf(fn, "%s.%02d.%c.tab", index_prefix, table_no, IS_COLOR_SPACE?'c':'b');
				SUBREADfflush(stdout);

				if(!dump_res)dump_res |= gehash_dump(&table, fn);
				if(VALUE_ARRAY_INDEX)
				{
					sprintf(fn, "%s.%02d.%c.array", index_prefix, table_no, IS_COLOR_SPACE?'c':'b');
					if(!dump_res)dump_res |= gvindex_dump(&value_array_index, fn);
				}

				table_no ++;

				gehash_destory(&table);
				if(VALUE_ARRAY_INDEX)
					gvindex_destory(&value_array_index);

				read_len -= seek_back_reads;
				read_len += 16;

				while(seek_back_reads)
				{
					fseek(ginp.input_fp, -1, SEEK_CUR);
					char bnch = fgetc(ginp.input_fp);
					if ((bnch >='A' && bnch <= 'Z' ) || (bnch >='a' && bnch <= 'z' ) || bnch == '-' || bnch == 'N' || bnch=='.')seek_back_reads--;
					fseek(ginp.input_fp, -1, SEEK_CUR);
				}

				for (i=0; i<16; i++)
				{
					char nch = geinput_next_char(&ginp);
					if (nch == 'N' ) skips = 16;
					else if (skips>0) skips--;
					window[i] = nch;
				}
			
				if(IS_COLOR_SPACE)
				{
					int_key = genekey2color('A',window);
					if(VALUE_ARRAY_INDEX)
						array_int_key = genekey2int(window,GENE_SPACE_BASE);
					last_last_color_base = -1;
					last_color_base = window[15];
				}
				else
					array_int_key = int_key = genekey2int(window, GENE_SPACE_BASE);
				
				if(gehash_create_ex(&table, segment_size, 0, SUBINDEX_VER2, GENE_SLIDING_STEP, padding_around_contigs)) return 1;
				if(MARK_NONINFORMATIVE_SUBREADS)
					copy_non_informative_subread(&table, huge_table);
				if(VALUE_ARRAY_INDEX)
					if(gvindex_init(&value_array_index, offset)) return 1;
			}
	
			status = 0;

			if(skips || (IS_COLOR_SPACE && last_last_color_base<0))
				all_skips ++;
			else
			{
				int is_no_info = gehash_exist(huge_table, int_key);
				if(!is_no_info)
				{
					//SUBREADprintf("INSERT KEY=%u AT %u\n", int_key, offset);
					if(gehash_insert(&table, int_key, offset - (IS_COLOR_SPACE?1:0))) return 1;
				}
				if(VALUE_ARRAY_INDEX)
				{
					gvindex_set(&value_array_index, offset - (IS_COLOR_SPACE?0:0), array_int_key, padding_around_contigs);
				}
			}

			if((!IS_FORCED_ONE_BLOCK) && table.current_items >= segment_size && (read_len > MIN_READ_SPLICING || read_len < 32))
			{
				status = FULL_PARTITION;
				continue;
			}

			for (i=0; i<GENE_SLIDING_STEP; i++)
			{
				next_char = geinput_next_char(&ginp);
				if(next_char < 0)
				{
					gvindex_set(&value_array_index, offset - (IS_COLOR_SPACE?0:0), array_int_key, padding_around_contigs);

					if (next_char == -1) status = NEXT_READ;
					if (next_char == -2) status = NEXT_FILE;
					if (next_char == -3) return 0;
					break;
				}
				//SUBREADprintf("NEXT_CH=%c\n", next_char);

				if (next_char == 'N' )skips = 16;
				else if (skips>0){
					skips--;
					last_color_base = -1;
				}

				int_key = int_key << 2;


				if (IS_COLOR_SPACE)
				{
					if(last_color_base>0)
						int_key += chars2color(last_color_base, next_char);
					if(VALUE_ARRAY_INDEX)
					{
						array_int_key = array_int_key << 2;
						array_int_key += base2int (next_char);
					}

					last_last_color_base = last_color_base;
					last_color_base = next_char;
				}
				else
				{
					int_key += base2int (next_char); 
					array_int_key = int_key;
				}

				if(all_bases > 12){
					if (offset>0 && offset % (all_bases / 12) == 0)
					{
						double finished_rate = offset*1.0/all_bases;
						double base_per_second = offset / (miltime()-local_begin_ftime);
						double ETA = (all_bases-offset) / base_per_second;
						print_build_log(finished_rate,base_per_second,ETA, all_bases);
						SUBREADfflush(stdout) ;
					}
				}

				offset ++;
				read_len ++;

				if(offset > 0xFFFFFFFDU)	
				{
					SUBREADprintf("ERROR: The chromosome data contains too many bases. The size of the input FASTA file(s) should be less than 4G Bytes\n") ;
					return -1;
				}

			}

		}
	}
	free(read_names);
	free(read_offsets);
	if(dump_res){
		SUBREADprintf("No index was built.\n");
		sprintf(fn, "%s.files", index_prefix);
		unlink(fn);
		sprintf(fn, "%s.reads", index_prefix);
		unlink(fn);
		int index_i;
		for(index_i = 0; index_i <= 99; index_i++){
			sprintf(fn, "%s.%02d.b.tab", index_prefix, index_i);
			unlink(fn);
			sprintf(fn, "%s.%02d.c.tab", index_prefix, index_i);
			unlink(fn);
			sprintf(fn, "%s.%02d.b.array", index_prefix, index_i);
			unlink(fn);
			sprintf(fn, "%s.%02d.c.array", index_prefix, index_i);
			unlink(fn);
		}
	}
	free(fn);
	return 0;
}

int add_repeated_subread(gehash_t * tab , unsigned int subr, unsigned char ** huge_index)
{
	unsigned int times;

	int huge_byte = (subr>>2) &0x3fffffff;
	int huge_offset = (subr % 4) * 2;
	unsigned int byte_value = huge_index[ (huge_byte >> 20) & 1023 ][huge_byte&0xfffff] ;

	int huge_value = (byte_value>> huge_offset) & 0x3;
	if(huge_value <3)
	{
		huge_value ++;
		huge_index[ (huge_byte >> 20) & 1023 ][huge_byte&0xfffff] = (byte_value & (~(0x3 << huge_offset))) | (huge_value << huge_offset);
		return 0;
	}

	int matched = gehash_get(tab, subr, &times, 1);
	if(matched)
	{
		gehash_update(tab, subr, times+1);
	}
	else
		if(gehash_insert(tab, subr,4)) return 1;
	return 0;
}


int scan_gene_index(const char index_prefix [], char ** chro_files, int chro_file_number, int threshold, gehash_t *huge_table)
{
	int file_number, i ,j;
	int status = NEXT_FILE;
	unsigned int offset, read_no;
	double local_begin_ftime = miltime();
	long long int all_bases = guess_gene_bases(chro_files,chro_file_number);

	gehash_t occurance_table;
	unsigned char * huge_index[1024];

	for(i=0;i<1024;i++)
	{
		huge_index[i] = (unsigned char *)malloc(1024*1024); 
		if(!huge_index[i])
		{
			for(j=0;j<i;j++) free(huge_index[j]);
			SUBREADprintf("ERROR: You need at least one point five gigabytes of memory for building the index.\n");
			return -1;
		}
		memset(huge_index[i], 0 , 1024*1024);
	}


	if(gehash_create(&occurance_table , 100000000, 0)) return 1;


	gene_input_t ginp;

	print_in_box(80,0,0,"Scan uninformative subreads in reference sequences ...");

	if (chro_file_number > 199)
	{
		SUBREADprintf("ERROR: There are too many FASTA files. You may merge them into one FASTA file.\n");
		return -1;
	}

	if (strlen (index_prefix) > 290)
	{
		SUBREADprintf("ERROR: The path is too long. It should not be longer than 290 chars.\n");
		return -1;
	}
	if(all_bases<1)
	{
		SUBREADprintf("ERROR: File '%s' is inaccessible.\n", chro_files[-all_bases-1]);
		return -1;
	}


	file_number = 0;
	offset = 0;
	read_no = 0;

	char * fn = malloc(3100);
	sprintf(fn, "%s.files", index_prefix);
	unlink(fn);

	status = NEXT_FILE;


	{
		char window [16], last_color_base=-1, last_last_color_base=-1;
		int i, read_len = 0;
		unsigned int int_key = 0, array_int_key = 0;
		int skips=0, all_skips = 0;

		//Pre-fill
		while(1)
		{
			//Subread Cycle
			char next_char;

			if (status == NEXT_FILE)
			{
				if(file_number == chro_file_number)
				{
					geinput_close(&ginp);

					break;
				}
				else
				{
					if (file_number)
						geinput_close(&ginp);
					geinput_open(chro_files[file_number++], &ginp);
					status = NEXT_READ;
				}
			}
			if (status == NEXT_READ)
			{

				geinput_readline(&ginp, fn, 0);

				for (i=0; i<16; i++)
				{
					char nch = geinput_next_char(&ginp);
					if (nch == 'N') skips = 16;
					else if (skips>0) skips--;
					window[i] = nch;
				}
				read_len = 16;
				read_no ++;

				if(IS_COLOR_SPACE)
				{
					int_key = genekey2color('A',window);
					if(VALUE_ARRAY_INDEX)
						array_int_key = genekey2int(window,GENE_SPACE_BASE);
					last_last_color_base = -1;
					last_color_base = window[15];
				}
				else
					array_int_key = int_key = genekey2int(window,GENE_SPACE_BASE);
	
			}
	
			status = 0;

			if(skips || (IS_COLOR_SPACE && last_last_color_base<0))
				all_skips ++;
			else
			{
				add_repeated_subread(&occurance_table, int_key, huge_index);
			}


			for (i=0; i<GENE_SLIDING_STEP; i++)
			{
				next_char = geinput_next_char(&ginp);
				if(next_char < 0)
				{
					if (next_char == -1) status = NEXT_READ;
					if (next_char == -2) status = NEXT_FILE;
					if (next_char == -3) return 0;
					break;
				}

				if (next_char == 'N' )skips = 16;
				else if (skips>0){
					skips--;
					last_color_base = -1;
				}

				int_key = int_key << 2;


				if (IS_COLOR_SPACE)
				{
					if(last_color_base>0)
						int_key += chars2color(last_color_base, next_char);
					if(VALUE_ARRAY_INDEX)
					{
						array_int_key = array_int_key << 2;
						array_int_key += base2int (next_char);
					}

					last_last_color_base = last_color_base;
					last_color_base = next_char;
				}
				else
				{
					int_key += base2int (next_char); 
					array_int_key = int_key;
				}

				offset ++;
				read_len ++;


				if (all_bases>12){
					if (offset % (all_bases / 12) == 0)
					{
						double finished_rate = offset*1.0/all_bases;
						double base_per_second = offset / (miltime()-local_begin_ftime);
						double ETA = (all_bases-offset) / base_per_second;
						print_build_log(finished_rate,base_per_second,ETA, all_bases);
						SUBREADfflush(stdout) ;
					}
				}

				if(offset > 0xFFFFFFFDU)	
				{
					SUBREADprintf("ERROR: The chromosome data contains too many bases. The size of the input FASTA files should be less than 4G Bytes\n") ;
					return -1;
				}

			}

		}
	}



	free(fn);

	for (i=0; i<occurance_table.buckets_number; i++)
	{
		struct gehash_bucket * current_bucket = &(occurance_table.buckets[i]);

		if(current_bucket -> current_items>=1)
		{
			for (j=0; j<current_bucket -> current_items; j++)
			{
				if(current_bucket -> item_values [j] > threshold)
				{
					if(gehash_insert(huge_table, current_bucket -> item_keys[j], 1)) return 1;
				}
			}
		}
	}

	for(i=0;i<1024;i++)
		free(huge_index[i]);
	gehash_destory(&occurance_table);


	if(huge_table -> current_items)
	{
		print_in_box(80,0,0,"%llu uninformative subreads were found.", huge_table -> current_items);
		print_in_box(80,0,0,"These subreads were excluded from index building.");
	}

	return 0;
}


char *rtrim(char *s)
{
	char* back = s + strlen(s);
	while(isspace(*--back));
	*(back+1) = '\0';
	return s;
}

int ERROR_FOUND_IN_FASTA = 0;
#define CHAR_ESC 27
void check_and_convert_warn(char * FN, long long int fpos_line_head, unsigned line_no, int line_pos, char * msg, FILE * log_fp)
{
	int x1,brs=0;
	long long int back_search_ptr;
	char * line_buf = malloc(MAX_READ_LENGTH+1);

	ERROR_FOUND_IN_FASTA += 1;

	fprintf(log_fp,"\n");

	//fprintf(log_fp,"%c[33m", CHAR_ESC);
	for(x1=0;x1<81;x1++)
		fprintf(log_fp,"=");
	fprintf(log_fp,"\n");
	//fprintf(log_fp,"%c[32m", CHAR_ESC);
	fprintf(log_fp,"Input file '%s':\n", FN);
	//fprintf(log_fp,"%c[31m", CHAR_ESC);
	fprintf(log_fp,"%s\n", msg);
	//fprintf(log_fp,"%c[33m", CHAR_ESC);
	for(x1=0;x1<81;x1++)
		fprintf(log_fp,".");
	fprintf(log_fp,"\n");
	//fprintf(log_fp,"%c[37m", CHAR_ESC);
	

	FILE * warn_fp = f_subr_open(FN, "r");

	for(back_search_ptr = fpos_line_head - 1; back_search_ptr>=0; back_search_ptr--)
	{
		int nch;
		fseeko(warn_fp, back_search_ptr, SEEK_SET);
		nch = fgetc(warn_fp);
		if(nch == '\n') brs++;
		if(brs >2) break;
		fseeko(warn_fp, back_search_ptr, SEEK_SET);
	}

	if(back_search_ptr<=0)brs++;

	int print_line_no = line_no - brs + 1;
	while(1)
	{
		char * ret = fgets(line_buf, MAX_READ_LENGTH, warn_fp);
		if(!ret)break;

		/*
		if(ftello(warn_fp) > fpos_line_head)
			fprintf(log_fp,"%c[9m%c[31m", CHAR_ESC, CHAR_ESC);
		else
			fprintf(log_fp,"%c[29m%c[37m", CHAR_ESC, CHAR_ESC);
		*/
		fprintf(log_fp," % 9d ", print_line_no++);

		rtrim(line_buf);
		
		fprintf(log_fp,"%s%s\n",line_buf,strlen(line_buf)<16?"              ":"");
		if(ftello(warn_fp) > fpos_line_head)
			break;
	}
	for(x1=0;x1<line_pos+11;x1++)
		fprintf(log_fp," ");
	fprintf(log_fp,"^\n");

	//fprintf(log_fp,"%c[29m%c[37m", CHAR_ESC, CHAR_ESC);
	for(x1=0;x1<2;x1++)
	{
		char * ret = fgets(line_buf, MAX_READ_LENGTH, warn_fp);
		if(!ret)break;
		fprintf(log_fp," % 9d ", print_line_no++);
		fprintf(log_fp,"%s",line_buf);
	}
	fclose(warn_fp);
	//fprintf(log_fp,"%c[33m", CHAR_ESC);
	for(x1=0;x1<81;x1++)
		fprintf(log_fp,"=");
	fprintf(log_fp,"\n");
	fprintf(log_fp,"\n");

	//fprintf(log_fp,"%c[0m", CHAR_ESC);
	free(line_buf);
}

int check_and_convert_FastA(char ** input_fas, int fa_number, char * out_fa, unsigned int ** chrom_lens, FILE * log_fp, char * log_fn)
{
	int is_R_warnned = 0, is_repeated_chro= 0;
	char * line_buf = malloc(MAX_READ_LENGTH);
	char * read_head_buf = malloc(MAX_READ_LENGTH * 3);
	unsigned int inp_file_no, line_no;
	int written_chrs = 0, is_disk_full = 0;
	int chrom_lens_max_len = 100;
	int chrom_lens_len = 0;
	ERROR_FOUND_IN_FASTA = 0;
	FILE * out_fp = f_subr_open(out_fa,"w");

	if(!out_fp)
	{
		SUBREADprintf("ERROR: the output directory is not writable, but the index builder needs to create temporary files in the current directory. Please change the working directory and rerun the index builder.\n");
		return -1;
	}

	(*chrom_lens) = malloc(chrom_lens_max_len*sizeof(unsigned int));
	memset((*chrom_lens), 0, chrom_lens_max_len*sizeof(unsigned int));


	HashTable * rep_name_table = HashTableCreate(9999);
	HashTableSetDeallocationFunctions(rep_name_table, free, NULL);
	HashTableSetKeyComparisonFunction(rep_name_table, (int (*) (const void *, const void *)) strcmp);
	HashTableSetHashFunction(rep_name_table , fc_chro_hash);
	
	print_in_box( 80,0,0,"Check the integrity of provided reference sequences ...");
	for(inp_file_no = 0; inp_file_no < fa_number; inp_file_no++)
	{
		FILE * in_fp = f_subr_open(input_fas[inp_file_no],"r");
		long long int last_read_head_pos = 0;
		unsigned int last_read_line_no = 1;

		if(!in_fp)
		{
			SUBREADprintf("ERROR: Input file '%s' is not found or is not accessible. No index was built.\n", input_fas[inp_file_no]);
			HashTableDestroy(rep_name_table);
			return -1;
		}

		line_no = 0;
		int is_head_written=0;
		read_head_buf[0]=0;
		while(!feof(in_fp))
		{
			long long int line_head_pos = ftello(in_fp);
			unsigned int read_len = 0;
			char * ret = fgets(line_buf, MAX_READ_LENGTH-1 , in_fp);
			if(!ret) break;
			line_no ++;
			int line_buf_len = strlen(line_buf);

			for(; line_buf[line_buf_len-1] == '\r' || line_buf[line_buf_len-1]=='\n' ;line_buf_len--)
			{
				if(line_buf[line_buf_len-1]=='\r')
				{
					if(!is_R_warnned)
					{
						is_R_warnned=1;
						check_and_convert_warn(input_fas[inp_file_no], line_head_pos, line_no, line_buf_len -1 ,"This line ends with '\\r\\n'. It is not a problem for building the index but we suggest to use Unix-styled line breaks.", log_fp);
					}	
				}
				line_buf[line_buf_len-1] =0;
			}


			if(line_buf_len<1)
			{
				check_and_convert_warn(input_fas[inp_file_no], line_head_pos, line_no ,0 ,"This line is empty. This is not allowed in the FASTA file.", log_fp);
				continue;
			}

			if(line_buf[0]=='>')
			{
				if(line_no>1 &&!is_head_written)
				{
					check_and_convert_warn(input_fas[inp_file_no], last_read_head_pos, last_read_line_no, 0,"This sequence has less than 16 bases. It is ignored in the index because no subreads can be extracted.", log_fp);
				}
				is_head_written = 0;
				last_read_line_no = line_no;
				last_read_head_pos = line_head_pos;
				read_len = 0;
				read_head_buf[0]=0;

				strcat(read_head_buf, line_buf);
				strcat(read_head_buf, "\n");

				int chro_name_end=1;
				while(line_buf[chro_name_end])
				{
					if(line_buf[chro_name_end]==' '||line_buf[chro_name_end]=='|' || line_buf[chro_name_end]=='\t')break;
					chro_name_end++;
				}

				line_buf[chro_name_end] = 0;
				int is_exist = HashTableGet(rep_name_table , line_buf+1) - NULL;
				if(is_exist)
				{
					SUBREADprintf("ERROR: repeated chromosome name '%s' is observed in the FASTA file(s).\nThe index was NOT built.\n", line_buf+1);
					is_repeated_chro=1;
					break;
				}

				char * keymem = malloc(chro_name_end);
				strcpy(keymem, line_buf+1);
				HashTablePut(rep_name_table, keymem, NULL+1);

			}
			else if(line_head_pos<1)
			{
				check_and_convert_warn(input_fas[inp_file_no], 0, 0, 0 ,"This file is not started with a header line. It seems not to be a FASTA file.", log_fp);
			}
			else
			{
				int xk2;
				for(xk2=0; xk2 < line_buf_len; xk2++)
				{
					int nextch = line_buf[xk2];
					int lowerch = tolower(nextch);
					if(!( lowerch == 'a' || lowerch == 't' || lowerch == 'g' || lowerch == 'c' || nextch == '.' || nextch=='-' || lowerch=='n'))
					{	
						check_and_convert_warn(input_fas[inp_file_no], line_head_pos, line_no, xk2, "The pointed base was converted to an 'A'.", log_fp);
						line_buf[xk2] = 'A';
					}
					else if(nextch == '.' || nextch=='-' ||  lowerch=='n')
						line_buf[xk2] = 'A';
					else
						line_buf[xk2] = toupper(nextch);
				}

				read_len += line_buf_len;

				if(read_len > 16 && !is_head_written)
				{
					fputs(read_head_buf, out_fp);
					written_chrs++;
					is_head_written = 1;

					chrom_lens_len++;
					if((chrom_lens_max_len-1) <= chrom_lens_len)
					{
						(*chrom_lens) = realloc((*chrom_lens), 2*chrom_lens_max_len*sizeof(unsigned int));
						chrom_lens_max_len*=2;
					}
				}

				if(is_head_written)
				{
					int line_buf_len = strlen(line_buf);
					int writen_len = fprintf(out_fp,"%s\n", line_buf);
					if(writen_len < line_buf_len){
						SUBREADprintf("ERROR: unable to write into the temporary file. Please check the free space of the output directory.\n");
						is_disk_full = 1;
						break;
					}
					(*chrom_lens)[chrom_lens_len-1] = read_len;
					(*chrom_lens)[chrom_lens_len] = 0;
				}
				else
				{
					strcat(read_head_buf, line_buf);
					strcat(read_head_buf, "\n");
				}
				
			}
		}

		fclose(in_fp);
		if(is_disk_full) break;
	}


	HashTableDestroy(rep_name_table);
	free(line_buf);
	free(read_head_buf);
	fclose(out_fp);

	if(!written_chrs){
		SUBREADprintf("ERROR: No index was built because there were no subreads extracted. A chromosome needs at least 16 bases to be indexed.");
		return 1;
	}

	if(is_repeated_chro|| is_disk_full)
		return 1;

	if(ERROR_FOUND_IN_FASTA)
	{
		print_in_box( 80,0,0,"There were %d notes for reference sequences.", ERROR_FOUND_IN_FASTA);
		print_in_box( 89,0,0,"The notes can be found in the log file, %c[36m'%s'%c[0m.", CHAR_ESC, log_fn, CHAR_ESC);
	}
	else	print_in_box( 80,0,0,"No format issues were found");

	return 0;
}

char * tmp_file_for_signal;

void SIGINT_hook(int param)
{
	#ifdef MAKE_STANDALONE
	if(tmp_file_for_signal[0])
	{
		unlink(tmp_file_for_signal);
		SUBREADprintf("\n\nReceived a terminal signal. The temporary file was removed. The index was NOT built successfully. Please DO NOT use the new index until they are rebuilt.\n\n");
	}

	exit(param);
	#endif
}

static struct option ib_long_options[]={
	{0, 0, 0, 0}
};

#ifdef MAKE_STANDALONE
int main(int argc,char ** argv)
#else
int main_buildindex(int argc,char ** argv)
#endif
{
	int threshold = 100, optindex=0;
	int memory_limit;	// 8000 MBytes
	char output_file[330], c, tmp_fa_file[300], log_file_name[350];
	char *ptr_tmp_fa_file[1];
	unsigned int * chromosome_lengths;

	if(sizeof(char *)>4) memory_limit=8000;
	else memory_limit=3000;

	ptr_tmp_fa_file[0]=tmp_fa_file;
	output_file[0] = 0;
	tmp_fa_file[0] = 0;
	tmp_file_for_signal = tmp_fa_file;

	IS_FORCED_ONE_BLOCK = 0;
	GENE_SLIDING_STEP = 3;
	IS_COLOR_SPACE = 0;

	SUBREADprintf("\n");

	optind = 0;
	
	while ((c = getopt_long (argc, argv, "kvcBFM:o:f:D?", ib_long_options, &optindex)) != -1)
		switch(c)
		{
			case 'B':
				IS_FORCED_ONE_BLOCK =1;
				break;
			case 'F':
				GENE_SLIDING_STEP =1;
				break;
			case 'v':
				core_version_number("Subread-buildindex");
				return 0;
			case 'c':
				IS_COLOR_SPACE = 1;
				break;
			case 'M':
				memory_limit = atoi(optarg);
				break;
			case 'f':
				threshold = atoi(optarg);
				break;
			case 'o':
				strncpy(output_file, optarg, 299);
				break;
			case 'k':
				if(memcmp(SUBREAD_VERSION , "1.3.",4)==0)
				{
					SUBREADprintf("The \"-k\" option is not supported in version 1.3.x. ");
					return -1;
				}
				MARK_NONINFORMATIVE_SUBREADS = 1;
				break;	
			case '?':
				return -1 ;
		}

	if (argc == optind || !output_file[0])
	{
		SUBREADprintf("Version %s\n\n", SUBREAD_VERSION);

	/*
		SUBREADputs("Usage:");
		SUBREADputs("");
		SUBREADputs(" ./subread-buildindex [options] -o <basename> {FASTA file1} [FASTA file2] ...");
		SUBREADputs("");
		SUBREADputs("Required arguments:");
		SUBREADputs("");
		SUBREADputs("    -o <basename>   base name of the index to be created");
		SUBREADputs("");
		SUBREADputs("Optional arguments:");
		SUBREADputs("");
		SUBREADputs("    -G              build a gapped index for the reference genome. 16bp subreads");
		SUBREADputs("");
		SUBREADputs("    -M <int>        size of requested memory(RAM) in megabytes. Index is split into blocks if necessary.");
		SUBREADputs("");
		SUBREADputs("    -f <int>        specify the threshold for removing uninformative subreads");
		SUBREADputs("                    (highly repetitive 16mers in the reference). 24 by default.");
		SUBREADputs("");
		SUBREADputs("    -c              build a color-space index.");
		SUBREADputs("");
		SUBREADputs("    -v              output version of the program.");
		SUBREADputs("");
		SUBREADputs("For more information about these arguments, please refer to the User Manual.\n");
	*/
		 SUBREADputs("Usage:");
		 SUBREADputs("");
		 SUBREADputs(" ./subread-buildindex [options] -o <basename> {FASTA file1} [FASTA file2] ...");
		 SUBREADputs("");
		 SUBREADputs("Required arguments:");
		 SUBREADputs("");
		 SUBREADputs("    -o <basename>   base name of the index to be created");
		 SUBREADputs("");
		 SUBREADputs("Optional arguments:");
		 SUBREADputs("");
		 SUBREADputs("    -F              build a full index for the reference genome. 16bp subreads");
		 SUBREADputs("                    will be extracted from every position of the reference");
		 SUBREADputs("                    genome. Size of the index is typically 3 times the size of");
		 SUBREADputs("                    index built from using the default setting.");
		 SUBREADputs("");
		 SUBREADputs("    -B              create one block of index. The built index will not be split");
		 SUBREADputs("                    into multiple pieces. This makes the largest amount of");
		 SUBREADputs("                    memory be requested when running alignments, but it enables");
		 SUBREADputs("                    the maximum mapping speed to be achieved. This option");
		 SUBREADputs("                    overrides -M when it is provided as well.");
		 SUBREADputs("");
		 SUBREADputs("    -M <int>        size of requested memory(RAM) in megabytes, 8000 by default.");
		 SUBREADputs("");
		 SUBREADputs("    -f <int>        specify the threshold for removing uninformative subreads");
		 SUBREADputs("                    (highly repetitive 16mers in the reference). 100 by default.");
		 SUBREADputs("");
		 SUBREADputs("    -c              build a color-space index.");
		 SUBREADputs("");
		 SUBREADputs("    -v              output version of the program.");
		 SUBREADputs("");
		 SUBREADputs("For more information about these arguments, please refer to the User Manual.\n");

		return -1 ;
	}


	// **********************************************
	//	print config summary
	// **********************************************

	print_subread_logo();

	SUBREADputs("");
	print_in_box(80, 1, 1, "indexBuilder setting");
	print_in_box(80, 0, 1, "");
	print_in_box(80, 0, 0, "               Index name : %s", get_short_fname(output_file));
	print_in_box(80, 0, 0, "              Index space : %s", IS_COLOR_SPACE?"color-space":"base-space");

	if(IS_FORCED_ONE_BLOCK)
	{
		print_in_box(80, 0, 0, "          One block index : yes");
		memory_limit = GENE_SLIDING_STEP==1?22000:11500;
	}
	else
	{
		if(memory_limit > 12000 && GENE_SLIDING_STEP>2)
		{
			print_in_box(80, 0, 0, "                   Memory : %u -> %u Mbytes", memory_limit, 12000);
			memory_limit = 12000;
		}
		else
			print_in_box(80, 0, 0, "                   Memory : %u Mbytes", memory_limit);
	}
	print_in_box(80, 0, 0, "         Repeat threshold : %d repeats", threshold);
	print_in_box(80, 0, 0, " Distance to next subread : %d", GENE_SLIDING_STEP);
	print_in_box(80, 0, 0, "");
	print_in_box(80, 0, 0, "              Input files : %d file%s in total",  argc - optind, (argc - optind>1)?"s":"");

	int x1;
	for(x1=0;x1< argc - optind; x1++)
	{
		char * fasta_fn = *(argv+optind+x1);
		int f_type = probe_file_type_fast(fasta_fn);
		char o_char = 'o';
		if(f_type != FILE_TYPE_FASTA){
			o_char = '?';
		}
		print_in_box(94, 0, 0, "                            %c[32m%c%c[36m %s%c[0m", CHAR_ESC, o_char, CHAR_ESC,  get_short_fname(fasta_fn) , CHAR_ESC);
	}
	print_in_box(80, 0, 0, "");
	print_in_box(80, 2, 1, "http://subread.sourceforge.net/");
	SUBREADputs("");

	print_in_box(80, 1, 1, "Running");
	print_in_box(80, 0, 0, "");

	for(x1=0;x1< argc - optind; x1++)
	{
		char * fasta_fn = *(argv+optind+x1);
		int f_type = probe_file_type_fast(fasta_fn);
		if(f_type != FILE_TYPE_FASTA && f_type != FILE_TYPE_NONEXIST){
			SUBREADprintf("ERROR: '%s' is not a Fasta file.\n", fasta_fn);
			if(f_type == FILE_TYPE_GZIP_FASTA){
				SUBREADprintf("The index builder does not accept gzipped Fasta files.\nPlease decompress the gzipped Fasza files before building the index.\n");
			}
			STANDALONE_exit(-1);
		}
	}


	begin_ftime = miltime();

	for(x1 = strlen(output_file); x1 >=0; x1--){
		if(output_file[x1]=='/'){
			memcpy(tmp_fa_file, output_file, x1);
			tmp_fa_file[x1]=0;
			break;
		}
	}
	if(tmp_fa_file[0]==0)strcpy(tmp_fa_file, "./");

	sprintf(tmp_fa_file+strlen(tmp_fa_file), "/subread-index-sam-%06u-XXXXXX", getpid());
	int tmpfdd = mkstemp(tmp_fa_file);
	if(tmpfdd == -1){
		SUBREADprintf("ERROR: cannot create temp file\n");
		return -1;
	}

	sprintf(log_file_name, "%s.log", output_file);
	FILE * log_fp = f_subr_open(log_file_name,"w");

	signal (SIGTERM, SIGINT_hook);
	signal (SIGINT, SIGINT_hook);


	int ret = check_and_convert_FastA(argv+optind , argc - optind, tmp_fa_file, &chromosome_lengths, log_fp, log_file_name);
	if(log_fp)
		fclose(log_fp);

	if(!ret)
	{
		gehash_t huge_table;
		gehash_create(& huge_table, 50000000 * (GENE_SLIDING_STEP==1?3:1), 0);
		ret = ret || scan_gene_index(output_file, ptr_tmp_fa_file , 1, threshold, &huge_table);
		ret = ret || build_gene_index(output_file, ptr_tmp_fa_file , 1,  memory_limit, threshold, &huge_table, chromosome_lengths);
		if(!ret){
			print_in_box(80, 0, 1, "Total running time: %.1f minutes.", (miltime()-begin_ftime)/60);
			print_in_box(89, 0, 1, "Index %c[36m%s%c[0m was successfully built!", CHAR_ESC, output_file, CHAR_ESC);
		}
		gehash_destory(& huge_table);
		//     ^^^^^^^ should be destroy
		free(chromosome_lengths);
	}

	unlink(tmp_fa_file);
	tmp_fa_file[0]=0;

	print_in_box(80, 0, 0, "");
	print_in_box(80, 2, 1, "http://subread.sourceforge.net/");
	SUBREADputs("");
	return ret;
}

