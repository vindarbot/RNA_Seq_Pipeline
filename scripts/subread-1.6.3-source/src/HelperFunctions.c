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
#include <string.h>
#include <assert.h>
#include <math.h>


#ifdef MACOS

#include <sys/types.h>
#include <sys/socket.h>
#include <sys/ioctl.h>
#include <sys/sysctl.h>
#include <net/if.h>
#include <net/if_dl.h>
#include <netinet/in.h>
#include <arpa/inet.h>

#else

#ifndef __MINGW32__
#include <sys/ioctl.h>
#include <netinet/in.h>
#include <net/if.h>
#endif
#include <unistd.h>
#endif


#include "subread.h"
#include "input-files.h"
#include "seek-zlib.h"
#include "gene-algorithms.h"
#include "HelperFunctions.h"


char * get_short_fname(char * lname){
	char * ret = lname;

	int x1;
	for(x1 = strlen(lname)-1; x1>=0; x1--){
		if(lname [x1] == '/'){
			ret = lname + x1 + 1;
			break;
		}
	}
	return ret;
}



// This assumes the first part of Cigar has differet strandness to the main part of the cigar.
// Pos is the LAST WANTED BASE location before the first strand jump (split by 'b' or 'n').
// The first base in the read actually has a larger coordinate than Pos. 
// new_cigar has to be at least 100 bytes.
unsigned int reverse_cigar(unsigned int pos, char * cigar, char * new_cigar) {
	int cigar_cursor = 0;
	new_cigar[0]=0;
	unsigned int tmpi=0;
	int last_piece_end = 0;
	int last_sec_start = 0;
	unsigned int chro_pos = pos, this_section_start = pos, ret = pos;
	int is_positive_dir = 0;
	int read_cursor = 0;
	int section_no = 0;

	for(cigar_cursor = 0 ;  ; cigar_cursor++)
	{
		if( cigar [cigar_cursor] == 'n' ||  cigar [cigar_cursor] == 'b' ||  cigar [cigar_cursor] == 0)
		{
			int xk1, jmlen=0, nclen=strlen(new_cigar);
			char jump_mode [13];

			if(cigar [cigar_cursor] !=0)
			{
				sprintf(jump_mode, "%u%c", tmpi,  cigar [cigar_cursor] == 'b'?'n':'b');
				jmlen = strlen(jump_mode);
			}

			for(xk1=nclen-1;xk1>=0; xk1--)
				new_cigar[ xk1 +  last_piece_end + jmlen - last_sec_start ] = new_cigar[ xk1 ];
			new_cigar [nclen + jmlen + last_piece_end - last_sec_start ] = 0;

			memcpy(new_cigar , jump_mode, jmlen);
			memcpy(new_cigar + jmlen , cigar + last_sec_start, last_piece_end - last_sec_start);

			last_sec_start = cigar_cursor+1;

			if(is_positive_dir && cigar [cigar_cursor] !=0)
			{
				if(cigar [cigar_cursor] == 'b') chro_pos -= tmpi - read_cursor - 1;
				else	chro_pos += tmpi - read_cursor - 1;
			}
			if((!is_positive_dir) && cigar [cigar_cursor] !=0)
			{
				if(cigar [cigar_cursor] == 'b') chro_pos = this_section_start - tmpi - read_cursor - 1;
				else	chro_pos = this_section_start + tmpi - read_cursor - 1;
			}

			this_section_start = chro_pos;

			if(section_no == 0)
				ret = chro_pos;

			is_positive_dir = ! is_positive_dir;
			section_no++;
			tmpi=0;
		}
		else if(isalpha(cigar [cigar_cursor]))
		{
			if(cigar [cigar_cursor]=='M' || cigar [cigar_cursor] == 'S')
				read_cursor += tmpi;
			tmpi=0;
			last_piece_end = cigar_cursor+1;
		}
		else tmpi = tmpi*10 + (cigar [cigar_cursor] - '0');

		if(cigar [cigar_cursor] == 0)break;
	}

	SUBREADprintf("REV CIGAR: %s  =>  %s\n", cigar, new_cigar);
	return ret;
}

unsigned int find_left_end_cigar(unsigned int right_pos, char * cigar){
	int delta_from_right = 0;
	int cigar_cursor = 0;
	unsigned int tmpi = 0;
	while(1){
		int nch = cigar[cigar_cursor++];
		if(nch == 0) break;
		if(isdigit(nch)){
			tmpi = tmpi * 10 + nch - '0';
		}else{
			if(nch == 'M'||nch == 'D' || nch == 'N'){
				delta_from_right +=tmpi;
			}
			tmpi = 0;
		}
	}
	return right_pos - delta_from_right;
}


char contig_fasta_int2base(int v){
	if(v == 1) return 'A';
	if(v == 2) return 'T';
	if(v == 3) return 'G';
	if(v == 4) return 'C';
	return 'N';
}

int contig_fasta_base2int(char base){
	base = tolower(base);
	if((base) == 'a'){ return 1;}
	else if((base) == 't' || (base) == 'u'){ return 2;}
	else if((base) == 'g'){ return 3;}
	else if((base) == 'c'){ return 4;}
	else return 15 ;
}

int get_contig_fasta(fasta_contigs_t * tab, char * chro, unsigned int pos, int len, char * out_bases){
	unsigned int this_size = HashTableGet( tab -> size_table, chro ) - NULL;
	if(this_size > 0){
		if(this_size >= len && pos <= this_size - len){
			char * bin_block = HashTableGet(tab -> contig_table, chro );
			unsigned int bin_byte = pos / 2;
			int bin_bit = 4*(pos % 2), x1;

			for(x1 = 0 ;x1 < len; x1++)
			{
				int bin_int = (bin_block[bin_byte] >> bin_bit) & 0xf;
				if(bin_bit == 4) bin_byte++;
				bin_bit = (bin_bit == 4)?0:4;
				out_bases[x1] = contig_fasta_int2base(bin_int);
			}

			return 0;
		}
	} 
	return 1;
}

void destroy_contig_fasta(fasta_contigs_t * tab){
	HashTableDestroy( tab -> size_table );
	HashTableDestroy( tab -> contig_table );
}
int read_contig_fasta(fasta_contigs_t * tab, char * fname){
	FILE * fp = f_subr_open(fname, "r");
	if(fp != NULL){
		tab -> contig_table = HashTableCreate(3943);
		tab -> size_table = HashTableCreate(3943);

		HashTableSetDeallocationFunctions(tab -> contig_table, free, free);
		HashTableSetDeallocationFunctions(tab -> size_table, NULL, NULL);

		HashTableSetKeyComparisonFunction(tab -> contig_table, fc_strcmp_chro);
		HashTableSetKeyComparisonFunction(tab -> size_table, fc_strcmp_chro);

		HashTableSetHashFunction(tab -> contig_table, fc_chro_hash);
		HashTableSetHashFunction(tab -> size_table, fc_chro_hash);

		char chro_name[MAX_CHROMOSOME_NAME_LEN];
		unsigned int inner_cursor = 0, current_bin_space = 0;
		int status = 0;
		char * bin_block = NULL;
		chro_name[0]=0;

		while(1){
			char nch = fgetc(fp);
			if(status == 0){
				assert(nch == '>');
				status = 1;
			}else if(status == 1){
				if(inner_cursor == 0){
					bin_block = calloc(sizeof(char),10000);
					current_bin_space = 10000;
				}
				if(nch == '|' || nch == ' ') status = 2;
				else if(nch == '\n'){
					status = 3;
					inner_cursor = 0;
				}else{
					chro_name[inner_cursor++] = nch;
					chro_name[inner_cursor] = 0;
				}
			}else if(status == 2){
				if(nch == '\n'){
					status = 3;
					inner_cursor = 0;
				}
			}else if(status == 3){
				if(nch == '>' || nch <= 0){
					char * mem_chro = malloc(strlen(chro_name)+1);
					strcpy(mem_chro, chro_name);
					HashTablePut(tab -> size_table , mem_chro, NULL + inner_cursor);
					HashTablePut(tab -> contig_table , mem_chro, bin_block);
		//			SUBREADprintf("Read '%s' : %u bases\n", chro_name, inner_cursor);
					inner_cursor = 0;
					status = 1;
					if(nch <= 0) break;
				}else if(nch != '\n'){
					int bin_bytes = inner_cursor / 2;
					int bin_bits = 4*(inner_cursor % 2);
					int base_int = contig_fasta_base2int(nch);
					if(bin_bytes >= current_bin_space){
						unsigned int new_bin_space = current_bin_space / 4 * 5;
						if(current_bin_space > 0xffff0000 /5 * 4){
							assert(0);
						}
						bin_block = realloc(bin_block, new_bin_space);
						memset(bin_block + current_bin_space, 0, new_bin_space - current_bin_space);
						current_bin_space = new_bin_space;
					}
					bin_block[bin_bytes] |= (base_int << bin_bits);
					inner_cursor++;
				}
			}
		}

		fclose(fp);
		return 0;
	}
	return 1;
}

int RSubread_parse_CIGAR_Extra_string(int FLAG, char * MainChro, unsigned int MainPos, const char * CIGAR_Str, const char * Extra_Tags, int max_M, char ** Chros, unsigned int * Staring_Chro_Points, unsigned short * Section_Start_Read_Pos, unsigned short * Section_Length, int * is_junction_read){
	int ret = RSubread_parse_CIGAR_string(MainChro, MainPos, CIGAR_Str, max_M, Chros, Staring_Chro_Points, Section_Start_Read_Pos, Section_Length, is_junction_read);

	char read_main_strand = (((FLAG & 0x40)==0x40) == ((FLAG & 0x10) == 0x10 ))?'-':'+';
	int tag_cursor=0;
	//SUBREADprintf("EXTRA=%s\n", Extra_Tags);
	int status = PARSE_STATUS_TAGNAME;
	char tag_name[2], typechar=0;
	int tag_inner_cursor=0;

	char current_fusion_char[MAX_CHROMOSOME_NAME_LEN];
	unsigned int current_fusion_pos = 0;
	char current_fusion_strand = 0;
	char current_fusion_cigar[max_M * 15];
	current_fusion_cigar [0] =0;
	current_fusion_char [0]=0;

	while(1){
		int nch = Extra_Tags[tag_cursor];
		if(status == PARSE_STATUS_TAGNAME){
			tag_name[tag_inner_cursor++] = nch;
			if(tag_inner_cursor == 2){
				status = PARSE_STATUS_TAGTYPE;
				tag_cursor += 1;
				assert(Extra_Tags[tag_cursor] == ':');
			}
		}else if(status == PARSE_STATUS_TAGTYPE){
			typechar = nch;
			tag_cursor +=1;
			assert(Extra_Tags[tag_cursor] == ':');
			tag_inner_cursor = 0;
			status = PARSE_STATUS_TAGVALUE;
		}else if(status == PARSE_STATUS_TAGVALUE){
			if(nch == '\t' || nch == 0 || nch == '\n'){
				if(current_fusion_cigar[0] && current_fusion_char[0] && current_fusion_pos && current_fusion_strand){
					//SUBREADprintf("ENTER CALC:%s\n", current_fusion_char );
					unsigned int left_pos = current_fusion_pos;
					if(current_fusion_strand!=read_main_strand)
						left_pos = find_left_end_cigar(current_fusion_pos, current_fusion_cigar);
					ret += RSubread_parse_CIGAR_string(current_fusion_char, left_pos, current_fusion_cigar, max_M - ret, Chros + ret, Staring_Chro_Points+ ret, Section_Start_Read_Pos+ ret, Section_Length + ret, is_junction_read);

					current_fusion_pos = 0;
					current_fusion_strand = 0;
					current_fusion_cigar [0] =0;
					current_fusion_char [0]=0;
					//SUBREADprintf("EXIT CALC:%s\n", current_fusion_char );
				}

				tag_inner_cursor = 0;
				status = PARSE_STATUS_TAGNAME;
			}else{
				if(tag_name[0]=='C' && tag_name[1]=='C' && typechar == 'Z'){
					current_fusion_char[tag_inner_cursor++]=nch;
					current_fusion_char[tag_inner_cursor]=0;
				}else if(tag_name[0]=='C' && tag_name[1]=='G' && typechar == 'Z'){
					current_fusion_cigar[tag_inner_cursor++]=nch;
					current_fusion_cigar[tag_inner_cursor]=0;
				}else if(tag_name[0]=='C' && tag_name[1]=='P' && typechar == 'i'){
					current_fusion_pos = current_fusion_pos * 10 + (nch - '0');
				}else if(tag_name[0]=='C' && tag_name[1]=='T' && typechar == 'Z'){
					//SUBREADprintf("pos=%d %c -> %c\n", tag_cursor, current_fusion_strand, nch);
					current_fusion_strand = nch;
					//SUBREADprintf("spo=%d %c -> %c\n", tag_cursor, current_fusion_strand, nch);
				}
			}
		}

		if(nch == 0 || nch == '\n'){
			assert(status == PARSE_STATUS_TAGNAME);
			break;
		}

		tag_cursor++;
		//SUBREADprintf("CUR=%d [%s], c=%d\n", tag_cursor, Extra_Tags, Extra_Tags[tag_cursor]);
	}
	return ret;
}

int RSubread_parse_CIGAR_string(char * chro , unsigned int first_pos, const char * CIGAR_Str, int max_M, char ** Section_Chromosomes, unsigned int * Section_Start_Chro_Pos,unsigned short * Section_Start_Read_Pos, unsigned short * Section_Chro_Length, int * is_junction_read)
{
	unsigned int tmp_int=0;
	int cigar_cursor=0;
	unsigned short current_section_chro_len=0, current_section_start_read_pos = 0, read_cursor = 0;
	unsigned int chromosome_cursor=first_pos;
	int ret=0, is_first_S = 1;

	for(cigar_cursor=0; ; cigar_cursor++)
	{
		char ch = CIGAR_Str[cigar_cursor];

		if(ch >='0' && ch <= '9')
		{
			tmp_int=tmp_int*10+(ch - '0');
		}
		else
		{
			if(ch == 'S'){
				if(is_first_S) current_section_start_read_pos = tmp_int;
				read_cursor += tmp_int;
			}
			else if(ch == 'M' || ch == 'X' || ch == '=') {
				read_cursor += tmp_int;
				current_section_chro_len += tmp_int;
				chromosome_cursor += tmp_int;
			} else if(ch == 'N' || ch == 'D' || ch=='I' || ch == 0) {
				if('N' == ch)(*is_junction_read)=1;
				if(ret < max_M)
				{
					if(current_section_chro_len>0)
					{
						Section_Chromosomes[ret] = chro;
						Section_Start_Chro_Pos[ret] = chromosome_cursor - current_section_chro_len;
						Section_Start_Read_Pos[ret] = current_section_start_read_pos;
						Section_Chro_Length[ret] = current_section_chro_len;
						ret ++;
					}
				}else assert(0);
				current_section_chro_len = 0;
				if(ch == 'I') read_cursor += tmp_int;
				else if(ch == 'N' || ch == 'D') chromosome_cursor += tmp_int;
				current_section_start_read_pos = read_cursor;

				if(ch == 0) break;
			}
			//printf("C=%c, TV=%d, CC=%d, RC=%d\n", ch, tmp_int, chromosome_cursor, current_section_chro_len);
			tmp_int = 0;
			is_first_S = 0;
		}
		if(cigar_cursor>100) return -1;
	}

	return ret;
}

void display_sections(char * CIGAR_Str)
{
	//int is_junc=0;
	int Section_Start_Chro_Pos[FC_CIGAR_PARSER_ITEMS];
	unsigned short Section_Start_Read_Pos[FC_CIGAR_PARSER_ITEMS];
	unsigned short Section_Chro_Length[FC_CIGAR_PARSER_ITEMS];

	int retv = 0;//RSubread_parse_CIGAR_string(CIGAR_Str, Section_Start_Chro_Pos, Section_Start_Read_Pos, Section_Chro_Length, &is_junc);

	int x1;
	SUBREADprintf("Cigar=%s ; Sections=%d\n", CIGAR_Str, retv);
	for(x1=0; x1<retv; x1++)
	{
		SUBREADprintf("   Section #%d: chro_offset=%d, read_offset=%u  length=%u\n",x1, Section_Start_Chro_Pos[x1], Section_Start_Read_Pos[x1], Section_Chro_Length[x1]);
	}
	SUBREADprintf("\n");
	
}


#define GECV_STATE_BEFORE 10
#define GECV_STATE_NAME 20
#define GECV_STATE_GAP 30
#define GECV_STATE_VALUE 40
#define GECV_STATE_QVALUE 50
#define GECV_STATE_QV_END 60
#define GECV_STATE_ERROR 9999

int GTF_extra_column_istoken_chr(char c)
{
	return (isalpha(c)||isdigit(c)||c=='_');
}

int GTF_extra_column_value(const char * Extra_Col, const char * Target_Name, char * Target_Value, int TargVal_Size)
{
	int state = GECV_STATE_BEFORE;
	int col_cursor = 0, is_correct_name=0;
	char name_buffer[200];
	int name_cursor = 0, value_cursor=-1;

	while(1)
	{
		if(name_cursor>190) return -1;
		char nch = Extra_Col[col_cursor];
		if(nch == '\n' || nch == '\r') nch = 0;
		if(state == GECV_STATE_BEFORE)
		{
			if(GTF_extra_column_istoken_chr(nch))
			{
				name_buffer[0] = nch;
				name_cursor = 1;
				state = GECV_STATE_NAME;
			}
			else if(nch != ' ' && nch != 0)
			{
				state = GECV_STATE_ERROR;
			}
		}
		else if(state == GECV_STATE_NAME)
		{
			if(nch == ' ' || nch == '=')
			{
				state = GECV_STATE_GAP;
				name_buffer[name_cursor] = 0;
				is_correct_name = (strcmp(name_buffer , Target_Name) == 0);
				//printf("CORR=%d : '%s'\n", is_correct_name, name_buffer);
			}
			else if(nch == '"')
			{
				name_buffer[name_cursor] = 0;
				is_correct_name = (strcmp(name_buffer , Target_Name) == 0);
				state = GECV_STATE_QVALUE;
				if(is_correct_name)
					value_cursor = 0;
			}
			else if(GTF_extra_column_istoken_chr(nch))
				name_buffer[name_cursor++] = nch;
			else
			{
				state = GECV_STATE_ERROR;
				//printf("ERR2  : '%c'\n", nch);
			}
			
		}
		else if(state == GECV_STATE_GAP)
		{
			if(nch == '"')
			{
				state = GECV_STATE_QVALUE;
				if(is_correct_name)
					value_cursor = 0;
			}
			else if(nch != '=' && isgraph(nch))
			{
				state = GECV_STATE_VALUE;
				if(is_correct_name)
				{
					Target_Value[0]=nch;
					value_cursor = 1;
				}
			}
			else if(nch != ' ' && nch != '=')
				state = GECV_STATE_ERROR;
		}
		else if(state == GECV_STATE_VALUE)
		{
			if(nch == ';' || nch == 0)
			{
				state = GECV_STATE_BEFORE;
				if(is_correct_name)
				{
					Target_Value[value_cursor] = 0;
				}
				is_correct_name = 0;
			}
			else{
				if(value_cursor < TargVal_Size-1 && is_correct_name)
					Target_Value[value_cursor++] = nch;
			}
		}
		else if(state == GECV_STATE_QVALUE)
		{
			if(nch == '"')
			{
				state = GECV_STATE_QV_END;
				if(is_correct_name)
					Target_Value[value_cursor] = 0;
				is_correct_name = 0;
			}
			else
			{
				if(value_cursor < TargVal_Size-1 && is_correct_name)
				{
					if(nch !=' ' || value_cursor>0)
						Target_Value[value_cursor++] = nch;
				}
			}
		}
		else if(state == GECV_STATE_QV_END)
		{
			if(nch == ';' || nch == 0)
				state = GECV_STATE_BEFORE;
			else if(nch != ' ')
				state = GECV_STATE_ERROR;
				
		}

		if (GECV_STATE_ERROR == state){
			Target_Value[0]=0;
			return -1;
		}
		if (nch == 0)
		{
			if(state == GECV_STATE_BEFORE && value_cursor>0)
			{
				int x1;
				for(x1 = value_cursor-1; x1>=0; x1--)
				{
					if(Target_Value[x1] == ' '){
						value_cursor --;
						Target_Value[x1]=0;
					}
					else break;
				}

				if(value_cursor>0)
					return value_cursor;
			}
			Target_Value[0]=0;
			return -1;
		}
		col_cursor++;
	}
}


void hpl_test2_func()
{
	char * extra_column = " gene_id \"PC4-013  \"; 013=ABCD  ; PC4 =  CCXX  ";
	char * col_name = "gene_id";
	char col_val[100];

	int col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);

	col_name = "013";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);

	col_name = "PC4";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);


	col_name = "XXX";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);

	extra_column = "gene_id =   \"PC4-013  ;=\"  ;013 = AXXD ; PC4=x";
	col_name = "013";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);

	col_name = "gene_id";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);


	col_name = "PC4";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);




	extra_column = " gene_id\"  PC4-013  ;=  \"; XXX='123' ;013 :ABCD  ; PC4 =  CCXX=  ;";
	col_name = "013";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);


	col_name = "XXX";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);


	col_name = "PC4";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);

	col_name = "gene_id";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);



	extra_column = "gene_id \"653635\"; transcript_id \"TR:653635\";";
	col_name = "gene_id";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);





}

void hpl_test1_func()
{
	display_sections("");
	display_sections("*");
	display_sections("5S10M2D10M800N12M3I12M450N12M12D99M6S");
	display_sections("110M2I10M800N32M3I12M6S");
	display_sections("200S110M2I10M800N32M3I12M200N40M");
	display_sections("3M1663N61M1045N36M3D20M66N10M2D10M77N3M1663N61M1045N36M3D20M66N103M1663N61M1045N36M3D20M66N9M");
}

#ifdef RSUBREAD_TEST_HELPER_FUNCTIONS
void main()
#else
void testi_helper_1_main()
#endif
{
	hpl_test1_func();
}

char *str_replace(char *orig, char *rep, char *with) {
    char *result; // the return string
    char *ins;    // the next insert point
    char *tmp;    // varies
    int len_rep;  // length of rep
    int len_with; // length of with
    int len_front; // distance between rep and end of last rep
    int count;    // number of replacements

    if (!orig)
        return NULL;
    if (!rep)
        rep = "";
    len_rep = strlen(rep);
    if (!with)
        with = "";
    len_with = strlen(with);

    ins = orig;
    for (count = 0; NULL != (tmp = strstr(ins, rep)); ++count) {
        ins = tmp + len_rep;
    }
    tmp = result = malloc(strlen(orig) + (len_with - len_rep) * count + 1);

    if (!result)
        return NULL;

    while (count--) {
        ins = strstr(orig, rep);
        len_front = ins - orig;
        tmp = strncpy(tmp, orig, len_front) + len_front;
        tmp = strcpy(tmp, with) + len_with;
        orig += len_front + len_rep; // move to next "end of rep"
    }
    strcpy(tmp, orig);
    return result;
}



// rule: the string is ABC123XXXXXX...
// This is the priroity:
// First, compare the letters part.
// Second, compare the pure numeric part.
// Third, compare the remainder.
int strcmp_number(char * s1, char * s2)
{
	int x1 = 0;
	int ret = 0;

	while(1)
	{
		char nch1 = s1[x1];
		char nch2 = s2[x1];

		if((!nch1) || !nch2){return nch2?1:(nch1?(-1):0);}
		if(isdigit(nch1) && isdigit(nch2))break;

		ret = nch1 - nch2;
		if(ret) return ret;
		x1++;
	}

	int v1 = 0, v2 = 0;
	while(1)
	{
		char nch1 = s1[x1];
		char nch2 = s2[x1];
		if((!nch1) || !nch2){
			if(nch1 || nch2)
				return nch2?(-1):1;
			break;
		}
		int is_chr1_digit = isdigit(nch1);
		int is_chr2_digit = isdigit(nch2);

		if(is_chr1_digit || is_chr2_digit)
		{
			if(is_chr1_digit && is_chr2_digit)
			{
				v1 = v1*10+(nch1-'0');
				v2 = v2*10+(nch2-'0');
			}
			else
			{
				ret = nch1 - nch2;
				return ret;
			}
		}
		else break;
		x1++;
	}

	if(v1==v2)
		return strcmp(s1+x1, s2+x1);	
	else
	{
		ret = v1 - v2;
		return ret;
	}
}



int mac_str(char * str_buff)
{
#if defined(FREEBSD) || defined(__MINGW32__)
	return 1;
#else
#ifdef MACOS
    int         mib[6], x1, ret = 1;
	size_t		len;
    char            *buf;
    unsigned char       *ptr;
    struct if_msghdr    *ifm;
    struct sockaddr_dl  *sdl;


	for(x1 = 0 ; x1 < 40; x1++)
    {
		mib[0] = CTL_NET;
		mib[1] = AF_ROUTE;
		mib[2] = 0;
		mib[3] = AF_LINK;
		mib[4] = NET_RT_IFLIST;
		mib[5] = x1;

		if (sysctl(mib, 6, NULL, &len, NULL, 0) >=0) {
			buf = malloc(len);
			if (sysctl(mib, 6, buf, &len, NULL, 0) >=0) {

				ifm = (struct if_msghdr *)buf;
				sdl = (struct sockaddr_dl *)(ifm + 1);
				ptr = (unsigned char *)LLADDR(sdl);

				if(sdl -> sdl_nlen < 1) continue;

				char * ifname = malloc(sdl -> sdl_nlen + 1);

				memcpy(ifname, sdl -> sdl_data, sdl -> sdl_nlen);
				ifname[sdl -> sdl_nlen] = 0;
				if(ifname[0]!='e'){
					free(ifname);
					continue;
				}
				free(ifname);

				sprintf(str_buff,"%02X%02X%02X%02X%02X%02X",  *ptr, *(ptr+1), *(ptr+2),
					*(ptr+3), *(ptr+4), *(ptr+5));
				ret = 0;
				break;
			}
			free(buf);
		}
	}
    return ret;
#else
#if defined(IFHWADDRLEN)
    struct ifreq ifr;
    struct ifconf ifc;
    char buf[1024];
    int success = 0;

    int sock = socket(AF_INET, SOCK_DGRAM, IPPROTO_IP);
    if (sock == -1) { /* handle error*/ };

    ifc.ifc_len = sizeof(buf);
    ifc.ifc_buf = buf;
    if (ioctl(sock, SIOCGIFCONF, &ifc) == -1) { /* handle error */ }

    struct ifreq* it = ifc.ifc_req;
    const struct ifreq* const end = it + (ifc.ifc_len / sizeof(struct ifreq));

    for (; it != end; ++it) {
        strcpy(ifr.ifr_name, it->ifr_name);
        if (ioctl(sock, SIOCGIFFLAGS, &ifr) == 0) {
            if (! (ifr.ifr_flags & IFF_LOOPBACK)) { // don't count loopback
                if (ioctl(sock, SIOCGIFHWADDR, &ifr) == 0) {
                      success = 1;
                      break;
                }
            }
        }
    }

    close(sock);

    unsigned char mac_address[6];

    if (success){
	memcpy(mac_address, ifr.ifr_hwaddr.sa_data, 6);
	    int x1;
	    for(x1 = 0; x1 < 6; x1++){
		 sprintf(str_buff+2*x1, "%02X",mac_address[x1]);
	    }
		return 0;
	}
#endif
	return 1;
#endif
#endif
}

int rand_str(char * str_buff){
	int ret = 1;
	FILE * fp = fopen("/dev/urandom","r");
	if(fp){
		int x1;
		for(x1=0; x1<6; x1++){
			sprintf(str_buff + 2*x1 , "%02X", fgetc(fp));
		}
		fclose(fp);
		ret = 0;
	}
	return ret;
}

int mathrand_str(char * str_buff){
	myrand_srand((int)(miltime()*100));
	int x1;
	for(x1 = 0; x1 < 6; x1++){
		sprintf(str_buff+2*x1, "%02X", myrand_rand() & 0xff );
	}
	return 0;
}

int mac_or_rand_str(char * str_buff){
	return mac_str(str_buff) && rand_str(str_buff) && mathrand_str(str_buff);
}

#define PI_LONG 3.1415926535897932384626434L

long double fast_fractorial(unsigned int x, long double * buff, int buflen){
	if(x<2) return 0;
	
	if(buff != NULL && x < buflen && buff[x]!=0){
		return buff[x];
	}
	long double ret;

	if(x<50){
		int x1;
		ret =0.L;
		for(x1 = 2 ; x1 <= x; x1++) ret += logl((long double)(x1));
	}else{
		ret = logl(x)*1.0L*x - 1.0L*x + 0.5L * logl(2.0L*PI_LONG* x) + 1.L/(12.L*x) - 1.L/(360.L* x*x*x) +  1.L/(1260.L* x*x*x*x*x) -  1.L/(1680.L*x*x*x*x*x*x*x);
	}
	if(buff != NULL && x < buflen) buff[x]=ret;
	return ret;
}


#define BUFF_4D 36

long double fast_freq( unsigned int tab[2][2] , long double * buff, int buflen){
	int x0 = -1;
	if(buff && buflen > BUFF_4D * BUFF_4D * BUFF_4D * BUFF_4D && tab[0][0] < BUFF_4D && tab[0][1] < BUFF_4D && tab[1][0] < BUFF_4D && tab[1][1] < BUFF_4D ){
		x0 = buflen + tab[0][0]* BUFF_4D * BUFF_4D * BUFF_4D + tab[0][1] * BUFF_4D * BUFF_4D + tab[1][0] * BUFF_4D + tab[1][1];
		if(buff[x0]!=0) return buff[x0];
		
	}
	long double ret = fast_fractorial(tab[0][0]+tab[0][1],buff,buflen)+fast_fractorial(tab[1][0]+tab[1][1],buff,buflen) + 
		fast_fractorial(tab[0][0]+tab[1][0],buff,buflen)+fast_fractorial(tab[0][1]+tab[1][1],buff,buflen) -
		fast_fractorial(tab[0][0],buff,buflen) - fast_fractorial(tab[0][1],buff,buflen) - 
		fast_fractorial(tab[1][0],buff,buflen) - fast_fractorial(tab[1][1],buff,buflen) - 
		fast_fractorial(tab[0][0] + tab[0][1] + tab[1][0] + tab[1][1],buff,buflen);
	if(x0>=0) buff[x0] = ret;
	return ret;
}



double fast_fisher_test_one_side(unsigned int a, unsigned int b, unsigned int c, unsigned int d, long double * buff, int buflen){
	unsigned int tab[2][2];
	long double P0, P_delta, ret;

	tab[0][0]=a;
	tab[0][1]=b;
	tab[1][0]=c;
	tab[1][1]=d;

	int x_a, y_a, x_min=-1, y_min=-1;
	unsigned int min_a = 0xffffffff;
	for(x_a = 0; x_a < 2; x_a++)
		for(y_a = 0; y_a < 2; y_a++){
			if(tab[x_a][y_a]<min_a){
				min_a = tab[x_a][y_a];
				x_min = x_a;
				y_min = y_a;
			}
		}
	P_delta = fast_freq(tab, buff, buflen);
	P0 = ret = exp(P_delta);
	//printf("P0=%LG\n", P0);
	if(min_a>0){
		unsigned int Qa = min_a;
		unsigned int Qb = tab[x_min][!y_min];
		unsigned int Qc = tab[!x_min][y_min];
		unsigned int Qd = tab[!x_min][!y_min];
		for(; ; ){
			min_a --;
			Qb++;Qc++;
			P_delta -= logl(Qb*Qc);
			P_delta += logl(Qa*Qd);
			Qa--;Qd--;
			ret += expl(P_delta);
		//	printf("%LG %LG %LG\n", ret, 1 - (ret - P0), expl(P_delta));
			if(min_a < 1) break;
		}
	}

	double ret1 = ret, ret2 = 1 - (ret - P0);

	if(min(ret1, ret2) < 0){
		return 0.0;
	}
	return  min(ret1, ret2) ;

}

int load_features_annotation(char * file_name, int file_type, char * gene_id_column, char * transcript_id_column, char * used_feature_type,
 void * context, int do_add_feature(char * gene_name, char * transcript_name, char * chro_name, unsigned int start, unsigned int end, int is_negative_strand, void * context)  ){
	char * file_line = malloc(MAX_LINE_LENGTH+1);
	int lineno = 0, is_GFF_txid_warned = 0, is_GFF_geneid_warned = 0, loaded_features = 0;
	autozip_fp afp;
	int aret = autozip_open(file_name, &afp);

	if(aret < 0){
		SUBREADprintf("Error: unable to open the annotation file : %s\n", file_name);
		return -1;
	}

	while(1){
		int is_tx_id_found = 0, is_gene_id_found = 0, is_negative_strand = -1;
		char * token_temp = NULL, * feature_name, *transcript_id = NULL, * chro_name = NULL;
		char feature_name_tmp[FEATURE_NAME_LENGTH], txid_tmp[FEATURE_NAME_LENGTH];
		feature_name = feature_name_tmp;
		
		unsigned int start = 0, end = 0;
		aret = autozip_gets(&afp, file_line, MAX_LINE_LENGTH);
		if(aret < 1) break;

		lineno++;
		if(is_comment_line(file_line, file_type, lineno-1))continue;

		if(file_type == FILE_TYPE_RSUBREAD)
		{
			feature_name = strtok_r(file_line,"\t",&token_temp);
			int feature_name_len = strlen(feature_name);
			if(feature_name_len > FEATURE_NAME_LENGTH) feature_name[FEATURE_NAME_LENGTH -1 ] = 0;
			
			chro_name = strtok_r(NULL,"\t", &token_temp);
			int chro_name_len = strlen(chro_name);
			if(chro_name_len > MAX_CHROMOSOME_NAME_LEN) chro_name[MAX_CHROMOSOME_NAME_LEN -1 ] = 0;

			char * start_ptr = strtok_r(NULL,"\t", &token_temp);
			char * end_ptr = strtok_r(NULL,"\t", &token_temp);

			if(start_ptr == NULL || end_ptr == NULL){
				SUBREADprintf("\nWarning: the format on the %d-th line is wrong.\n", lineno);
			}
			long long int tv1 = atoll(start_ptr);
			long long int tv2 = atoll(end_ptr);

			if( isdigit(start_ptr[0]) && isdigit(end_ptr[0]) ){
				if(strlen(start_ptr) > 10 || strlen(end_ptr) > 10 || tv1 > 0x7fffffff || tv2> 0x7fffffff){
					SUBREADprintf("\nError: Line %d contains a coordinate greater than 2^31!\n", lineno);
					return -2;
				}
			}else{
				SUBREADprintf("\nError: Line %d contains a format error. The expected annotation format is SAF.\n", lineno);
				return -2;
			}
			
			start = atoi(start_ptr);// start
			end = atoi(end_ptr);//end
			
			char * strand_str = strtok_r(NULL,"\t", &token_temp);
			if(strand_str == NULL)
				is_negative_strand = 0;
			else
				is_negative_strand = ('-' ==strand_str[0]);
				
			is_gene_id_found = 1;
			
		} else if(file_type == FILE_TYPE_GTF) {
			chro_name = strtok_r(file_line,"\t",&token_temp);
			strtok_r(NULL,"\t", &token_temp);// source
			char * feature_type = strtok_r(NULL,"\t", &token_temp);// feature_type
			
			if(strcmp(feature_type, used_feature_type)==0){
				char * start_ptr = strtok_r(NULL,"\t", &token_temp);
				char * end_ptr = strtok_r(NULL,"\t", &token_temp);
				

				if(start_ptr == NULL || end_ptr == NULL){
					SUBREADprintf("\nWarning: the format on the %d-th line is wrong.\n", lineno);
				}
				long long int tv1 = atoll(start_ptr);
				long long int tv2 = atoll(end_ptr);
				
				
				if( isdigit(start_ptr[0]) && isdigit(end_ptr[0]) ){
					if(strlen(start_ptr) > 10 || strlen(end_ptr) > 10 || tv1 > 0x7fffffff || tv2> 0x7fffffff){
						SUBREADprintf("\nError: Line %d contains a coordinate greater than 2^31!\n", lineno);
						return -2;
					}
				}else{
					SUBREADprintf("\nError: Line %d contains a format error. The expected annotation format is GTF/GFF.\n", lineno);
					return -2;
				}
				start = atoi(start_ptr);// start
				end = atoi(end_ptr);//end

				if(start < 1 || end<1 ||  start > 0x7fffffff || end > 0x7fffffff || start > end)
					SUBREADprintf("\nWarning: the feature on the %d-th line has zero coordinate or zero lengths\n\n", lineno);

					
				strtok_r(NULL,"\t", &token_temp);// score
				is_negative_strand = ('-' == (strtok_r(NULL,"\t", &token_temp)[0]));//strand
				strtok_r(NULL,"\t",&token_temp);	// "frame"
				char * extra_attrs = strtok_r(NULL,"\t",&token_temp);   // name_1 "val1"; name_2 "val2"; ...
				if(extra_attrs && (strlen(extra_attrs)>2)){
					int attr_val_len = GTF_extra_column_value(extra_attrs , gene_id_column , feature_name_tmp, FEATURE_NAME_LENGTH);
					if(attr_val_len>0) is_gene_id_found=1;

					if(transcript_id_column){
							transcript_id = txid_tmp;
							attr_val_len = GTF_extra_column_value(extra_attrs , transcript_id_column , txid_tmp, FEATURE_NAME_LENGTH);
							if(attr_val_len>0) is_tx_id_found=1;
							else transcript_id = NULL;
					}
				}

				if(!is_gene_id_found){
					if(!is_GFF_geneid_warned){
						int ext_att_len = strlen(extra_attrs);
						if(extra_attrs[ext_att_len-1] == '\n') extra_attrs[ext_att_len-1] =0;
						SUBREADprintf("\nERROR: failed to find the gene identifier attribute in the 9th column of the provided GTF file.\nThe specified gene identifier attribute is '%s'.\nAn example of attributes included in your GTF annotation is '%s'.\nThe program has to terminate.\n\n",  gene_id_column, extra_attrs);
					}
					is_GFF_geneid_warned++;
				}
					
				if(transcript_id_column && !is_tx_id_found){
					if(!is_GFF_txid_warned){
						int ext_att_len = strlen(extra_attrs);
						if(extra_attrs[ext_att_len-1] == '\n') extra_attrs[ext_att_len-1] =0;
						SUBREADprintf("\nERROR: failed to find the transcript identifier attribute in the 9th column of the provided GTF file.\nThe specified transcript identifier attribute is '%s'.\nAn example of attributes included in your GTF annotation is '%s'.\nThe program has to terminate\n\n", transcript_id_column, extra_attrs);
					}
					is_GFF_txid_warned++;
				}
			}
		}
		
		if(is_gene_id_found){
			do_add_feature(feature_name, transcript_id, chro_name, start, end, is_negative_strand, context);
			loaded_features++;
		}
		
	}
	autozip_close(&afp);
	free(file_line);

	if(is_GFF_txid_warned || is_GFF_geneid_warned)return -2;
	if(loaded_features<1){
		SUBREADprintf("\nERROR: No feature was loaded from the annotation file. Please check if the annotation format was correctly specified, and also if the feature type was correctly specified if the annotation is in the GTF format.\n\n");
		return -2;
	}
	return loaded_features;
}

HashTable * load_alias_table(char * fname) {
	FILE * fp = f_subr_open(fname, "r");
	if(!fp)
	{
		print_in_box(80,0,0,"WARNING unable to open alias file '%s'", fname);
		return NULL;
	}

	char * fl = malloc(2000);

	HashTable * ret = HashTableCreate(1013);
	HashTableSetDeallocationFunctions(ret, free, free);
	HashTableSetKeyComparisonFunction(ret, fc_strcmp_chro);
	HashTableSetHashFunction(ret, fc_chro_hash);
	
	while (1)
	{
		char *ret_fl = fgets(fl, 1999, fp);
		if(!ret_fl) break;
		if(fl[0]=='#') continue;
		char * sam_chr = NULL;
		char * anno_chr = strtok_r(fl, ",", &sam_chr);
		if((!sam_chr)||(!anno_chr)) continue;

		sam_chr[strlen(sam_chr)-1]=0;
		char * anno_chr_buf = malloc(strlen(anno_chr)+1);
		strcpy(anno_chr_buf, anno_chr);
		char * sam_chr_buf = malloc(strlen(sam_chr)+1);
		strcpy(sam_chr_buf, sam_chr);
		HashTablePut(ret, sam_chr_buf, anno_chr_buf);
	}

	fclose(fp);

	free(fl);
	return ret;
}

