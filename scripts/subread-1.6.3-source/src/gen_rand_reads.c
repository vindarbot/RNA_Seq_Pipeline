#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<assert.h>

#define gvindex_baseno2offset_m(base_number, index, offset_byte, offset_bit)    {offset_byte =  (base_number - index -> start_base_offset) >>2; offset_bit = base_number % 4 * 2;}
#define int2base(c) (1413695297 >> (8*(c))&0xff);



typedef struct{
        unsigned int start_base_offset;
        unsigned int start_point;
        unsigned int length;
        unsigned char * values;
        unsigned int values_bytes;
} gene_value_index_t;


void gvindex_baseno2offset(unsigned int base_number, gene_value_index_t * index, unsigned int * offset_byte, unsigned int * offset_bit)
{
        // the base number corrsponding to the 0-th bit in the whole value array;

        unsigned int offset = (base_number - index -> start_base_offset);

        * offset_byte = offset >>2 ;
        * offset_bit = base_number % 4 * 2;
}





int gvindex_load(gene_value_index_t * index, const char filename [])
{
        FILE * fp = f_subr_open(filename, "rb");
        int read_length;
        read_length = fread(&index->start_point,4,1, fp);
        assert(read_length>0);
        read_length = fread(&index->length,4,1, fp);
        assert(read_length>0);

        //SUBREADprintf ("\nBINDEX %s : %u ~ +%u\n",filename, index->start_point, index->length );

        unsigned int useful_bytes, useful_bits;
        index -> start_base_offset = index -> start_point - index -> start_point%4;
        gvindex_baseno2offset (index -> length+ index -> start_point, index ,&useful_bytes,&useful_bits);
        index -> values = malloc(useful_bytes);
        index -> values_bytes = useful_bytes;
	//printf("GVLD=%u\n", index-> values_bytes);
        if(!index->values)
        {
                return 1;
        }


        read_length =fread(index->values, 1, useful_bytes, fp);
        assert(read_length>0);

        fclose(fp);
        return 0;

}


int gvindex_get(gene_value_index_t * index, int offset)
{
        unsigned int offset_byte, offset_bit;
        gvindex_baseno2offset_m(offset, index , offset_byte, offset_bit);

        if(offset_byte >= index-> values_bytes){
		char nch_i = rand() % 4 ;
		char nch = nch_i?(nch_i<2?'T':(nch_i < 3?'C':'G')):'A';
		return nch;
	}

	//printf("GV=%u\n", index-> values_bytes);

        unsigned int one_base_value = (index->values [offset_byte]) >> (offset_bit);

        return int2base(one_base_value & 3);
}


char * __converting_char_table = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTNGNNNCNNNNNNNNNNNNAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTNGNNNCNNNNNNNNNNNNAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN  ";


void reverse_read(int read_len, char * InBuff)
{
	int i;
	for (i=0; i<read_len/2; i++)
	{
		int rll1 = read_len - 1 - i;
		unsigned char tmp = InBuff[rll1];

		InBuff[rll1] = *(__converting_char_table+InBuff[i]);
		InBuff[i] = *(__converting_char_table+tmp);

	}
	if(i*2 == read_len-1)
	{
		InBuff[i] = *(__converting_char_table+InBuff[i]);
	}

}

int main(int argc, char ** argv)
{
	gene_value_index_t bv;

	gvindex_load(&bv, argv[1]);

	int reads = atoi(argv[2]);
	int x1, x2;
	int read_len = 110;

	char read_txt[500];
	char qual_txt[500];

	for(x1=0; x1<reads; x1++)
	{
		unsigned long long int read_seed = 0;
		unsigned long long int var_seed = 0;
		for(x2=0; x2<8; x2++)
		{
			read_seed << 8;
			read_seed ^= rand();
			var_seed <<= 8;
			var_seed ^= rand();
		}

		unsigned long long int linear_pos = read_seed%bv.length;
		unsigned long long int linear_pos2 = var_seed%bv.length;
		unsigned long long int spl_point = ((read_seed^var_seed)% 9191919191) % 110;
			
		int wptr=0, x3;
		for(x2=0; x2<read_len; x2++)
		{
			if(wptr > 110)break;
			if(var_seed %(91919191-x2) < 119191)
			{
				for(x3 = rand()%5; x3>=0; x3--)
				{
					char nch_i = rand() % 4 ;
					char nch = nch_i?(nch_i<2?'T':(nch_i < 3?'C':'G')):'A';
					read_txt[wptr++]=nch;
				}
			}
			else
				if(var_seed %(61919191-x2) < 0*219191)
				{
					linear_pos+=rand() % 10 ;
					linear_pos2+=rand() % 10 ;
				}
			if(wptr<spl_point)
				read_txt[wptr++]=gvindex_get(&bv, linear_pos+x2);
			else
				read_txt[wptr++]=gvindex_get(&bv, linear_pos2+x2-spl_point);
		}

		read_txt[wptr]=0;

		if(rand()% 199 < 60)
		{
			if(rand()% 199 < 99)reverse_read(spl_point, read_txt);
			else reverse_read(wptr-spl_point, read_txt+spl_point);
		}

		for(x2=0; read_txt[x2]; x2++) qual_txt[x2]='I';
		qual_txt[x2]=0;

		printf("@read%d\n%s\n+\n%s\n", x1+1, read_txt, qual_txt);
		assert(strlen(read_txt) == strlen(qual_txt));
	}

	return 0;
}
