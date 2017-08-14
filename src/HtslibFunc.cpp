//
//
//	#include "PlinkBed.h" should be ahead of  #include <Rdefines.h>. There will be error otherwise. 
#define __STDC_LIMIT_MACROS
//#include <cstdint> 
#include "int_info.h"
 #include <stdint.h> 
 #include <limits>
#include "zlib.h"
#include "zconf.h"
#include "hts.h"
#include "vcf.h"
#include "kstring.h"
#include "kseq.h"
#include "khash_str2int.h"
#include "PlinkBed.h"
#include <algorithm>
#include <sstream>
#include "SNP_info.h" 

#include <fstream>
#include <iterator>
#include <vector>
#include <iostream>
#include <cstring>
#include <math.h>




int CountLines(char *fam)  
{  
    std::ifstream ReadFile;  
    int n=0;  
    std::string tmp;  
    ReadFile.open(fam,ios::in);//ios::in 
    if(ReadFile.fail())// 
    {  
        return 0;  
    }  
    else
    {  
        while(getline(ReadFile,tmp,'\n'))  
        {  
            n++;  
        }  
        ReadFile.close();  
        return n;  
    }  
} 
 

static inline void bcf_format_gt_new(bcf_fmt_t *fmt, int isample, int * genotype)
{
#define BRANCH(type_t, missing, vector_end) { \
    type_t *ptr = (type_t*) (fmt->p + isample*fmt->size); \
    int i; *genotype=0;\
    for (i=0; i<fmt->n && ptr[i]!=vector_end; i++) \
    { \
        if ( !(ptr[i]>>1) ) *genotype=9; \
        else *genotype +=(ptr[i]>>1) - 1; \
    } \
    if (i == 0) *genotype=9;\
}
    switch (fmt->type) {
        case BCF_BT_INT8:  BRANCH(int8_t,  bcf_int8_missing, bcf_int8_vector_end); break;
        case BCF_BT_INT16: BRANCH(int16_t, bcf_int16_missing, bcf_int16_vector_end); break;
        case BCF_BT_INT32: BRANCH(int32_t, bcf_int32_missing, bcf_int32_vector_end); break;
        default: fprintf(stderr,"FIXME: type %d in bcf_format_gt?\n", fmt->type); abort(); break;
    }
#undef BRANCH
}


/* Test read VCF file */
int Get_Variant_Info(const char *fname){

    int i, re;
    htsFile *fp = hts_open(fname, "rb");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    kstring_t s1;
    kstring_t * s = &s1;
    
    int count=0;
    bcf1_t *v    = bcf_init1();

    while ( bcf_read1(fp, hdr, v)>=0 ){
        
        
        bcf_unpack(v, BCF_UN_FLT);
        count++;
        if(count % 1000 ==  0){
            printf("[%d]\t%s\t%d\n", count, hdr->id[BCF_DT_CTG][v->rid].key, v->pos + 1);
            fflush(NULL);
           // break;

        }
    }
    printf("[%d]\t%s\t%d\n", count, hdr->id[BCF_DT_CTG][v->rid].key, v->pos + 1);
    bcf_destroy1(v);
    bcf_hdr_destroy(hdr);
    
  	  
    if ( (re=hts_close(fp)) )
    {
        fprintf(stderr,"hts_close(%s): non-zero status %d\n",fname,re);
        exit(re);
    }
    return 1;

    
}



/* Test read VCF file */
int Get_VCF_Header(const char *fname){
    
    int i, re;
    htsFile *fp = hts_open(fname, "r");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    hts_idx_t *idx;
    hts_itr_t *iter;
    kstring_t s1;
    kstring_t * s = &s1;
    
    
    bcf1_t *v    = bcf_init1();
    re = bcf_read1(fp, hdr, v);
    bcf_unpack(v, BCF_UN_ALL);
    int gen=0;
    if (v->n_sample)
    {
        int i,j;
        if ( v->n_fmt)
        {
            int gt_i = -1;
            bcf_fmt_t *fmt = v->d.fmt;
            int first = 1;
            for (i = 0; i < (int)v->n_fmt; ++i) {
                if ( !fmt[i].p ) continue;
                kputc(!first ? ':' : '\t', s); first = 0;
                if ( fmt[i].id<0 ) //!bcf_hdr_idinfo_exists(h,BCF_HL_FMT,fmt[i].id) )
                {
                    fprintf(stderr, "[E::%s] invalid BCF, the FORMAT tag id=%d not present in the header.\n", __func__, fmt[i].id);
                    abort();
                }
                kputs(hdr->id[BCF_DT_ID][fmt[i].id].key, s);
                if (strcmp(hdr->id[BCF_DT_ID][fmt[i].id].key, "GT") == 0) gt_i = i;
            }
            if ( first ) kputs("\t.", s);
            for (j = 0; j < v->n_sample; ++j) {
                kputc('\t', s);
                first = 1;
                for (i = 0; i < (int)v->n_fmt; ++i) {
                    bcf_fmt_t *f = &fmt[i];
                    if ( !f->p ) continue;
                    if (!first) kputc(':', s); first = 0;
                    if (gt_i == i){
                        bcf_format_gt_new(f,j,&gen);
                        printf("[%d]",gen);
                    } else {
                        bcf_fmt_array(s, f->n, f->type, f->p + j * f->size);
                    }
                }
                if ( first ) kputc('.', s);
            }
        }
        printf("\n");
        fflush(NULL);
    }

    bcf_destroy1(v);
    bcf_hdr_destroy(hdr);
    

    if ( (re=hts_close(fp)) )
    {
        fprintf(stderr,"hts_close(%s): non-zero status %d\n",fname,re);
        exit(re);
    }
    return 1;
    
}



/* Test read VCF file */
int Get_Genotype(bcf_hdr_t *hdr, bcf1_t *v, int * genotype ){


    int i,j;
    int n_geno_count=0;
    if ( v->n_fmt){
        
        int gt_i = -1;
        bcf_fmt_t *fmt = v->d.fmt;
        int first = 1;
        for (i = 0; i < (int)v->n_fmt; ++i) {
            if ( !fmt[i].p ) continue;
            if ( fmt[i].id<0 ) //!bcf_hdr_idinfo_exists(h,BCF_HL_FMT,fmt[i].id) )
            {
                fprintf(stderr, "[E::%s] invalid BCF, the FORMAT tag id=%d not present in the header.\n", __func__, fmt[i].id);
                abort();
            }
            if (strcmp(hdr->id[BCF_DT_ID][fmt[i].id].key, "GT") == 0){
                gt_i = i;
                break;
            }
        }
        i = gt_i;
        for (j = 0; j < v->n_sample; ++j) {
            bcf_fmt_t *f = &fmt[i];
            if ( !f->p ) continue;
                
            bcf_format_gt_new(f,j,genotype +j);
            n_geno_count += genotype[j];
        }
    }


    return n_geno_count;
    
}





void decode_byte(int m_line_counter, std::vector<SNP_info1> & m_snp_sets, int* bits_val,size_t * individuals_counter, int* temp_snp_info0, int* temp_snp_info1, size_t snp_set_ind) 
//m_line_counter,  m_snp_sets
{


	int flag = 0;
	for (int i = 0; i < 4; ++i)
	{
		if(bits_val[i*2] == 0 && bits_val[i*2+1] == 0)
		{
			*individuals_counter += 1;
			if (*individuals_counter > m_line_counter)
				return;
			m_snp_sets[snp_set_ind].total_counter_per_letter[0] += 2;
			temp_snp_info0[*individuals_counter - 1] = 2;
			temp_snp_info1[*individuals_counter - 1] = 0;
			flag = 1 ;//Homozegote 1 for example GG      write 20 ; 00 - will point to [0] letter   + 2 to [0]

		}
		else if(bits_val[i*2] == 1 && bits_val[i*2+1] == 1)
		{
			*individuals_counter += 1;
			if (*individuals_counter > m_line_counter)
				return;
			m_snp_sets[snp_set_ind].total_counter_per_letter[1] += 2;
			temp_snp_info0[*individuals_counter - 1] = 0;
			temp_snp_info1[*individuals_counter - 1] = 2;
			flag = 2 ;//Homozegote 2 for example AA      write 02 ; 11 - will point to [1] letter   + 2 to [1]
		}
		else if(bits_val[i*2] == 1 && bits_val[i*2+1] == 0)
		{
			*individuals_counter += 1;		
			if (*individuals_counter > m_line_counter)
				return;
			temp_snp_info0[*individuals_counter - 1] = 9;
			temp_snp_info1[*individuals_counter - 1] = 9;
			flag = 3 ; //Missing value                   nothing to add - write 99 ;
		}
		else if(bits_val[i*2] == 0 && bits_val[i*2+1] == 1)
		{
			*individuals_counter += 1;
			if (*individuals_counter > m_line_counter)
				return;
			m_snp_sets[snp_set_ind].total_counter_per_letter[0] ++;
			m_snp_sets[snp_set_ind].total_counter_per_letter[1] ++;
			temp_snp_info0[*individuals_counter - 1] = 1;
			temp_snp_info1[*individuals_counter - 1] = 1;

			flag = 4 ; //Heterozegote for example AG  or GA     write 11 ; 01 - will point to [0] and [1] letter   +1 +1 to [1]
		}
		else
			flag = 5 ; //Error
	}

    //printf("Dec[%d]%d-%d\n",snp_set_ind,  m_snp_sets[snp_set_ind].total_counter_per_letter[0], m_snp_sets[snp_set_ind].total_counter_per_letter[1]);
}


void encode(int* temp_snp_info,char* encoded_snp_info, int m_line_counter , int m_size_of_esi)
{

	size_t i = 0;
	unsigned int j = 0;
	int number = 0;
	size_t ind4enc = 0;
	int a[8]; 

	//=======================================================================================
	// converting every 4 bytes of temp_snp_info to "a" 
	// "a"- array of 8 integers ("kind of" 8 bits) that will be converted to one encoded number
	// possible values for temp_snp_info = 9,1,2,0
	while (i < m_line_counter)
	{
		memset(a, 0, sizeof(a)); 
		for (j = 0; j < 4; ++j)
		{	
			if (temp_snp_info[i] == 9) //missing value
			{
				a[j*2] = 1;
				a[j*2+1] = 0;
			}
			else if (temp_snp_info[i] == 1)//AG
			{
				a[j*2] = 0;
				a[j*2+1] = 1;
			}
			else if (temp_snp_info[i] == 0)//GG
			{
				a[j*2] = 0;
				a[j*2+1] = 0;
			}
			else if (temp_snp_info[i] == 2)//AA
			{
				a[j*2] = 1;
				a[j*2+1] = 1;
			}
			i ++;
		}

		ind4enc ++;
		if (ind4enc == m_size_of_esi+1)
			break; 
		else
		{	
			number = 0;
			//==============================================
			//converting 8 ints of "a" to one encoded number
			//for example a= 00100010 "number" will be: 34
			for (int  ii = 0; ii < 8; ++ii)
				number += a[ii] * (int)pow(2.0,(7-ii));
			//saving this encoded number to array that will be written into *.mwa file.
			encoded_snp_info[ind4enc-1] = (char)number;
			//=============================================
		}
	}
}


int Convert_BCF_to_SSD_step2(const char *ssd, const char *PlinkPrefix,  std::vector<std::string> & SetID_name){

    char bim[256]; // <- danger, only storage for 256 characters.
    strncpy(bim, PlinkPrefix, sizeof(bim));
    strncat(bim, ".bim", sizeof(bim));
    char fam[256]; // <- danger, only storage for 256 characters.
    strncpy(fam, PlinkPrefix, sizeof(fam));
    strncat(fam, ".fam", sizeof(fam));
    char bed[256]; // <- danger, only storage for 256 characters.
    strncpy(bed, PlinkPrefix, sizeof(bed));
    strncat(bed, ".bed", sizeof(bed));
    std::string delimiter=" ";
    std::vector<std::string> token;
    ifstream in2(bim);
    std::string line;
    std::vector<SNP_info1> m_snp_sets;
    int index=0;
    std::vector<std::string> snp_exist;
    std::map<string,int> maplive;  

    if(in2)   
    {   int i=0; 

        while (getline (in2, line))   
        {   
            
            gTokenize(line, token,delimiter);
            std::string snp1 =token[1];
            SNP_info1 temp(token[1],token[4],token[5],token[0],token[3],atoi(token[2].c_str()),i);
            m_snp_sets.push_back(temp);
            snp_exist.push_back(snp1);
            maplive.insert(pair<string,int>(token[1],index));
            index++;
            token.clear();  
            i++;
        }  
    }  

    //std::sort(m_snp_sets.begin(), m_snp_sets.end() );///
    //std::sort(snp_exist.begin(), snp_exist.end() );///


    ifstream infile; 
    infile.open(bed, ios::binary); 
    ofstream out_ssd;
    out_ssd.open(ssd,ios::out | ios::binary);
    char ssdinfo[256]; // <- danger, only storage for 256 characters.
    strncpy(ssdinfo, ssd, sizeof(ssdinfo));
    strncat(ssdinfo, "infotemp", sizeof(ssdinfo));

    char ssdinfo_re[256]; // <- danger, only storage for 256 characters.
    strncpy(ssdinfo_re, ssd, sizeof(ssdinfo_re));
    strncat(ssdinfo_re, "info", sizeof(ssdinfo_re));


    ofstream out_info(ssdinfo);
    ofstream out_info_re(ssdinfo_re);


    int m_line_counter=CountLines(fam);
    int esi=(m_line_counter+3)/4;
    //cout<<"esi="<<esi<<endl;
    char* buff= new char [esi]; 
    int m_num_of_snps_insetid= CountLines(bim);


    int m_MAFConvert=1;

    out_info<< "-999" << "\tWindowSize" << std::endl;
    out_info<< m_MAFConvert << "\tMAFConvert" << std::endl;
    out_info<<  m_num_of_snps_insetid << "\tNumberOfSNPs" << std::endl;
    out_info<< m_line_counter<< "\tNumberOfIndividuals" << std::endl;

    out_info_re<< "-999" << "\tWindowSize" << std::endl;
    out_info_re<< m_MAFConvert << "\tMAFConvert" << std::endl;
    out_info_re<< m_num_of_snps_insetid << "\tNumberOfSNPs" << std::endl;
    out_info_re<< m_line_counter << "\tNumberOfIndividuals" << std::endl;

	
	size_t ii;
	//int count_win_size = 1;
	int MY_CHAR_BIT=8;
	int bits_val[MY_CHAR_BIT];
	int* temp_snp_info0 = new int[m_line_counter];
	int* temp_snp_info1 = new int[m_line_counter];
	for (ii = 0; ii < m_line_counter; ++ii)
	{
		temp_snp_info0[ii] = 0; 
		temp_snp_info1[ii] = 0; 
	}

	size_t individuals_counter = 0;
	char* encoded_snp_info = new char[(m_line_counter+3)/4] ;
	int m_size_of_esi = (m_line_counter+3)/4;  // !!!Number of bytes per one snp = one line length!!!
//	char* buff = new char [m_size_of_esi]; 
	char tmpbuf[3]; memset(tmpbuf, '\0', sizeof(tmpbuf));
	infile.read(tmpbuf,sizeof(tmpbuf)); // three first bytes - permanent in bed file

	//=========OFFSET===================
	size_t set_counter = 1; //   Will be changed based on "m_setidf_setid" if it's change - new set!
	size_t begin, current; 
	begin = out_ssd.tellp();
	current = out_ssd.tellp();
	//this->m_info << "#=================================================#" << std::endl;
	out_info <<"SET#\tOFFSET\tSET_ID\tSET_SIZE" << std::endl;
	size_t m_begin4rw = out_info.tellp();



	//=========================================
	// 
	//Go in loop "i = 0 ... ht->m_num_of_snps_insetid-1" toward "ht->m_setidf_setid" &&
	//based on ht->m_hash_table go by OFFSET in this->m_file , read relevant line: "this->m_size_of_esi" bytes ! 
	//put the readed bytes into "encoded_snp_info"
	//

	size_t setSize = 0;

	for (size_t j = 0; j < m_num_of_snps_insetid; ++ j)
	{

 		 
		if (j == 0)
		{
			//start new set
			out_info << set_counter << "\t" << (current - begin) << "\t" << SetID_name[j]; 
			set_counter ++;
		}
		else if ( SetID_name[j]!= SetID_name[j - 1])
		{	//start new set

			out_ssd << std::endl; //ENTER - empty line between every two snp sets.
			current = out_ssd.tellp();
			out_info << "\t" << setSize << std::endl;
			out_info <<  set_counter << "\t" << (current - begin) << "\t" << SetID_name[j]; // << std::endl;
			set_counter ++;
			setSize = 0;
                        
            
		}

		buff[0] = '\0';
			//moving inside of *.bed to reach specified location of specified snp - based on lookup table: ht->m_hash_table[j]
		//then read from there exactly one line 
    	       

        	int pos_bim=m_snp_sets[j].flag;

	        //int pos_bim=m_snp_sets[pos_temp].flag;
		infile.seekg(m_size_of_esi *pos_bim + 3 ,std::ios::beg); // +3 because of first 3 bytes in the file
		infile.read(buff,m_size_of_esi); 
		setSize ++;

        	std::string temp_id = snp_exist[j];
		out_ssd << temp_id << " ";

		for(size_t i=0; i<m_size_of_esi; i++)   //process the buff info
		{	//===============================================================					
			//=== This part converts Byte "buff[i]" to bits values "bits_val"
			//=== for example byte buff[0] = "w" ->  bits_val = 11101110
			memset(bits_val, NULL, sizeof(bits_val));
			int k = MY_CHAR_BIT;  //8
			while (k > 0)
			{
				-- k;
				bits_val[k] = (buff[i]&(1 << k) ? 1 : 0);
			}
			//==========================================================
			//=== here interpret Bit information "bits_val" to snps and count it - decode it
			decode_byte(m_line_counter, m_snp_sets,bits_val, &individuals_counter, temp_snp_info0, temp_snp_info1, pos_bim);

		} //end of for(int i=0; i<this->m_size_of_esi; i++)   //process the buff info

		{
			//Check who is MAGORITY/MINORITY
			//write data to output file from temp_snp_info0, temp_snp_info1)
			//encode temp_snp_info0, temp_snp_info1;

			// ENCODE OUTPUT
			//if (this->m_encode_output == 1)
			{
				memset(encoded_snp_info,0,sizeof(encoded_snp_info));
				
                //printf("%d:%d-%d\n",ht->m_hash_table[j],  m_snp_sets[ht->m_hash_table[j]].total_counter_per_letter[0], this->m_snp_sets[ht->m_hash_table[j]].total_counter_per_letter[1]);
                /* If m_MAFConvert==1 && 0 allele is the major allele */
                /* Very important!*/
                /* Be careful, total_counter_per_letter[0] < total_counter_per_letter[1] indicates 0 is the minor allele, so
                 we need to encode for 0 (add minor allele vector for encode) */
                
                if (m_snp_sets[pos_bim].total_counter_per_letter[0] < m_snp_sets[pos_bim].total_counter_per_letter[1])
				{
					//write snp information as encoded
                    // Add minor allele 
					encode(temp_snp_info0,encoded_snp_info, m_line_counter ,  m_size_of_esi );
					for (ii = 0; ii < m_size_of_esi; ++ii)
					{
						out_ssd << encoded_snp_info [ii] ; 
					}
				} else
				{
					encode(temp_snp_info1,encoded_snp_info, m_line_counter , m_size_of_esi );
					for (ii = 0; ii < m_size_of_esi; ++ii)
						out_ssd << encoded_snp_info [ii] ; 
				}
			}// END OF ENCODE OUTPUT

			out_ssd << std::endl; //ENTER at the end of every line
    			individuals_counter = 0;
			for (ii = 0; ii < m_line_counter; ++ii)
			{
				temp_snp_info0[ii] = 0; 
				temp_snp_info1[ii] = 0; 
			}
			//snp_set_ind ++;
		}  //end of if (individuals_counter >= this->m_line_counter)
	} // end of for (int j = 0; j < ht->m_num_of_snps_insetid; ++ j)

	delete[] temp_snp_info0;
	delete[] temp_snp_info1;
	delete[] encoded_snp_info;

	out_ssd << std::endl; //ENTER at the end of every line
	out_ssd << '\0'; 
	out_ssd.close();

	out_info << "\t" << setSize  << std::endl;  // print last set size
	out_info << "#=================================================#" << std::endl;
	out_info << m_size_of_esi + 1 << 
		"\tDECODED NumberOfDECODEDbytesPerSNP(NotIncludesSnpID&SpaceAfter_Includes\\n) " << std::endl;
	out_info << set_counter - 1 << "\tTotalNumberOfSets" << std::endl;

	out_info_re << m_size_of_esi + 1 << 
		"\tDECODED NumberOfDECODEDbytesPerSNP(NotIncludesSnpID&SpaceAfter_Includes\\n) " << std::endl;
	out_info_re << set_counter - 1 << "\tTotalNumberOfSets" << std::endl;


	size_t m_set_counter = set_counter - 1;

	out_info << '\0';
	out_info.close();

	std::vector<std::string>().swap(snp_exist);
	std::vector<SNP_info1>().swap(m_snp_sets);
 	maplive.clear();

	//==========================================================
	//  REWRITE INFO FILE
	//==========================================================
	ifstream ssd_info; 
	ssd_info.open(ssdinfo);
	/*if (!ssd_info)
	{
		printf("Error: !\n");
		return;
	}*/

	ssd_info.seekg (m_begin4rw, std::ios::beg);
	//this->m_info_rewr << "#=================================================#" << std::endl;
	out_info_re <<"SET#\tOFFSET\tSET_ID\tSET_SIZE" << std::endl;


	size_t kk = 0;
	while (kk < m_set_counter ) 
	{
		getline(ssd_info, line);
		kk ++;
		out_info_re << line << std::endl;
	}
	out_info_re.close();

		



	return 1;
}

bool ValueCmp(refgene const & a, refgene const & b)
{
    return a.chr < b.chr;
}



int Convert_BCF_to_SSD(const char *ssd, CPlinkBed_Write * PlinkBed, const char *PlinkPrefix, const char *BCF_file, const char *SetID, int SetIDtype, int nmax){
    
    std::ifstream in(SetID);
    std::string line;  
    std::vector<std::string> snp;
    std::vector<std::string> SetID_name;
    std::map<std::string,std::string> map_setid;  
    std::vector<std::string> token;
    std::vector<std::string> tokens2;
    std::vector<std::string> tokens3;
    std::string delimiter=" ";
    std::string temp2;


    if (SetIDtype==1 || SetIDtype==2){	
    	while (getline (in, line))   
    	{   
          	token.clear();
    		gTokenize(line, token,"\t\n\r");
		if (token.size() < 2)
		{
			token.clear();
			gTokenize(line, token, " ");
		}
		if (token.size() >= 2)
		{
			tokens2.clear();
			gTokenize(token.at(1).c_str(), tokens2, " ");
			tokens3.clear();
			gTokenize(tokens2.at(0).c_str(), tokens3, "\r");
			
		       	std::string snp_temp=tokens3.at(0).c_str();
        		snp.push_back(snp_temp);
			map_setid.insert(pair<string,string>(snp_temp,token.at(0).c_str()));
		}
        		        
    	}
    }


    if (SetIDtype==3){
  
    	int pos_int1=0;
    	int pos_int2=0; 
	
    	while (getline (in, line))   
    	{   
                token.clear();
    		gTokenize(line, token,"\t\n\r");
		if (token.size() < 4)
		{
			token.clear();
			gTokenize(line, token, " ");
		}
		if (token.size() >= 4)
		{
			
			pos_int1 = atoi(token[2].c_str());
			pos_int2 = atoi(token[3].c_str());
	//printf("wulala1");
			for (int i=pos_int1;i<=pos_int2;i++){
	
				std::ostringstream temp1; 
				temp1<<i;
				std::string temp11=temp1.str();
				std::string snp_temp=token[1]+":"+temp11;
	        	       	snp.push_back(snp_temp);
//printf("wulala2");
				map_setid.insert(pair<std::string,std::string>(snp_temp,token[0]));
//printf("wulala3");
			}
		}
 
        	         
    	}
    }
    std::vector<refgene> refarray;
    if (SetIDtype==4){
  
    	int pos_int1=0;
    	int pos_int2=0; 

	//token[2]  [9] [10] [12];
	
	int countii=0;
    	while (getline (in, line))   
    	{   
                token.clear();
		tokens2.clear();
		tokens3.clear();
    		gTokenize(line, token,"\t\n\r");
	
		if (token.size() < 12)
		{
			token.clear();
			gTokenize(line, token, " ");
		}
		if (token.size() >= 12)
		{
			temp2=token[2].substr(3);
			//std::cout<<temp2<<"  ";
			// need to minus "chr";
			gTokenize(token[9], tokens2,",");
			gTokenize(token[10], tokens3,",");
			int tokensize1=tokens2.size();
			int tokensize2=tokens3.size();
			int tokensize= std::min(tokensize1,tokensize2); 
			//std::cout<<tokensize<<endl;
			for (int ii=0;ii<tokensize; ii++){
				pos_int1 = atoi(tokens2[ii].c_str());
				pos_int2 = atoi(tokens3[ii].c_str());
				refgene temp_ref(temp2,pos_int1,pos_int2,token[12]);
				refarray.push_back(temp_ref); 

				//std::cout<<pos_int1<<" "<<pos_int2<<endl;	
			/*	for (int i=pos_int1;i<=pos_int2;i++){
	
					std::ostringstream temp1; 
					temp1<<i;
					std::string temp11=temp1.str();
					std::string snp_temp=temp2+":"+temp11;
					if (map_setid.count(snp_temp) == 0){
	        	       			snp.push_back(snp_temp);
	
//printf("wulala2");
				 		map_setid.insert(pair<std::string,std::string>(snp_temp,token[12]));
					}
					countii++;
					if (countii%100000==0) {std::cout<<countii;}
				} */
		//		std::cout<<"I am done"<<endl;
//printf("wulala3");
			}//for ii
		//	std::cout<<"I am done1"<<endl;
		}//if >12
 		//std::cout<<"I am done2"<<endl;
        	         
    }//getline
   
   }//4		

  
    //std::cout<<"I am done"<<endl;


    int i,re;
    htsFile *fp = hts_open(BCF_file, "rb");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    hts_idx_t *idx;
    hts_itr_t *iter;
    kstring_t s1;
    kstring_t * s = &s1;
    
    size_t nSample = bcf_hdr_nsamples(hdr);//nSample == lines of fam file? ;
    
    
   	
    /* Init(const char * FilePrefix, bool AddLocalTime=true, size_t nSampleSize=0); */
    if(PlinkBed->Init(nSample) == 0){
        printf("Error: CPlinkBed Init!\n");
        return 0;
    }

    PlinkBed->Write_FamFile(hdr->samples);
    int * genotypes = new int[nSample];

    
    unsigned long count=0;
    bcf1_t *v    = bcf_init1();

    typedef std::string     kstring_t;

if (SetIDtype<4) {
    std::sort(snp.begin(), snp.end() );
    


}



std::vector<int> ref_num;



if (SetIDtype==4){
	int nn=1;
	std::sort(refarray.begin(), refarray.end(), ValueCmp);
	for (int i4=1; i4<refarray.size();i4++){
		if (refarray[i4].chr==refarray[i4-1].chr) {nn++;} else {
			ref_num.push_back(nn);
			//std::cout<<nn<<endl;///
			nn=1;
		}

	
	}


}

   	//char setid_default[256]; /// <- danger, only storage for 256 characters.
    	//strncpy(setid_default, PlinkPrefix, sizeof(setid_default));///
    	//strncat(setid_default, "setid.txt", sizeof(setid_default));///
    	//ofstream out(setid_default);///


    while ( bcf_read1(fp, hdr, v)>=0 ){

        std::string chr;
        size_t pos;
        std::string A1="";
        std::string A2="";
        
        bcf_unpack(v, BCF_UN_ALL);
        chr = hdr->id[BCF_DT_CTG][v->rid].key;
        pos = v->pos + 1;
        std::string snp1;
	if (SetIDtype>1 ){
		std::ostringstream ostr; //output string stream    	
		ostr << pos; //use the string stream just like cout,
    		std::string snp_pos = ostr.str();
		snp1=chr+":"+snp_pos;
	} else {
	
		snp1 = v->d.id;
		
	}
		
if (SetIDtype<4) {	

	if (std::binary_search(snp.begin(), snp.end(), snp1 )){
	    	map<std::string ,std::string >::const_iterator map_temp;
	        map_temp=map_setid.find(snp1);
	        std::string setid_temp=map_temp->second;
		SetID_name.push_back(setid_temp);
        	
        	if (v->n_allele > 0) A1=v->d.allele[0];
        	else A1=".";
        	if (v->n_allele > 1) {
	            	for (i = 1; i < v->n_allele; ++i) {
	                	A2 +=v->d.allele[i];
	            	}
        	} else A2=".";
        
	        memset(genotypes, 0, sizeof(int)* nSample);
	        int n_geno_count = Get_Genotype(hdr, v, genotypes );
	       // printf("[%d:%d:%d:%d:%d:%d]\n",count, n_geno_count, genotypes[0],genotypes[1],genotypes[2],genotypes[3]);
        
	        PlinkBed->Write_OneSNP1(snp1.c_str(),chr.c_str(), pos, A1.c_str(), A2.c_str(), genotypes);
        
	        count++;
	        if( count % 100000 == 0){
	            printf("%lu variants were read from the input file\n", count);
	        }
	        if(nmax > 0 && count >= nmax){
	            break;
	      	}
	}
}
	
	
if (SetIDtype==4) {	


 	int i_4=0;
	int j_4=0;
	int findit=0;
	//std::cout<<refarray.size()<<"size"<<endl;///
	while ( i_4 < refarray.size()){
		if (chr==refarray[i_4].chr ){
			for (int k_4=0;k_4<ref_num[j_4];k_4++){
				if (pos<=refarray[i_4+k_4].pos_end && pos>=refarray[i_4+k_4].pos_start){				
					//out<<SetID4[i_4]<<" "<<chr<<":"<<pos<<std::endl;///
					SetID_name.push_back(refarray[i_4+k_4].setid4);

					if (v->n_allele > 0) A1=v->d.allele[0];
        				else A1=".";
        				if (v->n_allele > 1) {
	            				for (i = 1; i < v->n_allele; ++i) {
	                				A2 +=v->d.allele[i];
	            				}
        				} else A2=".";
        
			        	memset(genotypes, 0, sizeof(int)* nSample);
			        	int n_geno_count = Get_Genotype(hdr, v, genotypes );
	        //printf("[%d:%d:%d]\n",count, n_geno_count, genotypes[0]);
        
				        PlinkBed->Write_OneSNP1(snp1.c_str(),chr.c_str(), pos, A1.c_str(), A2.c_str(), genotypes);
       			 		count++;
				        if( count % 100000 == 0){
					         printf("%lu variants were read from the input file\n", count);
				        }
					findit=1;
					break;
				} 
			}
		}
		if (findit==1){break;}
		i_4=i_4+ref_num[j_4];
		j_4++;
	}

	if(nmax > 0 && count >= nmax){
            break;
      	}


	
} 


}


   // out.close();///
   
    PlinkBed->Close_BedBim_Files();
    
    delete [] genotypes;
    

   
    bcf_destroy1(v);
    bcf_hdr_destroy(hdr);
    
    std::vector<std::string>().swap(snp);


    map_setid.clear();
    

    
     if ( (re=hts_close(fp)) )
    {
        exit(re);
   }


 


    std::vector<refgene>().swap(refarray);


    Convert_BCF_to_SSD_step2(ssd, PlinkPrefix, SetID_name);
    std::vector<std::string>().swap(SetID_name);
    return 1;
}





int Convert_BCF_to_PlinkBed(CPlinkBed_Write * PlinkBed, const char *BCF_file, int nmax){
   
    int i, re;
    htsFile *fp = hts_open(BCF_file, "rb");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    hts_idx_t *idx;
    hts_itr_t *iter;
    kstring_t s1;
    kstring_t * s = &s1;
    
    size_t nSample = bcf_hdr_nsamples(hdr);
    

   
    /* Init(const char * FilePrefix, bool AddLocalTime=true, size_t nSampleSize=0); */
    if(PlinkBed->Init(nSample) == 0){
        printf("Error: CPlinkBed Init!\n");
        return 0;
    }
    PlinkBed->Write_FamFile(hdr->samples);
    int * genotypes = new int[nSample];

    
    unsigned long count=0;
    bcf1_t *v    = bcf_init1();
    while ( bcf_read1(fp, hdr, v)>=0 ){

        std::string chr;
        size_t pos;
        std::string A1="";
        std::string A2="";
        
        bcf_unpack(v, BCF_UN_ALL);
        chr = hdr->id[BCF_DT_CTG][v->rid].key;
        pos = v->pos + 1;
        
        if (v->n_allele > 0) A1=v->d.allele[0];
        else A1=".";
        if (v->n_allele > 1) {
            for (i = 1; i < v->n_allele; ++i) {
                A2 +=v->d.allele[i];
            }
        } else A2=".";
        
        memset(genotypes, 0, sizeof(int)* nSample);
        int n_geno_count = Get_Genotype(hdr, v, genotypes );
        //printf("[%d:%d:%d]\n",count, n_geno_count, genotypes[0]);
        std::string snp1 = v->d.id;

        PlinkBed->Write_OneSNP(snp1.c_str(), chr.c_str(), pos, A1.c_str(), A2.c_str(), genotypes);
        
        count++;
        if( count % 100000 == 0){
            printf("%lu variants were read from the input file\n", count);
        }
        if(nmax > 0 && count >= nmax){
            break;
        }
    }
    PlinkBed->Close_BedBim_Files();
    
    
    delete [] genotypes;
    


   
    bcf_destroy1(v);
    bcf_hdr_destroy(hdr);
    
    
    if ( (re=hts_close(fp)) )
    {
        exit(re);
    }
    return 1;

    
}

extern 	"C" {
int Convert_BCF_to_PlinkBed_Work(char *PlinkPrefix, char *BCF_file, int nmax){

	
	CPlinkBed_Write PlinkBed;
	PlinkBed.Generate_FileName(PlinkPrefix, false, false);
    
    int re = Convert_BCF_to_PlinkBed(&PlinkBed, BCF_file,  nmax);
    
    return(re);

}
}









extern "C" {
int Convert_BCF_to_SSD_Work(char *ssd, char *PlinkPrefix, char *BCF_file, char *SetID, int SetIDtype, int nmax){

	
	CPlinkBed_Write PlinkBed;
	PlinkBed.Generate_FileName(PlinkPrefix, false, false);
    
    int re = Convert_BCF_to_SSD(ssd, &PlinkBed, PlinkPrefix, BCF_file, SetID, SetIDtype, nmax);

 
    return(re);

}
}

#include "mwo_reader.h"


extern "C" {

static MwoFileReader* MWA_FILE_ID = NULL;



void Open_MWA(char* MWA_File, char* Info, int* myerror)
{
	//it will take also "CEU.bed.INFO.txt" from same place - both of them shoul be existed
	//MWA_FILE_ID = new (MWA_File); 
	MWA_FILE_ID = new  MwoFileReader(MWA_File, myerror, Info); 

}


void Close_MWA() 
{
	delete MWA_FILE_ID;
}



void Get_Genotypes( int Set_number, int* Z, int size, int Is_MakeFile, int* myerror) // set_number base on INFO file. The result will be printed to file.
{
	MWA_FILE_ID->get_set(Set_number, Z, size, myerror, Is_MakeFile); //some integer - enter the value base on CEU.bed.INFO.txt
}	

void Get_Genotypes_withID( int Set_number, int* Z, char * SNPID, int size, int Is_MakeFile, int* myerror) // set_number base on INFO file. The result will be printed to file.
{
	MWA_FILE_ID->get_set(Set_number, Z, size, myerror, Is_MakeFile, SNPID); //some integer - enter the value base on CEU.bed.INFO.txt
}

}



