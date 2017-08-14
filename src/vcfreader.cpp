#include "vcfreader.h"
#define __STDC_LIMIT_MACROS
#include "int_info.h"
//#include <cstdint> 
 #include <stdint.h> 

static VCFFileReader* VCF_FILE_ID = NULL;

void Open_VCF(char* VCF_File, int* myerror, size_t* size, int* form)
{
	
	VCF_FILE_ID = new  VCFFileReader(VCF_File, myerror, size, form); 

}


void Close_VCF() 
{
	delete VCF_FILE_ID;
}

void BCF_oneline_work( int *Z, int*err, int  size, char * chr, int * pos, char* snpid)
{
	VCF_FILE_ID->BCF_oneline( Z, err, size, chr, pos , snpid); //??

	
}

void BCF_oneline_work1( double *Z, int*err, int  size, char * chr, int * pos, char* snpid)
{
	VCF_FILE_ID->BCF_oneline1( Z, err, size, chr, pos , snpid); //??

	
}

extern "C" {


void R_BCF_oneline( int *Z, int * err, int * size, char* chr, int * pos,   char * snpid){

	BCF_oneline_work( Z, err, * size, chr, pos, snpid);
}

void R_BCF_oneline1( double *Z, int * err, int * size, char* chr, int * pos,   char * snpid){

	BCF_oneline_work1( Z, err, * size, chr, pos, snpid);
}

void R_Open_VCF(char** VCF_File, int * err, size_t * size, int * form)
{
	Open_VCF(VCF_File[0], err, size, form);

}

void R_Close_VCF() 
{
	Close_VCF();
}

/*
void BCF_search( char *regions, int * err) 
{
	VCF_FILE_ID->BCF_region(regions); 
}
*/	

}

/*

int vcfSetRegion( const char *regions) 
{
	return(VCF_FILE_ID->BCF_region(regions)); 
}

*/

void VCFFileReader::bcf_format_gt_new1(bcf_fmt_t *fmt, int isample, int * genotype)
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



void VCFFileReader::bcf_format_float(bcf_fmt_t *fmt, int isample, float* genotype)
{
#define BRANCH(type_t, missing, vector_end) { \
    type_t *ptr = (type_t*) (fmt->p + isample*fmt->size); \
    int i; *genotype=*ptr; \
}
    switch (fmt->type) {
        case BCF_BT_FLOAT:  BRANCH(float,  bcf_float_missing, bcf_float_vector_end); break;
        default: fprintf(stderr,"FIXME: type %d in bcf_format_gt?\n", fmt->type); abort(); break;
    }
#undef BRANCH
   
}



/* Test read VCF file */
int VCFFileReader::GetGenotype(bcf_hdr_t *hdr, bcf1_t *v, int * genotype ){


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
              
            this->bcf_format_gt_new1(f,j,genotype +j);

            n_geno_count += genotype[j];
        }
    }


    return n_geno_count;
    
}


int VCFFileReader::GetGenotype1(bcf_hdr_t *hdr, bcf1_t *v, double * genotype ){


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
            if (strcmp(hdr->id[BCF_DT_ID][fmt[i].id].key, "DS") == 0){
                gt_i = i;
                break;
            }
        }
        i = gt_i;
        for (j = 0; j < v->n_sample; ++j) {
            bcf_fmt_t *f = &fmt[i];
            if ( !f->p ) continue;
            float genotype_tmp;  
            this->bcf_format_float(f,j,&genotype_tmp);
	    genotype[j] = (double) genotype_tmp;
		
            n_geno_count += genotype[j];
        }
    
	}
	
	
    return n_geno_count;
    
}


/*************************************************************
 *
 * Moving window project
 * File: mwo_reader.h
 * Date: Dec 9, 2010
 * Author: Larissa Miropolsky
 *
 * Description:
 *   Deal with *.mwa file
 *
 * ver 1.1
 **************************************************************/


#include <bitset>
#include <fstream>
#include <iostream>
#include <cstring>
#include "mwo_reader.h"
#include "error_messages.h"






/*
int VCFFileReader::BCF_region(const char *regions_list){
/*     	Tabix file(this->m_filename);
	std::string reg(regions);
	file.setRegion(reg);
	std::string line;
        while (file.getNextLine(line)) {
                    cout << line << endl;

	}
	return ;


}


int SetRegion(char *fname, char * regions_list){
*/
/*
	
	//Close_bcf();
	bcf_srs_t * g_files;
	bcf_hdr_t * g_hdr;
	int g_nSample;
	g_files   = bcf_sr_init();

	if ( bcf_sr_set_regions(g_files, regions_list, 0)<0 ){
            printf("Failed to read the regions: %s\n", regions_list);
            exit(1);
    }

    if ( !bcf_sr_add_reader(g_files, this->m_filename) ){
    	printf("Failed to open %s: %s\n", this->m_filename, bcf_sr_strerror(g_files->errnum));
    	exit(1);
    }

 	g_hdr = g_files->readers[0].header;
 	bcf_hdr_sync(g_hdr);
 	
 	g_nSample = bcf_hdr_nsamples(g_hdr);
 	

 	printf("5\n");
	fflush(NULL);
 	return g_nSample;
 	
}
*/





//==========================================================
//Constructor
//Reads INFO file and takes all relevant info from there.
//Also based on INFO creates Offset table to move easily in ".mwa" file.
//Inputs:
// filename - path to ".mwa" file"
// myerror - allocated memory to put there the final information
//		if some error happen during the run.
// info - path to "INFO" file
//=======================================================================
VCFFileReader::VCFFileReader(char* filename, int* myerror, size_t * size, int * form)
{
	*myerror = NO_ERRORS;

	this->m_filename = filename;

	this->fp = hts_open(filename, "rb");

	if (!this->fp)	{
		*myerror = CANT_OPEN_MWA_FILE4READ;
		return;
	}

    	this->hdr= bcf_hdr_read(this->fp);

	this-> nSample = bcf_hdr_nsamples(this->hdr);//nSample == lines of fam file? ; 
	this-> v = bcf_init1();
	this->flag=0;
	this->form=form;
	size_t ss=this-> nSample;


	*size=ss;

	return;
    //===============================

}




//==========================================================
//Destructor - free all dynamically allocated memory
//==========================================================
VCFFileReader::~VCFFileReader()	{	
	int re;
	bcf_destroy1(v);
	bcf_hdr_destroy(hdr);
    
		//delete [] this->m_offsetarr;
		//delete [] this->m_set_size;
		//this->m_file.close();
   	if ( ( re=hts_close(this->fp)) )
    	{
        	fprintf(stderr,"hts_close: non-zero status %d\n",re);
        	exit(re);
    	}
}


/*
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

extern "C" {
SEXP	R_SetRegion_search( SEXP regions_list_){
	
		SEXP re_R;
		int re;
		const char * regions_list;
		
		//PROTECT(fname_ = AS_CHARACTER(fname_));
		PROTECT(regions_list_ = AS_CHARACTER(regions_list_));
		//fname = CHAR(STRING_ELT(fname_, 0));
		regions_list = CHAR(STRING_ELT(regions_list_, 0));
		
		//printf("[%s][%s]", fname, regions_list);
		//return R_NilValue;
		
		re =  vcfSetRegion( (char *)regions_list);
		
		
		re_R = PROTECT(ScalarInteger(re));
		
		UNPROTECT(3);
		return re_R;
}
}
*/





//
//
//	#include "PlinkBed.h" should be ahead of  #include <Rdefines.h>. There will be error otherwise. 
#include "zlib.h"
#include "zconf.h"
#include "hts.h"
#include "vcf.h"
#include "kstring.h"
#include "kseq.h"
#include "khash_str2int.h"
#include <algorithm>
#include <sstream>

#include <fstream>
#include <iterator>
#include <vector>
#include <iostream>
#include <cstring>
#include <math.h>

int Get_Genotype(bcf_hdr_t *hdr, bcf1_t *v, int * genotype );

// Added by SLEE
#include "synced_bcf_reader.h"

bcf_srs_t * g_files;
bcf_hdr_t * g_hdr;
int	g_nSample;

int Close_bcf(){

	if(g_hdr != NULL){
	
		//printf("1\n"); 
		fflush(NULL);
		bcf_hdr_destroy(g_hdr);
		g_hdr=NULL;
	}
	
	if(g_files != NULL){
		//printf("2\n"); 
		fflush(NULL);
		//bcf_sr_destroy(g_files); // this produces an error but don't know why.
		g_files=NULL;
	}
		//printf("3\n"); 
		fflush(NULL);	
}


int VCFFileReader::SetRegion(char *fname, char * regions_list){


	
	Close_bcf();
	
	g_files   = bcf_sr_init();

	if ( bcf_sr_set_regions(g_files, regions_list, 0)<0 ){
            printf("Failed to read the regions: %s\n", regions_list);
            exit(1);
    }

    if ( !bcf_sr_add_reader(g_files, fname) ){
    	printf("Failed to open %s: %s\n", fname, bcf_sr_strerror(g_files->errnum));
    	exit(1);
    }

 	g_hdr = g_files->readers[0].header;
 	bcf_hdr_sync(g_hdr);
 	
 	g_nSample = bcf_hdr_nsamples(g_hdr);
 	

 	//printf("5\n");
	fflush(NULL);

	this->flag=1;

 	return g_nSample;
 	
}

int	Get_Next_Genotypes(int * genotype, int * pos){
	
	// v->n_sample
	int re=-1;
	int nSample = 0;
	//printf("RunGenotype1", re, nSample); fflush(NULL);
	
	if(bcf_sr_next_line(g_files)){
  		bcf1_t *line = g_files->readers[0].buffer[0];
  		bcf_unpack(line, BCF_UN_ALL);
  		//printf("RunGenotype2", re, nSample); fflush(NULL);
	
        if ( line->errcode ){
        
        	printf("Undefined tags in the header, cannot proceed in the sample subset mode.\n");
        	exit(1);
        }

		*pos = line->pos + 1;	
		re = Get_Genotype(g_hdr, line, genotype );
		
		printf("pos [%d]", *pos); fflush(NULL);
	} 
	
	//printf("RunGenotype %d, %d\n", re, nSample); fflush(NULL);
	return nSample;

}




void VCFFileReader::BCF_oneline(int* Z, int* myerror, int  size , char * chr, int * pos, char * snpid){
     
    	int nsample= size;
	if (this->flag==0){
	if (bcf_read1(this->fp, this->hdr, this->v)>=0 ) { 
		bcf_unpack(this->v, BCF_UN_ALL);
	
		
		strncpy(chr, this->hdr->id[BCF_DT_CTG][this->v->rid].key, 100-1);		
        	*pos = v->pos + 1;
		
		strncpy(snpid , this->v->d.id, SNP_ID_SIZE_MAX-1);

		//std::cout <<snpid << endl;
		memset(Z, 0, sizeof(int)* nsample);
		
		int n_geno_count = this->GetGenotype(this->hdr, this->v, Z );

	}
		
	} 
	if (this->flag==1){
	// v->n_sample
	int re=-1;
	
	//printf("RunGenotype1", re, nSample); fflush(NULL);
	
	if(bcf_sr_next_line(g_files)){
  		bcf1_t *line = g_files->readers[0].buffer[0];
  		bcf_unpack(line, BCF_UN_ALL);
  		//printf("RunGenotype2", re, nSample); fflush(NULL);
	
        if ( line->errcode ){
        
        	printf("Undefined tags in the header, cannot proceed in the sample subset mode.\n");
        	exit(1);
        }
		strncpy(chr, g_hdr->id[BCF_DT_CTG][line->rid].key, 100-1);
		*pos = line->pos + 1;	
		strncpy(snpid , line->d.id, SNP_ID_SIZE_MAX-1);
		memset(Z, 0, sizeof(int)* nsample);
		
		int n_geno_count = this->GetGenotype(this->hdr, this->v, Z );


	} 
	



	}

	return ;


}



void VCFFileReader::BCF_oneline1(double* Z, int* myerror, int  size , char * chr, int * pos, char * snpid){
     
    	int nsample= size;
	if (this->flag==0){
	if (bcf_read1(this->fp, this->hdr, this->v)>=0 ) { 
		bcf_unpack(this->v, BCF_UN_ALL);
	
		
		strncpy(chr, this->hdr->id[BCF_DT_CTG][this->v->rid].key, 100-1);		
        	*pos = v->pos + 1;
		
		strncpy(snpid , this->v->d.id, SNP_ID_SIZE_MAX-1);

		//std::cout <<snpid << endl;
		memset(Z, 0, sizeof(int)* nsample);

		int n_geno_count = this->GetGenotype1(this->hdr, this->v, Z );
	}
		
	} 
	if (this->flag==1){
	// v->n_sample
	int re=-1;
	
	//printf("RunGenotype1", re, nSample); fflush(NULL);
	
	if(bcf_sr_next_line(g_files)){
  		bcf1_t *line = g_files->readers[0].buffer[0];
  		bcf_unpack(line, BCF_UN_ALL);
  		//printf("RunGenotype2", re, nSample); fflush(NULL);
	
        if ( line->errcode ){
        
        	printf("Undefined tags in the header, cannot proceed in the sample subset mode.\n");
        	exit(1);
        }
		strncpy(chr, g_hdr->id[BCF_DT_CTG][line->rid].key, 100-1);
		*pos = line->pos + 1;	
		strncpy(snpid , line->d.id, SNP_ID_SIZE_MAX-1);
		memset(Z, 0, sizeof(int)* nsample);

		int n_geno_count = this->GetGenotype1(this->hdr, this->v, Z );

	} 
	



	}
	
	return ;


}

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <string.h>

extern "C" {
	SEXP	R_SetRegion(SEXP fname_, SEXP regions_list_){
	
		SEXP re_R;
		int re;
		const char * fname, * regions_list;
		
		PROTECT(fname_ = AS_CHARACTER(fname_));
		PROTECT(regions_list_ = AS_CHARACTER(regions_list_));
		fname = CHAR(STRING_ELT(fname_, 0));
		regions_list = CHAR(STRING_ELT(regions_list_, 0));
		
		//printf("[%s][%s]", fname, regions_list);
		//return R_NilValue;
		
		re =  VCF_FILE_ID->SetRegion((char *) fname, (char *)regions_list);
		
		
		re_R = PROTECT(ScalarInteger(re));
		
		UNPROTECT(3);
		return re_R;
	}
	
	SEXP R_Close_bcf(){
	
		Close_bcf();
		return R_NilValue;
	}
	
	SEXP Get_Next_Genotypes(){
		
		int i;
		int * genotype_, * pos_;
		SEXP genotype, pos, vec;
		
		PROTECT(genotype = NEW_INTEGER(g_nSample));
		PROTECT(pos = NEW_INTEGER(1));
		vec = PROTECT(allocVector(VECSXP, 2));
		
		genotype_ = INTEGER_POINTER(genotype);
		pos_ = INTEGER_POINTER(pos);
		Get_Next_Genotypes(genotype_, pos_);
		
		SET_VECTOR_ELT(vec, 0, genotype);
  		SET_VECTOR_ELT(vec, 1, pos);

		UNPROTECT(3);
		return vec;
		
	}





}

