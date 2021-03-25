
#ifndef _VCF_READER_H        
#define _VCF_READER_H 


#include "vcfreader.h"
#define __STDC_LIMIT_MACROS
#include "int_info.h"
 #include <stdint.h> 

#include "zlib.h"
#include "zconf.h"
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "htslib/khash_str2int.h"
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
#include "htslib/synced_bcf_reader.h"



using namespace std;


//int vcfSetRegion( const char *regions); 

//===============================================================
//===============================================================


class VCFFileReader
{	
public:

	VCFFileReader(char* filename, int* myerror, size_t * size,int * form );
	~VCFFileReader();


	void bcf_format_gt_new1(bcf_fmt_t *fmt, int isample, int * genotype);
	int GetGenotype(bcf_hdr_t *hdr, bcf1_t *v, int * genotype );
	int GetGenotype1(bcf_hdr_t *hdr, bcf1_t *v, double * genotype );
	void BCF_oneline(  int* Z, int* myerror, int  size, char * chr, int * pos, char * snpid,char * A1, char *A2);
	void BCF_oneline1(  double* Z, int* myerror, int  size, char * chr, int * pos, char * snpid,char * A1, char *A2);
	//int BCF_region(const char *regions);
	int SetRegion(char *fname, char * regions_list);
	void bcf_format_float(bcf_fmt_t *fmt, int isample, float* genotype);
	void getid(char * samplesid);
	
private:
    	char* m_filename;
	htsFile *fp;
	bcf_hdr_t *hdr;
	size_t nSample; 
	bcf1_t *v ;
	int flag;
	int* form;


};




#endif //_MWO_READER_H

