
#ifndef _BGEN_READER_H        
#define _BGEN_READER_H 



#define __STDC_LIMIT_MACROS
#include "int_info.h"
 #include <stdint.h> 

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
#include "synced_bcf_reader.h"


#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <memory>

//#include "genfile/View.hpp"
//#include "genfile/IndexQuery.hpp"
#include <sstream>
#include <time.h>
#include "bgen.hpp"



using namespace std;
using namespace genfile::bgen;

typedef uint8_t byte_t ;
//int vcfSetRegion( const char *regions); 

//===============================================================
//===============================================================


class BGENFileReader
{	

public:

	BGENFileReader(char* filename, int* myerror, size_t * size );
	~BGENFileReader();



	void BGEN_oneline(  double* Z, int* myerror, int  size, char * chr, int * pos, char * snpid);
	//int BCF_region(const char *regions);

private:
	
    	char* m_filename;

	size_t nSample; 

	
	std::auto_ptr< std::istream > gm_stream;
	uint32_t  gm_offset ;
	//genfile::bgen::Context gm_context ;

	double markerInfo;
};



	

#endif //
