#include "bgen_reader.h"
#define __STDC_LIMIT_MACROS
#include "int_info.h"
//#include <cstdint> 
 #include <stdint.h> 

/*
static BGENFileReader* BGEN_FILE_ID = NULL;///////////////

void Open_BGEN(char* BGEN_File, int* myerror, size_t* size)
{
	
	BGEN_FILE_ID = new BGENFileReader(BGEN_File, myerror, size); 

}


void Close_BGEN() 
{
	delete BGEN_FILE_ID;
}

void BGEN_oneline_work( double *Z, int*err, int  size, char * chr, int * pos, char* snpid)
{::
	BGEN_FILE_ID->BGEN_oneline( Z, err, size, chr, pos , snpid); //??

	
}



extern "C" {


void R_BGEN_oneline( double *Z, int * err, int * size, char* chr, int * pos,   char * snpid){

	BGEN_oneline_work( Z, err, * size, chr, pos, snpid);
}



void R_Open_BGEN(char** VCF_File, int * err, size_t * size)
{
	Open_BGEN(VCF_File[0], err, size);

}

void R_Close_BGEN() 
{
	Close_BGEN();
}

}


*/

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
#include "error_messages.h"

namespace {
        template< typename Integer >
        std::string atoi( Integer const value ) {
                std::ostringstream stream ;
                stream << value ;
                return stream.str() ;
        }

}


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
/*
BGENFileReader::BGENFileReader(char* filename, int* myerror, size_t * size)
{}

 	int numMarkers;
   	
	using namespace genfile::bgen ;
    	using namespace Rcpp;

     	


	*myerror = NO_ERRORS;

	this->m_filename = filename;

	this->gm_stream.reset(
            new std::ifstream( filename, std::ifstream::binary )
          ) ;/net/hunt/zhowei/script/tools/bgen


          if( !*gm_stream ) {
            throw std::invalid_argument( filename ) ;
          }

          //printf("1\n");fflush(NULL);
          this->gm_stream->seekg( 0, std::ios::beg ) ;
          genfile::bgen::read_offset( *this->gm_stream, &this->gm_offset ) ;
          //printf("2\n");fflush(NULL);   
          genfile::bgen::read_header_block( *this->gm_stream, &this->gm_context ) ;
          // Jump to the first variant data block.
          this->gm_stream->seekg(this-> gm_offset + 4 ) ;
          //printf("4\n");fflush(NULL);
          uint Mbgen = this->gm_context.number_of_variants;
          *size = size_t(Mbgen);

          std::cout << "All " << size << " markers will be analyzed " << std::endl;
	this->nSample=size_t(Mbgen);
	return;          

}

*/

/*
//==========================================================
//Destructor - free all dynamically allocated memory
//==========================================================
BGENFileReader::~BGENFileReader()	{	
	 free (m_filename);
	//free (nSample);
	//free (gm_sample_ids );
	//free( gm_stream);
	//free(gm_offset );
	//free(gm_context);
	//free(gm_have_sample_ids );
	//free(genoToTest_bgenDosage);
	//free(markerInfo);
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





double  Parse(unsigned char * buf, size_t bufLen,  std::string & snpName, uint Nbgen,std::vector< double > & dosages){

    size_t destLen = bufLen;

    unsigned char * bufAt = buf;
    uint N = bufAt[0]|(bufAt[1]<<8)|(bufAt[2]<<16)|(bufAt[3]<<24); bufAt += 4;
    if (N != Nbgen) {
      std::cerr << "ERROR: " << snpName << " has N = " << N << " (mismatch with header block)" << std::endl;
      exit(1);
    }
    uint K = bufAt[0]|(bufAt[1]<<8); bufAt += 2;
    if (K != 2U) {
      std::cerr << "ERROR: " << snpName << " has K = " << K << " (non-bi-allelic)" << std::endl;
      exit(1);
    }
    uint Pmin = *bufAt; bufAt++;
    if (Pmin != 2U) {
      std::cerr << "ERROR: " << snpName << " has minimum ploidy = " << Pmin << " (not 2)" << std::endl;
      exit(1);
    }
    uint Pmax = *bufAt; bufAt++;
    if (Pmax != 2U) {
      std::cerr << "ERROR: " << snpName << " has maximum ploidy = " << Pmax << " (not 2)" << std::endl;
      exit(1);
    }
    for (uint i = 0; i < N; i++) {
      uint ploidyMiss = *bufAt; bufAt++;
      if (ploidyMiss != 2U) {
        std::cerr << "ERROR: " << snpName << " has ploidy/missingness byte = " << ploidyMiss
             << " (not 2)" << std::endl;
        exit(1);
      }
    }
    uint Phased = *bufAt; bufAt++;
    if (Phased != 0U) {
      std::cerr << "ERROR: " << snpName << " has Phased = " << Pmax << " (not 0)" << std::endl;
      exit(1);
    }
    uint B = *bufAt; bufAt++;
    if (B != 8U) {
      std::cerr << "ERROR: " << snpName << " has B = " << B << " (not 8)" << std::endl;
      exit(1);
    }

        // Parse 
    double lut[256];
    for (int i = 0; i <= 255; i++)
      lut[i] = i/255.0;

    double sum_eij = 0, sum_fij_minus_eij2 = 0; // for INFO
    double p11,p10,dosage,eij,fij;
    for (uint i = 0; i < N; i++) {
      p11 = lut[*bufAt]; bufAt++;
      p10 = lut[*bufAt]; bufAt++;
      	
      dosage = 2*p11 + p10;
      dosages.push_back(2 - dosage);
/*
    if(i < 3){
        std::cout << "p11: " << p11 << std::endl;
        std::cout << "p10: " << p10 << std::endl;
        std::cout << "dosage: " << 2-dosage << std::endl;
    }
*/
      eij = dosage;
      fij = 4*p11 + p10;
      sum_eij += eij;
      sum_fij_minus_eij2 += fij - eij*eij;
    }
    double thetaHat = sum_eij / (2*N);
    double info = thetaHat==0 || thetaHat==1 ? 1 :
      1 - sum_fij_minus_eij2 / (2*N*thetaHat*(1-thetaHat));
//    std::cout << "info: " << info << std::endl;



    return(info);
}

namespace bgen1 {
		// Read an integer stored in little-endian format into an integer stored in memory.
		template< typename IntegerType >
		byte_t const* read_little_endian_integer( byte_t const* buffer, byte_t const* const end, IntegerType* integer_ptr ) {
			assert( end >= buffer + sizeof( IntegerType )) ;
			*integer_ptr = 0 ;
			for( std::size_t byte_i = 0; byte_i < sizeof( IntegerType ); ++byte_i ) {
				(*integer_ptr) |= IntegerType( *reinterpret_cast< byte_t const* >( buffer++ )) << ( 8 * byte_i ) ;
			}
			return buffer ;
		}
		

		// Read an integer stored in little-endian format into an integer stored in memory.
		// The stream is assumed to have sizeof( Integertype ) readable bytes.
		template< typename IntegerType >
		void read_little_endian_integer( std::istream& in_stream, IntegerType* integer_ptr ) {
			byte_t buffer[ sizeof( IntegerType ) ] ;
			in_stream.read( reinterpret_cast< char* >( buffer ), sizeof( IntegerType )) ;
			if( !in_stream ) {
				throw genfile::bgen::BGenError() ;
			}
			read_little_endian_integer( buffer, buffer + sizeof( IntegerType ), integer_ptr ) ;
		}


		 

}



void genfile::bgen::read_offset( std::istream& iStream, uint32_t* offset ) {
	bgen1::read_little_endian_integer( iStream, offset ) ;
}

BGENFileReader::BGENFileReader(char* filename, int* myerror, size_t * size)
{

	*myerror = NO_ERRORS;

	this->m_filename = filename;

	this->gm_stream.reset(
            new std::ifstream( filename, std::ifstream::binary )
          ) ;

          if( !*this->gm_stream ) {
            throw std::invalid_argument( filename ) ;
          }

          //printf("1\n");fflush(NULL);
          this->gm_stream->seekg( 0, std::ios::beg ) ;

         read_offset( *this->gm_stream, &this->gm_offset ) ;
/*
          //printf("2\n");fflush(NULL);   
          genfile::bgen::read_header_block( *this->gm_stream, &this->gm_context ) ;
          // Jump to the first variant data block.
          this->gm_stream->seekg( this->gm_offset + 4 ) ;
          //printf("4\n");fflush(NULL);
          uint Mbgen = this->gm_context.number_of_variants;
          *size = size_t(Mbgen);
          std::cout << "All " << size << " markers will be analyzed " << std::endl;
*/
          return;

}



//==========================================================
//Destructor - free all dynamically allocated memory
//==========================================================
BGENFileReader::~BGENFileReader()	{	
	 free (m_filename);
	//free (nSample);
	//free (gm_sample_ids );
	//free( gm_stream);
	//free(gm_offset );
	//free(gm_context);
	//free(gm_have_sample_ids );
	//free(genoToTest_bgenDosage);
	//free(markerInfo);
}




/*


void BGENFileReader::BGEN_oneline(double* Z, int* myerror, int  size , char * chr, int * pos, char * snpid){
  
  using namespace genfile;
  //using namespace Rcpp;

  bool temp;
  std::string SNPID, RSID, chromosome, first_allele,second_allele ;
  genfile::bgen::uint32_t position ;
  std::vector< std::string > alleles ;
  std::vector< std::vector< double > > probs ;
  std::vector< byte_t > buffer1;
  std::vector< byte_t > buffer2;
  std::vector< double > dosages;

  //clock_t t1,t2;
  //t1=clock();
  temp = genfile::bgen::read_snp_identifying_data(
                        *this->gm_stream,
                        this->gm_context,
                        &SNPID,
                        &RSID,
                        &chromosome,
                        &position,
                        &first_allele,
                        &second_allele
                ) ;
  //t2=clock();
  //float diff = ((float)t2-(float)t1);
  //float seconds = diff / CLOCKS_PER_SEC;
  //std::cout << seconds << std::endl;



  //t1=clock();
  genfile::bgen::read_genotype_data_block(
                        *this->gm_stream,
                        this->gm_context,
                        &buffer1
                ) ;
  //t2=clock();
  //diff = ((float)t2-(float)t1);
  //seconds = diff / CLOCKS_PER_SEC;
  //std::cout << seconds << std::endl;



  //t1=clock();
  genfile::bgen::uncompress_probability_data(
                        this->gm_context,
                        buffer1,
                        &buffer2
                ) ;
  //t2=clock();	
  //diff = ((float)t2-(float)t1);
  //seconds = diff / CLOCKS_PER_SEC;
  //std::cout << seconds << std::endl;



  unsigned char * buf  = (unsigned char *) buffer2.data();
  uint Nbgen = this->gm_context.number_of_samples;

  //t1=clock();
  this->markerInfo = Parse(buf, buffer2.size(), SNPID, Nbgen, dosages);
  //t2=clock();
  //diff =  ((float)t2-(float)t1);
  //seconds = diff / CLOCKS_PER_SEC;
  //std::cout << seconds << std::endl;  

  chr= const_cast<char*>(chromosome.c_str());
  int temppos=  static_cast<int>(position);
  pos=&temppos;
  snpid=const_cast<char*>(RSID.c_str());
  Z=&dosages[0];

  dosages.clear();

  return;


}

*/




