#define __STDC_LIMIT_MACROS
#include "int_info.h"
 #include <stdint.h> 
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "hts.h"
//#include "zlib.h"
//#include "zconf.h"

#include "HtslibFunc.h"


SEXP localfunc() {
    Rprintf("Hello from localfunc\n");
    return R_NilValue;
}

SEXP extvers() {
    return Rf_mkString(hts_version());
}

SEXP zlibVers() {
	return Rf_mkString(zlibVersion());
}


void Convert_BCFtoPlink(char **  PlinkPrefix,  char ** BCF_file, int * nmax, int * err) {
    
	int re = Convert_BCF_to_PlinkBed_Work(PlinkPrefix[0], BCF_file[0], nmax[0]);
	err[0] = re;
	
}

void Convert_BCFtoSSD( char **  ssd,char ** PlinkPrefix, char ** BCF_file, char ** SetID, int * SetIDtype, int * nmax, int * err) {
    
	int re = Convert_BCF_to_SSD_Work(ssd[0], PlinkPrefix[0],BCF_file[0], SetID[0], SetIDtype[0], nmax[0]);
	err[0] = re;
	
}




void Open_MWA(char* MWA_File, char* Info,int* myerror);
void Close_MWA() ;
void Get_Genotypes( int Set_number, int* Z, int size, int Is_MakeFile,int* myerror);
void Get_Genotypes_withID( int Set_number, int* Z, char * SNPID,int size, int Is_MakeFile, int* myerror);





void R_Open_MWA(char** MWA_File, char** Info, int * err)
{
	Open_MWA(MWA_File[0], Info[0],err);

}

void R_Close_MWA() 
{
	Close_MWA();
}





void R_Get_Genotypes( int *Set_number, int * Z , int * size, int *Is_MakeFile, int * err) // set_number base on INFO file. The result will be printed to file.
{
	Get_Genotypes( *Set_number, Z, *size,  * Is_MakeFile, err);
}	


void R_Get_Genotypes_withID( int *Set_number, int * Z , char * SNPID, int * size, int *Is_MakeFile, int * err) // set_number base on INFO file. The result will be printed to file.
{
	Get_Genotypes_withID( *Set_number, Z, SNPID,  *size,  * Is_MakeFile, err);
}	











