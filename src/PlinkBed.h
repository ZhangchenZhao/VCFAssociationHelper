#ifndef __PlinkBed_H_
#define __PlinkBed_H_




#include <map>
#include <string>
#include <utility>
#include <iostream>
#include <fstream>
#include <vector>

class CPlinkBed_Write {

public:
    CPlinkBed_Write();
    ~CPlinkBed_Write();

    // Set file names and open Bed and Bim files for write
    int     Init(size_t nSampleSize=0);
    
    // Generate Bed, Bim and Fam file names
    int     Generate_FileName(const char * FilePrefix, bool AddTemp, bool AddLocalTime);
    
    // Write one SNP info
    int     Write_OneSNP(const char * chr, size_t pos, const char * A1, const char * A2, int * genotype);
    int     Write_OneSNP1(const char * snp1, const char * chr, size_t pos, const char * A1, const char * A2, int * genotype);
    int     Write_FamFile(char ** SampleID);
    int     Close_BedBim_Files();
    int     Delete_All_Files();
    int     Delete_BedBim_Files();
    
    std::string Get_BedFileName(){ return m_BedFileName;}
    std::string Get_BimFileName(){ return m_BimFileName;}
    std::string Get_FamFileName(){ return m_FamFileName;}
    
protected:
    
    std::string  m_BedFileName;
    std::string  m_BimFileName;
    std::string  m_FamFileName;

    size_t    m_nSNP;
    size_t    m_nSample;
    
     // Number of bytes that hold all genotype info per one snp = one line length!!!
    size_t    m_size_of_esi;
    
    std::ofstream m_BedFile;
    std::ofstream m_BimFile;
    
    char * m_buffer;
    char m_magic_num[3];
    
    bool Is_FileName;
    
};



void BedFile_encode(int* temp_snp_info,char* encoded_snp_info, size_t nsample, size_t size_of_esi);
void gTokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters );

#endif
