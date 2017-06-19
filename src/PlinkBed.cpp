#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <math.h>
#include <ctime>
#include <cstring>

#include "PlinkBed.h"

CPlinkBed_Write::CPlinkBed_Write(){
    
    Is_FileName=false;
    m_buffer = NULL;
    
    // magic number 01101100 00011011 00000001
    m_magic_num[0]= 0x6C;
    m_magic_num[1]= 0x1B;
    m_magic_num[2]= 0x1;
    
    
    
}

CPlinkBed_Write::~CPlinkBed_Write(){
    
    if(m_buffer){
        delete [] m_buffer;
    }
    
}

int     CPlinkBed_Write::Init(size_t nSampleSize){
    
    if(!Is_FileName){
        return 0;
    }

    /* Open Bed and Bim files */
    m_BedFile.open(m_BedFileName.c_str(), std::ofstream::out | std::ofstream::binary);
	if (!m_BedFile)
	{
        printf("Cannot open %s\n", m_BedFileName.c_str());
		return 0;
	}
    // write magic number
    m_BedFile.write(m_magic_num, 3 * sizeof(char));
    
    m_BimFile.open(m_BimFileName.c_str(), std::ofstream::out );
	if (!m_BimFile)
	{
        printf("Cannot open %s\n", m_BimFileName.c_str());
		return 0;
	}
    
    // Set other values
    m_nSNP=0;
    m_nSample = nSampleSize;
    m_size_of_esi = (nSampleSize+3)/4;
    m_buffer = new char[m_size_of_esi*10 +1];
    
    
    return 1;
}

int     CPlinkBed_Write::Generate_FileName(const char * FilePrefix, bool AddTemp, bool AddLocalTime){
    
 
    time_t rawtime;
    struct tm * timeinfo;
    
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    char time_buff[1000];
    
    strftime (time_buff,80,"%y_%j_%H_%M_%S",timeinfo);
    
    
    m_BedFileName = FilePrefix;
    m_BimFileName = FilePrefix;
    m_FamFileName = FilePrefix;
    
    if(AddLocalTime==true){
        m_BedFileName += time_buff;
        m_BimFileName += time_buff;
        //m_FamFileName += time_buff;
    }
    
    if(AddTemp==true){
        m_BedFileName += ".temp";
        m_BimFileName += ".temp";
        //m_FamFileName += ".temp";
    } 
    
    m_BedFileName += ".bed";
    m_BimFileName += ".bim";
    m_FamFileName += ".fam";        
    
    
    Is_FileName=true;
    return 1;
}
    

int     CPlinkBed_Write::Write_OneSNP(const char * chr, size_t pos, const char * A1, const char * A2, int * genotype){
    
    char str[10000];
    
    // Bim file
    sprintf(str,"%s %s:%zu 0 %zu %s %s\n", chr, chr, pos, pos, A1, A2);
    m_BimFile << str;

    // Bed file
    BedFile_encode(genotype, m_buffer, m_nSample, m_size_of_esi);
    m_BedFile.write(m_buffer, m_size_of_esi* sizeof(char));

    
   	if (!m_BedFile)
	{
		return 0;
	}
    
    
    m_nSNP++;
    return 1;
    
    
}


int     CPlinkBed_Write::Write_OneSNP1(const char * snp1, const char * chr, size_t pos, const char * A1, const char * A2, int * genotype){
    
    char str[10000];
    
    // Bim file
    sprintf(str,"%s %s 0 %zu %s %s\n", chr,snp1, pos, A1, A2);
    m_BimFile << str;

    // Bed file
    BedFile_encode(genotype, m_buffer, m_nSample, m_size_of_esi);
    m_BedFile.write(m_buffer, m_size_of_esi* sizeof(char));

    
   	if (!m_BedFile)
	{
		return 0;
	}
    
    
    m_nSNP++;
    return 1;
    
    
}




int     CPlinkBed_Write::Write_FamFile(char ** SampleID){
    
    /* Open Bed and Bim files */
    std::ofstream m_FamFile;
    m_FamFile.open(m_FamFileName.c_str(), std::ofstream::out);
	if (!m_FamFile)
	{
		return 0;
	}

    char str[10000];
    
    for(int i=0; i<m_nSample;i++){
        
        // FAM file
        sprintf(str,"%s\t%s\t0\t0\t0\t0\n", SampleID[i], SampleID[i]);
        m_FamFile << str;
    }
    m_FamFile.close();
	if (!m_FamFile)
	{
		return 0;
	}
    
    return 1;
}

int     CPlinkBed_Write::Close_BedBim_Files(){
    
    m_BimFile.close();
	if (!m_BimFile)
	{
		return 0;
	}
    m_BedFile.close();
	if (!m_BedFile)
	{
		return 0;
	}
    
    return 1;
}

int     CPlinkBed_Write::Delete_All_Files(){
    
    std::remove(m_BedFileName.c_str());
    std::remove(m_BimFileName.c_str());
    std::remove(m_FamFileName.c_str());
    
    return 1;
}


int     CPlinkBed_Write::Delete_BedBim_Files(){
    
    std::remove(m_BedFileName.c_str());
    std::remove(m_BimFileName.c_str());
    
    return 1;
}





//==========================================================
// This function converting(encoding) every 4 bytes of temp_snp_info to "a"
// "a"- array of 8 integers ("kind of" 8 bits) that will be converted to one encoded number
// Compataible to plink bed file encode
// This encode is different from encode in bed_reader.cpp
//==========================================================
void BedFile_encode(int* temp_snp_info,char* encoded_snp_info, size_t nsample, size_t size_of_esi)
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
	while (i < nsample)
	{
		memset(a, 0, sizeof(a));
		for (j = 0; j < 4; ++j)
		{
			if (temp_snp_info[i] == 9) //missing value
			{
                // different from bed_reader.cpp encode
				a[j*2] = 0;
				a[j*2+1] = 1;
			}
			else if (temp_snp_info[i] == 1)//AG
			{
                // different from bed_reader.cpp encode
				a[j*2] = 1;
				a[j*2+1] = 0;
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
		if (ind4enc == size_of_esi+1)
			break;
		else
		{
			number = 0;
			//==============================================
			//converting 8 ints of "a" to one encoded number
			//for example a= 00100010 "number" will be: 34
			for (int  ii = 0; ii < 8; ++ii)
				number += a[ii] * (int)pow(2.0,(7-ii));
			//saving this encoded number.
			encoded_snp_info[ind4enc-1] = (char)number;
			//=============================================
		}
	}
}

//=======================================================================
// This function split "str" by "delimiters" and put the result to "tokens"
// Inputs:
// str - source string
// tokens - target vector of strings
// delimeters - string separator
// Notes:
// "tokens" - vector of strings, so to move inside use:
// tokens.at(0).c_str(); tokens.at(1).c_str().
// important clear "tokens" after each iteration: tokens.clear();
// otherwise it will append new results to existing.
//=======================================================================

void gTokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters ){
    
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
    
    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}




