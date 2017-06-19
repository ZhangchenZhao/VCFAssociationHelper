/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   SNP_info.h
 * Author: zczhao
 *
 * Created on October 8, 2016, 9:04 AM
 */

#ifndef SNP_INFO_H
#define SNP_INFO_H



#endif /* SNP_INFO_H */

#include <fstream>  
#include <string>  
#include <iostream>  
#include <stdlib.h> 
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <map>
#include <iterator>

using namespace std; 


class SNP_info1 
{
    
public:
    std::string snp_id;
    std::string A1;
    std::string A2;
       std::string Chr;
   
       std::string  Pos;
    int dis;
    int flag;
    int total_counter_per_letter[2]; //= {0,0};
    int line_counter_per_letter[2]; //= {0,0};
	//int flag;
    
    SNP_info1(   std::string a,    std::string b,    std::string c,    std::string d,   std::string e, int f, int g){
        snp_id=a;
        A1=b;
        A2=c;
        Chr=d;
        Pos=e;
        dis=f;
        flag=g;
        total_counter_per_letter[0] = 0;
        total_counter_per_letter[1] = 0;
        line_counter_per_letter[0] = 0;
        line_counter_per_letter[1] = 0; 
 


    }


};



class refgene
{
    
public:
    std::string chr;
    size_t pos_start;
    size_t pos_end;
    std::string setid4;
    
    refgene(   std::string a,    size_t b,    size_t c, std::string d){
        chr=a;
        pos_start=b;
        pos_end=c;
	setid4=d;
    }


};

/*
class generate_SSD
{
public:
    string SetID; 
    string bim; 
    string bed;
    string bim_id;
    string ssd;
    int nsample;
          
    generate_SSD(char *a, char *b, char* c, char* d, char* e, int f){
        SetID= a;
        bim=b;
        bed=c;
        bim_id=d;
        ssd=e;
        nsample=f;  
     }

    //void gTokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters );     
    void Checkandcopy(char *SetID, char *bim, char *bim_id,int nsample, char *bed, char *ssd);
    
};
    

 */   
        
/*       return *snp;
    }        
    



};
*/
   
//void create_ssd(int nsample, char *bed, char *ssd)


    
