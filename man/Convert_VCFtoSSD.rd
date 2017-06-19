 \name{Convert_VCFtoSSD }
 \alias{Convert_VCFtoSSD }
 \title{Generate SNP set data file (SSD)}
 \description{
	Generate a SNP set data file (SSD) from VCF files using user specified SNP sets.
 }
 \usage{

	Convert_VCFtoSSD(VCF_file, SSD_file, SetID_file=NULL, SetIDFormat, nmax=-1)
   
 }
\arguments{
      \item{VCF_file}{name of VCF, bzipped VCF or BCF file.}
      \item{SSD_file}{name of SSD file to be created.}
      \item{SetID_file}{name of the Set ID file. If NULL, the function will use refGene (GRCh37) to identify genes.}
      \item{SetIDFormat}{the format of SetID file. See details.}
      \item{nmax}{maximum snps that the function can read. If nmax=-1, the function will read all SNPs in the VCF file.}
}
\details{
 	The SetID file has three formats: SNPID, POS and Region. SNPID contains two columns: SetID and SNP_ID. POS contains two columns. The first column is SetID and the second one is chromosome and position delimitered with “:”, such as 22:16050075. Region SetID file contains four columns: SetID, chromosome, starting position and ending position. SetID file is a white-space (space or tab) delimitered file.
 
Please keep in mind that there should be no header! The SNP_IDs and SetIDs should be less than 50 characters, otherwise, it will return an error message.

The SSD file is a binary formatted file with genotypes. The SSD info file is a text file with general information on data and SNP sets (first 6 rows), and information on each set (after 8th row). In addition, SSDinfo file is generated with SSD file.   

}

\author{Zhangchen Zhao, Seunggeun Lee}





