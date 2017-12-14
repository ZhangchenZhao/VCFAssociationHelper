 \name{read.VCF.SetID}
 \alias{read.VCF.SetID}
 \title{read VCF file based on SetID file}
 \description{
	 Read the genotype infomation wanted to read in the VCF file.  
 }
 \usage{
	read.VCF.SetID(VCF_file, format="GT",SetID_file=NULL, SetIDFormat) 
 }
\arguments{
      \item{VCF_file}{Name of VCF, bzipped VCF or BCF file.}
      \item{format}{There are two formats: "GT" and "DS". }
      \item{SetID_file}{Name of SetID file. There are three formats.}
      \item{SetIDFormat}{The format of SetID file. See details.}
}

\value{
 	\item{Genotype Information}{The genotype information: including SNPID, Chromosome, Position, etc.}
	\item{Genotype matrix}{The genotype information of the reading line.}
}


\details{
The SetID file has three formats: SNPID, POS and Region. SNPID contains two columns: SetID and SNP_ID. POS contains two columns. The first column is SetID and the second one is chromosome and position delimitered with “:”, such as 22:16050075. Region SetID file contains four columns: SetID, chromosome, starting position and ending position. SetID file is a white-space (space or tab) delimitered file.
 
Please keep in mind that there should be no header! The SNP_IDs and SetIDs should be less than 50 characters, otherwise, it will return an error message.
}


\author{Seunggeun Lee, Zhangchen Zhao}

