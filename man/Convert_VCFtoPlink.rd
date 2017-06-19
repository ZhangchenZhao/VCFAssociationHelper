 \name{Convert_VCFtoPlink }
 \alias{Convert_VCFtoPlink }
 \title{Generate Plink file}
 \description{
	Generate the Plink file from VCF file. 
 }
 \usage{
	Convert_VCFtoPlink( VCF_file, Plink_file, nmax=-1)
 
 }
\arguments{

	\item{VCF_file}{name of VCF, bzipped VCF or BCF file.}
	\item{Plink_file}{name of Plink file prefix. } 
	\item{nmax}{maximum snps that the function can read. If nmax=-1, the function will read all SNPs in the VCF file.}


}
\details{
 	This function will create Plink binary format files:  Bim file, Fam file and Bed file. 
}
\author{Zhangchen Zhao, Seunggeun Lee}




