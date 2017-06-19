 \name{GetGenotypesRegionVCF}
 \alias{GetGenotypesRegionVCF}
 \title{Select the genotype regions}
 \description{
	 Select the genotype regions wanted to read in the VCF file.  
 }
 \usage{
	GetGenotypesRegionVCF(VCF_info, chromosome=1, pos_start=NA, pos_end=NA)
 }
\arguments{
      \item{VCF_info}{VCF_info object returned from Open_VCF.}
      \item{chromosome}{The chromosome the user wants to read. }
      \item{pos_start}{The starting position the user wants to read in the specific chormosome.}
      \item{pos_end}{The ending position the user wants to read in the specific chormosome.}
}

\details{
	If pos_start and pos_end is missing, the function will select the whole chromosome. Then the user could call GetGenotypesVCF to get the genotype information line by line.
}


\author{Seunggeun Lee, Zhangchen Zhao}

