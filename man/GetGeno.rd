 \name{GetGeno}
 \alias{GetGeno}
 \title{Directly read the selected genotype regions}
 \description{
	 read selected genotype regions in the VCF file.  
 }
 \usage{
	GetGeno(VCF_info, chromosome=1, pos_start=NA, pos_end=NA)
 }
\arguments{
      \item{VCF_info}{VCF_info object returned from Open_VCF.}
      \item{chromosome}{The chromosome the user wants to read. }
      \item{pos_start}{The starting position the user wants to read in the specific chormosome.}
      \item{pos_end}{The ending position the user wants to read in the specific chormosome.}
}

\value{
 	\item{SNP Info}{The data frame contains the SNP information: SNPID, chromosome and position.}
	\item{genotype}{A genotype matrix with n rows and m columns, where n is the number of samples and m is the number of SNPs in the selected region.}
}

\details{
	If pos_start and pos_end is missing, the function will select the whole chromosome. 

	The function only read the SNP information that is included in 'Refgene' file. That is to say, the SNP information would be discarded if it doesn't belong to 'Refgene' file.
}


\author{Zhangchen Zhao}

