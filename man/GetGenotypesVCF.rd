 \name{GetGenotypesVCF}
 \alias{GetGenotypesVCF}
 \title{Get Genotype data from VCF file}
 \description{
	Get one line's genotype information from the VCF file. 
 }
 \usage{
	GetGenotypesVCF(VCF_info)
 }
\arguments{
      \item{VCF_info}{VCF_info object returned from OpenVCF.}
 
}
\value{
 	\item{Genotype Information}{The genotype information: including SNPID, Chromosome, Position, etc.}
	\item{Genotype matrix}{The genotype information of the reading line.}
}
\details{
	The user can use GetGenotypesRegionVCF to select the reading region. If the user keep call the function, it will read the line by line. If the user directly calls this fucnction, it will read from the first line of the VCF file.  
}

\author{Seunggeun Lee, Zhangchen Zhao}

