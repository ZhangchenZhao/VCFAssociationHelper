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
 	\item{SNPID}{The SNPID information from the reading line in VCF file.}
	\item{chromosome}{The chromosome that the reading line belongs to in VCF file.}
	\item{position}{The position of the reading line in the specific chromosome.}
	\item{genotype}{The genotype information of the reading line.}
}
\details{
	The user can use GetGenotypesRegionVCF to select the reading region. If the user keep call the function, it will read the line by line. If the user directly calls this fucnction, it will read from the first line of the VCF file.  
}

\author{Seunggeun Lee, Zhangchen Zhao}

