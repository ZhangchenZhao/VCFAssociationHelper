 \name{read.VCF.filter}
 \alias{read.VCF.filter}
 \title{Read VCF file with filter}
 \description{
	 Get the genotype information based on requirement of MAF and MAC.  
 }
 \usage{
	read.VCF.filter(VCF_file, format="GT",upperMAC=NA, lowerMAC=0, upperMAF=1,lowerMAF=0,nmax=1000) 
 }
\arguments{
      \item{VCF_file}{Name of VCF, bzipped VCF or BCF file.}
      \item{format}{There are two formats: "GT" and "DS". }
      \item{upperMAC}{The maximum MAC that the function will read.}
      \item{lowerMAC}{The minimum MAC that the function will read.}
      \item{upperMAC}{The maximum MAF that the function will read.}
      \item{lowerMAC}{The minimum MAF that the function will read.}
      \item{nmax}{The maximum genotypes that the function will read.}
}

\value{
 	\item{Genotype Information}{The genotype information: including SNPID, Chromosome, Position, etc.}
	\item{Genotype matrix}{The genotype information of the reading line.}
}

\author{Seunggeun Lee, Zhangchen Zhao}

