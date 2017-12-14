 \name{testyourfunction}
 \alias{testyourfunction}
 \title{Test your function in the specific region}
 \description{
	You can create your own single-variant-test function and this function can test it in the specific region.
 }
 \usage{
	testyourfunction(fun=NA,VCF_info=NA,Chr=NA,pos_start=NA,pos_end=NA,Y=NA,numcores=1)
 }
\arguments{
      \item{fun}{Your single-variant-test function.}
      \item{VCF_info}{VCF_info object returned from Open_VCF.}
      \item{Chr}{The chromosome the user wants to test. }
      \item{pos_start}{The starting position the user wants to test in the specific chormosome.}
      \item{pos_end}{The ending position the user wants to test in the specific chormosome.}
      \item{Y}{The phenotype information of the individuals in the VCF file.}
      \item{numcores}{The number of cores that the user wants to use.}
}

\value{
 	\item{Matrix Information}{The matrix contains four columns: SNPID, Chromosome, Position, Test Value.}
}

\author{Seunggeun Lee, Zhangchen Zhao}

