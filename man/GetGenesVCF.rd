 \name{GetGenesVCF}
 \alias{GetGenesVCF}
 \title{Get Genotype data of a typical gene}
 \description{
	Get the genotype information of a typical gene from VCF file
 }
 \usage{
	GetGenesVCF(VCF_info,Gene_name)
 }
\arguments{
      \item{VCF_info}{VCF_info object returned from OpenVCF.}
      \item{Gene_name}{the name of the gene}
 
}
\value{
 	\item{Genotype Information}{The genotype information: including SNPID, Chromosome, Position, etc.}
	\item{Genotype matrix}{The genotype information of the gene}
}
\details{
	We use Refgene file in default to find the position of the gene, further find the genotype information of this gene.
}

\author{Seunggeun Lee, Zhangchen Zhao}

