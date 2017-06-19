 \name{GetGenotypesSSD}
 \alias{GetGenotypesSSD}
 \title{Get Genotype data from SSD file}
 \description{
	Read a SSD file and return a genotype matrix.
 }
 \usage{
	GetGenotypesSSD(SSD_INFO, Set_Index, is_ID = FALSE)
 }
\arguments{
      \item{SSD_INFO}{SSD_INFO object returned from Open_SSD.}
      \item{Set_Index}{a numeric value of Set index. The set index of each set can be found from SetInfo object of SSD.INFO. }
      \item{is_ID}{a logical value indicating whether to read SNP ID (default=FALSE). If TRUE, it reads SNP IDs and use them as column names.}
}
\value{
 	A genotype matrix with n rows and m columns, where n is the number of samples and m is the number of SNPs.
}
\author{Seunggeun Lee, Larisa Miropolsky}

