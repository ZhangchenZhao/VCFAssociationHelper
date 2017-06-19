 \name{OpenVCF}
 \alias{OpenVCF}
 \title{Open VCF file}
 \description{
 Open the VCF file. After finishing using the VCF file, you must close the file by callinsg CloseVCF function. 
 }
 \usage{
	OpenVCF(File.VCF, format="GT")
 }
\arguments{
      \item{File.VCF}{ name of the VCF file . }
      \item{format}{ the genotype format to read. See details. }
}
\value{
 	 a list object of VCF_info.
}
\details{
 	We can read two genotype formats in VCF file: GT and SD. The default is GT. 

}

\author{Seunggeun Lee, Zhangchen Zhao}

