### R code from vignette source 'VCFAssociationHelper.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Plink
###################################################
library(VCFAssociationHelper)
VCF_file =paste(path.package("VCFAssociationHelper"),"/extdata/Example.vcf.gz",sep="")  
PlinkPrefix = paste(path.package("VCFAssociationHelper"),"/extdata/Example.plink",sep="")  
nmax=1000
system.time({Convert_VCFtoPlink ( VCF_file, PlinkPrefix, nmax)})


###################################################
### code chunk number 2: SSD
###################################################
# SetID with SNPID
VCF_file =paste(path.package("VCFAssociationHelper"),"/extdata/Example.vcf.gz",sep="")  
SetID = paste(path.package("VCFAssociationHelper"),"/extdata/SetID1.txt",sep="")  
SSD=paste(path.package("VCFAssociationHelper"),"/extdata/SSD1",sep="") 
SetIDFormat="SNPID"
nmax=100000
system.time({Convert_VCFtoSSD ( VCF_file,  format="GT",SSD, SetID, SetIDFormat, nmax)})

# SetID with POS
VCF_file =paste(path.package("VCFAssociationHelper"),"/extdata/Example.vcf.gz",sep="")  
SetID = paste(path.package("VCFAssociationHelper"),"/extdata/SetID2.txt",sep="")  
SSD=paste(path.package("VCFAssociationHelper"),"/extdata/SSD2",sep="") 
SetIDFormat="POS"
nmax=100000
system.time({Convert_VCFtoSSD ( VCF_file, format="GT", SSD, SetID, SetIDFormat, nmax)})

# SetID with SNPID
VCF_file =paste(path.package("VCFAssociationHelper"),"/extdata/Example.vcf.gz",sep="")  
SetID = paste(path.package("VCFAssociationHelper"),"/extdata/SetID3.txt",sep="")  
SSD=paste(path.package("VCFAssociationHelper"),"/extdata/SSD3",sep="") 
SetIDFormat="Region"
nmax=100000
system.time({Convert_VCFtoSSD ( VCF_file,  format="GT",SSD, SetID, SetIDFormat, nmax)})


###################################################	
### code chunk number 3: geno
###################################################
#Open SSD file first.
SSD=paste(path.package("VCFAssociationHelper"),"/extdata/SSD3",sep="") 
SSDinfo=paste(SSD,"info",sep="") 
ssdtemp=OpenSSD(SSD, SSDinfo)
Set_Index=1

#Get the genotypes from SSD.
geno=GetGenotypesSSD(ssdtemp, Set_Index, is_ID = FALSE)
dim(geno)
head(geno)
#SSD must be closed after finishing the test!
CloseSSD()


###################################################
### code chunk number 4: VCFread
###################################################
VCF_file =paste(path.package("VCFAssociationHelper"),"/extdata/Example.vcf.gz",sep="")  

#Open VCF file with GT format.
vcf_info=OpenVCF(VCF_file,"GT")

#Get the genotype information of the first line in VCF file.
geno1=GetGenotypesVCF( vcf_info)
geno1[[1]]
sum(geno1[[2]])
#Get the genotype information of the second line in VCF file.
geno2=GetGenotypesVCF( vcf_info)
geno2[[1]]
sum(geno2[[2]])
CloseVCF()


###################################################
### code chunk number 5: VCFregion_read
###################################################
VCF_file =paste(path.package("VCFAssociationHelper"),"/extdata/Example.vcf.gz",sep="")  
vcf_info=OpenVCF(VCF_file,"GT")

#Select the region: chromosome 22, position 16000000~17000000.
#Method 1: read genotype information one line by one line
GetGenotypesRegionVCF( vcf_info,22,16000000,17000000)
geno1=GetGenotypesVCF( vcf_info)
geno1[[1]]
sum(geno1[[2]])
geno2=GetGenotypesVCF( vcf_info)
geno2[[1]]
sum(geno2[[2]])

#Method 2: read all selected genotype information
geno_all=GetGeno( vcf_info,22,16000000,17000000)
head(geno_all[[1]])
CloseVCF()


