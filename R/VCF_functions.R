SSD1.env <- new.env()

callLocal <- function() {
    .Call("localfunc")
}

htsVersion <- function() {
    .Call("extvers")
}



Print_Error_SSD<-function(code){

	if(code == 0 ){
		return(0);
	} else if(code == 1){
		stop("Error Can't open BIM file")
	} else if(code == 2){
		stop("Error Can't open FAM file")
	} else if(code == 3){
		stop("Error Can't open BED file")
	} else if(code == 4){
		stop("Error Can't open SETID file")
	} else if(code == 5){
		stop("Error Can't write SSD file")
	} else if(code == 6){
		stop("Error Can't read SSD file")
	} else if(code == 7){
		stop("Error Can't write INFO file")
	} else if(code == 8){
		stop("Error Can't read INFO file")
	} else if(code == 9){
		stop("Error Can't write INFO file")
	} else if(code == 13){
		stop("Error Wrong SNP or Individual sizes")
	} else if(code == 14){
		stop("Error SetID not found")
	} else {
		MSG<-sprintf("Error [%d]\n",code)
		stop(MSG)
	}
	
	return(1)
}



Check_File_Exists<-function(FileName){
	
	if(!file.exists(FileName)){
		Msg<-sprintf("File %s does not exist\n",FileName)
		stop(Msg)
	}

}

Print_File_Info<-function(INFO){
	
	MSG<-sprintf("%d Samples, %d Sets, %d Total SNPs\n",INFO$nSample, INFO$nSets, INFO$nSNPs)
	cat(MSG)

}

Read_File_Info<-function(File.Info){

	Check_File_Exists(File.Info)
	
	info1<-read.table(File.Info, header=FALSE, nrows= 7, sep='\t')
	info2<-read.table(File.Info, header=FALSE, skip=8, sep='\t', stringsAsFactors=FALSE,comment.char = "")

	INFO<-list()
	INFO$WindowSize<-as.numeric(as.character(info1[1,1]))
	INFO$MAFConvert<-as.numeric(as.character(info1[2,1]))
	INFO$nSNPs<-as.numeric(as.character(info1[3,1]))	
	INFO$nSample<-as.numeric(as.character(info1[4,1]))
	INFO$format<-as.character(info1[5,1])
	INFO$nDecodeSize<-as.numeric(as.character(info1[6,1]))
	INFO$nSets<-as.numeric(as.character(info1[7,1]))
	INFO$SetInfo<-data.frame(SetIndex= info2[,1], SetID = as.character(info2[,3])
	, SetSize = info2[,4], Offset = info2[,2],stringsAsFactors=FALSE)

	Print_File_Info(INFO)
	
	return(INFO)


}









# SEXP Convert_BCFtoPlink(const char * BCF_file, const char *PlinkPrefix,  const char *SetID, nmax) 


# SEXP Convert_BCFtoPlink(const char * BCF_file, const char *PlinkPrefix, nmax) 
Convert_VCFtoPlink <- function( VCF_file, Plink_file, nmax=-1) {
    
    err_code<-0
	VCF_file<-normalizePath(VCF_file ,mustWork =FALSE)
	Plink_file<-normalizePath(Plink_file ,mustWork =FALSE)

	Check_File_Exists(VCF_file)


    temp<-.C("Convert_BCFtoPlink", as.character(Plink_file), as.character(VCF_file), as.integer(nmax), as.integer(err_code))
	error_code<-temp[[4]]
	cat(error_code)
	
}



Convert_VCFtoSSD <- function(VCF_file, format="GT",SSD_file,  SetID_file=NULL, SetIDFormat, nmax=-1) {
    if (format!="GT"){
	if (format!="DS"){
		Msg<-sprintf("VCF_file Reading Format is wrong!\n")
		stop(Msg)
	}
    }
    if (format=="GT"){formatcode=0}
    if (format=="DS"){formatcode=1}
   
    
    if (is.null(SetID_file)) {
    	SetID_file=paste(path.package("VCFAssociationHelper"),"/extdata/refGene.txt",sep="")
		print(sprintf("The SetID is missing. refGene (GRCh37) will be used to define genes.\n"))
		print(sprintf("refGene file can be found at: %s\n", SetID_file))
		SetIDFormat="Ref"
	}
	if (is.na(SetIDFormat)) {
		Msg<-sprintf("SetID Format is missing!\n")
		stop(Msg)
	}
	err_code<-0
	VCF_file<-normalizePath(VCF_file ,mustWork =FALSE)
	
##PlinkPrefix<-SSD_file
	SetID_file<-normalizePath(SetID_file ,mustWork =FALSE)
	SSD_file<-normalizePath(SSD_file ,mustWork =FALSE)
	PlinkPrefix<-paste(SSD_file,"plink",sep=".")
	ssdinfo<-paste(SSD_file,"info",sep=".")
	Check_File_Exists(VCF_file)	
	Check_File_Exists(SetID_file)
	SetIDtype<-0
	if (SetIDFormat=="SNPID") { SetIDtype=1} else {
		if (SetIDFormat=="POS") { SetIDtype=2} else {
			if (SetIDFormat=="Region") { SetIDtype=3} else {
				if (SetIDFormat=="Ref"){SetIDtype=4} else {
					Msg<-sprintf("SetID Format isn't correct!\n")
					stop(Msg)
				}
			}	 
		}
	}

		
	print(PlinkPrefix)
    temp<-.C("Convert_BCFtoSSD", as.character(SSD_file), as.character(PlinkPrefix),as.character(VCF_file),as.character(SetID_file), as.integer(SetIDtype),as.integer(nmax), as.integer(err_code),as.integer(formatcode))
	error_code<-temp[[7]]
	cat(error_code)
	


	
}





Convert_PlinktoSSD <- function(File.bed,File.bim,File.fam, SSD_file, nmax=-1) {
#    if (format!="GT"){
#	if (format!="DS"){
#		Msg<-sprintf("VCF_file Reading Format is wrong!\n")
#		stop(Msg)
#	}
#   }
#    if (format=="GT"){formatcode=0}
#    if (format=="DS"){formatcode=1}
   
    
#    if (is.null(SetID_file)) {
#    	SetID_file=paste(path.package("VCFAssociationHelper"),"/extdata/refGene.txt",sep="")
#		print(sprintf("The SetID is missing. refGene (GRCh37) will be used to define genes.\n"))
#		print(sprintf("refGene file can be found at: %s\n", SetID_file))
#		SetIDFormat="Ref"
#	}
#	if (is.na(SetIDFormat)) {
#		Msg<-sprintf("SetID Format is missing!\n")
#		stop(Msg)
#	}
	err_code<-0
	File.bed<-normalizePath(File.bed ,mustWork =FALSE)
	File.bim<-normalizePath(File.bim ,mustWork =FALSE)
	File.fam<-normalizePath(File.fam ,mustWork =FALSE)



##PlinkPrefix<-SSD_file
#	SetID_file<-normalizePath(SetID_file ,mustWork =FALSE)
	SSD_file<-normalizePath(SSD_file ,mustWork =FALSE)
	PlinkPrefix<-paste(SSD_file,"plink",sep=".")
	ssdinfo<-paste(SSD_file,"info",sep=".")
	Check_File_Exists(File.bed)	
	Check_File_Exists(File.bim)
	Check_File_Exists(File.fam)

	print(PlinkPrefix)
	SetID_file=paste(path.package("VCFAssociationHelper"),"/extdata/refGene.txt",sep="")

    temp<-.C("Convert_PlinktoSSD", as.character(SSD_file), as.character(PlinkPrefix),as.character(File.bed),as.character(File.bim), 
as.character(File.fam), as.character(SetID_file), as.integer(nmax), as.integer(err_code))
	error_code<-temp[[7]]
	cat(error_code)
	


	
}



##################################################################
#
#	Open and Close the SSD Files


assign("Helper_SSD_FILE_OPEN.isOpen", 0, envir=SSD1.env)
assign("Helper_SSD_FILE_OPEN.FileName","", envir=SSD1.env)

CloseSSD<-function(){

	if(get("Helper_SSD_FILE_OPEN.isOpen", envir=SSD1.env) == 1){
		temp<-.C("R_Close_MWA", PACKAGE="VCFAssociationHelper")
		Msg<-sprintf("Close the opened SSD file: %s\n"
		,get("Helper_SSD_FILE_OPEN.FileName", envir=SSD1.env));
		cat(Msg)
		assign("Helper_SSD_FILE_OPEN.isOpen", 0, envir=SSD1.env);
	} else{
		Msg<-sprintf("No opened SSD file!\n");
		cat(Msg)		
	}
}

OpenSSD<-function(File.SSD, File.Info){

	err_code<-0
	File.SSD<-normalizePath(File.SSD ,mustWork =FALSE)
	File.Info<-normalizePath(File.Info ,mustWork =FALSE)

	Check_File_Exists(File.SSD)
	Check_File_Exists(File.Info)

	if(get("Helper_SSD_FILE_OPEN.isOpen", envir=SSD1.env) == 1){
		CloseSSD();
	}

	# Read Info File
	INFO<-Read_File_Info(File.Info)
	##format_code=0
	##if (INFO$format=="DS"){ format_code=1}


	# Read SSD File
	temp<-.C("R_Open_MWA", as.character(File.SSD), as.character(File.Info)
	, as.integer(err_code), PACKAGE="VCFAssociationHelper")

	error_code<-temp[[3]]
	Print_Error_SSD(error_code)


	Msg<-sprintf("Open the SSD file\n");
	cat(Msg)

	#SSD_FILE_OPEN.isOpen<<-1
	#SSD_FILE_OPEN.FileName<<-File.SSD

	assign("Helper_SSD_FILE_OPEN.isOpen", 1, envir=SSD1.env)
	assign("Helper_SSD_FILE_OPEN.FileName",File.SSD, envir=SSD1.env)


	return(INFO)
	
}

#######################################################


##################################################################
#
#	Get Genotype Matrix

GetGenotypesSSD<-function(SSD_INFO, Set_Index, is_ID = FALSE){

	SNP_ID_SIZE=1024 # it should be the same as SNP_ID_SIZE_MAX in error_messages.h 
	
	Is_MakeFile=0
	if(get("Helper_SSD_FILE_OPEN.isOpen", envir=SSD1.env) == 0){
		stop("SSD file is not opened. Please open it first!")
	}

	id1<-which(SSD_INFO$SetInfo$SetIndex == Set_Index)
	if(length(id1) == 0){
		MSG<-sprintf("Error: cannot find set index [%d] from SSD!\n", Set_Index)
		stop(MSG)
	}	
	Set_Index<-SSD_INFO$SetInfo$SetIndex[id1]

	err_code<-0
	N.SNP<-SSD_INFO$SetInfo$SetSize[id1]
	N.Sample<-SSD_INFO$nSample
	size<-N.SNP * N.Sample

	Z.out.t=NULL
	if (SSD_INFO$format=="GT"){
	Z<-rep(9,size)

	if(!is_ID){
		temp<-.C("R_Get_Genotypes",as.integer(Set_Index),as.integer(Z),as.integer(size)
		,as.integer(Is_MakeFile), as.integer(err_code), PACKAGE="VCFAssociationHelper")

		error_code<-temp[[5]]
		Print_Error_SSD(error_code)
		
		Z.out.t<-matrix(temp[[2]],byrow=TRUE, nrow=N.SNP)
		
	} else {
		SNPID=raw(N.SNP* SNP_ID_SIZE)
	
		temp<-.C("R_Get_Genotypes_withID",as.integer(Set_Index),as.integer(Z), SNPID, as.integer(size)
		,as.integer(Is_MakeFile), as.integer(err_code), PACKAGE="VCFAssociationHelper")

		error_code<-temp[[6]]
		Print_Error_SSD(error_code)
		
		SNPID.m<-matrix(temp[[3]], byrow=TRUE, nrow=N.SNP)
		SNPID.c<-apply(SNPID.m, 1, rawToChar)

		Z.out.t<-matrix(temp[[2]],byrow=TRUE, nrow=N.SNP)
		rownames(Z.out.t)<-SNPID.c
		
		#SNPID.c1<<-SNPID.c
		#SNPID.1<<-SNPID	
	}
	}


	if (SSD_INFO$format=="DS"){
	Z<-rep(9,size)

	if(!is_ID){
		temp<-.C("R_Get_Genotypes_ds",as.integer(Set_Index),as.double(Z),as.integer(size)
		,as.integer(Is_MakeFile), as.integer(err_code), PACKAGE="VCFAssociationHelper")

		error_code<-temp[[5]]
		Print_Error_SSD(error_code)
		
		Z.out.t<-matrix(temp[[2]],byrow=TRUE, nrow=N.SNP)
		
	} else {
		SNPID=raw(N.SNP* SNP_ID_SIZE)
	
		temp<-.C("R_Get_Genotypes_withID_ds",as.integer(Set_Index),as.double(Z), SNPID, as.integer(size)
		,as.integer(Is_MakeFile), as.integer(err_code), PACKAGE="VCFAssociationHelper")

		error_code<-temp[[6]]
		Print_Error_SSD(error_code)
		
		SNPID.m<-matrix(temp[[3]], byrow=TRUE, nrow=N.SNP)
		SNPID.c<-apply(SNPID.m, 1, rawToChar)

		Z.out.t<-matrix(temp[[2]],byrow=TRUE, nrow=N.SNP)
		rownames(Z.out.t)<-SNPID.c
		
		#SNPID.c1<<-SNPID.c
		#SNPID.1<<-SNPID	
	}
	}
	
	return(t(Z.out.t))
}



#######################################################


##################################################################
#
#	Get Genotype Matrix

GetGenotypesSSD_New<-function(SSD_INFO, Set_Index){

	SNP_ID_SIZE=1024 # it should be the same as SNP_ID_SIZE_MAX in error_messages.h 
	
	Is_MakeFile=0
	if(get("Helper_SSD_FILE_OPEN.isOpen", envir=SSD1.env) == 0){
		stop("SSD file is not opened. Please open it first!")
	}

	id1<-which(SSD_INFO$SetInfo$SetIndex == Set_Index)
	if(length(id1) == 0){
		MSG<-sprintf("Error: cannot find set index [%d] from SSD!\n", Set_Index)
		stop(MSG)
	}	
	Set_Index<-SSD_INFO$SetInfo$SetIndex[id1]

	err_code<-0
	N.SNP_total<-SSD_INFO$SetInfo$SetSize[id1]
	N.Sample<-SSD_INFO$nSample
	
	N.SNP_left=N.SNP_total 
	Z.out.t=NULL
	Z.out=NULL
	Pos=SSD_INFO$SetInfo$Offset[id1]
	
	if (SSD_INFO$format=="GT"){
		flag=FALSE
		i=1
		while (flag==FALSE){
			if (N.SNP_left>10 ){
				N.SNP=10
				N.SNP_left=N.SNP_left-10
			} else {
				flag=TRUE
				N.SNP=N.SNP_left	
			}
			size<-N.SNP * N.Sample
			Z<-rep(9,size)

		
			SNPID=raw(N.SNP* SNP_ID_SIZE)
			
		

			temp<-.C("R_Get_Genotypes_withID_new",as.integer(Set_Index),as.integer(Z), SNPID, as.integer(size)
			,as.integer(Is_MakeFile), as.integer(err_code), as.integer(Pos),as.integer(N.SNP),PACKAGE="VCFAssociationHelper")

	
			
			error_code<-temp[[6]]
			Print_Error_SSD(error_code)
		
			SNPID.m<-matrix(temp[[3]], byrow=TRUE, nrow=N.SNP)
			SNPID.c<-apply(SNPID.m, 1, rawToChar)

			Z.out.t<-Matrix(temp[[2]],byrow=TRUE, nrow=N.SNP,sparse=TRUE)
			rownames(Z.out.t)<-SNPID.c
			

			Pos=temp[[7]]-1
			rm(temp)
			gc()
			if (i==1){Z.out<-Z.out.t} else {Z.out=Matrix(rbind(Z.out,Z.out.t), sparse=TRUE)}
			i=i+1;		
			#Z.out=Matrix(rbind(Z.out,Z.out.t), sparse=TRUE)

		}	
	
	}

	return(Z.out)
}




	

