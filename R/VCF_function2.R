lines.env <- new.env()

callLocal <- function() {
    .Call("localfunc")
}

htsVersion <- function() {
    .Call("extvers")
}





GetGenotypesVCF<-function(VCF_info){
######

	
############
	
	
	if(get("Helper_Lines_OPEN.isOpen", envir=lines.env) == 0){
		stop("VCF file is not opened. Please open it first!")
	}

##	id1<-which(SSD_INFO$SetInfo$SetIndex == Set_Index)
##	if(length(id1) == 0){
##		MSG<-sprintf("Error: cannot find set index [%d] from SSD!\n", Set_Index)
##		stop(MSG)
##	}	
##	Set_Index<-SSD_INFO$SetInfo$SetIndex[id1]
	size<-VCF_info$samplesize
	CHR_SIZE=100	
	chr=raw(CHR_SIZE)
	pos=0
	snpid="null"
	Z<-rep(9,size)
	err_code<-0
	SNP_ID_SIZE=1024
	snpid=raw(SNP_ID_SIZE)
	if (VCF_info$format=="GT"){
	temp<-.C("R_BCF_oneline",as.integer(Z), as.integer(err_code),as.integer(size),chr,as.integer(pos),snpid,PACKAGE="VCFAssociationHelper")
	}
	if (VCF_info$format=="DS"){
	temp<-.C("R_BCF_oneline1",as.double(Z), as.integer(err_code),as.integer(size),chr,as.integer(pos),snpid,PACKAGE="VCFAssociationHelper")
	}
	error_code<-temp[[2]]
		
	geno=list()
	##print(temp[[1]])	
	##geno$temp=	
	Z.out.t<-matrix(temp[[1]],byrow=TRUE, nrow=1)
	
	geno$SNPID=rawToChar(temp[[6]])
	geno$chromosome=rawToChar(temp[[4]])
	geno$position=temp[[5]]
	geno$genotypes=Z.out.t

	return(geno)

}

############
	




assign("Helper_Lines_OPEN.isOpen", 0, envir=lines.env)
assign("Helper_Lines_OPEN.FileName","", envir=lines.env)

CloseVCF<-function(){

	if(get("Helper_Lines_OPEN.isOpen", envir=lines.env) == 1){
		temp<-.C("R_Close_VCF", PACKAGE="VCFAssociationHelper")
		Msg<-sprintf("Close the opened VCF file: %s\n"
		,get("Helper_Lines_OPEN.FileName", envir=lines.env));
		cat(Msg)
		assign("Helper_Lines_OPEN.isOpen", 0, envir=lines.env);
	} else{
		Msg<-sprintf("No opened VCF file!\n");
		cat(Msg)		
	}
}

OpenVCF<-function(File.VCF,format="GT"){

	err_code<-0
	File.VCF<-normalizePath(File.VCF ,mustWork =FALSE)

	Check_File_Exists(File.VCF)

	if(get("Helper_Lines_OPEN.isOpen", envir=lines.env) == 1){
		CloseVCF();
	}
	
	size<-0
	if (format=="GT") {format_1=0} else {if (format=="DS"){format_1=1} else {stop("Wrong genotype format setting!")}}
	# Read VCF File
	temp<-.C("R_Open_VCF", as.character(File.VCF)
	, as.integer(err_code), as.integer(size),as.integer(format_1), PACKAGE="VCFAssociationHelper")

	error_code<-temp[[2]]
	Print_Error_SSD(error_code)

	if (format=="GT"){
	Msg<-sprintf("Open the VCF file with GT format\n");
	cat(Msg)} else { 	
		if (format=="DS"){
			Msg<-sprintf("Open the VCF file with DS format\n");
			cat(Msg)} } 
	

	#SSD_FILE_OPEN.isOpen<<-1
	#SSD_FILE_OPEN.FileName<<-File.SSD

	assign("Helper_Lines_OPEN.isOpen", 1, envir=lines.env)
	assign("Helper_Lines_OPEN.FileName",File.VCF, envir=lines.env)
	INFO=list()
	INFO$file=File.VCF
	INFO$format=format
	INFO$samplesize=temp[[3]]
	MSG<-sprintf("Sample size of VCF is \n", INFO$samplesize)
	cat(Msg)
	
	return(INFO)
	
}




GetGenotypesRegionVCF<-function(VCF_info,chromosome=1,pos_start=NA,pos_end=NA){
 	options("scipen"=100, "digits"=4)
	region=chromosome
	if (is.na(pos_start)) {
		print(sprintf("The start position is missing.\n"))
		if (is.na(pos_end)) {
			print(sprintf("The end position is also missing.The whole chromosome will be read.\n"))
			regions=region
		} else {
			print(sprintf("The program will read from the position 1.\n"))
			pos_start=1
			regions=paste(region,":",pos_start,"-",pos_end,sep="")
		}
	} else {
		if (is.na(pos_end)){
			stop("Please provide the end position!")} else {
		if (pos_end<pos_start){stop("Wrong position setting!")}
		regions=paste(region,":",pos_start,"-",pos_end,sep="")
		}
	}

	if(get("Helper_Lines_OPEN.isOpen", envir=lines.env) == 0){
		stop("VCF file is not opened. Please open it first!")
	}
	

	err_code<-0

	

	##temp1<-.C("R_SetRegion_search",as.character(regions), PACKAGE="VCFAssociationHelper")

	out = .Call("R_SetRegion", VCF_info$file, regions)
	return(out)



}
######




Get_Next_Genotypes<-function(){

	re = .Call("Get_Next_Genotypes")
	names(re)<-c("Genotype", "Pos")
	return(re)
}
