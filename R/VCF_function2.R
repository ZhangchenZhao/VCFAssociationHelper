callLocal <- function() {
    .Call("localfunc")
}

htsVersion <- function() {
    .Call("extvers")
}


lines.env <- new.env()


assign("Helper_Lines_OPEN.isOpen", 0, envir=lines.env)
assign("Helper_Lines_OPEN.FileName","", envir=lines.env)


	
CloseVCF<-function(){

	if(get("Helper_Lines_OPEN.isOpen", envir=lines.env) == 1){
		temp<-.C("R_Close_VCF", PACKAGE="VCFAssociationHelper")
		Msg<-sprintf("Close the opened VCF file: %s\n"
		,get("Helper_Lines_OPEN.FileName", envir=lines.env));
		cat(Msg)
		assign("Helper_Lines_OPEN.isOpen", 0, envir=lines.env);
		assign("Helper_Lines_OPEN.FileName","", envir=lines.env)
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
	sample=c("null")
	temp<-.C("R_Open_VCF", as.character(File.VCF)
	, as.integer(err_code), as.integer(size),as.integer(format_1),PACKAGE="VCFAssociationHelper")

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

	SIZE=1024 
	samples=raw(INFO$samplesize*SIZE)
	samplesid<-.C("R_getid",samples,PACKAGE="VCFAssociationHelper")
	
	ID.m<-matrix(samplesid[[1]], byrow=TRUE, nrow=INFO$samplesize)
	ID.c<-apply(ID.m, 1, rawToChar)

	INFO$samples=ID.c
	MSG<-sprintf("Sample size of VCF is %i.\n", INFO$samplesize)
	cat(MSG)
	
	return(INFO)
	
}


GetGenotypesVCF<-function(VCF_info){
######

	
############
	#print("R1")
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
	A1=raw(CHR_SIZE)
	A2=raw(CHR_SIZE)
	#print("R2")
	if (VCF_info$format=="GT"){
	temp<-.C("R_BCF_oneline",as.integer(Z), as.integer(err_code),as.integer(size),chr,as.integer(pos),snpid,A1,A2,PACKAGE="VCFAssociationHelper")
	}
	if (VCF_info$format=="DS"){
	temp<-.C("R_BCF_oneline1",as.double(Z), as.integer(err_code),as.integer(size),chr,as.integer(pos),snpid,A1,A2,PACKAGE="VCFAssociationHelper")
	}
	##print(temp[[5]])
	error_code<-temp[[2]]
	##print(temp[[1]][1])	
	geno=list()
	##print(temp[[1]])	
	##geno$temp=	
	Z.out.t<-matrix(temp[[1]],byrow=TRUE, nrow=1)
	##print(temp[[1]][1])

	##print(temp[[6]])
	##print(temp[[4]])
	SNPID=rawToChar(temp[[6]])
	chromosome=rawToChar(temp[[4]])
	##print(SNPID)
	##print(chromosome)
	position=temp[[5]]
	allele1=rawToChar(temp[[7]])
	allele2=rawToChar(temp[[8]])
	##print(position)
	##print(allele1)
	##print(allele2)
	geno$info=cbind(SNPID,chromosome,position,allele1,allele2)
	geno$genotypes=Z.out.t
	##print(geno$info)
	return(geno)

}

############



GetGenesVCF<-function(VCF_info,Gene_name){
############
	#print("R1")
	if(get("Helper_Lines_OPEN.isOpen", envir=lines.env) == 0){
		stop("VCF file is not opened. Please open it first!")
	}

##	id1<-which(SSD_INFO$SetInfo$SetIndex == Set_Index)
##	if(length(id1) == 0){
##		MSG<-sprintf("Error: cannot find set index [%d] from SSD!\n", Set_Index)
##		stop(MSG)
##	}	
##	Set_Index<-SSD_INFO$SetInfo$SetIndex[id1]
	err_code<-0

	SetID_file=paste(path.package("VCFAssociationHelper"),"/extdata/refGene.txt",sep="")
	
	con <- file(SetID_file, "r")
	flag2=0
	line=1
	out=list("SNP info"="","Genotype Matrix"="")
	while( length(line) != 0 ) {
     		line=readLines(con,n=1)
		a=strsplit(line,"\t")
		if (length(a)>=1){
			if (a[[1]][13]==Gene_name){
				chromosome=substr(a[[1]][3],4,20)
				b1=strsplit(a[[1]][10],",",fixed=TRUE)
				b2= strsplit(a[[1]][11],",",fixed=TRUE)
				n1=length(b1[[1]])
				n2=length(b2[[1]])
				for (i in 1:min(n1,n2)){
					
					pos_int1=as.numeric(b1[[1]][i])
					pos_int2=as.numeric(b2[[1]][i])
					geno_temp=GetGenotypesRegionVCF( VCF_info,chromosome,pos_int1,pos_int2);
					geno1=GetGenotypesVCF(VCF_info)
					if (geno1[[1]][1]==""){
						geno_temp=GetGenotypesRegionVCF( VCF_info,paste("chr",chromosome,sep=""),pos_int1,pos_int2);
						geno1=GetGenotypesVCF(VCF_info)
					}			
					while (geno1[[1]][1]!=""){
						if (flag2==0){out=geno1;flag2=1;} else {
							out[[1]]=rbind(out[[1]],geno1[[1]])
							out[[2]]=rbind(out[[2]],geno1[[2]])
						}
						geno1=GetGenotypesVCF(VCF_info)
	
					

					}

				}

			}
		
		}	
	}

	close(con)
	out1=list("SNP info"="","Genotype Matrix"="")
	if (out[[1]][1]==""){out1=out} else {
		kk=duplicated(out[[1]])
		out1[[1]]=out[[1]][!kk, ]
		out1[[2]]=out[[2]][!kk, ]
	}



	return(out1)

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










GetGeno<-function(VCF_info,chromosome=1,pos_start=NA,pos_end=NA){
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

	SetID_file=paste(path.package("VCFAssociationHelper"),"/extdata/refGene.txt",sep="")
	pos_status=rep(FALSE, pos_end-pos_start+1)
	con <- file(SetID_file, "r")
	flag2=0
	line=1
	out=list("SNP info"="","Genotype Matrix"="")
	while( length(line) != 0 ) {
     		line=readLines(con,n=1)
		a=strsplit(line,"\t")
		if (length(a)>=1){
		if (a[[1]][3]==paste("chr",chromosome,sep="")){
			b1=strsplit(a[[1]][10],",",fixed=TRUE)
			b2= strsplit(a[[1]][11],",",fixed=TRUE)
			n1=length(b1[[1]])
			n2=length(b2[[1]])
			for (i in 1:min(n1,n2)){
				pos_int1=as.numeric(b1[[1]][i])
				pos_int2=as.numeric(b2[[1]][i])
				flag1=0;
				if (pos_int1<=pos_end & pos_int2>= pos_start){
					if (pos_start<=pos_int1 ) {
						if (pos_end>=pos_int2){ flag1=1;} else{ pos_int2=pos_end;flag1=1;}} else { 
						pos_int1=pos_start;
						if (pos_end>=pos_int2){flag1=1;} else {pos_int2=pos_end;flag1=1;}
						 
					}
				}
				if (flag1==1){
					geno_temp=GetGenotypesRegionVCF( VCF_info,chromosome,pos_int1,pos_int2);
					geno1=GetGenotypesVCF(VCF_info)
								
					while (geno1[[1]][1]!=""){
						if (flag2==0){out=geno1;flag2=1;pos_status[as.numeric(geno1[[1]][3])-pos_start+1]=TRUE;} else {
							if (pos_status[as.numeric(geno1[[1]][3])-pos_start+1]==FALSE){
								out[[1]]=rbind(out[[1]],geno1[[1]])
								out[[2]]=rbind(out[[2]],geno1[[2]])
								pos_status[as.numeric(geno1[[1]][3])-pos_start+1]=TRUE;
							}
						}
						geno1=GetGenotypesVCF(VCF_info)

					}

				}

			}

		}
		}
	
	}

	close(con)
	GetGenotypesRegionVCF( VCF_info,chromosome,pos_end+1,10E15);
	return(out)



}



