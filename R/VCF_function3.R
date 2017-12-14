######read setid line by line##############
read.VCF.SetID<-function(VCF_file, format="GT",SetID_file=NULL, SetIDFormat) {
	if (is.null(SetID_file)) {
    		SetID_file=paste(path.package("VCFAssociationHelper"),"/extdata/refGene.txt",sep="")
		stop(sprintf("The SetID file is missing.\n"))
	}
	if (is.na(SetIDFormat)) {
		Msg<-sprintf("SetID Format is missing!\n")
		stop(Msg)
	}
	if (format!="GT"){
		if (format!="DS"){
			Msg<-sprintf("VCF_file Reading Format is wrong!\n")
			stop(Msg)
		}
	}
	SetID_file<-normalizePath(SetID_file ,mustWork =FALSE)
	VCF_file<-normalizePath(VCF_file ,mustWork =FALSE)
	line=1
	temp=OpenVCF(VCF_file,format)
	out=NULL
	if (SetIDFormat=="SNPID"){
		snppool=read.table(SetID_file)[,2]
		tpgeno=GetGenotypesVCF(temp)
		while (tpgeno[[1]][1]!="" & length(snppool)>0){
			if (tpgeno[[1]][1] %in% snppool){
				if (is.null(out)){
					out=tpgeno
				} else {
					out[[1]]=rbind(out[[1]],tpgeno[[1]])
					out[[2]]=rbind(out[[2]],tpgeno[[2]])
				}
			} 
			tpgeno=GetGenotypesVCF(temp)
		}
	} 
	if (SetIDFormat=="POS"){
		con <- file(SetID_file, "r")
		line=readLines(con,n=1)
		while( length(line) !=0 ) {
			if (line !=""){
				a=strsplit(line," ")
				b=strsplit(a[[1]][2],":",fixed=TRUE)
				chr_temp=as.numeric(b[[1]][1])
				pos_temp=as.numeric(b[[1]][2])
				GetGenotypesRegionVCF( temp,chr_temp,pos_temp,pos_temp)
				tpgeno=GetGenotypesVCF(temp)
				if (tpgeno[[1]][1]!=""){
					if (is.null(out)){
						out=tpgeno
					} else {
						out[[1]]=rbind(out[[1]],tpgeno[[1]])
						out[[2]]=rbind(out[[2]],tpgeno[[2]])
					}
				}
			}
			line=readLines(con,n=1)
		}
		close(con)
		
	} 
	if (SetIDFormat=="Region"){
		con <- file(SetID_file, "r")
		line=readLines(con,n=1)
		while( length(line) !=0) {
			if (line!=""){
				line1=line
				a=strsplit(line," ")
				GetGenotypesRegionVCF( temp,as.numeric(a[[1]][2]),as.numeric(a[[1]][3]),as.numeric(a[[1]][4]))
				tpgeno=GetGenotypesVCF(temp)
				while (tpgeno[[1]][1]!=""){
					if (is.null(out)){
						out=tpgeno
					} else {
						out[[1]]=rbind(out[[1]],tpgeno[[1]])
						out[[2]]=rbind(out[[2]],tpgeno[[2]])
					}
					tpgeno=GetGenotypesVCF(temp)			
				
				}
			}
			line=readLines(con,n=1)
			
		}
		close(con)
		
	} 
	CloseVCF()
	return(out)

}


read.VCF.filter<-function(VCF_file, format="GT",upperMAC=NA, lowerMAC=0, upperMAF=1,lowerMAF=0,nmax=1000) {
	if (upperMAC<lowerMAC & !is.na(upperMAC)){stop(sprintf("Wrong MAC setting!\n"))}
	if (upperMAF<lowerMAF){stop(sprintf("Wrong MAF setting!\n"))}
	VCF_file<-normalizePath(VCF_file ,mustWork =FALSE)
	line=1
	if (format!="GT"){
		if (format!="DS"){
			Msg<-sprintf("VCF_file Reading Format is wrong!\n")
			stop(Msg)
		}
	}
	temp=OpenVCF(VCF_file,format)
	out=NULL
	tempmin=max(lowerMAC/temp$samplesize,lowerMAF)
	if (is.na(upperMAC)) {tempmax=upperMAF} else {tempmax=min(upperMAC/temp$samplesize,upperMAF)}
	while (nmax>0){
		tpgeno=GetGenotypesVCF(temp)
		meantpgeno=mean(tpgeno[[2]])	
		if (meantpgeno>=tempmin & meantpgeno<=tempmax){
			if (is.null(out)){
				out=tpgeno
			} else {
				out[[1]]=rbind(out[[1]],tpgeno[[1]])
				out[[2]]=rbind(out[[2]],tpgeno[[2]])
			}
			nmax=nmax-1
		}
		if (tpgeno[[1]][1]==""){CloseVCF;return (out);}
	}
	CloseVCF()
	
	return(out)				

}


testyourfunction<-function(fun=NA,VCF_info=NA,Chr=NA,pos_start=NA,pos_end=NA,Y=NA,numcores=1){
	library( parallel)
	if (class(fun)!="function"){stop(sprintf("The testing funciton is missing!\n"))}
	if (is.na(VCF_info[[1]])){stop(sprintf("The VCF info is missing!\n"))}
	if (is.na(Chr)){stop(sprintf("The chromosome is missing!\n"))}
	if (is.na(pos_start) & !is.na(pos_end)){print(sprintf("The program will read from the position 1.\n"));pos_start=1}
	if (is.na(pos_end)){stop(sprintf("The end position is missing!\n"))}
	if (sum(is.na(Y))>0){stop(sprintf("The outcome Y is missing!\n"))}
	if (pos_end<pos_start){stop(sprintf("Wrong postion setting!\n"))}
	if (numcores!=round(numcores) ){stop(sprintf("Wrong number of cores!\n"))}
	if (numcores<1 ){stop(sprintf("Wrong number of cores!\n"))}	
	##affinity <- c(1:numcores)
	##mclapply(X = list(A,B,C), FUN = varFilter,
## mc.preschedule = FALSE, affinity.list = affinity))
	temp_sum=list()
	dif=round((pos_end-pos_start)/numcores)
	for (i in 1:numcores){
		temp=list()
		temp[[1]]=VCF_info
		temp[[2]]=Y
		if (i == numcores){
			temp[[3]]=c(Chr,pos_start+(i-1)*dif,pos_end)
		} else {
			temp[[3]]=c(Chr,pos_start+(i-1)*dif,pos_start+i*dif-1)
		}
		temp[[4]]=fun
		temp_sum[[i]]=temp
	}
	out_temp=mclapply(X = temp_sum, FUN = test_core, mc.cores = numcores,
 mc.preschedule = FALSE)
	out=out_temp[[1]]
	if (i >1){
		for (i in 2:numcores){
			out=rbind(out,out_temp[[i]])
		}
	}
	return(out)
	
}


test_core<-function(X){
	VCF_info=X[[1]]
	Y=X[[2]]
	Chr=X[[3]][1]
	pos_start=X[[3]][2]
	pos_end=X[[3]][3]
	fun=X[[4]]
	GetGenotypesRegionVCF( VCF_info,Chr,pos_start,pos_end)
	temp=GetGenotypesVCF(VCF_info)
	out=NULL
	while (temp[[1]][1]!=""){
		if (as.numeric(temp[[1]][3])<=pos_end & as.numeric(temp[[1]][3])>=pos_start){
			a=fun(t(temp[[2]]),Y)
			if (is.null(out)){
				out=t(c(temp[[1]][1],temp[[1]][2],temp[[1]][3],a))
			} else {
				out=rbind(out,t(c(temp[[1]][1],temp[[1]][2],temp[[1]][3],a)))
			}
		}
		temp=GetGenotypesVCF(VCF_info)
	}
	rownames(out) <- c()
	return(out)
}

