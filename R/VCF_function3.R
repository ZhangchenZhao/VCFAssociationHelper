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



