######read setid line by line##############
read.VCF.SetID<-function(VCF_file, format="GT",SetID_file=NULL, SetIDFormat) {
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
	if (format!="GT"){
		if (format!="DS"){
			Msg<-sprintf("VCF_file Reading Format is wrong!\n")
			stop(Msg)
		}
	}
	SetID_file<-normalizePath(SetID_file ,mustWork =FALSE)
	SetID<-read.table(SetID_file)
	if (SetIDFormat=="SNPID"){
		snppool=SetID[,2]
		temp=OpenVCF(VCF_file,format)
		tpgeno=GetGenotypesVCF(temp)
		out=NULL
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
		CloseVCF()
	} 
	if (SetIDFormat=="POS"){
		out=NULL
		for (i in 1:nrow(SetID)){
			temp=OpenVCF(VCF_file,format)
			tp=strsplit(as.character(SetID[i,2]), "\\:")
			tp1=tp[[1]][1]
			tp2=tp[[1]][2]
			GetGenotypesRegionVCF(temp,tp1,tp2,tp2)
			tpgeno=GetGenotypesVCF(temp)
			if (is.null(out)){
				out=tpgeno
			} else {
				out[[1]]=rbind(out[[1]],tpgeno[[1]])
				out[[2]]=rbind(out[[2]],tpgeno[[2]])
			}
			CloseVCF()

		}
	} 
	if (SetIDFormat=="Region"){
		out=NULL
		for (i in 1:nrow(SetID)){
			temp=OpenVCF(VCF_file,format)
			
			GetGenotypesRegionVCF(temp,SetID[i,2],SetID[i,3],SetID[i,4])
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
			CloseVCF()

		}
	} 
	return(out)

}

