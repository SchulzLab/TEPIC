args<-commandArgs(TRUE)

library("tRap")
library("parallel")
library("Biostrings")

sequences<-readDNAStringSet(args[1])
seql<-length(sequences)

#Load Data
jaspar<-read.jaspar(args[4]);
lj<-length(jaspar)

#Generate matrix column names containing TF names
cnames<-c()
for (m in jaspar)
{
    cnames<-c(cnames,attr(m,"accession"))
}


#Preallocation output
pS<-c(1:seql)

trapAF<-function(j){
	s<-toString(sequences[j])
	lenS<-nchar(s)
	l<-c(1:(lj))
	counter=1
	for (mat in jaspar){
		if (lenS > ncol(mat)){
			af<-affinity(mat,s,gc.content=0.41)
			}
		else {	
			af<-0.0
			}
		l[counter]<-round(af,10)
		counter=counter+1
	}
	l
}

#Run tRap
cl<-makeCluster(args[3],"FORK")
m<-parSapply(cl,pS,trapAF)
stopCluster(cl)

#Divide result into affinity and p-value matrices
Affinity=m
Affinity<-t(Affinity)

if (lj > 1) {
	dimnames(Affinity)<-list(names(sequences),cnames)
	write.table(Affinity,file=args[2],sep='\t',col.names=NA,quote=FALSE,na="0.0")
	} else {
	dimnames(Affinity)<-list(cnames,names(sequences))
	Affinity<-as.matrix(Affinity)
	Affinity<-t(Affinity)
	write.table(Affinity[,1],file=args[2],sep='\t',quote=FALSE,col.names=NA)
}
