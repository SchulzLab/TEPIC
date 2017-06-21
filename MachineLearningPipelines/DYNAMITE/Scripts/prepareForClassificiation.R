args<-commandArgs(TRUE)

data<-read.table(args[1],header=TRUE,stringsAsFactors=FALSE)
u<-which(data$Expression > 0)
n<-which(data$Expression < 0)
dataN<-data[,-1]
dataN$Expression[n]<-0
dataN$Expression[u]<-1
dataN$Expression<-as.factor(dataN$Expression)
colnames(dataN)<-c(colnames(dataN))
write.table(dataN,args[2],quote=FALSE,sep="\t",col.names=NA,row.names=T)
