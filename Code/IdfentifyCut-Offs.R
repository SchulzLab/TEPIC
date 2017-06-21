library(data.table)

#arg1 file containing TF affinities computed on the randomised data
#arg2 file containing TF affinities computed on the actual data
#arg3 output file
#arg4 cut off used
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4){
	print("Usage: <Affinities computed on random set> <Affinities computed on the actual data> <Output file> <Cut off (The quantile used is 1-value)>")
	q(save="no")
}

#Function defintions
getLength<-function(x){as.numeric(x[3])-as.numeric(x[2])}
normalise<-function(x,y){x/y}
determineQuantile<-function(x,y){quantile(x,c(y))}

#Setting file names
randomRegionsFile=args[1]
originalRegionsFile=args[2]

##Computing quantiles
#Loading affinity data based on random regions
randomData<-fread(randomRegionsFile,stringsAsFactors = F)
sLength<-strsplit(randomData$V1,"[[:punct:][:punct:]]")

#Computing length of oc region
lengthV<-rep(0,length(sLength))
for (i in c(1:length(sLength))){
  lengthV[i]<-getLength(sLength[[i]])
}

#Normalise affinities
randomDataNormalised<-randomData[,-1]
randomDataNormalised<-apply(randomData[,-1],2,normalise,lengthV)

#Compute quantiles
quantiles<-apply(randomDataNormalised,2,determineQuantile,1-as.numeric(args[4]))

#Remove objects that are not required anymore
rm(randomData,randomDataNormalised,sLength)

#Load affinity data computed on the actual data
originalData<-fread(originalRegionsFile,stringsAsFactors = F)
sLength<-strsplit(originalData$V1,"[[:punct:][:punct:]]")
lengthV<-rep(0,length(sLength))
for (i in c(1:length(sLength))){
  lengthV[i]<-getLength(sLength[[i]])
}

#Normalise the data
originalDataNormalised<-originalData[,-1]
originalDataNormalised<-apply(originalData[,-1],2,normalise,lengthV)

#Filter the data
reprocessedData<-originalData
for (i in c(2:dim(originalData)[2])){
  reprocessedData[which(originalDataNormalised[,i-1]<quantiles[i-1]),i]<-0.0
}
colnames(reprocessedData)[1]<-"\t"
write.table(reprocessedData,args[3],quote=F,sep="\t",col.names=TRUE,row.names=FALSE)
