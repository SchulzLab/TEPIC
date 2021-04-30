args <- commandArgs(TRUE)
library('methods')
ggplotAvailable<-require("ggplot2")
gplotsAvailable<-require("gplots") 
fontSize=.9

if(length(args) < 1) {
  args <- c("--help")
}
 
if("--help" %in% args) {
  cat("
      Arguments:
      --outDir= Output path (required)
      --dataDir= Data directory (required)
      --out_var= Name of the response variable (required)
      --cores= Number of the cores to use (1 as default)
      --alpha= Alpha parameter stepsize (0.1 as default)
      --testsize= Size of test data (0.2 as default)
      --Ofolds= Number of outer folds for model validation (3 as default)
      --Ifolds= Number of inner cross validation folds (6 as default)
      --balanced= Flag indicating whether the data should be balanced through downsampling (default TRUE)
      --performance= Flag indicating whether performance measures should be computed (default TRUE)
      --randomise= Flag indicating whether a model should be learned on randomised data (default FALSE)	
")      
  q(save="no")
}

#Processing command line arguments
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

#Checking required arguments
if(is.null(argsL$outDir)) {
 cat("Error!\n\n  Please specify an output directory \n\n ")
q(save="no")
}
argsL$outDir<-paste0(argsL$outDir,"/")

if(is.null(argsL$dataDir)) {
 cat("Error!\n\n  Please specify the data directory \n\n ")
q(save="no")
}
argsL$dataDir<-paste0(argsL$dataDir,"/")

if(is.null(argsL$out_var)) {
 cat("Error!\n\n Please enter the name of the response variable \n\n ")
q(save="no")
}

#Checking for optional arguments
if(is.null(argsL$cores)) {
argsL$cores<-20
}
if(is.null(argsL$testsize)){
argsL$testsize<-0.2
}

if(is.null(argsL$Ofolds)){
argsL$Ofolds<-3
}

if(is.null(argsL$alpha)) {
argsL$alpha <-0.1
}

if(is.null(argsL$balanced)) {
argsL$balanced <-TRUE
}

if(is.null(argsL$performance)) {
argsL$performance <-TRUE
}


if (is.null(argsL$Ifolds)){
argsL$Ifolds <-6
}

if (is.null(argsL$randomise)){
argsL$randomise <- FALSE
}

#Creating output directory if necessary
dir.create(argsL$outDir,showWarning=FALSE,recursive=TRUE)

permute<-function(x){
s<-sample(length(x))
x[s]
}

#Loading required packages and initialising index variables
require(glmnet)
library("doMC")
i<-0
k<-0
counter<-0
registerDoMC(argsL$cores)
#Checking the total number of samples
FileList<-list.files(path=argsL$dataDir)
print(paste("Samples to analyse:",as.character(length(FileList)),sep=" "))
for(Sample in FileList){
	counter<-counter+1
	print(paste(as.character(counter),unlist(unlist(strsplit(Sample,".txt")))))
	}

#Creating the header for the overview table
SampleOverview<-c("Name","Mean Test Accuracy","Var Test Accuracy","Mean F1_1","Var F1_1","Mean F1_2","Var F1_2","Mean Training Accuracy","Var Train Accuracy")

#Declare elastic net functions
elaBinomial<-function(a,x,y,i){
	print(paste0("Learning model for alpha = ",a))
	elasticnet<-cv.glmnet(data.matrix(x), y, alpha=a,family="binomial",type.measure="class",parallel=TRUE,nfolds=i)
	min(elasticnet$cvm)
	}
elaMultinomial<-function(a,x,y,i){
	print(paste0("Learning model for alpha = ",a))
	elasticnet<-cv.glmnet(data.matrix(x), y, alpha=a,family="multinomial",type.measure="class",parallel=TRUE,nfolds=i)
	min(elasticnet$cvm)
	}

coefficients<-vector("list",length(FileList))

# Data Preprocessing 
for(Sample in FileList){
	# Process file names for plot title
	name<-strsplit(Sample, ".",fixed=TRUE)[[1]][1]
	print(paste0("Learning model for ",name))
	i<-i+1
	#Loading the data matrix
	M<-read.table(paste(argsL$dataDir,Sample,sep="/"),header=TRUE,stringsAsFactors=FALSE,row.names=1)
	#Removing features with standard deviation zero
	SD<-apply(M,2,sd)
	Feature_zero_SD<-as.vector(which(SD==0))
	if(length(Feature_zero_SD)>0){
		print("Warning! The following features are constant and will be removed from the feature matrix")
		print(colnames(M)[Feature_zero_SD])
		M<-M[,-c(Feature_zero_SD)]
	}
	#Initialising FeatureNames, Test and Training Accuracy vectors and retrieving the column containing the response variable
	FeatureName<-colnames(M)
	Response_Variable_location<- grep(argsL$out_var,FeatureName)
	TestAcc<-c(1:argsL$Ofolds)
	TrainAcc<-c(1:argsL$Ofolds)
	F1_1<-c(1:argsL$Ofolds)
	F1_2<-c(1:argsL$Ofolds)
	alphas=seq(0.0,1.0,as.numeric(argsL$alpha))
	coefficients[[i]]<-vector("list",as.numeric(argsL$Ofolds))
	#Normalising the data matrix
	M[,-Response_Variable_location]<-log2(M[,-Response_Variable_location]+1)
	M[,-Response_Variable_location]<-scale(M[,-Response_Variable_location],center=TRUE, scale=TRUE)
	if (argsL$randomise == TRUE){
          MP<-apply(M[,-Response_Variable_location],2,permute)
          M<-cbind(data.frame(scale(MP,center=TRUE, scale=TRUE)),M[,Response_Variable_location])
	      colnames(M)<-FeatureName
	}

	#Balancing the data set
	if (argsL$balanced==TRUE){
		subM<-list(list())
		bM<-c()
		index=0
		for (j in unique(M[,Response_Variable_location])){
			index<-index+1
			subM[[index]]<-M[which(M[,Response_Variable_location]==j),]
		}
		mSize=min(sapply(subM,dim)[1,])
		for (l in 1:length(subM)){
			rndselect=sample(x=nrow(subM[[l]]),size=mSize)
			subM[[l]]=subM[[l]][rndselect,]
			bM<-rbind(bM,subM[[l]])
		}
	}
	test_size<-1/as.numeric(argsL$testsize)

	if (argsL$performance){
	#Loop through the outer folds
	for (k in 1:argsL$Ofolds){
		print(paste("Outer CV fold ",as.character(k),sep=" "))
		if (argsL$balanced==TRUE){
			#Balanced selection of test and training data
			Test_Data<-c()
			Train_Data<-c()
			for (j in 1:length(subM)){
				rndselect=sample(x=nrow(subM[[j]]), size=mSize/test_size)
				Test_Data<-rbind(Test_Data,subM[[j]][rndselect,])
				Train_Data<-rbind(Train_Data,subM[[j]][-rndselect,])
			}
		}else{
			rndselect=sample(x=nrow(M),size=as.numeric(argsL$testsize)*nrow(M))
			Test_Data<-M[rndselect,]
			Train_Data<-M[-rndselect,]
		}

		#Split the features from response
		x_train<-data.matrix(Train_Data[,-Response_Variable_location])
		x_test<-data.matrix(Test_Data[,-Response_Variable_location])
		y_train<-as.vector(unlist(Train_Data[,Response_Variable_location,drop=FALSE]))
		y_test<-as.vector(unlist(Test_Data[,Response_Variable_location]))
		write.table(x_train,paste0(argsL$outDir,"/tmp1.txt"),sep="\t",quote=F)
		#Training model parameters in the inner cross validation
		print(argsL$Ifolds)
		if (length(unique(M[,Response_Variable_location]))==2){
			aError<-sapply(alphas,elaBinomial,x_train,y_train,as.numeric(argsL$Ifolds))
		}else{
			aError<-sapply(alphas,elaMultinomial,x_train,y_train,as.numeric(argsL$Ifolds))
		}

		#Determine best model and retrain it
		index<-which(aError==min(aError))
		if (length(unique(M[,Response_Variable_location]))==2){
			elasticnet<-cv.glmnet(x_train, y_train,alpha=alphas[index],family="binomial",type.measure="class",parallel=TRUE,nfolds=as.numeric(argsL$Ifolds))
		}else{
			elasticnet<-cv.glmnet(x_train, y_train,alpha=alphas[index],family="multinomial",type.measure="class",parallel=TRUE,nfolds=as.numeric(argsL$Ifolds))
		}

		#Generating a plot visualising the model selection
		svg(paste(paste0(argsL$outDir,"/Misclassification_vs_Lambda_Fold_",k,"_",name),"svg",sep="."))
		plot(elasticnet)
		dev.off()

		#Applying the selected  model to the test data
		predict_test<-predict(elasticnet, x_test, s="lambda.min",type="class")
		predict_train<-predict(elasticnet, x_train, s="lambda.min",type="class")

		#Printing and storing the performance table
		tTest<-table(y_test,predict_test)
		write.table(tTest,file=paste(argsL$outDir,"/Confusion-Matrix_",k,"_",name,".txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
		tTrain<-table(y_train,predict_train)

		#Storing Feature values
		if (length(unique(M[,Response_Variable_location]))==2){
			coefficients[[i]][[k]]<-coef(elasticnet, s = "lambda.min")
		}else{
			coefficients[[i]][[k]]<-c()
		}

		#Calculating Accuracy measures
		TestAcc[k]<-0
		TrainAcc[k]<-0
		F1_1[k]<-0
		F1_2[k]<-0
		for (index in 1:dim(tTest)[1]){
			TestAcc[k]<-TestAcc[k]+tTest[(index-1)*dim(tTest)[1]+index]
			TrainAcc[k]<-TrainAcc[k]+tTrain[(index-1)*dim(tTest)[1]+index]
			}
		F1_1[k]<-(2*tTest[1])/(2*tTest[1]+tTest[2]+tTest[3])
		F1_2[k]<-(2*tTest[4])/(2*tTest[4]+tTest[2]+tTest[3])
		TestAcc[k]<-TestAcc[k]/sum(tTest)
		TrainAcc[k]<-TrainAcc[k]/sum(tTrain)
		}
	
	#Storing the mean performance values
	sampleResult<-c(name,mean(TestAcc),var(TestAcc),mean(F1_1),var(F1_1),mean(F1_2),var(F1_2),mean(TrainAcc),var(TrainAcc))
	SampleOverview<-rbind(SampleOverview,sampleResult)
	}

	#Learning one model on the full data set for feature analysis	
	print("Learning model on the entire data set")
	if (argsL$balanced==TRUE){
		x_com<-as.matrix(bM[,-Response_Variable_location])
		y_com<-as.vector(unlist(bM[,Response_Variable_location,drop=FALSE]))
	}else{
		x_com<-as.matrix(M[,-Response_Variable_location])
		y_com<-as.vector(unlist(M[,Response_Variable_location,drop=FALSE]))
	}
	aError=0
	if (length(unique(M[,Response_Variable_location]))==2){
		aError<-sapply(alphas,elaBinomial,x_com,y_com,as.numeric(argsL$Ifolds))
	}else{
		aError<-sapply(alphas,elaMultinomial,x_com,y_com,as.numeric(argsL$Ifolds))
	}
	index<-which(aError==min(aError))
	if (length(unique(M[,Response_Variable_location]))==2){
		elasticnet<-cv.glmnet(x_com, y_com,alpha=alphas[index],family="binomial",type.measure="class",parallel=TRUE,nfolds=as.numeric(argsL$Ifolds))
	}else{
		elasticnet<-cv.glmnet(x_com, y_com,alpha=alphas[index],family="multinomial",type.measure="class",parallel=TRUE,nfolds=as.numeric(argsL$Ifolds))
	}

	#Store the feature values
	if (length(unique(M[,Response_Variable_location]))>2){
		coefO<-coef(elasticnet,s="lambda.min")
		for (j in names(coefO)){
			nf<-(coefO[j][[1]][,1])[-1]
			nf2<-nf/max(abs(nf))
			nf3<-t(nf2)
			nf3<-as.data.frame(nf3)
			nf4<-t(rbind(colnames(nf3),nf3))
			colnames(nf4)<-c("TF","value")
			nf4[,2]<-as.numeric(as.character(nf4[,2]))
			nf4<-as.data.frame(nf4)
			nf4[,2]<-as.numeric(as.character(nf4[,2]))
			write.table(nf4,file=paste(argsL$outDir,paste("Class",j,"Regression_Coefficients_Entire_Data_Set.txt",sep="-"),sep='/'),quote=FALSE,sep="\t",row.names=FALSE)
			nf4<-nf4[which(nf4$value >0.025),]
			np<-c(1:length((nf4[,2])))
			np[which(nf4[,2]>0)]<-1
			np[which(nf4[,2]<=0)]<-0
			if (ggplotAvailable){
				ggplot2::ggplot(nf4,aes(x=reorder(TF,value),y=value,width=0.8,fill=np))+
				geom_bar(stat="identity")+
				theme_bw(10)+ylab("Normalised coefficient")+xlab("TF")+
				theme(axis.text.x=element_text(angle=45,hjust=1))+
				theme(strip.background  = element_blank())+
				theme(legend.position="none")
				ggsave(paste0(argsL$outDir,"Regression_Coefficients_Entire_Data_Set_Class",j,"_",name,".png"),width=min(50,5+(.3*length(nf4$TF))),height=5,limitsize=F)
			}
		}
	}else{
			nf<-coef(elasticnet,s="lambda.min")[,1][-1]
			nf2<-nf/max(abs(nf))
			nf3<-t(nf2)
			nf3<-as.data.frame(nf3)
			nf4<-t(rbind(colnames(nf3),nf3))
			colnames(nf4)<-c("TF","value")
			nf4[,2]<-as.numeric(as.character(nf4[,2]))
			nf4<-as.data.frame(nf4)
			nf4[,2]<-as.numeric(as.character(nf4[,2]))
			write.table(nf4,file=paste(argsL$outDir,paste0("Regression_Coefficients_Entire_Data_Set_",name,".txt"),sep='/'),quote=FALSE,sep="\t",row.names=FALSE)
			nf4<-nf4[which(nf4$value !=0.0),]
			np<-c(1:length((nf4[,2])))
			np[which(nf4[,2]>0)]<-1
			np[which(nf4[,2]<0)]<-0
			if (ggplotAvailable){
				ggplot2::ggplot(nf4,aes(x=reorder(TF,value),y=value,width=0.8,fill=np))+
				geom_bar(stat="identity")+
				theme_bw(10)+ylab("Normalised coefficient")+xlab("TF")+
				theme(axis.text.x=element_text(angle=45,hjust=1))+
				theme(strip.background  = element_blank())+
				theme(legend.position="none")
				ggsave(paste0(argsL$outDir,"Regression_Coefficients_Entire_Data_Set",name,".png"),width=min(50,5+(.3*length(nf4$TF))),height=5,limitsize=F)
		}
	}	
}

if (argsL$performance){
	print("Performance measures:")
	print(paste("Mean Test Accuracy",mean(TestAcc),sep=" "))
	print(paste("Mean Training Accuracy",mean(TrainAcc),sep=" "))
	print(paste("Mean F1_1",mean(F1_1),sep=" "))
	print(paste("Mean F1_2",mean(F1_2),sep=" "))	

	for (i in 1:length(FileList)){
		featureMatrix<-c()
		for (j in 1:length(coefficients[[i]])){
			if (length(coefficients[[i]][[j]]>1)){
				featureMatrix<-rbind(featureMatrix,coefficients[[i]][[j]][,1])
			}
		}
		if (length(featureMatrix > 1)){
			write.table(featureMatrix,paste(argsL$outDir,"Regression_Coefficients_Cross_Validation_",FileList[i],sep=""),quote=F,sep="\t",col.names=NA)
			meanFeature<-apply(featureMatrix,2,median)
			featureMatrix<-featureMatrix[,-1]
			meanFeature<-meanFeature[-1]
			featureMatrix<-featureMatrix[,which(meanFeature!=0.0)]
			meanFeature<-meanFeature[which(meanFeature!=0.0)]
			if(length(meanFeature) > 1){
				allFeatures<-featureMatrix[,order(meanFeature,decreasing=TRUE)]
				meanFeatures<-meanFeature[order(meanFeature,decreasing=TRUE)]
				all<-rbind(allFeatures,meanFeatures)
				all<-all[,order(all[dim(all)[1],],decreasing=TRUE)]
				row.names(all)<-c(paste("Fold ",c(1:(dim(all)[1]-1))),"Median")
				if (gplotsAvailable){
					library("gplots")
					svg(paste(argsL$outDir,"Regression_Coefficients_Cross_Validation_Heatmap_",unlist(unlist(strsplit(FileList[i],".txt")))[1],".svg",sep=""),width=min(54,7+(.3*length(meanFeature))),height=min(50,5+(.5*as.numeric(argsL$Ofolds))))
					if(any(allFeatures < 0)){
						allFeatures<-featureMatrix[,order(meanFeature,decreasing=TRUE)]
						meanFeature<-meanFeature[order(meanFeature,decreasing=TRUE)]							
						all<-rbind(allFeatures,meanFeatures)
						all<-all[,order(all[dim(all)[1],],decreasing=TRUE)]
						row.names(all)<-c(paste("Fold ",c(1:(dim(all)[1]-1))),"Median")
						heatmap.2(all,trace="none",col=bluered(250),srtCol=45,cexRow=fontSize,cexCol=fontSize,density.info="none",distfun = distF, dendrogram="none", margins=c(6,6),Colv=FALSE,Rowv=FALSE,key.xlab="Regression coefficient",keysize=0.9) 
					}else{
						allFeatures<-featureMatrix[,order(meanFeature,decreasing=TRUE)]
						meanFeatures<-meanFeature[order(meanFeature,decreasing=TRUE)]
						all<-rbind(allFeatures,meanFeatures)
						all<-all[,order(all[dim(all)[1],],decreasing=TRUE)]
						row.names(all)<-c(paste("Fold ",c(1:(dim(all)[1]-1))),"Median")
						heatmap.2(all,trace="none",col=heat.colors(250),srtCol=45,cexRow=fontSize,cexCol=fontSize,density.info="none",distfun = distF, dendrogram="none", margins=c(6,6),Colv=FALSE,Rowv=FALSE,key.xlab="Regression coefficient",keysize=0.9) 
					}
					dev.off()
				}
			}else{
				if (ggplotAvailable){
					library("ggplot2")
					text<-"Mean regression coefficients are to small, a heatmap can not be shown."
					ggplot()+annotate("text",x=4,y=25,size=8,label=text)+theme_void() 
					ggsave(filename=paste(argsL$outDir,"Regression_Coefficients_Cross_Validation_Heatmap_",unlist(unlist(strsplit(FileList[i],".txt")))[1],".png",sep=""),width=11,height=4)
				}
			}
		}
	}
	#Generating summary output table
	colnames(SampleOverview)<-SampleOverview[1,]
	SampleOverview<-SampleOverview[-1,]
	if(class(SampleOverview)=="matrix"){
		ggplotSampleOverview<-as.data.frame(rbind(cbind(SampleOverview[,1:3],rep("Test",length(SampleOverview[,1]))),cbind(SampleOverview[,c(1,8,9)],rep("Training",length(SampleOverview[,1]))),cbind(SampleOverview[,c(1,4,5)],rep("F1_1",length(SampleOverview[,1]))),cbind(SampleOverview[,c(1,6,7)],rep("F1_2",length(SampleOverview[,1])))))
	}else{
		ggplotSampleOverview<-as.data.frame(rbind(c(SampleOverview[1:3],"Test Accuracy"),c(SampleOverview[c(1,8,9)],"Trainings Accuracy"),c(SampleOverview[c(1,4,5)],"F1-Up"),c(SampleOverview[c(1,6,7)],"F1-Down")))
	}

	colnames(ggplotSampleOverview)<-c("Name","Mean","Variance","Measure")
	row.names(ggplotSampleOverview)<-c(1:length(row.names(ggplotSampleOverview)))
	ggplotSampleOverview[,2]<-as.numeric(as.character(ggplotSampleOverview[,2]))
	ggplotSampleOverview[,3]<-as.numeric(as.character(ggplotSampleOverview[,3]))
	write.table(ggplotSampleOverview,file=paste0(argsL$outDir,"/Performance_overview.txt"),quote=FALSE,sep='\t',row.names=FALSE)
	if (ggplotAvailable){
		ggplot2::ggplot(ggplotSampleOverview,aes(x=Name,y=Mean,width=0.8,fill=Measure,group=Measure))+
		geom_bar(stat="identity",position="dodge")+
		theme_bw(22)+ylab("Value")+xlab("Sample")+
		theme(axis.text.x=element_text(angle=45,hjust=1))+
		scale_y_continuous(breaks=seq(0.0,1.0,0.05))+
		scale_fill_manual(values=c("red","blue","orange","black"))+
		coord_cartesian(ylim=c(0.0,1.0))+
		geom_errorbar(aes(ymin=Mean-sqrt(Variance),ymax=Mean+sqrt(Variance)),position=position_dodge(0.8),width=.4,size=.5)+
		theme(strip.background  = element_blank())
		ggsave(paste0(argsL$outDir,"/Performance_Barplots.png"),width=20,height=20,units=c("cm"))
	}
}
