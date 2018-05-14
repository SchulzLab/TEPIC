args <- commandArgs(TRUE)
library("methods")
library("glmnet")
library("doMC")
library("parallel")
if(length(args) < 1) {
	args <- c("--help")
}
 
## Help 
if("--help" %in% args) {
	cat("
	INVOKE offers linear regression with Lasso, Ridge, and Elastic Net regularisation.
	Arguments:
	--outDir Output directory (will be created if it does not exist)
	--dataDir Directory containing the data
	--response Name of the response variable
	--cores Number of cores to be use (default 1)
	--fixedAlpha Use a fixed value for the alpha parameter in elastic net regulatisation, do not perform a grid search
	--alpha Stepsize to optimise the alpha parameter in elastic net regularisation (default 0.05)
	--testsize Size of test data[%] (default 0.2)
	--regularisation L for Lasso, R for Ridge, and E for Elastic net (default E)
	--innerCV Number of folds for inner cross-validation (default 6)
	--outerCV Number of iterations of outer cross-validation to determine test error (default 3)
	--constraint Specifies a constraint on the coefficent sign, enter N for negative and P for positive constraint
	--performance Flag indiciating whether the performance of the model should be assessed (default TRUE)
	--seed Random seed used for random number generation (default random)
	--leaveOneOutCV Flag indicating whether a leave one out cross-validation should be used (default FALSE)
	--asRData Store feature coefficients as RData files (default FALSE)
	--randomise Randomise the feature matrix (default FALSE) 
	--logResponse Flag indicating whether the response variable should be log transformed (default TRUE)
	--ftest Flag indicating whether partial F-test should be computed to assess the significance of each feature (default FALSE)
	--help=print this text
")
	q(save="no")
}

# Process command arguments
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

if(is.null(argsL$outDir)) {
	cat("No output directory specified. Use the --outDir option to specify an output directory.")
	q(save="no")
}
argsL$outDir<-paste0(argsL$outDir,"/")

if(is.null(argsL$dataDir)) {
	cat("No data directory specified. Use the --dataDir option to specify a data directory.")
	q(save="no")
}
argsL$dataDir<-paste0(argsL$dataDir,"/")
Data_Directory <- argsL$dataDir

if(is.null(argsL$response)) {
	cat("No response variable name specified. Use the --response option to specify a response variable.")
	q(save="no")
}

if(is.null(argsL$ftest)){
	argsL$ftest<-FALSE
}

if(is.null(argsL$testsize)){
	argsL$testsize <- 0.2
}

if(is.null(argsL$innerCV)){
	argsL$innerCV <- 6
}


if(is.null(argsL$outerCV)){
	argsL$outerCV<- 3
}

if (as.numeric(argsL$outerCV) < 2){
	cat("Number of outer cross validation folds must be at least 2")
	q(save="no")
}

if(is.null(argsL$alpha)) {
	argsL$alpha <- 0.05
}

if(is.null(argsL$cores)) {
	argsL$cores <- 1
}

if(is.null(argsL$regularisation)){
	argsL$regularisation<-c("E")
}

if(is.null(argsL$constraint)){
	lower_bound <- NULL
	upper_bound <- NULL
}else if(argsL$constraint=="P"){
	lower_bound <- 0
}else if(argsL$constraint=="N"){
	upper_bound <- 0
}

if(is.null(argsL$fixedAlpha)){
	argsL$fixedAlpha <- -1
}

if (is.null(argsL$performance)){
	argsL$performance <- TRUE
}

if (! is.null(argsL$seed)){
	set.seed(as.numeric(argsL$seed))
}

if(is.null(argsL$leaveOneOutCV)){
     argsL$leaveOneOutCV <- FALSE
}

if (is.null(argsL$asRData)){
	argsL$asRData <- FALSE
}

if (is.null(argsL$randomise)){
	argsL$randomise <- FALSE
}

if (is.null(argsL$logResponse)){
	argsL$logResponse <- TRUE
}

if ((argsL$leaveOneOutCV==TRUE) & (argsL$ftest==TRUE)){
     cat("The F-Test can not be combined with leave one out cross validation.")
     q(save="no")
	}

registerDoMC(cores = argsL$cores)

permute<-function(x,resPos){
s<-sample(length(x))
s<-s[which(s != resPos)]
c(x[s],x[resPos])
}

#Check output directory, create it if necessary
dir.create(argsL$outDir,showWarning=FALSE)

#Initilaise lists for storage of intermediate results
FileList<-list.files(path=Data_Directory)
numFiles=length(FileList)
pearson_correlation<-vector("list",numFiles)
spearman_correlation<-vector("list",numFiles)
test_error<-vector("list",numFiles)
rss_error<-vector("list",numFiles)
ftest_result<-vector("list",numFiles)
coefficients<-vector("list",numFiles)
coefficientsF<-vector("list",numFiles)
Sample_View<-vector("list",numFiles)
validSamples<-vector("logical",numFiles)

#Print sample names
print(paste("Total number of samples:",as.character(length(FileList)),sep=" "))
if (length(FileList)==0){
	print("No samples available! Aborting")
	exit()
}
counter<-0
for(Sample in FileList){
	counter<-counter+1
	print(paste(as.character(counter),unlist(unlist(strsplit(Sample,".txt")))))
	}

#Loop through sample files
i<-0
for(Sample in FileList){
	i<-i+1
	print(paste("Learning sample ",as.character(i),sep=" "))	
	#Loading and preprocessing data
	print("Processing sample matrix. This can take a few minutes. Please wait.")
	M<-read.table(paste(Data_Directory,Sample,sep=""),header=TRUE,sep="",row.names=1)
	
	M<-unique(M)
	M<-data.frame(M)
	FeatureNames_temp<-colnames(M)
	Response_Variable_location_temp <- grep(argsL$response,FeatureNames_temp)
	if (min(M[,Response_Variable_location_temp]) >= 0){
		if (argsL$logResponse == TRUE){
			M<-log2(M+1)
		}
	}else{
		print("Applying log2 transformation only to features")
		Response_Variable_location_temp <- grep(argsL$response,FeatureNames_temp)
		M[,-Response_Variable_location_temp]<-log2(M[,-Response_Variable_location_temp]+1)
	}
	SD<-apply(M,2,sd)
	Feature_zero_SD<-as.vector(which(SD==0))
	if(length(Feature_zero_SD)>0){
		print("Warning, there are constant features. These are not considered for further analysis.")
		if (Response_Variable_location_temp %in% Feature_zero_SD){
			print("Warning, response is constant, this sample is excluded")
			validSamples[i]=FALSE
			next;
			}
		M<-M[,-c(Feature_zero_SD)]
	}
	if (is.null(dim(M))){
		validSamples[i]=FALSE
		print("Warning, sample matrix is null, this sample is excluded")
		next;
	}
	if (dim(M)[2] < 2){
		print("Warning, no data included")
		validSamples[i]=FALSE
		next;
	}
	print(length(which(M==0)))
	if (length(which(M==0))>(dim(M)[1]*dim(M)[2]*0.5)){
		validSamples[i]=FALSE
		print("Warning, insufficient data coverage")
		next;
	}

	FeatureNames<-colnames(M)
	M<-data.frame(scale(M,center=TRUE, scale=TRUE))
     if (dim(M)[1] < 30){
          print("Warning, less then 30 samples available. This file is not processed.")
          validSamples[i]=FALSE
          next;
     }else{
          validSamples[i]=TRUE
     }

	if (argsL$leaveOneOutCV==FALSE){
		vectorLength=as.numeric(argsL$outerCV)
		pearson_correlation[[i]]<-vector("list",vectorLength)
		spearman_correlation[[i]]<-vector("list",vectorLength)
		test_error[[i]]<-vector("list",vectorLength)
		coefficients[[i]]<-vector("list",vectorLength)
	}else{
		vectorLength<-nrow(M)    
		pearson_correlation[[i]]<-cbind(vector("list",vectorLength),vector("list",vectorLength))
		spearman_correlation[[i]]<-cbind(vector("list",vectorLength),vector("list",vectorLength))
		test_error[[i]]<-vector("list",vectorLength)
		coefficients[[i]]<-vector("list",vectorLength)
	}
	if (argsL$ftest ==TRUE){
		ftest_result[[i]]<-vector("list",dim(M)[2]-1)
	}

	name<-unlist(unlist(strsplit(Sample, ".txt")))
	Response_Variable_location<- grep(argsL$response,FeatureNames)
	#Randomise the data
	if (argsL$randomise == TRUE){
		MP<-t(apply(M,1,permute,Response_Variable_location))
		colnames(MP)<-colnames(M)
		M<-data.frame(scale(MP,center=TRUE, scale=TRUE))  
	}

	if (argsL$performance == TRUE){
		if (argsL$leaveOneOutCV==FALSE){
		#Looping through the outer folds
		for (k in 1:argsL$outerCV){
			print(paste("Outer cross validation fold: ",as.character(k),sep=" "))
			# Partition data into test and training data sets
			Test_size<-round(nrow(M)/(1/as.numeric(argsL$testsize)))
			rndselect<-sample(x=nrow(M), size=Test_size)
			Test_Data<-M[rndselect,]
			Train_Data<-M[-rndselect,]

			# Split the features from response
			x_train<-as.matrix(Train_Data[,-Response_Variable_location])
			x_test<-as.matrix(Test_Data[,-Response_Variable_location])
			y_train<-as.vector(unlist(Train_Data[,Response_Variable_location,drop=FALSE]))
			y_test<-as.vector(unlist(Test_Data[,Response_Variable_location]))

			#Creating alpha vector
		 	A<-c()
			if(argsL$regularisation=="L"){
				alphaslist <- c(1.0)
				print("The value of alpha is set to 1.0 (Lasso penalty)")
			}else{
				if(argsL$regularisation=="R"){
				alphaslist <- c(0.0)
				print("The value of alpha is set to 0.0 (Ridge penalty)")
				}else{
					alphaslist<-seq(0,1,by=as.numeric(argsL$alpha))
				}
			}
			#Learning model on training data
			if(argsL$regularisation=="E"){   
				if(argsL$fixedAlpha==-1){
				print("Elastic net with grid search")
					if(is.null(argsL$constraint)){
						elasticnet<-mclapply(alphaslist, function(x){print(paste0("Outer cross validation fold=",k,", Testing alpha=",x));cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV))}, mc.cores=argsL$cores)
					}else{ 
						if(argsL$constraint=="P"){
							elasticnet<-mclapply(alphaslist, function(x){print(paste0("Outer cross validation fold=",k,",Testing alpha=",x));cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV))}, mc.cores=argsL$cores)
							}else{
								if(argsL$constraint=="N"){
								elasticnet<-mclapply(alphaslist, function(x){print(paste0("Outer cross validation fold=",k,",Testing alpha=",x));cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV))}, mc.cores=argsL$cores)
								}
							}
			      		}
				}else{
				x=argsL$fixedAlpha
				print(paste("Elastic net with fixed alpha",x,sep=" "))
				if(is.null(argsL$constraint)){
					elasticnet<-cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
				}else{ 
					if(argsL$constraint=="P"){
							elasticnet<-cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
					}else{
						if(argsL$constraint=="N"){
							elasticnet<-cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
								}
						}
		      	}
				}   
			}else{
			x=alphaslist[1]
			if(argsL$regularisation=="L"){
				print(paste0("Outer cross validation fold=",k,", Learning model using Lasso penalty"))
			}else{
					print(paste0("Outer cross valdiation fold=",k,", Learning model using Ridge penalty"))
			}
		
			if(is.null(argsL$constraint)){
					elasticnet<-cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
					}else{ 
					if(argsL$constraint=="P"){
						elasticnet<-cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
					}else{
						if(argsL$constraint=="N"){
							elasticnet<-cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
							}
						}
	     			}
				}
			if(length(elasticnet[[1]]) > 1){
				if (argsL$regularisation=="E"){
						if(argsL$fixedAlpha==-1){
							for (j in 1:length(alphaslist)) {
								A[j]<-min(elasticnet[[j]]$cvm)
							}
							#Determine best alpha value from training data
							index<-which(A==min(A), arr.ind=TRUE)
							model<-elasticnet[[index]]
					}
				}else{
					model<-elasticnet
					}
			}
			if (length(elasticnet[[1]]) > 1){ 
				#Determine error of the best alpha model on hold out data and on training data
				predict_fit<-predict(model, x_test, s="lambda.min")
				predict_fit_train<-predict(model, x_train, s="lambda.min")
				coefficients[[i]][[k]]<-coef(model, s = "lambda.min")
				png(paste0(argsL$outDir,"Scatter_Plot_",name,"_outer_fold_",k,".png"),width=800,height=800)
				plot(predict_fit,y_test,xlab="Predicted gene expression",ylab="Measured gene expression")
				dev.off()
				if (var(predict_fit,y_test) > 0){
					pearson_correlation[[i]][k]<-cor(predict_fit,y_test)
					spearman_correlation[[i]][k]<-cor(predict_fit,y_test,method='spearman')
				}else{
					pearson_correlation[[i]][k]<-0.0
					spearman_correlation[[i]][k]<-0.0
				}
				test_error[[i]][k]<-sum((y_test-predict_fit)^2)/length(y_test)
				rss_error[k]<-sum((y_test-predict_fit)^2)
			}else{
				coefficients[[i]][k]<-c()
				pearson_correlation[[i]][k]<-0.0
				spearman_correlation[[i]][k]<-0.0
				test_error[[i]][k]<-1.0
				rss_error[k]<-1.0
				validSamples[i]=FALSE
			}
		}
		if (argsL$ftest==TRUE){
			numFeatures=dim(M)[2]
			if (numFeatures-1 > 1){
				featureMatrixF<-c()
		          for (j in 1:length(coefficients[[i]])){
		               if (length(coefficients[[i]][[j]]>1)){
		                    featureMatrixF<-rbind(featureMatrixF,coefficients[[i]][[j]][,1])
		               }
		          }
				featureMatrixF<-featureMatrixF[,-1]
		          if (length(featureMatrixF > 1)){
						ftest_result[[i]]<-rep(1.0,numFeatures-1)
		                    meanFeature<-apply(featureMatrixF,2,median)
						nonZeroFeatures<-which(meanFeature != 0)
						print(nonZeroFeatures)
						print(paste0("Running F-test using ",length(nonZeroFeatures)," features"))
						featureIndex<-nonZeroFeatures
						for (m in featureIndex){
							print(paste0("Excluding feature ",m," ",FeatureNames[m]))
							MF<-M[,-m]
								error<-c(1:argsL$outerCV)
							for (k in 1:argsL$outerCV){
									Test_size<-round(nrow(MF)/(1/as.numeric(argsL$testsize)))
									rndselect<-sample(x=nrow(MF), size=Test_size)
									Test_Data<-MF[rndselect,]
									Train_Data<-MF[-rndselect,]
									Response_Variable_location_MF<- grep(argsL$response,colnames(MF))
									x_train<-as.matrix(Train_Data[,-Response_Variable_location_MF])
									x_test<-as.matrix(Test_Data[,-Response_Variable_location_MF])
										y_train<-as.vector(unlist(Train_Data[,Response_Variable_location_MF,drop=FALSE]))
									y_test<-as.vector(unlist(Test_Data[,Response_Variable_location_MF]))
							 		A<-c()
									if(argsL$regularisation=="L"){
										alphaslist <- c(1.0)
									}else{
										if(argsL$regularisation=="R"){
											alphaslist <- c(0.0)
										}else{
											alphaslist<-seq(0,1,by=as.numeric(argsL$alpha))
										}
									}
								if(argsL$regularisation=="E"){   
									if(argsL$fixedAlpha==-1){
										if(is.null(argsL$constraint)){
											elasticnet<-mclapply(alphaslist, function(x){cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV))}, mc.cores=argsL$cores)
										}else{ 
											if(argsL$constraint=="P"){
												elasticnet<-mclapply(alphaslist, function(x){cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV))}, mc.cores=argsL$cores)
											}else{
												if(argsL$constraint=="N"){
													elasticnet<-mclapply(alphaslist, function(x){cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV))}, mc.cores=argsL$cores)
												}
											}
							      		}
									}else{
										x=argsL$fixedAlpha
										if(is.null(argsL$constraint)){
											elasticnet<-cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
										}else{ 
										if(argsL$constraint=="P"){
											elasticnet<-cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
										}else{
											if(argsL$constraint=="N"){
												elasticnet<-cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
													}
												}
								     	 	}
										}			   
								}else{
								x=alphaslist[1]
								if(is.null(argsL$constraint)){
									elasticnet<-cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
								}else{ 
									if(argsL$constraint=="P"){
										elasticnet<-cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
									}else{
										if(argsL$constraint=="N"){
											elasticnet<-cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
											}
										}
					     			}
								}
							if(length(elasticnet[[1]]) > 1){
								if (argsL$regularisation=="E"){
									if(argsL$fixedAlpha==-1){
										for (j in 1:length(alphaslist)) {
											A[j]<-min(elasticnet[[j]]$cvm)
										}
									index<-which(A==min(A), arr.ind=TRUE)
									model<-elasticnet[[index]]
									}
								}else{
									model<-elasticnet
									}
								}
							if (length(elasticnet[[1]]) > 1){ 
								predict_fit<-predict(model, x_test, s="lambda.min")
								predict_fit_train<-predict(model, x_train, s="lambda.min")
								error[k]<-sum((y_test-predict_fit)^2)	
							}else{
								error[k]<-1.0
							}
						}
						mRSS<-mean(error)
						print(mRSS)
						mORSS<-mean(unlist(rss_error[[i]]),na.rm=TRUE)
						print(mORSS)
						MSe<-mORSS/(dim(MF)[1]-(length(nonZeroFeatures)))
						print(MSe)
						partialRSS<-mRSS-mORSS
						print(partialRSS)
						fvalue<-partialRSS/MSe
						print(fvalue)
						pValue<-1.0-pf(as.numeric(fvalue),length(nonZeroFeatures)-1,(dim(MF)[1]-length(nonZeroFeatures)))
						ftest_result[[i]][m]<-pValue
						}
					}
				###Generate output
			}else{
				print(paste0("Computing significance for feature ",FeatureNames[1]))
				mORSS<-mean(unlist(rss_error[[i]]),na.rm=TRUE)
				MSe<-mORSS/(dim(MF)[1]-2)
				fvalue<-mORSS/MSe
				pValue<-1.0-pf(as.numeric(fvalue),1,dim(MF)[1]-2)
				ftest_result[[i]][1]<-pValue
			}
	}
	}else{
	#Leave one out cross validation
	for (k in 1:nrow(M)){                                                                                                                                                                                                                                                                                                                      
		#print(paste("Leave one out cross validation fold: ",as.character(k),sep=" "))
		# Partition data into test and training data sets
		Test_Data<-M[k,]
		Train_Data<-M[-k,]
		# Split the features from response
		x_train<-as.matrix(Train_Data[,-Response_Variable_location])
		x_test<-as.matrix(Test_Data[,-Response_Variable_location])
		y_train<-as.vector(unlist(Train_Data[,Response_Variable_location,drop=FALSE]))
		y_test<-as.vector(unlist(Test_Data[,Response_Variable_location]))
		#Creating alpha vector
		A<-c()
		if(argsL$regularisation=="L"){
			print("The value of alpha is set to 1.0 (Lasso penalty)")
			alphaslist <- c(1.0)
		}else{
		if(argsL$regularisation=="R"){
			print("The value of alpha is set to 0.0 (Ridge penalty)")
			alphaslist <- c(0.0)
			}else{
				alphaslist<-seq(0,1,by=as.numeric(argsL$alpha))
			}
		}
	
		#Learning model on training data
		if(argsL$regularisation=="E"){   
			if(argsL$fixedAlpha==-1){
		if(is.null(argsL$constraint)){
			elasticnet<-mclapply(alphaslist, function(x){cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV))}, mc.cores=argsL$cores)
		}else{ 
			if(argsL$constraint=="P"){
			elasticnet<-mclapply(alphaslist, function(x){cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV))}, mc.cores=argsL$cores)
			}else{
				if(argsL$constraint=="N"){
					elasticnet<-mclapply(alphaslist, function(x){cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV))}, mc.cores=argsL$cores)
					}
				}
			}
		}else{
			x=argsL$fixedAlpha
			if(is.null(argsL$constraint)){
				elasticnet<-cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
				}else{ 
					if(argsL$constraint=="P"){
						elasticnet<-cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
					}else{
						if(argsL$constraint=="N"){
							elasticnet<-cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
						}
					}
				}
			}   
	}else{
	x=alphaslist[1]
	if(is.null(argsL$constraint)){
		elasticnet<-cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
		}else{ 
			if(argsL$constraint=="P"){
				elasticnet<-cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
				}else{
				if(argsL$constraint=="N"){
					elasticnet<-cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
					}
				}
			}
		}
	
	if(length(elasticnet[[1]]) > 1){
		if (argsL$regularisation=="E"){
			if(argsL$fixedAlpha==-1){
				for (j in 1:length(alphaslist)) {
					A[j]<-min(elasticnet[[j]]$cvm)
					}
				#Determine best alpha value from training data
				index<-which(A==min(A), arr.ind=TRUE)
				model<-elasticnet[[index]]
				}
			}else{
			model<-elasticnet
		}
	}  
	if (length(elasticnet[[1]]) > 1){ 
	#Determine error of the best alpha model on hold out data and on training data
		predict_fit<-predict(model, x_test, s="lambda.min")
		predict_fit_train<-predict(model, x_train, s="lambda.min")
		coefficients[[i]][[k]]<-coef(model, s = "lambda.min")
		pearson_correlation[[i]][k,1]<-y_test
		pearson_correlation[[i]][k,2]<-predict_fit
		spearman_correlation[[i]][k,1]<-y_test
		spearman_correlation[[i]][k,2]<-predict_fit
		test_error[[i]][k]<-(y_test-predict_fit)^2
	}else{
		coefficients[[i]][k]<-c()
		pearson_correlation[[i]][k]<-0.0
		spearman_correlation[[i]][k]<-0.0
		test_error[[i]][k]<-1.0
		validSamples[i]<-FALSE
			}
		}		     
	}
}

if (validSamples[i]){
	#Learning the model once on the full data set
	if (! is.null(argsL$seed)){
		set.seed(as.numeric(argsL$seed))
	}
	print("Learning model on entire data set for feature interpretation")
	# Partition data into test and training data sets
	Train_Data<-M

	# Split the features from response
	x_train<-as.matrix(Train_Data[,-Response_Variable_location])
	y_train<-as.vector(unlist(Train_Data[,Response_Variable_location,drop=FALSE]))
	#Creating alpha vector
 	A<-c()
	if(argsL$regularisation=="L"){
		print("The value of alpha is set to 1.0 (Lasso penalty)")
		alphaslist <- c(1.0)
	}else{
	if(argsL$regularisation=="R"){
		print("The value of alpha is set to 0.0 (Ridge penalty)")
		alphaslist <- c(0.0)
		}else{
			alphaslist<-seq(0,1,by=as.numeric(argsL$alpha))
		}
	}
	#Learning model on training data
	if(argsL$regularisation=="E"){   
		if(argsL$fixedAlpha==-1){
			if(is.null(argsL$constraint)){
				elasticnet<-mclapply(alphaslist, function(x){print(paste0("Testing alpha=",x));cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV))}, mc.cores=argsL$cores)
			}else{ 
				if(argsL$constraint=="P"){
					elasticnet<-mclapply(alphaslist, function(x){print(paste0("Testing alpha=",x));cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV))}, mc.cores=argsL$cores)
					}else{
					if(argsL$constraint=="N"){
						elasticnet<-mclapply(alphaslist, function(x){print(paste0("Testing alpha=",x));cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV))}, mc.cores=argsL$cores)
						}
					}
      		}
		}else{
		x=argsL$fixedAlpha
		if(is.null(argsL$constraint)){
			elasticnet<-cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
		}else{ 
			if(argsL$constraint=="P"){
				elasticnet<-cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
			}else{
				if(argsL$constraint=="N"){
					elasticnet<-cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
					}
				}
      	}
		}   
	}else{
	x=alphaslist[1]
	if(is.null(argsL$constraint)){
			elasticnet<-cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
		}else{ 
			if(argsL$constraint=="P"){
				elasticnet<-cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
			}else{
				if(argsL$constraint=="N"){
					elasticnet<-cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
					}
				}
		}
	}
	#Finished training the model
	if(length(elasticnet[[1]]) > 1){
		if (argsL$regularisation=="E"){
			if(argsL$fixedAlpha==-1){
				for (j in 1:length(alphaslist)) {
					A[j]<-min(elasticnet[[j]]$cvm)
				}					
				#Determine best alpha value from training data
				index<-which(A==min(A), arr.ind=TRUE)
				model<-elasticnet[[index]]
				}
			}else{
			model<-elasticnet
		}
	}
	if (length(elasticnet[[1]]) > 1){ 
		coefficientsF[[i]]<-coef(model, s = "lambda.min")
	}else{
		coefficientsF[[i]]<-c()
	}
	}	
}

if (argsL$performance == TRUE){
	if (argsL$leaveOneOutCV == FALSE){
	     for (i in 1:length(FileList)){
	          if (validSamples[i]==FALSE){
     	          next;
	          }
	          cm<-mean(unlist(pearson_correlation[[i]]),na.rm=TRUE)
	          csd<-var(unlist(pearson_correlation[[i]]),na.rm=TRUE)
	          cms<-mean(unlist(spearman_correlation[[i]]),na.rm=TRUE)
     	     csds<-var(unlist(spearman_correlation[[i]]),na.rm=TRUE)
	          erm<-mean(unlist(test_error[[i]]),na.rm=TRUE)
	          ersd<-var(unlist(test_error[[i]]),na.rm=TRUE)
	          Sample_View[[i]]<-data.frame(Sample_Name=FileList[i],Pearson=cm,PearsonVar=csd,Spearman=cms,SpearmanVar=csds,MSE=erm,MSEVar=ersd)
	     }
	}else{
	     for (i in 1:length(FileList)){
	          if (validSamples[i]==FALSE){
	               next;
     	     }
	          cm<-cor(unlist(pearson_correlation[[i]][,1]),unlist(pearson_correlation[[i]][,2]))
	          csd<-0.0
	          cms<-cor(unlist(spearman_correlation[[i]][,1]),unlist(spearman_correlation[[i]][,2]),method="spearman")
	          csds<-0.0
	          erm<-mean(unlist(test_error[[i]]),na.rm=TRUE)
	          ersd<-var(unlist(test_error[[i]]),na.rm=TRUE)
	          Sample_View[[i]]<-data.frame(Sample_Name=FileList[i],Pearson=cm,PearsonVar=csd,Spearman=cms,SpearmanVar=csds,MSE=erm,MSEVar=ersd)
			write.table(cbind(pearson_correlation[[i]][,1],pearson_correlation[[i]][,2]),paste0(argsL$outDir,"Leave_One_Out_Predicitions_",FileList[i]),quote=FALSE,sep="\t",row.names=F)
     	}
	}
	Sample_ViewF<-do.call("rbind",Sample_View)
	if (argsL$asRData==TRUE){
		save(coefficients,file=paste(argsL$outDir,"Coefficients.RData",sep=""))
	}
	write.table(Sample_ViewF,paste(argsL$outDir,"Performance_Overview.txt",sep=""),quote=FALSE,sep="\t",row.names=F)
}
fontSize=1.3
ggplotAvailable<-require("ggplot2")
gplotsAvailable<-require("gplots")

if (argsL$performance == TRUE){
	#Generating plots and interpretable output data
	for (i in 1:length(FileList)){
		if (validSamples[i]==FALSE){
			next;
		}
		#Storing Feature coefficients
		featureMatrix<-c()
		for (j in 1:length(coefficients[[i]])){
			if (length(coefficients[[i]][[j]]>1)){
				featureMatrix<-rbind(featureMatrix,coefficients[[i]][[j]][,1])
			}
		}
		if (length(featureMatrix > 1)){
			if (argsL$ftest==TRUE){
				meanFeature<-apply(featureMatrix,2,median)
				print(dim(featureMatrix))
				print(length(meanFeature))
				print(length(unlist(ftest_result[[i]])))
				featureMatrixTemp<-rbind(featureMatrix,meanFeature,c(1.0,unlist(ftest_result[[i]])))
				row.names(featureMatrixTemp)<-c(paste("Fold ",c(1:(dim(featureMatrixTemp)[1]-2))),"Median","p-value")
				write.table(featureMatrixTemp[,-1],paste(argsL$outDir,"Regression_Coefficients_",FileList[i],sep=""),quote=F,sep="\t",col.names=NA)
			}else{
				write.table(featureMatrix[,-1],paste(argsL$outDir,"Regression_Coefficients_",FileList[i],sep=""),quote=F,sep="\t",col.names=NA)
			}
			meanFeature<-apply(featureMatrix,2,median)
			featureMatrix<-featureMatrix[,-which(meanFeature==0)]
			meanFeature<-meanFeature[-which(meanFeature==0)]
			if(length(meanFeature) > 1){
				limit<-floor(min(10,length(meanFeature)/2))
				allFeatures<-cbind(featureMatrix[,order(meanFeature,decreasing=TRUE)[1:limit]],(featureMatrix[,order(meanFeature,decreasing=FALSE)[1:limit]]))
				meanFeatures<-c(meanFeature[order(meanFeature,decreasing=TRUE)[1:limit]],(meanFeature[order(meanFeature,decreasing=FALSE)[1:limit]]))
				all<-rbind(allFeatures,meanFeatures)
				all<-all[,order(all[dim(all)[1],],decreasing=TRUE)]
				row.names(all)<-c(paste("Fold ",c(1:(dim(all)[1]-1))),"Median")
				if (gplotsAvailable){
					library("gplots")
					png(paste(argsL$outDir,"Coefficients_Heatmap_",unlist(unlist(strsplit(FileList[i],".txt")))[1],".png",sep=""),width=800,height=800)
					heatmap.2(all,trace="none",col=bluered(250),srtCol=45,cexRow=fontSize,cexCol=fontSize,density.info="none",distfun = distF, dendrogram="none", margins=c(8,12),Colv=FALSE,Rowv=FALSE,key.xlab="Regression coefficient")
					dev.off()
				}
			}else{
				if (ggplotAvailable){
					library("ggplot2")
					text<-"Median regression coefficients are to small, a heatmap can not be shown."
					ggplot()+annotate("text",x=4,y=25,size=8,label=text)+theme_void()	
					ggsave(filename=paste(argsL$outDir,"Coefficients_Heatmap_",unlist(unlist(strsplit(FileList[i],".txt")))[1],".png",sep=""),width=11,height=4)
				}
			}
		}
	}
}
	


if (argsL$performance==TRUE){
	if (argsL$leaveOneOutCV==FALSE){
		for (i in 1:length(FileList)){
			if (validSamples[i]){
				#Generating Correlation boxplot
				pearson_correlation_data<-c()
				spearman_correlation_data<-c()
				for (j in 1:length(pearson_correlation[[i]])){
					pearson_correlation_data<-rbind(pearson_correlation_data,pearson_correlation[[i]][[j]])
					spearman_correlation_data<-rbind(spearman_correlation_data,spearman_correlation[[i]][[j]])
				}
				correlationMatrix<-cbind(pearson_correlation_data,spearman_correlation_data)
				colnames(correlationMatrix)<-c("Pearson Correlation","Spearman Correlation")
				write.table(correlationMatrix,file=paste(argsL$outDir,"Correlation_Matrix_",unlist(unlist(strsplit(FileList[i],".txt")))[1],".txt",sep=""),quote=F,row.names=F,sep="\t")
				png(paste(argsL$outDir,"Correlation_Boxplot",unlist(unlist(strsplit(FileList[i],".txt")))[1],".png",sep=""),width=800,height=800)
				par(cex.lab=1.5)
				par(cex.axis=1.5)
				boxplot(correlationMatrix,ylab="Value",main="Correlation on Test data")
				dev.off()
			}
		}
	}
}

for (i in 1:length(FileList)){
	if (validSamples[i]){
	#Generating perSample features
		features<-coefficientsF[[i]][which(abs(coefficientsF[[i]])!=0)]
		if (length(features > 1)){
			names(features)<-row.names(coefficientsF[[i]])[which(abs(coefficientsF[[i]])!=0)]
			coefG<-cbind(names(features),features)
			colnames(coefG)<-c("name", "value")
			row.names(coefG)<-c(1:dim(coefG)[1])
			coefG<-as.data.frame(coefG)
			write.table(coefG,paste(argsL$outDir,"Regression_Coefficients_Entire_DataSet_",FileList[i],sep=""),sep="\t",quote=F,row.names=F)
			if (ggplotAvailable){
				library("ggplot2")
				features<-coefficientsF[[i]][which(abs(coefficientsF[[i]])>0.025)]
				if (length(features>1)){
					names(features)<-row.names(coefficientsF[[i]])[which(abs(coefficientsF[[i]])>0.025)]
					coefG<-cbind(names(features),features)
					if ((round(dim(coefG)[1]*1.5) < 40) & (dim(coefG)[1] > 0)){
						colnames(coefG)<-c("name", "value")
						row.names(coefG)<-c(1:dim(coefG)[1])
						coefG<-as.data.frame(coefG)
						coefG$value<-as.numeric(as.character(coefG$value))	 	
						ggplot2::ggplot(coefG,ggplot2::aes(x=reorder(name,-value),y=value))+ggplot2::geom_bar(stat="identity")+ggplot2::theme_bw(20)+ggplot2::xlab("Feature name")+ggplot2::ylab("Regression coefficient")+ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
						ggsave(filename=paste(argsL$outDir,"Coefficients_Barplot_",unlist(unlist(strsplit(FileList[i],".txt")))[1],".png",sep=""),dpi=600,width=round(dim(coefG)[1]*1.5),limitsize=F)
					} else {
						text<-"No regression coefficients > 0.025"
						ggplot()+annotate("text",x=4,y=25,size=8,label=text)+theme_void()	
						ggsave(filename=paste(argsL$outDir,"Coefficients_Barplot_",unlist(unlist(strsplit(FileList[i],".txt")))[1],".png",sep=""),width=10,height=4)
					}
				}
			}
		}
	}
}
