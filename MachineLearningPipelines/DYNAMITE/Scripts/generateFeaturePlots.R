#This script can be used to generate figures to interpret the results of a DYNAMITE analysis.
#It includes a density plot of TF affinities for all differentially expressed genes, a scatter plot
#of the affinity values coloured according to the expression changes of the genes, a miniature heatmap
#showing the regression coefficient of the selected TF as well as a density and a scatter plot on a reduced
#version of the dataset (0.9 quantile), that is more stable against outliers.
#
#
#arg1: Output directory specified in the DYNAMITE configuration file
#arg2: Name of the feature for which the analysis plots should be generated
#arg3: Identifier that should be used for group1
#arg4: Identifier that should be used for group2

library("reshape2")
library("ggplot2")
library("gridExtra")
args<-commandArgs(T)
if (length(args) < 4){
cat("Arguments:\n")
cat("arg1: Output directory specified in the DYNAMITE configuration file\n")
cat("arg2: Name of the feature for which the analysis plots should be generated\n")
cat("arg3: Identifier that should be used for group1\n")
cat("arg4: Identifier that should be used for group2\n")
q(save="no")
}

print("Reading affinities for group 1")
group1<-read.table(paste0(args[1],"/Affinities/Mean/Mean_Affinities_group1.txt"),header=T,stringsAsFactors = F)
print("Reading affinities for group 2")
group2<-read.table(paste0(args[1],"/Affinities/Mean/Mean_Affinities_group2.txt"),header=T,stringsAsFactors = F)
print("Reading expression data")
expression<-read.table(paste0(args[1],"/IntegratedData/Log2/Integrated_Data_Log2_Quotient.txt"),header=T,stringsAsFactors = F)
print("Reading coefficients")
regCoef<-read.table(paste0(args[1],"/Learning_Results/Regression_Coefficients_Cross_Validation_Integrated_Data_For_Classification.txt"),header=T,stringsAsFactors = F)
feature=args[2]
group1L=args[3]
group2L=args[4]

fontsize=8.5
if ((feature != "Peak_Counts") & (feature !="Peak_Length")){
  prefix<-"Affinity of TF "
}else{
  prefix<-"Value of "
}

commonIDs=intersect(intersect(group1$geneID,group2$geneID),expression$geneID)
g1M<-match(commonIDs,group1$geneID)
g2M<-match(commonIDs,group2$geneID)
index=which(colnames(group1)==feature)
selection=cbind(group1[g1M,index],group2[g2M,index])
selection<-as.data.frame(selection)
colnames(selection)=c(group1L,group2L)
mSelection<-melt(selection)
colnames(mSelection)<-c("Type","value")

#Density plot
pValue<-ks.test(selection[,1],selection[,2])$p.value
den<-density(mSelection$value)
l<-paste0("p-value of a KS test: ",signif(pValue,5))
plot1<-ggplot2::ggplot(mSelection,aes(value,fill=Type,guide="legend"))+
  geom_density(adjust=.8,alpha=0.2)+
  scale_fill_manual(values = c("black","green"))+
  theme_bw(fontsize)+
  labs(fill="Cell type")+
  ggtitle("Full data set")+
  annotate("text",x=max(mSelection$value)*.6,y=max(den$y)*.9,label=l,size=0.35*fontsize)+
  xlab(paste0(prefix,feature))+
  theme(legend.position="none")

#Density plot without outliers
outliers<-union(which(selection[,1] > quantile(selection[,1],probs=0.9)),which(selection[,2] > quantile(selection[,2],probs=0.9)))
selectionOR<-selection[-outliers,]
mSelectionOR<-melt(selectionOR)
colnames(mSelectionOR)<-c("Type","value")
plot3<-ggplot2::ggplot(mSelectionOR,aes(value,fill=Type,guide="legend"))+
  geom_density(adjust=.8,alpha=0.2)+
  scale_fill_manual(values = c("black","green"))+
  theme_bw(fontsize)+
  labs(fill="Cell type")+
  ggtitle("Reduced data set (0.9 quantile)")+
  xlab(paste0(prefix,feature))

#Scatter plot 
selectionC<-cbind(selection,expression$Expression[match(commonIDs,expression$geneID)])
colnames(selectionC)<-c(group1L,group2L,"log2(FoldChange)")
plot2<-ggplot2::ggplot(selectionC,aes(x=selectionC[,1],y=selectionC[,2],colour=selectionC[,3]))+
      geom_point()+
      xlab(paste0(prefix,feature," in ",group1L))+
      ylab(paste0(prefix,feature," in ",group2L))+
      theme_bw(fontsize)+
      geom_abline()+
      labs(colour=paste0("log2(\nexp(",group1L,")/\n         exp(",group2L,"))"))+
      ggtitle(" ")+
      scale_colour_gradient2(low = "blue",high="red",mid="white",midpoint=0)+
      geom_density2d()+
      scale_x_continuous(limits=c(min(selectionC[,1],selectionC[,2]),max(selectionC[,1],selectionC[,2])))+
      scale_y_continuous(limits=c(min(selectionC[,1],selectionC[,2]),max(selectionC[,1],selectionC[,2])))+
      theme(legend.position="none")

#Scatter plot without outliers
selectionCOR<-selectionC[-outliers,]
plot4<-ggplot2::ggplot(selectionCOR,aes(x=selectionCOR[,1],y=selectionCOR[,2],colour=selectionCOR[,3]))+
  geom_point()+
  xlab(paste0(prefix,feature," in ",group1L))+
  ylab(paste0(prefix,feature," in ",group2L))+
  theme_bw(fontsize)+
  geom_abline()+
  labs(colour=paste0("log2(\nexp(",group1L,")/\nexp(",group2L,"))"))+
  ggtitle(" ")+
  scale_colour_gradient2(low = "blue",high="red",mid="white",midpoint=0)+
  geom_density2d()+
  scale_x_continuous(limits=c(min(selectionCOR[,1],selectionCOR[,2]),max(selectionCOR[,1],selectionCOR[,2])))+
  scale_y_continuous(limits=c(min(selectionCOR[,1],selectionCOR[,2]),max(selectionCOR[,1],selectionCOR[,2])))


#Geom tile plot showing the regression coefficient
coefs<-as.data.frame(cbind(rep((feature),dim(regCoef)[1]+1),c(regCoef[,which(colnames(regCoef)==feature)],median(regCoef[,which(colnames(regCoef)==feature)])),c(paste0("Fold ",seq(1:dim(regCoef)[1])),"Median")))
colnames(coefs)<-c("Feature","Value","Fold")
coefs$Value<-as.numeric(as.character(coefs$Value))
coefs$Fold<-factor(coefs$Fold,levels=rev(coefs$Fold))
plot5<-ggplot2::ggplot(coefs,aes(x=Feature,y=Fold,fill=Value))+
        geom_tile(colour="black")+
        theme_bw(fontsize)+
        xlab("")+
        scale_fill_gradient2(high="red",low="blue",mid="white",midpoint=0)+
        labs(fill="Regression\n coefficient")
  
svg(paste0(args[1],"/Feature_overview_",group1L,"vs",group2L,"_",feature,".svg"),width=10,height=6)
grid.arrange(plot1, plot2, plot3, plot4, plot5,layout_matrix=rbind(c(1,1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3,3,3,3,5,5,5,5,5),
                                                                   c(1,1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3,3,3,3,5,5,5,5,5),
                                                                   c(2,2,2,2,2,2,2,2,2,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5),
                                                                   c(2,2,2,2,2,2,2,2,2,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5)))
dev.off()

