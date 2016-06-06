library(ggplot2)
library(scales)


data <- read.table('~/Documents/TEPIC/Code/2016-05-24-15-33_ResolutionPerWindowToSize.txt', header=TRUE, stringsAsFactors=FALSE)
data$reso<-as.factor(data$reso)
data$window<-as.factor(data$window)

ggplot(data, aes(window, count)) +
  geom_boxplot(aes(fill=reso,position="dodge")) +
  ggtitle('Loops in TSS-windows according to different Hi-C resolutions') +
  xlab('TSS window-size') +
  ylab('Number of loops') +
  theme_bw() +
  scale_fill_discrete(name="Resolution") +
  scale_x_discrete(labels=c("1kb", "2kb", "3kb", "4kb", "5kb", "10kb", "15kb", "25kb", "50kb", "100kb", "250kb", "500kb", "1mb"))

#coveragegrob <- grobTree(textGrob(paste("Total coverage: ",overall$Coverage," %"), x=0.77,  y=0.93, hjust=0,
                                  #gp=gpar(fontsize=12)))

#ggplot(st, aes(x=as.factor(st$window), y=st$count, group=st$window)) + geom_boxplot(aes(fill=st$resolution)) + xlab('TSS Window-size') + ylab('Number of loops')
  #ggtitle('test') + 
  #scale_fill_discrete(name="Classification", breaks=c("Remaining", "Count"), labels=c("Not in window", "In window")) + 
  #annotation_custom(coveragegrob) + geom_text(data=cst, aes(y=stcopy[,c("Count")], label=paste(cst$value, "%")), color="black", size=3, position = position_dodge(0.9), vjust = -0.25)
