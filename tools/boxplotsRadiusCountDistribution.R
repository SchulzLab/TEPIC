library(ggplot2)
library(scales)


data <- read.table('~/Documents/TEPIC/Code/2016-05-24-15-33_ResolutionPerWindowToSize.txt', header=TRUE, stringsAsFactors=FALSE)
data$reso<-as.factor(data$reso)
data$window<-as.factor(data$window)

ggplot(data, aes(window, count)) +
  geom_boxplot(aes(fill=reso,position="dodge")) +
  ggtitle('Number of loops in loop-windows according to different Hi-C resolutions') +
  xlab('Loop-window') +
  ylab('Number of loops') +
  theme_bw() +
  scale_fill_discrete(name="Hi-C Resolution") +
  scale_x_discrete(labels=c("1kb", "2kb", "3kb", "4kb", "5kb", "10kb", "15kb", "25kb", "50kb", "100kb", "250kb", "500kb", "1mb"))
