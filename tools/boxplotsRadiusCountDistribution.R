library(ggplot2)
library(scales)


data <- read.table('~/Documents/TEPIC/Code/2016-05-24-15-33_ResolutionPerWindowToSize.txt', header=TRUE, stringsAsFactors=FALSE)
data$reso<-as.factor(data$reso)
data$window<-as.factor(data$window)

# Define text sizes
ts <- 20
texs <- 18
titles <- 22
ltexs <- 16

ggplot(data, aes(window, count)) +
  geom_boxplot(aes(fill=reso,position="dodge")) +
  #ggtitle('Histogram on number of loops in loop-windows of different size\nfor K562') +
  xlab('Loop-window size') +
  ylab('Number of loops') +
  theme_bw() +
  scale_fill_discrete(name="Hi-C Resolution") +
  scale_x_discrete(labels=c("1kb", "2kb", "3kb", "4kb", "5kb", "10kb", "15kb", "25kb", "50kb", "100kb", "250kb", "500kb", "1mb")) +
  theme(axis.title = element_text(size=ts), plot.title = element_text(size=titles), axis.text = element_text(size=texs), legend.title=element_text(size=texs), legend.text=element_text(size=ltexs), legend.position="bottom")
