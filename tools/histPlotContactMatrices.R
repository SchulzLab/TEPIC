library(ggplot2)
library(grid)
library(plyr)
library(scales)

options(scipen=999)

df <- read.table('~/Downloads/chr1_5kb.RAWobserved', header=FALSE, stringsAsFactors=FALSE, col.names=c("x", "y", "observations")) # 748 und 25 sowie 550 und 10
#df <- read.table('~/Downloads/chr2_5kb.RAWobserved', header=FALSE, stringsAsFactors=FALSE, col.names=c("x", "y", "observations")) # 400 10
dfred <- df[df$observations < 550 & df$x != df$y,]
#dfred <- df[df$x != df$y,]

gr <- grobTree(textGrob('binwidth: 10', x=0.77, y=0.93, hjust=0))

#dfred <- df[df$observations<400 & df$observations>0,]
#histinfo <-hist(log(dfred$observations), col="red", freq=TRUE, breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))
ggplot(dfred, aes(x = observations)) +
  geom_histogram(binwidth=10) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  ggtitle('Histogram of observed contacts for Chr.1 in K562') +
  annotation_custom(gr)

#quantile(df$observations, c(.99999952))
test <- df[df$observations>=550,]
