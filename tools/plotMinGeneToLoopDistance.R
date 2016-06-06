#Import libs and set global options
library(ggplot2)
library(plyr)
library(scales)
options(scipen=999)
xlabels <- c(-1, "1kb", "5kb", "10kb", "15kb", "20kb", "25kb", "30kb", "35kb", "40kb", "45kb", "50kb", "75kb", "100kb", "200kb")

# Read raw data
data5kb <- read.table('~/Documents/TEPIC/Code/2016-06-06-13-51_MinDistances_Res5000.txt', header=TRUE, stringsAsFactors=FALSE)
data10kb <- read.table('~/Documents/TEPIC/Code/2016-06-06-13-51_MinDistances_Res10000.txt', header=TRUE, stringsAsFactors=FALSE)
data25kb <- read.table('~/Documents/TEPIC/Code/2016-06-06-13-52_MinDistances_Res25000.txt', header=TRUE, stringsAsFactors=FALSE)

# Transform data, prepare it for plotting
data <- rbind(data5kb, data10kb, data25kb)
data[,2] <- apply(subset(data, select=c("distance")), 2, cut, c(-Inf,-1, 1000, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 75000, 100000, 200000,Inf), labels=xlabels<-c(xlabels, Inf))
#data[data=="Inf"] <- "100kb+" # activate this to include value greater than the last interval cut, and disable the next line!
data <- data[data$distance != Inf , ]
data <- data[data$distance != -1, ]
uni <- subset(data, !duplicated(geneID))
data <- data[,c("distance", "resolution")]
uni <- uni[,c("distance", "resolution")]

# Set labels properly
xlabelsreduced <- head(xlabels, -1)
xlabelsreduced <- tail(xlabelsreduced, -1)
data$distance <- factor(data$distance, levels=xlabelsreduced)
uni$distance <- factor(uni$distance, levels=xlabelsreduced)

# Apply sums and merge data in one dataframe
countdata <- ddply(data, .(distance,resolution), nrow)
countunidata <- ddply(uni, .(distance,resolution), nrow)
temp <- ddply(countunidata, .(distance), summarize, V1 = sum(V1))
temp$resolution <- "All"
cdata <- rbind(countdata, temp)
cdata <- cdata[order(cdata$distance,cdata$resolution),]


# only barplot
ggplot(cdata, aes(distance), ordered=TRUE) +
  geom_bar(aes(y=cdata$V1, fill=factor(resolution)), colour="black", position = "dodge", stat="identity") +
  scale_fill_discrete(name="Hi-C Resolutions") +
  ggtitle('Binned minimum window size around TSS containing nearest loop-edge in K562 cell-line\nHistogram') +
  xlab('Window size') +
  ylab('#genes') +
  theme_bw()

#only line and point plot
ggplot(cdata, aes(distance, cumsum(cdata$V1)), ordered=TRUE) +
  geom_line(aes(group=factor(resolution), colour=factor(resolution))) + 
  geom_point(colour="darkgrey") +
  scale_colour_discrete(name = "Hi-C Resolutions") +
  ggtitle('Binned minimum window size around TSS containing nearest loop-edge in K562 cell-line\nCumulative histogram') +
  xlab('Window size') +
  ylab('#genes') +
  theme_bw()

# combined plot
ggplot(cdata, aes(distance, y=cumsum(cdata$V1)), ordered=TRUE) +
  geom_bar(aes(y=cdata$V1, fill=factor(resolution)), colour="black", position = "dodge", stat="identity") +
  geom_line(aes(group=factor(resolution), colour=factor(resolution))) + 
  geom_point(colour="darkgrey") +
  scale_fill_discrete(name="Hi-C Resolutions") +
  scale_colour_discrete(name = "Hi-C Resolutions\ncumulative #genes") +
  ggtitle('Binned minimum window size around TSS containing nearest loop-edge in K562 cell-line\nCombined histogram') +
  xlab('Window size') +
  ylab('#genes') +
  theme_bw()
