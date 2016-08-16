#Import libs and set global options
library(ggplot2)
library(plyr)
library(data.table)
library(scales)
options(scipen=999)
xlabels <- c(-1, "1kb", "5kb", "10kb", "15kb", "20kb", "25kb", "30kb", "35kb", "40kb", "45kb", "50kb")

# Read raw data
data5kb <- read.table('~/Documents/TEPIC/Code/2016-06-06-13-51_MinDistances_Res5000.txt', header=TRUE, stringsAsFactors=FALSE)
data10kb <- read.table('~/Documents/TEPIC/Code/2016-06-06-13-51_MinDistances_Res10000.txt', header=TRUE, stringsAsFactors=FALSE)
data25kb <- read.table('~/Documents/TEPIC/Code/2016-06-06-13-52_MinDistances_Res25000.txt', header=TRUE, stringsAsFactors=FALSE)

# Transform data, prepare it for plotting
dat <- rbind(data5kb, data10kb, data25kb)
dat[,2] <- apply(subset(dat, select=c("distance")), 2, cut, c(-Inf,-1, 1000, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000,Inf), labels=xlabels<-c(xlabels, Inf))
#data[data=="Inf"] <- "50kb+" # activate this to include value greater than the last interval cut, and disable the next line!
dat <- dat[dat$distance != Inf , ]
dat <- dat[dat$distance != -1, ]
uni <- subset(dat, !duplicated(geneID))
dat <- dat[,c("distance", "resolution")]
uni <- uni[,c("distance", "resolution")]

# Set labels properly and apply ordering
xlabelsreduced <- head(xlabels, -1)
xlabelsreduced <- tail(xlabelsreduced, -1)
dat$distance <- factor(dat$distance, levels=xlabelsreduced)
uni$distance <- factor(uni$distance, levels=xlabelsreduced)

# Count rows and apply cumsums for distinct resolutions
countdata <- ddply(dat, .(distance,resolution), nrow)
dt <- data.table(countdata)
dt <- dt[,list(distance, V1, cumsum = cumsum(V1)),by=list(resolution)]
dt <- as.data.frame.matrix(dt)

# Count rows and apply cumsums for all resolutions combined
countunidata <- ddply(uni, .(distance,resolution), nrow)
countunidata <- ddply(countunidata, .(distance), summarize, V1 = sum(V1))
temp <- data.table(countunidata)
temp <- temp[,list(distance, V1, cumsum = cumsum(V1))]
temp <- as.data.frame.matrix(temp)
temp$resolution <- "All"

# Merge dataframes and apply ordering
cdata <- rbind(dt, temp)
cdata$distance <- factor(cdata$distance, levels=xlabelsreduced)

# only barplot
ggplot(cdata, aes(distance), ordered=TRUE) +
  geom_bar(aes(y=cdata$V1, fill=factor(resolution)), colour="black", position = "dodge", stat="identity") +
  scale_fill_discrete(name="Hi-C Resolution") +
  ggtitle('Binned minimum window size around TSS containing nearest loop-edge in K562 cell-line\nHistogram') +
  xlab('Loop-window') +
  ylab('Number of genes') +
  theme_bw()

#only line and point plot
ggplot(cdata, aes(distance, y=cdata$cumsum), ordered=TRUE) +
  geom_line(aes(group=factor(resolution, levels=c("5000","10000","25000","All")), colour=factor(resolution,levels=c("5000","10000","25000","All"))),size=1.5) + 
  geom_point(colour="darkgrey") +
  scale_colour_discrete(name = "Hi-C Resolution") +
  ggtitle('Binned minimum window size around TSS containing nearest loop-edge in K562 cell-line\nCumulative histogram') +
  xlab('Loop-window') +
  ylab('Number of genes') +
  theme_bw()

# combined plot
ggplot(cdata, aes(distance, y=cdata$cumsum), ordered=TRUE) +
  geom_bar(aes(y=cdata$V1, fill=factor(resolution)), colour="black", position = "dodge", stat="identity") +
  geom_line(aes(group=factor(resolution), colour=factor(resolution))) + 
  geom_point(colour="darkgrey") +
  scale_fill_discrete(name="Hi-C Resolution") +
  scale_colour_discrete(name = "Hi-C Resolution\ncumulative #genes") +
  ggtitle('Binned minimum window size around TSS containing nearest loop-edge in K562 cell-line\nCombined histogram') +
  xlab('Loop-window') +
  ylab('Number of genes') +
  theme_bw()
