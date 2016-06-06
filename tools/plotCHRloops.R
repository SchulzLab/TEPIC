library(ggplot2)
library(reshape2)
library(scales)
library(grid)


st <- read.table('~/Documents/TEPIC/Code/2016-05-23-16-47_3kb.txt', header=TRUE, col.names=c("ChrNo", "Counts", "Remaining", "Total", "Coverage"))
overall <- tail(st, 1)
stcopy <- head(st,-1)
stcopy <- stcopy[(stcopy$ChrNo !="X" & stcopy$ChrNo !='Y'),]
stcopy$ChrNo <- as.numeric(stcopy$ChrNo)
stcopy <- stcopy[order(stcopy$ChrNo),]
stcopy <- rbind(stcopy, st[(st$ChrNo =="X"),])
stcopy <- rbind(stcopy, st[(st$ChrNo =="Y"),])
stcopy$ChrNo <- factor(stcopy$ChrNo, levels=stcopy$ChrNo, ordered=TRUE)

mst <- melt(stcopy[,c('ChrNo', 'Counts', 'Remaining')], id.vars = "ChrNo")
cst <- melt(stcopy[,c('ChrNo', 'Coverage')], id.vars = "ChrNo")

coveragegrob <- grobTree(textGrob(paste("Total coverage: ",overall$Coverage," %"), x=0.77,  y=0.93, hjust=0,
                            gp=gpar(fontsize=12)))

ggplot(mst, aes(ChrNo, value, fill=variable), ordered=TRUE) + geom_bar( stat="identity" ) + xlab('Chromosome') + ylab('Number of loops') +
  ggtitle('Loop distribution over Chromosomes in the K562 cell-line for a 3kb TSS-window') + 
  scale_fill_discrete(name="Classification", breaks=c("Remaining", "Counts"), labels=c("Not in window", "In window")) + 
  annotation_custom(coveragegrob) + geom_text(data=cst, aes(y=stcopy[,c("Counts")], label=paste(cst$value, "%")), color="black", size=3, position = position_dodge(0.9), vjust = -0.25) + theme_bw()
