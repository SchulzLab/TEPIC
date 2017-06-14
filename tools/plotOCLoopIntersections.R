library(ggplot2)
library(reshape2)

df <- read.table('~/Downloads/OC_HiC_Intersections.txt', header=TRUE, stringsAsFactors=FALSE)

dfcopy = df
dfcopy$Coverage = dfcopy$Intersections/dfcopy$Loops

df$Loops = df$Loops - df$Intersections
df <- melt(df, id.var="Cell.line")
dfcopymelt <- melt(dfcopy[,c('Cell.line', 'Coverage')], id.vars = "Cell.line")

# Define text sizes
ts <- 20
texs <- 18
titles <- 22
ltexs <- 18

ggplot(df, aes(Cell.line, value, fill=variable), ordered=TRUE) + geom_bar( stat="identity" ) + xlab('Cell line') + ylab('Number of loops') +
  scale_fill_discrete(name="Classification", breaks=c("Loops", "Intersections"), labels=c("Total number of loops", "Number of loops intersecting with open chromatin peaks")) +
  geom_text(data=dfcopymelt, aes(y=dfcopy$Intersections, label=paste(round(dfcopymelt$value*100, 0), "%"), fill=NA), color="black", size=6, position = position_dodge(0.9), vjust = -0.25) +
  theme_bw() + theme(axis.title = element_text(size=ts), axis.text = element_text(size=texs), legend.title=element_text(size=ts), legend.text=element_text(size=ltexs), legend.position="bottom")