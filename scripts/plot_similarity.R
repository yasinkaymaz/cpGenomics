library(tidyverse)
library(reshape2)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

data1 <- read.delim(args[1], header = TRUE)
data <- data.frame(Pos=data1$start, Sim2Ref=data1$Sim2Ref)
mdata <- melt(data, id=c("Pos"))
plot.1 <- ggplot(mdata, aes(Pos/1000, y=value, colour=variable))+
  geom_line(size=1,show.legend = FALSE)+
  theme(axis.text.x = element_text(angle = 0, hjust = .5,colour = "black", vjust=0.5))+
  xlab("Genomic Position (kb)")+
  ylab("% Dissimilarity")+
  ylim(0.0, 2.0)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf(args[2], width = 9, height = 5)
plot.1 
dev.off()
