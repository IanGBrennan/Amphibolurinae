library(ggplot2)

setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

astats <- read.csv("Data/Agamidae_SampleInfo.csv")

ingroup <- dplyr::filter(astats, Source == "AusARG_Agamidae")

testo <- tidyr::pivot_longer(data=ingroup, cols=c("AHE","gene","uce"))

ggplot(testo, aes(fill=name, y=value, x=lineage)) + 
  geom_bar(position="stack", stat="identity") + scale_x_discrete(limits=rev) + 
  coord_flip() + theme_classic()
