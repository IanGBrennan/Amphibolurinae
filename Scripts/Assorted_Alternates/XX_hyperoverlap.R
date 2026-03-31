
# set the working directory to the GitHub repo
setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

############################################################################

# Load necessary packages
library(phytools)

# load the agamid data
load("Data/Amphibolurinae_Data.RData")
load("Data/Ancestral_Trait_Estimates.RData")

############################################################################

eco.df <- read.csv("Data/Amphibolurinae_Ecology.csv", h=T)


# I should try this process with the PCs

LSR.samples <- read.csv("Data/Amphibolurinae_SampleLSR.csv")
LSR.samples$Genus_species <- sapply(LSR.samples$Genus_species, function(x) stringr::str_replace(x," ","_"))
LSR.samples <- LSR.samples[which(LSR.samples$Genus_species %in% eco.df$Genus_species),]
LSR.samples$Regime <- unlist(sapply(LSR.samples$Genus_species, function(x) eco.df[which(eco.df$Genus_species == x),"rpt"]))
LSR.samples$Ecology <- unlist(sapply(LSR.samples$Genus_species, function(x) eco.df[which(eco.df$Genus_species == x),"Ecology"]))


samples.regime <- hyperoverlap::hyperoverlap_set(LSR.samples[,2:20], as.factor(LSR.samples$Regime))
ho.plot(samples.regime)
ho.plot(dplyr::filter(samples.regime, entity1=="Amph")) # you can subset states for plotting


ho.plot <- function (x, cols = pal) 
{
  x2 <- x
  x2$entity2 <- x$entity1; x2$entity1 <- x$entity2
  x <- rbind(x, x2)
  
  entity1 <- entity2 <- NULL
  pal = c("red", "blue", "lightgrey")
  overlap.plot = ggplot2::ggplot(data = x, aes(entity1, entity2, 
                                               fill = as.factor(x$result))) + geom_tile(color = "white", 
                                                                                        na.rm = TRUE, lwd = 3) + theme_void() + scale_fill_manual(values = cols) + 
    scale_x_discrete(position = "top") + theme(axis.text.x = element_text(angle = 90, 
                                                                          vjust = 0, size = 8, hjust = 0), axis.text.y = element_text(vjust = 0, 
                                                                                                                                      size = 8, hjust = 0)) + coord_fixed()
  overlap.plot$labels$fill <- "Result"
  return(overlap.plot)
}

# A LOGICAL STEP WOULD BE TO TAKE N SAMPLES OF THE ANCESTRAL NODES AND
# THEN DO HYPEROVERLAP AGAINST THE KNOWN REGIMES TO DETERMINE WHICH THEY
# OVERLAP WITH




