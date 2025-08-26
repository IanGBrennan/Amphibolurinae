setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

#######################################################################

load("Data/Amphibolurinae_Data.RData")
load("Data/Ancestral_Trait_Estimates.RData")

#######################################################################


# if we want to show a heatmap of the euclidean distances among species
dist.mat <- as.matrix(dist(allLSR[,1:19], method="euclidean", diag=T))

# organize distance matrix rows and columns to match tree tip label order
dmr <- dist.mat[match(agam.tree$tip.label, rownames(dist.mat)), ] 
dmrc <- dmr[,match(agam.tree$tip.label, colnames(dmr))]

# plot the heatmap
heatmap(dmrc, Rowv=NA, Colv=NA, symm=T)
