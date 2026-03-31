library(ggplot2)

setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

############################################################################

# required scripts
source("/Users/ianbrennan/Documents/GitHub/Amphibolurinae/Scripts/fabric_Functions.R")
source("/Users/ianbrennan/Documents/GitHub/Amphibolurinae/Scripts/trait.at.time.R")


############################################################################

# load the data
load("Data/Amphibolurinae_Data.RData")
load("Data/Ancestral_Trait_Estimates.RData")

############################################################################

# combine the observed and estimated trait values
LSR.anc <- allLSR[,1:19]
# combine observed and ancestral traits into a single dataframe
LSR.anc <- dplyr::bind_rows(LSR.anc, ancestors)

############################################################################


# practice
pca.lsr <- prcomp(LSR.anc[,1:19])$x


lineage.density <- function(phy, trait, tip.spread){
  # build a dataframe of basic tree statistics
  ndel <- data.frame(edge = 1:nrow(phy$edge),
                     node.parent = phy$edge[,1],
                     node.child  = phy$edge[,2],
                     length = phy$edge.length,
                     timestart = sapply(phy$edge[,1], function(x) nodeheight(phy, x))-max(nodeHeights(phy)),
                     timestop  = sapply(phy$edge[,2], function(x) nodeheight(phy, x))-max(nodeHeights(phy)),
                     parent.name = paste0("n",phy$edge[,1]))
  ndel$child.name <- sapply(ndel$node.child, function(x) ifelse(x <= Ntip(phy), phy$tip[[x]], paste0("n",x)))
  
  # get the most recent common node of our focal tips
  mrcn <- getMRCA(phy, tip.spread)
  # extract all the descendant edges from our mrcn
  edge.list <- getEdges(phy, node=mrcn)
  # get all the descendant tips from our mrcn
  desc.all <- getDescendants(phy, mrcn)
  desc.tip <- desc.all[which(desc.all<=Ntip(phy))]
  desc.tip.names <- phy$tip.label[desc.tip]
  
  
  # determine the multivariate euclidean distance along each edge (between pairs of parent/child nodes)
  edist <- NULL
  for (k in 1:nrow(ndel)){edist <- append(edist, euclidean(trait[ndel[k,"node.parent"],], trait[ndel[k,"node.child"],]))}
  ndel$edist <- edist
  
  # subset the tree data down to the focal edges
  focal.ndel <- ndel[ndel$edge%in%edge.list,]
  
  # estimate the euclidean distance travelled across all the focal edges
  raw.L <- sum(focal.ndel$edist)
  L <- pracma::nthroot(raw.L,ncol(trait))
  
  # estimate the volume of the focal trait data
  trait.hv <- hypervolume::hypervolume(subset(trait, rownames(trait) %in% desc.tip.names))
  V <- trait.hv@Volume
  
  # calculate the lineage density
  D <- L/V
  
  # report stats
  return(list("Raw L" = raw.L, 
              "L" = L,
              "V" = V,
              "D" = D,
              "number of edges" = length(edge.list),
              "number of tips" = length(desc.tip)))
}

D.cten <- lineage.density(phy = agam.tree,
                trait = pca.lsr[,1:3],
                tip.spread = c("Ctenophorus_clayi", "Ctenophorus_fordi"))
D.tymp <- lineage.density(phy = agam.tree,
                          trait = pca.lsr[,1:3],
                          tip.spread = c("Tympanocryptis_einasleighensis", "Tympanocryptis_lineata"))
D.pogo <- lineage.density(phy = agam.tree,
                          trait = pca.lsr[,1:3],
                          tip.spread = c("Rankinia_diemensis", "Pogona_barbata"))
D.dipo <- lineage.density(phy = agam.tree,
                          trait = pca.lsr[,1:3],
                          tip.spread = c("Diporiphora_superba", "Diporiphora_sobria"))
D.amph <- lineage.density(phy = agam.tree,
                          trait = pca.lsr[,1:3],
                          tip.spread = c("Chlamydosaurus_kingii", "Amphibolurus_muricatus"))
D.hyps <- lineage.density(phy = agam.tree,
                          trait = pca.lsr[,1:3],
                          tip.spread = c("Hypsilurus_godeffroyi", "Hypsilurus_magnus"))
D.aust <- lineage.density(phy = agam.tree,
                          trait = pca.lsr[,1:3],
                          tip.spread = c("Moloch_horridus", "Diporiphora_sobria"))

res.df <- data.frame(t(data.frame(Ctenophorus = t(data.frame(D.cten)),
                     Tympanocryptis = t(data.frame(D.tymp)),
                     Pogona = t(data.frame(D.pogo)),
                     Diporiphora = t(data.frame(D.dipo)),
                     Amphibolurus = t(data.frame(D.amph)),
                     Hypsilurus = t(data.frame(D.hyps)),
                     Australia = t(data.frame(D.aust))
                     )))
res.df$group <- rownames(res.df)

ggplot(res.df) +
  geom_point(aes(x=L, y=V, color=group, size=D)) +
  theme_classic()




################################################################################

# if we wanted to compare to null values from under BM

rate.pc1 <- geiger::fitContinuous(agam.tree, pca.lsr[1:119,1], model="BM")$opt$sigsq
rate.pc2 <- geiger::fitContinuous(agam.tree, pca.lsr[1:119,2], model="BM")$opt$sigsq
rate.pc3 <- geiger::fitContinuous(agam.tree, pca.lsr[1:119,3], model="BM")$opt$sigsq

sim.pc1 <- data.frame(fastBM(agam.tree, sig2=rate.pc1, nsim=10, internal=T))
sim.pc2 <- data.frame(fastBM(agam.tree, sig2=rate.pc2, nsim=10, internal=T))
sim.pc3 <- data.frame(fastBM(agam.tree, sig2=rate.pc3, nsim=10, internal=T))


ld.sims <- NULL
for(k in 1:10){
  curr.df <- data.frame(PC1=sim.pc1[,k], PC2=sim.pc2[,k], PC3=sim.pc3[,k])
  rownames(curr.df) <- rownames(LSR.anc)
  
  curr.ld <- lineage.density(phy = agam.tree,
                             trait = curr.df,
                             tip.spread = c("Diporiphora_superba", "Diporiphora_sobria"))
                             #tip.spread = c("Tympanocryptis_einasleighensis", "Tympanocryptis_lineata"))
                             #tip.spread = c("Rankinia_diemensis","Pogona_barbata"))
                             #tip.spread = c("Chlamydosaurus_kingii","Amphibolurus_muricatus"))
  ld.sims <- rbind(ld.sims, data.frame(curr.ld, sim=k))
}
D.tymp <- lineage.density(phy = agam.tree,
                          trait = pca.lsr[,1:3],
                          tip.spread = c("Tympanocryptis_einasleighensis", "Tympanocryptis_lineata"))

D.tymp$sim <- "EMP"
D.tymp.df <- data.frame(D.tymp)
ld.sims <- rbind(ld.sims, D.tymp.df)
ld.sims$color <- c(rep("black",10),"red")
plot(ld.sims$V ~ ld.sims$L)

t.tree <- extract.clade(agam.tree, node=191)
t.data <- pca.lsr[97:119,1:3]
t.sig2 <- geiger::fitContinuous(t.tree, t.data)



################################################################################

# NEXT STEP WOULD BE TO BUILD A METHOD THAT ESTIMATES THE LINEAGE DENSITY
# MEASURE FROM NON-CLADES, E.G. BY ECOLOGY, WHICH COULD BE GENERALIZED TO
# WORK ON A WHOLE TREE, INCLUDING ANCESTRAL EDGES

################################################################################

lineage.density.sim <- function(phy, trait, tip.spread, nsim){
  
  # get the mrca of focal tips and extract the subtree
  mrca <- getMRCA(phy, tip.spread)
  new.tree <- extract.clade(phy, mrca)
  
  # downsample the trait data to just the focal taxa
  new.trait <- trait[which(rownames(trait)%in%new.tree$tip.label),]
  
  # fit some BM models and extract the estimated rate
  rate.pc1 <- geiger::fitContinuous(new.tree, new.trait[,1], model="BM")$opt$sigsq
  rate.pc2 <- geiger::fitContinuous(new.tree, new.trait[,2], model="BM")$opt$sigsq
  rate.pc3 <- geiger::fitContinuous(new.tree, new.trait[,3], model="BM")$opt$sigsq
  
  # simulate some data under the estimated parameters
  sim.pc1 <- data.frame(fastBM(new.tree, sig2=rate.pc1, nsim=10, internal=T))
  sim.pc2 <- data.frame(fastBM(new.tree, sig2=rate.pc2, nsim=10, internal=T))
  sim.pc3 <- data.frame(fastBM(new.tree, sig2=rate.pc3, nsim=10, internal=T))
  
  # build a dataframe of basic tree statistics
  ndel <- data.frame(edge = 1:nrow(new.tree$edge),
                     node.parent = new.tree$edge[,1],
                     node.child  = new.tree$edge[,2],
                     length = new.tree$edge.length,
                     timestart = sapply(new.tree$edge[,1], function(x) nodeheight(new.tree, x))-max(nodeHeights(new.tree)),
                     timestop  = sapply(new.tree$edge[,2], function(x) nodeheight(new.tree, x))-max(nodeHeights(new.tree)),
                     parent.name = paste0("n",new.tree$edge[,1]))
  ndel$child.name <- sapply(ndel$node.child, function(x) ifelse(x <= Ntip(new.tree), new.tree$tip[[x]], paste0("n",x)))
  
  # get the most recent common node of our focal tips
  mrcn <- getMRCA(new.tree, tip.spread)
  # extract all the descendant edges from our mrcn
  edge.list <- getEdges(new.tree, node=mrcn)
  # get all the descendant tips from our mrcn
  desc.all <- getDescendants(new.tree, mrcn)
  desc.tip <- desc.all[which(desc.all<=Ntip(new.tree))]
  desc.tip.names <- new.tree$tip.label[desc.tip]
  
  
  # determine the multivariate euclidean distance along each edge (between pairs of parent/child nodes)
  edist <- NULL
  for (k in 1:nrow(ndel)){edist <- append(edist, euclidean(trait[ndel[k,"node.parent"],], trait[ndel[k,"node.child"],]))}
  ndel$edist <- edist
  
  # subset the tree data down to the focal edges
  focal.ndel <- ndel[ndel$edge%in%edge.list,]
  
  # estimate the euclidean distance travelled across all the focal edges
  raw.L <- sum(focal.ndel$edist)
  L <- pracma::nthroot(raw.L,ncol(trait))
  
  # estimate the volume of the focal trait data
  trait.hv <- hypervolume::hypervolume(subset(trait, rownames(trait) %in% desc.tip.names))
  V <- trait.hv@Volume
  
  # calculate the lineage density
  D <- L/V
  
  # report stats
  return(list("Raw L" = raw.L, 
              "L" = L,
              "V" = V,
              "D" = D,
              "number of edges" = length(edge.list),
              "number of tips" = length(desc.tip)))
}


############################################################################




