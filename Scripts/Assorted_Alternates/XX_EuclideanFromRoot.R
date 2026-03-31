







source("/Users/ianbrennan/Documents/GitHub/Amphibolurinae/Scripts/fabric_Functions.R")

ll.log <- inspectBT.log(log.file = "/Applications/BayesTraitsV4.1.2/Amphibolurinae_Fabric_Traits/UpperLeg_TRAIT.txt.Log.txt",
                         single.trait = T, skip.lines = 284, model = "fabric", ESS.threshold = 200)



ul.res <- fabricProcess(fpp.path = "/Applications/BayesTraitsV4.1.2/FabricPostProcessor/Amphibolurinae/UpperLeg_TRAIT.txt.VarRates.txt.csv",
                        phy = agam.tree, trait.name = "Upper Leg")
fabricPlot(process.obj = ul.res, phy = agam.tree, trait.name = "Upper Leg")

ll.res <- fabricProcess(fpp.path = "/Applications/BayesTraitsV4.1.2/FabricPostProcessor/Amphibolurinae/LowerLeg_TRAIT.txt.VarRates.txt.csv",
                        phy = agam.tree, trait.name = "LowerLeg")
fabricPlot(process.obj = ll.res, phy = agam.tree, trait.name = "Lower Leg")

ft.res <- fabricProcess(fpp.path = "/Applications/BayesTraitsV4.1.2/FabricPostProcessor/Amphibolurinae/Foot_TRAIT.txt.VarRates.txt.csv",
                        phy = agam.tree, trait.name = "Foot")
fabricPlot(process.obj = ft.res, phy = agam.tree, trait.name = "Foot")

ua.res <- fabricProcess(fpp.path = "/Applications/BayesTraitsV4.1.2/FabricPostProcessor/Amphibolurinae/UpperArm_TRAIT.txt.VarRates.txt.csv",
                        phy = agam.tree, trait.name = "Upper Leg")
fabricPlot(process.obj = ua.res, phy = agam.tree, trait.name = "Upper Arm")

la.res <- fabricProcess(fpp.path = "/Applications/BayesTraitsV4.1.2/FabricPostProcessor/Amphibolurinae/LowerArm_TRAIT.txt.VarRates.txt.csv",
                        phy = agam.tree, trait.name = "LowerLeg")
fabricPlot(process.obj = la.res, phy = agam.tree, trait.name = "Lower Arm")

ha.res <- fabricProcess(fpp.path = "/Applications/BayesTraitsV4.1.2/FabricPostProcessor/Amphibolurinae/Hand_TRAIT.txt.VarRates.txt.csv",
                        phy = agam.tree, trait.name = "Hand")
fabricPlot(process.obj = ha.res, phy = agam.tree, trait.name = "Hand")




tail.res <- fabricProcess(fpp.path = "/Applications/BayesTraitsV4.1.2/FabricPostProcessor/Amphibolurinae/TailLength_TRAIT.txt.VarRates.txt.csv",
                          phy = agam.tree, trait.name = "TailLength")
fabricPlot(process.obj = tail.res, phy = agam.tree, trait.name = "Tail Length")

size.res <- fabricProcess(fpp.path = "/Applications/BayesTraitsV4.1.2/FabricPostProcessor/Amphibolurinae/Size_TRAIT.txt.VarRates.txt.csv",
                          phy = agam.tree, trait.name = "Size")
fabricPlot(process.obj = size.res, phy = agam.tree, trait.name = "Size")

hw.res <- fabricProcess(fpp.path = "/Applications/BayesTraitsV4.1.2/FabricPostProcessor/Amphibolurinae/HeadWidth_TRAIT.txt.VarRates.txt.csv",
                          phy = agam.tree, trait.name = "Head Width")
fabricPlot(process.obj = hw.res, phy = agam.tree, trait.name = "Head Width")


fab.files <- dir("/Applications/BayesTraitsV4.1.2/FabricPostProcessor/Amphibolurinae", ".csv", full.names = T) 
trait.names <- sapply(sapply(fab.files, function(x) strex::str_after_last(x, "/")), function(y) strex::str_before_first(y, "_"))
par(mfrow=c(5,4))
full.shifts <- NULL
for (j in 1:length(fab.files)){
  cur.res <- fabricProcess(fpp.path = fab.files[[j]],
                           phy = agam.tree, trait.name = trait.names[[j]])
  full.shifts <- rbind(full.shifts, cur.res)
  fabricPlot(process.obj = cur.res, phy = agam.tree, trait.name = trait.names[[j]], font.size = 0.2, line.width =  1)
}

fabricPlot(process.obj = full.shifts, phy = agam.tree, trait.name = "All Traits", font.size = 0.3)

all.b <- dplyr::filter(full.shifts, parameter == "b")
nrow(dplyr::filter(all.b, scale < 0))
nrow(dplyr::filter(all.b, scale > 0))

all.v <- dplyr::filter(full.shifts, parameter == "v")
nrow(dplyr::filter(all.v, scale < 1))
nrow(dplyr::filter(all.v, scale > 1))



# satisfying explanation of difference between estimating multivariate
# distance via euclidean and mahalanobis measures
# https://pjbartlein.github.io/GeogDataAnalysis/lec18.html

# mahalanobis distances from log-shape ratios
mdist <- mahalanobis(x=LSR.anc[1:119,], center=unlist(LSR.anc[120,]), cov=cov(LSR.anc[1:119,]))
mah.df <- data.frame(taxon = names(mdist), distance = mdist)
View(mah.df[order(mah.df$distance),])

# mahalanobis distances from PCA of log-shape ratios
pca.res <- princomp(LSR.anc[1:120,], cor=F)
mdist.pca <- mahalanobis(x=pca.res$scores[1:119,], center=unlist(pca.res$scores[120,]), cov=cov(pca.res$scores[1:119,]))
mah.pca.df <- data.frame(taxon = names(mdist.pca), distance = mdist.pca)
View(mah.pca.df[order(mah.pca.df$distance),])

# euclidean distances from log-shape ratios
edist <- apply(LSR.anc[1:119,], 1, function(x) euclidean(x, LSR.anc[120,]))
ozdist <- apply(LSR.anc[1:119,], 1, function(x) euclidean(x, LSR.anc[129,]))
euc.df <- data.frame(taxon = names(edist), distance = edist, OzMRCA = -ozdist)
#View(euc.df[order(euc.df$distance),])
euc.df$genus <- sapply(euc.df$taxon, function(x) strsplit(x,"_")[[1]][1])

# median distance
#euc.df$distance <- euc.df$distance - median(euc.df$distance)
# rescale to mean
#euc.df$distance <- euc.df$distance - min(euc.df$distance)

ggplot(euc.df) +
  geom_jitter(aes(x=distance, y=genus, fill=genus), width=0.25, shape=21, size=3) +
  scale_y_discrete(limits=rev) +
  xlab("Euclidean Multivariate Distance from Amphibolurine TMRCA") +
  theme_bw() + theme(legend.position="none")

ggplot(euc.df) +
  #geom_jitter(aes(x=distance, y=genus, fill=genus), width=0.25, shape=21, size=3) +
  geom_segment(aes(x=0, xend=distance, y=taxon, yend=taxon,color=genus)) +
  geom_point(aes(x=distance, y=taxon, fill=genus), shape=21, size=3) +
#  geom_segment(aes(x=0, xend=OzMRCA, y=taxon, yend=taxon,color=genus)) +
#  geom_point(aes(x=OzMRCA, y=taxon, fill=genus), shape=21, size=3) +
  scale_y_discrete(limits=rev) +
  geom_vline(xintercept=0, linetype="dotted") +
  #xlim(0,3) +
  xlab("Euclidean Multivariate Distance from Amphibolurine MRCA") +
  theme_classic() + theme(legend.position="none")


table(all.b$trait)
table(all.v$trait)

pathDistance <- function(phy, trait.df){
  
  node.edge <- data.frame(edge = 1:nrow(phy$edge),
                          parent.node = phy$edge[,1],
                          child.node  = phy$edge[,2],
                          parent      = paste0("n",phy$edge[,1]))
  child <- NULL
  for(k in 1:nrow(node.edge)){
    if(node.edge[k,"child.node"] > Ntip(phy)){child <- append(child, paste0("n",node.edge[k,"child.node"]))}
    else if(node.edge[k,"child.node"] <= Ntip(phy)){child <- append(child, phy$tip.label[node.edge[k,"child.node"]])}
  }
  node.edge$child <- child
  
  node.edge$euc.dist <- apply(node.edge, 1, function(x) euclidean(LSR.anc[x["parent"],], LSR.anc[x["child"],]))
  
  phy$edge.length <- node.edge$euc.dist
  
  plot.phylo(phy, cex=0.5, type="fan")
}





############################################################################

# Next we'll focus on the species mean data

# isolate traits (exclude Genus, Genus_species)
LSR.anc <- allLSR[,1:19]
# combine observed and ancestral traits into a single dataframe
LSR.anc <- dplyr::bind_rows(LSR.anc, ancestors)

############################################################################





