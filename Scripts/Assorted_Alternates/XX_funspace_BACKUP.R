setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")
############################################################################

library(funspace)

############################################################################

source("Scripts/trait.at.time.R")
source("Scripts/innovate_elaborate.R")

############################################################################

load("Data/Amphibolurinae_Data.RData")
load("Data/Ancestral_Trait_Estimates.RData")

############################################################################

sampleLSR <- read.csv("Data/Amphibolurinae_SampleLSR.csv", header=T, row.names=1)
regimes <- read.csv("Data/Amphibolurinae_Ecology.csv", header=T)

############################################################################
# start by generating a functional space from all the samples included (450+)
sample.pca <- princomp(sampleLSR, cor=F)
plot.rotations(pca.obj = sample.pca, plot = T, which.pca = "princomp")
fs.sample <- funspace(x = sample.pca, PCs = c(2,1), n_divisions = 300)
plot(x = fs.sample,
     type = "global",
     quant.plot = T,
     globalContour = T,
     #arrows = T, arrows.length = 0.5, arrows.label.cex = 0.5, 
     pnt = T,
     colors = brewer.pal(5, "YlOrRd"))

############################################################################

# Next we'll focus on the species mean data

# isolate traits (exclude Genus, Genus_species)
LSR.anc <- allLSR[,1:19]
# combine observed and ancestral traits into a single dataframe
LSR.anc <- dplyr::bind_rows(LSR.anc, ancestors)

############################################################################

# estimate and visualize the contemporary functional space
curr.pca <- princomp(allLSR[,1:19], cor=F)
plot.rotations(pca.obj = curr.pca, plot = T, which.pca = "princomp")
fs.curr <- funspace::funspace(x = curr.pca, PCs = c(1,2), n_divisions = 300)
#summary(tsg)
plot(x = fs.curr,
     type = "global",
     quant.plot = T, quant = c(0.25,0.5,0.95,0.999),
     globalContour = T,
     arrows = T, arrows.length = 3, arrows.label.cex = 0.5, 
     pnt = T,
     threshold = 0.95, 
     colors = brewer.pal(5, "YlOrRd"))
# note, the above trait dataframe is just contemporary trait values (observed)

############################################################################

# estimate and visualize the contemporary functional space per genus
genera <- allLSR$Genus
fs.genera <- funspace::funspace(x = curr.pca, PCs = c(1,2), n_divisions = 300, group.vec = genera)
#summary(tsg)
plot(x = fs.genera,
     type = "groups",
     quant.plot = T, quant = c(0.25,0.5,0.95,0.999),
     globalContour = T,
#     arrows = T, arrows.length = 3, arrows.label.cex = 0.5, 
     pnt = T,
     threshold = 0.95, 
     colors = brewer.pal(5, "YlOrRd"))
# note, the above trait dataframe is just contemporary trait values (observed)

############################################################################

# estimate and visualize the functional space through time

# this first option works, but isn't ideal
fs.anc <- funspace.at.time(timeslices = 5, phy = agam.tree, trait.df = LSR.anc[1,19], 
                                    pcs = c(1,2), plot = T, corr = F, smooth = 100)
# note the above trait dataframe includes ancestral trait values (observed + estimated)

# this second version is preferrable
fs.anc <- funspace.through.time(timeslices = 0.5, phy = agam.tree, trait.df = LSR.anc[,1:19], 
                           pcs = c(1,2), plot = T, corr = F, smooth = 300, save.img = T, quantile = 0.999, tpd = F)

# we can also plot just the global contour to show how the adaptive landscape might have moved through time
funspace.through.time(timeslices = 0.5, phy = agam.tree, trait.df = LSR.anc[,1:19], 
                      pcs = c(1,2), plot = T, corr = F, smooth = 100, save.img = F, quantile = 0.999, tpd = F)

############################################################################

# We can plot the functional spaces as defined by PhyloEM regimes

curr.pca <- princomp(allLSR[,1:19], cor=F)
LSR.eco <- dplyr::full_join(allLSR, regimes)

fs.regimes <- funspace(x = curr.pca, PCs = c(1,2), n_divisions = 300, group.vec = LSR.eco$rpt)
plot(x = fs.regimes,
     type = "groups",
     quant.plot = T,
     quant = 0.999,
     quant.lwd = 3,
     quant.col = "grey",
     #globalContour = T,
     #globalContour.quant = 0.999,
     colors = brewer.pal(5, "YlOrRd"))


############################################################################

# We can also investigate the transitions between the different regimes

anc.regimes <- read.csv("Data/Amphibolurinae_Ecology_Ancestors.csv", header=T)

cur.pca <- princomp(LSR.anc[1:119,], cor=F)
clades <- anc.regimes$rppc[1:119]

cur.group <- funspace(x = cur.pca, PCs = c(1,2), n_divisions = 300, group.vec = clades, threshold = 0.90)
#summary(anc.group)
plot(x = cur.group,
     type = "groups",
     quant.plot = T,
     quant = 0.9,
     quant.lwd = 3,
     quant.col = "grey",
     globalContour = F,
     colors = brewer.pal(5, "YlOrRd"),
     xlim = c(-2,2), ylim = c(-1.5,1.5))


anc.pca <- princomp(LSR.anc, cor=F)
clades <- anc.regimes$rppc_anc

anc.group <- funspace(x = anc.pca, PCs = c(1,2), n_divisions = 300, group.vec = clades, threshold = 0.9)
#summary(anc.group)
plot(x = anc.group,
     type = "groups",
     quant.plot = T,
     quant = 0.9,
     quant.lwd = 3,
     quant.col = "grey",
     globalContour = F,
     colors = brewer.pal(5, "YlOrRd"),
     xlim = c(-2,2), ylim = c(-1.5,1.5))


############################################################################

require(patchwork)

# plot the transition between two nodes as lollipops
node.to.node <- function(trait.data, parent.node, child.node, palette, ret=c("plot","df")){
  new.df <- data.frame(t(trait.data[child.node,] - trait.data[parent.node,]))
  colnames(new.df) <- "value"
  new.df$trait <- rownames(new.df)
  
  out.plot <- ggplot(new.df) + 
    geom_hline(yintercept=0) +
    geom_segment(aes(x=trait, xend=trait, y=0, yend=value)) +
    geom_point(aes(x=trait, y=value, fill=value), shape=21, size=4) +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    scale_fill_distiller(palette=palette) +
    scale_y_continuous(position="right") +
    theme_classic() + theme(legend.position = "none")
  
  if(ret=="plot"){return(out.plot)}
  if(ret=="df"){return(new.df)}
}

root.to.loph <- node.to.node(trait.data=LSR.anc, parent.node=120, child.node=127, palette="Spectral")
root.to.oz   <- node.to.node(trait.data=LSR.anc, parent.node=120, child.node=131, palette="Spectral")
oz.to.moloch <- node.to.node(trait.data=LSR.anc, parent.node=129, child.node="Moloch_horridus", palette="Spectral")
oz.to.cten2  <- node.to.node(trait.data=LSR.anc, parent.node=131, child.node=133, palette="Spectral")
cten2.to.crypt  <- node.to.node(trait.data=LSR.anc, parent.node=133, child.node="Cryptagama_aurita", palette="Spectral")
cten2.to.cten1  <- node.to.node(trait.data=LSR.anc, parent.node=133, child.node=142, palette="Spectral")
oz.to.chlam  <- node.to.node(trait.data=LSR.anc, parent.node=173, child.node="Chlamydosaurus_kingii", palette="Spectral")
oz.to.pogona <- node.to.node(trait.data=LSR.anc, parent.node=180, child.node=183, palette="Spectral")
oz.to.tymp   <- node.to.node(trait.data=LSR.anc, parent.node=172, child.node=191, palette="Spectral")
oz.to.dip    <- node.to.node(trait.data=LSR.anc, parent.node=172, child.node=213, palette="Spectral")
dip.to.sup   <- node.to.node(trait.data=LSR.anc, parent.node=213, child.node="Diporiphora_superba", palette="Spectral")

# all transitions between optima
(root.to.loph /
  root.to.oz /
  oz.to.moloch /
  oz.to.cten2 /
  cten2.to.crypt /
  cten2.to.cten1 /
  oz.to.chlam /
  oz.to.pogona /
  oz.to.tymp /
  oz.to.dip /
  dip.to.sup) + plot_layout(axes="collect_x")

# selected transitions between optima
selected <- (root.to.oz /
    oz.to.moloch /
    oz.to.pogona /
    oz.to.tymp /
    cten2.to.crypt /
    dip.to.sup) + plot_layout(axes="collect_x")

plot_spacer() | plot_spacer() | selected

############################################################################

# identify the magnitude of change along each axis on each branch
node.to.node.tree <- function(trait.data, phy, palette){
  
  # build a dataframe of basic tree statistics
  ndel <- data.frame(edge = 1:nrow(phy$edge),
                     node.parent = phy$edge[,1],
                     node.child  = phy$edge[,2],
                     length = phy$edge.length,
                     name.parent = paste0("n",phy$edge[,1]))
  ndel$name.child <- sapply(ndel$node.child, function(x) ifelse(x <= Ntip(phy), phy$tip[[x]], paste0("n",x)))
  
  change.results <- NULL
  for(k in 1:nrow(ndel)){
    parent.node <- ndel[k,"node.parent"]
    child.node  <- ndel[k,"node.child"]
    change.results <- rbind(change.results, data.frame(t(node.to.node(trait.data=LSR.anc, parent.node=parent.node, child.node=child.node, palette="Spectral", ret="df")[1])))
  }
  #change.results <- abs(change.results)
  max.change <- apply(abs(change.results),1,max)
  change.results$primary <- colnames(change.results)[apply(abs(change.results),1,which.max)]
  rownames(change.results) <- NULL
  ndel <- cbind(ndel, change.results); ndel <- cbind(ndel, max.change)                                                                 
}

ndel$color <- sapply(ndel$primary, function(x) ifelse(x%in%c("Size","Tail_Length"), "lightblue","orangered"))

plot(agam.tree, edge.col=ndel$color, cex=0.3, edge.width=2); axisPhylo()

############################################################################









############################################################################

# CURRENT + ANCESTORS
clades <- NULL
for(k in 1:nrow(LSR.anc)){
  curr.taxon <- rownames(LSR.anc)[[k]]
  if(curr.taxon %in% c(extract.clade(agam.tree, node=213)$tip.label, extract.clade(agam.tree, node=213)$node.label)){clades <- append(clades, "Diporiphora"); next}
  if(curr.taxon %in% c(extract.clade(agam.tree, node=191)$tip.label, extract.clade(agam.tree, node=191)$node.label)){clades <- append(clades, "Tympanocryptis"); next}
  if(curr.taxon %in% c(extract.clade(agam.tree, node=133)$tip.label, extract.clade(agam.tree, node=133)$node.label)){clades <- append(clades, "Ctenophorus"); next}
  if(curr.taxon %in% c(extract.clade(agam.tree, node=122)$tip.label, extract.clade(agam.tree, node=122)$node.label)){clades <- append(clades, "Tree"); next}
  if(curr.taxon %in% c(extract.clade(agam.tree, node=127)$tip.label, extract.clade(agam.tree, node=127)$node.label)){clades <- append(clades, "Tree"); next}
  if(curr.taxon %in% c(extract.clade(agam.tree, node=182)$tip.label, extract.clade(agam.tree, node=182)$node.label)){clades <- append(clades, "Pogona"); next}
  if(curr.taxon %in% c(extract.clade(agam.tree, node=173)$tip.label, extract.clade(agam.tree, node=173)$node.label)){clades <- append(clades, "Lophognathus"); next}
#  if(curr.taxon == "Moloch_horridus"){clades <- append(clades, "Moloch"); next}
#  if(curr.taxon == "Chelosania_brunnea"){clades <- append(clades, "Chelosania"); next}
#  if(curr.taxon == "Physignathus_cocincinus"){clades <- append(clades, "Tree"); next}
#  if(curr.taxon == "Gowidon_longirostris"){clades <- append(clades, "Gowidon"); next}
#  if(curr.taxon == "Tropicagama_temporalis"){clades <- append(clades, "Tropicagama"); next}
#  if(curr.taxon == "Intellagama_lesueurii"){clades <- append(clades, "Tree"); next}
  else{clades <- append(clades, "Background")}
}

anc.pca <- princomp(LSR.anc, cor=F)
plot.rotations(pca.obj = anc.pca, plot = T, which.pca = "princomp")
ts.group <- funspace(x = anc.pca, PCs = c(1,3), n_divisions = 300, group.vec = clades)
#summary(tsgf)
plot(x = ts.group,
     type = "groups",
     quant.plot = T,
     globalContour = T,
     colors = brewer.pal(5, "YlOrRd"))


# CURRENT
clades <- NULL
for(k in 1:nrow(allLSR)){
  curr.taxon <- allLSR$Genus_species[[k]]
  if(curr.taxon %in% extract.clade(agam.tree, node=213)$tip.label){clades <- append(clades, "Diporiphora"); next}
  if(curr.taxon %in% extract.clade(agam.tree, node=191)$tip.label){clades <- append(clades, "Tympanocryptis"); next}
  if(curr.taxon %in% extract.clade(agam.tree, node=140)$tip.label){clades <- append(clades, "Ctenophorus"); next}
  if(curr.taxon %in% extract.clade(agam.tree, node=122)$tip.label){clades <- append(clades, "Tree"); next}
  if(curr.taxon %in% extract.clade(agam.tree, node=129)$tip.label){clades <- append(clades, "Tree"); next}
  if(curr.taxon %in% extract.clade(agam.tree, node=183)$tip.label){clades <- append(clades, "Pogona"); next}
  if(curr.taxon %in% extract.clade(agam.tree, node=173)$tip.label){clades <- append(clades, "General"); next}
  #  if(curr.taxon == "Moloch_horridus"){clades <- append(clades, "Moloch"); next}
  #  if(curr.taxon == "Chelosania_brunnea"){clades <- append(clades, "Chelosania"); next}
    if(curr.taxon == "Physignathus_cocincinus"){clades <- append(clades, "Tree"); next}
    if(curr.taxon == "Gowidon_longirostris"){clades <- append(clades, "General"); next}
    if(curr.taxon == "Tropicagama_temporalis"){clades <- append(clades, "General"); next}
  #  if(curr.taxon == "Intellagama_lesueurii"){clades <- append(clades, "Tree"); next}
    if(curr.taxon == "Cryptagama_aurita"){clades <- append(clades, "Cryptagama"); next}
    if(curr.taxon == "Diporiphora_superba"){clades <- append(clades, "Dsuperba"); next}
  else{clades <- append(clades, "Background")}
}

# estimate and visualize the contemporary functional space
curr.pca <- princomp(allLSR[,1:19], cor=F)
plot.rotations(curr.pca, plot=T, which.pca="princomp")
fs.group <- funspace(x = curr.pca, PCs = c(1,3), n_divisions = 300, group.vec = clades, threshold=0.9)
#fs.group <- funspace(x = curr.pca, PCs = c(1,3), n_divisions = 300, group.vec = allLSR$Genus)
#summary(tsg)
plot(x = fs.group,
     type = "groups",
     quant.plot = T, quant=0.9,
     globalContour = F,
     #arrows = T, arrows.length = 0.5, arrows.label.cex = 0.5, 
     pnt = T,
     colors = brewer.pal(5, "YlOrRd"))



c# CURRENT ALL SAMPLES
genera <- sapply(rownames(sampleLSR), function(x) strsplit(x," ")[[1]][1])

# estimate and visualize the contemporary functional space
fs.sample.group <- funspace(x = sample.pca, PCs = c(2,1), n_divisions = 100, group.vec = genera)
#summary(tsg)
plot(x = fs.sample.group,
     type = "groups",
     quant.plot = T, 
     globalContour = T,
     #arrows = T, arrows.length = 0.5, arrows.label.cex = 0.5, 
     #pnt = T,
     colors = brewer.pal(5, "YlOrRd"))



km <- factoextra::eclust(anc.pca$scores, "kmeans", hc_metric="euclidean", k=2)
pam <- factoextra::eclust(anc.pca$scores, "pam", k=2)


# estimate and visualize the contemporary functional space
curr.raw.pca <- princomp(sp.means, cor=F)
fs.raw.curr <- funspace(x = curr.raw.pca, PCs = c(1,2), n_divisions = 100)
#summary(tsg)
plot(x = fs.raw.curr,
     type = "global",
     quant.plot = T,
     globalContour = T,
     #arrows = T, arrows.length = 0.5, arrows.label.cex = 0.5, 
     #pnt = T,
     colors = brewer.pal(5, "YlOrRd"))


# edge 25  connects nodes 132/133

n132 <- LSR.anc[which(rownames(LSR.anc)=="n132"),]
n133 <- LSR.anc[which(rownames(LSR.anc)=="n133"),]
e25 <- abs(n132 - n133)

primary.axis <- function(phy, trait.df){
  # stop if there aren't ancestors and contemporary taxa in the dataframe
  if(!nrow(trait.df) == (Ntip(phy) + Nnode(phy))){stop("your trait dataframe does not include both tip and ancestral trait values")}
  
  edge.key <- data.frame(edge = 1:nrow(phy$edge),
                         parent.node = phy$edge[,1],
                         child.node  = phy$edge[,2])
  primaries <- NULL
  for (j in 1:nrow(edge.key)){
    pn <- edge.key[j,"parent.node"]
    cn <- edge.key[j,"child.node"]
    parent <- paste0("n",pn)
    if(cn <= Ntip(phy)){child <- phy$tip.label[[cn]]}
    else(child <- paste0("n",cn))
    
    edge.diff <- abs(
      trait.df[which(rownames(trait.df)==parent),] -
        trait.df[which(rownames(trait.df)==child),]
    )
    primaries <- append(primaries, names(edge.diff)[[which(edge.diff == max(edge.diff))]])
  }
  edge.key$primary <- primaries
  return(edge.key)
}
testo <- primary.axis(phy = agam.tree, trait.df = LSR.anc)

testa <- paintBranches(tree = agam.tree, edge = c(1,119), state = "a")
plot(testa)


paintEdges <- function (tree, edge, state, anc.state = "1") 
{
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  if (is.null(tree$maps)) 
    maps <- lapply(tree$edge.length, function(x) setNames(x, 
                                                          anc.state))
  else maps <- tree$maps
  ii <- edge
  for (i in 1:length(ii)) maps[[ii[i]]] <- setNames(tree$edge.length[[ii[i]]], 
                                                    state)
  s <- vector()
  for (i in 1:nrow(tree$edge)) s <- c(s, names(maps[[i]]))
  s <- unique(s)
  mapped.edge <- matrix(0, length(tree$edge.length), length(s), 
                        dimnames = list(edge = apply(tree$edge, 1, function(x) paste(x, 
                                                                                     collapse = ",")), state = s))
  for (i in 1:length(maps)) for (j in 1:length(maps[[i]])) mapped.edge[i, 
                                                                       names(maps[[i]])[j]] <- mapped.edge[i, names(maps[[i]])[j]] + 
    maps[[i]][j]
  tree$mapped.edge <- mapped.edge
  tree$maps <- maps
  class(tree) <- c("simmap", setdiff(class(tree), "simmap"))
  tree
}

paint.axes <- function(phy, edge.key){
  states <- unique(edge.key$primary)
  
  for(k in 1:length(states)){
    curr.state <- dplyr::filter(edge.key, primary == states[[k]])
    phy <- paintEdges(tree = phy, edge = curr.state$edge, state = curr.state[1,"primary"])
  }
  plot(phy)
  return(phy)
}
testq <- paint.axes(phy = agam.tree, edge.key = testo)

# "black"   "#DF536B"   "#61D04F"   "#2297E6"   "#28E2E5"   "#CD0BBC"   "#F5C710"

############################################################################
# start by generating a functional space from all the samples included (450+)
raw.pca <- princomp(all.raw.traits[,2:ncol(all.raw.traits)], cor=F)
plot.rotations(pca.obj = raw.pca, plot = T, which.pca = "princomp")
fs.raw <- funspace(x = raw.pca, PCs = c(1,2), n_divisions = 300)
plot(x = fs.raw,
     type = "global",
     quant.plot = T,
     globalContour = T,
     #arrows = T, arrows.length = 0.5, arrows.label.cex = 0.5, 
     pnt = T,
     colors = brewer.pal(5, "YlOrRd"))

exp(LSR.anc["n120",]*(exp(LSR.anc["n120","Size"])))
exp(LSR.anc["Intellagama_lesueurii",]*LSR.anc["Intellagama_lesueurii","Size"])


############################################################################

# function to get the euclidean distance between two points
euclidean <- function(a, b) sqrt(sum((a - b)^2))


