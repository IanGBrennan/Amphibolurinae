setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

############################################################################

library(funspace)

############################################################################

source("Scripts/Analysis_Helper_Scripts/trait.at.time.R")
source("Scripts/Assorted_Alternates/innovate_elaborate.R")

############################################################################

load("Data/Amphibolurinae_Data.RData")
load("Data/Ancestral_Trait_Estimates.RData")

############################################################################

regimes <- read.csv("Data/Amphibolurinae_Ecology.csv")
regimes <- dplyr::select(regimes, Genus_species, opt5); colnames(regimes) <- c("Genus_species","regimes")
allLSR <- left_join(allLSR, regimes)

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

# estimate and visualize the contemporary functional space per regime (k=5)
regime <- allLSR$regimes
fs.regime <- funspace::funspace(x = curr.pca, PCs = c(1,2), n_divisions = 300, group.vec = regime)
#summary(tsg)
plot(x = fs.regime,
     type = "groups",
     quant.plot = T, quant = c(0.25,0.5,0.95,0.999),
     globalContour = T,
     #     arrows = T, arrows.length = 3, arrows.label.cex = 0.5, 
     pnt = T,
     threshold = 0.95, 
     colors = brewer.pal(5, "YlOrRd"))
# note, the above trait dataframe is just contemporary trait values (observed)
# a = Ctenophorus/Cryptagama; b = Hypsilurus/Lophosaurus + arboreal allies; c = Amphibolurus + semiarboreal allies; D = Moloch; e = Tympanocryptis

############################################################################

# estimate and visualize the functional space through time
fs.anc <- funspace.through.time(timeslices = 0.5, phy = agam.tree, trait.df = LSR.anc[,1:19], 
                           pcs = c(1,2), plot = T, corr = F, smooth = 300, save.img = T, quantile = 0.999, tpd = F)

# we can also plot just the global contour to show how the adaptive landscape might have moved through time
funspace.through.time(timeslices = 0.5, phy = agam.tree, trait.df = LSR.anc[,1:19], 
                      pcs = c(1,2), plot = T, corr = F, smooth = 100, save.img = F, quantile = 0.999, tpd = F)


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

root.to.loph <- node.to.node(trait.data=LSR.anc, parent.node=120, child.node=127, palette="Spectral",ret="plot")
root.to.oz   <- node.to.node(trait.data=LSR.anc, parent.node=120, child.node=131, palette="Spectral",ret="plot")
oz.to.moloch <- node.to.node(trait.data=LSR.anc, parent.node=129, child.node="Moloch_horridus", palette="Spectral",ret="plot")
oz.to.cten2  <- node.to.node(trait.data=LSR.anc, parent.node=131, child.node=133, palette="Spectral",ret="plot")
cten2.to.crypt  <- node.to.node(trait.data=LSR.anc, parent.node=133, child.node="Cryptagama_aurita", palette="Spectral",ret="plot")
cten2.to.cten1  <- node.to.node(trait.data=LSR.anc, parent.node=133, child.node=142, palette="Spectral",ret="plot")
oz.to.chlam  <- node.to.node(trait.data=LSR.anc, parent.node=173, child.node="Chlamydosaurus_kingii", palette="Spectral",ret="plot")
oz.to.pogona <- node.to.node(trait.data=LSR.anc, parent.node=180, child.node=183, palette="Spectral",ret="plot")
oz.to.tymp   <- node.to.node(trait.data=LSR.anc, parent.node=172, child.node=191, palette="Spectral",ret="plot")
oz.to.dip    <- node.to.node(trait.data=LSR.anc, parent.node=172, child.node=213, palette="Spectral",ret="plot")
dip.to.sup   <- node.to.node(trait.data=LSR.anc, parent.node=213, child.node="Diporiphora_superba", palette="Spectral",ret="plot")

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
