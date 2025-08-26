# Code below incorporates 3 extant Hypsilurus species for which
# we have morphological trait data.
# Phylogenetic placement and divergence times for these taxa
# are based on results from Title et al. (2024)


setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")
t1 <- read.tree("Trees/Amphibolurinae_MCMCTree.tre")
t1 <- ladderize(t1)
t1$edge.length <- t1$edge.length*100
plot(t1, cex=0.3); 
t1 <- drop.tip(t1, c("Sphenodon_punctatus", "Phrynocephalus_putjatai", "Diporiphora_linga2", "Pogona_minor_minor2", "Moloch_sp."))
plot(t1, cex=0.4); axisPhylo()

t2 <- bind.tip(tree = t1,
               tip.label = "Hypsilurus_godeffroyi", 
               where = getMRCA(t1, c("Hypsilurus_modestus", "Hypsilurus_magnus")),
               #where = 126,
               position = 20.35 - (max(nodeHeights(t1)) - nodeheight(t1, getMRCA(t1, c("Hypsilurus_modestus", "Hypsilurus_magnus")))),
               edge.length = 20.35)
plot(t2, cex=0.3); axisPhylo()

t3 <- bind.tip(tree = t2,
               tip.label = "Hypsilurus_papuensis", 
               where = which(t2$tip.label=="Hypsilurus_magnus"),
               position = 8.78,
               edge.length = 8.78)
plot(t3, cex=0.3); axisPhylo()

t4 <- bind.tip(tree = t3,
               tip.label = "Hypsilurus_nigrigularis", 
               where = which(t3$tip.label=="Hypsilurus_magnus"),
               position = 7.01,
               edge.length = 7.01)
plot(t4, cex=0.3); axisPhylo()

t4 <- ladderize(t4)
t4 <- ladderize(t4, right=F)
plot(t4, cex=0.3); axisPhylo()
write.tree(t4, "/Users/ianbrennan/Documents/GitHub/Amphibolurinae/Trees/Amphibolurinae_forR.tre")

