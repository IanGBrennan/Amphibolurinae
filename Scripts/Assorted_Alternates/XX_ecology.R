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
rownames(eco.df) <- eco.df$Genus_species
eco.df <- eco.df[,c(1:6,10:13)]
eco.df$regime_l1ou_PC <- sapply(eco.df$regime_l1ou_PC, function(x) letters[x])
eco.df$regime_PhyloEM_PC <- sapply(eco.df$regime_PhyloEM_PC, function(x) letters[x])
eco.df$regime_PhyloEM_Trait <- sapply(eco.df$regime_PhyloEM_Trait, function(x) letters[x])
eco.df$regime_PhyloEM_opt <- sapply(eco.df$regime_PhyloEM_opt, function(x) letters[x])

############################################################################

pca.res <- data.frame(prcomp(allLSR[,1:19])$x)
pca.res$Genus_species <- rownames(pca.res)
eco.res <- dplyr::full_join(pca.res, eco.df, by="Genus_species")

ggplot(eco.res) +
  geom_point(data = transform(eco.res, regime_l1ou_PC = NULL), aes(x=PC1, y=PC2),color="lightgrey", size=2) +
  geom_point(aes(x=PC1, y=PC2, color=regime_l1ou_PC), size=3) +
  theme_bw() + theme(legend.position="none") + facet_wrap(~regime_l1ou_PC)
ggplot(eco.res) +
  geom_point(data = transform(eco.res, regime_PhyloEM_PC = NULL), aes(x=PC1, y=PC2),color="lightgrey", size=2) +
  geom_point(aes(x=PC1, y=PC2, color=regime_PhyloEM_PC), size=3) +
  theme_bw() + theme(legend.position="none") + facet_wrap(~regime_PhyloEM_PC)
ggplot(eco.res) +
  geom_point(data = transform(eco.res, regime_PhyloEM_raw = NULL), aes(x=PC1, y=PC2),color="lightgrey", size=2) +
  geom_point(aes(x=PC1, y=PC2, color=regime_PhyloEM_raw), size=3) +
  theme_bw() + theme(legend.position="none") + facet_wrap(~regime_PhyloEM_raw)
ggplot(eco.res) +
  geom_point(data = transform(eco.res, Ecology = NULL), aes(x=PC1, y=PC2),color="lightgrey", size=2) +
  geom_point(aes(x=PC1, y=PC2, color=Ecology), size=3) +
  theme_bw() + theme(legend.position="none") + facet_wrap(~Ecology)

ggplot(eco.res) +
  geom_point(data = transform(eco.res, rlp = NULL), aes(x=PC1, y=PC2),color="lightgrey", size=2) +
  geom_point(aes(x=PC1, y=PC2, color=rlp), size=3) +
  theme_bw() + theme(legend.position="none") + facet_wrap(~rlp)

# PhyloEM K=8 PC1/PC2
ggplot(eco.res) +
  geom_point(data = transform(eco.res, rpt = NULL), aes(x=PC1, y=PC2),color="lightgrey", size=2) +
  geom_point(aes(x=PC1, y=PC2, color=rpt), size=3) +
  theme_bw() + theme(legend.position="none") + facet_wrap(~rpt)

# PhyloEM K=4 PC1/PC2
ggplot(eco.res) +
  geom_point(data = transform(eco.res, rpo = NULL), aes(x=PC1, y=PC2),color="lightgrey", size=2) +
  geom_point(aes(x=PC1, y=PC2, color=rpo), size=3) +
  theme_bw() + theme(legend.position="none") + facet_wrap(~rpo)

# PhyloEM K=4 PC1/PC3
ggplot(eco.res) +
  geom_point(data = transform(eco.res, rpo = NULL), aes(x=PC1, y=PC3),color="lightgrey", size=2) +
  geom_point(aes(x=PC1, y=PC3, color=rpo), size=3) +
  theme_bw() + theme(legend.position="none") + facet_wrap(~rpo)

############################################################################

library(ggplot2)
library(ggsankey)

eco.res$Genus <- sapply(eco.res$Genus_species, function(x) strsplit(x, "_")[[1]][1])

# Old Genera --> Old Names --> New Genera
optima <- eco.res %>% make_long(regime_l1ou_PC, Ecology, regime_PhyloEM_PC)
optima <- eco.res %>% make_long(Ecology, regime_l1ou_PC)
optima <- eco.res %>% make_long(Ecology, regime_PhyloEM_PC)
optima <- eco.res %>% make_long(Ecology, Genus, regime_l1ou_PC)
optima <- eco.res %>% make_long(Ecology, Genus, rlp)
optima <- eco.res %>% make_long(Ecology, rlp)
optima <- eco.res %>% make_long(Genus, rpo, rppc)
optima <- eco.res %>% make_long(Genus, k5, k12)
optima <- eco.res %>% make_long(k5, Genus, k12)


#frogs.df$node <- factor(frogs.df$node, levels = rev(levels(factor(frogs.df$node))))
#frogs.df$next_node < - factor(frogs.df$next_node, levels = rev(levels(factor(frogs.df$next_node))))
#frogs.df$
#frogs$genus <- factor(frogs$genus, levels = c("Cyclorana","Litoria","Nyctimystes"))
ggplot(optima,aes(x = x, 
                    next_x = next_x,
                    node = node, 
                    next_node = next_node, 
                    #color = node,
                    fill = node,
                    label = node)) +
  geom_sankey(show.legend = F, linetype = "solid", flow.alpha=0.75, node.color=1) + 
  #geom_sankey_label(show.legend = F) +
  geom_sankey_label(size = 3.5, color = 1, fill = "white", show.legend=F) +
  theme_sankey(base_size = 7) +
  scale_fill_viridis_d(alpha = 1)


############################################################################


eco.df <- dplyr::select(eco.df, ecology, genus_species)
#eco.df <- dplyr::select(eco.df, ecology)

# isolate traits (exclude Genus, Genus_species)
LSR.anc <- allLSR[,1:19]
# combine observed and ancestral traits into a single dataframe
LSR.anc <- dplyr::bind_rows(LSR.anc, ancestors)
LSR.anc$genus_species <- rownames(LSR.anc)

trait.eco <- full_join(LSR.anc, eco.df)
trait.eco$ecology[is.na(trait.eco$ecology)] <- "ancestors"
rownames(trait.eco) <- rownames(LSR.anc)

te.pca <- prcomp(trait.eco[,1:19])
te.pca.x <- data.frame(te.pca$x)
te.pca.x$ecology <- trait.eco$ecology

pca.l1ou <- prcomp(allLSR[,1:19])$x


cols <- setNames(c(hcl.colors(8)[c(1,3,5,7)],"grey"),
                 c("arboreal","saxicolous","semiarboreal","terrestrial","ancestors"))

ggplot(te.pca.x) +
  geom_point(aes(x=PC1, y=PC2, color=ecology), size=3) +
  scale_color_manual(values = cols) +
  theme_classic()
ggplot(te.pca.x) +
  geom_point(aes(x=PC1, y=PC3, color=ecology), size=3) +
  scale_color_manual(values = cols) +
  theme_classic()
ggplot(te.pca.x) +
  geom_point(aes(x=PC2, y=PC3, color=ecology), size=3) +
  scale_color_manual(values = cols) +
  theme_classic()

ggplot() +
  geom_point(data = transform(te.pca.x, ecology = NULL), 
             color = "lightGrey", aes(x=PC1, y=PC2)) +
  geom_point(data=te.pca.x, aes(x=PC1, y=PC2, 
                             color=ecology), alpha=1) +
  scale_color_manual(values=cols) +
  theme_classic() + theme(legend.position = "none") + facet_wrap(~ecology)


library(hyperoverlap)
test1 <- iris[which(iris$Species!="versicolor"),c(1:3,5)]
setosa_virginica3d <- hyperoverlap_detect(test1[,1:3], test1$Species)
#hyperoverlap_plot(setosa_virginica3d)

ho.test <- hyperoverlap_set(te.pca.x[,1:3],te.pca.x$ecology)

avt <- dplyr::filter(te.pca.x, ecology %in% c("arboreal", "terrestrial"))
res.avt <- hyperoverlap_detect(avt[,1:3], avt$ecology)
res.avt.n <- res.avt <- hyperoverlap_detect(avt[,1:19], avt$ecology)
hyperoverlap_lda(res.avt.n)

avs <- dplyr::filter(te.pca.x, ecology %in% c("arboreal", "semiarboreal"))
res.avs <- hyperoverlap_detect(avs[,1:3], avs$ecology)
res.avs.n <- hyperoverlap_detect(avs[,1:19], avs$ecology)

tvr <- dplyr::filter(te.pca.x, ecology %in% c("terrestrial", "saxicolous"))
res.tvr <- hyperoverlap_detect(tvr[,1:3], tvr$ecology)
res.tvr@result
res.tvr.n <- hyperoverlap_detect(tvr[,1:19], tvr$ecology)
res.tvr.n@result

allLSR$Ecology <- eco.df$ecology

ArTe <- dplyr::filter(allLSR, Ecology %in% c("arboreal", "terrestrial"))
res.ArTe <- hyperoverlap_detect(ArTe[,1:19], ArTe$Ecology)
res.ArTe@result
#hyperoverlap_lda(res.ArTe)

ArSe <- dplyr::filter(allLSR, Ecology %in% c("arboreal", "semiarboreal"))
res.ArSe <- hyperoverlap_detect(ArSe[,1:19], ArSe$Ecology)
res.ArSe@result
#hyperoverlap_lda(res.ArTe)

TeSa <- dplyr::filter(allLSR, Ecology %in% c("terrestrial", "saxicolous"))
res.TeSa <- hyperoverlap_detect(TeSa[,1:19], TeSa$Ecology)
res.TeSa@result
#hyperoverlap_lda(res.ArTe)


test2 <- dplyr::filter(te.pca.x, !ecology == "ancestors")
all <- hyperoverlap_set(test2[,1:3], test2$ecology)

test5 <- iris
all_spp <- hyperoverlap_set(test5[,1:4],test5$Species)
hyperoverlap_pairs_plot(all_spp)





plot.treePath(phy=agam.tree, focus="clade", tip.spread=c("Diporiphora_adductus","Diporiphora_vescus"),
              lwd=5, plot.type="phylogram", font.size=0.5)



