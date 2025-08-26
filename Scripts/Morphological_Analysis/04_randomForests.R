# set the working directory to the GitHub repo
setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

############################################################################

# Load necessary packages
library(phytools)
library(randomForest)
library(ggRandomForests)
library(randomForestExplainer)
# useful information here: https://htmlpreview.github.io/?https://github.com/geneticsMiNIng/BlackBoxOpener/blob/master/randomForestExplainer/inst/doc/randomForestExplainer.html

# helper functions for plotting
source("Scripts/innovate_elaborate.R")

# load the agamid data
load("Data/Amphibolurinae_Data.RData")
load("Data/Ancestral_Trait_Estimates.RData")

############################################################################

# read in the ecology data
eco.df <- read.csv("Data/Amphibolurinae_Ecology.csv", h=T)
eco.df <- dplyr::select(eco.df, Genus_species, opt6, opt11)

# bind the ecology data to the trait shape data
rfdata <- dplyr::full_join(allLSR, eco.df)
rownames(rfdata) <- rfdata$Genus_species

############################################################################

# choose the data we're interested in
rfdata.eco <- rfdata[,c(1:19,22)]

# add in the regime states estimated from PhyloEM
#rfdata.reg <- rfdata[,c(1:19,28)] # opt6
rfdata.reg <- rfdata[,c(1:19,29)] # opt11
colnames(rfdata.reg)[[20]] <- "Regime"
rfdata.reg$Regime <- as.factor(rfdata.reg$Regime)

############################################################################

# set a seed to make this reproducible
set.seed(2024)

# run the RandomForest and visualize it
# forest <- randomForest(Ecology ~ ., data = rfdata.eco, ntree = 10000, localImp = TRUE)
forest <- randomForest(Regime ~ ., data = rfdata.reg, ntree = 10000, localImp = TRUE)
plot(forest)

# extract the categorizations for each taxon, map to regime designation
rf.votes <- data.frame(forest$votes)
rf.votes$Regime <- rfdata.reg$Regime

# now use the existing categorizations to estimate regimes for ancestral taxa
rf.predictions <- predict(forest, ancestors, type="prob")
rfp <- predict(forest, ancestors, type="prob")
rownames(rfp) <- 120:(120+(nrow(rfp)-1))

# match up the regimes iwth the data
reg <- setNames(rfdata.reg$Regime, rownames(rfdata.reg))


############################################################################

# Visualize the resulting output as pie charts at nodes

# we'll use the phytools architecture to plot the RF likelihoods for ancestral states
#fit.STP.full <- fitMk(agam.tree, reg, model=STP.full, pi="fitzjohn")
fit.ER <- fitMk(agam.tree, reg, model="ER", pi="fitzjohn")


# extract the marginal ancestral state estimates
anc <- ancr(fit.ER, type="marginal")
# replace the ancr ace matrix with our RF predictions
anc$ace <- rfp


# Plot the marginal ancestral states on the phylogeny

# set colors
c1 <- colorRampPalette(RColorBrewer::brewer.pal(9,"Spectral"))
cols <- c1(length(unique(reg)))

# set the pies to be bigger if there is greater uncertainty (no state with >50%)
node.cex<-apply(anc$ace,1,function(x) if(any(x>0.5)) 0.3 else 0.6)


# plot tree
plot(anc,
     args.plotTree=list(type="arc", arc_height=0.5, fsize=0.5, offset=3, color="grey"),
     args.nodelabels=list(cex=node.cex,piecol=cols),
     args.tiplabels=list(cex=0.2,piecol=cols),
     legend=FALSE)

# plot legend
legend(x=-70, y=0,
       levels(reg),pch=16,col=cols,
       horiz=T,cex=0.8,bty="n",pt.cex=2,
       x.intersp=0.5)




























# TESTING OUT LDA (VERY DIFFERENT RESULTS!)
lda.test <- MASS::lda(Ecology~., rfdata)
ggord::ggord(lda.test, lda.test$Genus)
plot.rotations(pca.obj = lda.test, which.pca = "lda")

lda.sp <- MASS::lda(Genus_species~.,sampleLSR[,c(1:20)])
ggord::ggord(lda.sp, lda.sp$Genus)
plot.rotations(pca.obj = lda.sp, which.pca = "lda")
lda.sp$scores <- data.frame(predict(lda.sp)$x)
lda.sp$scores$Genus <- sampleLSR$Genus
lda.sp$scores$Genus_species <- sampleLSR$Genus_species

ggplot() +
  geom_point(data = transform(lda.sp$scores, Genus = NULL), 
             colour = "grey85", aes(x=LD1, y=LD2), size=2) +
  #geom_polygon(data=x.hull, aes(x=PC1, y=PC2, fill=genus), alpha=0.5) +
  geom_point(data=lda.sp$scores, aes(x=LD1, y=LD2, fill=Genus), size=3, pch=21) +
  theme_classic() + theme(legend.position = "none") + facet_wrap(~Genus)

lda.tpd <- TPD::TPDs(species = rep("Amphibolurinae",119), lda.sp$scores[,1:2], alpha=0.999)
TPD::plotTPD(TPD = lda.tpd)
fun.lda <- funspace(x=lda.tpd, PCs=c(1,2), n_divisions=100)



testo <- princomp(allLSR[,1:19])


rfdata.anc <- trait.eco[,c(1:19,21)]
colnames(rfdata.anc)[[20]] <- "Ecology"



min_depth_frame <- min_depth_distribution(forest)
head(min_depth_frame, n = 10)
plot_min_depth_distribution(min_depth_frame)

p1 <- predict(forest, rfdata)
caret::confusionMatrix(p1, rfdata$Ecology)

p2 <- predict(forest, rfdata.anc)
caret::confusionMatrix(p2, rfdata.anc$Ecology)


data<-iris

set.seed(222)
ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.7, 0.3))
train <- data[ind==1,]
test <- data[ind==2,]


plotTree(agam.tree,
     args.plotTree=list(type="arc", arc_height=0.5, fsize=0.5, offset=3, color="grey"),
     #     args.plotTree=list(type="fan",part=0.75, color="grey", fsize=0.5,offset=3),
     args.nodelabels=list(piecol=cols),
     args.tiplabels=list(cex=0.2,piecol=cols),
     legend=FALSE)

plotTree(agam.tree, type="arc")

