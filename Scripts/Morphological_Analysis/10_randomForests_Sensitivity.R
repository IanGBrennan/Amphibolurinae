setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

#######################################################################

library(phytools)
library(randomForest)
source()

load("Data/Amphibolurinae_Data.RData")
load("Data/Ancestral_Trait_Estimates.RData")

############################################################################

# read in the ecology data
eco.df <- read.csv("Data/Amphibolurinae_Ecology.csv", h=T)
eco.df <- dplyr::select(eco.df, Genus_species, opt5)

# bind the ecology data to the trait shape data
rfdata <- dplyr::full_join(allLSR, eco.df)
rownames(rfdata) <- rfdata$Genus_species

############################################################################

# choose the data we're interested in
rfdata.eco <- rfdata[,c(1:19,22)]

# add in the regime states estimated from PhyloEM
rfdata.reg <- rfdata[,c(1:19,22)] # opt5
#rfdata.reg <- rfdata[,c(1:19,23)] # opt11
colnames(rfdata.reg)[[20]] <- "Regime"
rfdata.reg$Regime <- as.factor(rfdata.reg$Regime)

# remind ourselves of the number and names of the traits
traits <- colnames(allLSR)[1:19]

############################################################################

# We need to get rate and trait estimates from the empirical data

# turn our tree/data into a simmap object
agam.simmap <- make.simmap(agam.tree, setNames(rfdata.reg$Regime, rownames(rfdata.reg)))
# fit a multi OU model to each trait
oum.fits <- lapply(traits, function(x) mvOU(agam.simmap, dplyr::select(allLSR,x), model="OUM"))
names(oum.fits) <- traits

############################################################################

# Now we're ready to simulate and fit data

# choose the number of simulations we'd like
nsim <- 500

# generate some random taxon names
new.names <- stringi::stri_rand_strings(119, 8, '[A-Z]')

# rename the tips of the tree
new.tree <- agam.tree; new.tree$tip.label <- new.names

# build a reasonable equal-rates Q matrix for transitions among regimes
q <- 0.01 # rate
Q <- q*matrix(c(
  -1,1,1,1,1,
  1,-1,1,1,1,
  1,1,-1,1,1,
  1,1,1,-1,1,
  1,1,1,1,-1),5,5,
  dimnames=list(c("a","b","c","d","e"),
                c("a","b","c","d","e")))

# Make a loop that will simulate data, fit a random forest model
# to that simulated data, estimate ancestral trait values
# and fit those ancestors to the modern regimes. This should give
# us some idea of the sensitivity/accuracy of the RandomForest approach
# and tell us if there is a bias in inferring a given state

# create a loop of length nsim to get results
rfp.res <- NULL; sim.res <- NULL; anc.res <- NULL; sim.trees <- NULL
for (p in 1:nsim){
  
  # simulate a discrete trait history on the tree under the established rate matrix
  tt <- sim.history(new.tree,Q)
  sim.trees[[p]] <- tt
  
  # keep the information about what ancestral node states are
  edge.states <- data.frame(cbind(tt$edge, tt$node.states)[,c(1,3)])
  edge.states <- edge.states[!duplicated(edge.states),]
  colnames(edge.states) <- c("node","state")
  rownames(edge.states) <- paste0("n",edge.states$node)
  # we want to compare this to our RandomForest later
  anc.res[[p]] <- edge.states
  
  # visualize the distribution of new traits
#  cols<-setNames(hcl.colors(n=5),c("a","b","c","d","e"))
#  plot(tt,cols,ftype="off",type="arc",
#       arc_height=0.2,lwd=3)
#  legend("bottomright",c("a","b","c","d","e"),lwd=3,
#         col=cols,bty="n",horiz=TRUE)
  
  # extract the states for each tip
  regimes <- tt$states
  
  # let's simulate data based on empirical fits 
  sim.data <- data.frame(sapply(traits, function(x) mvSIM(tt, nsim=1, model="OUM", param=oum.fits[[x]])))
  rownames(sim.data) <- tt$tip.label
  sim.data$Regime <- as.factor(tt$states)
  sim.res[[p]] <- sim.data
  
  # estimate the internal node values
  anc.data <- data.frame(sapply(traits, function(x) fastAnc(tt, setNames(sim.data[,x],rownames(sim.data)))))
  
  # run the RandomForest classifier
  forest.new <- randomForest(Regime ~ ., data = sim.data, ntree = 10000, localImp = TRUE)
  
  # now use the existing categorizations to estimate regimes for ancestral taxa
  rfp <- predict(forest.new, anc.data, type="prob")
  rownames(rfp) <- 120:(120+(nrow(rfp)-1))
  
  # add the results to a lists we can inspect later
  rfp.res[[p]] <- rfp
}

############################################################################

# Summarize the results by comparing the simulated and estimated fits

# extract the estimated states at each node
# if all state estimates <0.5 we call it equivocal
all.res <- NULL
for (j in 1:length(rfp.res)){
  curr <- data.frame(rfp.res[[j]])
  for(t in 1:nrow(curr)){
    if(any(curr[t,]>0.5)){all.res <- append(all.res, names(curr[t,])[which(curr[t,]>0.5)])}
    if(!any(curr[t,]>0.5)){all.res <- append(all.res,"equivocal")}
  }
}

# extract the true states of each node
all.true <- NULL; for(k in 1:nsim){all.true <- append(all.true, anc.res[[k]]$state)}
# get the node numbers the states apply to
all.node <- NULL; for(k in 1:nsim){all.node <- append(all.node, rownames(anc.res[[k]]))}

# bind these data into a single dataframe
sim.estim <- data.frame(sim.state=all.true, estim.state=all.res, node=all.node)
# compare the estimated and true states to see where they agree/disagree
sim.estim$match <- sim.estim$sim.state==sim.estim$estim.state
# determine how frequently the RandomForest approach is wrong/right
table(sim.estim$match)/(nsim*Nnode(agam.tree))
# 80% accuracy isn't bad.

table(sim.estim$estim.state)[1:5]/table(sim.estim$sim.state)

res.a <- dplyr::filter(sim.estim, sim.state=="a")

by.state <- sim.estim %>%
  group_by(sim.state) %>%
  select(estim.state) %>%
  table()
sim.numbers <- table(sim.estim$sim.state)


# if we want to assess accuracy by node, we might expect older to be less accurate
node.acc <- as.data.frame.matrix(sim.estim %>%
                select(node, match) %>%
                table()/nsim)
node.acc$height <- unlist(node.depth.edgelength(agam.tree)[120:237])
colnames(node.acc) <- c("incorrect","correct","height")
library(ggrepel)
ggplot() +
  #geom_text(data=node.acc, aes(x=height,y=correct,label=node),nudge_x=2,nudge_y=0.01,size=3) +
  geom_label_repel(data=node.acc,aes(x=height,y=correct,label=node),
                   box.padding   = 0.15, 
                   #point.padding = 0.5,
                   label.padding = 0.1,
                   label.size=0.1,
                   max.overlaps=20,
                   segment.color = 'grey50',
                   size=2) +  
  geom_point(data=node.acc, aes(x=height,y=correct)) + theme_bw()
# older nodes are less accurate, younger nodes more


############################################################################

# Plot the marginal ancestral states on the phylogeny
# for both the truth and RandomForest

# choose which simulation you want to see
state <- 500

# visualize the TRUE tree
plot(sim.trees[[state]],cols,ftype="off",type="arc",
     arc_height=0.2,lwd=3)

# prepare the data
# match up the regimes iwth the data
sim.reg <- sim.trees[[state]]$states

# we'll use the phytools architecture to plot the RF likelihoods for ancestral states
#fit.STP.full <- fitMk(agam.tree, reg, model=STP.full, pi="fitzjohn")
fit.ER <- fitMk(sim.trees[[state]], sim.reg, model="ER", pi="fitzjohn")

# extract the marginal ancestral state estimates
anc <- ancr(fit.ER, type="marginal")

# set colors
cols<-setNames(hcl.colors(n=5),c("a","b","c","d","e"))

# set the pies to be bigger if there is greater uncertainty (no state with >50%)
node.cex<-apply(anc$ace,1,function(x) if(any(x>0.5)) 0.3 else 0.6)

# plot the tree with the ER result
plot(anc,
     args.plotTree=list(type="arc", arc_height=0.5, fsize=0.5, offset=3, color="grey"),
     args.nodelabels=list(cex=node.cex,piecol=cols),
     args.tiplabels=list(cex=0.2,piecol=cols),
     legend=F)
legend("bottom",c("a","b","c","d","e"),
       pch=16,col=cols,cex=1,pt.cex=2,bty="n",horiz=T)

# Now plot the RandomForest result
# replace the ancr ace matrix with our RF predictions
anc$ace <- rfp.res[[state]]

# plot the RandomForest prediction
plot(anc,
     args.plotTree=list(type="arc", arc_height=0.5, fsize=0.5, offset=3, color="grey"),
     args.nodelabels=list(cex=node.cex,piecol=cols),
     args.tiplabels=list(cex=0.2,piecol=cols),
     legend=F)
legend("bottom",c("a","b","c","d","e"),
       pch=16,col=cols,cex=1,pt.cex=2,bty="n",horiz=T)

# mostly looks like increased uncertainty at the oldest nodes

############################################################################

pca1 <- prcomp(sim.res[[5]][,1:19]); summary(pca1)
pca1 <- data.frame(pca1$x)
pca1$Regime <- sim.res[[1]]$Regime

library(ggplot2)
ggplot() +
  geom_point(data=pca1, aes(x=PC1,y=PC2,fill=Regime),shape=21,size=3) +
  theme_bw()

pca.emp <- prcomp(rfdata.reg[,1:19]); summary(pca.emp)
pca.emp <- data.frame(pca.emp$x)
pca.emp$Regime <- rfdata.reg$Regime

ggplot() +
  geom_point(data=pca.emp, aes(x=PC1,y=PC2,fill=Regime),shape=21,size=3)

lapply(sim.trees, function(x) plot(x,cols,ftype="off",type="arc",
                                   arc_height=0.2,lwd=3))
sapply(1:10, function(x) plot(sim.trees[[x]],cols,ftype="off",type="arc",
                              arc_height=0.2,lwd=3))

emp.reg <- setNames(rfdata.reg$Regime, rownames(rfdata.reg))
fit.ER <- fitMk(agam.tree, emp.reg, model="ER", pi="fitzjohn")
fit.SY <- fitMk(agam.tree, emp.reg, model="SYM", pi="fitzjohn")

