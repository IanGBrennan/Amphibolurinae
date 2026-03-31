
# set the working directory to the GitHub repo
setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

############################################################################

# Load necessary packages
library(phytools)

# load the agamid data
load("Data/Amphibolurinae_Data.RData")

############################################################################

eco.df <- read.csv("Data/Amphibolurinae_Ecology.csv", h=T)
eco <- setNames(eco.df$Ecology, eco.df$Genus_species)
levels(eco) <- sort(unique(eco))

############################################################################

# fit Markov models with phytools under different evolutionary scenarios

# Equal Rates among all states
fit.ER  <- fitMk(agam.tree, eco, model="ER", pi="fitzjohn")

# All-Rates Different
fit.ARD <- fitMk(agam.tree, eco, model="ARD", pi="fitzjohn")

# Symmetric rates among states
fit.SYM <- fitMk(agam.tree, eco, model="SYM", pi="fitzjohn")

# This is the hidden rates model, but maybe not appropriate here
# it also provides a poor fit to this dataset
#fit.HRM <- fitHRM(agam.tree, eco, ncat=2, pi="fitzjohn", parallel=TRUE, niter=20)

# Now we can build some more complex models by specifying the
# rates we're interested in estimating

# A Stepped model with Symmetric rates
# Saxicolous <--> Terrestrial <--> Semiarboreal <--> Arboreal
STP <- matrix(c(
  0,0,1,0,
  0,0,0,3,
  1,0,0,2,
  0,3,2,0),
  4,4,
  byrow=TRUE,
  dimnames=list(sort(unique(eco)),sort(unique(eco)))
)
# fit the Stepped model with symmetric rates
fit.STP.symm <- fitMk(agam.tree, eco, model=STP, pi="fitzjohn")

# A Stepped model with independent rates
# Saxicolous <- -> Terrestrial <- -> Semiarboreal <- -> Arboreal
STP.full <- matrix(c(
  0,0,1,0,
  0,0,0,3,
  2,0,0,5,
  0,4,6,0),
  4,4,
  byrow=TRUE,
  dimnames=list(sort(unique(eco)),sort(unique(eco)))
)
# fit the Stepped model with independent rates
fit.STP.full <- fitMk(agam.tree, eco, model=STP.full, pi="fitzjohn")
plot(fit.STP.full, show.zeros=T, color=T, width=T)
as.Qmatrix(fit.STP.full)

############################################################################

# compare the models and plot the transition matrix of the best one
anova(fit.ER, fit.ARD, fit.SYM, fit.STP.symm, fit.STP.full)

# best model fit
plot(fit.STP.symm, show.zeros=T, width=T, color=T)
plot(fit.STP.full, show.zeros=T, width=T, color=T)

############################################################################

# extract the marginal ancestral state estimates
eco.anc <- ancr(fit.STP.full, type="marginal")

# Plot the marginal ancestral states on the phylogeny
# set colors
cols <- setNames(c(hcl.colors(8)[c(1,3,5,7)]),
               levels(eco))
node.cex<-apply(eco.anc$ace,1,
                function(x) if(any(x>0.90)) 0.3 else 0.6)
# plot tree
plot(eco.anc,
     args.plotTree=list(type="arc", arc_height=0.5, fsize=0.5, offset=3, color="grey"),
#     args.plotTree=list(type="fan",part=0.75, color="grey", fsize=0.5,offset=3),
     args.nodelabels=list(cex=node.cex,piecol=cols),
     args.tiplabels=list(cex=0.2,piecol=cols),
     legend=FALSE)
# plot legend
legend(x=-40, y=0,
       levels(eco),pch=16,col=cols,
       horiz=T,cex=0.8,bty="n",pt.cex=2,
       x.intersp=1)
















############################################################################
# You can also fit a model with a custom Qmatrix but the 
# behavior is unreliable!

# A Stepped model where 
test <- matrix(c(
  0,0,1,0,
  0,0,0,3,
  2,0,0,5,
  0,4,0,0),
  4,4,
  byrow=TRUE,
  dimnames=list(sort(unique(eco)),sort(unique(eco)))
)
test.qmat <- as.Qmatrix(fit.STP.full)
test.qmat[1,3] <- 0.01888531
test.qmat[3,1] <- 0.01340638
#test.qmat[4,3] <- 0.1

fit.TEST <- fitMk(agam.tree, eco, model=test, pi="fitzjohn", fixedQ=test.qmat)

#fit.TESTO <- fitMk(agam.tree, eco, model=test, pi="fitzjohn")
rates <- c(as.Qmatrix(fit.TEST)[1,3],
           as.Qmatrix(fit.TEST)[3,1],
           as.Qmatrix(fit.TEST)[2,4],
           as.Qmatrix(fit.TEST)[4,2],
           as.Qmatrix(fit.TEST)[3,4],
           as.Qmatrix(fit.TEST)[4,3])
fit.TEST$rates <- rates
fit.TEST$index.matrix <- test; diag(fit.TEST$index.matrix) <- NA
