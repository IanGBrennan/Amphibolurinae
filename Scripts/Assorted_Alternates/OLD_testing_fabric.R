setwd("/Users/ianbrennan/Documents/GitHub/EARtH")

anc.fabric <- function (tree, x, maxit = 2000) 
{
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  if (is.ultrametric(tree)) 
    cat("Warning: the trend model is generally non-identifiable for ultrametric trees.\n")
  tol <- 1e-08
  if (any(tree$edge.length <= (10 * .Machine$double.eps))) 
    stop("some branch lengths are 0 or nearly zero")
  D <- dist.nodes(tree)
  ntips <- length(tree$tip.label)
  Cii <- D[ntips + 1, ]
  C <- D
  C[, ] <- 0
  counts <- vector()
  for (i in 1:nrow(D)) for (j in 1:ncol(D)) C[i, j] <- (Cii[i] + 
                                                          Cii[j] - D[i, j])/2
  dimnames(C)[[1]][1:length(tree$tip)] <- tree$tip.label
  dimnames(C)[[2]][1:length(tree$tip)] <- tree$tip.label
  C <- C[c(1:ntips, (ntips + 2):nrow(C)), c(1:ntips, (ntips + 
                                                        2):ncol(C))]
  x <- x[tree$tip.label]
  likelihood <- function(theta, x, C) {
    a <- theta[1]
    u <- 0
    sig2 <- theta[2]
    y <- theta[3:length(theta)]
    logLik <- mnormt::dmnorm(x = c(x, y), mean = (a + diag(C) * u), 
                     varcov = sig2 * C, log = TRUE)
    return(-logLik)
  }
  a <- mean(x)
  sig2 <- var(x)/max(C)
  result <- optim(par = c(a, sig2, rep(a, tree$Nnode - 1)), 
                  likelihood, x = x, C = C, method = "L-BFGS-B", 
                  lower = c(-Inf, tol, rep(-Inf, tree$Nnode - 1)), control = list(maxit = maxit))
  ace <- c(result$par[c(1, 3:length(result$par))])
  names(ace) <- c(as.character(tree$edge[1, 1]), rownames(C)[(length(tree$tip.label) + 
                                                                1):nrow(C)])
  obj <- list(ace = ace, sig2 = result$par[2], 
              logL = -result$value, convergence = result$convergence, 
              message = result$message)
  class(obj) <- "anc.trend"
  obj
}

#testo <- anc.fabric(tree.trend, x=tip.traits)

########################################################################

# this method isn't right because the C matrix represents the 
# full height above the root for each tip/node and not just the
# length of each individual branch. So adjusting via the input
# beta vector over-corrects
anc.fabric2 <- function (tree, x, beta, maxit = 2000) 
{
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  if (is.ultrametric(tree)) 
    cat("Warning: the trend model is generally non-identifiable for ultrametric trees.\n")
  tol <- 1e-08
  if (any(tree$edge.length <= (10 * .Machine$double.eps))) 
    stop("some branch lengths are 0 or nearly zero")
  D <- dist.nodes(tree)
  ntips <- length(tree$tip.label)
  Cii <- D[ntips + 1, ]
  C <- D
  C[, ] <- 0
  counts <- vector()
  for (i in 1:nrow(D)) for (j in 1:ncol(D)) C[i, j] <- (Cii[i] + 
                                                          Cii[j] - D[i, j])/2
  dimnames(C)[[1]][1:length(tree$tip)] <- tree$tip.label
  dimnames(C)[[2]][1:length(tree$tip)] <- tree$tip.label
  C <- C[c(1:ntips, (ntips + 2):nrow(C)), c(1:ntips, (ntips + 
                                                        2):ncol(C))]
  x <- x[tree$tip.label]
  likelihood <- function(theta, x, C) {
    a <- theta[1]
    u <- beta
    sig2 <- theta[2]
    y <- theta[3:length(theta)]
    logLik <- mnormt::dmnorm(x = c(x, y), mean = (a + diag(C) * u), 
                             varcov = sig2 * C, log = TRUE)
    return(-logLik)
  }
  a <- mean(x)
  sig2 <- var(x)/max(C)
  result <- optim(par = c(a, sig2, rep(a, tree$Nnode - 1)), 
                  likelihood, x = x, C = C, method = "L-BFGS-B", lower = c(-Inf, 
                                                                           tol, rep(-Inf, tree$Nnode - 1)), control = list(maxit = maxit))
  ace <- c(result$par[c(1, 3:length(result$par))])
  names(ace) <- c(as.character(tree$edge[1, 1]), rownames(C)[(length(tree$tip.label) + 
                                                                1):nrow(C)])
  obj <- list(ace = ace, sig2 = result$par[2], 
              logL = -result$value, convergence = result$convergence, 
              message = result$message)
  class(obj) <- "anc.trend"
  obj
}

#testa <- anc.fabric2(tree=tree, x=tip.traits, beta=test.beta)

########################################################################

# Testing the effect of scaling the tree by Beta
# actually: (Mean (Beta * BL) NZ)

# read in an example tree I've made
tree <- read.tree("fabric_test/test.tre")
# visualize the tree
plot(tree); nodelabels(fr="c", cex=1); axisPhylo()

# generate some traits that have evolved under BM
traits <- fastBM(tree, sig2=0.15, a=5, internal=T)
# isolate just the traits at the tips
tip.traits <- traits[1:5]
phenogram(tree, x=tip.traits)
# save the tree and data together
save(tree, traits, file="fabric_test/testdata1.RData")

# rescale the tree by a Beta of -0.5 (decrease of 0.5 per unit time)
tree.trend <- tree
tree.trend$edge.length[3] <- 7 + (7*0.11) # scale this edge by Beta -0.5
plot(tree.trend); nodelabels(fr="c", cex=1); axisPhylo()
# I have to say, I don't think this is right, as it suggests there's less change
# happening on that branch, when in actually it should just be a decrease
# in the trait value

# rescale the tree by a 

# let's estimate ancestral states under a couple different scenarios
# (1) BM
traits.anc.BM <- fastAnc(tree, x=tip.traits)
# (2) BM but on the transformed tree
traits.anc.trend <- fastAnc(tree.trend, x=tip.traits)
# (3) Trend
trend.fit <- anc.trend(tree, x=tip.traits)
# (4) Fabric
fabric.fit <- anc.fabric(tree.trend, x=tip.traits)
# (5) Fabric v2
test.beta <- c(0,0,0,0,0,0,0.11,0)
fabric2.fit <- anc.fabric2(tree, x=tip.traits, beta=test.beta)

phenogram(tree, x=c(tip.traits,traits.anc.BM), sub="anc under BM")
phenogram(tree, x=c(tip.traits,traits.anc.trend), sub="anc under Trend")
phenogram(tree, x=c(tip.traits,trend.fit$ace), sub="anc under anc.Trend")
phenogram(tree, x=c(tip.traits,fabric.fit$ace), sub="anc under anc.fabric")
phenogram(tree, x=c(tip.traits,fabric2.fit$ace), sub="anc under anc.fabric2")


trend.fit <- anc.trend(tree, x=tip.traits)
phenogram(tree, x=traits)
phenogram(tree, x=c(tip.traits,trend.fit$ace))
phenogram(tree, x=tip.traits)






# Ideas:
# use anc.fabric with a provided 

# Can't scale tree/edges by Beta because 



########################################################################
# TROUBLESHOOTING WITH THE MICROHYLIDS

micro.tree <- read.nexus("/Applications/BayesTraitsV4/Microhylidae/Microhylidae_BT.tre")
plot(micro.tree); nodelabels(fr="c", cex=0.5); axisPhylo()
micro.data <- read.csv("/Applications/BayesTraitsV4/Microhylidae/Microhylidae_SIZE.txt", header=F, sep="\ ")
SIZE <- micro.data[,2]; names(SIZE) <- micro.data[,1]
SIZE.df <- data.frame(SIZE); rownames(SIZE.df) <- names(SIZE)
phenogram(micro.tree, SIZE)

micro.beta <- rep(0, nrow(micro.tree$edge))
micro.tree$edge
micro.beta[[27]] <- 0.62 #/ 3.4967   
micro.beta[[23]] <- 0.83 #/ 1.2804 
micro.beta[[10]] <- 0.82 #/ 4.6392 

micro.tree.trend <- micro.tree
micro.tree.trend$edge.length <- micro.tree$edge.length + (micro.tree$edge.length * micro.beta)
plot(micro.tree.trend, cex=0.5); axisPhylo()

# estimate ancestral states under BM
micro.bm <- fastAnc(micro.tree, x = SIZE)
# under fabric
micro.fabric <- anc.fabric(tree = micro.tree.trend, x = SIZE)
# under the bad fabric2
micro.fabric2 <- anc.fabric2(tree = micro.tree, x = SIZE, beta = micro.beta)
# lastly test what happens if we fix some nodes following the fabric results and do the rest under BM
micro.bmtrend <- fastAnc(micro.tree, x = SIZE, anc.states = micro.fabric$ace[which(names(micro.fabric$ace)%in%c(33,32,31,30,24))])
# this last one is a precursor to adding rate shifts

phenogram(tree = micro.tree, x = SIZE, sub="anc under BM")
phenogram(tree = micro.tree, x = c(SIZE, micro.fabric$ace), sub="anc under fabric")
phenogram(tree = micro.tree, x = c(SIZE, micro.fabric2$ace), sub="anc under fabric2")
phenogram(tree = micro.tree, x = c(SIZE, micro.bmtrend), sub="anc under BM with trend")

# saxatilis
# kula & paka
# pet & zweif

########################################################################
# next steps:

# simulate small tree (e.g. 4 tips) and known age
# simulate bigger tree (e.g. 50 tips)

# combine trees to form a single tree (54 tips)
# with a known clade of divergent evolutionary mode

# simulate trait under Trend model on 4 tip tree with fastBM with mu = 0.8 and
# simulate trait under BM on 20 tip tree with 

