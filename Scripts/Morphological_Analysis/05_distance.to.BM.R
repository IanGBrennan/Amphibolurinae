
setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

#######################################################################

library(ape)
library(phytools)
library(mvMORPH)
library(RColorBrewer)

load("Data/Amphibolurinae_Data.RData")
load("Data/Ancestral_Trait_Estimates.RData")

#######################################################################

# combine the extant and estimated ancestral traits together

# isolate traits (exclude Genus, Genus_species)
LSR.anc <- allLSR[,1:19]
# combine observed and ancestral traits into a single dataframe
LSR.anc <- dplyr::bind_rows(LSR.anc, ancestors)

#######################################################################

# function to get the euclidean distance between two points
euclidean <- function(a, b) sqrt(sum((a - b)^2))

# function to calculate the multivariate distances traveled across a tree
dist.nodes <- function(phy, trait.df, metric=c("RootToChild", "ParentToChild", "RootParentChild"), distance=c("euclidean","mahalanobis")){
  
  # build a dataframe of basic tree statistics
  ndel <- data.frame(edge = 1:nrow(phy$edge),
                     node.parent = phy$edge[,1],
                     node.child  = phy$edge[,2],
                     length = phy$edge.length,
                     name.parent = paste0("n",phy$edge[,1]))
  ndel$name.child <- sapply(ndel$node.child, function(x) ifelse(x <= Ntip(phy), phy$tip[[x]], paste0("n",x)))
  
  # estimate the multivariate distance between each parent/child node (along each edge)
  if(metric == "ParentToChild"){
    edist <- NULL
    for (k in 1:nrow(ndel)){
      parent <- trait.df[which(rownames(trait.df)==ndel$name.parent[[k]]),]
      child  <- trait.df[which(rownames(trait.df)==ndel$name.child[[k]]),]
      if(distance=="euclidean"){cdist <- euclidean(parent, child)}
      if(distance=="mahalanobis"){cdist <- mahalanobis(x=unlist(child), center=unlist(parent), cov=cov(trait.df))}
      edist <- append(edist, cdist)
    }
  }
  
  # estimate the multivariate distance between the root and each node/tip
  if(metric == "RootToChild"){
    rdist <- NULL
    if(distance=="euclidean"){for(k in 1:nrow(ndel)){rdist <- append(rdist, euclidean(trait.df[Ntip(phy)+1,], trait.df[which(rownames(trait.df)==ndel$name.child[[k]]),]))}}
    if(distance=="mahalanobis"){for(k in 1:nrow(ndel)){rdist <- append(rdist, mahalanobis(x=unlist(trait.df[Ntip(phy)+1,]), 
                                                                                          center=unlist(trait.df[which(rownames(trait.df)==ndel$name.child[[k]]),]),
                                                                                          cov=cov(trait.df)))}}
  }

  # estimate the difference in multivariate distance between the root and a pair of parent and child nodes
  if(metric == "RootParentChild"){
    adist <- NULL
    if(distance=="euclidean"){for(k in 1:nrow(ndel)){adist <- append(adist, (euclidean(trait.df[Ntip(phy)+1,], trait.df[which(rownames(trait.df)==ndel$name.child[[k]]),]) - 
                                                                               euclidean(trait.df[Ntip(phy)+1,], trait.df[which(rownames(trait.df)==ndel$name.parent[[k]]),])))}}
  }

  # specify which metric we're using
  if(metric == "ParentToChild"){ndel$distance <- edist}
  if(metric == "RootToChild"){ndel$distance <- rdist}
  if(metric == "RootParentChild"){ndel$distance <- adist}
  
  # rescale the distance metric for ease of plotting
  ndel$distance.scaled <- round((ndel$distance - min(ndel$distance))/diff(range(ndel$distance)) * 99) + 1
  
  # tell us which metric was used
  ndel$metric <- metric

  # kick out the ndel dataframe for use with other methods (e.g. lollipop plotting the distance)
  return(ndel)
}

# function to compare the empirical morphological distances to simulated (BM) data
distance.to.BM <- function(phy, trait.df, nsim=10, tree.type=c("fan","phylogram")){
  
  # subset the trait data to just the tips
  trait.obs <- trait.df[which(rownames(trait.df)%in%phy$tip.label),]
  
  # establish the multivariate distance travelled between parent and child nodes (compared to the root)
  dist.emp <- dist.nodes(phy=phy, trait.df=trait.df, metric="RootParentChild", distance="euclidean")
  
  # fit a Brownian Motion model to each individual trait and extract the theta and sigma estimates
  sigmas <- NULL; thetas <- NULL
  for(k in 1:ncol(all.LSR.traits)){
    curr.fit <- mvBM(tree=phy, data=dplyr::select(trait.obs, colnames(trait.obs)[[k]]), model="BM1", echo=F, diagnostic=F)
    sigmas <- append(sigmas, curr.fit$sigma)
    thetas <- append(thetas, curr.fit$theta)
  }
  
  # simulate nsim number of traits based on the sigma/theta estimated from the empirical data
  # this step is slow, so just be patient
  distances <- NULL
  for(j in 1:nsim){
    # simulate a set of data
    sim.df <- data.frame(sapply(1:ncol(trait.df), function(x) fastBM(tree=phy, a=thetas[x], sig2=sigmas[x], internal=T)))
    # fix the rownames for the internal nodes
    rownames(sim.df)[(Ntip(phy)+1):nrow(trait.df)] <- paste0("n", rownames(sim.df)[(Ntip(phy)+1):nrow(trait.df)])
    # estimate the multivariate trait change
    dist.sim <- dist.nodes(phy=phy, trait.df=sim.df, metric="RootParentChild", distance="euclidean")
    # keep just the distance measured
    distances <- cbind(distances, dist.sim$distance)
    # let me know where we are
    print(j)
  }
  distances <- data.frame(distances)
  # collect the distances as absolute values for comparison
  abs.distances <- abs(distances)
  
  # get the 95% quantiles for the simulated data (needs a function first)
  apply.quantile <- function(x){quantile(x,probs=c(0.05,0.95))}
  sim.ci <- data.frame(t(apply(abs.distances, 1, apply.quantile)))
  colnames(sim.ci) <- c("lower","upper")
  # get the mean of the simulated data (needs to be from the absolute values)
  sim.ci$mean <- apply(sim.ci, 1, mean)
  # combine the empirical with the quantiles from the simulated data
  full.ci <- cbind(dist.emp, sim.ci)
  # make sure to take the absolute values of the empirical data too
  full.ci$distance.abs <- abs(full.ci$distance)
  
  results <- NULL; color <- NULL
  for(x in 1:nrow(full.ci)){
    #if(full.ci[x,"distance.abs"] > full.ci[x,"upper"]){results<-append(results,"greater")}
    #if(full.ci[x,"distance.abs"] < full.ci[x,"lower"]){results<-append(results, "less")}
    if(all(full.ci[x,"distance.abs"] > full.ci[x,c("upper","lower")])){results<-append(results,"greater"); color<-append(color,"blue"); next}
    if(all(full.ci[x,"distance.abs"] < full.ci[x,c("upper","lower")])){results<-append(results,"less"); color<-append(color,"red"); next}
    else{results<-append(results,"expected");color<-append(color,"grey")}
  }
  # append the results (greater, less, expected), color, and fold-changes
  full.ci$results <- results
  full.ci$edge.color <- color
  full.ci$change.fold <- full.ci$distance.abs / full.ci$mean
  
  # make a second dataframe to plot the results with less change than expected
  full.less <- full.ci
  full.less$change.fold <- full.less$mean / full.less$distance.abs
  
  # this is how you'd scale it normally, but we have an extreme outlier
  #full.ci$change.scaled <- round((full.ci$change.fold - 0)/diff(c(0,max(full.ci$change.fold))) * 99) + 1
  
  # so instead we'l do ti like this
  full.ci$change.scaled <- round((full.ci$change.fold - 0)/diff(c(0,sort(full.ci$change.fold,partial=nrow(full.ci)-1)[nrow(full.ci)-1])) * 99) + 1
  full.ci$change.scaled[which(full.ci$change.scaled>100)] <- 100
  
  full.less$change.scaled <- round((full.less$change.fold - 0)/diff(c(0,max(full.less$change.fold))) * 99) + 1
  
  
  colors.change <- ((colorRampPalette(brewer.pal(9, "YlOrRd"))(100)))
  #colors.change <- (colorRampPalette(viridis::viridis(n=12)[3:12])(100))
  full.ci$color.scaled <- colors.change[full.ci$change.scaled]
  full.less$color.scaled <- colors.change[full.less$change.scaled]
  
  
  for(p in 1:nrow(full.ci)){if(full.ci$results[[p]]%in%c("less","expected")){full.ci$color.scaled[[p]]<-"lightgrey"}}
  for(p in 1:nrow(full.ci)){if(full.ci$results[[p]]%in%c("less","expected")){full.ci$change.fold[[p]]<-1}}
  
  for(p in 1:nrow(full.less)){if(full.less$results[[p]]%in%c("greater","expected")){full.less$color.scaled[[p]]<-"lightgrey"}}
  for(p in 1:nrow(full.less)){if(full.less$results[[p]]%in%c("greater","expected")){full.less$change.fold[[p]]<-1}}
  
  
  # if you want to scale the edge lengths to be proportional to the fold change
  phy2 <- phy
  phy2$edge.length <- full.ci$length * full.ci$change.fold
  full.ci$length.scaled <- full.ci$length * full.ci$change.fold
  
  # plot the results on the tree
  plot.phylo(phy, edge.color = unlist(full.less$color.scaled), edge.width=4, type=tree.type, cex=0.3)
  title("Less than Expected Change")
  
  plot.phylo(phy, edge.color = unlist(full.ci$color.scaled), edge.width=4, type=tree.type, cex=0.3)
  title("Greater than Expected Change")

  return(full.ci)
  
}
# distance is measured as:
  # RootToChild: the euclidean distance between the root and each node/tip
  # ParentToChild: the euclidean distance between each parent node and the focal descendant node
  # RootParentChild: the difference in euclidean distance between the root-to-parent and root-to-child distances
# the function fits a constant rate BM to each trait and extracts fitted sigma/theta values
# simulates nsim datasets using estimated values, then compares the empirical morphological distances to the simulated
# differences between empirical/simulated distances are scaled and plotted
# plots two trees (1) Less than Expected Change, and (2) Greater than Expected Change


#######################################################################

# fit the comparison function to our empirical data
distRes <- distance.to.BM(agam.tree, LSR.anc, nsim=100, tree.typ="fan")



