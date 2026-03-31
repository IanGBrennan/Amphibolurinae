setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")
############################################################################

library(ggplot2)

############################################################################

load("Data/Amphibolurinae_Data.RData")
load("Data/Ancestral_Trait_Estimates.RData")

############################################################################

source("Scripts/morphotrajectory.R")

############################################################################

# function to get the euclidean distance between two points
euclidean <- function(a, b) sqrt(sum((a - b)^2))

############################################################################

# We'll focus on the species mean data

# isolate traits (exclude Genus, Genus_species)
LSR.anc <- allLSR[,1:19]
# combine observed and ancestral traits into a single dataframe
LSR.anc <- dplyr::bind_rows(LSR.anc, ancestors)
# correct some weird ancestral trait value for n120's foot
LSR.anc[120,"Foot"] <- 0.4977063

############################################################################

# we might want to identify which traits are the most variable 
trait.variance <- apply(LSR.anc, 2, function(x) max(x) - min(x))
trait.variance[order(trait.variance, decreasing = T)]

############################################################################

# transform Log-Shape Ratios back to real measurements
unLSR.anc <- exp(LSR.anc)
unLSR.anc <- unLSR.anc[,1:18] * unLSR.anc$Size

# euclidean distances from log-shape ratios
edist <- apply(unLSR.anc[1:119,], 1, function(x) euclidean(x, unLSR.anc[120,]))
euc.df <- data.frame(taxon = names(edist), distance = edist)
#View(euc.df[order(euc.df$distance),])
euc.df$genus <- sapply(euc.df$taxon, function(x) strsplit(x,"_")[[1]][1])

# plot the multivariate euclidean distance of each species from the MRCA
ggplot(euc.df) +
  geom_segment(aes(x=0, xend=distance, y=taxon, yend=taxon,color=genus)) +
  geom_point(aes(x=distance, y=taxon, fill=genus), shape=21, size=3) +
  scale_y_discrete(limits=rev) +
  xlab("Euclidean Multivariate Distance from Amphibolurine MRCA") +
  theme_classic() + theme(legend.position="none")

############################################################################

# mahalanobis distances from log-shape ratios
mdist <- mahalanobis(x=unLSR.anc[1:119,], center=unlist(unLSR.anc[120,]), cov=cov(unLSR.anc[1:119,]))
mah.df <- data.frame(taxon = names(mdist), distance = mdist)
mah.df$distance.min <- mah.df$distance - min(mah.df$distance)

############################################################################

# plot the multivariate euclidean distance
plotTree.lollipop(agam.tree, dplyr::select(euc.df, distance),
                  args.plotTree=list(fsize=0.5,ftype="i"))

# plot the multivariate mahalanobis distance
plotTree.lollipop(agam.tree, dplyr::select(mah.df, distance.min),
                  args.plotTree=list(fsize=0.5,ftype="i"))


############################################################################


# correct some weird ancestral trait value for n120's foot
LSR.anc[120,"Foot"] <- 0.4977063

# extract the distance between each species and the MRCA (n120) for each trait using shape ratios
LSR.mrca <- apply(LSR.anc[1:119,], 1, function(x) (x - LSR.anc[120,]))
LSR.mrca <- do.call(rbind.data.frame, LSR.mrca)

# plot the leg traits
plotTree.lollipop(agam.tree, LSR.mrca[,c("UpperArm","LowerArm","Hand",
                                         "UpperLeg","LowerLeg","Foot")],
                  args.plotTree=list(fsize=0.5,ftype="i"))

# plot the body traits
plotTree.lollipop(agam.tree, LSR.mrca[,c("Interlimb","BodyWidth",
                                         "PelvicWidth","PelvicHeight")],
                  args.plotTree=list(fsize=0.5,ftype="i"))

# plot the head traits
plotTree.lollipop(agam.tree, LSR.mrca[,c("HeadWidth","HeadDepth",
                                         "EyeDiameter","SnoutEye")],
                  args.plotTree=list(fsize=0.5,ftype="i"))

# plto the tail and size traits
plotTree.lollipop(agam.tree, LSR.mrca[,c("TailLength","TailWidth",
                                         "Size")],
                  args.plotTree=list(fsize=0.5,ftype="i"))


LSR.anc.center <- LSR.anc - LSR.anc[rep(129,237),]
plotTree.lollipop(agam.tree, LSR.anc.center[1:119,c("BodyWidth","TailLength")],
                  args.plotTree=list(fsize=0.5,ftype="i"))



# this is the lollipop plot method provided in new versions of phytools
plotTree.lollipop<-function(tree,x,
                            args.plotTree=list(),args.lollipop=list(),...){
  if(!inherits(x,c("matrix","data.frame"))) x<-as.matrix(x)
  h<-max(nodeHeights(tree))
  if(hasArg(panel_height)) panel_height<-list(...)$panel_height
  else panel_height<-1.0
  panel_height<-panel_height*h
  args.plotTree$tree<-tree
  args.plotTree$direction<-"upwards"
  if(is.null(args.plotTree$mar)) 
    args.plotTree$mar<-c(0.1,5.1,0.1,0.1)
  if(is.null(args.plotTree$ylim)) 
    args.plotTree$ylim<-c(0,h+ncol(x)*panel_height)
  if(is.null(args.plotTree$ftype)) 
    args.plotTree$ftype<-"off"
  if(is.null(args.plotTree$lwd)) args.plotTree$lwd<-1
  do.call(plotTree,args.plotTree)
  pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  if(pp$font){
    dx<-abs(diff(pp$x.lim))
    pdin<-par()$din[2]
    sh<-(pp$cex*strwidth(paste(" ",tree$tip.label,sep=""))+
           0.3*pp$cex*strwidth("W"))*(par()$din[1]/par()$din[2])*
      (diff(par()$usr[3:4])/diff(par()$usr[1:2]))
    new_h<-h+max(sh)
    panel_height<-(h-new_h+ncol(x)*panel_height)/ncol(x)
    h<-new_h
  }
  if(hasArg(ylab)) ylab<-list(...)$ylab
  else ylab<-if(!is.null(colnames(x))) 
    colnames(x) else rep("",ncol(x))
  for(i in ncol(x):1){
    d<-max(c(diff(range(x[,i])),max(x[,i])))
    y<-setNames(x[,i]/d*0.8*panel_height,rownames(x))
    lower<-h+(i-1)*panel_height+panel_height*0.05
    upper<-h+(i-1)*panel_height+panel_height*0.95
    polygon(c(0,max(pp$xx)+1,max(pp$xx)+1,0),
            c(lower,lower,upper,upper),
            border=FALSE,col="#F2F2F2")
    hh<-lower-min(c(0,min(y)))+0.05*panel_height
    lines(range(pp$xx),rep(hh,2),col="black",lty="dotted")
    segments(x0=pp$xx[1:Ntip(tree)],y0=rep(hh,Ntip(tree)),
             x1=pp$xx[1:Ntip(tree)],y1=y[tree$tip.label]+hh)
    labs<-pretty(c(min(c(0,min(x[,i]))),x[,i]),n=4)
    labs[!(labs>max(x[,i]))]->labs
    labs[!(labs<min(c(0,min(x[,i]))))]->labs
    axis(2,at=hh+max(y)/max(x[,i])*labs,
         labels=labs,las=1,cex.axis=0.6)
    args.lollipop$bg<-setNames(
      hcl.colors(n=100)[ceiling(99*((y-
                                       min(y))/diff(range(y))))+1],
      names(y))
    args.lollipop$bg<-args.lollipop$bg[tree$tip.label]
    if(is.null(args.lollipop$pch)) args.lollipop$pch<-21
    if(is.null(args.lollipop$cex)) args.lollipop$cex<-1.2
    args.lollipop$x<-pp$xx[1:Ntip(tree)]
    args.lollipop$y<-y[tree$tip.label]+hh
    do.call(points,args.lollipop)
    mtext(ylab[i],2,line=3,at=mean(hh+max(y)/
                                     max(x[,i])*labs),cex=0.8)
  }
}


# create a plotting function that visualizes the multivariate phenotypic change
# either along individual branches (ParentToChild)
# or between the root and each tip/node (RootToChild)
edist.btwnodes <- function(phy, trait.df, lwd=2, palette, tree.type=c("fan","phylogram"), 
                           metric=c("RootToChild", "ParentToChild", "RootParentChild"), distance=c("euclidean","mahalanobis")){
  
  # build a dataframe of basic tree statistics
  ndel <- data.frame(edge = 1:nrow(phy$edge),
                     node.parent = phy$edge[,1],
                     node.child  = phy$edge[,2],
                     length = phy$edge.length,
                     edge.color = "lightGrey",
                     edge.width = lwd,
                     name.parent = paste0("n",phy$edge[,1]))
  ndel$name.child <- sapply(ndel$node.child, function(x) ifelse(x <= Ntip(phy), phy$tip[[x]], paste0("n",x)))
  
  # estimate the multivariate distance between each parent/child node (along each edge)
  edist <- NULL
  for (k in 1:nrow(ndel)){
    parent <- trait.df[which(rownames(trait.df)==ndel$name.parent[[k]]),]
    child  <- trait.df[which(rownames(trait.df)==ndel$name.child[[k]]),]
    if(distance=="euclidean"){cdist <- euclidean(parent, child)}
    if(distance=="mahalanobis"){cdist <- mahalanobis(x=unlist(child), center=unlist(parent), cov=cov(trait.df))}
    edist <- append(edist, cdist)
  }
  
  # estimate the multivariate distance between the root and each node/tip
  rdist <- NULL
  if(distance=="euclidean"){for(k in 1:nrow(ndel)){rdist <- append(rdist, euclidean(trait.df[Ntip(phy)+1,], trait.df[which(rownames(trait.df)==ndel$name.child[[k]]),]))}}
  if(distance=="mahalanobis"){for(k in 1:nrow(ndel)){rdist <- append(rdist, mahalanobis(x=unlist(trait.df[Ntip(phy)+1,]), 
                                                                                       center=unlist(trait.df[which(rownames(trait.df)==ndel$name.child[[k]]),]),
                                                                                       cov=cov(trait.df)))}}
  
  # estimate the difference in multivariate distance between the root and a pair of parent and child nodes
  adist <- NULL
  if(distance=="euclidean"){for(k in 1:nrow(ndel)){adist <- append(adist, (euclidean(trait.df[Ntip(phy)+1,], trait.df[which(rownames(trait.df)==ndel$name.child[[k]]),]) - 
                                                                           euclidean(trait.df[Ntip(phy)+1,], trait.df[which(rownames(trait.df)==ndel$name.parent[[k]]),])))}}
  
  # specify which metric we're using
  if(metric == "ParentToChild"){ndel$distance <- edist}
  if(metric == "RootToChild"){ndel$distance <- rdist}
  if(metric == "RootParentChild"){ndel$distance <- adist}
  
  # rescale the distance metric for ease of plotting
  ndel$distance.scaled <- round((ndel$distance - min(ndel$distance))/diff(range(ndel$distance)) * 99) + 1
  
  # make a color ramp that suits the data
  if(palette %in% c("magma", "inferno", "plasma", "viridis","cividis", "rocket", "mako")){
    new.cols <- viridis::viridis(n=100, option=palette)
  }else{
    col.ramp <- colorRampPalette(RColorBrewer::brewer.pal(9, palette))
    new.cols <- (col.ramp(100))
    if(palette %in% c("Spectral","RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn","PiYG","BrBG","YlOrRd","YlGnBu")){new.cols <- rev(new.cols)}
  }
  
  # match the scaled distances to the available colors (1-100)
  ndel$edge.color <- unlist(sapply(ndel$distance.scaled, function(x) new.cols[x]))
  
  # plot the tree
  plot.phylo(phy, edge.color = unlist(ndel$edge.color), edge.width = 4, type = tree.type, cex = 0.3)
  
  # provide information about which metric the user applied
  if(metric == "ParentToChild" & distance == "euclidean"){message("The distance variable represents euclidean distance between each parent and child node pair along each edge")}
  if(metric == "RootToChild" & distance == "euclidean"){message("The distance variable represents euclidean distance between the root node and each child node and not along the tree path")}
  if(metric == "ParentToChild" & distance == "mahalanobis"){message("The distance variable represents mahalanobis distance between each parent and child node pair along each edge")}
  if(metric == "RootToChild" & distance == "mahalanobis"){message("The distance variable represents mahalanobis distance between the root node and each child node and not along the tree path")}
  if(metric == "RootParentChild" & distance == "euclidean"){message("The distance variable represents the difference euclidean distance between the root and each parent and child node pair")}
  
  # kick out the ndel dataframe for use with other methods (e.g. lollipop plotting the distance)
  return(ndel)
}

# plot the ParentToChild distances (edgewise)
dist.edges <- edist.btwnodes(phy = agam.tree, trait.df = LSR.anc, tree.type = "fan", 
                             palette = "RdYlBu", metric="ParentToChild", distance="euclidean")
# plot the RootToChild distances (root-to-node/tip-wise)
dist.2root <- edist.btwnodes(phy = agam.tree, trait.df = LSR.anc, tree.type = "fan",
                             palette = "RdYlBu", metric="RootToChild", distance="euclidean")
# plot the RootParentChild distances (edgewise)
dist.rpc <- edist.btwnodes(phy = agam.tree, trait.df = LSR.anc, tree.type = "fan", 
                             palette = "RdYlBu", metric="RootParentChild", distance="euclidean")

plot(dist.edges$distance ~ dist.rpc$distance)

# if we wanted to extract the estimated distances to plot a different way (lollipop)
r2c.dists <- data.frame(setNames(dist.2root$distance, dist.2root$name.child))
r2c.dists <- subset(r2c.dists, rownames(r2c.dists) %in% agam.tree$tip.label)
names(r2c.dists) <- "distance"

plotTree.lollipop(agam.tree, r2c.dists,
                  args.plotTree=list(fsize=0.5,ftype="i"))


pca.x


