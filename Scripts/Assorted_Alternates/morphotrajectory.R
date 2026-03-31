require(ggplot2)
require(RColorBrewer)
require(dplyr)
require(tibble)
require(phytools)
require(ggdensity)
require(mvMORPH)

# the 'morphotrajectory.RR' and 'morphotrajectory.VR' functions take empirical traits and plot 
# a phylomorphospace 

#####################

load("Data/Amphibolurinae_Data.RData")

#####################

# isolate traits (exclude Genus, Genus_species)
LSR.anc <- allLSR[,1:19]
# combine observed and ancestral traits into a single dataframe
LSR.anc <- dplyr::bind_rows(LSR.anc, ancestors)



morphotrajectory <- function(phy, trait, tip.spread, 
                                focus = c("tip", "clade"),
                                PCs = c(1,2),
                                psize=3, lsize=2){

  # check that the tip names are correct
  if(!all(tip.spread %in% phy$tip.label)){stop("WARNING: check tip names")}
  # check that all the tips are in the trait dataframe
#  if(nrow(trait) != length(phy$tip))
#    stop("your trait dataframe must have the same number of rows as species in your tree")
  
#  trait <- trait %>%
#    rownames_to_column('taxon') %>%
#    arrange(match(taxon, phy$tip.label)) %>%
#    column_to_rownames('taxon')
  
  # reduce the dimensionality of your data if >2 traits are included
  if(ncol(trait) > 2){pcs <- data.frame(prcomp(trait)$x)}
  if(ncol(trait) == 2){pcs <- trait}
  
  # separate the tip (empirical) and internal node (ancestral) values
  # and reorder the tip data so that it matches the tree order
  pcs.emp <- pcs[1:Ntip(phy),]
  pcs.emp <- pcs.emp %>%
    tibble::rownames_to_column('taxon') %>%
    arrange(match(taxon, phy$tip.label)) %>%
    tibble::column_to_rownames('taxon')
  pcs.anc <- pcs[(Ntip(phy)+1):nrow(pcs),]
  
  t1 <- c(pcs.emp[,PCs[[1]]], pcs.anc[,PCs[[1]]])
  t2 <- c(pcs.emp[,PCs[[2]]], pcs.anc[,PCs[[2]]])
  
  # specify the root trait data
  root.data <- data.frame(pcs[Ntip(phy)+1,c(PCs[[1]],PCs[[2]])])
  colnames(root.data) <- c("T1","T2")
  
  # extract the xy coordinates of the nodes and segments
  phydat <- data.frame(xstart=t1[phy$edge[,1]],
                       ystart=t2[phy$edge[,1]],
                       xstop=t1[phy$edge[,2]],
                       ystop=t2[phy$edge[,2]],
                       nodestart=phy$edge[,1],
                       nodestop =phy$edge[,2])
  
  # if you're working with a clade
  if(focus == "clade"){
    # determine the MRCA node of your target tips
    mrcn <- getMRCA(phy, tip = tip.spread)
    # get all the descendant nodes of your MRCA including internals
    all.desc <- getDescendants(phy, mrcn)    
  }
  # if you're working with a single tip
  if(focus == "tip"){all.desc <- which(phy$tip.label==tip.spread[[1]])}
  
  # identify which descendant nodes are tips
  tip.desc <- all.desc[which(all.desc <= Ntip(phy))]
  # get the node path from the root to each tip, and make a vector of the unique nodes
  focal.nodes <- unique(unlist(sapply(tip.desc, function(x) nodepath(phy, from=Ntip(phy)+1, to=x))))
  # subset the phydat dataframe to just the nodes of interest
  focal.dat <- dplyr::filter(phydat, nodestart %in% focal.nodes & nodestop %in% focal.nodes)
  focal.dat <- focal.dat[order(focal.dat$nodestart),]
  
  # select the points for just the tips
  tip.dat <- dplyr::filter(focal.dat, nodestop <= Ntip(phy))
  
  # plot the branches from root to tips colored by evolutionary rate, with simulated data
  ggplot()+
#    geom_segment(data=phydat, aes(x=xstart,y=ystart,xend=xstop,yend=ystop), lwd=lsize-1, color="grey") +
    geom_point(data=pcs.anc, aes(x=pcs.anc[,PCs[[1]]], y=pcs.anc[,PCs[[2]]]), size=psize-1, shape=20, fill="lightGrey", color="lightGrey") +
    geom_point(data=pcs.emp, aes(x=pcs.emp[,PCs[[1]]], y=pcs.emp[,PCs[[2]]]), size=psize, shape=21, fill="lightGrey", color="#757574") +
    geom_point(data=root.data, aes(x=T1,y=T2), size=psize+2, shape=3, color="black") +
    geom_point(data=root.data, aes(x=T1,y=T2), size=psize+2, shape=4, color="black") +
#    geom_segment(data=focal.dat,aes(x=xstart,y=ystart,xend=xstop,yend=ystop),color=focal.dat$color, lwd = lsize) +
    geom_segment(data=focal.dat,aes(x=xstart,y=ystart,xend=xstop,yend=ystop),color="black", lwd = lsize) +
    geom_point(data=focal.dat, aes(x=xstart, y=ystart), size=psize, shape=21, fill="white", color="black") +
    geom_point(data=focal.dat, aes(x=xstop, y=ystop),   size=psize, shape=21, fill="white", color="black") +
    geom_point(data=tip.dat, aes(x=xstop, y=ystop),   size=psize+1, shape=21, fill="#66C2A5", color="black") +
    labs(x=colnames(pcs)[[PCs[[1]]]], y=colnames(pcs)[[PCs[[2]]]], subtitle=paste("Path to",tip.spread)) +
    theme_classic()
}

morphotrajectory(phy = agam.tree,
                 trait = LSR.anc,
                 tip.spread = "Cryptagama_aurita",
                 focus = "tip", PCs = c(1,2),
                 psize = 3, lsize = 2)

morphotrajectory(phy = agam.tree,
                 trait = LSR.anc,
                 tip.spread = "Diporiphora_superba",
                 focus = "tip", PCs = c(1,2),
                 psize = 3, lsize = 2)

morphotrajectory(phy = agam.tree,
                 trait = LSR.anc,
                 tip.spread = "Chlamydosaurus_kingii",
                 focus = "tip", PCs = c(1,2),
                 psize = 3, lsize = 2)

morphotrajectory(phy = agam.tree,
                 trait = LSR.anc,
                 tip.spread = c("Lophognathus_horneri","Amphibolurus_norrisi"),
                 focus = "clade",
                 psize = 3, lsize = 2)

morphotrajectory(phy = agam.tree,
                 trait = LSR.anc[,c("TailLength","TailWidth")],
                 tip.spread = c("Lophognathus_horneri","Amphibolurus_norrisi"),
                 focus = "clade",
                 psize = 3, lsize = 2)


# Representative Directional Trends for Optima/Clades

morphotrajectory(phy = agam.tree,
                 trait = LSR.anc[,c("BodyWidth","TailLength")],
                 tip.spread = c("Hypsilurus_godeffroyi","Hypsilurus_magnus"),
                 focus = "clade",
                 psize = 3, lsize = 2)

morphotrajectory(phy = agam.tree,
                 trait = LSR.anc[,c("SnoutEye","EyeDiameter")],
                 tip.spread = c("Lophosaurus_spinipes","Lophosaurus_boydii"),
                 focus = "clade",
                 psize = 3, lsize = 2)

morphotrajectory(phy = agam.tree,
                 trait = LSR.anc[,c("BodyWidth","TailLength")],
                 tip.spread = c("Moloch_horridus"),
                 focus = "tip",
                 psize = 3, lsize = 2)

morphotrajectory(phy = agam.tree,
                 trait = LSR.anc[,c("PelvicWidth","TailLength")],
                 tip.spread = c("Ctenophorus_maculosus","Ctenophorus_chapmani"),
                 focus = "clade",
                 psize = 3, lsize = 2)

morphotrajectory(phy = agam.tree,
                 trait = LSR.anc[,c("Size","HeadWidth")],
                 tip.spread = c("Cryptagama_aurita"),
                 focus = "tip",
                 psize = 3, lsize = 2)

morphotrajectory(phy = agam.tree,
                 trait = LSR.anc[,c("Foot","LowerLeg")],
                 tip.spread = c("Ctenophorus_gibba","Ctenophorus_infans"),
                 focus = "clade",
                 psize = 3, lsize = 2)

morphotrajectory(phy = agam.tree,
                 trait = LSR.anc[,c("Size","Neck")],
                 tip.spread = c("Chlamydosaurus_kingii"),
                 focus = "tip",
                 psize = 3, lsize = 2)

morphotrajectory(phy = agam.tree,
                 trait = LSR.anc[,c("Size","SnoutEye")],
                 tip.spread = c("Lophognathus_horneri","Amphibolurus_norrisi"),
                 focus = "clade",
                 psize = 3, lsize = 2)

morphotrajectory(phy = agam.tree,
                 trait = LSR.anc[,c("Size","Foot")],
                 tip.spread = c("Rankinia_diemensis","Pogona_barbata"),
                 focus = "clade",
                 psize = 3, lsize = 2)

morphotrajectory(phy = agam.tree,
                 trait = LSR.anc[,c("BodyWidth","Foot")],
                 tip.spread = c("Tympanocryptis_einasleighensis","Tympanocryptis_gigas"),
                 focus = "clade",
                 psize = 3, lsize = 2)

morphotrajectory(phy = agam.tree,
                 trait = LSR.anc[,c("HeadWidth","TailLength")],
                 tip.spread = c("Diporiphora_superba"),
                 focus = "tip",
                 psize = 3, lsize = 2)

morphotrajectory(phy = agam.tree,
                 trait = LSR.anc[,c("BodyWidth","Size")],
                 tip.spread = c("Diporiphora_perplexa","Diporiphora_gracilis"),
                 focus = "clade",
                 psize = 3, lsize = 2)









