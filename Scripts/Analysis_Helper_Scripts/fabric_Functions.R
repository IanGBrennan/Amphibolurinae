
# setwd("/Applications/BayesTraitsV4.1.2/Tiliquini_Fabric_Traits")

require(ape)
require(phytools)

###############################################################################

# Generate a control file for fabric analysis that tracks the ancestral trait values
# for all nodes of a tree

fabricControl <- function(phy, burnin, iterations, beta.prior, model, traits=NULL, sigma=NULL, stones=NULL, file.name="fabric_Control.txt"){
  
  # establish the control file and basic parameters
  writeLines("7", file.name) # indepnedent contrasts
  control <- file(file.name, "a")    #Open connection to append
  writeLines("2", control) # MCMC
  writeLines(model, control) # model of interest (e.g. 'Fabric', 'VarRates')
  writeLines(paste("Prior", beta.prior), control)

  # if we have traits and want to specify priors based on that information
  if(!is.null(traits)){
    # fit a constant rate BM model
    bm.fit <- geiger::fitContinuous(phy, traits, model="BM")
    # write the Alpha prior to file
    writeLines(paste("Prior", "Alpha-1", "normal", round(bm.fit$opt$z0,6), 0.75), control)
    if(sigma=="estimate"){
      writeLines(paste("Prior", "Sigma-1", "gamma", 2, round(bm.fit$opt$sigsq,6)), control) # visualize here: https://homepage.divms.uiowa.edu/~mbognar/applets/gamma.html a=2 b=ML fit
    }
  }
  if(length(sigma)==2){writeLines(paste("Prior", paste0("Sigma-",j), "uniform", sigma[1], sigma[2]), control)}

  # add 'tags' to keep track of ancestral trait values
  for(k in (Ntip(phy)+1):(Nnode(phy)+Ntip(phy))){
    tip.nos <- getDescendants(phy, k)
    tip.nos <- tip.nos[tip.nos<=Ntip(phy)]
    tip.nos <- phy$tip.label[tip.nos]
    writeLines(paste("AddTag", paste0("Tag0",k), paste(tip.nos, collapse=" ")), control)
    writeLines(paste("AddMRCA", paste0("n",k), paste0("Tag0",k)), control)
  }

  # add stepping stone sampler if you want to estimate the marginal likelihood
  if(!is.null(stones)){writeLines(paste("Stones", stones[1], stones[2]), control)}
  
  # specify the burnin, iterations, and run commands
  writeLines(paste("Burnin", burnin), control)
  writeLines(paste("Iterations", iterations), control)
  writeLines("Run", control)
  close(control)                        #Close connection
}
# NOTE: iterations and burnin must be in quotes to avoid being translated to scientific notation

# e.g. 
# phy <- read.nexus("/Applications/BayesTraitsV4.1.2/Tiliquini_Fabric_Traits/Tiliquini_BT.tre")
# prior <- "VR_LS_BL weibull 1.1 1.5"
# burnin <- "1000000"
# iterations <- "10000000"
# model <- "Fabric"

# I've removed the following prior for Alpha because it's uniform and I should be using a better shaped prior
#   range.df <- data.frame(t(apply(traits, 2, range)))
#   range.df$mean <- apply(range.df, 1, mean)
#   range.df$diff <- apply(range.df[,1:2], 1, function(x) diff(x)) # 
#   range.df$min <- range.df$mean - (range.df$diff)
#   range.df$max <- range.df$mean + (range.df$diff)
#   print(range.df)
#   for(j in 1:nrow(range.df)){
#     writeLines(paste("Prior", paste0("Alpha-",j), "uniform", round(range.df[j,"min"],6), round(range.df[j,"max"],6)), control)
#     if(sigma=="estimate" & !is.null(traits)){
#       bm.fit <- geiger::fitContinuous(phy, traits, model="BM")
#       #writeLines(paste("Prior", paste0("Sigma-",j), "uniform", 0, bm.fit$opt$sigsq*2), control) # range: 0 to twice the ML fit under constant rate BM
#       writeLines(paste("Prior", paste0("Sigma-",j), "gamma", 2, round(bm.fit$opt$sigsq,6)), control) # visualize here: https://homepage.divms.uiowa.edu/~mbognar/applets/gamma.html a=2 b=ML fit
#     }
#     if(length(sigma)==2){writeLines(paste("Prior", paste0("Sigma-",j), "uniform", sigma[1], sigma[2]), control)}


###############################################################################

# Extract the ancestral trait values for all nodes tracked in the control file

fabricAncestors <- function(log.path, phy, probs=c(0.25,0.5,0.75)){
  # read in the log
  log.in <- read.delim(log.path, header=T, skip = (48+(2*Nnode(phy))))
  # clean colnames
  colnames(log.in) <- gsub(".Alpha", "", colnames(log.in))
  # summarize the columns of interest
  trait.values <- data.frame(t(apply(log.in[,6:(5+Nnode(phy))], 2, function(x) quantile(x, probs = c(0.25, 0.5, 0.75)))))
  colnames(trait.values) <- c("lower","mean","upper")
  return(trait.values)
}
# NOTE: if this isn't working, you might have to manually adjust the number of
# lines to skip when reading in the log file

# compare results against BM
# fabric.anc <- trait.values$mean
# names(fabric.anc) <- gsub("n", "", rownames(trait.values))
# 
# il.dat <- read.csv("/Applications/BayesTraitsV4.1.2/Tiliquini_Fabric_Traits/Tiliquini_Interlimb.txt", header=F, sep="\ ")
# interlimb <- il.dat[,2]; names(interlimb) <- il.dat[,1]
# 
# phenogram(phy, interlimb, sub="anc under BM")
# phenogram(phy, c(interlimb, fabric.anc), sub="anc under Fabric")

###############################################################################

# Make a shell script to run lots of BayesTraits runs

fabric.shell <- function(in.files, tree.name){
  curr.dir <- getwd()
  #dir.create(paste0(curr.dir,"/MAFFT")) # make a directory for the new files
  
  for (k in 1:length(in.files)){
    curr.call <- paste("/Applications/BayesTraitsV4.1.2/BayesTraitsV4", 
                       paste0(getwd(),"/",tree.name), 
                       paste0(getwd(),"/",in.files[[k]]),
                       "<", 
                       paste0(getwd(),"/fabric_Control.txt"))
    write.table(curr.call, file=paste0(curr.dir,"/","fabric_Shell.sh"), append=T, row.names=F, col.names=F, quote=F)
  }
}

# Run the fabric_Shell.sh file in parallel with GNU Parallel
# $ parallel -j [#parallel] --bar :::: [shell_file.sh]

###############################################################################

# Inspect the BayesTraits log files for appropriate ESS values

inspectBT.log <- function(log.file, single.trait=T, skip.lines=43, 
  model=c("VarRates","fabric"), ESS.threshold=200){
  # establish how many lines to skip
  if(single.trait==T && is.null(skip.lines)){
      if(model=="VarRates"){skip.lines=43}
      if(model=="fabric"){skip.lines=45}
  }
  if(single.trait==F){skip.lines <- skip.lines}

  # read in the log file, skipping all the way down to the actual mcmc
  bt.log <- read.delim(log.file, skip=skip.lines)
  # remove the empty column at the end
  bt.log <- bt.log[,1:(ncol(bt.log)-1)]
  # Inspect ESS values
  ESS <- apply(bt.log, 2, coda::effectiveSize)[4:ncol(bt.log)]
  bad.ESS <- ESS[which(ESS<ESS.threshold)]
  if(length(bad.ESS)>0){
    baddies <- stringr::str_flatten(names(bad.ESS), " ")
    warning(paste("inspect variables with ESS <", ESS.threshold,":",baddies))
  }else{message(paste("all ESS >", ESS.threshold))}


  var.mean <- apply(bt.log, 2, mean)[4:ncol(bt.log)]
  if(single.trait==T){sigma.qts <- quantile(bt.log$Sigma.2.1, probs = c(0.001,0.05,0.5,0.95,0.999))}
  if(single.trait==F){
    sigmas <- dplyr::select(bt.log, starts_with("Sigma"))
    sigma.qts <- apply(sigmas, 2, function(x) quantile(x, probs=c(0.001,0.05,0.5,0.95,0.999)))
  }
  
  return(list(ESS=ESS, Mean=var.mean, Sigma_Quantiles=sigma.qts, Log=bt.log))
}

###############################################################################

# Generate a reference key to match the MD5SUM code to node numbers

fabricMD5 <- function(fpp.path, phy){
  # read in the post-processed file
  fpp <- read.csv(fpp.path)
  
  # determine the start/end columns for taxon assignments
  tax.first <- which(colnames(fpp)=="Taxa.1")
  tax.last  <- ncol(fpp)

  # 
  node.no <- NULL; node.name <- NULL
  for(k in 1:nrow(fpp)){
    curr.taxa <- unlist(c(fpp[k,tax.first:tax.last])); curr.taxa <- curr.taxa[nzchar(curr.taxa)]
    if(length(curr.taxa)==1){
      node.no <- append(node.no, which(phy$tip.label==curr.taxa))
      node.name <- append(node.name, curr.taxa[[1]])
    }
    if(length(curr.taxa)>1){
      node.no  <- append(node.no, getMRCA(phy=phy, tip=curr.taxa))
      node.name <- append(node.name, paste0("n",getMRCA(phy=phy, tip=curr.taxa)))
    }
  }

  # clean it up into a dataframe and return the object
  key <- data.frame(md5sum = fpp$Md5.Sum, node = node.no, node.name = node.name)
  return(key)
}
# fpp.path is the full path to the fabric Post Processed file (.csv usually)
# phy is the appropriate phylo object
# here, the output "node" corresponds to the MRCA, so is the same as what I call the 'child.node' elsewhere

###############################################################################

# Process the Post-Processed Fabric file to extract significant b/v values

fabricProcess <- function(fpp.path, phy, trait.name){
  # read in the post-processed file
  fpp <- read.csv(fpp.path)
  
  # determine the start/end columns for taxon assignments
  tax.first <- which(colnames(fpp)=="Taxa.1")
  tax.last  <- ncol(fpp)
  
  # filter significant beta scores (directional trends)
  sig.b <- dplyr::filter(fpp, Sig..Beta...BL. > 0)
  
  # filter significant v scores (evolvability)
  sig.v <- dplyr::filter(fpp, Sig.Scalar > 0)
  
  # summarize the significant b shifts
  sig.b.info <- NULL
  if(nrow(sig.b)>0){
    for(k in 1:nrow(sig.b)){
      curr.taxa <- unlist(c(sig.b[k,tax.first:tax.last])); curr.taxa <- curr.taxa[nzchar(curr.taxa)]
      if(length(curr.taxa)==1){node.no <- which(phy$tip.label==curr.taxa)}
      if(length(curr.taxa)>1){node.no  <- getMRCA(phy=phy, tip=curr.taxa)}
      all.desc <- phytools::getDescendants(phy, node.no); no.desc <- length(which(all.desc<=Ntip(phy)))
      curr.res <- data.frame(node = node.no,
                             scale = sig.b[k, "Mean..Beta...BL..NZ"],
                             parameter = "b",
                             trait = trait.name,
                             descendants = no.desc,
                             edge.length = phy$edge.length[which(phy$edge[,2]==node.no)])
      sig.b.info <- rbind(sig.b.info, curr.res)
    }    
  }
  
  # summarize the significant v shifts
  sig.v.info <- NULL
  if(nrow(sig.v)>0){
    for(k in 1:nrow(sig.v)){
      curr.taxa <- unlist(c(sig.v[k,tax.first:tax.last])); curr.taxa <- curr.taxa[nzchar(curr.taxa)]
      if(length(curr.taxa)==1){node.no <- which(phy$tip.label==curr.taxa)}
      if(length(curr.taxa)>1){node.no  <- getMRCA(phy=phy, tip=curr.taxa)}
      all.desc <- phytools::getDescendants(phy, node.no); no.desc <- length(which(all.desc<=Ntip(phy)))
      curr.res <- data.frame(node = node.no,
                             scale = sig.v[k, "Mean.Non.1.Scalar"],
                             parameter = "v",
                             trait = trait.name,
                             descendants = no.desc,
                             edge.length = phy$edge.length[which(phy$edge[,2]==node.no)])
      sig.v.info <- rbind(sig.v.info, curr.res)
    }
  }

  # combine the info from both shifts
  combo.info <- rbind(sig.b.info, sig.v.info)

  return(combo.info)
}
# fpp.path is the full path to the fabric Post Processed file (.csv usually)
# phy is the appropriate phylo object
# trait.name is the name of the trait to be stored in the output dataframe


###############################################################################

# Process the Merged Post-Processed Fabric file to extract significant b/v values

fabricProcessMerge <- function(fppm.path, phy, md5key, trait.name, threshold){
  # read in the post-processed file
  fppm <- read.csv(fppm.path)
  
  # filter significant beta scores (directional trends)
  sig.b <- dplyr::filter(fppm, No.Sig.Beta > threshold)
  
  # filter significant v scores (evolvability)
  sig.v <- dplyr::filter(fppm, No.Sig.Nodes > threshold)
  
  # summarize the significant b shifts
  sig.b.info <- NULL
  for(k in 1:nrow(sig.b)){
    curr.res <- data.frame(node = md5key[which(md5key$md5sum == sig.b[k,"Md5.Sum"]),"node"],
                           #scale = sig.b[k, "Mean..Beta...BL..NZ"],
                           no.sig = sig.b[k,"No.Sig.Beta"],
                           parameter = "b",
                           trait = trait.name)
    sig.b.info <- rbind(sig.b.info, curr.res)
  }
  
  # summarize the significant v shifts
  sig.v.info <- NULL
  for(k in 1:nrow(sig.v)){
    curr.res <- data.frame(node = md5key[which(md5key$md5sum == sig.v[k,"Md5.Sum"]),"node"],
                           #scale = sig.v[k, "Mean.Non.1.Scalar"],
                           no.sig = sig.v[k,"No.Sig.Nodes"],
                           parameter = "v",
                           trait = trait.name)
    sig.v.info <- rbind(sig.v.info, curr.res)
  }
  
  # combine the info from both shifts
  combo.info <- rbind(sig.b.info, sig.v.info)
  return(combo.info)
}
# fppm.path is the full path to the fabric Post Processed file (.csv usually)
# phy is the appropriate phylo object
# trait.name is the name of the trait to be stored in the output dataframe
# md5key is the output dataframe of 'fabricMD5'
# threshold is the number of times a significant effect needs to be logged to
#     be kept in the output (if you ran 6 chains and consider 5 sufficient)

###############################################################################

# plot fabric shift results quickly

fabricPlot <- function(process.obj, phy, trait.name, font.size=0.5, line.width=1, plot.type=c("fan","phylogram")){
  
  # isolate the edges with beta shifts and color them using simmap
  b.shift <- dplyr::filter(process.obj, parameter == "b")
  if(nrow(b.shift)>0){
    for (j in 1:nrow(b.shift)){
      if(b.shift$scale[j] > 0){phy <- paintBranches(tree=phy, edge=b.shift$node[j], state="increase")}
      if(b.shift$scale[j] < 0){phy <- paintBranches(tree=phy, edge=b.shift$node[j], state="decrease")}
    }
    colorz <- setNames(c("black", "#f0554d", "#4fbcf7"), c("1","increase","decrease"))
    plotSimmap(phy, colors=colorz, type=plot.type, fsize=font.size, lend=0, lwd=line.width)    
  }
  else if(nrow(b.shift)==0){plotTree(phy, type=plot.type, fsize=font.size, lend=0, lwd=line.width)}
  
  # isolate the nodes with evolvability shifts and plot them as circles
  v.shift <- dplyr::filter(process.obj, parameter == "v")
  if(nrow(v.shift)>0){
    for (j in 1:nrow(v.shift)){
      if(v.shift$scale[j] > 1){nodelabels(text="", node=noquote(v.shift$node[j]), frame="circle", bg="#f0554d", cex=0.3, col=0.5)}
      if(v.shift$scale[j] < 1){nodelabels(text="", node=noquote(v.shift$node[j]), frame="circle", bg="#4fbcf7", cex=0.3, col=0.5)}
    }
  }

  # add title to plot
  #axisPhylo()
  title(main = trait.name)
  #title(main = trait.name, sub = paste("b shift =",nrow(b.shift),":","v shift =",nrow(v.shift)))
}
# process.obj is the output of the 'fabricProcess' function
# phy is the appropriate phylo object
# trait.name is the name of the trait to be plotted as the title


###############################################################################

# plot fabric shift results in a pretty way

fabricPlotVB <- function(process.obj, phy, trait.name, font.size=0.5, line.width=1, plot.type=c("fan","phylogram"), legend=F, trait=NULL){
  require(RColorBrewer)
  require(dplyr)

  # Set the layout
  if(legend==T){layout(
    matrix(c(1,2,1,3), ncol=2, byrow=T), 
    widths=c(4,1), 
    heights=c(1,1))
  }
  
  # generate a dataframe of basic tree stats
  ndel <- data.frame(edge = 1:nrow(phy$edge),
                     parent = phy$edge[,1],
                     child  = phy$edge[,2],
                     length = phy$edge.length)
  ndel$rate <- 1
  ndel$edge.color <- "lightGrey"
  ndel$edge.lty <- 1

  # make a new tree for the purpose
  ts.tree <- phy
  
  # subset the directional trends (b)
  ts.b <- dplyr::filter(process.obj, parameter == "b")
  # process directional trends if there are any
  if(nrow(ts.b) > 0){
      ts.b$rescale <- exp(ts.b$scale)
      ts.b$edge <- sapply(ts.b$node, function(x) which(phy$edge[,2]==x))
      
      # rescale the length of trend branches by the directional trend parameter
      # NOTE: this assumes your data were natural-log transformed (not log10, and not raw)
      for(j in 1:nrow(ts.b)){
        curr.edge <- which(ts.tree$edge[,2]==ts.b$node[[j]])
        ts.tree$edge.length[[curr.edge]] <- ts.tree$edge.length[[curr.edge]] * ts.b$rescale[[j]]
        ndel$edge.lty[[curr.edge]] <- 3
      }
  }

  # subset the evolvability trends (v)
  ts.v <- dplyr::filter(process.obj, parameter == "v")
  # process the evolvability shifts if there are any
  if(nrow(ts.v) > 0){
      # rescale color of the evolvability branches by rate
      for(k in 1:nrow(ts.v)){
        curr.edge <- ts.v[k,]
        descendants <- getEdges(phy, curr.edge$node)
        ndel[which(ndel$edge %in% descendants),"rate"] <- curr.edge$scale
      }
      ndel$rate.round <- round((ndel$rate - 1)/diff(range(ndel$rate)) * 99) + 100
      
      # generate the color ramp (red faster, blue slower, grey background)
      col.ramp.red <- colorRampPalette(brewer.pal(9,"Reds")[1:8])
      col.ramp.blu <- colorRampPalette(rev(brewer.pal(9,"Blues")[3:8]))
      new.cols <- c(col.ramp.blu(99), "lightGrey", col.ramp.red(99))
      #scales::show_col(new.cols)
      ndel$edge.color <- new.cols[ndel$rate.round]
  }

  # plot the tree
  plot.phylo(ts.tree, edge.color = unlist(ndel$edge.color), edge.width=line.width, edge.lty=unlist(ndel$edge.lty), type=plot.type, cex=font.size, open.angle=5, no.margin=T)
  if(!is.null(trait)){title(trait)}

  # add colored circles to highlight branches with directional trends (red bigger, blue smaller)
  if(nrow(ts.b) > 0){
      for (j in 1:nrow(ts.b)){
        if(ts.b[j,"scale"] < 0){edgelabels(text="", edge=noquote(ts.b$edge[j]), frame="circle", bg="#4292C6", cex=0.4, col=0.5)}
        if(ts.b[j,"scale"] > 0){edgelabels(text="", edge=noquote(ts.b$edge[j]), frame="circle", bg="#EF3B2C", cex=0.4, col=0.5)}
      } 
  }
  
  # if you chose to add a legend, do that now
  min.max <- c(min(ndel$rate), max(ndel$rate))
  extreme <- max(abs(min.max))
  if(legend == T){
    plot(0,type='n',axes=FALSE,ann=FALSE);
    color.bar(new.cols, min = round(-extreme,3), max = round(extreme,3), title="v")
  }
  # return the summary dataframe 'ndel'
  return(ndel)
}


###############################################################################

# plot a summary of all shifts across multiple traits 
# will show b and v on separate phylogenies

fabricSuper <- function(all.process.obj, phy, col.palette, legend=T, no.traits=NULL, scores=c("all","up","down")){
  
  node.edge <- data.frame(edge = 1:nrow(phy$edge),
                          parent.node = phy$edge[,1],
                          child.node  = phy$edge[,2],
                          sig.b = 0, sig.v = 0)
  
  # isolate b effects and report the outcomes
  all.b <- dplyr::filter(all.process.obj, parameter == "b")
  if(scores=="up"){all.b <- dplyr::filter(all.b, scale > 0)}
  if(scores=="down"){all.b <- dplyr::filter(all.b, scale < 0)}
  message(paste("number of significant trend (b) scores < 0 :",nrow(dplyr::filter(all.b, scale < 0))))
  message(paste("number of significant trend (b) scores > 0 :",nrow(dplyr::filter(all.b, scale > 0))))
  table.b <- table(all.b$node)

  
  # isolate v effects and report the outcomes
  all.v <- dplyr::filter(all.process.obj, parameter == "v")
  if(scores=="up"){all.v <- dplyr::filter(all.v, scale > 1)}
  if(scores=="down"){all.v <- dplyr::filter(all.v, scale < 1)}
  message(paste("number of significant evolvability (v) scores < 1 :",nrow(dplyr::filter(all.v, scale < 1))))
  message(paste("number of significant evolvability (v) scores > 1 :",nrow(dplyr::filter(all.v, scale > 1))))
  table.v <- table(all.v$node)
  
  # summarize how many times an edge has had b or v shifts
  for(j in 1:nrow(node.edge)){
    if(j %in% names(table.b)){node.edge[which(node.edge$child.node==j),"sig.b"] <- table.b[which(names(table.b)==j)]}
    if(j %in% names(table.v)){node.edge[which(node.edge$child.node==j),"sig.v"] <- table.v[which(names(table.v)==j)]}
  }

  # identify how many possible shifts could occur on a single branch (# traits)
  #if(is.null(no.traits)){no.poss <- length(unique(all.process.obj$trait))}
  if(is.null(no.traits)){
    no.poss <- max(c(max(table.b),max(table.v)))
    max.b <- max(table.b); max.v <- max(table.v)
  }
  message(paste("the maximum number of estimated b shifts on a single branch is:", max.b))
  message(paste("the maximum number of estimated v shifts on a single branch is:", max.v))
  message(paste("the maximum number of shifts on a single branch is:", no.poss))
  
  # make a color ramp that suits the data
  if(col.palette %in% c("magma", "inferno", "plasma", "viridis","cividis", "rocket", "mako")){
    new.cols <- viridis::viridis(n=no.poss, option=col.palette)
  }else{
    col.ramp <- colorRampPalette(RColorBrewer::brewer.pal(9, col.palette)[2:8]) # if you need to limit upper/lower colors
    new.cols <- (col.ramp(no.poss))
    if(col.palette %in% c("Spectral","RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn","PiYG","BrBG")){new.cols <- rev(new.cols)}
  }
  new.cols <- c("black",new.cols)
  #scales::show_col(new.cols)

  node.edge$color.b <- unlist(sapply(node.edge$sig.b, function(x) new.cols[x+1]))
  node.edge$color.v <- unlist(sapply(node.edge$sig.v, function(x) new.cols[x+1]))
  node.edge$lwd.b <- sapply(node.edge$sig.b, function(x) ifelse(x==0,2,4))
  node.edge$lwd.v <- sapply(node.edge$sig.v, function(x) ifelse(x==0,2,4))
  
  par(mfrow=c(1,2))
  plot.phylo(phy, edge.color = unlist(node.edge$color.b), edge.width = node.edge$lwd.b, type = "fan", cex = 0.3)
  title("All Directional Trends (b)")
  plot.phylo(phy, edge.color = unlist(node.edge$color.v), edge.width = node.edge$lwd.v, type = "fan", cex = 0.3)
  title("All Evolvability Effects (v)")
}
# all.process.obj is a combined dataframe of multiple dataframes output from fabricProcess
# phy is the appropriate phylo object
# col.palette is a color palette, preferrably from viridis or RColorBrewer
# no.traits is the number of separate traits represented in the all.process.obj, which
#     represents the maximum number of times a single branch can be shifted
# scores indicates if you'd like to plot increasing shifts ('up'), decreasing shifts ('down'), or both ('all')

###############################################################################

fabricTraitgram <- function(process.obj, phy, trait.name, 
  lsize=2, trait=NULL, tip.spread=NULL,
  focus=c("tip","clade","all"), background.color="lightGrey"){
  require(RColorBrewer)
  require(dplyr)

  if(is.null(trait) | nrow(trait) <= Ntip(phy)){stop("trait dataframe must include ancestral trait estimates")}

  # Set the layout
#  if(legend==T){layout(
#    matrix(c(1,2,1,3), ncol=2, byrow=T), 
#    widths=c(4,1), 
#    heights=c(1,1))
#  }
  
  # build a dataframe of basic tree statistics
  ndel <- data.frame(edge = 1:nrow(phy$edge),
                     node.parent = phy$edge[,1],
                     node.child  = phy$edge[,2],
                     length = phy$edge.length,
                     timestart = sapply(phy$edge[,1], function(x) nodeheight(phy, x))-max(nodeHeights(phy)),
                     timestop  = sapply(phy$edge[,2], function(x) nodeheight(phy, x))-max(nodeHeights(phy)),
                     parent.name = paste0("n",phy$edge[,1]))
  ndel$child.name <- sapply(ndel$node.child, function(x) ifelse(x <= Ntip(phy), phy$tip[[x]], paste0("n",x)))
  ndel$rate <- 1
  
  # select and rename the parameters we're interested in
  focus.params <- dplyr::select(process.obj, scale.mean, parameter, node.parent, node)
  colnames(focus.params) <- c("scale","parameter","node.parent","node.child")
  
  # combine the dataframes to get all the info we're after
  ndel2 <- dplyr::full_join(ndel, focus.params)
  
  # provide the trait start and stop times (for the y axis)
  tdf <- dplyr::select(trait, trait.name)
  ndel2$traitstart <- sapply(ndel2$parent.name, function(x) tdf[x,])
  ndel2$traitstop  <- sapply(ndel2$child.name,  function(x) tdf[x,])
  
  # subset the evolvability trends (v)
  ts.v <- dplyr::filter(process.obj, parameter == "v")
  # process the evolvability shifts if there are any
  if(nrow(ts.v) > 0){
    # rescale color of the evolvability branches by rate
    for(k in 1:nrow(ts.v)){
      curr.edge <- ts.v[k,]
      descendants <- getEdges(phy, curr.edge$node)
      ndel2[which(ndel2$edge %in% descendants),"rate"] <- curr.edge$scale.mean
    }
    ndel2$rate.round <- round((ndel2$rate - 1)/diff(range(ndel2$rate)) * 99) + 100
    
    # generate the color ramp (red faster, blue slower, grey background)
    col.ramp.red <- colorRampPalette(RColorBrewer::brewer.pal(9,"Reds")[1:8])
    col.ramp.blu <- colorRampPalette(rev(RColorBrewer::brewer.pal(9,"Blues")[3:8]))
    new.cols <- c(col.ramp.blu(99), "lightGrey", col.ramp.red(99))
    #scales::show_col(new.cols)
    ndel2$edge.color <- new.cols[ndel2$rate.round]
  }
  
  # if you're working with a clade
  if(focus == "clade"){
    # determine the MRCA node of your target tips
    mrcn <- getMRCA(phy, tip = tip.spread)
    # get all the descendant nodes of your MRCA including internals
    all.desc <- getDescendants(phy, mrcn)    
  }
  # if you're working with a single tip
  if(focus == "tip"){all.desc <- which(phy$tip.label==tip.spread[[1]])}
  # or if you want to plot the whole tree (it will look like a mess)
  if(focus == "all"){all.desc <- 1:((Ntip(phy)+1)+Nnode(phy))}
  
  # identify which descendant nodes are tips
  tip.desc <- all.desc[which(all.desc <= Ntip(phy))]
  # get the node path from the root to each tip, and make a vector of the unique nodes
  focal.nodes <- unique(unlist(sapply(tip.desc, function(x) nodepath(phy, from=Ntip(phy)+1, to=x))))
  # subset the phydat dataframe to just the nodes of interest
  focal.dat <- dplyr::filter(ndel2, node.parent %in% focal.nodes & node.child %in% focal.nodes)
  
  # plot the traitgram with annotations
  tgram <- ggplot() +
    # plot the background lines
    {if(background.color=="lightGrey")geom_segment(data=ndel2, aes(x=timestart,y=traitstart,xend=timestop,yend=traitstop), lwd=lsize-1, color="lightGrey", alpha=0.5, lineend="round")} +
    # plot all the background (outline color) for focal branches
    geom_segment(data=focal.dat,
                 aes(x=timestart,y=traitstart,xend=timestop,yend=traitstop), lwd=lsize+1, alpha=1, 
                 color="#706f6f", lineend="round") +
    # an older version fo the evolvability effects
    geom_segment(data=focal.dat,
                 aes(x=timestart,y=traitstart,xend=timestop,yend=traitstop), lwd=lsize, alpha=1, 
                 color=focal.dat$edge.color, lineend="round") +
    # plot the directional effects (if any) as arrows
    {if(nrow(dplyr::filter(focal.dat, parameter=="b"))>0)geom_segment(data=dplyr::filter(focal.dat, parameter=="b"),
                                                                      aes(x=timestart,xend=timestop,y=traitstart,yend=traitstop),
                                                                      lineend="round",linejoin="round",size=lsize-1,arrow=arrow(length = unit(0.25, "cm")),colour="black")} +
    # plot the position of evolvability effects as a white circle
    {if(nrow(dplyr::filter(focal.dat, parameter=="v"))>0)geom_point(data=dplyr::filter(focal.dat, parameter=="v"),
                                                                    aes(x=timestop,y=traitstop), fill="white",pch=21,size=lsize+1)} +
    xlab("Million Years Ago") + ylab(trait.name) +
    theme_classic()

  # return the summary dataframe 'focal.data'
  return(tgram)
  #return(focal.dat)
}

###############################################################################

# Produce a biplot of the path from the root to a specified tip
# and annotate directional trends (arrows) and evolvability changes (colors)

fabricPathgram <- function(process.obj, phy, trait.names, 
                         lsize=2, trait=NULL, tip.spread=NULL,
                         focus=c("tip","clade","all"), background.color="lightGrey"){
  
  # downsample the obj to the focal traits
  process.obj <- dplyr::filter(process.obj, trait %in% trait.names)
  
  # build a dataframe of basic tree statistics
  ndel <- data.frame(edge = 1:nrow(phy$edge),
                     node.parent = phy$edge[,1],
                     node.child  = phy$edge[,2],
                     length = phy$edge.length,
                     timestart = sapply(phy$edge[,1], function(x) nodeheight(phy, x))-max(nodeHeights(phy)),
                     timestop  = sapply(phy$edge[,2], function(x) nodeheight(phy, x))-max(nodeHeights(phy)),
                     parent.name = paste0("n",phy$edge[,1]))
  ndel$child.name <- sapply(ndel$node.child, function(x) ifelse(x <= Ntip(phy), phy$tip[[x]], paste0("n",x)))
  ndel$rate <- 1
  ndel$edge.color <- "lightGrey"
  
  # select and rename the parameters we're interested in
  focus.params <- dplyr::select(process.obj, scale.mean, parameter, node.parent, node)
  colnames(focus.params) <- c("scale","parameter","node.parent","node.child")
  
  # combine the dataframes to get all the info we're after
  ndel2 <- dplyr::full_join(ndel, focus.params)
  
  # provide the trait start and stop times (for the y axis)
  # TRAIT 1
  t1df <- dplyr::select(trait, trait.names[[1]])
  ndel2$trait1start <- sapply(ndel2$parent.name, function(x) t1df[x,])
  ndel2$trait1stop  <- sapply(ndel2$child.name,  function(x) t1df[x,])
  # TRAIT 2
  t2df <- dplyr::select(trait, trait.names[[2]])
  ndel2$trait2start <- sapply(ndel2$parent.name, function(x) t2df[x,])
  ndel2$trait2stop  <- sapply(ndel2$child.name,  function(x) t2df[x,])
  
  # subset the evolvability trends (v)
  ts.v <- dplyr::filter(process.obj, parameter == "v")
  # process the evolvability shifts if there are any
  if(nrow(ts.v) > 0){
    # rescale color of the evolvability branches by rate
    for(k in 1:nrow(ts.v)){
      curr.edge <- ts.v[k,]
      descendants <- getEdges(phy, curr.edge$node)
      ndel2[which(ndel2$edge %in% descendants),"rate"] <- curr.edge$scale.mean
    }
    ndel2$rate.round <- round((ndel2$rate - 1)/diff(range(ndel2$rate)) * 99) + 100
    
    # generate the color ramp (red faster, blue slower, grey background)
    col.ramp.red <- colorRampPalette(RColorBrewer::brewer.pal(9,"Reds")[1:8])
    col.ramp.blu <- colorRampPalette(rev(RColorBrewer::brewer.pal(9,"Blues")[3:8]))
    new.cols <- c(col.ramp.blu(99), "lightGrey", col.ramp.red(99))
    #scales::show_col(new.cols)
    ndel2$edge.color <- new.cols[ndel2$rate.round]
  }
  
  # if you're working with a clade
  if(focus == "clade"){
    # determine the MRCA node of your target tips
    mrcn <- getMRCA(phy, tip = tip.spread)
    # get all the descendant nodes of your MRCA including internals
    all.desc <- getDescendants(phy, mrcn)    
  }
  # if you're working with a single tip
  if(focus == "tip"){all.desc <- which(phy$tip.label==tip.spread[[1]])}
  # or if you want to plot the whole tree (it will look like a mess)
  if(focus == "all"){all.desc <- 1:((Ntip(phy)+1)+Nnode(phy))}
  
  # identify which descendant nodes are tips
  tip.desc <- all.desc[which(all.desc <= Ntip(phy))]
  # get the node path from the root to each tip, and make a vector of the unique nodes
  focal.nodes <- unique(unlist(sapply(tip.desc, function(x) nodepath(phy, from=Ntip(phy)+1, to=x))))
  # subset the phydat dataframe to just the nodes of interest
  focal.dat <- dplyr::filter(ndel2, node.parent %in% focal.nodes & node.child %in% focal.nodes)
  focal.dat <- focal.dat[order(focal.dat$node.parent),]
  
  # plot the traitgram with annotations
  pathgram <- ggplot() +
    # plot the background lines
    #{if(background.color=="lightGrey")geom_segment(data=ndel2, aes(x=timestart,y=traitstart,xend=timestop,yend=traitstop), lwd=lsize-1, color="lightGrey", alpha=0.5, lineend="round")} +
    # plot all the background (outline color) for focal branches
    geom_segment(data=focal.dat,
                 aes(x=trait1start,y=trait2start,xend=trait1stop,yend=trait2stop), lwd=lsize+1, alpha=1, 
                 color="#706f6f", lineend="round") +
    # an older version fo the evolvability effects
    geom_segment(data=focal.dat,
                 aes(x=trait1start,y=trait2start,xend=trait1stop,yend=trait2stop), lwd=lsize, alpha=1, 
                 color=focal.dat$edge.color, lineend="round") +
    geom_point(data=focal.dat,
               aes(x=trait1start,y=trait2start), size=lsize+1, fill="white", pch=21) +
    # plot the directional effects (if any) as arrows
    {if(nrow(dplyr::filter(focal.dat, parameter=="b"))>0)geom_segment(data=dplyr::filter(focal.dat, parameter=="b"),
                                                                      aes(x=trait1start,xend=trait1stop,y=trait2start,yend=trait2stop),
                                                                      lineend="round",linejoin="round",size=lsize-1,arrow=arrow(length = unit(0.25, "cm")),colour="black")} +
    xlab(trait.names[[1]]) + ylab(trait.names[[2]]) +
    theme_classic()
  
  return(pathgram)
}




###############################################################################

# plot a phylogeny with the path to a given tip or clade highlighted

plot.treePath <- function(phy, focus=c("tip","clade"), tip.spread=NULL,
  lwd=1, plot.type=c("phylogram","fan"), font.size=0.5){

  # build a dataframe of basic tree statistics
  ndel <- data.frame(edge = 1:nrow(phy$edge),
                     node.parent = phy$edge[,1],
                     node.child  = phy$edge[,2],
                     length = phy$edge.length,
                     edge.color = "lightGrey",
                     edge.width = lwd)
  
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
  focal.dat <- dplyr::filter(ndel, node.parent %in% focal.nodes & node.child %in% focal.nodes)
  # specify which edges are the focus with color
  ndel[which(ndel$edge %in% focal.dat$edge),"edge.color"] <- "black"
  # specify which edges are the focus with linewidth
  ndel[which(ndel$edge %in% focal.dat$edge),"edge.width"] <- lwd+1
  
  # plot the resulting tree
  plot.phylo(phy, edge.color = unlist(ndel$edge.color), edge.width=unlist(ndel$edge.width), type=plot.type, cex=font.size, open.angle=5, no.margin=T)
}
# phy is the appropriate phylo object
# focus is either "clade" for a group of species, or "tip" for a single species
# tip.spread is a vector the summarizes the extent of the focal clade, or just a focal species name
# plot.type can be either "phylogram" or "fan" 


###############################################################################

## This function pulls out all descendant edges from either a node or edge number

getEdges <- function(phy, node=NULL, edge=NULL){
  if(is.null(edge)){node <- node}
  if(is.null(node)){node <- tree$edge[edge,2]}
  child.nodes <- phytools::getDescendants(phy, node)
  edges <- sapply(child.nodes, function(x) which(phy$edge[,2]==x))
  return(edges)
}

#########################################################################

# this function is just to plot the scale bar, which is dumb, but it works, so hey.
color.bar <- function(lut, min=0, max=100, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  #dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}
