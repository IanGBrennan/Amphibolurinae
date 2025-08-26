setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

source("Scripts/fabric_Functions.R")
load("Data/Amphibolurinae_Data.RData")


############################################################################

# identify all the fabric directories
fab.dir <- dir("/Applications/BayesTraitsV4.1.2/Amphibolurinae_Fabric_Traits/fabric_Runs_Priors", 
               pattern="Run", include.dirs=T, full.names=T)

############################################################################

all.ancestors <- NULL
for (p in 1:length(fab.dir)){
  
  # get names of log files for each morphological trait
  log.files <- dir(fab.dir[[p]], "Log.txt", full.names = T)
  
  ############################################################################
  
#  # make sure to inspect each log file for stationarity
#  log.res.list <- NULL
#  for(k in 1:length(log.files)){
#    log.res <- inspectBT.log(log.file = log.files[[k]], single.trait = T, 
#                             skip.lines = 284, model = "fabric", 
#                             ESS.threshold = 200)
#    log.res.list[[k]] <- log.res
#  }
#  trait.name <- sapply(sapply(log.files, function(y) strex::str_after_last(y,"/")), function(x) strex::str_before_first(x,"_"))
#  names(log.res.list) <- trait.name
#  
#  # save summary of log files to RData
#  save(log.res.list, file = "Data/Fabric_LogESS_Summary.RData")
  
  ############################################################################
  
  # process the fabric outputs to extract the ancestral state values
  message(paste("processing chain", p))
  anc.data <- NULL
  for(k in 1:length(log.files)){
    anc.res <- fabricAncestors(log.path = log.files[[k]], phy = agam.tree)
    anc.vec <- anc.res$mean; names(anc.vec) <- rownames(anc.res)
    #  names(anc.vec) <- gsub("n", "", names(anc.vec))
    anc.data[[k]] <- anc.vec
  }
  # get the trait names and apply to the list
  names(anc.data) <- sapply(log.files, function(x) strex::str_before_first(strex::str_after_last(x, "/"), "_"))
  # convert these data into a dataframe in case that's useful
  anc.df <- data.frame(t(data.frame(Reduce(rbind, anc.data)))); 
  colnames(anc.df) <- names(anc.data)
  anc.df$chain <- p
  anc.df$node <- rownames(anc.df)
  
  # add the current results to the results across chains
  all.ancestors <- rbind(all.ancestors, anc.df)
}

# Summarize results from across chains 
unique.nodes <- unique(all.ancestors$node)
ancestors <- NULL
for (t in 1:length(unique.nodes)){
  curr.node <- dplyr::filter(all.ancestors, node == unique.nodes[[t]])
  curr.ancs <- apply(curr.node[1:19], 2, mean)
  ancestors <- rbind(ancestors, curr.ancs)
}
ancestors <- data.frame(ancestors)
rownames(ancestors) <- unique.nodes

############################################################################

# save to file
save(ancestors, file = "Data/Ancestral_Trait_Estimates.RData")
