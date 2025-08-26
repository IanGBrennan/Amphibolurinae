setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

source("Scripts/fabric_Functions.R")
load("Data/Amphibolurinae_Data.RData")


############################################################################

# get names of log files for each morphological trait
log.files <- dir("/Applications/BayesTraitsV4.1.2/Amphibolurinae_Fabric_Traits", "Log.txt", full.names = T)

############################################################################

# make sure to inspect each log file for stationarity
log.res.list <- NULL
for(k in 1:length(log.files)){
  log.res <- inspectBT.log(log.file = log.files[[k]], single.trait = T, 
                           skip.lines = 284, model = "fabric", 
                           ESS.threshold = 200)
  log.res.list[[k]] <- log.res
}
trait.name <- sapply(sapply(log.files, function(y) strex::str_after_last(y,"/")), function(x) strex::str_before_first(x,"_"))
names(log.res.list) <- trait.name

# save summary of log files to RData
save(log.res.list, file = "Data/Fabric_LogESS_Summary.RData")

############################################################################

# process the fabric outputs to extract the ancestral state values
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
anc.df <- data.frame(t(data.frame(Reduce(rbind, anc.data))); colnames(anc.df) <- names(anc.data))
# remove the 'n' from node names in anc.vec
for(k in 1:length(anc.data)){names(anc.data[[k]]) <- gsub("n", "", names(anc.data[[k]]))}

############################################################################

# save to file
save(anc.data, anc.df, file = "Data/Ancestral_Trait_Estimates.RData")
