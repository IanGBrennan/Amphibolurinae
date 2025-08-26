setwd("/Applications/BayesTraitsV4.1.2/Amphibolurinae_Fabric_Traits/fabric_Runs_Priors/Run02")
stfiles <- dir(getwd(), "Stones.txt", full.names = T)

trait.names <- sapply(sapply(stfiles, function(x) strex::str_after_last(x,"/")), function(y) strex::str_before_last(y,"_TRAIT"))

margLik.df <- NULL
for(k in 1:length(stfiles)){
 stone <- read.delim(stfiles[[k]]) 
 margLik <- as.numeric(stone[nrow(stone),])
 margLik.df <- rbind(margLik.df, data.frame(trait = trait.names[[k]],
                                margLik = margLik))
}
write.csv(margLik.df, file="/Users/ianbrennan/Documents/GitHub/Amphibolurinae/Data/BayesTraits_MarginalLikelihoods_fabric_Priors_R02.csv", row.names=F, quote=F)
