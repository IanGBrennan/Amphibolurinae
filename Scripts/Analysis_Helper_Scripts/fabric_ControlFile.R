source("/Users/ianbrennan/Documents/GitHub/Amphibolurinae/Scripts/fabric_Functions.R")

setwd("/Applications/BayesTraitsV4.1.2/Amphibolurinae_Fabric_Traits")

tree <- read.tree("/Users/ianbrennan/Documents/GitHub/Amphibolurinae/Trees/Amphibolurinae_forR.tre")
#write.nexus(tree, file="/Applications/BayesTraitsV4.1.2/Amphibolurinae_Fabric_Traits/Amphibolurinae_BT.tre")

# basic control file for all individual traits (no sigma or alpha priors)
fabricControl(phy = tree,
              beta.prior = "VR_LS_BL weibull 1.1 1.5",
              burnin = "1000000",
              iterations = "10000000",
              sigma = c(0,0.25),
              model = "Fabric",
              stones = c("100","100000"),
              file.name = "fabric_Control_noPriors.txt")


# shell script for running fabric across lots of traits
trait.files <- dir(".", "_TRAIT.txt")
fabric.shell(in.files = trait.files, tree.name = "Amphibolurinae_BT.tre")


# let's make control files with trait specific sigma/alpha priors
trait.names <- colnames(allLSR[,1:19])
for(j in 1:19){
  fabricControl(phy = tree,
                beta.prior = "VR_LS_BL weibull 1.1 1.5",
                burnin = "1000000",
                iterations = "10000000",
                sigma = "estimate",
                model = "Fabric",
                stones = c("100","100000"),
                file.name = paste0("fabric_Control_",trait.names[[j]],".txt"),
                traits = dplyr::select(allLSR, trait.names[[j]]))
}



