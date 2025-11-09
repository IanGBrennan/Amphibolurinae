library(phytools)

setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

#######################################################################

load("Data/Amphibolurinae_Data.RData")
load("Data/Ancestral_Trait_Estimates.RData")

source("Scripts/Analysis_Helper_Scripts/trait.at.time.R")

#######################################################################

# combine the extant and estimated ancestral traits together

# isolate traits (exclude Genus, Genus_species)
LSR.anc <- allLSR[,1:19]
# combine observed and ancestral traits into a single dataframe
LSR.anc <- dplyr::bind_rows(LSR.anc, ancestors)

ttt <- trait.at.time.multi(timeslices = 0.25, phy = agam.tree, trait.matrix = LSR.anc, plot = F)
rownames(ttt) <- NULL

testy <- extract.volume(ttt, PCs=3, parallel=8, plot=T)

#################################################################

# THIS WAS MY SYSTEM FOR SIMULATING TRAITS UNDER REASONABLE BM MODELS
# BUT IT'S CLEAR THAT THE MODEL FIT IS POOR AND SO THE RATE ESTIMATES
# ARE POOR AND THEN THE SIMULATED TRAITS ARE NOT VERY GOOD
# INSTEAD I'M GOING TO TRY USING fastAnc TO GET THE INTERNAL NODE
# VALUES UNDER BM SO THAT THE RESULTING CONTEMPORARY VOLUME IS EQUIVALENT TO 
# THE EMPIRICAL VOLUME
# THIS MIGHT BE SOME WORK

# get evolutionary rate estimates for each trait using mvMORPH
emp.rates <- NULL; emp.theta <- NULL
for(j in 1:19){
  curr.res <- geiger::fitContinuous(agam.tree, 
                                    dplyr::select(allLSR, colnames(allLSR)[[j]]),
                                    model="BM")
  emp.rates[[j]] <- curr.res$opt$sigsq
  emp.theta[[j]] <- curr.res$opt$z0
}
emp.rates <- unlist(emp.rates); names(emp.rates) <- colnames(allLSR)[1:19]
emp.theta <- unlist(emp.theta); names(emp.theta) <- colnames(allLSR)[1:19]

# simulate traits under the empirical BM estimates
nsim <- 100; ntrait <- 19
sim.list <- lapply(1:19, function(x) fastBM(agam.tree, a=emp.theta[[x]], sig2=emp.rates[[x]], internal=T, nsim=nsim))
names(sim.list) <- names(emp.rates)

sapply(rownames(sim.list[[1]])[(Ntip(agam.tree)+1):nrow(sim.list[[1]])], function(x) paste("n",x,sep=""))

# do a stupid loop to rename the internal nodes to match the tree node.labels
for(y in 1:ntrait){rownames(sim.list[[y]])[(Ntip(agam.tree)+1):nrow(sim.list[[y]])] <- sapply(rownames(sim.list[[y]])[(Ntip(agam.tree)+1):nrow(sim.list[[y]])], function(x) paste("n",x,sep=""))}


sim.dfs <- NULL
for(k in 1:nsim){
  curr.df <- NULL
  for(x in 1:ntrait){
    curr.df <- cbind(curr.df, sim.list[[x]][,k])
  }
  colnames(curr.df) <- names(sim.list)
  sim.dfs[[k]] <- data.frame(curr.df)
}

# extract the volume-through-time for each simulated dataset
sim.vols <- NULL
for(j in 1:nsim){
  sss <- trait.at.time.multi(timeslices = 0.25, phy = agam.tree, trait.matrix = sim.dfs[[j]], plot = F)
  curr.res <- extract.volume(sss, PCs=3, parallel=8, plot=F)
  rownames(curr.res) <- NULL
  sim.vols <- cbind(sim.vols, curr.res$volume)
}
j <- 12

sim.vols2 <- sim.vols

sim.vols2 <- cbind(sim.vols2,data.frame(t(apply(sim.vols2, 1, quantile, probs=c(0.05,0.95)))))

sim.vols2$time <- curr.res$time.rev
sim.vols2$emp <- testy$volume

library(ggplot2)
ggplot() + 
  geom_ribbon(data=sim.vols2, aes(x=time, ymin=X5., ymax=X95.)) + 
  geom_line(data=sim.vols2, aes(x=time, y=emp), color="green")


rhy <- reduce.hypervolume(dplyr::filter(sss, time==31.274),3)

sim1 <- setNames(sim.list[[1]][1:119,1], rownames(sim.list[[1]])[1:119])
emp1 <- setNames(allLSR$Interlimb, rownames(allLSR))

phytools::phenogram(agam.tree, emp1, cex=0.3)
phytools::phenogram(agam.tree, sim1, cex=0.3)

dplyr::filter(sss, time==31.274)
reduce.hypervolume(dplyr::filter(ttt,time==31.274),3)
reduce.hypervolume(dplyr::filter(sss,time==31.274),3)

base.tree <- pbtree(n=100)
base.sim <- fastBM(base.tree, sig=0.002)
base.fit <- geiger::fitContinuous(base.tree, base.sim, model="BM")
sim.rate <- base.fit$opt$sigsq

# NEED TO GET fastAnc FOR EACH TRAIT
# PUT INTO A DATAFRAME
# THEN RUNIF FROM ON EACH 95%CI TO GET 100 TRAIT DATASETS

# make the empirical dataframes into separate named vectors
allLSR.vec <- apply(allLSR[,1:19], 2, setNames, nm=rownames(allLSR))
# fit fastAnc to each trait vector
all.anc <- apply(allLSR.vec, 2, fastAnc, tree=agam.tree, CI=T)
# extract the 95%CI for each trait vector 
anc.CI.df <- lapply(all.anc, function(x) data.frame(x$CI95))
# run a loop that generates trait values and fits the hypervolume for each
# simulated dataset
sim.vols <- NULL
for(j in 1:nsim){
  # random values from within 95% CIs
  anc.SIM <- data.frame(lapply(anc.CI.df, function(x) with(x,runif(nrow(x),X1,X2))))
  # name the rows
  rownames(anc.SIM) <- paste("n",120:237,sep="")
  # bind to the empirical dataset
  curr.df <- rbind(allLSR[,1:19], anc.SIM) 
  # generate the trait through time object
  curr.ttt <- trait.at.time.multi(timeslices = 0.25, phy = agam.tree, trait.matrix = curr.df, plot = F)
  # estimate the volume through time
  curr.vol <- extract.volume(curr.ttt, PCs=3, parallel=8, plot=T)
  # add the volume data to a dataframe
  sim.vols <- cbind(sim.vols, curr.vol$volume)
}

#################################################################
sim.vols <- data.frame(sim.vols)
sim.vols$time <- curr.vol$time.rev
sim.vols.all <- tidyr::pivot_longer(sim.vols, cols=starts_with("X"))

sim.vols2 <- data.frame(t(apply(sim.vols, 1, quantile, probs=c(0.05,0.95))))
sim.vols2$time <- curr.vol$time.rev
ggplot() +
  geom_ribbon(data=sim.vols2, aes(x=time,ymin=X5.,ymax=X95.)) +
  geom_line(data=testy, aes(x=time.rev,y=volume),color="green") +
  theme_classic()

ggplot() +
#  geom_ribbon(data=sim.vols2, aes(x=time,ymin=X5.,ymax=X95.)) +
  geom_line(data=sim.vols.all, aes(x=time,y=value,group=name),color="black",alpha=0.25,lwd=1) +
  geom_line(data=testy, aes(x=time.rev,y=volume),color="orange") +
  theme_classic()



#################################################################

# COULD DO DIFFERENCE OF SLOPES?

#################################################################




sim.vols <- NULL
for(j in 1:nsim){
  sss <- trait.at.time.multi(timeslices = 0.25, phy = agam.tree, trait.matrix = sim.dfs[[j]], plot = F)
  curr.res <- extract.volume(sss, PCs=3, parallel=8, plot=F)
  rownames(curr.res) <- NULL
  sim.vols <- cbind(sim.vols, curr.res$volume)
}




lapply(all.anc, function(x) with(x$CI95))


anc <- fastAnc(agam.tree, setNames(allLSR$Interlimb,rownames(allLSR)), CI=T)
cis <- data.frame(anc$CI95)
#

apply(cis, 1, runif(100,X1,X2))
with(cis,runif(nrow(cis),X1,X2))

f <- Vectorize(function(x,y) runif(1,min=x,max=y),vectorize.args = c("x","y"))
testo <- data.frame(t(with(cis,f(X1,X2))))
rownames(testo) <- paste("n",120:237,sep="")

test.dataset <- c(allLSR[1:19])
sss1 <- trait.at.time.multi(timeslices = 0.25, phy = agam.tree, trait.matrix = sim.dfs[[j]], plot = F)
curr.res <- extract.volume(sss, PCs=3, parallel=8, plot=F)


