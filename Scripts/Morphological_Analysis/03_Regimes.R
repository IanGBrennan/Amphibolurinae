setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

#######################################################################

library(l1ou)
source("/Users/ianbrennan/Documents/GitHub/Amphibolurinae/Scripts/shifts.to.simmap.l1ou.R")
load("Data/Amphibolurinae_Data.RData")
load("Data/Ancestral_Trait_Estimates.RData")

#######################################################################

# make sure the tree is ultrametric
agam.tree <- phytools::force.ultrametric(agam.tree)

#######################################################################

# Set up data and run l1ou on all the trait dimensions
trait.data <- adjust_data(agam.tree, allLSR[,1:19]) # c(1:6,8:16,19)
trait.data <- adjust_data(agam.tree, prcomp(allLSR[,1:19])$x[,1:6])

#######################################################################

# Run PhylogeneticEM which estimates multi-OU models
# and corrects for correlations among traits

library(PhylogeneticEM)

# Set up data and run l1ou on all the trait dimensions
trait.data <- adjust_data(agam.tree, allLSR[,c(1:6,8:16,19)]) # 1:19 # the full trait set
trait.data <- adjust_data(agam.tree, prcomp(allLSR[,1:19])$x[,1:6]) # PCAs

## Run algorithm
full.em <- PhyloEM(phylo = trait.data$tree,
               #Y_data = t(trait.data$Y[,c(1:6,8:16,19)]),
               Y_data = t(trait.data$Y),
               process = "scOU",                   ## scalar OU model
               random.root = TRUE,                 ## Root is stationary (true model)
               stationary.root = TRUE,
               #independent = TRUE,
               #alpha = alpha_grid,                 ## On a grid of alpha
               nbr_alpha = 30,
               K_max = 15,                         ## Maximal number of shifts
               parallel_alpha = TRUE,              ## This can be set to TRUE for
               Ncores = 6)                         ## parallel computations
# plot the preferred result
plot(full.em)
# plot the model selection criterion
plot_criterion(full.em)

# isolate the results for k=4 & k=8 & k=9 & k=10 & k=...
res.k4 <- params_process(full.em, K=4)
res.k9 <- params_process(full.em, K=9)
res.k10 <- params_process(full.em, K=10)
res.kX <- params_process(full.em, K=11)

# plot the preferred result 
plot(full.em, params = res.k4)
plot(full.em, params = res.k9)
plot(full.em, params = res.k10)
plot(full.em, params = res.kX)

# we can also plot equivalent shifts for any number
plot(equivalent_shifts(trait.data$tree, res.k4),
     show_shifts_values = FALSE, shifts_cex = 0.5)

# compare the loglikelihoods
log_likelihood(x=res.k4, Y_data=t(trait.data$Y[,c(1:6,8:16,19)]), phylo=trait.data$tree)
log_likelihood(x=res.k9, Y_data=t(trait.data$Y[,c(1:6,8:16,19)]), phylo=trait.data$tree)
log_likelihood(x=res.kX, Y_data=t(trait.data$Y[,c(1:6,8:16,19)]), phylo=trait.data$tree)
log_likelihood(x=res.kX, Y_data=t(trait.data$Y[,c(1:6,8:16,19)]), phylo=trait.data$tree)

log_likelihood(x=res.k4, Y_data=t(trait.data$Y), phylo=trait.data$tree)
log_likelihood(x=res.k9, Y_data=t(trait.data$Y), phylo=trait.data$tree)
log_likelihood(x=res.k10, Y_data=t(trait.data$Y), phylo=trait.data$tree)
log_likelihood(x=res.kX, Y_data=t(trait.data$Y), phylo=trait.data$tree)


# notice how the k=4 shifts are nested inside the k=8 shifts
params_process(full.em, K=4)$shifts
params_process(full.em, K=11)$shifts

