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

l1ou.fit <- estimate_shift_configuration(trait.data$tree, trait.data$Y,
                                         nCores = 8, quietly = F, criterion = "pBIC",
                                         max.nShifts = 12)
# Extract the shift config and plot it
trait.shift <- shifts.to.simmap.l1ou(l1ou.fit)
plotColorSimmap(trait.shift)

#######################################################################

# Set up data and run l1ou on a subset of the trait dimensions
l1ou.fit.limit <- estimate_shift_configuration(trait.data$tree, trait.data$Y[,c(1:6,8:16,19)],
                                         nCores = 8, quietly = F, criterion = "pBIC",
                                         max.nShifts = 12)
# Extract the shift config and plot it
trait.shift.limit <- shifts.to.simmap.l1ou(l1ou.fit.limit)
plotColorSimmap(trait.shift.limit)


#######################################################################

# Set up data and run l1ou on the most important PCA axes
trait.data.pca <- adjust_data(agam.tree, prcomp(allLSR[,1:19])$x[,1:6])
# l1ou on PCA axes: 1-6=90% variance; 1-8=95% variance
l1ou.fit.pca <- estimate_shift_configuration(trait.data.pca$tree, trait.data.pca$Y,
                                             nCores = 6, quietly = F, criterion = "AICc",
                                             max.nShifts = 15)
# Extract the shift config and plot it
pca.shift <- shifts.to.simmap.l1ou(l1ou.fit.pca)
plotColorSimmap(pca.shift)

#######################################################################

# We can identify if any of the shifts are convergent
l1ou.con.pca <- estimate_convergent_regimes(l1ou.fit.pca, criterion="pBIC", nCores=6)
plot(l1ou.con.pca, plot.bar=F)
pca.conv <- shifts.to.simmap.l1ou(l1ou.con.pca)
plotColorSimmap(pca.conv)

#######################################################################

# Extract the tip states for each species for plotting
tip.states <- data.frame(Genus_species = pca.conv$tip.label,
                   state = getStates(pca.conv, type="tips"))


# combine the PCA data with the tip states to visualize the groups
pcax <- data.frame(prcomp(allLSR[,1:19])$x)
#pcax$Genus <- sampleLSR$Genus
pcax$Genus_species <- rownames(pcax)
pcax2 <- dplyr::full_join(pcax, tip.states)
pcax2$state <- letters[as.numeric(pcax2$state)]

# Plot PC 1 vs 2
PC.1.2 <- ggplot() +
  geom_point(data = transform(pcax2, state = NULL), 
             colour = "grey85", aes(x=PC1, y=PC2), size=2) +
  #geom_polygon(data=x.hull, aes(x=PC1, y=PC2, fill=genus), alpha=0.5) +
  geom_point(data=pcax2, aes(x=PC1, y=PC2, fill=state), size=3, pch=21) +
  theme_classic() + theme(legend.position = "none") + facet_wrap(~state)
# Plot PC 1 vs 3
PC.1.3 <- ggplot() +
  geom_point(data = transform(pcax2, state = NULL), 
             colour = "grey85", aes(x=PC1, y=PC3), size=2) +
  #geom_polygon(data=x.hull, aes(x=PC1, y=PC2, fill=genus), alpha=0.5) +
  geom_point(data=pcax2, aes(x=PC1, y=PC3, fill=state), size=3, pch=21) +
  theme_classic() + theme(legend.position = "none") + facet_wrap(~state)
# Plot PC 2 vs 3
PC.2.3 <- ggplot() +
  geom_point(data = transform(pcax2, state = NULL), 
             colour = "grey85", aes(x=PC2, y=PC3), size=2) +
  #geom_polygon(data=x.hull, aes(x=PC1, y=PC2, fill=genus), alpha=0.5) +
  geom_point(data=pcax2, aes(x=PC2, y=PC3, fill=state), size=3, pch=21) +
  theme_classic() + theme(legend.position = "none") + facet_wrap(~state)

PC.1.2 / PC.1.3 / PC.2.3



#######################################################################

# For comparison we can run PhylogeneticEM
# which similarly estimates multi-OU models

library(PhylogeneticEM)

# Set up data and run l1ou on all the trait dimensions
trait.data <- adjust_data(agam.tree, allLSR[,c(1:6,8:16,19)]) # 1:19
trait.data <- adjust_data(agam.tree, prcomp(allLSR[,1:19])$x[,1:6])

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
# notice there are best scores at more than one point

# isolate the results for k=4 & k=8 & k=9 & k=10 & k=...
res.k4 <- params_process(full.em, K=4)
res.k5 <- params_process(full.em, K=5)
res.k8 <- params_process(full.em, K=8)
res.k9 <- params_process(full.em, K=9)
res.k10 <- params_process(full.em, K=10)
res.kX <- params_process(full.em, K=11)

# plot the preferred result 
plot(full.em, params = res.k4)
plot(full.em, params = res.k8)
plot(full.em, params = res.k9)
plot(full.em, params = res.k10)
plot(full.em, params = res.kX)

# we can also plot equivalent shifts for any number
plot(equivalent_shifts(trait.data$tree, res.k8),
     show_shifts_values = FALSE, shifts_cex = 0.5)

# compare the loglikelihoods
log_likelihood(x=res.k4, Y_data=t(trait.data$Y[,c(1:6,8:16,19)]), phylo=trait.data$tree)
log_likelihood(x=res.k8, Y_data=t(trait.data$Y[,c(1:6,8:16,19)]), phylo=trait.data$tree)
log_likelihood(x=res.kX, Y_data=t(trait.data$Y[,c(1:6,8:16,19)]), phylo=trait.data$tree)
log_likelihood(x=res.kX, Y_data=t(trait.data$Y[,c(1:6,8:16,19)]), phylo=trait.data$tree)

log_likelihood(x=res.k4, Y_data=t(trait.data$Y), phylo=trait.data$tree)
log_likelihood(x=res.k8, Y_data=t(trait.data$Y), phylo=trait.data$tree)
log_likelihood(x=res.k10, Y_data=t(trait.data$Y), phylo=trait.data$tree)
log_likelihood(x=res.kX, Y_data=t(trait.data$Y), phylo=trait.data$tree)


# notice how the k=4 shifts are nested inside the k=8 shifts
params_process(full.em, K=4)$shifts
params_process(full.em, K=8)$shifts

k11 <- shifts_to_simmap(trait.data$tree, res.kX$shifts$edges)
plot(k11)

k8 <- shifts_to_simmap(trait.data$tree, res.k8$shifts$edges)
plot(k8)

pc.ou11 <- mvOU(k11, trait.data$Y, model="OUM")

