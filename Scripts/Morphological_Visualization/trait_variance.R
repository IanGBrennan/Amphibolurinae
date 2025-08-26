setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")
############################################################################

library(ggplot2)
library(ggridges)

############################################################################

source("Scripts/morphotrajectory.R")

############################################################################

load("Data/Amphibolurinae_Data.RData")
load("Data/Ancestral_Trait_Estimates.RData")

############################################################################

# transform Log-Shape Ratios back to real measurements
unLSR.anc <- exp(LSR.anc)
unLSR.anc <- unLSR.anc[,1:18] * unLSR.anc$Size

MRCA <- unLSR.anc[120,]
MRCA$Interlimb + MRCA$SnoutEye + MRCA$EyeDiameter + MRCA$PosSkull + MRCA$Neck

ozMRCA <- unLSR.anc[131,]
ozMRCA$Interlimb + ozMRCA$SnoutEye + ozMRCA$EyeDiameter + ozMRCA$PosSkull + ozMRCA$Neck



# we might want to identify which traits are the most variable 
trait.variance <- apply(allLSR[,1:19], 2, function(x) max(x) - min(x))
#trait.variance[order(trait.variance, decreasing = T)]


test.traits <- allLSR[,1:19]
#test.variance <- apply(test.traits, 2, function(x) max(x) - min(x))
#test.traits <- rbind(test.traits, test.variance)

testo <- test.traits %>%
  tidyr::pivot_longer(cols = 1:19, names_to = "trait")

testo$variance <- sapply(testo$trait, function(x) trait.variance[which(names(trait.variance)==x)])

ggplot(testo, aes(x = value, y = trait, fill = variance)) + 
  geom_density_ridges(jittered_points = TRUE, scale = 3,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_size = 1, point_alpha = 0.5, alpha = 0.9) +
  scale_fill_distiller(palette = "Spectral") +
  scale_y_discrete(limits=rev) +
  #scale_fill_viridis_c(direction=-1) +
  theme_classic()

# trait disparity
disp.stat <- function(indat){geiger::disparity(data=as.matrix(c(indat)))}
trait.disparity <- apply(allLSR[,1:19], 2, disp.stat)
trait.disparity[order(trait.disparity)]

# trait functional diversity
div.stat <- function(indat){dispRity::func.div(as.matrix(c(indat)))}
trait.funcdiv <- apply(allLSR[,1:19], 2, div.stat)

# trait functional evenness
even.stat <- function(indat){dispRity::func.eve(as.matrix(c(indat)))}
trait.funceve <- apply(allLSR[,1:19], 2, even.stat)

# Brownian Motion evolutionary Rate
fit.BM <- function(indat){geiger::fitContinuous(phy=agam.tree, dat=indat, model="BM")}
all.res <- apply(allLSR[,1:19], 2, fit.BM)
all.sig <- unlist(lapply(all.res, function(x) x$opt$sig))

stats <- data.frame(variance = trait.variance,
                    disparity = trait.disparity,
                    func_div = trait.funcdiv,
                    func_even = trait.funceve,
                    BM_rate = all.sig)

stats <- stats[sort(rownames(stats)),]
write.csv(stats, file="~/Desktop/trait_stats.csv")

# plot the difference between fabric_full and BM (delta) as a function of some trait statistics
trait.stats <- read.csv("Data/BayesTraits_MarginalLikelihoods.csv", row.names = 1)
trait.stats$delta <- trait.stats$fabric_full_prior - trait.stats$BM_noPrior
plot(trait.stats$delta ~ trait.stats$BM_rate)  
plot(trait.stats$delta ~ trait.stats$variance)
plot(trait.stats$delta ~ trait.stats$disparity)
