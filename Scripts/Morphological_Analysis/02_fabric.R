setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

source("Scripts/fabric_Functions.R")

############################################################################

load("Data/Amphibolurinae_Data.RData")
# if you're here just for plotting, load here:
# load("Data/Amphibolurinae_fabric_AllChainsSummary.RData")

############################################################################

# generate an md5sum key to match against nodes
# any single POST-PROCESSED (FPP) fabric file (using the same tree) will 
# be appropriate and will provide corresponding md5 and taxa lists
nodekey <- fabricMD5("/Applications/BayesTraitsV4.1.2/Amphibolurinae_Fabric_Traits/fabric_Runs_NoPriors/Run01/UpperLeg_Run01_TRAIT.txt.VarRates.txt.csv", agam.tree)

############################################################################

# process a single trait and plot the fabric results
ul.res <- fabricProcess(fpp.path = "/Applications/BayesTraitsV4.1.2/Amphibolurinae_Fabric_Traits/Run01/UpperLeg_Run01_TRAIT.txt.VarRates.txt.csv",
                        phy = agam.tree, trait.name = "Upper Leg")
fabricPlot(process.obj = ul.res, phy = agam.tree, trait.name = "Upper Leg", plot.type="phylogram")
# circles denote shifts in the evolvability (v), red increases, blue decreases
# colored branches denote directional trends (b) in traits, red increases, blue decreases

# or plot it differently:
ul.fab <- fabricPlotVB(process.obj = ul.res, phy = agam.tree, plot.type = "phylogram", 
             line.width = 1, font.size = 0.3, trait = "Upper Leg")

############################################################################

# process and plot all the individual trait files together

# identify all the fabric files
fab.files <- dir("/Applications/BayesTraitsV4.1.2/Amphibolurinae_Fabric_Traits/Run01", ".csv", full.names = T) 

# extract the trait names
trait.names <- sapply(sapply(fab.files, function(x) strex::str_after_last(x, "/")), function(y) strex::str_before_first(y, "_"))

# set the plotting frame (20 spots for 19 traits)
par(mfrow=c(5,4))

# loop through each but keep a joint output of all the b & v shifts
full.shifts <- NULL
for (j in 1:length(fab.files)){
  cur.res <- fabricProcess(fpp.path = fab.files[[j]],
                           phy = agam.tree, trait.name = trait.names[[j]])
  full.shifts <- rbind(full.shifts, cur.res)
  print(paste("plotting", trait.names[[j]]))
  fabricPlot(process.obj = cur.res, phy = agam.tree, trait.name = trait.names[[j]], font.size = 0.2, line.width =  1)
}

############################################################################

# plot the frequency and strength of trend (b) and evolvability (v) parameters
# from a single MCMC run (e.g. Run01)

plot.b <- ggplot(dplyr::filter(full.shifts, parameter == "b")) +
  geom_jitter(aes(x=trait, y=scale, color=trait), width=0.25) +
  #scale_color_continuous(palette="Greens") +
  geom_hline(yintercept=0, linetype="dotted") +
  scale_x_discrete(limits=rev) + coord_flip() +
  ylab("beta (b x t)") +
  theme_bw() + theme(legend.position = "none")

plot.v <- ggplot(dplyr::filter(full.shifts, parameter == "v")) +
  geom_jitter(aes(x=trait, y=scale, color=trait), width=0.25) +
  #scale_color_continuous(palette="Greens") +
  geom_hline(yintercept=1, linetype="dotted") +
  scale_x_discrete(limits=rev) + coord_flip() +
  ylab("v (evolvability)") +
  theme_bw() + theme(legend.position = "none")

plot.b /
  plot.v

############################################################################

# identify all the fabric directories
fab.dir <- dir("/Applications/BayesTraitsV4.1.2/Amphibolurinae_Fabric_Traits", 
               pattern="Run", include.dirs=T, full.names=T)

# loop across each directory and summarize the run information
all.chains <- NULL
for (p in 1:length(fab.dir)){
  
  # identify all the fabric files
  fab.files <- dir(fab.dir[[p]], ".csv", full.names = T) 
  
  # extract the trait names
  trait.names <- sapply(sapply(fab.files, function(x) strex::str_after_last(x, "/")), function(y) strex::str_before_first(y, "_"))
  
  # loop through each but keep a joint output of all the b & v shifts
  full.shifts <- NULL; ancestors <- NULL
  for (j in 1:length(fab.files)){
    cur.res <- fabricProcess(fpp.path = fab.files[[j]],
                             phy = agam.tree, trait.name = trait.names[[j]])
    full.shifts <- rbind(full.shifts, cur.res)
  }
  # specify which chain the result comes from
  full.shifts$chain <- p
  # combine all shifts results with other chains
  all.chains <- rbind(all.chains, full.shifts)
}

# summarize the shifts into data about individual branches
chains.summary <- NULL
unique.combo <- expand(all.chains, nesting(node, parameter, trait))
for (h in 1:nrow(unique.combo)){
  curr.combo <- dplyr::filter(all.chains, 
                       node==unique.combo$node[[h]],
                       parameter==unique.combo$parameter[[h]],
                       trait==unique.combo$trait[[h]])
  curr.res <- data.frame(node = curr.combo$node[[1]],
                         scale.mean = mean(curr.combo$scale),
                         scale.sd   = sd(curr.combo$scale),
                         parameter = curr.combo$parameter[[1]],
                         trait = curr.combo$trait[[1]],
                         #chains = paste0(min(curr.combo$chain),":",max(curr.combo$chain))
                         chains = nrow(curr.combo))
  chains.summary <- rbind(chains.summary, curr.res)
}
# each row of 'chains.summary' will correspond to one branch for one trait

# exclude shifts from the posterior that weren't found in all chains
chain.sum <- dplyr::filter(chains.summary, chains == 6)

# save these data to file
save(all.chains, chains.summary, chain.sum,
     file = "Data/Amphibolurinae_fabric_AllChainsSummary.RData")


############################################################################

# plot the frequency and strength of trend (b) and evolvability (v) parameters
# across all chains. We'll first filter it down to shifts in all chains

plot.b <- ggplot(dplyr::filter(chain.sum, parameter == "b")) +
  geom_jitter(aes(x=trait, y=scale.mean, fill=trait), 
              width=0.25, shape=21, size=3) +
  #scale_color_continuous(palette="Greens") +
  geom_hline(yintercept=0, linetype="dotted") +
  scale_x_discrete(limits=rev) + coord_flip() +
  ylab("beta (b x t)") +
  theme_bw() + theme(legend.position = "none")

plot.v <- ggplot(dplyr::filter(chain.sum, parameter == "v")) +
  geom_jitter(aes(x=trait, y=scale.mean, fill=trait), 
              width=0.25, shape=21, size=3) +
  #scale_color_continuous(palette="Greens") +
  geom_hline(yintercept=1, linetype="dotted") +
  scale_x_discrete(limits=rev) + coord_flip() +
  ylab("v (evolvability)") +
  theme_bw() + theme(legend.position = "none")

plot.b /
  plot.v

plot.b.all <- ggplot(chain.sum) +
  geom_jitter(aes(scale.mean, parameter), width=0.25, shape=21) +
  #scale_color_continuous(palette="Greens") +
  geom_hline(yintercept=0, linetype="dotted") +
  scale_x_discrete(limits=rev) + coord_flip() +
  ylab("beta (b x t)") +
  theme_bw() + theme(legend.position = "none")

ggplot(chain.sum) + 
  geom_jitter(aes(x=parameter, y=scale.mean, fill=trait,), shape=21, size=3, alpha=0.75) + 
  facet_wrap(~parameter, scales="free") +
  theme_bw()

############################################################################

# plot the summary of all b and v effects per branch
# will show how frequently a branch shows either b or v effects

# this will plot for a single chain (all traits)
fabricSuper(all.process.obj = full.shifts, phy = agam.tree, col.palette = "Spectral", scores = "all")

# this will plot for all chains combined and summarized (all traits)
colnames(chain.sum)[[2]] <- "scale"
fabricSuper(all.process.obj = chain.sum, phy = agam.tree, col.palette = "YlGnBu", scores = "all")

# downsample the chain.sum object to include only traits for which the fabric model fit best
# see Data/BayesTraits_MarginalLikelihoods.csv
chain.best <- chain.sum %>%
  dplyr::filter(., trait %in% c("BodyWidth","Foot","HeadWidth","Neck","PelvicWidth","PosSkull","Size","SnoutEye","TailLength","TailWidth")) %>%
  dplyr::filter(., parameter=="b" & abs(scale) >= 0.25 | parameter=="v")

fabricSuper(all.process.obj = chain.best, phy = agam.tree, col.palette = "YlOrRd", scores = "all")

# plot only increases in trends and rates 
fabricSuper(all.process.obj = chain.sum, phy = agam.tree, col.palette = "Blues", scores = "down")
# plot only decreases in trends and rates
fabricSuper(all.process.obj = chain.sum, phy = agam.tree, col.palette = "Reds", scores = "up")

# Gridded plot of all the traits separately
# set the plotting frame (20 spots for 19 traits)
par(mfrow=c(5,4))
# loop through each but keep a joint output of all the b & v shifts
for (j in 1:length(unique(chain.sum$trait))){
  cur.res <- dplyr::filter(chain.sum, trait == sort(unique(chain.sum$trait))[[j]])
  fabricPlot(process.obj = cur.res, phy = agam.tree, trait.name = cur.res$trait[1], font.size = 0.2, line.width =  1)
}

# set the plotting frame (20 spots for 19 traits)
par(mfrow=c(5,4))
all.traits <- sort(unique(chain.sum$trait))
for (k in 1:length(all.traits)){
  ts <- dplyr::filter(chain.sum, trait == all.traits[[k]])
  ts$scale <- ts$scale.mean
  fabricPlotVB(process.obj=ts, phy=agam.tree, plot.type="phylogram", 
               line.width=1, font.size=0.3, trait=all.traits[[k]])
}


############################################################################


# we can plot the trajectory of individual traits using traitgrams

# we can also plot traitgrams for clades of interest
fabricTraitgram(process.obj=dplyr::filter(chain.sum, trait == "TailWidth"),
                phy=agam.tree, trait.name="TailWidth", trait=LSR.anc, focus="clade", 
                tip.spread=c("Ctenophorus_badius","Ctenophorus_fordi"))

fabricTraitgram(process.obj=dplyr::filter(chain.sum, trait == "Size"),
                phy=agam.tree, trait.name="Size", trait=LSR.anc, focus="clade", 
                tip.spread=c("Pogona_barbata","Rankinia_diemensis"))

fabricTraitgram(process.obj=dplyr::filter(chain.sum, trait == "BodyWidth"),
                phy=agam.tree, trait.name="BodyWidth", trait=LSR.anc, focus="tip", 
                tip.spread=c("Moloch_horridus"))

fabricTraitgram(process.obj=dplyr::filter(chain.sum, trait == "TailLength"),
                phy=agam.tree, trait.name="TailLength", trait=LSR.anc, focus="tip", 
                tip.spread=c("Diporiphora_superba"))

fabricTraitgram(process.obj=dplyr::filter(chain.sum, trait == "UpperArm"),
                phy=agam.tree, trait.name="UpperArm", trait=LSR.anc, focus="tip", 
                tip.spread=c("Moloch_horridus"))

fabricTraitgram(process.obj=dplyr::filter(chain.sum, trait == "PelvicHeight"),
                phy=agam.tree, trait.name="PelvicHeight", trait=LSR.anc, focus="clade", 
                tip.spread=c("Ctenophorus_ornatus","Ctenophorus_slateri"))

# now imagine we might want to look at the evolution of a trait in real space (non-log shape ratio'd)

fabricTraitgram(process.obj=dplyr::filter(chain.sum, trait == "PelvicHeight"),
                phy=agam.tree, trait.name="PelvicHeight", trait=unLSR.anc, focus="clade", 
                tip.spread=c("Ctenophorus_ornatus","Ctenophorus_slateri"))

fabricTraitgram(process.obj=dplyr::filter(chain.sum, trait == "BodyWidth"),
                phy=agam.tree, trait.name="BodyWidth", trait=unLSR.anc, focus="tip", 
                tip.spread=c("Moloch_horridus"))

############################################################################

# we can also plot the path of a tip from the root to observed states in two
# dimensions, and annotate it with directional (arrows) and evolvability (color) shifts

fabricPathgram(process.obj = chain.sum,
               phy = agam.tree,
               trait.names = c("Size","TailLength"),
               lsize = 2,
               trait = LSR.anc,
               tip.spread = "Pogona_barbata",
               focus = "tip")

fabricPathgram(process.obj = chain.sum,
               phy = agam.tree,
               trait.names = c("TailWidth","TailLength"),
               lsize = 2,
               trait = LSR.anc,
               tip.spread = "Diporiphora_superba",
               focus = "tip")



############################################################################

# plot the directional effects against the edge length on which they occur
b.edge <- ggplot(dplyr::filter(chain.sum, parameter == "b")) +
  geom_point(aes(x=edge.length, y=abs(scale.mean))) + 
  geom_smooth(aes(x=edge.length, y=abs(scale.mean))) + 
  theme_bw() + theme(legend.position = "none")

# plot the evolvability effects against the number of descendant taxa
v.desc <- ggplot(dplyr::filter(chain.sum, parameter == "v")) +
  geom_point(aes(x=descendants, y=scale.mean)) + 
  geom_smooth(aes(x=descendants, y=scale.mean)) + 
  theme_bw() + theme(legend.position = "none")

# we can do the same with linear regressions if we'd like
lm.b <- lm(abs(res.b$scale.mean) ~ res.b$edge.length)
reg.b <- ggplotRegression(lm.b)

lm.v <- lm(res.v$scale.mean ~ res.v$descendants)
reg.v <- ggplotRegression(lm.v)

(b.edge | v.desc) / (reg.b | reg.v)

