
setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

source("Scripts/ClaDS_Plotting.R")

#################################################################

# Summarize the ClaDS run for the Amphibolurinae
load("Data/Amphibolurinae_ClaDS.RData")
clads.out <- CladsOutput; rm(CladsOutput)
# plot the DTT
a.DTT <- clads_DTT(clads.out, log.sp=T, return.df=F)
# plot the RTT
a.RTT <- clads_RTT(clads.out, log.rate=T, return.df=F, ci=T)
# plot the branch rates
clads_TREE(clads.out, bpal="RdYlBu", show.tip.label=T, 
           cex=0.1, log=T, tree.type="phylogram")


#################################################################

# Extract the tip speciation rates for different groups
all.rates <- data.frame(sp_rates = clads.out$lambdatip_map[2:6], group = "Hypsilurus")
all.rates <- add_row(all.rates, sp_rates = clads.out$lambdatip_map[7:9], group = "Lophosaurus")
all.rates <- add_row(all.rates, sp_rates = clads.out$lambdatip_map[13:52], group = "Ctenophorus")
all.rates <- add_row(all.rates, sp_rates = clads.out$lambdatip_map[62:70], group = "Pogona")
all.rates <- add_row(all.rates, sp_rates = clads.out$lambdatip_map[71:93], group = "Tympanocryptis")
all.rates <- add_row(all.rates, sp_rates = clads.out$lambdatip_map[94:119], group = "Diporiphora")

# Plot the tip speciation rates among genera
library(ggridges)
library(ggplot2)
ggplot(all.rates, aes(x=sp_rates, y=group, 
                      point_color=sp_rates)) +
  geom_density_ridges(jittered_points = T, scale = 2,
                      point_shape = "|", alpha = 0.75,
                      point_size = 2, size = 0.5,
                      position = position_points_jitter(height=0),
                      aes(fill = group)) +
  scale_fill_brewer(palette = "RdYlBu") +
  theme_bw() +
  theme(legend.position="none")

#################################################################

# extract the speciation rates for each branch
a.lambda <- clads_LAMBDA(clads.out)
# identify the branch type (internal or terminal)
a.lambda$edge.type <- sapply(a.lambda$child.node, function(x) ifelse(x <= Ntip(agam.tree),"terminal","internal"))
# identify the ages of the branches (when they started)
a.lambda$parent.age <- sapply(a.lambda$parent.node, function(x) max(nodeHeights(agam.tree))-nodeheight(agam.tree, x))

# plot the speciation rate against the edge length
lambda.edge <- ggplot(a.lambda, aes(x=branch.length, y=lambda, color=edge.type)) +
  geom_point() + 
  theme_classic() + 
  theme(legend.position="bottom") +
  geom_smooth(method="lm", aes(fill=edge.type))

# plot the speciation rate against the edge age
lambda.age <- ggplot(a.lambda, aes(x=parent.age, y=lambda, color=edge.type)) +
  geom_point() + 
  theme_classic() + 
  theme(legend.position="bottom") +
  geom_smooth(method="lm", aes(fill=edge.type)) +
  scale_x_reverse()

# plot the two together
lambda.edge + lambda.age + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

