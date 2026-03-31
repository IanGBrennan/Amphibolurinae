setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

library(phytools)
library(glmmTMB)
library(dplyr)

#######################################################################

load("Data/Amphibolurinae_Data.RData")
load("Data/Ancestral_Trait_Estimates.RData")

# function to get the euclidean distance between two points
euclidean <- function(a, b) sqrt(sum((a - b)^2))

#######################################################################

# combine the extant and estimated ancestral traits together

# isolate traits (exclude Genus, Genus_species)
LSR.anc <- allLSR[,1:19]
# combine observed and ancestral traits into a single dataframe
LSR.anc <- dplyr::bind_rows(LSR.anc, ancestors)

#######################################################################
# determine the geometric center of the morphospace
all.mid <- apply(LSR.anc[,1:19],2,function(x) (max(x)+min(x))/2)

# get the Euclidean distance to the middle of the morphospace
mid.dist <- apply(LSR.anc[,1:19], 1, function(x) euclidean(x,all.mid))

# get the Mahalanobis distance to the middle of the morphospace
mid.mahl <- apply(LSR.anc[,1:19], 1, function(x) mahalanobis(x=x, center=all.mid, cov=solve(cov(LSR.anc[,1:19]))))

# make a distances dataframe
distances <- data.frame(euc = mid.dist, mah = mid.mahl)

#######################################################################

# read in the niche breadth data
niche <- read.csv("Data/Amphibolurinae_Ecology.csv")

# select the niche breadth variable
breadth <- niche %>%
  dplyr::filter(Genus_species %in% agam.tree$tip.label) %>%
  #  tibble::column_to_rownames(var="Genus_species") %>%
  dplyr::select(Genus_species, breadth)
breadth <- setNames(breadth$breadth, breadth$Genus_species)

# specify a stepwise jump model where transitiosn are allowed
# only between n+2 states (e.g. 1-->2, 1-->3), but those transitions rates are all the same
nj   <- matrix(c(0,1,1,0,0,
                 1,0,1,1,0,
                 1,1,0,1,1,
                 0,1,1,0,1,
                 0,0,1,1,0),5)
fit.NJ <- phytools::fitMk(tree=agam.tree, x=breadth, model=nj); plot(fit.NJ,width=T); AIC(fit.NJ)

# extract the fit of the best model (here: NJ)
anc.fit <- phytools::ancr(fit.NJ, type="marginal")

# extract the ancestral states for substrate breadth
aces <- anc.fit$ace # this comes from 09_Modelling_NicheBreadth.R
best <- apply(aces, 1, function(x) which(x==max(x)))
names(best) <- paste0("n",names(best))
breadth.anc <- c(breadth, best)

# combine the distance estimates and the breadth data
breadth.dist <- data.frame(breadth = breadth.anc, euc = mid.dist, mah = mid.mahl)

##################################################################################
##################################################################################
##################################################################################

# WE'LL START BY MODELLING BREADTH AS A FUNCTION OF DISTANCE TO THE MORPH CENTER

# fit a GLM of substrate breadth as a function of Mahalanobis distances to the center
bm.glm <- glm(breadth ~ mah, data=breadth.dist, family="poisson"); summary(bm.glm); AIC(bm.glm)
d.weight <- seq(from=min(breadth.dist$mah), to=max(breadth.dist$mah), length.out=237)
b.weight <- predict(bm.glm, list(mah = d.weight),type="response")
plot(breadth.dist$mah, breadth.dist$breadth, pch = 21, cex = 2, xlab = "Mahalanobis distance to center", ylab = "breadth")
lines(d.weight, b.weight)

##################################################################################

# now we'll do the same but include the phylogeny via a vcv matrix
# we'll do this with glmmTMB and rerun the basic glm at the same time

# get the variance covariance matrix for the tree including ancestors
agam.mat <- phytools::vcvPhylo(agam.tree, anc.nodes=T)
rownames(agam.mat)[120:236] <- paste0("n",120:236) # rename rows
colnames(agam.mat)[120:236] <- paste0("n",120:236) # rename columns
agam.mat <- agam.mat[order(rownames(agam.mat)),] # alphabetically order rows
agam.mat <- agam.mat[,order(colnames(agam.mat))] # alphabetically order columns

# downsample the data to match the vcv matrix
bdist <- breadth.dist[which(rownames(breadth.dist)%in%rownames(agam.mat)),]
bdist$species <- rownames(bdist)
bdist$g <- 1
bdist <- bdist[order(rownames(bdist)),] # alphabetically order columns
bdist$tipnode <- ifelse(bdist$species%in%agam.tree$tip.label,"tip","node")

# let's fit a phylogenetically corrected GLM (requires data and propto to both be alphabetical order)
gmm <- glmmTMB(breadth ~ mah + (1|species) + propto(0 + species|g, agam.mat),
               data = bdist, family = "poisson",
               REML = F); summary(gmm); AIC(gmm);
bd.gmm.phy.res <- data.frame(summary(gmm)$coef$cond)

# we can plot the predictions with ggeffects
library(ggeffects)
# extract predictions
pred.df <- ggpredict(gmm, terms = "mah")
# plot
breadth.mah <- ggplot() +
  geom_ribbon(data=pred.df, aes(x=x,ymin=conf.low,ymax=conf.high),fill="lightGrey", alpha=0.5) +
  #geom_jitter(data=bdist, aes(x=mah,y=breadth, color=tipnode), size=3, alpha=0.5, height=0.1) +
  geom_line(data=pred.df, aes(x=x,y=predicted), color="darkgrey", lwd=2, lineend="round", alpha=0.5) +
  theme_bw()

# compare against the same model without phylogenetic correction
gmm3 <- glmmTMB(breadth ~ mah + (1|species),
                data = bdist, family = "poisson",
                REML = F); summary(gmm3); AIC(gmm3)
bd.gmm.nophy.res <- data.frame(summary(gmm3)$coef$cond)

# compare against an intercept-only model
gmm.int <- glmmTMB(breadth ~ 1 + (1|species),
                   data = bdist, family = "poisson",
                   REML = F); summary(gmm.int); AIC(gmm.int)
bd.gmm.int.res <- data.frame(summary(gmm.int)$coef$cond)

bd.AICs <- data.frame(phy = summary(gmm)$AICtab,
                   nonphy = summary(gmm3)$AICtab,
                   intercept = summary(gmm.int)$AICtab)

save(gmm, gmm3, gmm.int,
     bd.gmm.phy.res,
     bd.gmm.nophy.res,
     bd.gmm.int.res, 
     bd.AICs, file = "Data/GLMs/Breadth_Distance.RData")



##################################################################################
##################################################################################
##################################################################################

# NOW LET'S MODEL BREADTH AS A FUNCTION OF PC SCORES

# create a data frame of the breadth data and PC scores
lp <- data.frame(prcomp(LSR.anc[,1:19])$x)
lp$breadth <- breadth.anc[which(names(breadth.anc)==rownames(lp))]

# fit a GLM of substrate breadth as a function of PC scores
ex.glm <- glm(breadth ~ PC1 + I(PC1^2) + PC2 + I(PC2^2) + PC3 + I(PC3^2), family="poisson", data=lp); summary(ex.glm); AIC(ex.glm)
ex.int <- glm(breadth ~ 1, family="poisson", data=lp); summary(ex.int); AIC(ex.int)

# extrapolate the data to plot the fit (curves)
x1weight <- seq(from=min(lp$PC1), to=max(lp$PC1), length.out=500)
x2weight <- seq(from=min(lp$PC2), to=max(lp$PC2), length.out=500)
x3weight <- seq(from=min(lp$PC3), to=max(lp$PC3), length.out=500)
yweight <- predict(ex.glm, list(PC1 = x1weight, 
                                PC2 = x2weight,
                                PC3 = x3weight), type="response")
plot(lp$PC1, lp$breadth, pch = 16, xlab = "PC1", ylab = "breadth"); lines(x1weight, yweight)
plot(lp$PC2, lp$breadth, pch = 16, xlab = "PC2", ylab = "breadth"); lines(x2weight, yweight)
plot(lp$PC3, lp$breadth, pch = 16, xlab = "PC3", ylab = "breadth"); lines(x3weight, yweight)

pred.pc1 <- predict(ex.glm, list(PC1=x1weight))
pred.pc2 <- predict(ex.glm, PC2=x2weight)
testo <- data.frame(PC1=x1weight, pred.pc1)

##################################################################################

lps <- lp[which(rownames(lp)%in%rownames(agam.mat)),]
lps$species <- rownames(lps)
lps$g <- 1
lps <- lps[order(rownames(lps)),]

# let's fit a phylogenetically corrected GLM (requires data and propto to both be alphabetical order)
gmm2 <- glmmTMB(breadth ~ I(PC1^2) + PC2 + I(PC2^2) + PC3 + I(PC3^2) + (1|species) + propto(0 + species|g, agam.mat),
               data = lps, family = "poisson",
               REML = F); summary(gmm2); AIC(gmm2)
gmm.pc.phy.res <- data.frame(summary(gmm2)$coef$cond)

# predict the data from the model
pred.df2 <- ggpredict(gmm2, terms = c("PC1 [all]"))
pred.df3 <- ggpredict(gmm2, terms = c("PC2 [all]"))

# plot
breadth.pc <- ggplot() +
  geom_ribbon(data=pred.df2, aes(x=x,ymin=conf.low,ymax=conf.high),fill="lightGrey", alpha=0.5) +
  #geom_jitter(data=lps, aes(x=PC1,y=breadth), size=3, alpha=0.5, height=0.1) +
  geom_smooth(data=pred.df2, aes(x=x,y=predicted), method="loess", se=F) +
  theme_bw()

# let's fit a the non-phylo glmm
gmm4 <- glmmTMB(breadth ~ PC1 + I(PC1^2) + PC2 + I(PC2^2) + PC3 + I(PC3^2) + (1|species),
                data = lps, family = "poisson",
                REML = F); summary(gmm4); AIC(gmm4)
gmm.pc.nophy.res <- data.frame(summary(gmm4)$coef$cond)

# and the intercept model
gmm5 <- glmmTMB(breadth ~ 1,
                data = lps, family = "poisson",
                REML = F); summary(gmm5); AIC(gmm5)
gmm.pc.int.res <- data.frame(summary(gmm5)$coef$cond)

pc.AICs <- data.frame(phy = summary(gmm2)$AICtab,
                      nonphy = summary(gmm4)$AICtab,
                      intercept = summary(gmm5)$AICtab)

save(gmm2, gmm4, gmm5,
     gmm.pc.phy.res,
     gmm.pc.nophy.res,
     gmm.pc.int.res, 
     pc.AICs, file = "Data/GLMs/Breadth_PCs.RData")


##################################################################################

library(patchwork)

breadth.mah / breadth.pc &
  ylim(0, 4.25)

##################################################################################
##################################################################################
##################################################################################


##################################################################################

# if we want to visualize how morphologies move towards or away from the 
# center of the morphospace we can do that on a branch-by-branch basis


phy <- agam.tree
trait.df <- LSR.anc[,1:19]
metric <- "ParentToChild"
distance <- "euclidean"

# this function is just to plot the scale bar, which is dumb, but it works, so hey.
color.bar <- function(lut, min=0, max=100, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  #dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

# Function to calculate angles of a triangle given side lengths
calculate_angles <- function(a, b, c) {
  # Check if the sides form a valid triangle using the triangle inequality
  if (a + b <= c || a + c <= b || b + c <= a) {
    return("Error: The provided side lengths do not form a valid triangle.")
  }
  
  # Calculate angles using the Law of Cosines
  # cos(A) = (b^2 + c^2 - a^2) / (2 * b * c)
  # R uses acos() for inverse cosine, which returns the angle in radians
  angle_A_rad <- acos((b^2 + c^2 - a^2) / (2 * b * c))
  angle_B_rad <- acos((a^2 + c^2 - b^2) / (2 * a * c))
  
  # The third angle can be found by subtracting the first two from pi (180 degrees)
  angle_C_rad <- pi - angle_A_rad - angle_B_rad
  
  # Convert angles from radians to degrees (1 radian = 180/pi degrees)
  angle_A_deg <- angle_A_rad * 180 / pi
  angle_B_deg <- angle_B_rad * 180 / pi
  angle_C_deg <- angle_C_rad * 180 / pi
  
  # Return the angles in a named vector
  angles_deg <- c(A = angle_A_deg, B = angle_B_deg, C = angle_C_deg)
  return(angles_deg)
}

color.angles <- function(phy, trait.df, metric, distance, plot=T){
  require(RColorBrewer)
  
  # build a dataframe of basic tree statistics
  ndel <- data.frame(edge = 1:nrow(phy$edge),
                     node.parent = phy$edge[,1],
                     node.child  = phy$edge[,2],
                     length = phy$edge.length,
                     name.parent = paste0("n",phy$edge[,1]))
  ndel$name.child <- sapply(ndel$node.child, function(x) ifelse(x <= Ntip(phy), phy$tip[[x]], paste0("n",x)))
  
  # estimate the multivariate distance between each parent/child node (along each edge)
  if(metric == "ParentToChild"){
    edist <- NULL; mdist <- NULL; ndist <- NULL
    for (k in 1:nrow(ndel)){
      parent <- trait.df[which(rownames(trait.df)==ndel$name.parent[[k]]),]
      child  <- trait.df[which(rownames(trait.df)==ndel$name.child[[k]]),]
      if(distance=="euclidean"){cdist <- euclidean(parent, child); # distance from parent to child
      c2dist <- euclidean(parent, all.mid); # distance from parent to morphospace center
      c3dist <- euclidean(child, all.mid)} # distance from child to morphospace center
      if(distance=="mahalanobis"){cdist <- mahalanobis(x=unlist(child), center=unlist(parent), cov=cov(trait.df))}
      edist <- append(edist, cdist)
      mdist <- append(mdist, c2dist)
      ndist <- append(ndist, c3dist)
    }
  }
  
  
  
  ndel$edist <- edist
  ndel$mdist <- mdist
  ndel$ndist <- ndist
  ndel$movement <- mdist - ndist
  ndel$angle <- sapply(1:nrow(ndel), function(x) calculate_angles(ndel$edist[[x]],
                                                                  ndel$mdist[[x]],
                                                                  ndel$ndist[[x]])[[3]])
  
  
  ndel$color <- ndel$color <- ifelse(ndel$angle<90,"blue","red")
  ndel$change <- ifelse(ndel$angle<90,"towards","away")
  
  #ndel$color <- ndel$color <- ifelse(ndel$angle<=45,"blue","red")
  #ndel$change <- ifelse(ndel$angle<=45,"towards","away")
  
  
  # if movement is > 0 then trait change is moving towards the center
  # if movement is < 0 then trait change is moving away from the center
  
  plot(agam.tree, show.tip.label = F, edge.color = ndel$color)
  
  output <- ndel
  
  # scale and color the amount of change
  # output$movement.abs <- abs(output$movement)
  # output$distance.scaled <- round((output$movement.abs - min(output$movement.abs))/diff(range(output$movement.abs)) * 99) + 1
  
  # split the dataframes, we'll bring them back together at the end
  out.tow <- dplyr::filter(output, change=="towards")
  out.awa <- dplyr::filter(output, change=="away")
  
  out.tow$distance.scaled <- round((out.tow$angle - min(out.tow$angle))/diff(range(out.tow$angle)) * 99) + 1
  out.tow$distance.scaled <- 101 - out.tow$distance.scaled
  out.awa$distance.scaled <- round((out.awa$angle - min(out.awa$angle))/diff(range(out.awa$angle)) * 99) + 1
  
  
  #if("neither" %in% output$change){out.neither <- dplyr::filter(output, change=="neither")}
  # make appropriate color ramps
  colors.tow <- (colorRampPalette(brewer.pal(9, "Blues")[2:9])(100))
  colors.awa <- (colorRampPalette(brewer.pal(9, "Oranges")[2:9])(100))
  # get scaled colors from our ramps
  out.tow$color.scaled <- colors.tow[out.tow$distance.scaled]
  out.awa$color.scaled <- colors.awa[out.awa$distance.scaled]
  #if("neither" %in% output$change){out.neither$color.scaled <- "lightgrey"}
  # combine the dataframe back together
  #if("neither" %in% output$change){output <- rbind(out.tow, out.awa, out.neither)}
  #else{output <- rbind(out.tow, out.awa)}
  output <- rbind(out.tow, out.awa)
  # extract the colors
  color.df <- output[order(output$edge),]
  #color.df <- color.df[-1,] # the root
  # plot if you'd like
  if(plot==T){
    layout(
      matrix(c(1,2,1,3), nrow=2, ncol=2, byrow=TRUE), 
      widths=c(5,1), 
      heights=c(1,1)
    )
    plot(phy, edge.col=color.df$color.scaled, cex=0.3, edge.width=2); axisPhylo()
    color.bar(colors.tow, min=89, max=0, title="toward center")
    color.bar(colors.awa, min=90, max=round(max(out.awa$angle),2), title="away from center")
  }
  
}

color.angles(phy = agam.tree, 
             trait.df = LSR.anc[,1:19],
             metric = "ParentToChild",
             distance = "euclidean")




#####################


### Visualize the surface

# Assume fit is your model and you want to predict over ranges of x1 and x2
pred_grid <- expand.grid(
#  PC1 = seq(min(lp$PC1), max(lp$PC1), length.out = 1000),
  PC1 = seq(-2,2,length.out=100),
#  PC2 = seq(min(lp$PC2), max(lp$PC2), length.out = 1000),
  PC3 = seq(-2,2,length.out=100)
#  PC3 = seq(min(lp$PC3), max(lp$PC3), length.out = 100)
)

ex.glm2 <- glm(breadth ~ PC1 + I(PC1^2) + PC3 + I(PC3^2), family="poisson", data=lp)
pred_grid$predicted_value <- predict(ex.glm2, newdata = pred_grid)

library(ggplot2)
ggplot(pred_grid, aes(x = PC1, y = PC3, fill = predicted_value)) +
  geom_tile() +
#  scale_color_brewer(palette="Spectral") +
  scale_fill_gradient(low = "#10B1E7", high = "#D32427") + # Optional: customize colors
#  scale_fill_gradient(low = "#10B1E7", high = "white") + # Optional: customize colors
  theme_minimal()







#
