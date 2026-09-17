library(phytools)

setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

#######################################################################

load("Data/Amphibolurinae_Data.RData")

#######################################################################

# read in the niche breadth data
niche <- read.csv("Data/Amphibolurinae_Ecology.csv")

# select the niche breadth variable
breadth <- niche %>%
  dplyr::filter(Genus_species %in% agam.tree$tip.label) %>%
#  tibble::column_to_rownames(var="Genus_species") %>%
  dplyr::select(Genus_species, breadth)
breadth <- setNames(breadth$breadth, breadth$Genus_species)

# REMEMBER THAT CUSTOM MATRICES ARE READ AS COLUMN-->ROW

#######################################################################

# fit a series of models of niche breadth evolution

# start with the most basic, an equal rates model for all transitions
fit.ER <- fitMk(tree=agam.tree, x=breadth, model="ER"); plot(fit.ER,width=T); AIC(fit.ER)

# next allow all the transition rates to be differe
fit.ARD <- fitMk(tree=agam.tree, x=breadth, model="ARD"); plot(fit.ARD,width=T,color=T); AIC(fit.ARD)

# specify a stepwise model where transitions are allowed
# only between adjacent states, but those transition rates vary
stepwise <- matrix(c(0,1,0,0,0,
                     2,0,3,0,0,
                     0,4,0,5,0,
                     0,0,6,0,7,
                     0,0,0,8,0),5)
fit.STP <- fitMk(tree=agam.tree, x=breadth, model=stepwise); plot(fit.STP,width=T,color=T); AIC(fit.STP)

# specify a stepwise model where transitions are allowed
# only between adjacent states, but those transition rates are all the same
stp.er  <- matrix(c(0,1,0,0,0,
                     1,0,1,0,0,
                     0,1,0,1,0,
                     0,0,1,0,1,
                     0,0,0,1,0),5)
fit.STPER <- fitMk(tree=agam.tree, x=breadth, model=stp.er); plot(fit.STPER,width=T); AIC(fit.STPER)

# specify a stepwise model where transitions are allowed
# only between adjacent states, but those transition rates are all the same
stp.sym <- matrix(c(0,1,0,0,0,
                    2,0,1,0,0,
                    0,2,0,1,0,
                    0,0,2,0,1,
                    0,0,0,2,0),5)
fit.STPSYM <- fitMk(tree=agam.tree, x=breadth, model=stp.sym); plot(fit.STPSYM,width=T); AIC(fit.STPSYM)

# specify a stepwise jump model where transitiosn are allowed
# only between n+2 states (e.g. 1-->2, 1-->3), but those transitions rates are all the same
nj   <- matrix(c(0,1,1,0,0,
                 1,0,1,1,0,
                 1,1,0,1,1,
                 0,1,1,0,1,
                 0,0,1,1,0),5)
fit.NJER <- fitMk(tree=agam.tree, x=breadth, model=nj); plot(fit.NJER,width=T); AIC(fit.NJER)

# specify the stepwise jump model as above, but allowing 
# different rates for increasing and decreasing specialization/generalism
nj2  <- matrix(c(0,1,1,0,0,
                 2,0,1,1,0,
                 2,2,0,1,1,
                 0,2,2,0,1,
                 0,0,2,2,0),5)
fit.NJ <- fitMk(tree=agam.tree, x=breadth, model=nj2); plot(fit.NJ,width=T); AIC(fit.NJ)

# specify a model where transitions from generalists to specialists
# are favored over the reverse (all rates equal)
gen  <- matrix(c(0,1,0,0,0,
                 1,0,1,0,0,
                 1,1,0,1,1,
                 1,1,1,0,1,
                 1,1,1,1,0),5)
fit.GEN <- fitMk(tree=agam.tree, x=breadth, model=gen); plot(fit.GEN,width=T); AIC(fit.GEN)

# specify a model where increasing/decreasing specialization have different rates
inc  <- matrix(c(0,1,0,0,0,
                 2,0,1,0,0,
                 2,2,0,1,1,
                 2,2,2,0,1,
                 2,2,2,2,0),5)
fit.INC <- fitMk(tree=agam.tree, x=breadth, model=inc); plot(fit.INC,width=T,color=T); AIC(fit.INC)

# specify a model where transitions can only happen from generalism towards
# increasing specialization
spc  <- matrix(c(0,1,1,1,1,
                 0,0,1,1,1,
                 0,0,0,1,1,
                 0,0,0,0,1,
                 0,0,0,0,0),5)
fit.SPC <- fitMk(tree=agam.tree, x=breadth, model=spc); plot(fit.SPC,width=T,color=T); AIC(fit.SPC)

# Compare all models and save the object
anova.nb <- anova(fit.ER, fit.ARD, fit.STP, fit.STPER, fit.NJER, fit.NJ, fit.INC, fit.SPC)
anova.nb$delta <- anova.nb$AIC - min(anova.nb$AIC)
anova.nb <- anova.nb[order(anova.nb$delta, decreasing=T),]

# Save the results to file
save(fit.ER, fit.ARD, fit.STP, fit.STPER, fit.NJER, fit.NJ, fit.INC, fit.SPC,
  anova.nb, file="Data/ModellingResults_NicheBreadth.RData")

# Estimate Ancestral States under a model-averaging approach
anc.fit <- ancr(anova.nb, type="marginal")
plot(anc.fit)


# Create a likelihood ratio test nested models
phytools.LRT <- function(m1, m2){
  # m1 is complex model, m2 is simpler model
  lr.stat <- 2*(m1$logLik - m2$logLik)
  df <- length(m1$rates) - length(m2$rates)
  p.val <- pchisq(lr.stat, df=df, lower.tail=F)
  return(p.val)
}

##########################################################################

# Plot the structure of all the competing models

## create plot
par(mfrow=c(3,4))
plot(fit.ER,cex.rates=0.25,width=T)
legend("topleft",legend=paste("AIC =",round(AIC(fit.ER),1)),bty="n",cex=1)
mtext("(g) ER",line=0,adj=0,cex=1)

plot(fit.ARD,cex.rates=0.25,show.zeros=FALSE,width=T)
legend("topleft",legend=paste("AIC =",round(AIC(fit.ARD),1)),bty="n",cex=1)
mtext("(h) ARD",line=0,adj=0,cex=1)

plot(fit.STP,cex.rates=0.25,show.zeros=F,width=T,color=FALSE)
legend("topleft",legend=paste("AIC =",round(AIC(fit.STP),1)),bty="n",cex=1)
mtext("(i) STP",line=0,adj=0,cex=1)

plot(fit.STPER,cex.rates=0.25,show.zeros=F,width=T,color=FALSE)
legend("topleft",legend=paste("AIC =",round(AIC(fit.STP),1)),bty="n",cex=1)
mtext("(j) STP-ER",line=0,adj=0,cex=1)

plot(fit.NJER,cex.rates=0.25,show.zeros=FALSE,width=T,color=FALSE)
legend("topleft",legend=paste("AIC =",round(AIC(fit.NJER),1)),bty="n",cex=1)
mtext("(k) NJ-ER",line=0,adj=0,cex=1)

plot(fit.NJ,cex.rates=0.25,show.zeros=FALSE,width=T,color=FALSE)
legend("topleft",legend=paste("AIC =",round(AIC(fit.NJ),1)),bty="n",cex=1)
mtext("(l) NJ",line=0,adj=0,cex=1)

plot(fit.INC,cex.rates=0.25,show.zeros=FALSE,width=T,color=FALSE)
legend("topleft",legend=paste("AIC =",round(AIC(fit.INC),1)),bty="n",cex=1)
mtext("(m) INC",line=0,adj=0,cex=1)

plot(fit.SPC,cex.rates=0.25,show.zeros=FALSE,width=T,color=FALSE)
legend("topleft",legend=paste("AIC =",round(AIC(fit.SPC),1)),bty="n",cex=1)
mtext("(n) SPC",line=0,adj=0,cex=1)

#######################################################################

# plot the model-averaged result

# set colors
cols <- RColorBrewer::brewer.pal(9, "Reds")[c(9,7,5,3,1)]
node.cex<-apply(anc.fit$ace,1,
                function(x) if(any(x>0.7)) 0.3 else 0.8)
# plot tree
plot(anc.fit,
     args.plotTree=list(type="arc", arc_height=0.5, fsize=0.3, offset=1, color="grey"), # type="arc"
     args.nodelabels=list(cex=node.cex,piecol=cols),
     args.tiplabels=list(cex=0.2,piecol=cols),
     legend=FALSE)
# plot legend
legend(x=0, y=80,
       fit.NJ$states,pch=16,col=cols,
       horiz=T,cex=0.8,bty="n",pt.cex=2,
       x.intersp=0.5)


#######################################################################

library(phytools)

setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

#######################################################################

load("Data/Amphibolurinae_Data.RData")

#######################################################################

# read in the niche breadth data
niche <- read.csv("Data/Amphibolurinae_Ecology.csv")

# select the niche breadth variable
genspec <- niche %>%
  dplyr::filter(Genus_species %in% agam.tree$tip.label) %>%
  #  tibble::column_to_rownames(var="Genus_species") %>%
  dplyr::select(Genus_species, gen_spec)
genspec <- setNames(genspec$gen_spec, genspec$Genus_species)

#######################################################################

# start with the most basic, an equal rates model for all transitions
fit.ER <- fitMk(tree=agam.tree, x=genspec, model="ER"); plot(fit.ER,width=T)

# next allow all the transition rates to be differe
fit.ARD <- fitMk(tree=agam.tree, x=genspec, model="ARD"); plot(fit.ARD,width=T,color=T);AIC(fit.ARD)

# next allow all the transition rates to be differe
fit.SYM <- fitMk(tree=agam.tree, x=genspec, model="SYM"); plot(fit.SYM,width=T,color=T)


# specify a model where transitions from generalists to specialists
# are favored over the reverse (all rates equal)
gen  <- matrix(c(0,1,1,1,
                 1,0,0,0,
                 1,0,0,0,
                 1,0,0,0),4)
fit.GEN <- fitMk(tree=agam.tree, x=genspec, model=gen); plot(fit.GEN,width=T)

gen.sym <- matrix(c(0,1,2,3,
                    1,0,0,0,
                    2,0,0,0,
                    3,0,0,0),4)
fit.GENSYM <- fitMk(tree=agam.tree, x=genspec, model=gen.sym); plot(fit.GENSYM,width=T); AIC(fit.GENSYM)

gen.spc <- matrix(c(0,1,1,1,
                    2,0,2,2,
                    2,2,0,2,
                    2,2,2,0),4)
fit.GENSPC <- fitMk(tree=agam.tree, x=genspec, model=gen.spc); plot(fit.GENSPC,width=T,color=T); AIC(fit.GENSPC)

# The GENSPC
gen.spc2 <-matrix(c(0,1,1,1,
                    1,0,2,2,
                    1,2,0,2,
                    1,2,2,0),4)
fit.GENSPC2 <- fitMk(tree=agam.tree, x=genspec, model=gen.spc2); plot(fit.GENSPC2,width=T,color=T); AIC(fit.GENSPC2)

# gen.spc3 <-matrix(c(0,1,2,3,
#                     1,0,4,4,
#                     2,4,0,4,
#                     3,4,4,0),4)
# fit.GENSPC3 <- fitMk(tree=agam.tree, x=genspec, model=gen.spc3); plot(fit.GENSPC3,width=T,color=T); AIC(fit.GENSPC3)
# 
# gen.spc4<- matrix(c(0,2,2,2,
#                     1,0,2,2,
#                     1,2,0,2,
#                     1,2,2,0),4)
# fit.GENSPC4 <- fitMk(tree=agam.tree, x=genspec, model=gen.spc4); plot(fit.GENSPC4,width=T,color=T); AIC(fit.GENSPC4)
# 
# # this model gives rate estimates equivalent to gen.spc2, so is redundant
# gen.spc5<- matrix(c(0,3,3,3,
#                     1,0,2,2,
#                     1,2,0,2,
#                     1,2,2,0),4)
# fit.GENSPC5 <- fitMk(tree=agam.tree, x=genspec, model=gen.spc5); plot(fit.GENSPC5,width=T,color=T); AIC(fit.GENSPC5)
# 
# # this matches exactly what we see in the data but with 2 rates
# gen.cust <-matrix(c(0,1,1,1,
#                     1,0,0,0,
#                     1,0,0,0,
#                     1,2,0,0),4)
# fit.CUST <- fitMk(tree=agam.tree, x=genspec, model=gen.cust); plot(fit.CUST,width=T,color=T); AIC(fit.CUST)

# this matches exactly what we see in the data (single rate)
# only between generalists and specialists, except terrestrial can go to rock
gen.rock <-matrix(c(0,1,1,1,
                    1,0,0,0,
                    1,0,0,0,
                    1,1,0,0),4)
fit.ROCK <- fitMk(tree=agam.tree, x=genspec, model=gen.rock); plot(fit.ROCK,width=T,color=T); AIC(fit.ROCK)


# Compare all models and save the object
anova.gs <- anova(fit.ER, fit.ARD, fit.SYM, fit.GEN, fit.GENSPC2, fit.ROCK)
anova.gs$delta <- anova.gs$AIC - min(anova.gs$AIC)
anova.gs <- anova.gs[order(anova.gs$delta, decreasing=T),]

# Save the results to file
save(fit.ER, fit.ARD, fit.SYM, fit.GEN, fit.GENSPC2, fit.ROCK,
     anova.gs, file="Data/ModellingResults_GeneralistSpecialist.RData")

# Estimate Ancestral States under a model-averaging approach
anc.fit <- ancr(anova.gs, type="marginal")
plot(anc.fit)

#######################################################################

par(mfrow=c(2,3))

## create plot
plot(fit.ER,cex.rates=0.25,width=T)
legend("topleft",legend=paste("AIC =",round(AIC(fit.ER),1)),bty="n",cex=1)
mtext("(a) ER",line=0,adj=0,cex=1)

plot(fit.ARD,cex.rates=0.25,show.zeros=FALSE,width=T)
legend("topleft",legend=paste("AIC =",round(AIC(fit.ARD),1)),bty="n",cex=1)
mtext("(b) ARD",line=0,adj=0,cex=1)

plot(fit.SYM,cex.rates=0.25,show.zeros=FALSE,width=T,color=FALSE)
legend("topleft",legend=paste("AIC =",round(AIC(fit.SYM),1)),bty="n",cex=1)
mtext("(c) SYM",line=0,adj=0,cex=1)

plot(fit.GEN,cex.rates=0.25,show.zeros=FALSE,width=T,color=FALSE)
legend("topleft",legend=paste("AIC =",round(AIC(fit.GEN),1)),bty="n",cex=1)
mtext("(d) GEN",line=0,adj=0,cex=1)

plot(fit.GENSPC2,cex.rates=0.25,show.zeros=FALSE,width=T,color=FALSE)
legend("topleft",legend=paste("AIC =",round(AIC(fit.GENSPC2),1)),bty="n",cex=1)
mtext("(e) GENSPC",line=0,adj=0,cex=1)

plot(fit.ROCK,cex.rates=0.25,show.zeros=FALSE,width=T,color=FALSE)
legend("topleft",legend=paste("AIC =",round(AIC(fit.ROCK),1)),bty="n",cex=1)
mtext("(f) ROCK",line=0,adj=0,cex=1)





#######################################################################

# plot the best fitting model

# set colors
cols <- RColorBrewer::brewer.pal(9, "Spectral")[rev(c(1,3,7,9))]
node.cex<-apply(anc.fit$ace,1,
                function(x) if(any(x>0.8)) 0.3 else 0.8)
# plot tree
plot(anc.fit,
     args.plotTree=list(type="arc", arc_height=0.5, fsize=0.3, offset=1, color="grey"), # type="arc"
     args.nodelabels=list(cex=node.cex,piecol=cols),
     args.tiplabels=list(cex=0.2,piecol=cols),
     legend=FALSE)
# plot legend
legend(x=0, y=80,
       fit.GENSPC2$states,pch=16,col=cols,
       horiz=T,cex=0.8,bty="n",pt.cex=2,
       x.intersp=0.5)





#######################################################################

# NOW LET'S WORK WITH POLYMORPHIC DATA

#######################################################################

load("Data/Amphibolurinae_Data.RData")

#######################################################################

# read in the niche breadth data
niche <- read.csv("Data/Amphibolurinae_Ecology.csv")

# select the niche breadth variable
npoly <- niche %>%
  dplyr::filter(Genus_species %in% agam.tree$tip.label) %>%
  #  tibble::column_to_rownames(var="Genus_species") %>%
  dplyr::select(Genus_species, polyalph)
npoly <- setNames(npoly$polyalph, npoly$Genus_species)

# switch letters to numbers
npoly <- gsub("a",0,npoly)
npoly <- gsub("b",1,npoly)
npoly <- gsub("c",2,npoly)
npoly <- gsub("d",3,npoly)
npoly <- gsub("e",4,npoly)

###################################################################

eqr.unorder <- fitpolyMk(agam.tree, npoly, model="ER"); plot(eqr.unorder)
eqr.ordered <- fitpolyMk(agam.tree, npoly, model="ER", ordered=T); plot(eqr.ordered)
two.unorder <- fitpolyMk(agam.tree, npoly, model="transient", ordered=F); plot(two.unorder)
two.ordered <- fitpolyMk(agam.tree, npoly, model="transient", ordered=T, max.states=5,pi="fitzjohn"); plot(two.ordered)

jump.mat1 <- eqr.unorder$index.matrix
jump.mat1[1,2:5] <- 1; jump.mat1[2,3:5] <- 1; jump.mat1[3,4:5] <- 1; jump.mat1[4,5] <- 1
jump.mat1[2:5,1] <- 1; jump.mat1[3:5,2] <- 1; jump.mat1[4:5,3] <- 1; jump.mat1[5,4] <- 1
eqr.unorder.jump <- fitMk(agam.tree, to.matrix(npoly,colnames(jump.mat1)), model=jump.mat1)

jump.mat2 <- two.unorder$index.matrix
jump.mat2[1,2:5] <- 2; jump.mat2[2,3:5] <- 2; jump.mat2[3,4:5] <- 2; jump.mat2[4,5] <- 2
jump.mat2[2:5,1] <- 1; jump.mat2[3:5,2] <- 1; jump.mat2[4:5,3] <- 1; jump.mat2[5,4] <- 1
two.unorder.jump <- fitMk(agam.tree, to.matrix(npoly,colnames(jump.mat2)), model=jump.mat2)

jump.mat3 <- eqr.ordered$index.matrix
jump.mat3[1,5] <- 1; jump.mat3[5,c(1,9)] <- 1; jump.mat3[9,c(5,12)] <- 1; jump.mat3[12,c(9,14)] <- 1; jump.mat3[14,12] <- 1
eqr.order.jump <- fitMk(agam.tree, to.matrix(npoly,colnames(jump.mat3)), model=jump.mat3)

jump.mat4 <- two.ordered$index.matrix
jump.mat4[1,6]<-2; jump.mat4[6,c(1,10)]<-2; jump.mat4[10,c(6,13)]<-2; jump.mat4[13,c(10,15)]<-2 
two.order.jump <- fitMk(agam.tree, to.matrix(npoly,colnames(jump.mat4)), model=jump.mat4)

jump.mat5 <- two.ordered$index.matrix
jump.mat5[6,1] <- 2
rock.order.jump <- fitMk(agam.tree, to.matrix(npoly,colnames(jump.mat5)), model=jump.mat5)

#jump.mat6 <- two.ordered$index.matrix
#jump.mat6[1,c(6,10,13,15)] <- 2
#jump.mat6[c(10,13,15),1] <- 2
#two.order.norock.jump <- fitMk(agam.tree, to.matrix(npoly,colnames(jump.mat6)), model=jump.mat6)
#
#direct.mat <- eqr.ordered$index.matrix
#direct.mat[1:nrow(direct.mat),1:ncol(direct.mat)] <- 0
#direct.mat[1,2:5] <- 1
#direct.mat[2:5,1] <- 1
#eqr.direct <- fitMk(agam.tree, to.matrix(npoly,colnames(direct.mat)), model=direct.mat)


poly.res <- anova(eqr.unorder,
                  eqr.ordered,
                  two.unorder,
                  two.ordered,
                  eqr.unorder.jump,
                  two.unorder.jump,
                  eqr.order.jump,
                  two.order.jump,
                  rock.order.jump)
poly.res$delta <- poly.res$AIC - min(poly.res$AIC)
poly.res <- poly.res[order(poly.res$delta, decreasing=T),]

rownames(poly.res[which(poly.res$weight > 0.01),])
poly.results <- anova(two.ordered, two.order.jump, rock.order.jump)
poly.results$delta <- poly.results$AIC - min(poly.results$AIC)
anova.poly <- poly.results[order(poly.results$delta, decreasing=T),]

# Save the results to file
save(eqr.unorder, eqr.ordered, two.unorder, two.ordered,
     eqr.unorder.jump, two.unorder.jump,
     eqr.order.jump, two.order.jump,
     rock.order.jump,
     poly.res, anova.poly, file="Data/ModellingResults_Poly.RData")

# Estimate Ancestral States under a model-averaging approach
anc.fit <- ancr(anova.poly, type="marginal")
plot(anc.fit)

#######################################################################

# plot the best fitting model

# set colors
cols <- RColorBrewer::brewer.pal(9, "Spectral")[rev(c(1,3,7,9))]
node.cex<-apply(anc.fit$ace,1,
                function(x) if(any(x>0.8)) 0.3 else 0.8)
# plot tree
plot(anc.fit,
     args.plotTree=list(type="arc", arc_height=0.5, fsize=0.3, offset=1, color="grey"), # type="arc"
     args.nodelabels=list(cex=node.cex),
     args.tiplabels=list(cex=0.2),
     legend=F)

#######################################################################

par(mfrow=c(2,4))

## create plot
plot(eqr.unorder,cex.rates=0.2,width=T)
legend("topleft",legend=paste("AIC =",round(AIC(eqr.unorder),1)),bty="n",cex=1)
mtext("(a) eqr.unorder",line=0,adj=0,cex=1)

plot(two.unorder,cex.rates=0.2,show.zeros=FALSE,width=T)
legend("topleft",legend=paste("AIC =",round(AIC(two.unorder),1)),bty="n",cex=1)
mtext("(b) two.unorder",line=0,adj=0,cex=1)

plot(eqr.unorder,cex.rates=0.2,show.zeros=FALSE,width=T,color=FALSE)
legend("topleft",legend=paste("AIC =",round(AIC(eqr.unorder.jump),1)),bty="n",cex=1)
mtext("(c) eqr.unorder.jump",line=0,adj=0,cex=1)

plot(eqr.ordered,cex.rates=0.2,show.zeros=FALSE,width=T,color=FALSE)
legend("topleft",legend=paste("AIC =",round(AIC(eqr.ordered),1)),bty="n",cex=1)
mtext("(d) eqr.ordered",line=0,adj=0,cex=1)

plot(two.ordered,cex.rates=0.2,show.zeros=FALSE,width=T,color=FALSE)
legend("topleft",legend=paste("AIC =",round(AIC(two.ordered),1)),bty="n",cex=1)
mtext("(e) two.ordered",line=0,adj=0,cex=1)

plot(two.ordered,cex.rates=0.2,show.zeros=FALSE,width=T,color=FALSE)
legend("topleft",legend=paste("AIC =",round(AIC(two.order.jump),1)),bty="n",cex=1)
mtext("(f) two.order.jump",line=0,adj=0,cex=1)

plot(two.ordered,cex.rates=0.2,show.zeros=FALSE,width=T,color=FALSE)
legend("topleft",legend=paste("AIC =",round(AIC(rock.order.jump),1)),bty="n",cex=1)
mtext("(f) rock.order.jump",line=0,adj=0,cex=1)


###################################################################

# Plot the traits at the tips of the tree as split pies

plotTree(agam.tree,ftype="off",lwd=1,type="arc")
X<-strsplit(setNames(as.character(npoly),names(npoly)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(agam.tree),5,dimnames=list(agam.tree$tip.label,
                                          0:4))
pie.size <- unlist(lapply(X,function(y) length(y)))
pie.size <- pie.size[match(rownames(pies), names(pie.size))]

for(i in 1:Ntip(agam.tree)) 
  pies[agam.tree$tip.label[i],X[[agam.tree$tip.label[i]]]]<-
  rep(1/length(X[[agam.tree$tip.label[i]]]),
      length(X[[agam.tree$tip.label[i]]]))
tiplabels(pie=pies,piecol=palette()[1:5],
          cex=pie.size/6)

legend(x="topleft",legend=0:4,pt.cex=2,pch=21,
       pt.bg=palette()[1:5])



###################################################################

# Plot the traits at the tips of the tree in concentric circles

npoly.mat <- strsplit(npoly, "+", fixed = TRUE)
nm <- matrix(nrow=119,ncol=6,0)
rownames(nm) <- names(npoly.mat)
colnames(nm) <- c("0","1","2","3","4","5")
for(k in 1:length(npoly.mat)){
  if(length(npoly.mat[k][[1]])==1){nm[k,npoly.mat[k][[1]]]<-1;next}
  if(length(npoly.mat[k][[1]])>1){
    for(j in npoly.mat[k][[1]]){
      nm[k,j] <- 1
    }
  }
  
}
nm[c(79,86),6] <- 1 # add the aquatic state in for the Physignathus/Intellagama

plotFanTree.wTraits(agam.tree, nm, type="arc", part=0.5, ftype="off")

##
