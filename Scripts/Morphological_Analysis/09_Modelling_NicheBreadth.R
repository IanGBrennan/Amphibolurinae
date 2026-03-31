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
fit.NJ <- fitMk(tree=agam.tree, x=breadth, model=nj); plot(fit.NJ,width=T); AIC(fit.NJ)

# specify the stepwise jump model as above, but allowing 
# different rates for increasing and decreasing specialization/generalism
nj2  <- matrix(c(0,1,1,0,0,
                 2,0,1,1,0,
                 2,2,0,1,1,
                 0,2,2,0,1,
                 0,0,2,2,0),5)
fit.NJ2 <- fitMk(tree=agam.tree, x=breadth, model=nj2); plot(fit.NJ2,width=T); AIC(fit.NJ2)

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



# compare AIC values across model fits
AIC.df <- data.frame(AIC(fit.ER, fit.ARD, fit.STP, fit.STPER, fit.NJ, fit.NJ2, fit.INC, fit.SPC))
# calculate the AIC weights
AIC.df$W <- aic.w(AIC(fit.ER, fit.ARD, fit.STP, fit.STPER, fit.NJ, fit.NJ2, fit.INC, fit.SPC)$AIC)
# order models by AICw
AIC.df[order(AIC.df$W,decreasing=T),]

# extract the fit of the best model (here: NJ)
anc.fit <- ancr(fit.NJ, type="marginal")


## create plot
par(mfrow=c(3,4))
plot(fit.ER,cex.rates=0.9,width=T)
legend("topleft",legend=paste("AIC =",round(AIC(fit.ER),1)),bty="n",cex=1)
mtext("(e) ER",line=0,adj=0,cex=1)

plot(fit.ARD,cex.rates=0.9,show.zeros=FALSE,width=T)
legend("topleft",legend=paste("AIC =",round(AIC(fit.ARD),1)),bty="n",cex=1)
mtext("(f) ARD",line=0,adj=0,cex=1)

plot(fit.STP,cex.rates=0.9,show.zeros=F,width=T,color=FALSE)
legend("topleft",legend=paste("AIC =",round(AIC(fit.STP),1)),bty="n",cex=1)
mtext("(g) STP",line=0,adj=0,cex=1)

plot(fit.STPER,cex.rates=0.9,show.zeros=F,width=T,color=FALSE)
legend("topleft",legend=paste("AIC =",round(AIC(fit.STP),1)),bty="n",cex=1)
mtext("(h) STP-ER",line=0,adj=0,cex=1)

plot(fit.NJ,cex.rates=0.9,show.zeros=FALSE,width=T,color=FALSE)
legend("topleft",legend=paste("AIC =",round(AIC(fit.NJ2),1)),bty="n",cex=1)
mtext("(i) NJ-ER",line=0,adj=0,cex=1)

plot(fit.NJ2,cex.rates=0.9,show.zeros=FALSE,width=T,color=FALSE)
legend("topleft",legend=paste("AIC =",round(AIC(fit.NJ),1)),bty="n",cex=1)
mtext("(j) NJ",line=0,adj=0,cex=1)

plot(fit.INC,cex.rates=0.9,show.zeros=FALSE,width=T,color=FALSE)
legend("topleft",legend=paste("AIC =",round(AIC(fit.INC),1)),bty="n",cex=1)
mtext("(k) INC",line=0,adj=0,cex=1)

plot(fit.SPC,cex.rates=0.9,show.zeros=FALSE,width=T,color=FALSE)
legend("topleft",legend=paste("AIC =",round(AIC(fit.SPC),1)),bty="n",cex=1)
mtext("(l) SPC",line=0,adj=0,cex=1)

#######################################################################

# plot the best fitting model

# set colors
cols <- RColorBrewer::brewer.pal(9, "Reds")[c(9,7,5,3,1)]
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

gen.spc2 <-matrix(c(0,1,1,1,
                    1,0,2,2,
                    1,2,0,2,
                    1,2,2,0),4)
fit.GENSPC2 <- fitMk(tree=agam.tree, x=genspec, model=gen.spc2); plot(fit.GENSPC2,width=T,color=T); AIC(fit.GENSPC2)

gen.spc3 <-matrix(c(0,1,2,3,
                    1,0,4,4,
                    2,4,0,4,
                    3,4,4,0),4)
fit.GENSPC3 <- fitMk(tree=agam.tree, x=genspec, model=gen.spc3); plot(fit.GENSPC3,width=T,color=T); AIC(fit.GENSPC3)

gen.spc4<- matrix(c(0,2,2,2,
                    1,0,2,2,
                    1,2,0,2,
                    1,2,2,0),4)
fit.GENSPC4 <- fitMk(tree=agam.tree, x=genspec, model=gen.spc4); plot(fit.GENSPC4,width=T,color=T); AIC(fit.GENSPC4)

# this model gives rate estimates equivalent to gen.spc2, so is redundant
gen.spc5<- matrix(c(0,3,3,3,
                    1,0,2,2,
                    1,2,0,2,
                    1,2,2,0),4)
fit.GENSPC5 <- fitMk(tree=agam.tree, x=genspec, model=gen.spc5); plot(fit.GENSPC5,width=T,color=T); AIC(fit.GENSPC5)

# this matches exactly what we see in the data
gen.cust <-matrix(c(0,1,1,1,
                    1,0,0,0,
                    1,0,0,0,
                    1,2,0,0),4)
fit.CUST <- fitMk(tree=agam.tree, x=genspec, model=gen.cust); plot(fit.CUST,width=T,color=T); AIC(fit.CUST)

AIC.df <- AIC(fit.ER, fit.ARD, fit.SYM, fit.GEN, fit.GENSPC2)
#AIC.df <- AIC(fit.ER, fit.SYM, fit.GEN, fit.GENSYM, fit.GENSPC, fit.GENSPC2, fit.GENSPC3, fit.GENSPC4, fit.GENSPC5)
AIC.df$W <- aic.w(AIC.df$AIC)
AIC.df[order(AIC.df$W,decreasing=T),]

### THIS IS EXACTLY WHAT I NEEDED!
# extract the fit of the best model (here: NJ)
anc.fit <- ancr(fit.GENSPC2, type="marginal")
plot(anc.fit)

#######################################################################

## create plot
plot(fit.ER,cex.rates=0.9,width=T)
legend("topleft",legend=paste("AIC =",round(AIC(fit.ER),1)),bty="n",cex=1)
mtext("(a) ER",line=0,adj=0,cex=1)

plot(fit.ARD,cex.rates=0.9,show.zeros=FALSE,width=T)
legend("topleft",legend=paste("AIC =",round(AIC(fit.ARD),1)),bty="n",cex=1)
mtext("(b) ARD",line=0,adj=0,cex=1)

plot(fit.SYM,cex.rates=0.9,show.zeros=FALSE,width=T,color=FALSE)
legend("topleft",legend=paste("AIC =",round(AIC(fit.SYM),1)),bty="n",cex=1)
mtext("(c) SYM",line=0,adj=0,cex=1)

plot(fit.GENSPC2,cex.rates=0.9,show.zeros=FALSE,width=T,color=FALSE)
legend("topleft",legend=paste("AIC =",round(AIC(fit.GENSPC2),1)),bty="n",cex=1)
mtext("(d) GENSPC",line=0,adj=0,cex=1)

#######################################################################



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




load("Data/Amphibolurinae_Data.RData")

# read in the niche breadth data
niche <- read.csv("Data/Amphibolurinae_Ecology.csv")

# select the niche breadth variable
poly.state <- niche %>%
  dplyr::filter(Genus_species %in% agam.tree$tip.label) %>%
  #  tibble::column_to_rownames(var="Genus_species") %>%
  dplyr::select(Genus_species, poly)
poly.state <- setNames(poly.state$poly, poly.state$Genus_species)

# start with the most basic, an equal rates model for all transitions
fit.ER <- fitpolyMk(tree=agam.tree, x=poly.state, model="ER"); plot(fit.ER,width=T,zeros=T)

simmap.trees <- make.simmap(agam.tree,fit.ER$data,model=fit.ER$index.matrix,nsim=10)

cols<-setNames(colorRampPalette(c("blue","grey","brown","yellow","green","ForestGreen"))(11),
               colnames(fit$data))

cols <- setNames(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(63),colnames(fit.ER$data))
plot(summary(simmap.trees),colors=cols,type="fan",ftype="off")
legend(x="topleft",legend=colnames(fit.ER$data),pt.cex=2.4,pch=21,
       pt.bg=cols)

spc.trees<-mergeMappedStates(simmap.trees,fit.ER$states[1:21],"specialist")
spc.trees<-mergeMappedStates(spc.trees,fit.ER$states[22:63],"generalist")
cols<-setNames(c("#f06f24","#5acadb"),c("generalist","specialist"))
plot(summary(spc.trees),colors=cols,type="fan",ftype="off")
legend(x="topleft",legend=names(cols),pt.cex=2.4,pch=21,
       pt.bg=cols)


spc2spc <- fit.ER$index.matrix
spc2spc[1:63,1:63] <- 1
spc2spc[1:6, 1:6]  <- 2
diag(spc2spc) <- 0


testo <- fit.TRA$index.matrix
testa <- fitpolyMk2(tree=agam.tree, x=poly.state, model=spc2spc)

# start with the most basic, an equal rates model for all transitions
fit.SPC2SPC <- fitpolyMk(tree=agam.tree, x=poly.state, model="ER", max.poly=5); plot(fit.SPC2SPC,width=T,zeros=T)


fit.TRA <- fitpolyMk(tree=agam.tree, x=poly.state, model="transient"); plot(fit.TRA,width=T)

tra2 <- fit.TRA$index.matrix
tra2[1:6,1:6] <- 3
diag(tra2) <- 0

fit.TRA2 <- fitpolyMk2(tree=agam.tree, x=poly.state, model=tra2); plot(fit.TRA2,width=T)


# next allow all the transition rates to be differe
fit.ARD <- fitMk(tree=agam.tree, x=poly.state, model="ARD"); plot(fit.ARD,width=T,color=T)

# specify a stepwise model where transitions are allowed
# only between adjacent states, but those transition rates vary
stepwise <- matrix(c(0,1,0,0,0,
                     2,0,3,0,0,
                     0,4,0,5,0,
                     0,0,6,0,7,
                     0,0,0,8,0),5)
fit.STP <- fitMk(tree=agam.tree, x=breadth, model=stepwise); plot(fit.STP,width=T,color=T)



plotTree(agam.tree,ftype="off",lwd=1,type="fan")
X<-strsplit(setNames(as.character(poly.state),names(poly.state)),"+",fixed=TRUE)
pies<-matrix(0,Ntip(agam.tree),6,dimnames=list(agam.tree$tip.label,
                                               c("w","r","g","l","h","t")))
for(i in 1:Ntip(agam.tree)) 
  pies[agam.tree$tip.label[i],X[[agam.tree$tip.label[i]]]]<-
  rep(1/length(X[[agam.tree$tip.label[i]]]),
      length(X[[agam.tree$tip.label[i]]]))
#tiplabels(pie=pies,piecol=c("black","yellow","red","blue","green","pink"),cex=0.35)
tiplabels(pie=pies,piecol=c("blue","grey","brown","yellow","green","ForestGreen"),cex=0.5)
legend(x="topleft",legend=c("water","rock","ground","low","high","tree"),pt.cex=2,pch=21,
       pt.bg=c("blue","grey","brown","yellow","green","ForestGreen"))
