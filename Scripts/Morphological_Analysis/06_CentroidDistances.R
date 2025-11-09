
setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

#######################################################################

load("Data/Amphibolurinae_Data.RData")
load("Data/Ancestral_Trait_Estimates.RData")

#######################################################################

# combine the extant and estimated ancestral traits together

# isolate traits (exclude Genus, Genus_species)
LSR.anc <- allLSR[,1:19]
# combine observed and ancestral traits into a single dataframe
LSR.anc <- dplyr::bind_rows(LSR.anc, ancestors)

#######################################################################


allLSR$Genus <- sapply(allLSR$Genus, function(x) strsplit(x,"_")[[1]][1])

tymp <- dplyr::filter(allLSR, Genus == "Tympanocryptis")
dipo <- dplyr::filter(allLSR, Genus == "Diporiphora")
pogo <- dplyr::filter(allLSR, Genus %in% c("Pogona","Rankinia"))
amph <- dplyr::filter(allLSR, Genus %in% c("Tropicagama","Amphibolurus",
                                           "Gowidon", "Lophognathus",
                                           "Chlamydosaurus", "Intellagama",
                                           "Chelosania"))
molo <- dplyr::filter(allLSR, Genus_species == "Moloch_horridus")
hyps <- dplyr::filter(allLSR, Genus == "Hypsilurus")
loph <- dplyr::filter(allLSR, Genus %in% c("Lophosaurus","Physignathus"))

centroids <- data.frame(t(data.frame(
  tymp = apply(tymp[,1:19], 2, mean),
  dipo = apply(dipo[,1:19], 2, mean),
  pogo = apply(pogo[,1:19], 2, mean),
  amph = apply(amph[,1:19], 2, mean),
  molo = apply(molo[,1:19], 2, mean),
  hyps = apply(hyps[,1:19], 2, mean),
  loph = apply(loph[,1:19], 2, mean)
)))

euclidean <- function(a, b) sqrt(sum((a - b)^2))
euclidean.centroid <- function(a) sqrt(sum((a - curr.taxon)^2))

centroid.distances <- NULL
for(k in 1:nrow(LSR.anc)){
  curr.taxon <- LSR.anc[k,1:19]
  curr.dist <- apply(centroids, 1, euclidean.centroid)
  centroid.distances <- rbind(centroid.distances, curr.dist)
}
centroid.distances <- data.frame(centroid.distances)
centroid.distances$Taxon <- rownames(LSR.anc)

testo <- tidyr::pivot_longer(centroid.distances,
                             !Taxon, names_to = "centroid")
testo <- dplyr::filter(testo, stringr::str_detect(Taxon, "^n"))

ggplot() +
  geom_point(data=testo, aes(x=Taxon, y=value, color=centroid)) +
  coord_flip()


##################################

acentroid = apply(allLSR[,1:19], 2, mean)
euclidean(allLSR[1,1:19], acentroid)

centroid.distances <- NULL
for(k in 1:nrow(allLSR)){
  curr.taxon <- allLSR[k,1:19]
  #curr.dist <- apply(acentroid, 1, euclidean.centroid)
  curr.dist <- euclidean(acentroid, curr.taxon)
  centroid.distances <- rbind(centroid.distances, curr.dist)
}
centroid.distances <- data.frame(centroid.distances)
centroid.distances$Taxon <- rownames(allLSR)

testo <- tidyr::pivot_longer(centroid.distances,
                             !Taxon, names_to = "centroid")
testo <- dplyr::filter(testo, stringr::str_detect(Taxon, "^n"))

ggplot() +
  geom_point(data=testo, aes(x=Taxon, y=value, color=centroid)) +
  coord_flip()

distances <- setNames(centroid.distances$centroid.distances, centroid.distances$Taxon)
phytools::plotTree.lollipop(agam.tree, distances, args.plotTree=list(fsize=0.5,ftype="i"))

#####################################

midpoint <- function(a){(max(a)-min(a))/2}
mids <- apply(allLSR[,1:19],2,midpoint)

mid.distances <- NULL
for(k in 1:nrow(allLSR)){
  curr.taxon <- allLSR[k,1:19]
  #curr.dist <- apply(acentroid, 1, euclidean.centroid)
  curr.dist <- euclidean(mids, curr.taxon)
  mid.distances <- rbind(mid.distances, curr.dist)
}
mid.distances <- data.frame(mid.distances)
mid.distances$Taxon <- rownames(allLSR)

mdistances <- setNames(mid.distances$mid.distances, mid.distances$Taxon)
mdistances <- mdistances - min(mdistances)
phytools::plotTree.lollipop(agam.tree, mdistances, args.plotTree=list(fsize=0.5,ftype="i"))

#####################################

LSR.pca <- princomp(allLSR[,1:19])
mpca <- apply(LSR.pca$scores[,1:19], 2, mean)
pca.scores <- data.frame(LSR.pca$scores)

pca.distances <- NULL
for(k in 1:nrow(pca.scores)){
  curr.taxon <- pca.scores[k,1:3]
  #curr.dist <- apply(acentroid, 1, euclidean.centroid)
  curr.dist <- euclidean(mpca[1:3], curr.taxon)
  pca.distances <- rbind(pca.distances, curr.dist)
}
pca.distances <- data.frame(pca.distances)
pca.distances$Taxon <- rownames(allLSR)
pca.distances[order(pca.distances$pca.distances),]

mdistances <- setNames(mid.distances$mid.distances, mid.distances$Taxon)
mdistances <- mdistances - min(mdistances)
phytools::plotTree.lollipop(agam.tree, mdistances, args.plotTree=list(fsize=0.5,ftype="i"))


small.pca <- dplyr::filter(pca.scores, rownames == "Ctenophorus_slater")

pca.scores["center_all",] <- apply(pca.scores, 2, mean)
small.pca <- pca.scores %>%
  tibble::rownames_to_column(var = "name") %>%
  dplyr::filter(name %in% c("Moloch_horridus",
                            "Lophosaurus_boydii",
                            "Chlamydosaurus_kingii",
                            "Cryptagama_aurita",
                            "Tympanocryptis_cephalus",
                            "Diporiphora_superba",
                            "center_all"))
ggplot() +
  geom_point(data=small.pca, aes(x=Comp.1,y=Comp.2, color=name))


####################################

unlist(allLSR[1,1:19])

acentroid = apply(allLSR[,1:19], 2, mean)
mahalanobis(x=unlist(allLSR[1,1:19]),
            center=unlist(ancestors[1,]),
            cov=cov(allLSR[,1:19]))

mahalanobis(x=unlist(acentroid[1:6]),
            center=unlist(allLSR[1,1:6]),
            cov=cov(allLSR[,1:6]))

mahalanobis(x=unlist(allLSR[1,1:16]),
            center=unlist(acentroid[1:16]),
            cov=cov(allLSR[,1:16]))
testo <- apply(allLSR[,1:16], 1, function(x) mahalanobis(x=unlist(x),
                                         center=unlist(acentroid[1:16]),
                                         cov=solve(cov(allLSR[,1:16]))))
testo[order(testo)]

if(distance=="mahalanobis"){for(k in 1:nrow(ndel)){rdist <- append(rdist, mahalanobis(x=unlist(trait.df[Ntip(phy)+1,]), 
                                                                                      center=unlist(trait.df[which(rownames(trait.df)==ndel$name.child[[k]]),]),
                                                                                      cov=cov(trait.df)))}}

allLSR2 <- allLSR[,1:19]
for(j in 1:length(acentroid)){
  curr.trait <- names(acentroid)[[j]]
  allLSR2[curr.trait] <- abs(allLSR2[curr.trait] - acentroid[[j]])
}
allLSR2$distance <- rowSums(allLSR2)

distances <- setNames(allLSR2$distance, rownames(allLSR2))
phytools::plotTree.lollipop(agam.tree, distances, args.plotTree=list(fsize=0.5,ftype="i"))
