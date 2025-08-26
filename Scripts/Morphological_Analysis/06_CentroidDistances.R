
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

