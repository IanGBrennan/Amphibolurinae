library(ggplot2)
#library(ggfortify)
#library(patchwork)

setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

source("Scripts/trait.at.time.R")
source("Scripts/innovate_elaborate.R")

load("Data/Amphibolurinae_Data.RData")
load("Data/Ancestral_Trait_Estimates.RData")


############################################################################

# isolate traits (exclude Genus, Genus_species)
LSR.anc <- allLSR[,1:19]
# combine observed and ancestral traits into a single dataframe
LSR.anc <- dplyr::bind_rows(LSR.anc, ancestors)


############################################################################

# Extrapolate trait values along branches to get traits at specified timeslices
ttobj <- trait.at.time.multi(0.5, agam.tree, LSR.anc, plot=F)

############################################################################

# Fit PCA/Hypervolume at timeslices and then extract the volume (volume-through-time)
vtt <- extract.volume(obj = ttobj, PCs = 3, parallel = 8)
save(vtt, file="Data/Amphibolurinae_VolumeThroughTime.RData")

# assess functional evenness across timeslices
evenness <- ttobj %>%
  group_by(time) %>%
  group_map(~func.eve(., method="euclidean"))
even.df <- data.frame(time = unique(ttobj$time), even = unlist(evenness), time.rev = unique(ttobj$time) - max(ttobj$time))
plot(even.df$even ~ even.df$time.rev, type="l")

# assess functional divergence across timeslices
divergence <- ttobj %>%
  group_by(time) %>%
  group_map(~func.div(., method="euclidean"))
div.df <- data.frame(time = unique(ttobj$time), even = unlist(divergence), time.rev = unique(ttobj$time) - max(ttobj$time))
plot(div.df$even ~ div.df$time.rev, type="l")


############################################################################

library(hypervolume)
library(hyperoverlap)

# compare the hypervolumes of different ecologies
eco.df <- read.csv("Data/Amphibolurinae_Ecology.csv", h=T)
rownames(eco.df) <- eco.df$Genus_species
eco.df <- dplyr::select(eco.df, Genus_species, Ecology)
pca.df <- data.frame(prcomp(allLSR[,1:19])$x)
pca.df$Genus_species <- rownames(pca.df)
volume.df <- dplyr::full_join(pca.df, eco.df)

arbo.vol <- hypervolume(dplyr::filter(volume.df, Ecology=="arboreal")[,1:3])
semi.vol <- hypervolume(dplyr::filter(volume.df, Ecology=="semiarboreal")[,1:3])
terr.vol <- hypervolume(dplyr::filter(volume.df, Ecology=="terrestrial")[,1:3])
saxi.vol <- hypervolume(dplyr::filter(volume.df, Ecology=="saxicolous")[,1:3])

arbo.vol@Volume; semi.vol@Volume; terr.vol@Volume; saxi.vol@Volume

arb.sem.perm <- hypervolume::hypervolume_permute("~/Desktop/Arboreal_Semiarboreal", 
                                                 arbo.vol, semi.vol, n=200, cores=8)

hypervolume::hypervolume_overlap_test(tym.hyp, dip.hyp)


# Using PCs

avt <- dplyr::filter(volume.df, Ecology %in% c("arboreal", "terrestrial"))
res.avt <- hyperoverlap_detect(avt[,1:3], avt$Ecology)
res.avt.n <- res.avt <- hyperoverlap_detect(avt[,1:19], avt$Ecology)
# Result: non-overlap

avs <- dplyr::filter(volume.df, Ecology %in% c("arboreal", "semiarboreal"))
res.avs <- hyperoverlap_detect(avs[,1:3], avs$Ecology)
res.avs.n <- hyperoverlap_detect(avs[,1:19], avs$Ecology)
hyperoverlap_lda(res.avs.n)
# Result: overlap

avr <- dplyr::filter(volume.df, Ecology %in% c("arboreal", "saxicolous"))
res.avr <- hyperoverlap_detect(avr[,1:3], avr$Ecology)
res.avr.n <- hyperoverlap_detect(avr[,1:19], avr$Ecology)
# Result: non-overlap

tvs <- dplyr::filter(volume.df, Ecology %in% c("terrestrial", "semiarboreal"))
res.tvs <- hyperoverlap_detect(tvs[,1:3], tvs$Ecology)
res.tvs.n <- hyperoverlap_detect(tvs[,1:19], tvs$Ecology)
# Result: overlap

tvr <- dplyr::filter(volume.df, Ecology %in% c("terrestrial", "saxicolous"))
res.tvr <- hyperoverlap_detect(tvr[,1:3], tvr$Ecology)
res.tvr.n <- hyperoverlap_detect(tvr[,1:19], tvr$Ecology)
# Result: overlap (3PCs); non-overlap (19 PCs)

# Using Raw Traits

all.eco <- dplyr::full_join(allLSR, eco.df)
all.eco <- all.eco[,c(1:19,22)]

avt <- dplyr::filter(all.eco, Ecology %in% c("arboreal", "terrestrial"))
res.avt.n <- hyperoverlap_detect(avt[,1:19], avt$Ecology); res.avt.n@result
# Result: non-overlap

avs <- dplyr::filter(all.eco, Ecology %in% c("arboreal", "semiarboreal"))
res.avs.n <- hyperoverlap_detect(avs[,1:19], avs$Ecology); res.avs.n@result
# Result: non-overlap

avr <- dplyr::filter(all.eco, Ecology %in% c("arboreal", "saxicolous"))
res.avr.n <- hyperoverlap_detect(avr[,1:19], avr$Ecology); res.avr.n@result
# Result: non-overlap

tvs <- dplyr::filter(all.eco, Ecology %in% c("terrestrial", "semiarboreal"))
res.tvs.n <- hyperoverlap_detect(tvs[,1:19], tvs$Ecology); res.tvs.n@result
# Result: overlap

tvr <- dplyr::filter(all.eco, Ecology %in% c("terrestrial", "saxicolous"))
res.tvr.n <- hyperoverlap_detect(tvr[,1:19], tvr$Ecology); res.tvr.n@result
# Result: non-overlap

rvs <- dplyr::filter(all.eco, Ecology %in% c("saxicolous", "semiarboreal"))
res.rvs.n <- hyperoverlap_detect(rvs[,1:19], rvs$Ecology); res.rvs.n@result
# Result: non-overlap

############################################################################

tym <- extract.clade(agam.tree, node=191)
dip <- extract.clade(agam.tree, node=213)

tym.data <- filter(allLSR, Genus == "Tympanocryptis")
dip.data <- filter(allLSR, Genus == "Diporiphora")

tym.hyp <- hypervolume::hypervolume(prcomp(tym.data[,1:19])$x[,1:3])
dip.hyp <- hypervolume::hypervolume(prcomp(dip.data[,1:19])$x[,1:3])

tym.dip.perm <- hypervolume::hypervolume_permute("~/Desktop/Tympanocryptis_Diporiphora", 
                                                 tym.hyp, dip.hyp, n=200, cores=8)

hypervolume::hypervolume_overlap_test(tym.hyp, dip.hyp)

############################################################################




















il.30 <- dplyr::filter(ttobj, time==30)
dispRity::func.eve(as.matrix(il.20[,2:20]))


testy %>%
  dplyr::group_by(time) %>%
  dispRity::func.eve(as.matrix(.[,2:4]), method="euclidean")

testi <- testy %>% dplyr::group_by(time)


mtcars %>%
  group_by(cyl) %>%
  group_map(~ head(.x, 2L))

testy %>%
  group_by(time) %>%
  group_map(~func.eve(., method="euclidean"))
  #group_map(mean(.$Interlimb))

func.eve(as.matrix(il.30))


# extract trait disparity at each point in time
emp.dis.trait <- lapply(2:length(testo), function(x) extract.variance(testo[,c(1,x)], plot = T, metric = "disparity"))
names(emp.dis.trait) <- names(testo)[2:20]
plot(emp.dis.trait$TailLength$measure ~ emp.dis.trait$TailLength$time, type="l")
# extract trait disparity at each point in time
emp.eve.trait <- lapply(2:length(testo), function(x) extract.variance(testo[,c(1,x)], plot = T, metric = "functional evenness"))
names(emp.eve.trait) <- names(testo)[2:20]
plot(emp.eve.trait$TailLength$measure ~ emp.eve.trait$TailLength$time, type="l")
# extract trait variance at each point in time
emp.var.trait <- lapply(2:length(testo), function(x) extract.variance(testo[,c(1,x)], plot = T, metric = "variance"))
names(emp.var.trait) <- names(testo)[2:20]
plot(emp.var.trait$TailLength$measure ~ emp.var.trait$TailLength$time, type="l")
plot(emp.var.trait$HeadWidth$measure ~ emp.var.trait$HeadWidth$time, type="l")

emp.eve <- extract.

extract.variance(testo[,c(1,11)], plot="sideXside", metric="disparity")


agg.pca <- function(x) prcomp(x)$x

PCs <- 3

testy <- prcomp(testo[1:10,])

aggregate(testo[,2], list(testo$time), agg.pca)

curr <- filter(testo, time==31.635)
curr.pca <- prcomp(curr[,2:20])
curr.hyp <- hypervolume(curr.pca$x[,1:3])

curr.minus <- filter(testo, time==31.000)
curr.pca.minus <- prcomp(curr.minus[,2:20])
curr.hyp.minus <- hypervolume(curr.pca.minus$x[,1:3])


utime <- unique(testo$time)
volume <- NULL
for(k in 1:length(utime)){
  curr.traits <- dplyr::filter(testo, time==utime[[k]])
  if(!nrow(curr.traits) >= 3){next}
  curr.pca <- prcomp(curr.traits[,2:20])
  curr.hyp <- hypervolume(curr.pca$x[,1:3])
  volume <- rbind(volume, data.frame(time=utime[[k]], volume=get_volume(curr.hyp)))
}

# split the trait.time.object into a list of time dataframes
tto.list <- NULL
utime <- unique(testo$time)
for(j in 1:length(utime)){
  curr.traits <- dplyr::filter(testo, time==utime[[j]])
  if(!nrow(curr.traits) >= 3){next}
  tto.list[[j]] <- curr.traits
  names(tto.list)[[j]] <- paste0("time_",utime[[j]])
}
tto.list <- tto.list[which(!names(tto.list)=="NULL")]


hypervolumes <- parallel::mclapply(tto.list, reduce.hypervolume, mc.cores = 8)

volumes <- NULL
for(k in 1:length(tto.list)){
  message(paste("running hypervolume on", names(tto.list)[[k]]))
  volumes[[k]] <- reduce.hypervolume(tto.list[[k]])
}

volumes2 <- lapply(volumes, get_volume)
vdf <- data.frame(time = sapply(names(tto.list)[1:20], function(x) strex::str_after_first(x, "_")), volume = unlist(volumes2))
plot(vdf$volume ~ vdf$time)










pca <- prcomp(LSR.anc)

pca.rot <- data.frame(pca$rotation)
# create a variable with the rownames
pca.rot$trait <- rownames(pca.rot)
# reshape the data from wide to long format
pca.rot <- reshape2::melt(pca.rot, id.vars="trait", variable.name="PC")
# plot the resulting data
rot.plot <- ggplot(pca.rot,aes(x=trait,y=value)) +
  geom_col() + theme_bw() +
  facet_wrap(~PC) + coord_flip() + scale_x_discrete(limits=rev)



  
  
  
  
  
  
  


