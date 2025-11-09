library(dplyr)
library(ggplot2)
library(patchwork)

setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

#######################################################################

ranges <- read.csv("Data/GARD_IBRA.csv"); rownames(ranges) <- NULL
niches <- read.csv("Data/Habitat_Use.csv")

#######################################################################

# remove unclassified species
ranges <- ranges %>%
  tibble::column_to_rownames(var="taxon") %>%
  dplyr::select(corrected_subregions) %>%
  tidyr::drop_na()

# select just the niche breadth
niches <- niches %>%
  tibble::column_to_rownames(var="taxon") %>%
  dplyr::select(breadth)

# combine the data frames together
range.niche <- merge(ranges, niches, by="row.names")
# give it a column for 'Genus'
range.niche$genus <- strex::str_before_first(range.niche$Row.names, "_")

#######################################################################

# VISUALIZATIONS

# we can plot niche breadth as a function of number of IBRA subregions
ggplot() +
  geom_point(data=range.niche, aes(x=corrected_subregions,y=breadth,fill=genus),shape=21,size=4)
#  geom_point(data=range.niche, aes(x=breadth,y=genus,fill=genus),shape=21,size=4)
# but this doesn't really tell us that much

# plot the number of subregions per species by genus
subregions.plot <- ggplot() +
  geom_jitter(data=range.niche, aes(x=corrected_subregions,y=genus,fill=genus),shape=21,size=3) +
  scale_y_discrete(limits=rev) + theme_bw()

#######################################################################

# Transform the data so that we can see how many
# species in each genus fit in each breadth category
newie <- range.niche %>%
  dplyr::count(genus,breadth)

# get a table of the number of species in each genus
spp.table <- table(range.niche$genus)

# get the percent of species in each genus in each category
newie$percent <- sapply(1:nrow(newie), function(x) (newie$n[x]/spp.table[which(names(spp.table)==newie$genus[x])])*100)

# plot the breadth by genus
breadth.plot <- ggplot() +
  #geom_point(data=range.niche, aes(x=corrected_subregions,y=breadth,fill=genus),shape=21,size=4)
  geom_label(data=newie, aes(x=breadth,y=genus,label=n),nudge_x=0.25,size=3) +
  geom_point(data=newie, aes(x=breadth,y=genus,fill=genus,size=percent),shape=21) +
  scale_y_discrete(limits=rev) + theme_bw()

# combine the plots together
subregions.plot + breadth.plot + plot_layout(guides = "collect") & theme(legend.position="bottom")

#######################################################################

# read in the multi-OU regimes
regimes <- read.csv("Data/Amphibolurinae_Ecology.csv")

regimes <- regimes %>%
  tibble::column_to_rownames(var="Genus_species") %>%
  dplyr::select(opt5)


rownames(range.niche) <- range.niche$Row.names; range.niche <- dplyr::select(range.niche, !Row.names)
testo <- merge(range.niche, regimes, by="row.names")

# plot the number of subregions per species by genus
regime.range <- ggplot() +
  geom_jitter(data=testo, aes(x=corrected_subregions,y=opt5,fill=opt5),shape=21,size=3) +
  scale_y_discrete(limits=rev) + theme_bw()

newie2 <- testo %>%
  dplyr::count(opt5,breadth)

regime.table <- table(testo$opt5)

# get the percent of species in each genus in each category
newie2$percent <- sapply(1:nrow(newie2), function(x) (newie2$n[x]/regime.table[which(names(regime.table)==newie2$opt5[x])])*100)

regime.breadth <- ggplot() +
  #geom_point(data=range.niche, aes(x=corrected_subregions,y=breadth,fill=genus),shape=21,size=4)
  geom_label(data=newie2, aes(x=breadth,y=opt5,label=n),nudge_x=0.25,size=3) +
  geom_point(data=newie2, aes(x=breadth,y=opt5,fill=opt5,size=percent),shape=21) +
  scale_y_discrete(limits=rev) + theme_bw()

regime.range + regime.breadth + plot_layout(guides = "collect") & theme(legend.position="bottom")

            