# set the working directory to the GitHub repo
setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

############################################################################

# Load necessary packages
library(dplyr)
library(ggplot2)
library(phytools)
library(RColorBrewer)
library(tidyr)
library(googlesheets4); gs4_deauth()

############################################################################

# Read in the raw morphological data
rdata <- read_sheet("https://docs.google.com/spreadsheets/d/1NbvnODOT1HEv_9ojYJr5h0XUAW7I1HNDPSI3SQT3z-8/edit#gid=0")

# Select the variables of interest
xdata <- rdata %>%
  dplyr::filter(Subfamily == "Amphibolurinae") %>%
  dplyr::select(Genus_species, Snout_Vent, Snout_Axilla, Interlimb,
                Body_Width, Pelvic_Width, Pelvic_Height, 
                Head_Width, Head_Length, Snout_Eye, Eye_Diameter, Head_Depth, 
                Tail_Width, Tail_Length, #Tail_Loss,
                Upper_Arm, Lower_Arm, Hand, 
                Upper_Leg, Lower_Leg, Foot); nrow(xdata)

# exclude incomplete samples
xdata <- xdata[complete.cases(xdata),]; nrow(xdata)

############################################################################

# Some measurements are can be split into multiple variables

# determine neck length
xdata <- dplyr::mutate(xdata, Neck = Snout_Axilla - Head_Length)

# split head_length into tripartite model: snout, eye, posterior skull
xdata <- dplyr::mutate(xdata, Pos_Skull = Head_Length - (Snout_Eye + Eye_Diameter))

# We might want to drop the variables that you split (only if you want)
samp.data <- dplyr::select(xdata, -Snout_Axilla, -Head_Length, -Snout_Vent)
write.csv(samp.data, row.names=F,
          file="Data/Amphibolurinae_Morphology_AllSamples.csv")

############################################################################

# Create a tibble to get the species means  
sp.means <- xdata %>%
  group_by(Genus_species) %>%
  summarise_at(vars(Snout_Vent, Interlimb, Snout_Axilla, 
                    Body_Width, Pelvic_Width, Pelvic_Height, 
                    Head_Width, Snout_Eye, Eye_Diameter, Head_Depth, Head_Length,
                    Tail_Width, Tail_Length, 
                    Upper_Arm, Lower_Arm, Hand,
                    Upper_Leg, Lower_Leg, Foot,
                    Neck, Pos_Skull), mean)
# make the species names machine readable
sp.means$Genus_species <- sapply(sp.means$Genus_species, function(x) gsub(" ", "_", x))

# Save the species mean file
write.csv(sp.means, row.names=FALSE,
          file="Data/Amphibolurinae_Morphology_spMEANS.csv")

############################################################################

# Drop some variables before we calculate our log-shape ratios (geometric mean of size)
sp.means <- dplyr::select(sp.means, -Snout_Vent, -Snout_Axilla, -Head_Length)

# make sure your variables are all numeric
all.raw.traits <- sp.means

# Compute the geometric mean for obtaining size
allsize <- apply(sp.means[2:19], 1, prod)^(1/ncol(sp.means[2:19]))

# Compute the log shape ratios
allLSR <- sp.means[2:19]/allsize
allLSR$Size <- allsize
allLSR <- log(allLSR)

# alternatively you can do the log-shape ratios with MorphoInd::lsr(x, transf=T)
# allLSR <- MorphoInd::lsr(sp.means[2:19, transf=T])
# you'll need to reassign the allLSR$lsr object and add allLSR$size, but the
# results should be identical

# fix the names of the LSR traits
colnames(allLSR)[1:19] <- sapply(colnames(allLSR)[1:19], function(x) gsub("_", "", x))
# make a custom DF for the lsr traits
all.LSR.traits <- allLSR;
rownames(all.LSR.traits) <- sp.means$Genus_species

# Add the taxon, group, and family, names back onto the data frame
rownames(allLSR) <- sp.means$Genus_species;
allLSR$Genus <- sapply(rownames(allLSR), function(x) strsplit(x, " ")[[1]][1])
allLSR$Genus_species <- sp.means$Genus_species
#allLSR$Genus_species <- sapply(allLSR$Genus_species, function(x) gsub(" ", "_", x))

# Save the log shape ratio file
write.csv(allLSR, row.names=FALSE,
          file="Data/Amphibolurinae_AllLSR.csv")

############################################################################

# Trim our preferred tree down to match our data
atree <- read.tree("Trees/Amphibolurinae_ALL.tre")
agam.tree <- drop.tip(atree, setdiff(atree$tip.label, allLSR$Genus_species))
# Name the nodes
agam.tree$node.label <- paste("n", (Ntip(agam.tree)+1):(Ntip(agam.tree) + Nnode(agam.tree)), sep="")
# Save the tree
write.tree(agam.tree, file="Trees/Amphibolurinae_forR.tre")

############################################################################

# Trim our data down to match our tree
allLSR <- allLSR[which(allLSR$Genus_species %in% agam.tree$tip.label),]
all.LSR.traits <- all.LSR.traits[which(rownames(all.LSR.traits) %in% agam.tree$tip.label),]

############################################################################

# make individual BayesTraits input files for each trait
all.LSR.traits <- filter(all.LSR.traits, rownames(all.LSR.traits) %in% agam.tree$tip.label)
for (k in 1:ncol(all.LSR.traits)){
  trait <- select(all.LSR.traits, colnames(all.LSR.traits)[[k]])
  tname <- colnames(trait)
  filename <- paste0("~/Desktop/", tname, ".txt")
  write.table(trait, file=filename, sep=" ", quote=F, col.names=F)
}


############################################################################

# visualize the PCA data colored by genus
LSR.pca <- prcomp(allLSR[,1:19])
LSR.x <- data.frame(LSR.pca$x); rownames(LSR.x) <- rownames(allLSR)
LSR.x$genus <- sapply(rownames(LSR.x), function(y) strsplit(y, "_")[[1]][1])

ggplot() +
  geom_point(data = transform(LSR.x, genus = NULL), 
             colour = "grey85", aes(x=PC1, y=PC2), size=2) +
  geom_point(data=LSR.x, aes(x=PC1, y=PC2, fill=genus), size=3, pch=21) +
  theme_classic() + theme(legend.position = "none") + facet_wrap(~genus)

############################################################################

# Make an RData object we can easily load
save(agam.tree, all.raw.traits, all.LSR.traits, allLSR, LSR.pca, file="Data/Amphibolurinae_Data.RData")

############################################################################


