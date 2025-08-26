# set the working directory to the GitHub repo
setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

############################################################################

# Load necessary packages

library(phytools)
library(ggtree) # BiocManager::install("ggtree")

############################################################################

# read in the genus-level tree for amphibolurines
backbone <- read.tree('Trees/Amphibolurinae_Genera_hASTRAL.tre')

# plot the base tree
p + theme_tree2() + geom_nodepoint() + geom_tiplab()

# read in the concordance factor values
gcf <- read.csv('Data/Amphibolurinae_Genera_gCF.csv')
gcf$node <- gcf$ID

# select the columns of interest
gcf <- dplyr::select(gcf, node, gCF, gDF1, gDF2, gDFP)

# set colors for the pie charts and set names to the colors
gcf.colors <- c("#ffca38", "#4ec8de", "#ff8fb2", "#dbdbdb")
names(gcf.colors) <- c("gCF", "gDF1", "gDF2", "gDFP")

# establish the pie charts and add them to the tree
pies <- nodepie(gcf, cols=2:5, outline.color="black", outline.size=0.25)
pies <- lapply(pies, function(x) x+scale_fill_manual(values = gcf.colors))
inset(p+geom_tiplab(), pies, width=0.1, height=0.1, hjust=0)

############################################################################
