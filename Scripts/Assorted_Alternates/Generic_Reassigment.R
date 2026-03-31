# set the working directory to the GitHub repo
setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae")

############################################################################

library(ggplot2)
library(ggsankey)

############################################################################

taxo <- read.csv("Data/Generic_Reassignment.csv", h=T)

# Old Genera --> Old Names --> New Genera
assn <- taxo %>% make_long(Genus_Greer, Genus_2025)


ggplot(assn, aes(x = x, 
                  next_x = next_x,
                  node = node, 
                  next_node = next_node, 
                  #color = node,
                  fill = node,
                  label = node)) +
  geom_sankey(show.legend = F, linetype = "solid", flow.alpha=0.75, node.color=1) + 
  #geom_sankey_label(show.legend = F) +
  geom_sankey_label(size = 3.5, color = 1, fill = "white", show.legend=F) +
  theme_sankey(base_size = 7) +
  scale_fill_viridis_d(alpha = 1)
