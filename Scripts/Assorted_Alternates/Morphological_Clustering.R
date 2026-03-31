# Defining Morphological Clusters

# determine the pairwise euclidean distance among species
mdist <- vegan::vegdist(allLSR[,1:19], method="euclidean", diag=T, upper=T)

# cluster these distances into a dendrogram
top.clust <- hclust(dist(mdist), method="single")
plot(top.clust)

# make a more readable tree by converting the cluster object to use with ggplot
dendr <- ggdendro::dendro_data(top.clust, type="rectangle") 
#your own labels (now rownames) are supplied in geom_text() and label=label
ggplot() + 
  geom_segment(data=ggdendro::segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=ggdendro::label(dendr), aes(x=x, y=y, label=label, hjust=0), size=3) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())


# now we'll determine the clusters
library(cluster)
library(NbClust)
library(factoextra)

# chose just the traits we're interested in 
scaled <- scale(allLSR[,1:19])

# Quick assessment of the optimal number of morphological clusters
fviz_nbclust(scaled, kmeans, method="wss")
fviz_nbclust(scaled, kmeans, method="silhouette")
fviz_nbclust(scaled, pam, method="wss")
fviz_nbclust(scaled, pam, method="silhouette")

# we can do this more extensively
nb <- NbClust(scaled[,3:19], distance = "euclidean", min.nc = 2,
              max.nc = 10, method = "complete", index ="all")

# or we can use PAM with a predefined number of clusters
pam.res <- pam(allLSR[,1:19], metric="euclidean", 3)
output <- fviz_cluster(pam.res, stand = F, geom = "point",
                       ellipse.type = "norm", show.clust.cent=T)
# plot the clusters
output+theme_classic()

# save them to the original dataframe and export that
allLSR[,"k12"] <- pam.res$clustering
write.csv(allLSR, file="~/Desktop/Amphibolurinae_Morph_Cluster.csv")

# extra info here: 
# https://stats.stackexchange.com/questions/2597/what-stop-criteria-for-agglomerative-hierarchical-clustering-are-used-in-practice
# https://stats.stackexchange.com/questions/2717/clustering-with-a-distance-matrix