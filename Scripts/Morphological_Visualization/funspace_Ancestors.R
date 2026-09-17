ctenophorus1 <- c("Ctenophorus_maculosus",
                  "Ctenophorus_clayi",
                  "Ctenophorus_parviceps",
                  "Ctenophorus_butlerorum",
                  "Ctenophorus_adelaidensis",
                  "Ctenophorus_chapmani",
                  "n133","n134")
ctenophorus2 <- c(agam.tree$tip.label[20:50], 
                  paste("n",140:171,sep=""))
generalist1 <- c(agam.tree$tip.label[52:61],
                 agam.tree$tip.label[94:119],
                 paste("n",172:178,sep=""),
                 "n131","n132","n179","n180","n190")
tree1 <- c(agam.tree$tip.label[1:6],
          agam.tree$tip.label[9:11],
          paste("n",120:126,sep=""),
          "n127","n128","n129","n130",
          "Chelosania brunnea")
tree.anc <- c("n120","n121","n126","n127","n128")
anc.pogo <- c("n131","n132",
              paste("n",172:181,sep=""))
anc.cry <- c(anc.pogo, "n133","n134","n139")

group.c1 <- ifelse(rownames(curr.pca$scores) %in% ctenophorus1, "Ctenophorus 1", "Else")
group.c2 <- ifelse(rownames(curr.pca$scores) %in% ctenophorus2, "Ctenophorus 2", "Else")
group.g1 <- ifelse(rownames(curr.pca$scores) %in% generalist1, "Generalist 1", "Else")
group.t1 <- ifelse(rownames(curr.pca$scores) %in% tree1, "Tree 1", "Else")
group.ta <- ifelse(rownames(curr.pca$scores) %in% tree.anc, "Tree Ancestor", "Else")
group.ap <- ifelse(rownames(curr.pca$scores) %in% anc.pogo, "Generalist Ancestor to Pogona", "Else")
group.ac <- ifelse(rownames(curr.pca$scores) %in% anc.cry, "Generalist Ancestor to Cryptagama", "Else")


fs.group <- funspace::funspace(x = curr.pca, PCs = c(1,4), n_divisions = 300, group.vec = group.ac, threshold=0.95)

plot(x = fs.group,
     type = "groups",
     quant.plot = T, quant = 0.95,
     globalContour = F,
     #     arrows = T, arrows.length = 3, arrows.label.cex = 0.5, 
     pnt = F,
     threshold = 0.95, 
     xlim = c(-2,2), ylim = c(-0.7,0.7),
     colors = brewer.pal(5, "YlOrRd"))
