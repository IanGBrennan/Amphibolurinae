# following methodology: https://github.com/sabifo4/mammals_dating/blob/ed05638c2f8541ae5f8e1d5a50316b3f2165ccd2/01_SeqBayes_S1/03_Fit_ST_to_posteriors/00_Fit_skewT.R#L18

library(sn)

# Function to process MCMCTree mcmc file into long format for plotting
process_mcmc <- function(file.path, rescale=c(1,100), chain.name){
  chain <- read.csv(file.path, sep="\t")
  chain <- chain[-1,]
  chain <- dplyr::select(chain, -starts_with("sigma2"))
  chain <- dplyr::select(chain, -starts_with("mu"))
  chain <- dplyr::select(chain, -starts_with("lnL"))
  chain <- dplyr::select(chain, -starts_with("X"))
  
  chain <- reshape(chain, idvar=c("Gen"), varying=2:ncol(chain), direction="long", sep="_")
  rownames(chain) <- NULL
  colnames(chain) <- c("generation", "node", "time")
  chain$time <- chain$time*rescale
  chain$chain <- chain.name
  
  return(chain)
}
# e.g. process_mcmc(file.path="...", rescale=100, chain.name="posterior")

# process the MCMCTree combined chains to get the posteriors
posterior <- process_mcmc(file.path = "/Users/ianbrennan/Desktop/AusARG_MCMCTree/ILN/Agamidae/Agamidae_Genera/outgroup_Phrynocephalus/mcmctree_C5_GC1_R1/mcmc_C5_ALL_ILN_HKY_R1.txt",
                          rescale = 1, chain.name = "posterior")

#posterior$time <- posterior$time / 100


nodes <- unique(posterior$node)
st.fits <- NULL
for (i in 1:length(nodes)){
  print(paste("processing node", nodes[[i]]))
  curr.node <- dplyr::filter(posterior, node == nodes[[i]])
  curr.st <- st.mple(y = curr.node$time, opt.method = "BFGS")
  curr.dp <- data.frame(t(curr.st$dp))
  curr.dp$node <- nodes[[i]]
  curr.dp$norm.mean <- mean(curr.node$time*100)
  curr.dp$norm.sd <- sd(curr.node$time*100)
  st.fits <- rbind(st.fits, curr.dp)
}
saveRDS(st.fits, file="/Users/ianbrennan/Desktop/AusARG_MCMCTree/ILN/Agamidae/Agamidae_Genera/outgroup_Phrynocephalus/mcmctree_C5_GC1_R1/mcmctree_genera_C5_GC1_SkewTPriors.rds")
write.csv(st.fits, file="/Users/ianbrennan/Desktop/AusARG_MCMCTree/ILN/Agamidae/Agamidae_Genera/outgroup_Phrynocephalus/mcmctree_C5_GC1_R1/mcmctree_genera_C5_GC1_SkewTPriors.csv", quote=F, row.names=F)

par(mfrow=c(4,5))
#par(mar = c(4, 2, 0.8, 0.2))

for(k in 1:length(nodes)){
  curr.post <- dplyr::filter(posterior, node==nodes[[k]])
  min.post <- min(curr.post$time) - (min(curr.post$time)*0.1)
  max.post <- max(curr.post$time) + (max(curr.post$time)*0.1)
  curr.st <- st.fits[which(st.fits$node == nodes[[k]]),]
  plot(density(curr.post$time, adj=1),
       xlim = c(min.post,max.post), col="#ffc247", main=nodes[[k]],
       xlab = paste("mean =", round(mean(curr.post$time),2)), ylab = NA, lwd = 2)
  norm.val <- rnorm(100000, mean=curr.st$norm.mean/100, sd=curr.st$norm.sd/100)
  lines(density(norm.val), col = "lightblue", lwd = 2)
  curve(dst(x, xi = curr.st$xi, omega = curr.st$omega,
            alpha = curr.st$alpha, nu = curr.st$nu),
        from = min.post, to = max.post, add = T, col = "#db3125", lwd = 2)
  points(x=mean(curr.post$time), y=0, col = "lightblue", pch=19, cex=1)
}

intree <- read.tree("/Users/ianbrennan/Desktop/AusARG_MCMCTree/ILN/Agamidae/Agamidae_Genera/outgroup_Phrynocephalus/mcmctree_C3_GC1_R1/C3_R1.tre")
intree <- intree[[2]]
plot(intree, cex=0.4); nodelabels(fr="c", cex=0.6)


