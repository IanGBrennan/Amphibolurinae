require(sn)

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

# Function to sample from uniform MCMCTree priors for long format plotting
process_prior <- function(file.path, rescale=c(1,100), samples=20000, shape=c("uniform","skewT")){
  priors <- read.csv(file.path)
  
  if(shape=="uniform"){
    priors[,c("min","max")] <- priors[,c("min","max")]*rescale
    
    prior.df <- NULL
    for(k in 1:nrow(priors)){
      curr.prior <- priors[k,]
      if(is.na(curr.prior$min)){next}
      
      psamp <- runif(samples, min=curr.prior$min, max=curr.prior$max)
      psamp.df <- data.frame(generation=0, node=curr.prior$node,
                             time=psamp, chain="prior")
      prior.df <- rbind(prior.df, psamp.df)
    }
  }
  if(shape=="skewT"){
    prior.df <- NULL
    for(k in 1:nrow(priors)){
      curr.prior <- priors[k,]
      if(is.na(curr.prior$xi)){next}
      
      psamp <- rst(n=samples, xi=curr.prior$xi, omega=curr.prior$omega,
                   alpha=curr.prior$alpha, nu=curr.prior$nu)
      psamp.df <- data.frame(generation=0, node=curr.prior$node,
                             time=psamp*rescale, chain="prior")
      prior.df <- rbind(prior.df, psamp.df)
    }
  }

  rownames(prior.df) <- NULL
  return(prior.df)
}
# e.g. process_prior(file.path="...", rescale=100, samples=20000)

# process the MCMCTree analysis run on priors with no data (effective prior)
eff.prior <- process_mcmc(file.path = "/Users/ianbrennan/Desktop/AusARG_MCMCTree/ILN/Agamidae/outgroup_Phrynocephalus/mcmctree_C4osb_GC1_PRIOR/mcmc_C4_GC1_ILN_HKY_R2.txt",
                      rescale = 100, chain.name = "eff.prior")

# process the MCMCTree combined chains to get the posteriors
posterior <- process_mcmc(file.path = "/Users/ianbrennan/Desktop/AusARG_MCMCTree/ILN/Agamidae/outgroup_Phrynocephalus/mcmctree_C4osb_GC1_R1-3/mcmc_C4_GC1_ILN_HKY_R123.txt",
                          rescale = 100, chain.name = "posterior")

# sample from the full priors
prior <- process_prior(file.path = "/Users/ianbrennan/Documents/GitHub/Amphibolurinae/Data/Calibrations_Full_2.csv",
                       rescale = 100, samples = 20000)

# combine the applied priors with the effective priors and posteriors
p.p <- rbind(prior, rbind(eff.prior, posterior))

# plot the data
ggplot(p.p) +
  geom_density(aes(x=time, color=chain), trim=F) +
  facet_wrap(~node, scales="free") + 
  theme_classic() + theme(legend.position = "none")
