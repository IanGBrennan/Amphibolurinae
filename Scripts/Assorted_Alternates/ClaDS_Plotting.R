library(ggplot2)
library(RColorBrewer)
require(fields)

# load("/Users/ianbrennan/Documents/GitHub/Egernia_Evolution/ClaDS/Egerniinae_ClaDS.Rdata")
# egout <- CladsOutput

# create a function to plot the Diversity Through Time of a ClaDS object from RPANDA-Julia
clads_DTT <- function(CladsOutput, log.sp=F, return.df=F){
  cdf <- data.frame(time=CladsOutput$time_points-max(CladsOutput$time_points),
                    no.species = CladsOutput$DTT_mean)
  if(log.sp==T){
    dtt <- ggplot(cdf, aes(x=time, y=no.species)) +
      scale_y_log10() +
      geom_line(color="#4dd1b7") + theme_bw() +
      labs(x="Time Before Present", y="Number of Species (log10 axis)")
  }else{
    dtt <- ggplot(cdf, aes(x=time, y=no.species)) +
      geom_line(color="#4dd1b7") + theme_bw() +
      labs(x="Time Before Present", y="Number of Species")
  }
  if(return.df==T){print(dtt); print(cdf)}
  if(return.df==F){print(dtt)}
}

# clads_DTT(egout, log.sp=T, return.df=F)

# create a function to plot the mean Rate Through Time of a ClaDS object from RPANDA-Julia
clads_RTT <- function(CladsOutput, log.rate=F, return.df=F, ci=T){
  # check if there are missing rates in the data
  missing.rates <- (length(CladsOutput$time_points) - length(CladsOutput$RTT_map))
  
  
  # make a ClaDS data frame of the time and mean a posteriori rates
  cdf <- data.frame(time=CladsOutput$time_points-max(CladsOutput$time_points),
                    rate = c(rep(CladsOutput$RTT_map[[1]],missing.rates),
                             CladsOutput$RTT_map))
  
  # if you want to plot 90% confidence intervals around the MAP rate
  if(ci==T){
    # make a loop to iterate across the chains
    for(j in 1:length(CladsOutput$rtt_chains)){
      curr.chain <- CladsOutput$rtt_chains[[j]]
      # loop across each iteration within the chain, and remove the first 50% of samples
      for(k in (round(length(curr.chain)/2,0)):length(curr.chain)){
        chain.name <- paste("chain",j,k, sep=".")
        curr.curr <- data.frame(c(rep(curr.chain[[k]][1], missing.rates),curr.chain[[k]])); 
        colnames(curr.curr) <- chain.name
        cdf <- cbind(cdf, curr.curr)
      }
    }
    # estimate the quantiles of from the posterior
    quants <- t(apply(cdf[,3:ncol(cdf)], 1, function(x) quantile(x, probs=c(0.95, 0.5, 0.05))))
    cdf <- cbind(cdf, quants)
    qmax <- max(cdf$`95%`); qmin <- min(cdf$`5%`)
    
    if(log.rate==T){
      rtt <- ggplot(cdf) +
        scale_y_log10() +
        #scale_y_continuous(trans = "log") +
        geom_ribbon(aes(x=time, ymin=`5%`, ymax=`95%`), fill="#e6b035", alpha=0.25) +
        geom_line(aes(x=time, y=rate), color="#e6b035") + theme_bw() +
        labs(x="Time Before Present", y="Mean rate (log10 axis)")    
    }else{
      rtt <- ggplot(cdf) +
        geom_ribbon(aes(x=time, ymin=`5%`, ymax=`95%`), fill="#e6b035", alpha=0.25) +
        geom_line(aes(x=time, y=rate), color="#e6b035") + 
        theme_bw() + labs(x="Time Before Present", y="Mean Rate")     
    }
  }
  
  if(ci==F){
    if(log.rate==T){
      rtt <- ggplot(cdf, aes(x=time, y=rate)) +
        scale_y_log10() +
        #scale_y_continuous(trans = "log") +
        geom_line(color="#e6b035") + theme_bw() +
        labs(x="Time Before Present", y="Mean rate (log10 axis)")    
    }else{
      rtt <- ggplot(cdf, aes(x=time, y=rate)) +
        geom_line(color="#e6b035") + theme_bw() +
        labs(x="Time Before Present", y="Mean Rate")    
    }
  }

  if(return.df==T){print(rtt); return(cdf)}
  if(return.df==F){print(rtt)}
}

# clads_RTT(egout, log.rate=F, return.df=F, ci=T)

# create a function to plot the branch rates of a ClaDS object from RPANDA-Julia
clads_TREE <- function (CladsOutput, rates2 = NULL, same.scale = T, main = NULL, 
                        lwd = 2, log = T, show.tip.label = F, bpal = "Spectral", 
                        tree.type = c("phylogram","fan"), ...) 
{
  #Colors = colorRampPalette(c("steelblue2", "paleturquoise3", 
  #                            "palegreen2", "yellow2", "salmon1", "darkorange", "red", 
  #                            "red4"))(100)
  Colors = rev(colorRampPalette(brewer.pal(9, bpal))(100))
  phylo = CladsOutput$tree
  rates = CladsOutput$lambdai_map
  
  if (is.null(rates2)) {
    if (log) 
      rates = log(rates)
    if (isTRUE(all.equal(rep(as.numeric(rates[1]), length(rates)), 
                         as.numeric(rates)))) {
      col = rep(1, length(rates))
      plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
           main = main, edge.width = lwd, type = tree.type, ...)
      if (log) {
        image.plot(z = c(exp(rates[1]), 2 * exp(rates[1])), 
                   col = Colors, horizontal = T, legend.only = T)
      }
      else {
        image.plot(z = c(rates[1], 2 * rates[1]), col = Colors, 
                   horizontal = T, legend.only = T)
      }
    }
    else {
      col = round((rates - min(rates))/diff(range(rates)) * 
                    99) + 1
      plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
           main = main, edge.width = lwd, type = tree.type, ...)
      if (log) {
        min = min(rates)
        max = max(rates)
        m10 = floor(min/log(10))
        M10 = ceiling(max/log(10))
        if ((M10 - m10) < 4) {
          ticks = c(1, 2, 5)
        }
        else {
          ticks = 1
        }
        ticks = as.vector(sapply(m10:M10, function(k) {
          return(ticks * 10^k)
        }))
        lt = length(ticks[ticks > exp(min) & ticks < 
                            exp(max)])
        if (lt < 4) {
          ticks = c(round(exp(min), max(0, -1 * m10 + 
                                          (lt < 2))), ticks, round(exp(max), max(0, 
                                                                                 -1 * M10 + 1 + (lt < 2))))
        }
        image.plot(z = c(min, max), col = Colors, horizontal = T, 
                   legend.only = T, axis.args = list(at = log(ticks), 
                                                     labels = ticks))
      }
      else {
        image.plot(z = as.matrix(rates), col = Colors, 
                   horizontal = T, legend.only = T)
      }
    }
  }
  else {
    if (log) {
      rates = log(rates)
      rates2 = log(rates2)
    }
    if (same.scale) {
      min = min(min(rates), min(rates2))
      max = max(max(rates), max(rates2))
      par(mfrow = c(1, 2))
      col = round(((rates - min)/(max - min)) * 99) + 1
      plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
           edge.width = lwd, type = tree.type, ...)
      col = round(((rates2 - min)/(max - min)) * 99) + 
        1
      plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
           edge.width = lwd, type = tree.type, ...)
      par(mfrow = c(1, 1))
      if (log) {
        m10 = floor(min/log(10))
        M10 = ceiling(max/log(10))
        if ((M10 - m10) < 4) {
          ticks = c(1, 2, 5)
        }
        else {
          ticks = 1
        }
        ticks = as.vector(sapply(m10:M10, function(k) {
          return(ticks * 10^k)
        }))
        lt = length(ticks[ticks > exp(min) & ticks < 
                            exp(max)])
        if (lt < 4) {
          ticks = c(round(exp(min), max(0, -1 * m10 + 
                                          (lt < 2))), ticks, round(exp(max), max(0, 
                                                                                 -1 * M10 + 1 + (lt < 2))))
        }
        image.plot(z = c(min, max), col = Colors, horizontal = T, 
                   legend.only = T, axis.args = list(at = log(ticks), 
                                                     labels = ticks))
      }
      else {
        image.plot(z = c(min, max), col = Colors, horizontal = T, 
                   legend.only = T)
      }
    }
    else {
      par(mfrow = c(1, 2))
      if (isTRUE(all.equal(rep(rates[1], length(rates)), 
                           rates))) {
        col = rep(1, length(rates))
        plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
             edge.width = lwd, type = tree.type, ...)
        if (log) {
          image.plot(z = c(exp(rates[1]), 2 * exp(rates[1])), 
                     col = Colors, horizontal = T, legend.only = T)
        }
        else {
          image.plot(z = c(rates[1], 2 * rates[1]), col = Colors, 
                     horizontal = T, legend.only = T)
        }
      }
      else {
        col = round(((rates - min(rates))/(max(rates) - 
                                             min(rates))) * 99) + 1
        plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
             edge.width = lwd, type = tree.type, ...)
        if (log) {
          min = min(rates)
          max = max(rates)
          m10 = floor(min/log(10))
          M10 = ceiling(max/log(10))
          if ((M10 - m10) < 4) {
            ticks = c(1, 2, 5)
          }
          else {
            ticks = 1
          }
          ticks = as.vector(sapply(m10:M10, function(k) {
            return(ticks * 10^k)
          }))
          lt = length(ticks[ticks > exp(min) & ticks < 
                              exp(max)])
          if (lt < 4) {
            ticks = c(round(exp(min), max(0, -1 * m10 + 
                                            (lt < 2))), ticks, round(exp(max), max(0, 
                                                                                   -1 * M10 + 1 + (lt < 2))))
          }
          image.plot(z = c(min, max), col = Colors, horizontal = T, 
                     legend.only = T, axis.args = list(at = log(ticks), 
                                                       labels = ticks))
        }
        else {
          image.plot(z = as.matrix(rates), col = Colors, 
                     horizontal = T, legend.only = T)
        }
      }
      if (isTRUE(all.equal(rep(rates2[1], length(rates2)), 
                           rates2))) {
        col = rep(1, length(rates2))
        plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
             edge.width = lwd, type = tree.type, ...)
        if (log) {
          image.plot(z = c(exp(rates2[1]), 2 * exp(rates2[1])), 
                     col = Colors, horizontal = T, legend.only = T)
        }
        else {
          image.plot(z = c(rates2[1], 2 * rates2[1]), 
                     col = Colors, horizontal = T, legend.only = T)
        }
      }
      else {
        col = round(((rates2 - min(rates2))/(max(rates2) - 
                                               min(rates2))) * 99) + 1
        plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
             edge.width = lwd, type = tree.type, ...)
        if (log) {
          min = min(rates2)
          max = max(rates2)
          m10 = floor(min/log(10))
          M10 = ceiling(max/log(10))
          if ((M10 - m10) < 4) {
            ticks = c(1, 2, 5)
          }
          else {
            ticks = 1
          }
          ticks = as.vector(sapply(m10:M10, function(k) {
            return(ticks * 10^k)
          }))
          lt = length(ticks[ticks > exp(min) & ticks < 
                              exp(max)])
          if (lt < 4) {
            ticks = c(round(exp(min), max(0, -1 * m10 + 
                                            (lt < 2))), ticks, round(exp(max), max(0, 
                                                                                   -1 * M10 + 1 + (lt < 2))))
          }
          image.plot(z = c(min, max), col = Colors, horizontal = T, 
                     legend.only = T, axis.args = list(at = log(ticks), 
                                                       labels = ticks))
        }
        else {
          image.plot(z = as.matrix(rates2), col = Colors, 
                     horizontal = T, legend.only = T)
        }
      }
    }
    par(mfrow = c(1, 1))
  }
}

# clads_TREE(egout, bpal="Spectral", show.tip.label=T, cex=0.1, log=T, tree.type="phylogram")


# process the lambda values by edge
clads_LAMBDA <- function(CladsOutput){
  # identify the parent and child nodes
  edge.df <- data.frame(parent.node = CladsOutput$tree$edge[,1],
                        child.node  = CladsOutput$tree$edge[,2])
  # isolate the edge lengths
  edge.df$branch.length <- CladsOutput$tree$edge.length
  # isolate the lambda values
  edge.df$lambda <- CladsOutput$lambdai_map

  return(edge.df)
}



