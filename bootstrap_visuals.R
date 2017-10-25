################################################################################
# Loading packages
################################################################################
library("stringr")
library('xtable')
library('dirmult')
library('lattice')
source('/Users/bvegetabile/git/entropyRate/entropyRate.R')
setwd('/Users/bvegetabile/Dropbox/ucirvine/research/papers/2017_entropyrate/2017-10-JEBS/')

summary1 <- function(sims){
  n <- nrow(sims)
  sorted <- apply(sims, 2, sort)
  lower.values <- sorted[1,]
  upper.values <- sorted[n,]
  mean.values <- apply(sorted, 2, mean)
  median.values <- apply(sorted, 2, median)
  out.put <- cbind(mean.values, median.values, lower.values, upper.values)
  return(out.put)
}


summary2 <- function(sims){
  n <- nrow(sims)
  sorted <- apply(sims, 2, sort)
  lower.values <- unlist(lapply(sorted, min))
  upper.values <- unlist(lapply(sorted, max))
  mean.values <- apply(sims, 2, mean, na.rm=T)
  median.values <- apply(sims, 2, median, na.rm=T)
  out.put <- cbind(mean.values, median.values, lower.values, upper.values)
  return(out.put)
}

summarize.sims <- function(sims){
  if(sum(is.na(sims)) > 0){
    return(summary2(sims))
  } else{
    return(summary1(sims))
  }
}

plot.simresults <- function(i, titl){
  titles <- c('Empirical Stationary', 'Eigenvector Stationary', 'SWLZ')
  par(mfrow=c(1,3), oma=c(0,0,0,0), mar=c(5,4,3,0.5))
  fname <- paste('./simulation_oct2017.', i, '.RDS', sep='')
  sim.data <- readRDS(fname)
  tru <- sim.data[['True Entropy']]
  tab <- cbind(summarize.sims(sim.data[[7]]), 
               summarize.sims(sim.data[[8]]), 
               summarize.sims(sim.data[[9]]))
  emp_stderr <- apply(sim.data[[4]], 2, sd)
  for(j in 1:3){
    pts <- tab[,((j-1)*4 + 1):(j*4)]
    plot(1, xlim = c(0,0.5), ylim=c(0.75,5.25), type="n", axes=F,
         xlab="Est. Entropy", ylab="Chain Length",
         main=paste(titl, titles[j], '\nBootstrap Std. Errs', sep=" "))
    top <- nrow(pts) + 1
    for(k in 1:nrow(pts)){
      lines(c(pts[k,3], pts[k,4]), rep(top-k, 2))
      points(c(pts[k,1], pts[k,3], pts[k,4]), rep(top-k, 3), pch=4)
      points(emp_stderr[k], top-k, pch=19, col=rgb(0,0,0.75,0.5), cex=1.5)
    }
    # abline(v=tru, lty=3)
    axis(2, at=1:5, labels = c(5000,1000,500,250,50), las=2)
    axis(1, at=c(0,0.5,0.25))
    abline(v=0, lty=3)
  }
}


pdf('2017-10-24_lowent_stderr.pdf', height=2.5, width=10)
plot.simresults(4, 'Low')
dev.off()

pdf('2017-10-24_medent_stderr.pdf', height=2.5, width=10)
plot.simresults(5, 'Med.')
dev.off()

pdf('2017-10-24_highent_stderr.pdf', height=2.5, width=10)
plot.simresults(6, 'High')
dev.off()


plot.simresults(7, 'Periodic')


i <- 7
fname <- paste('./simulation_oct2017.', i, '.RDS', sep='')
sim.data <- readRDS(fname)
tru <- sim.data[['True Entropy']]
tab <- cbind(summarize.sims(sim.data[[7]]), 
             summarize.sims(sim.data[[8]]), 
             summarize.sims(sim.data[[9]]))
emp_stderr <- apply(sim.data[[4]], 2, sd)

plot(sim.data[[4]], sim.data[[7]])



plot_sims <- function(i,titl){
  titles <- c('Empirical Stationary', 'Eigenvector Stationary', 'SWLZ')
  par(mfrow=c(3,2), oma=c(0,0,3,0), mar=c(5,6,3,0.5))
  fname <- paste('./simulation_oct2017.', i, '.RDS', sep='')
  sim.data <- readRDS(fname)
  
  tru <- sim.data[['True Entropy']]
  
  tab_ests <- sim.data[[2]]
  
  tab_stderrs <- cbind(summarize.sims(sim.data[[7]]), 
               summarize.sims(sim.data[[8]]), 
               summarize.sims(sim.data[[9]]))
  
  emp_stderr <- apply(sim.data[[4]], 2, sd)
  
  
  for(j in 1:3){
    plot_ests(j, tab_ests, tru, titl, titles, sidelab = titles[j])
    plot_stderr(j, tab_stderrs, titl, titles, emp_stderr)
  }
  mtext(titl, side=3, outer=T)
  
}

plot_ests <- function(j, tab, tru, titl, titles, sidelab){
  pts <- tab[,((j-1)*4 + 1):(j*4)]
  plot(1, xlim = c(0,max(3,max(pts))), ylim=c(0.75,5.25), type="n", axes=F, 
       xlab="Est. Entropy Rate", ylab="Chain Length",
       main=paste('Entropy Rate Estimate', sep=" "))
  top <- nrow(pts) + 1
  for(k in 1:nrow(pts)){
    lines(c(pts[k,3], pts[k,4]), rep(top-k, 2))
    points(c(pts[k,1], pts[k,3], pts[k,4]), rep(top-k, 3), pch=4)
  }
  abline(v=tru, lty=3)
  axis(2, at=1:5, labels = c(5000,1000,500,250,50), las=2)
  axis(1, at=c(0,1.5,3))
  mtext(sidelab, side=2, line = 4)
}

plot_stderr <- function(j, tab, titl, titles, emp_stderr){
  pts <- tab[,((j-1)*4 + 1):(j*4)]
  plot(1, xlim = c(0,0.5), ylim=c(0.75,5.25), type="n", axes=F,
       xlab="Standard Error", ylab="Chain Length",
       main=paste('Bootstrap Std. Error', sep=" "))
  top <- nrow(pts) + 1
  for(k in 1:nrow(pts)){
    lines(c(pts[k,3], pts[k,4]), rep(top-k, 2))
    points(c(pts[k,1], pts[k,3], pts[k,4]), rep(top-k, 3), pch=4)
    points(emp_stderr[k], top-k, pch=19, col=rgb(0,0,0.75,0.5), cex=1.5)
  }
  # abline(v=tru, lty=3)
  axis(2, at=1:5, labels = c(5000,1000,500,250,50), las=2)
  axis(1, at=c(0,0.5,0.25))
  abline(v=0, lty=3)
}

pdf('./figures/2017-10-24_lowent.pdf', height=7, width=4.5)
plot_sims(4, 'Low Entropy Rate')
dev.off()

pdf('./figures/2017-10-24_medent.pdf', height=7, width=4.5)
plot_sims(5, 'Medium Entropy Rate')
dev.off()

pdf('./figures/2017-10-24_highent.pdf', height=7, width=4.5)
plot_sims(6, 'High Entropy Rate')
dev.off()