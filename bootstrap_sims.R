################################################################################
# 
# Simulations for the paper: Estimating the Entropy Rate of Finite Markov Chains 
#                            with Application to Behavior Studies
#                            
# Paper Authors            : Brian Vegetabile, Jenny Molet, Tallie Z. Baram, Hal Stern
# Code Written by:         : Brian Vegetabile (bvegetab [ATSYMBOL] uci [dot] edu)
# 
################################################################################

################################################################################
# Loading packages
################################################################################
library("stringr")
library('xtable')
library('dirmult')
library('lattice')
source('/Users/bvegetabile/git/entropyRate/entropyRate.R')
setwd('/Users/bvegetabile/Dropbox/ucirvine/research/papers/2017_entropyrate/2017-10-JEBS/')

################################################################################
# Visualization Functions for simulations
################################################################################

trellis.state.transition <- function(test.case, titl='True Entropy Rate:'){
  N <- nrow(test.case)
  row.names(test.case) <- paste("State", seq(1,N))
  colnames(test.case) <- paste("State", seq(1,N))
  sm <- CalcEigenStationary(test.case)
  ent <- CalcMarkovEntropyRate(test.case, sm)
  trellis.par.set(regions=list(col=colorRampPalette(c('white',
                                                      rgb(0.75,0,0,0.5),
                                                      'blue', 
                                                      "black"))(1000)))
  lattice.options(axis.padding=list(factor=0.5))
  print(levelplot(t(test.case[seq(N,1),]), cuts=1000, pretty=TRUE,
                  panel = function(...){
                    panel.levelplot(...)
                    panel.abline(h = seq(1.5,7.5, 1), col=rgb(0,0,0,.75))
                    panel.abline(v = seq(1.5,7.5, 1), col=rgb(0,0,0,.75))
                  },
                  at=seq(0,1,0.01),
                  xlab='', ylab='', #scales=list(x=list(rot=360)), 
                  main=paste(titl, round(ent,4))))
  
}

summarize.sims <- function(sims){
  n <- nrow(sims)
  sorted <- apply(sims, 2, sort)
  lower.values <- sorted[1,]
  upper.values <- sorted[n,]
  mean.values <- apply(sorted, 2, mean)
  median.values <- apply(sorted, 2, median)
  out.put <- cbind(mean.values, median.values, lower.values, upper.values)
  return(out.put)
}



################################################################################
#
# Creating simulations for paper
#
################################################################################

test.cases <- list()

################################################################################
## TWO STATES
################################################################################

state_space <- c(1,2)
###
# Two States - Low Entropy 
p10 <- 0.05
p01 <- 0.01 
test.cases[['two.state.lowentropy']] <- matrix(c(1-p10, p01, p10, 1-p01),2,2)

###
# Two States - Medium Entropy
p10 <- 0.1 
p01 <- 0.8
test.cases[['two.state.medentropy']] <- matrix(c(1-p10, p01, p10, 1-p01),2,2)

###
# Two States - High Entropy 
p10 <- 0.45
p01 <- 0.3 
test.cases[['two.state.highentropy']] <- matrix(c(1-p10, p01, p10, 1-p01),2,2)

################################################################################
## EIGHT STATES
################################################################################
# state_space <- seq(1,8)
# 
# ### 
# # Eight States - Low Entropy

test.cases[['eight.state.lowentropy']] <- matrix( c( .95, 0.05,   0,   0,   0,   0,   0,   0,
                                                     .025,  .95, .025,   0,   0,   0,   0,   0,
                                                     0, .025,  .95, .025,   0,   0,   0,   0,
                                                     0,   0, .025,  .95, .025,   0,   0,   0,
                                                     0,   0,   0, .025,  .95, .025,   0,   0,
                                                     0,   0,   0,   0, .025,  .95, .025,   0,
                                                     0,   0,   0,   0,   0, .025,  .95, .025,
                                                     0,   0,   0,   0,   0,   0,  .05,  .95), 
                                                  8,8, byrow=TRUE)
sample(x=seq(1,8), 8, replace=T)
# ###
# # Eight States - Medium Entropy 
set.seed(2015)
med.entropy <- matrix(NA, 8,8)
# for(i in 1:8){
#   med.entropy[i,] <- rdirichlet(n=1, alpha = sample(x=(1:8)^6, 8, replace=F))
# }
upper <- 6.3

med.entropy[1,] <- rdirichlet(n=1, alpha = sample(1:8, 8, replace = T)^upper)
med.entropy[2,] <- rdirichlet(n=1, alpha = sample(1:8, 8, replace = T)^upper)
med.entropy[3,] <- rdirichlet(n=1, alpha = sample(1:8, 8, replace = T)^upper)
med.entropy[4,] <- rdirichlet(n=1, alpha = sample(1:8, 8, replace = T)^upper)
med.entropy[5,] <- rdirichlet(n=1, alpha = sample(1:8, 8, replace = T)^upper)
med.entropy[6,] <- rdirichlet(n=1, alpha = sample(1:8, 8, replace = T)^upper)
med.entropy[7,] <- rdirichlet(n=1, alpha = sample(1:8, 8, replace = T)^upper)
med.entropy[8,] <- rdirichlet(n=1, alpha = sample(1:8, 8, replace = T)^upper)
test.cases[['eight.state.medentropy']] <- med.entropy

# ###
# # Eight States - High Entropy 
set.seed(2015)
test.cases[['eight.state.highentropy']] <- rdirichlet(n=8, alpha = c(rep(5,8)))

# ###
# # Eight States - Periodic Markov Chain
test.cases[['eight.state.periodic']] <- matrix( c( 0, 1/3, 1/3, 1/3, 0, 0, 0, 0,
                                                   1/3, 0, 0, 0, 0, 1/3, 1/3, 0,
                                                   1/3, 0, 0, 0, 1/3, 0, 1/3, 0,
                                                   1/3, 0, 0, 0, 1/3, 1/3, 0, 0,
                                                   0, 0, 1/3, 1/3, 0, 0, 0, 1/3,
                                                   0, 1/3, 0, 1/3, 0, 0, 0, 1/3,
                                                   0, 1/3, 1/3, 0, 0, 0, 0, 1/3,
                                                   0, 0, 0, 0, 1/3, 1/3, 1/3, 0), 
                                                8,8, byrow=TRUE)

# pdf('2017_sims/eightstate_low_entropy_plot.pdf', height=6, width=6)
trellis.state.transition(test.cases[[4]], 'Low Entropy Case, Entropy Rate:')
# dev.off()
# pdf('2017_sims/eightstate_med_entropy_plot.pdf', height=6, width=6)
trellis.state.transition(test.cases[[5]], 'Med. Entropy Case, Entropy Rate:')
# dev.off()
# pdf('2017_sims/eightstate_high_entropy_plot.pdf', height=6, width=6)
trellis.state.transition(test.cases[[6]], 'High Entropy Case, Entropy Rate:')
# dev.off()
# pdf('2017_sims/eightstate_periodic_plot.pdf', height=6, width=6)
trellis.state.transition(test.cases[[7]], 'Period Entropy Case, Entropy Rate:')
# dev.off()

par(mfrow=c(1,1))
set.seed(23124)
for(test in 6:7){
  print(paste('Simulation ', test, ' Started at: ', Sys.time()))
  nsimulations <- 100
  n_bootstrap <- 100
  
  trans_mat <- test.cases[[test]]
  nstates <- nrow(trans_mat)
  state_space <- seq(1:nstates)
  true.entropy.rate <- CalcMarkovEntropyRate(trans_mat, CalcEigenStationary(trans_mat))
  
  test.lengths <- c(50, 250, 500, 1000, 5000)
  n_tl <- length(test.lengths)
  
  max_sim_length = max(test.lengths)
  
  sample.chains <- matrix(0, nrow=nsimulations, ncol=max_sim_length)
  empirical.entropy.matrix <- matrix(NA, nrow=nsimulations, ncol=n_tl)
  lempel.ziv.matrix <- matrix(NA, nrow=nsimulations, ncol=n_tl)
  eigen.entropy.matrix <- matrix(0, nrow=nsimulations, ncol=n_tl)
  
  stderr.empirical.entropy.matrix <- matrix(NA, nrow=nsimulations, ncol=n_tl)
  stderr.lempel.ziv.matrix <- matrix(NA, nrow=nsimulations, ncol=n_tl)
  stderr.eigen.entropy.matrix <- matrix(NA, nrow=nsimulations, ncol=n_tl)
  
  coverage.empirical.entropy.matrix <- matrix(NA, nrow=nsimulations, ncol=n_tl)
  coverage.lempel.ziv.matrix <- matrix(NA, nrow=nsimulations, ncol=n_tl)
  coverage.eigen.entropy.matrix <- matrix(NA, nrow=nsimulations, ncol=n_tl)
  
  run_start <- Sys.time()
  for (i in 1:nsimulations){
    if (!(i %% 10)){
      message(paste("Simulation Number:", i))
      print(Sys.time()-run_start)
    }
    chain <- SimulateMarkovChain(trans_mat, n_sims = max_sim_length)
    sample.chains[i,] <- chain
    for(c in 1:length(test.lengths)){
      n <- test.lengths[c]
      ss <- chain[1:n]
      lz_entropyrate <- CalcEntropyRate(ss, 1:8, method = "SWLZ")
      tc <- CalcTransitionCounts(ss, nstates)
      tm <- CalcTransitionMatrix(tc)
      eig.sm <- CalcEigenStationary(tm)
      if(sum(is.na(eig.sm))>0){message(paste('Sim: ', i, ', Test Length:', n, sep=''))}
      emp.sm <- CalcEmpiricalStationary(ss, 1:nstates)
      emp.ent <- CalcMarkovEntropyRate(tm, emp.sm)
      eig.ent <- CalcMarkovEntropyRate(tm, eig.sm)
      
      
      bootstrap_ests <- matrix(NA, nrow=n_bootstrap, ncol=3)
      # Bootstrapping ----------------------------------------------------------
      for(b in 1:n_bootstrap){
        
        if(emp.ent != 0){
          p.emp <- emp.ent / log2(n)
          bs.emp <- stationary_bootstrap(ss, p.emp)
          bs.emp.tc <- CalcTransitionCounts(bs.emp, nstates)
          bs.emp.tm <- CalcTransitionMatrix(bs.emp.tc)
          bs.emp.sm <- CalcEmpiricalStationary(bs.emp, 1:nstates)
          bs.emp.ent <- CalcMarkovEntropyRate(bs.emp.tm, bs.emp.sm)  
        } else {
          bs.emp.ent <- NA
        }
        
        if(eig.ent != 0){
          p.eig <- eig.ent / log2(n)
          bs.eig <- stationary_bootstrap(ss, p.eig)
          bs.eig.tc <- CalcTransitionCounts(bs.eig, nstates)
          bs.eig.tm <- CalcTransitionMatrix(bs.eig.tc)
          bs.eig.sm <- CalcEigenStationary(bs.eig.tm)
          if(sum(is.na(bs.eig.sm))>0){message(paste('Sim: ', i, ', Bootstrap :', b, ', Test Length:', n, sep=''))}
          bs.eig.ent <- CalcMarkovEntropyRate(bs.eig.tm, bs.eig.sm)          
        } else {
          bs.eig.ent <- NA
        }
        
        bs.lz.ent <- CalcEntropyRate(bs.lz, 1:8, method = "SWLZ")
        p.lz <- lz_entropyrate / log2(n)
        bs.lz <- stationary_bootstrap(ss, p.lz)
        
        
        bootstrap_ests[b, ] <- c(bs.emp.ent, bs.eig.ent, bs.lz.ent)
      }
      bs_stderrs <- apply(bootstrap_ests, 2, sd, na.rm=T)
      
      emp_low <- emp.ent - 1.96 * bs_stderrs[1]
      eig_low <- eig.ent - 1.96 * bs_stderrs[2]
      lz_low <- lz_entropyrate - 1.96 * bs_stderrs[3]
      
      emp_upp <- emp.ent + 1.96 * bs_stderrs[1]
      eig_upp <- eig.ent + 1.96 * bs_stderrs[2]
      lz_upp <- lz_entropyrate + 1.96 * bs_stderrs[3]
      
      coverage.empirical.entropy.matrix[i,c] <- (emp_low < true.entropy.rate & emp_upp > true.entropy.rate)
      coverage.lempel.ziv.matrix[i,c] <- (lz_low < true.entropy.rate & lz_upp > true.entropy.rate)
      coverage.eigen.entropy.matrix[i,c] <- (eig_low < true.entropy.rate & eig_upp > true.entropy.rate)
      
      stderr.empirical.entropy.matrix[i,c] <- bs_stderrs[1]
      stderr.eigen.entropy.matrix[i,c] <- bs_stderrs[2]
      stderr.lempel.ziv.matrix[i,c] <- bs_stderrs[3]
      
      lempel.ziv.matrix[i, c] <- lz_entropyrate
      eigen.entropy.matrix[i, c] <- eig.ent
      empirical.entropy.matrix[i, c] <- emp.ent
    }  
  }
  run_dur <- Sys.time() - run_start
  message(paste('Script total run time: ', 
                round(as.numeric(run_dur, units='mins'),3), 'minutes'))
  
  simulation.table <- cbind(summarize.sims(empirical.entropy.matrix), 
                            summarize.sims(eigen.entropy.matrix), 
                            summarize.sims(lempel.ziv.matrix))
  
  latex.table <- xtable(simulation.table, digits = 4)
  simulation.output <- list()
  simulation.output[['True Entropy']] <- true.entropy.rate
  simulation.output[['Simulation Table']] <- simulation.table
  simulation.output[['LaTeX Table']] <- latex.table
  simulation.output[['Empirical Entropies']] <- empirical.entropy.matrix
  simulation.output[['Eigen Entropies']] <- eigen.entropy.matrix
  simulation.output[['Lempel Ziv']] <- lempel.ziv.matrix
  
  simulation.output[['Std.Err. Empirical Entropies']] <- stderr.empirical.entropy.matrix
  simulation.output[['Std.Err. Eigen Entropies']] <- stderr.eigen.entropy.matrix
  simulation.output[['Std.Err. Lempel Ziv']] <- stderr.lempel.ziv.matrix
  
  simulation.output[['Coverage Empirical Entropies']] <- coverage.empirical.entropy.matrix
  simulation.output[['Coverage Eigen Entropies']] <- coverage.eigen.entropy.matrix
  simulation.output[['Coverage Lempel Ziv']] <- coverage.lempel.ziv.matrix
  
  simulation.output[['Simulated Chains']] <- sample.chains
  saveRDS(simulation.output, file=paste('./simulation_oct2017.', test, '.RDS', sep=''))
  print(paste('Simulation ', test, ' Finished at: ', Sys.time()))
}
# 
plot.simresults <- function(i, titl){
  titles <- c('Empirical Stationary', 'Eigenvector Stationary', 'SWLZ')
  par(mfrow=c(1,3), oma=c(0,0,0,0), mar=c(5,4,3,0.5))
  fname <- paste('./simulation_may2017.', i, '.RDS', sep='')
  sim.data <- readRDS(fname)
  tru <- sim.data[['True Entropy']]
  tab <- sim.data[['Simulation Table']]
  for(j in 1:3){
    pts <- tab[,((j-1)*4 + 1):(j*4)]
    plot(1, xlim = c(0,max(3,max(pts))), ylim=c(0.75,6.25), type="n", axes=F,
         xlab="Est. Entropy", ylab="Chain Length",
         main=paste(titl, titles[j], sep=" "))
    top <- nrow(pts) + 1
    for(k in 1:nrow(pts)){
      lines(c(pts[k,3], pts[k,4]), rep(top-k, 2))
      points(c(pts[k,1], pts[k,3], pts[k,4]), rep(top-k, 3), pch=4)
    }
    abline(v=tru, lty=3)
    axis(2, at=1:6, labels = c(50000, 5000,1000,500,250,50), las=2)
    axis(1, at=c(0,1.5,3))
  }
}

plot_stderrs <- function(i, titl){
  titles <- c('Empirical Stationary', 'Eigenvector Stationary', 'SWLZ')
  par(mfrow=c(1,3), oma=c(0,0,0,0), mar=c(5,4,3,0.5))
  fname <- paste('./simulation_oct2017.', i, '.RDS', sep='')
  sim.data <- readRDS(fname)
  tru <- sim.data[['True Entropy']]
  tab <- sim.data[['Simulation Table']]
  for(j in 1:3){
    pts <- tab[,((j-1)*4 + 1):(j*4)]
    plot(1, xlim = c(0,max(3,max(pts))), ylim=c(0.75,6.25), type="n", axes=F,
         xlab="Est. Entropy", ylab="Chain Length",
         main=paste(titl, titles[j], sep=" "))
    top <- nrow(pts) + 1
    for(k in 1:nrow(pts)){
      lines(c(pts[k,3], pts[k,4]), rep(top-k, 2))
      points(c(pts[k,1], pts[k,3], pts[k,4]), rep(top-k, 3), pch=4)
    }
    abline(v=tru, lty=3)
    axis(2, at=1:6, labels = c(50000, 5000,1000,500,250,50), las=2)
    axis(1, at=c(0,1.5,3))
  }
}


