library("stringr")
library('xtable')
library('dirmult')
library('lattice')
source('/Users/bvegetabile/git/entropyRate/entropyRate.R')
source('/Users/bvegetabile/Dropbox/rCode/bvplot.R')
setwd('/Users/bvegetabile/Dropbox/ucirvine/research/papers/2017_entropyrate/methodpaper/')

################################################################################
#
# Visualization Functions
#
################################################################################

rat.trellis.state.transition <- function(test.case, titl='True Entropy Rate:'){
  N <- nrow(test.case)
  sm <- CalcEigenStationary(test.case)
  ent <- CalcEntropyRate(test.case, sm)
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
                  main=paste(titl)))
  
}


ratfile <- 'rodentdata//Rats.csv'
rats=read.csv(ratfile,header=T)

n.rats <- max(rats$Rat)

results.store <- matrix(NA, ncol=3, nrow=n.rats)
rat.results <- data.frame(matrix(vector(), n.rats, 10, 
                                 dimnames=list(c(), c('RodentID',
                                                      'Control',
                                                      'NumberOfTransitions',
                                                      'SWLZEntropy',
                                                      'M1_EigEntropy',
                                                      'M1_EmpEntropy',
                                                      'M2_EigEntropy',
                                                      'M2_EmpEntropy',
                                                      'MeanLickingGrooming',
                                                      'LGEpisodes'))),
                          stringsAsFactors = F)

for(rat in 1:n.rats){
  message(paste('Rat Number:', rat))
  rat.data <- rats[rats$Rat == rat,]
  rat.behaviors <- rat.data$Behavior #[1:290]
  name_map = levels(rat.behaviors)
  print(name_map)
  licking.grooming <- rat.data[rat.data[,4]=='LG',]
  mean.lg.time <- mean(licking.grooming[,3])
  
  emp_ent1 <- CalcEntropyRate(as.integer(rat.behaviors), 
                              state_space = 1:7, 
                              method = 'Markov', 
                              mc_order=1)
  emp_ent2 <- CalcEntropyRate(as.integer(rat.behaviors), 
                              state_space = 1:7, 
                              method = 'Markov', 
                              mc_order=2)
  eig_ent1 <- CalcEntropyRate(as.integer(rat.behaviors), 
                              state_space = 1:7,
                              method = 'Markov', 
                              stat_method = 'Eigen',  
                              mc_order=1)
  eig_ent2 <- CalcEntropyRate(as.integer(rat.behaviors), 
                              state_space = 1:7,
                              method = 'Markov', 
                              stat_method = 'Eigen', 
                              mc_order=2)
  swlz_ent <- CalcEntropyRate(as.integer(rat.behaviors),
                              state_space = 1:7,
                              method = 'SWLZ')
  
  rat.results$RodentID[rat] <- rat
  print(rat.data$Control[1])
  rat.results$Control[rat] <- as.character(rat.data$Control[1])
  tc <- CalcTransitionCounts(as.integer(rat.behaviors), 7)
  tm <- CalcTransitionMatrix(tc)
  colnames(tm) <- name_map
  rownames(tm) <- name_map
  rat.trellis.state.transition(tm, titl = paste('Estimated Transition Matrix, Rodent:', rat))
  rat.results$NumberOfTransitions[rat] <- length(rat.behaviors) - 1
  rat.results$SWLZEntropy[rat] <- swlz_ent
  rat.results$M1_EigEntropy[rat] <- eig_ent1
  rat.results$M1_EmpEntropy[rat] <- emp_ent1
  rat.results$M2_EigEntropy[rat] <- eig_ent2
  rat.results$M2_EmpEntropy[rat] <- emp_ent2
  rat.results$MeanLickingGrooming[rat] <- mean.lg.time
  rat.results$LGEpisodes[rat] <- nrow(licking.grooming)
}

rat.results$Control[rat.results$Control=='CES'] <- 'LBN'
par(mfrow=c(1,2))
boxplot(M1_EigEntropy ~ Control, data=rat.results, ylim=c(1.4,2.25))
boxplot(M2_EigEntropy ~ Control, data=rat.results, ylim=c(1,1.8))
t.test(SWLZEntropy ~ Control, data=rat.results)
t.test(SWLZEntropy ~ Control, data=rat.results, var.equal = T)
ctl.rats <- rat.results[rat.results$Control=='CTL',]
ces.rats <- rat.results[!rat.results$Control=='CTL',]


# Printing table for paper
print(xtable(rat.results[order(rat.results$Control),c(1,2,3,4,6,8,5,7)], 
             digits = c(0,0,0,0,4,4,4,4,4)), include.rownames=FALSE)

round(apply(tester[1:6,c(4, 5, 6, 7, 8)], 2, mean), 4)
round(apply(tester[7:12,c(4, 5, 6, 7, 8)], 2, mean), 4)

t.test(SWLZEntropy ~ Control, data=rat.results, var.equal = T)
t.test(M1_EmpEntropy ~ Control, data=rat.results, var.equal = T)
t.test(M2_EmpEntropy ~ Control, data=rat.results, var.equal = T)
t.test(M1_EigEntropy ~ Control, data=rat.results, var.equal = T)
t.test(M2_EigEntropy ~ Control, data=rat.results, var.equal = T)
