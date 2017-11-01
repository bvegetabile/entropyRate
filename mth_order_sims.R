setwd('/Users/bvegetabile/Dropbox/ucirvine/research/papers/2017_entropyrate/2017-10-JEBS/figures/')
source('/Users/bvegetabile/git/entropyRate/entropyRate.R')

#-------------------------------------------------------------------------------
# Setup and functions for second order process----------------------------------

TwoOrderTransition <- function(p,q,phi=0,gam=0){
  state_space <- c("11", "12", "21", "22")
  tm <- matrix(c((1 + phi)*(1/(1+phi) - p), p*(1+phi), 0,0,
                 0,0,(1+gam)*(q - gam/(1+gam)),(1-q)*(1+gam),
                 (1-p)*(1+phi),(1+phi)*(p-phi/(1+phi)),0,0,
                 0,0,q*(1+gam),(1+gam)*(1/(1+gam)-q)),
               4,4, T)
  colnames(tm) <- state_space
  rownames(tm) <- state_space
  return(tm)
}

SimulateTwoOrder <- function(tm, n_sims=100, ss=c("11", "12", "21", "22")){
  so_seq <- SimulateMarkovChain(tm, n_sims = n_sims)
  so_symb <- ss[so_seq]
  mc_chain <- c(as.integer(substr(so_symb[1], 1,1)), 
                as.integer(substr(so_symb[1], 2,2))) 
  for(i in 2:n_sims){
    mc_chain <- c(mc_chain, as.integer(substr(so_symb[i], 2,2)))
  }
  return(mc_chain)
}

#-------------------------------------------------------------------------------
# Plot to find valid pairs of p and a, as well as q and d ----------------------
n <- 500
a <- seq(0,1,length.out = n)
p <- seq(0,1,length.out = n)
phi_result <- matrix(NA, n, n)
c_result <- matrix(NA, n, n)
for(i in 1:n){
  for(j in 1:n){
    phi_result[i,j] <- a[i]/p[j] - 1
    c_result[i,j] <- a[i] - phi_result[i,j]
  }
}
par(mfrow=c(1,2))
image(a,p,phi_result, zlim = c(-1,1), col=heat.colors(1000))
contour(a,p,phi_result, zlim = c(-1,1), add=T, levels=seq(-1,1,0.1))
abline(h=seq(-.9,.9,0.1), lty=3)
abline(v=seq(-.9,.9,0.1), lty=3)
image(a,p,c_result, zlim = c(0,1), col=heat.colors(1000))
contour(a,p,c_result, zlim = c(0,1), add=T)
abline(h=seq(-.9,.9,0.1), lty=3)
abline(v=seq(-.9,.9,0.1), lty=3)

#-------------------------------------------------------------------------------
# Simulation 1 -----------------------------------------------------------------
p = 0.4
a = 0.1
phi <- a/p - 1
c = a - phi

q = 0.75
d = 0.2
gam = d/q - 1
b = d - gam

print(c(p,a,phi, c))
print(c(d,q,gam, b))

para_set1 <- c(phi, gam)

so_test_one_sim1 <- TwoOrderTransition(p, q, a/p - 1, d/q - 1)
fo_test_one_sim1 <- matrix(c(1-p, p, q, 1-q), 2,2, T)

true_ent_sim1 <- CalcMarkovEntropyRate(so_test_one_sim1, 
                                       CalcEigenStationary(so_test_one_sim1))

set.seed(124)
n_sims <- 1000
results_mat_sim1 <- matrix(NA, nrow=n_sims, ncol=4)
state_space <- c('1', '2')
for(i in 1:n_sims){
  mcchain <- SimulateTwoOrder(so_test_one_sim1, 1000)
  ent1 <- CalcEntropyRate(mcchain, state_space, mc_order = 1, stat_method = 'Empirical')
  ent2 <- CalcEntropyRate(mcchain, state_space, mc_order = 2, stat_method = 'Empirical')
  ent3 <- CalcEntropyRate(mcchain, state_space, mc_order = 3, stat_method = 'Empirical')
  ent4 <- CalcEntropyRate(mcchain, state_space, method = 'SWLZ')
  results_mat_sim1[i, ] <- c(ent1, ent2, ent3, ent4)[4:1]
}
colnames(results_mat_sim1) <- c('m=1', 'm=2', 'm=3', 'SWLZ')[4:1]
boxplot(results_mat_sim1, horizontal=T, ylim=c(0.25, 1))
abline(v=true_ent_sim1)
#-------------------------------------------------------------------------------
# Simulation 2 -----------------------------------------------------------------
p = 0.4
a = 0.52
phi <- a/p - 1
c = a - phi

q = 0.75
d = 0.95
gam = d/q - 1
b = d - gam

para_set2 <- c(phi, gam)

print(c(p,a,phi, c))
print(c(d,q,gam, b))

so_test_one_sim2 <- TwoOrderTransition(p, q, a/p - 1, d/q - 1)
fo_test_one_sim2 <- matrix(c(1-p, p, q, 1-q), 2,2, T)

true_ent_sim2 <- CalcMarkovEntropyRate(so_test_one_sim2, 
                                       CalcEigenStationary(so_test_one_sim2))

set.seed(124)
n_sims <- 1000
results_mat_sim2 <- matrix(NA, nrow=n_sims, ncol=4)
state_space <- c('1', '2')
for(i in 1:n_sims){
  mcchain <- SimulateTwoOrder(so_test_one_sim2, 1000)
  ent1 <- CalcEntropyRate(mcchain, state_space, mc_order = 1, stat_method = 'Empirical')
  ent2 <- CalcEntropyRate(mcchain, state_space, mc_order = 2, stat_method = 'Empirical')
  ent3 <- CalcEntropyRate(mcchain, state_space, mc_order = 3, stat_method = 'Empirical')
  ent4 <- CalcEntropyRate(mcchain, state_space, method = 'SWLZ')
  results_mat_sim2[i, ] <- c(ent1, ent2, ent3, ent4)[4:1]
}
colnames(results_mat_sim2) <- c('m=1', 'm=2', 'm=3', 'SWLZ')[4:1]
boxplot(results_mat_sim2, horizontal=T, ylim=c(0.25,1))
abline(v = true_ent_sim2, lty=3)
#-------------------------------------------------------------------------------
# Simulation 3 -----------------------------------------------------------------
p = 0.4
a = 0.01
phi <- a/p - 1
c = a - phi

q = 0.75
d = 0.01
gam = d/q - 1
b = d - gam

para_set3 <- c(phi, gam)

print(c(p,a,phi, c))
print(c(d,q,gam, b))

so_test_one_sim3 <- TwoOrderTransition(p, q, a/p - 1, d/q - 1)
fo_test_one_sim3 <- matrix(c(1-p, p, q, 1-q), 2,2, T)

true_ent_sim3 <- CalcMarkovEntropyRate(so_test_one_sim3, 
                                       CalcEigenStationary(so_test_one_sim3))

set.seed(124)
n_sims <- 1000
results_mat_sim3 <- matrix(NA, nrow=n_sims, ncol=4)
state_space <- c('1', '2')
for(i in 1:n_sims){
  mcchain <- SimulateTwoOrder(so_test_one_sim3, 1000)
  ent1 <- CalcEntropyRate(mcchain, state_space, mc_order = 1, stat_method = 'Empirical')
  ent2 <- CalcEntropyRate(mcchain, state_space, mc_order = 2, stat_method = 'Empirical')
  ent3 <- CalcEntropyRate(mcchain, state_space, mc_order = 3, stat_method = 'Empirical')
  ent4 <- CalcEntropyRate(mcchain, state_space, method = 'SWLZ')
  results_mat_sim3[i, ] <- c(ent1, ent2, ent3, ent4)[4:1]
}
colnames(results_mat_sim3) <- c('m=1', 'm=2', 'm=3', 'SWLZ')[4:1]
boxplot(results_mat_sim3, horizontal=T, ylim=c(0, 1))
abline(v=true_ent_sim3)

#-------------------------------------------------------------------------------
# Plotting Simulation Results --------------------------------------------------
# pdf('second_order_sims.pdf', height = 5, width=10)
par(mfrow=c(2,1), mar=c(3,4,2,1)+0.1)
boxplot(results_mat_sim1, horizontal=T, ylim=c(0.25, 1), axes=F, pch=19, 
        pars=list(outcol=rgb(0,0,0,0.25)),
        xlab='Entropy Rate Estimate\n',
        main=expression(paste('a = 0.1, c = 0.85, (a-c = -0.75) and d = 0.2, b = 0.933, (d-b = -.733)')))
abline(v=true_ent_sim1, lty=3, lwd=2, col=rgb(0.75,0,0,1))
axis(1, at=seq(0.3,1,0.1))
axis(2, 1:4, labels = colnames(results_mat_sim1), las=2)
boxplot(results_mat_sim2, horizontal=T, ylim=c(0.25,1), axes=F, pch=19,
        pars=list(outcol=rgb(0,0,0,0.25)),
        xlab='Entropy Rate Estimate\n',
        main=expression(paste('a = 0.52, c = 0.22 (a-c = 0.3) and d = 0.95, b = 0.6833, (d-b = 0.2667)')))
abline(v = true_ent_sim2, lty=3, lwd=2, col=rgb(0.75,0,0,1))
axis(1, at=seq(0.3,1,0.1))
axis(2, 1:4, labels = colnames(results_mat_sim2), las=2)
# dev.off()
#-------------------------------------------------------------------------------
# Plotting Analytical Differences ----------------------------------------------
plotEntropyDifference <- function(p, q,phi_res=100, gam_res=100, case_name=''){
  n_phi <- phi_res
  n_gam <- gam_res
  
  phi_max <- min(1, (1-p)/p, p/(1-p))
  gam_max <- min(1, (1-q)/q, q/(1-q))
  
  phi_test <- seq(-0.999, phi_max-0.001, length.out = n_phi)
  gam_test <- seq(-0.999, gam_max-0.001, length.out = n_gam)
  
  ent1_grid <- matrix(NA, nrow = n_phi, ncol=n_gam)
  ent2_grid <- matrix(NA, nrow = n_phi, ncol=n_gam)
  
  for(i in 1:length(phi_test)){
    for(j in 1:length(gam_test)){
      phi_i <- phi_test[i]
      gam_j <- gam_test[j]
      oneorder <- matrix(c(1-p, p,q, 1-q), 2,2,T)
      twoorder <- TwoOrderTransition(p,q,phi_i,gam_j)
      sm1 <- CalcEigenStationary(oneorder)
      sm2 <- CalcEigenStationary(twoorder)
      ent1 <- CalcMarkovEntropyRate(oneorder, sm1)
      ent2 <- CalcMarkovEntropyRate(twoorder, sm2)
      ent1_grid[i,j] <- ent1
      ent2_grid[i,j] <- ent2
    }
  }
  
  print(max(ent1_grid-ent2_grid))
  
  n_probs <- 1000
  p_tests <- seq(0,1,length.out = n_probs)
  prob_vals <- rep(NA, n_probs)
  for(k in 1:n_probs){
    prob_vals[k] <- sum((ent1_grid - ent2_grid) < p_tests[k]) / (n_phi*n_gam)
  }
  image(phi_test, gam_test, ent2_grid, 
        zlim = c(0,1),
        col=heat.colors(1000, alpha=0),
        xlab=expression(phi),
        ylab=expression(gamma), 
        main=paste('Entropy Rate, p = ', p, ', q = ', q,sep=''))
  contour(phi_test, gam_test, ent2_grid, add=T, col=rgb(0,0,0,0.5))
  abline(h = seq(-1,1,0.25), lty=3, col=rgb(0,0,0,0.5))
  abline(v = seq(-1,1,0.25), lty=3, col=rgb(0,0,0,0.5))
}

pdf('so_analytical_diff.pdf', height=5, width=6)
par(mfcol=c(1,1), mar=c(5,4,2,1)+0.1)
plotEntropyDifference(0.4, 0.75, 100, 100, '')
points(para_set1[1], para_set1[2], pch=18, lwd=2, cex=2)
text(para_set1[1], para_set1[2], "I", cex=1.5, pos = 2, family='Times')
points(para_set2[1], para_set2[2], pch=18, lwd=2, cex=2)
text(para_set2[1], para_set2[2], "II", cex=1.5, pos = 2, family='Times')
dev.off()
# 
#-------------------------------------------------------------------------------