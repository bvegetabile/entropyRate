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

p = 0.4
a = 0.05
phi <- a/p - 1
c = a - phi

q = 0.75
d = 0.01
gam = d/q - 1
b = d - gam

print(c(p,a,phi, c))
print(c(d,q,gam, b))

para_set1 <- c(phi, gam)

tm2 <- TwoOrderTransition(p, q, a/p - 1, d/q - 1)
ent2 <- CalcMarkovEntropyRate(tm2, CalcEigenStationary(tm2))
tm1 <- matrix(c(1-p, p, q, 1-q), 2,2, T)
ent1 <- CalcMarkovEntropyRate(tm1, CalcEigenStationary(tm1))

seq_len <- 1000
mc2 <- SimulateTwoOrder(tm2, seq_len)

samp_dist <- matrix(NA, nrow=5000, ncol=4)
for(i in 1:5000){
  mc2 <- SimulateTwoOrder(tm2, seq_len)
  samp_dist[i,] <- c(efficient_mc_er(mc2, 1),
                     efficient_mc_er(mc2, 2),
                     efficient_mc_er(mc2, 3),
                     SWLZEntRate(mc2))
  if(!(i %% 100) ){
    message(paste(i, ':', sep=""), appendLF=F)
  }
}

apply(samp_dist, 2, sd)

set.seed(8924)
sims <- 500
sim_study <- matrix(NA)
sd_res <- matrix(NA, nrow=sims, ncol=4)
mean_res <- matrix(NA, nrow=sims, ncol=4)
est_res <- matrix(NA, nrow=sims, ncol=4)
for(i in 1:sims){
  mc2 <- SimulateTwoOrder(tm2, seq_len)
  est1 <- efficient_mc_er(mc2, 1)
  est2 <- efficient_mc_er(mc2, 2)
  est3 <- efficient_mc_er(mc2, 3)
  est4 <- SWLZEntRate(mc2)
  p1 <- est1 / log2(seq_len)
  p2 <- est2 / log2(seq_len)
  p3 <- est3 / log2(seq_len)
  p4 <- est4 / log2(seq_len)

  B <- 100
  bs_res <- matrix(NA, nrow=B, ncol=4)
  # Bootstrap 1
  for(b in 1:B){
    mc_new1 <- stationary_bootstrap(mc2, p1)
    mc_new2 <- stationary_bootstrap(mc2, p2)
    mc_new3 <- stationary_bootstrap(mc2, p3)
    mc_new4 <- stationary_bootstrap(mc2, p4)

    bs1 <- efficient_mc_er(mc_new1)
    bs2 <- efficient_mc_er(mc_new2, 2)
    bs3 <- efficient_mc_er(mc_new3, 3)
    bs4 <- SWLZEntRate(mc_new4)
    bs_res[b,] <- c(bs1, bs2, bs3, bs4)
  }
  est_res[i, ] <- c(est1, est2, est3, est4)
  mean_res[i, ] <- apply(bs_res, 2, mean)
  sd_res[i, ] <- apply(bs_res, 2, sd)
  if(!(i %% 100) ){
    message(paste(i, ':', sep=""), appendLF=F)
  }
}

apply(samp_dist, 2, sd)
apply(sd_res, 2, median, na.rm=T)

low_res <- est_res  - 1.96*sd_res
upp_res <- est_res  + 1.96*sd_res

cov_res <- low_res  < ent2 & upp_res > ent2
coverage <- apply(cov_res, 2, mean)

pdf('bootstrap_stderr.pdf', height=12, width=12)
par(mfrow=c(2,2))
spot <- 1
plot(density(sd_res[,spot], na.rm=T, from = 0, to = max(sd_res)),
     xlab='Bootstrap Standard Error',
     ylab='Density',
     main=paste('Bootstrap Standard Error : First-Order Markov Chain\nCoverage = ', 
                round(coverage[spot],2), sep=''))
abline(v=mean(sd_res[,spot]), lty=2)
abline(v=apply(samp_dist, 2, sd)[2], lty=3, col=rgb(0.75,0,0,0.75))
legend('topright', c('Simulation Means',
                    'True Std. Err.'), 
       col=c('black', rgb(0.75,0,0,0.75)),
       lty=c(2,3))

spot <- 2
plot(density(sd_res[,spot], na.rm=T, from = 0, to = max(sd_res)),
     xlab='Bootstrap Standard Error',
     ylab='Density',
     main=paste('Bootstrap Standard Error : Second-Order Markov Chain\nCoverage = ', 
                round(coverage[spot],2), sep=''))
abline(v=mean(sd_res[,spot]), lty=2)
abline(v=apply(samp_dist, 2, sd)[2], lty=3, col=rgb(0.75,0,0,0.75))
legend('topright', c('Simulation Means',
                    'True Std. Err.'), 
       col=c('black', rgb(0.75,0,0,0.75)),
       lty=c(2,3))

spot <- 3
plot(density(sd_res[,spot], na.rm=T, from = 0, to = max(sd_res)),
     xlab='Bootstrap Standard Error',
     ylab='Density',
     main=paste('Bootstrap Standard Error : Third-Order Markov Chain\nCoverage = ', 
                round(coverage[spot],2), sep=''))
abline(v=mean(sd_res[,spot]), lty=2)
abline(v=apply(samp_dist, 2, sd)[2], lty=3, col=rgb(0.75,0,0,0.75))
legend('topright', c('Simulation Means',
                    'True Std. Err.'), 
       col=c('black', rgb(0.75,0,0,0.75)),
       lty=c(2,3))

spot <- 4
plot(density(sd_res[,spot], na.rm=T, from = 0, to = max(sd_res)),
     xlab='Bootstrap Standard Error',
     ylab='Density',
     main=paste('Bootstrap Standard Error : Lempel-Ziv\nCoverage = ', 
                round(coverage[spot],2), sep=''))
abline(v=mean(sd_res[,spot]), lty=2)
abline(v=apply(samp_dist, 2, sd)[2], lty=3, col=rgb(0.75,0,0,0.75))
legend('topright', c('Simulation Means',
                    'True Std. Err.'), 
       col=c('black', rgb(0.75,0,0,0.75)),
       lty=c(2,3))

dev.off()

