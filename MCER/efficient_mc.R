setwd('~/git/entropyRate/MCER/')
source('~/git/entropyRate/entropyRate.R')

# Functions --------------------------------------------------------------------

efficient_mc_er <- function(mc_seq, mc_order=1){
  unique_states <- as.character(unique(mc_seq))
  n_states <- length(unique_states)
  state_space <- data.frame(matrix(NA, nrow=n_states, ncol=mc_order))
  for(i in 1:mc_order){
    state_space[,i] <- unique_states
  }
  state_space <- c(expand.grid(state_space), sep="")
  state_space <- do.call(paste, state_space)
  char_seq <- paste(mc_seq, collapse = '')
  # print(unique_states)
  mcer(char_seq, mc_order, unique_states, state_space)
}

fixedwindow_lz77 <- function(mc_seq, window_size = 10){
  char_seq <- paste(mc_seq, collapse = '')
  lz77entropy(char_seq, window_size)
}

fixedquick_lz77 <- function(mc_seq, window_size = 10){
  char_seq <- paste(mc_seq, collapse = '')
  lz77entropy_quick(char_seq, window_size)
}

block_bootstrap <- function(mc_seq, block_size){
  n_blocks <- floor(length(mc_seq) / block_size)
  new_seq <- rep(NA, length = n_blocks * block_size)
  blocks <- sample(1:n_blocks, replace = T)
  for(i in 1:n_blocks){
    beg_ind <- (i-1)*block_size + 1
    end_ind <- beg_ind + block_size - 1

    block_start <- (blocks[i] - 1) * block_size + 1
    block_end <- block_start + block_size - 1
    new_seq[beg_ind:end_ind] <- mc_seq[block_start:block_end]
  }
  return(new_seq)
}

overlap_bootstrap <- function(mc_seq, block_size){
  seq_len <- length(mc_seq)
  n_blocks <- floor(seq_len / block_size)
  new_seq <- rep(NA, length = n_blocks * block_size)
  block_loc <- sample(1:(seq_len - block_size + 1), replace = T)
  for(i in 1:n_blocks){
    beg_ind <- (i-1)*block_size + 1
    end_ind <- beg_ind + block_size - 1
    block_start <- block_loc[i]
    block_end <- block_start + block_size - 1
    new_seq[beg_ind:end_ind] <- mc_seq[block_start:block_end]
  }
  return(new_seq)
}

stationary_bootstrap <- function(mc_seq, p){
  seq_len <- length(mc_seq)
  block_sizes <- rgeom(seq_len, p) + 1
  block_locs <- cumsum(block_sizes)
  n_blocks <- sum(block_locs <= seq_len) + 1
  block_starts <- block_locs - block_sizes + 1
  block_ind <- sample(1:seq_len, size = n_blocks, replace = T)
  new_seq <- rep(NA, length = cumsum(block_sizes)[n_blocks])
  for(i in 1:n_blocks){
    block_size <- block_sizes[i]
    # beg_ind <- (i-1)*block_size + 1
    beg_ind <- block_starts[i]
    end_ind <- beg_ind + block_size - 1
    block_start <- block_ind[i]
    block_end <- block_start + block_size - 1
    seq_pos <- ((block_start:block_end - 1) %% seq_len) + 1
    new_seq[beg_ind:end_ind] <- mc_seq[seq_pos]
  }
  new_seq[1:seq_len]
}

#-------------------------------------------------------------------------------

# Simulation Cases -------------------------------------------------------------

fo_sim_cases <- list('fo_low' = matrix(c(0.01,0.99,0.98,0.02),2,2,T),
                     'fo_med' = matrix(c(0.2,0.8,0.9,0.1),2,2,T),
                     'fo_hi' = matrix(c(0.45,0.55,0.55,0.45),2,2,T))



# Empirical Distribution of Entropy Rate ---------------------------------------

n_sims <- 1000
mc_lens <- c(100,250,500,1000,2500)
emp_results <- list()
for(sc in 1:length(fo_sim_cases)){
  tm <- fo_sim_cases[[sc]]
  true_er <- CalcMarkovEntropyRate(tm, CalcEigenStationary(tm))
  entropy_results <- matrix(NA, nrow = n_sims, ncol = length(mc_lens))
  for(j in 1:length(mc_lens)){
    for(i in 1:n_sims){
      mc <- SimulateMarkovChain(tm, mc_lens[j])
      entropy_results[i,j] <- efficient_mc_er(mc, 2)
    }
    message('.', appendLF = F)
  }
  emp_results[[sc]] <- entropy_results
}

pdf('empirical_sampling_dist.pdf', height=4, width=12)
par(mfrow=c(1,3))
for(sc in 1:length(fo_sim_cases)){
  boxplot.matrix(emp_results[[sc]], horizontal=T, axes = F, pch=19,
                 ylim=range(emp_results[[sc]], na.rm = T))
  axis(1, at=seq(0,1,0.05))
  axis(2, 1:5, mc_lens, las=2)
  abline(v=true_er, lty=3, col='red')
}
title('\n\nEmpirical Sampling Distribution of Entropy Rate - First Order MC', outer=T)
dev.off()

# ------------------------------------------------------------------------------

# Which LZ77 Method? ---------

tm <- fo_sim_cases[[2]]
true_er <- CalcMarkovEntropyRate(tm, CalcEigenStationary(tm))
seq_len = 1000
window_sizes <- seq(50,950,50)
n_sims <- 500
exp_res <- matrix(NA, nrow=n_sims, ncol=1)
fw_res <- matrix(NA, nrow=n_sims, ncol=length(window_sizes))
fwq_res <- matrix(NA, nrow=n_sims, ncol=length(window_sizes))
for(i in 1:n_sims){
  mc <- SimulateMarkovChain(tm, seq_len)
  exp_res[i, ] <- SWLZEntRate(mc, F)
  for(ws in 1:length(window_sizes)){
    fw_res[i, ws] <- fixedwindow_lz77(mc, window_sizes[ws])
    fwq_res[i, ws] <- fixedquick_lz77(mc, window_sizes[ws])
  }
  message('.', appendLF = 'F')
}


pdf('comparing_lz_methods.pdf', height=4, width=10)
par(mfrow=c(1,3))
plot(window_sizes, rep(quantile(exp_res)[3], length(window_sizes)),
     type='b', pch=19, ylim=c(0,1),
     main='Expanding Window LZ77')
lines(window_sizes, rep(quantile(exp_res)[1], length(window_sizes)), lty=3)
lines(window_sizes, rep(quantile(exp_res)[2], length(window_sizes)), lty=3)
lines(window_sizes, rep(quantile(exp_res)[4], length(window_sizes)), lty=3)
lines(window_sizes, rep(quantile(exp_res)[5], length(window_sizes)), lty=3)
abline(h=quantile(emp_results[[2]][,4])[1:5], lwd=2, col=rgb(0,0,0,0.5))
legend('bottomleft', c('IQR - EmpDist', 'LZ77 - Expanding'),
       lty=c(1,3))

plot(window_sizes, apply(fw_res, 2, quantile)[3,], ylim=c(0,1), type='b', pch=19,
     main='Fixed Width Window LZ77\nEvaluated at all Points')
lines(window_sizes, apply(fw_res, 2, quantile)[1,], lty=3)
lines(window_sizes, apply(fw_res, 2, quantile)[2,], lty=3)
lines(window_sizes, apply(fw_res, 2, quantile)[4,], lty=3)
lines(window_sizes, apply(fw_res, 2, quantile)[5,], lty=3)
abline(h=quantile(emp_results[[2]][,4])[1:5], lwd=2, col=rgb(0,0,0,0.5))
abline(h=true_er, col=rgb(1,0,0,0.75))
legend('bottomleft', c('IQR - EmpDist', 'LZ77 - Fixed Width'),
       lty=c(1,3))

plot(window_sizes, apply(fwq_res, 2, quantile)[3,], ylim=c(0,1), type='b', pch=19,
     main='Fixed Width Window LZ77\nLinear Runtime')
lines(window_sizes, apply(fwq_res, 2, quantile)[1,], lty=3)
lines(window_sizes, apply(fwq_res, 2, quantile)[2,], lty=3)
lines(window_sizes, apply(fwq_res, 2, quantile)[4,], lty=3)
lines(window_sizes, apply(fwq_res, 2, quantile)[5,], lty=3)
abline(h=quantile(emp_results[[2]][,4])[1:5], lwd=2, col=rgb(0,0,0,0.5))
abline(h=true_er, col=rgb(1,0,0,0.75))
legend('bottomright', c('IQR - EmpDist', 'LZ77 - Fixed Width Quick'),
       lty=c(1,3))
dev.off()


# Bootstrap simulations --------------------------------------------------------
tm <- fo_sim_cases[[2]]
true_er <- CalcMarkovEntropyRate(tm, CalcEigenStationary(tm))
seq_len = 1000
bootsamps <- 100
emp_sd <- sd(emp_results[[2]][,4])
n_sims <- 100
block_sizes <- c(10, 25, 50, 75, 100, 250)
sd1_res <- matrix(NA, nrow = n_sims, ncol=length(block_sizes))
sd2_res <- matrix(NA, nrow = n_sims, ncol=length(block_sizes))
for(i in 1:n_sims){
  mc <- SimulateMarkovChain(tm, seq_len)
  for(bs in 1:length(block_sizes)){
    block_size <- block_sizes[bs]
    bs1_res <- rep(NA, length(bootsamps))
    bs2_res <- rep(NA, length(bootsamps))
    for(B in 1:bootsamps){
      ns <- block_bootstrap(mc, block_size)
      bs1_res[B] <- SWLZEntRate(ns)
      bs2_res[B] <- efficient_mc_er(ns)
    }
    sd1_res[i, bs] <- sd(bs1_res)
    sd2_res[i, bs] <- sd(bs2_res)
  }
  message('.', appendLF = F)
}

pdf('simulation_stddev_block1.pdf', height=4, width=10)
par(mfrow=c(1,2), mar=c(4,4,6,4)+0.1)
boxplot.matrix(sd1_res, horizontal=T, axes = F, pch=19,
               ylim=range(sd1_res, na.rm = T))
axis(1, at=seq(0,1,0.05))
axis(2, 1:6, block_sizes, las=2)
abline(v=emp_sd, lty=3, col='red')
title('Expanding Window LZ77')
boxplot.matrix(sd2_res, horizontal=T, axes = F, pch=19,
               ylim=range(sd2_res, na.rm = T))
axis(1, at=seq(0,1,0.05))
axis(2, 1:6, block_sizes, las=2)
abline(v=emp_sd, lty=3, col='red')
title('Markov Chain')
title('\n\nNon-Overlapping Block Bootstrap', outer = T)
dev.off()

# Overlapping Block Bootstrap --------------------------------------------------

tm <- fo_sim_cases[[2]]
true_er <- CalcMarkovEntropyRate(tm, CalcEigenStationary(tm))
seq_len = 1000
bootsamps <- 100
emp_sd <- sd(emp_results[[2]][,4])
n_sims <- 100
block_sizes <- c(10, 25, 50, 75, 100, 250)
sd1_res <- matrix(NA, nrow = n_sims, ncol=length(block_sizes))
sd2_res <- matrix(NA, nrow = n_sims, ncol=length(block_sizes))
for(i in 1:n_sims){
  mc <- SimulateMarkovChain(tm, seq_len)
  for(bs in 1:length(block_sizes)){
    block_size <- block_sizes[bs]
    bs1_res <- rep(NA, length(bootsamps))
    bs2_res <- rep(NA, length(bootsamps))
    for(B in 1:bootsamps){
      ns <- overlap_bootstrap(mc, block_size)
      bs1_res[B] <- SWLZEntRate(ns)
      bs2_res[B] <- efficient_mc_er(ns)
    }
    sd1_res[i, bs] <- sd(bs1_res)
    sd2_res[i, bs] <- sd(bs2_res)
  }
  message('.', appendLF = F)
}

pdf('simulation_stddev_block2.pdf', height=4, width=10)
par(mfrow=c(1,2), mar=c(4,4,6,4)+0.1)
boxplot.matrix(sd1_res, horizontal=T, axes = F, pch=19,
               ylim=range(sd1_res, na.rm = T))
axis(1, at=seq(0,1,0.05))
axis(2, 1:6, block_sizes, las=2)
abline(v=emp_sd, lty=3, col='red')
title('Expanding Window LZ77')
boxplot.matrix(sd2_res, horizontal=T, axes = F, pch=19,
               ylim=range(sd2_res, na.rm = T))
axis(1, at=seq(0,1,0.05))
axis(2, 1:6, block_sizes, las=2)
abline(v=emp_sd, lty=3, col='red')
title('Markov Chain')
title('\n\nOverlapping Block Bootstrap', outer = T)
dev.off()


# Stationary Block Bootstrap --------------------------------------------------

tm <- fo_sim_cases[[2]]
true_er <- CalcMarkovEntropyRate(tm, CalcEigenStationary(tm))
seq_len = 1000
bootsamps <- 100
emp_sd <- sd(emp_results[[2]][,4])
n_sims <- 100
block_sizes <- c(0.005, 0.01, 0.05, 0.1, 0.2, 0.5)
sd1_res <- matrix(NA, nrow = n_sims, ncol=length(block_sizes))
sd2_res <- matrix(NA, nrow = n_sims, ncol=length(block_sizes))
pdf('stat_bootstrap_dist.pdf', height=8, width=14)
par(mfrow=c(2,3))
for(i in 1:n_sims){
  mc <- SimulateMarkovChain(tm, seq_len)
  for(bs in 1:length(block_sizes)){
    block_size <- block_sizes[bs]
    bs1_res <- rep(NA, length(bootsamps))
    bs2_res <- rep(NA, length(bootsamps))
    for(B in 1:bootsamps){
      ns <- stationary_bootstrap(mc, block_size)
      bs1_res[B] <- SWLZEntRate(ns)
      bs2_res[B] <- efficient_mc_er(ns)
    }
    sd1_res[i, bs] <- sd(bs1_res)
    sd2_res[i, bs] <- sd(bs2_res)
    if(i == 1){
      ymax <- max(max(density(emp_results[[2]])$y),
                  max(density(bs1_res)$y),
                  max(density(bs2_res)$y))
      plot(density(emp_results[[2]][,4]), xlim=c(0,1), ylim=c(0,ymax),
           xlab='', main=paste('Effective Block Size =',1/block_size))
      lines(density(bs1_res), col='red')
      lines(density(bs2_res), col='blue')
      legend('topleft',
             c('EmpDist', 'SWLZ', 'MC'),
             lty=1,
             col=c('black', 'red', 'blue'))
    }
  }
  message('.', appendLF = F)
}
dev.off()

pdf('simulation_stddev_block3.pdf', height=4, width=10)
par(mfrow=c(1,2), mar=c(4,4,6,4)+0.1)
boxplot.matrix(sd1_res, horizontal=T, axes = F, pch=19,
               ylim=range(sd1_res, na.rm = T))
axis(1, at=seq(0,1,0.05))
axis(2, 1:6, block_sizes, las=2)
abline(v=emp_sd, lty=3, col='red')
title('Expanding Window LZ77')
boxplot.matrix(sd2_res, horizontal=T, axes = F, pch=19,
               ylim=range(sd2_res, na.rm = T))
axis(1, at=seq(0,1,0.05))
axis(2, 1:6, block_sizes, las=2)
abline(v=emp_sd, lty=3, col='red')
title('Markov Chain')
title('\n\nStationary Block Bootstrap', outer = T)
dev.off()


tm <- fo_sim_cases[[2]]
true_er <- CalcMarkovEntropyRate(tm, CalcEigenStationary(tm))
seq_len = 1000
bootsamps <- 500
test_p <- true_er / log2(seq_len)
n_sims <- 100
bs1_res <- matrix(NA, nrow=n_sims, ncol = bootsamps)
bs2_res <- matrix(NA, nrow=n_sims, ncol = bootsamps)
for(i in 1:n_sims){
  mc <- SimulateMarkovChain(tm, seq_len)
  for(bs in 1:bootsamps){
    new_mc <- stationary_bootstrap(mc, test_p)
    lzent <- SWLZEntRate(new_mc)
    mcent <- efficient_mc_er(new_mc)
    bs1_res[i, bs] <- lzent
    bs2_res[i, bs] <- mcent
  }
  message('.', appendLF = F)
}

par(mfrow=c(1,2))
plot(density(emp_results[[2]][,4]), xlim=c(0,1), ylim=c(0,ymax),
     xlab='', main=paste('Effective Block Size =',1/test_p))
for(s in 1:nrow(bs1_res)){
  lines(density(bs1_res[s,]), col=rgb(0.75,0,0,0.25))
}



plot(density(emp_results[[2]][,4]), xlim=c(0,1), ylim=c(0,ymax),
     xlab='', main=paste('Effective Block Size =',1/test_p))
for(s in 1:nrow(bs1_res)){
  lines(density(bs2_res[s,]), col=rgb(0,0,0.75,0.25))
}





# set.seed(1278)
# mc <- SimulateMarkovChain(tm, 1000)
# bs_samps <- 500
# block_sizes <- seq(10, 100, 5)
# bs_results <- matrix(NA, nrow=bs_samps, ncol=length(block_sizes))
# bs_results2 <- matrix(NA, nrow=bs_samps, ncol=length(block_sizes))
# bs_results3 <- matrix(NA, nrow=bs_samps, ncol=length(block_sizes))
# bs_results4 <- matrix(NA, nrow=bs_samps, ncol=length(block_sizes))
# for(bs in 1:length(block_sizes)){
#   for(B in 1:bs_samps){
#     ns <- block_bootstrap(mc_seq = mc, block_size = block_sizes[bs])
#     ns2 <- overlap_bootstrap(mc_seq = mc, block_size = block_sizes[bs])
#     bs_results[B, bs] <- efficient_mc_er(ns, mc_order=1)
#     bs_results2[B, bs] <- SWLZEntRate(ns)
#     bs_results3[B, bs] <- efficient_mc_er(ns2, mc_order=1)
#     bs_results4[B, bs] <- SWLZEntRate(ns2)
#   }
# }
# par(mfrow=c(1,2))
# plot(block_sizes, apply(bs_results, 2, sd), type='l', ylim=c(0,0.25))
# lines(block_sizes, apply(bs_results2, 2, sd), type='l', col='blue')
# lines(block_sizes, apply(bs_results3, 2, sd), type='l', col='green')
# lines(block_sizes, apply(bs_results4, 2, sd), type='l', col='red')
# abline(h=apply(entropy_results,2,sd)[4])
# plot(block_sizes, apply(bs_results, 2, mean), type='l', ylim=c(0,1))
# lines(block_sizes, apply(bs_results2, 2, mean), type='l', col='blue')
# lines(block_sizes, apply(bs_results3, 2, mean), type='l', col='green')
# lines(block_sizes, apply(bs_results4, 2, mean), type='l', col='red')
# abline(h=apply(entropy_results,2,mean)[4])



# Stationary Bootstrap ------ Probs not steps------------------------
# set.seed(1278)
# mc <- SimulateMarkovChain(tm, seq_len)
# bs_samps <- 500
# block_sizes <- seq(0.005,0.5,length.out = 25)
# bs_results <- matrix(NA, nrow=bs_samps, ncol=length(block_sizes))
# bs_results2 <- matrix(NA, nrow=bs_samps, ncol=length(block_sizes))
# for(bs in 1:length(block_sizes)){
#     for(B in 1:bs_samps){
#         ns <- stationary_bootstrap(mc, block_sizes[bs])
#         bs_results[B, bs] <- efficient_mc_er(ns, mc_order=1)
#         bs_results2[B, bs] <- SWLZEntRate(ns)
#     }
# }
# par(mfrow=c(1,2))
# plot(block_sizes, apply(bs_results, 2, sd), type='l', ylim=c(0,0.25))
# lines(block_sizes, apply(bs_results2, 2, sd), type='l', col='blue')
# abline(h=apply(entropy_results,2,sd)[4])
# abline(v=true_er/log2(1000))
# plot(block_sizes, apply(bs_results, 2, mean), type='l', ylim=c(0,1))
# lines(block_sizes, apply(bs_results2, 2, mean), type='l', col='blue')
# abline(h=apply(entropy_results,2,mean)[4])
# abline(v=true_er/log2(1000))
#
#
#
#


#
# for(i in 1:500){
#   uh <- efficient_lz77(mc, i)
#   print(c(uh, CalcMarkovEntropyRate(tm, CalcEigenStationary(tm))))
# }
#
#
# mc <- SimulateMarkovChain(tm, 500)
# win_sizes <- seq(10,250, 10)
# win_res <- matrix(NA, nrow=length(win_sizes), ncol=2)
# for(i in 1:length(win_sizes)){
#   win_res[i,] <- c(win_sizes[i],
#                    efficient_lz77(mc, win_sizes[i]))
# }
# plot(win_res[,1], win_res[,2])
# abline(h=CalcMarkovEntropyRate(tm, CalcEigenStationary(tm)), lty=3, col='red')
# abline(h=SWLZEntRate(mc))
#
# n_sims <- 100
# sim_res <- matrix(NA, nrow=n_sims, ncol=3)
# for(s in 1:n_sims){
#   mc <- SimulateMarkovChain(tm, 5000)
#   sim_res[s, ] <- c(efficient_mc_er(mc),
#                     efficient_lz77(mc, 2500),
#                     SWLZEntRate(mc))
# }
# boxplot.matrix(sim_res)
# abline(h=CalcMarkovEntropyRate(tm, CalcEigenStationary(tm)), lty=3, col='red')
#
#
