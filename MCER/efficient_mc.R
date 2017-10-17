source('~/git/entropyRate/entropyRate.R')

tm <- matrix(c(0.2,0.8,0.9,0.1),2,2,T)
# tm <- matrix(c(0.49,0.51,0.51,0.49),2,2,T)
# tm <- matrix(c(0.01,0.99,0.98,0.02),2,2,T)
CalcEigenStationary(tm)
mc <- SimulateMarkovChain(tm, 1000)
CalcEmpiricalStationary(mc, 1:2)
true_er <- CalcMarkovEntropyRate(tm, CalcEmpiricalStationary(mc, 1:2))

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

efficient_lz77 <- function(mc_seq, window_size = 10){
  char_seq <- paste(mc_seq, collapse = '')
  lz77entropy(char_seq, window_size)
}


for(o in 1:10){
  print(efficient_mc_er(mc, mc_order=o))
}

n_sims <- 1000
mc_lens <- c(100,250,500,1000,2500)
entropy_results <- matrix(NA, nrow = n_sims, ncol = length(mc_lens))
for(j in 1:length(mc_lens)){
  for(i in 1:n_sims){
    mc <- SimulateMarkovChain(tm, mc_lens[j])
    entropy_results[i,j] <- efficient_mc_er(mc, 2)
  }
  message('.', appendLF = F)
}
plot(mc_lens, sqrt(apply(entropy_results, 2, var)))
CalcMarkovEntropyRate(tm, CalcEigenStationary(tm))
boxplot.matrix(entropy_results, ylim=c(0,1))
abline(h=CalcMarkovEntropyRate(tm, CalcEigenStationary(tm)), lty=3, col='red')




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
set.seed(1278)
mc <- SimulateMarkovChain(tm, 1000)
bs_samps <- 500
block_sizes <- seq(0.005,0.5,length.out = 25)
bs_results <- matrix(NA, nrow=bs_samps, ncol=length(block_sizes))
bs_results2 <- matrix(NA, nrow=bs_samps, ncol=length(block_sizes))
for(bs in 1:length(block_sizes)){
    for(B in 1:bs_samps){
        ns <- stationary_bootstrap(mc, block_sizes[bs])
        bs_results[B, bs] <- efficient_mc_er(ns, mc_order=1)
        bs_results2[B, bs] <- SWLZEntRate(ns)
    }
}
par(mfrow=c(1,2))
plot(block_sizes, apply(bs_results, 2, sd), type='l', ylim=c(0,0.25))
lines(block_sizes, apply(bs_results2, 2, sd), type='l', col='blue')
abline(h=apply(entropy_results,2,sd)[4])
abline(v=true_er/log2(1000))
plot(block_sizes, apply(bs_results, 2, mean), type='l', ylim=c(0,1))
lines(block_sizes, apply(bs_results2, 2, mean), type='l', col='blue')
abline(h=apply(entropy_results,2,mean)[4])
abline(v=true_er/log2(1000))







for(i in 1:500){
  uh <- efficient_lz77(mc, i)
  print(c(uh, CalcMarkovEntropyRate(tm, CalcEigenStationary(tm))))
}


mc <- SimulateMarkovChain(tm, 500)
win_sizes <- seq(10,250, 10)
win_res <- matrix(NA, nrow=length(win_sizes), ncol=2)
for(i in 1:length(win_sizes)){
  win_res[i,] <- c(win_sizes[i],
                   efficient_lz77(mc, win_sizes[i]))
}
plot(win_res[,1], win_res[,2])
abline(h=CalcMarkovEntropyRate(tm, CalcEigenStationary(tm)), lty=3, col='red')
abline(h=SWLZEntRate(mc))

n_sims <- 100
sim_res <- matrix(NA, nrow=n_sims, ncol=3)
for(s in 1:n_sims){
  mc <- SimulateMarkovChain(tm, 5000)
  sim_res[s, ] <- c(efficient_mc_er(mc),
                    efficient_lz77(mc, 2500),
                    SWLZEntRate(mc))
}
boxplot.matrix(sim_res)
abline(h=CalcMarkovEntropyRate(tm, CalcEigenStationary(tm)), lty=3, col='red')


