source('~/git/entropyRate/entropyRate.R')

tm <- matrix(c(0.2,0.8,0.9,0.1),2,2,T)
CalcEigenStationary(tm)
mc <- SimulateMarkovChain(tm, 1000)
CalcEmpiricalStationary(mc, 1:2)


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

set.seed(1278)
mc <- SimulateMarkovChain(tm, 1000)
bs_samps <- 100
block_sizes <- seq(10, 100, 5)
bs_results <- matrix(NA, nrow=bs_samps, ncol=length(block_sizes))
bs_results2 <- matrix(NA, nrow=bs_samps, ncol=length(block_sizes))
for(bs in 1:length(block_sizes)){
  for(B in 1:bs_samps){
    ns <- block_bootstrap(mc_seq = mc, block_size = block_sizes[bs])
    bs_results[B, bs] <- efficient_mc_er(ns, mc_order=1)
    bs_results2[B, bs] <- SWLZEntRate(ns)
  }
}
print()
plot(block_sizes, apply(bs_results, 2, sd), type='l', ylim=c(0,0.1))
lines(block_sizes, apply(bs_results2, 2, sd), type='l', col='blue')
abline(h=apply(entropy_results,2,sd)[4])

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
                    efficient_lz77(mc, 1000),
                    SWLZEntRate(mc))
}
boxplot.matrix(sim_res)
abline(h=CalcMarkovEntropyRate(tm, CalcEigenStationary(tm)), lty=3, col='red')


