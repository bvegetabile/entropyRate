


morder_transmat <- function(mc_order = 1, unique_states = 1:2){
  # mc_order <- 3
  # unique_states <- 1:8
  n_states <- length(unique_states)
  state_space <- data.frame(matrix(NA, nrow=n_states, ncol=mc_order))
  for(i in 1:mc_order){
    state_space[,i] <- unique_states
  }
  state_space <- c(expand.grid(state_space), sep="")
  state_space <- do.call(paste, state_space)
  trans_table <- matrix(NA, nrow=length(state_space), ncol=n_states)
  for(i in 1:length(state_space)){
    alphas <- rgamma(n_states, 500, 1000)
    p <- dirmult::rdirichlet(n=1, alphas)
    trans_table[i, ] <- as.vector(p)
  }
  
  trans_mat <- matrix(0, length(state_space), length(state_space))
  rownames(trans_mat) <- colnames(trans_mat) <- state_space
  for(i in 1:length(state_space)){
    for(j in 1:n_states){
      from_state <- state_space[i]
      to_state <- paste(state_space[i], unique_states[j], sep='')
      to_state <- substr(to_state, 2, 6)
      fs <- which(state_space == from_state)
      ts <- which(state_space == to_state)
      trans_mat[fs, ts] <- trans_table[i, j]
    }
  }
  list('trans_table' = trans_table,
       'trans_mat' = trans_mat,
       'true_er' = CalcMarkovEntropyRate(trans_mat, CalcEigenStationary(trans_mat)))
}

morder_seq <- function(mc_order, unique_states, trans_table, seq_len = 1000){
  n_states <- length(unique_states)
  state_space <- data.frame(matrix(NA, nrow=n_states, ncol=mc_order))
  for(i in 1:mc_order){
    state_space[,i] <- unique_states
  }
  state_space <- c(expand.grid(state_space), sep="")
  state_space <- do.call(paste, state_space)
  
  samp_seq <- rep(NA, seq_len)
  x0 <- sample(state_space, size = 1, replace = T)
  samp_seq[1:mc_order] <- as.integer(strsplit(x0, split = '')[[1]])
  start_spot <- mc_order + 1
  for(s in start_spot:seq_len){
    left_word <- paste(samp_seq[(s-mc_order):(s-1)], collapse='')
    tp <- trans_table[which(state_space == left_word), ]
    new_symbol <- sample(unique_states, size=1, replace = T, prob = tp)
    samp_seq[s] <- new_symbol
  }
  samp_seq
}

ent_ests <- rep(NA, 1000)
for(i in 1:1000){
  tm <- morder_transmat(1, 1:5)[[2]]
  ent_ests[i] <- CalcMarkovEntropyRate(tm, CalcEigenStationary(tm))
}
hist(ent_ests, breaks=seq(0, log2(5), length.out = 20))
abline(v=log2(5))

trans_table <- tm[[1]]


tm <-  morder_transmat(3, 1:5)
samp_seq <- morder_seq(3, 1:5, tm[[1]], seq_len = 10000)

for(i in 1:5){
  tc <- CalcTC_Mth_Order(samp_seq, 1:5, i)
  ugh <- CalcTransitionMatrix(tc)
  print(sum(tc * log(ugh), na.rm=T))
}

mco <- matrix(NA, nrow = 5, ncol=2)
for(o in 1:5){
  mco[o,] <- c(o, efficient_mc_er(samp_seq, o))
}
plot(mco, ylim=c(0,log2(5)))
abline(h=SWLZEntRate(samp_seq))
abline(h=tm[[3]], col='red')
abline(v=3, col='blue')

