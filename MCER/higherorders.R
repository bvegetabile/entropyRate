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

seq_len <- 5000
mc2 <- SimulateTwoOrder(tm2, seq_len)
fixed_window_ests <- matrix(NA, nrow=seq_len, ncol=2)
for(i in 1:seq_len){
    fixed_window_ests[i,] <- c(i, fixedquick_lz77(mc2, i))
}
plot(fixed_window_ests[,1], fixed_window_ests[,2], type='l')
abline(h=SWLZEntRate(mc2))
abline(h=ent2, col='red', lty=3)



SWLZEntRate(mc2)

mc_ents <- matrix(NA, nrow=15, ncol=2)
for(i in 1:15){
    mc_ents[i,] <- c(i, efficient_mc_er(mc2, i))
}
plot(mc_ents[,1], mc_ents[,2])
abline(h=SWLZEntRate(mc2))
abline(h=ent2, col='red', lty=3)
