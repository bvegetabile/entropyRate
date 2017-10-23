fotm <- function(p,q){matrix(c(1-p, p, q, 1-q), 2,2, T)}
# bigtm <- function(){CalcTransitionMatrix(matrix(abs(rgeom(64,0.4)), 8,8,T))}
# rgeom(64,0.7)
# hmm <- rep(NA, 1000)
# for(i in 1:1000){
#     tm <- bigtm()
#     hmm[i] <- CalcMarkovEntropyRate(tm, CalcEigenStationary(tm))
# }
# hist(hmm)


n_sims <- 1000
samp_sizes <- seq(100,1000,50)
results <- matrix(NA, nrow = length(samp_sizes), ncol = n_sims)
biases <- matrix(NA, nrow = length(samp_sizes), ncol = n_sims)
ents <- matrix(NA, nrow = length(samp_sizes), ncol = n_sims)
sizes <- matrix(NA, nrow = length(samp_sizes), ncol = n_sims)
for(ss in 1:length(samp_sizes)){
    for(i in 1:n_sims){
        pee <- runif(1,0,1)
        cue <- runif(1,0,1)
        tm <- fotm(pee, cue)
        ent <- CalcMarkovEntropyRate(tm, CalcEigenStationary(tm))
        mc <- SimulateMarkovChain(tm, samp_sizes[ss])
        est <- SWLZEntRate(mc)
        results[ss, i] <- est
        biases[ss, i] <- (ent - est)
        ents[ss, i] <- ent
        sizes[ss, i] <- samp_sizes[ss]
    }
    message('.', appendLF = F)
}
spot <- 19
samp_sizes[spot]
plot(results[spot,], biases[spot,], 
     xlim=c(0,1), pch=19, col=rgb(0,0,0,0.25),
     xlab='Estimated Entropy Rate',
     ylab='Bias')
abline(h=0, lwd=3, col=rgb(1,0,0,0.5))
smoothed <- loess.smooth(results[spot,], biases[spot,], span = 0.2)
lines(smoothed$x, smoothed$y, lwd=3, col=rgb(0,0,1,0.5))
# abline(h=mean(biases[18,]), lwd=3, col=rgb(0,0,1,0.5))

res <- matrix(NA, nrow=1000, ncol=3)
tm <- matrix(c(0.45,0.55,0.55,0.45),2,2,T)
for(i in 1:1000){
    ent <- CalcMarkovEntropyRate(tm, CalcEigenStationary(tm))
    mc <- SimulateMarkovChain(tm, 500)
    est <- SWLZEntRate(mc)
    res[i, 1] <- ent
    res[i, 2] <- est
    res[i, 3] <- ent - est
}
hist(res[,2])
