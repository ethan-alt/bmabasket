## Simulate data with 3 baskets
probs <- c(0.5, 0.5, 0.5)
n <- rep(350, length(probs))
y <- rbinom(length(probs), size = n, prob = probs)

datMat <- cbind(y, n-y)
pi0 <- probs
mu0 = 0.5
phi0 = 2
a0 <- mu0 * phi0
b0 <- (1 - mu0) * phi0
lbeta_a0b0 <- lbeta(a0, b0)


## Obtain posterior model probabilities
postModelProbs <- numeric(ncol(parts))
for ( j in 1:ncol(parts) ) {
  postModelProbs[j] <- logPostSurvProb(pi0, datMat, parts[, j], a0, b0, lbeta_a0b0)$sumLogBeta
}
postModelProbs <- round( exp(postModelProbs) / sum(exp(postModelProbs)), 3)
postModelProbs
# 
# 
# postProbs <- rep(0, length(probs))
# for ( j in 1:ncol(parts) ) {
#   postProbs <- postProbs + postModelProbs[j] * exp( logPostSurvProb(pi0, datMat, parts[, j], a0, b0, lbeta_a0b0)$logPostSurvProb )
# }
# postProbs



bma_cpp(pi0, datMat, parts, mu0 = mu0, phi0 = phi0, logModelPriors)

bma(.5, y, n)

