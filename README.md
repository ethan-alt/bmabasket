
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bmabasket

<!-- badges: start -->
[![R-CMD-check](https://github.com/ethan-alt/bmabasket/workflows/R-CMD-check/badge.svg)](https://github.com/ethan-alt/bmabasket/actions)
<!-- badges: end -->

The goal of bmabasket is to simulate basket trial data based on
hyperparameters and analyze things such as the family-wise error rate,
bias, and MSE. The package uses Bayesian model average (BMA) to compute
the posterior probability that response within a basket exceeds some
threshold.

## Installation

You can install the released version of bmabasket from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("bmabasket")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ethan-alt/bmabasket")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(bmabasket)
## REPEAT SIMS FROM BIOSTATISTICS JOURNAL PUBLICATION
nSims      <- 100             ## change to ~250000 to repeat journal results                   
meanTime   <- 0.01
sdTime     <- 0.0000000001
mu0        <- 0.45
phi0       <- 1.00
ppEffCrit  <- 0.985
ppFutCrit  <- 0.2750
pmp0       <- 2
n1         <- 7
n2         <- 16
targSSPer  <- c(n1, n2)
nInterim   <- 2
futOnly    <- 1
K0         <- 5
row        <- 0
mss        <- 4
minSSFut   <- mss  ## minimum number of subjects in basket to assess futility using BMA
minSSEff   <- mss  ## minimum number of subjects in basket to assess activity using BMA
rTarg      <- 0.45
rNull      <- 0.15
rRatesMod  <- matrix(rNull,(K0+1)+3,K0)
rRatesNull <- rep(rNull,K0)
rRatesMid  <- rep(rTarg,K0)
eRatesMod  <- rep(1, K0)

## min and max #' of new subjects per basket before next analysis (each row is interim)
minSSEnr <- matrix(rep(mss, K0), nrow=nInterim ,ncol=K0, byrow=TRUE) 
maxSSEnr <- matrix(rep(100, K0), nrow=nInterim, ncol=K0, byrow=TRUE) 

## construct matrix of rates
for (i in 1:K0)  
{
  rRatesMod[(i+1):(K0+1),i]= rTarg     
}
rRatesMod[(K0+2),] <- c(0.05,0.15,0.25,0.35,0.45)
rRatesMod[(K0+3),] <- c(0.15,0.30,0.30,0.30,0.45)
rRatesMod[(K0+4),] <- c(0.15,0.15,0.30,0.30,0.30)

## conduct simulation of trial data and analysis
x <- bma_design(
  nSims, K0, K0, eRatesMod, rRatesMod[i+1,], meanTime, sdTime, 
  ppEffCrit, ppFutCrit, as.logical(futOnly), rRatesNull, rRatesMid, 
  minSSFut, minSSEff, minSSEnr, maxSSEnr, targSSPer, nInterim, mu0, 
  phi0, priorModelProbs = NULL, pmp0 = pmp0
)
x
#> $hypothesis.testing
#> $hypothesis.testing$rr
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,] 0.86 0.94 0.88 0.86 0.88
#> 
#> $hypothesis.testing$fw.fpr
#> [1] 0
#> 
#> $hypothesis.testing$nerr
#> [1] 0
#> 
#> $hypothesis.testing$fut
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,] 0.06 0.03 0.03 0.06 0.08
#> 
#> 
#> $sample.size
#> $sample.size$basket.ave
#>       [,1]  [,2]  [,3]  [,4]  [,5]
#> [1,] 21.91 22.31 21.37 21.08 21.53
#> 
#> $sample.size$basket.med
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]   23 23.5   23   22   23
#> 
#> $sample.size$basket.min
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    4    4    4    6    4
#> 
#> $sample.size$basket.max
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]   33   32   33   33   34
#> 
#> $sample.size$overall.ave
#>       [,1]
#> [1,] 108.2
#> 
#> 
#> $point.estimation
#> $point.estimation$PM.ave
#>           [,1]   [,2]    [,3]      [,4]      [,5]
#> [1,] 0.4398896 0.4589 0.46378 0.4436683 0.4331298
#> 
#> $point.estimation$SP.ave
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 0.4334567 0.4601979 0.4698621 0.4414528 0.4266855
#> 
#> $point.estimation$PP.ave
#>           [,1]      [,2]      [,3]     [,4]      [,5]
#> [1,] 0.9603756 0.9787252 0.9701788 0.955995 0.9449077
#> 
#> $point.estimation$bias
#>            [,1]        [,2]       [,3]         [,4]       [,5]
#> [1,] -0.0101104 0.008900043 0.01378003 -0.006331726 -0.0168702
#> 
#> $point.estimation$mse
#>            [,1]        [,2]       [,3]       [,4]       [,5]
#> [1,] 0.01283483 0.009855582 0.01291377 0.01489975 0.01278448
#> 
#> 
#> $trial.duration
#> $trial.duration$average
#> [1] 58.41536
#> 
#> 
#> $early.stopping
#> $early.stopping$interim.stop.prob
#>      [,1] [,2]
#> [1,] 0.06 0.94
#> 
#> $early.stopping$baskets.continuing.ave
#>      [,1] [,2]
#> [1,] 4.46 0.32
```
