
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bmabasket

<!-- badges: start -->

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
#> [1,] 0.84 0.89 0.79 0.88 0.89
#> 
#> $hypothesis.testing$fw.fpr
#> [1] 0
#> 
#> $hypothesis.testing$nerr
#> [1] 0
#> 
#> $hypothesis.testing$fut
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,] 0.08 0.02 0.13 0.01 0.04
#> 
#> 
#> $sample.size
#> $sample.size$basket.ave
#>       [,1] [,2]  [,3]  [,4]  [,5]
#> [1,] 21.55 22.1 21.11 22.77 22.38
#> 
#> $sample.size$basket.med
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]   23   23   23   23   23
#> 
#> $sample.size$basket.min
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    5    5    4    5    4
#> 
#> $sample.size$basket.max
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]   34   34   34   34   34
#> 
#> $sample.size$overall.ave
#>        [,1]
#> [1,] 109.91
#> 
#> 
#> $point.estimation
#> $point.estimation$PM.ave
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 0.4264962 0.4448821 0.4146264 0.4470727 0.4450968
#> 
#> $point.estimation$SP.ave
#>           [,1]      [,2]      [,3]      [,4]     [,5]
#> [1,] 0.4218376 0.4441119 0.4068842 0.4458202 0.443398
#> 
#> $point.estimation$PP.ave
#>           [,1]      [,2]      [,3]     [,4]      [,5]
#> [1,] 0.9372382 0.9839219 0.9095636 0.987195 0.9706107
#> 
#> $point.estimation$bias
#>             [,1]         [,2]       [,3]         [,4]         [,5]
#> [1,] -0.02350383 -0.005117902 -0.0353736 -0.002927319 -0.004903175
#> 
#> $point.estimation$mse
#>            [,1]        [,2]       [,3]        [,4]      [,5]
#> [1,] 0.01630138 0.007935678 0.01983467 0.007982926 0.0122709
#> 
#> 
#> $trial.duration
#> $trial.duration$average
#> [1] 66.40398
#> 
#> 
#> $early.stopping
#> $early.stopping$interim.stop.prob
#>      [,1] [,2]
#> [1,] 0.03 0.97
#> 
#> $early.stopping$baskets.continuing.ave
#>      [,1] [,2]
#> [1,] 4.57 0.43
```
