rm(list = ls())
source("./R/SFET_extend.R")
source("./R/fast_SFET.R")
source("./R/SFET_tidy.R")
library(dplyr)
## scenario 1: small number of big strata, rare cases, basically no effect
n <- c(1200, 1210, 1100)
m <- c(670, 500, 490)
z <- c(43, 35, 44)
x <- c(21, 19, 20)

strata <- rep(1:3, each = 2)
treatment <- rep(1:2, 3)
nsubjects <- c(670, 530, 500, 710, 490, 610)
nresponders <- c(21, 22, 19, 16, 20, 24)
case <- 1

time11 <- system.time(pvalue11 <- SFET_extend2(strata = strata, treatment = treatment, nsubjects = nsubjects, nresponders = nresponders,
                                     case = case, side = "up"))
# p-value: 0.3947835
# user  system elapsed
# 0.55    0.00    0.55
time12 <- system.time(pvalue12 <- fast_SFET(strata = strata, treatment = treatment, nsubjects = nsubjects, nresponders = nresponders,
                                       case = case, side = "up"))

# p-value: 0.3947835
# user  system elapsed
# 0.19    0.00    0.18

## scenario 2: big number of small strata, rare cases, basiccaly no effect
n2 <- c(120, 121, 110,100, 132, 119, 106, 132, 140, 100)
m2 <- c(59, 69, 57, 66, 56, 66, 55, 59, 66, 63)
z2 <- c(4, 4, 4, 5, 3, 4, 4, 5, 4, 4)
x2 <- c(2, 1, 2, 2, 2, 2, 2, 2, 2, 2)

strata <- rep(1:10, each = 2)
treatment <- rep(1:2, 10)
nsubjects <- c(59, 61, 69, 52, 57, 53, 66, 34, 56, 76, 66, 53, 55, 51, 59, 73, 66, 74, 63, 37)
nresponders <- c(2, 2, 1, 3, 2, 2, 2, 3, 2, 1, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2)
case <- 1

time21 <- system.time(pvalue21 <- SFET_extend2(strata = strata2, treatment = treatment2, nsubjects = nsubjects2, nresponders = nresponders2,
                                     case = case2, side = "up"))
## p-value: 0.8573128
## user  system elapsed
## 132.63    0.76  135.52
time22 <- system.time(pvalue22 <- fast_SFET(strata = strata2, treatment = treatment2, nsubjects = nsubjects2, nresponders = nresponders2,
                                               case = case2, prec = 1e-3, side = "up"))
## p-value: 0.8573128
## user  system elapsed
## 136.61    0.02  136.79

## no sample points were reduced, thus take longer time. This may happen in balanced, rare cases

## scenario 3: medium number of strata
strata <- rep(1:3, each = 2)
treatment <- rep(1:2, 3)
nsubjects <- c(670, 530, 500, 710, 490, 610)
nresponders <- c(210, 170, 190, 167, 200, 240)
case <- 1

time31 <- system.time(pvalue31 <- SFET_extend2(strata = strata, treatment = treatment, nsubjects = nsubjects, nresponders = nresponders,
                                               case = case, side = "up"))
# p-value: 0.0007500332
# user  system elapsed
# 362.89    0.12  363.61
time32 <- system.time(pvalue32 <- fast_SFET(strata = strata, treatment = treatment, nsubjects = nsubjects, nresponders = nresponders,
                                            case = case, prec = 1e-5, side = "up"))
# p-value: 0.0007500156
# user  system elapsed
# 3.63    0.02    3.69

set.seed(1127)
test33 <- system.time(pvalue33 <- sfet_naive_mc(strata = strata, treatment = treatment, nsubjects = nsubjects, nresponders = nresponders,
                                          case = case, sample_size = 500000, side = "up"))

# p-value: 0.000702
# user  system elapsed
# 233.51    0.07  234.11
