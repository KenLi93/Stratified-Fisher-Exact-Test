#!/bin/Rscript
library(dplyr)
library(magrittr)

sfet_naive_mc <- function(strata, treatment, nsubjects, nresponders,
                          case, data = NULL, side = c("up"),
                          sample_size = 20000){
  nstrata <- strata %>% unique() %>% length()

  if (data %>% is.null()) {
    indata <- tibble(strata = strata, treatment = treatment, nsubjects = nsubjects,
                     nresponders = nresponders) %>% arrange(strata)
  } else {
    indata <- data %$% tibble(strata = strata, treatment = treatment, nsubjects = nsubjects,
                              nresponders = nresponders) %>% arrange(strata)
  }

  ss <- nresponders %>% subset(treatment == case) %>% sum()

  xx <- indata %>% filter(treatment == case) %>% .$nresponders

  zz <- indata %>%
    group_by(strata) %>%
    summarise(tot_responders = sum(nresponders)) %>%
    .$tot_responders

  mm <- indata %>%
    filter(treatment == case) %>%
    .$nsubjects

  nn <- indata %>%
    group_by(strata) %>%
    summarise(tot_subjects = sum(nsubjects)) %>%
    .$tot_subjects

  m_lower <- pmax(0, zz-nn+mm)
  m_upper <- pmin(zz, mm)

  ss_null <- replicate(sample_size, {
    (1:nstrata) %>%
      sapply(function(i){
        rhyper(1,
               m = zz[i],
               n = nn[i] - zz[i],
               k = mm[i])
      })
  }) %>% colSums()

  if (side == "up") {
    return(mean(ss_null > ss))
  } else if (side == "lower") {
    return(mean(ss_null < ss))
  } else if (side == "both") {
    return(min(1, 2 * min(mean(ss_null > ss), mean(ss_null < ss))))
  }
}

mod_top <- function(a, b) {
  res <- a %% b
  if (res == 0) {
    return(b)
  } else {
    return(res)
  }
}



SFET_extend <- function(strata, treatment, nsubjects, nresponders,
                        case, data = NULL, side = c("up")){

  nstrata <- strata %>% unique() %>% length()

  if (data %>% is.null()) {
    indata <- tibble(strata = strata, treatment = treatment, nsubjects = nsubjects,
                     nresponders = nresponders) %>% arrange(strata)
  } else {
    indata <- data %$% tibble(strata = strata, treatment = treatment, nsubjects = nsubjects,
                              nresponders = nresponders) %>% arrange(strata)
  }


  ss <- nresponders %>% subset(treatment == case) %>% sum()

  xx <- indata %>% filter(treatment == case) %>% select(nresponders)

  zz <- indata %>%
    group_by(strata) %>%
    summarise(tot_responders = sum(nresponders)) %>%
    .$tot_responders

  mm <- indata %>%
    filter(treatment == case) %>%
    .$nsubjects

  nn <- indata %>%
    group_by(strata) %>%
    summarise(tot_subjects = sum(nsubjects)) %>%
    .$tot_subjects

  m_lower <- pmax(0, zz-nn+mm)
  m_upper <- pmin(zz, mm)
  m_counts <- c(1, m_upper - m_lower + 1)  ## all plausible counts


  # x_range <- vector(mode = "list", length = nstrata)
  # for(i in 1:nstrata){
  #   x_range[[i]] <- m_lower[i]:m_upper[i]
  # }

  ## imaginary combination of underlying table: expand.grid(x_range)

  x_prob <- x_range <- vector(mode = "list", length = nstrata)
  for(i in 1:nstrata){
    x_range[[i]] <- m_lower[i]:m_upper[i]
  }

  for(i in 1:nstrata){
    x_prob[[i]] <- dhyper(x_range[[i]], m = mm[i], n = nn[i] - mm[i], k = zz[i], log=TRUE)
  }


  orep <- prod(m_counts)

  if (is.infinite(orep)) {
    stop("Sample space too large for this algorithm to handle")
  }

  if (side == "up") {

    pvalue <- 0
    i <- 1
    while(i <= orep) {
      index <- x_samp <- rep(0, nstrata)
      for(j in 1:nstrata) {
        index[j] <- ceiling(mod_top(i, prod(m_counts[1:(j+1)])) /
                              (prod(m_counts[1:j])))
        x_samp[j] <- x_range[[j]][index[j]]
      }

      if (sum(x_samp) >= ss) {
        x_samp_prob <- rep(0, nstrata)
        for (j in 1:nstrata) {
          x_samp_prob[j] <- x_prob[[j]][index[j]]
        }
        pvalue <- pvalue + exp(sum(x_samp_prob))
      }

      i <- i + 1
    }

    return(pvalue)

  } else if (side == "lower") {
    pvalue <- 0

    i <- 1

    while(i <= orep) {
      index <- x_samp <- rep(0, nstrata)
      for(j in 1:nstrata) {
        index[j] <- ceiling(mod_top(i, prod(m_counts[1:(j+1)])) /
                              (prod(m_counts[1:j])))
        x_samp[j] <- x_range[[j]][index[j]]
      }

      if (sum(x_samp) <= ss) {
        x_samp_prob <- rep(0, nstrata)
        for (j in 1:nstrata) {
          x_samp_prob[j] <- x_prob[[j]][index[j]]
        }
        pvalue <- pvalue + exp(sum(x_samp_prob))
      }

      i <- i + 1
    }

    return(pvalue)

  } else if (side == "both") {
    pvalue_up <- pvalue_lower <- 0

    i <- 1
    while(i <= orep) {
      index <- x_samp <- rep(0, nstrata)
      for(j in 1:nstrata) {
        index[j] <- ceiling(mod_top(i, prod(m_counts[1:(j+1)])) /
                              (prod(m_counts[1:j])))
        x_samp[j] <- x_range[[j]][index[j]]
      }
      if (sum(x_samp) >= ss) {
        x_samp_prob <- rep(0, nstrata)
        for (j in 1:nstrata) {
          x_samp_prob[j] <- x_prob[[j]][index[j]]
        }
        pvalue_up <- pvalue_up + exp(sum(x_samp_prob))
      }

      if (sum(x_samp) <= ss) {
        x_samp_prob <- rep(0, nstrata)
        for (j in 1:nstrata) {
          x_samp_prob[j] <- x_prob[[j]][index[j]]
        }
        pvalue_lower <- pvalue_lower + exp(sum(x_samp_prob))
      }

      i <- i + 1
    }

    pvalue <- min(1, 2 * min(pvalue_up, pvalue_lower))

    return(pvalue)
  } else {
    stop("side should be one of \"up\", \"lower\" or \"both\".")
  }
}

avandia.strata <- rep(1:48, each = 2)

avandia.treatment <- rep(1:2, 48)

avandia.nsubjects <- c(357, 176,
                      391, 207,
                      774, 185,
                      213, 109,
                      232, 116,
                      43, 47,
                      121, 124,
                      110, 114,
                      382, 384,
                      284, 135,
                      294, 302,
                      563, 142,
                      278, 279,
                      418, 212,
                      395, 198,
                      203, 106,
                      104, 99,
                      212, 107,
                      138, 139,
                      196, 96,
                      122, 120,
                      175, 173,
                      56, 58,
                      39, 38,
                      561, 276,
                      116, 111,
                      148, 143,
                      231, 242,
                      89, 88,
                      168, 172,
                      116, 61,
                      1172, 377,
                      706, 325,
                      204, 185,
                      288, 280,
                      254, 272,
                      314, 154,
                      162, 160,
                      442, 112,
                      394, 124,
                      2635, 2634,
                      1456, 2895,
                      70, 75,
                      25, 24,
                      232, 115,
                      101, 51,
                      196, 195,
                      676, 225)

avandia.mi <- c(2, 0,
                2, 1,
                1, 1,
                0, 1,
                1, 0,
                0, 1,
                1, 0,
                5, 2,
                1, 0,
                1, 0,
                0, 1,
                2, 0,
                2, 1,
                2, 0,
                2, 1,
                1, 1,
                1, 2,
                2, 0,
                3, 1,
                0, 0,
                0, 1,
                0, 1,
                1, 0,
                1, 0,
                0, 2,
                2, 3,
                1, 0,
                1, 0,
                1, 0,
                1, 0,
                0, 0,
                1, 0,
                0, 0,
                1, 2,
                1, 0,
                1, 0,
                1, 0,
                0, 0,
                1, 0,
                1, 0,
                15, 9,
                27, 41,
                0, 0,
                0, 0,
                0, 0,
                0, 0,
                0, 0,
                0, 0)

avandia.cvd <- c(1, 0,
                 0, 0,
                 0, 0,
                 0, 0,
                 1, 0,
                 0, 0,
                 0, 0,
                 3, 2,
                 0, 0,
                 0, 0,
                 2, 1,
                 0, 0,
                 0, 1,
                 0, 0,
                 2, 0,
                 1, 1,
                 0, 0,
                 1, 0,
                 1, 0,
                 1, 0,
                 0, 0,
                 0, 0,
                 0, 0,
                 0, 0,
                 1, 0,
                 2, 1,
                 2, 0,
                 1, 0,
                 0, 0,
                 1, 0,
                 0, 0,
                 1, 0,
                 1, 0,
                 0, 1,
                 1, 0,
                 0, 0,
                 0, 0,
                 0, 0,
                 1, 0,
                 1, 0,
                 12, 10,
                 2, 5,
                 0, 0,
                 0, 0,
                 0, 0,
                 0, 0,
                 0, 0,
                 0, 0)


set.seed(2542)
tt1 <- system.time(mi.pvalue <- sfet_naive_mc(strata = avandia.strata, treatment = avandia.treatment,
                  nsubjects = avandia.nsubjects, nresponders = avandia.mi,
                  case = 1, side = "both", sample_size = 50000))
## p-value: 0.028
## user.time: 44.77s
tt2 <- system.time(cvd.pvalue <- sfet_naive_mc(strata = avandia.strata, treatment = avandia.treatment,
                                              nsubjects = avandia.nsubjects, nresponders = avandia.cvd,
                                              case = 1, side = "both", sample_size = 50000))
## p-value: 0.04224
## user.time: 43.85s



# for MI, number of possible sample points is 1.95 * 10^18
# for CVD, number of possible sample points is 1.56 * 10^10
tt3 <- system.time(mi.pvalue <- SFET_extend(strata = avandia.strata, treatment = avandia.treatment,
                                            nsubjects = avandia.nsubjects, nresponders = avandia.mi,
                                            case = 1, side = "both"))
## p-value: 0.028
## user.time: 44.77s
tt4 <- system.time(cvd.pvalue <- SFET_extend(strata = avandia.strata, treatment = avandia.treatment,
                                             nsubjects = avandia.nsubjects, nresponders = avandia.cvd,
                                             case = 1, side = "both"))
## p-value: 0.04224
## user.time: 43.85s
save.image(file = "avandia_sfet.RData")
