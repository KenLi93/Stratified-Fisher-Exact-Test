
library(magrittr)
library(dplyr)

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


# strata <- rep(1:3, each = 2)
# treatment <- rep(1:2, 3)
# nsubjects <- c(670, 530, 500, 710, 490, 610)
# nresponders <- c(21, 22, 19, 16, 20, 24)
# case <- 1
#
#
# tt1 <- system.time(test.p.exact.up <- SFET_extend(strata, treatment, nsubjects, nresponders, case, side = "up"))
# # p-value: 0.395
# # user  system elapsed
# # 2.54     0.0     2.54
#
# tt2 <- system.time(test.p.exact.lower <- SFET_extend(strata, treatment, nsubjects, nresponders, case, side = "lower"))
# # p-value: 0.675
# # user  system elapsed
# # 2.39    0.00    2.39
#
# tt3 <- system.time(test.p.exact.both <- SFET_extend(strata, treatment, nsubjects, nresponders, case, side = "both"))
# # p-value: 0.79
# # user  system elapsed
# # 2.53    0.00    2.54
