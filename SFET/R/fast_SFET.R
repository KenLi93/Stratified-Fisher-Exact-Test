fast_SFET <- function(strata, treatment, nsubjects, nresponders,
                         case, data = NULL, prec = 1e-4, side = c("up")){

  nstrata <- strata %>% unique() %>% length()

  if (data %>% is.null()) {
    indata <- tibble(strata = strata, treatment = treatment, nsubjects = nsubjects,
                     nresponders = nresponders) %>% arrange(strata)
  } else {
    indata <- data %$% tibble(strata = strata, treatment = treatment, nsubjects = nsubjects,
                              nresponders = nresponders) %>% arrange(strata)
  }


  ss <- nresponders %>% subset(treatment == case) %>% sum()

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

  mlower <- pmax(0, zz-nn+mm)
  mupper <- pmin(zz, mm)
  mcounts <- mupper - mlower + 1  ## all plausible counts



  # x_range <- vector(mode = "list", length = nstrata)
  # for(i in 1:nstrata){
  #   x_range[[i]] <- m_lower[i]:m_upper[i]
  # }

  ## imaginary combination of underlying table: expand.grid(x_range)

  x_prob <- x_range <- vector(mode = "list", length = nstrata)

  for(i in 1:nstrata){
    x_range[[i]] <- mlower[i]:mupper[i]
  }

  ## calculate the hypergeometric density in each 2*2 table
  for(i in 1:nstrata){
    x_prob[[i]] <- dhyper(x_range[[i]], m = mm[i], n = nn[i] - mm[i], k = zz[i], log=TRUE)
  }


  ## acceleration, reduce the sample space for computation
  ## compute the maximum log probability that can be reduced from each stratum
  log_pthresh <- log(1 - (1-prec)^(1 / nstrata))

  ## for each stratum, remove the points that are unlikely
  for (i in 1:nstrata) {
    remained <- x_prob[[i]] >= log_pthresh - log(mcounts[i])
    mcounts[i] <- sum(remained)
    x_range[[i]] <- x_range[[i]][remained]
    x_prob[[i]] <- x_prob[[i]][remained]
    mlower[i] <- x_range[[i]][1]
    mupper[i] <- x_range[[i]][mcounts[i]]
  }






  ## initialize pvalue and the iterating machine
  pvalue <- pvalue_up <- pvalue_lower <- 0

  current <- rep(1, nstrata)

  while(1) {
    #    print(current)
    #    print(pvalue)
    ## iterating over all combinations of indices
    for (i in 1:(nstrata - 1)) {
      if (current[i] > mcounts[i]) {
        current[i] <- 1
        current[i+1] <- current[i+1] +1
      }
    }

    if (current[nstrata] > mcounts[nstrata]) {
      if (side == "up"| side =="lower") {
        return(pvalue)
      } else {
        pvalue <- min(1, 2 * min(pvalue_up, pvalue_lower))
        return(pvalue)
      }
    }

    if (side == "up") {
      if (sum(current) + sum(mlower) - nstrata >= ss) {
        x_samp_prob <- rep(NA, nstrata)

        for (j in 1:nstrata) {
          x_samp_prob[j] <- x_prob[[j]][current[j]]
        }
        pvalue <- pvalue + exp(sum(x_samp_prob))
      }
    } else if (side == "lower") {
      if (sum(current) + sum(mlower) - nstrata <= ss) {
        x_samp_prob <- rep(NA, nstrata)

        for (j in 1:nstrata) {
          x_samp_prob[j] <- x_prob[[j]][current[j]]
        }
        pvalue <- pvalue + exp(sum(x_samp_prob))

      }

    } else if (side == "both") {
      ## up-tail probability
      if (sum(current) + sum(mlower) - nstrata >= ss) {
        x_samp_prob <- rep(NA, nstrata)

        for (j in 1:nstrata) {
          x_samp_prob[j] <- x_prob[[j]][current[j]]
        }
        pvalue_up <- pvalue_up + exp(sum(x_samp_prob))
      }

      ## lower-tail probability
      if (sum(current) + sum(mlower) - nstrata <= ss) {
        x_samp_prob <- rep(NA, nstrata)

        for (j in 1:nstrata) {
          x_samp_prob[j] <- x_prob[[j]][current[j]]
        }
        pvalue_lower <- pvalue_lower + exp(sum(x_samp_prob))
      }


    } else {
      stop("side should be one of \"up\", \"lower\" or \"both\".")
    }

    current[1] <- current[1] + 1
  }
}


# strata <- rep(1:3, each = 2)
# treatment <- rep(1:2, 3)
# nsubjects <- c(670, 530, 500, 710, 490, 610)
# nresponders <- c(21, 22, 19, 16, 20, 24)
# case <- 1
# prec <- 1e-4
#
# tt1 <- system.time(test.p.exact.up <- fast_SFET(strata, treatment, nsubjects, nresponders, case, side = "up"))
# # p-value: 0.395
# # user  system elapsed
# # 0.43    0.00    0.44
#
# tt2 <- system.time(test.p.exact.lower <- fast_SFET(strata, treatment, nsubjects, nresponders, case, side = "lower"))
# # p-value: 0.675
# # user  system elapsed
# # 0.22    0.00    0.22
#
# tt3 <- system.time(test.p.exact.both <- fast_SFET(strata, treatment, nsubjects, nresponders, case, side = "both"))
# # p-value: 0.79
# # user  system elapsed
# # 0.36    0.00    0.36
