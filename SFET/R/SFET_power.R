library(dplyr)

library(magrittr)

## density of hypergeometric function when sampling propability of balls with different
## colors are unequal
dhyper_h1 <- function(x, m, n, k, log = FALSE, theta = 1) {
  if (x >= max(0, k-n)| x <= min(k, m)){
    denom <- 0
    rr <- max(0, k-n):min(k, m)
    p_x <- (choose(m, x) * choose(n, k-x) * theta^x) /
      sum(choose(m, rr) * choose(n, k-rr) * theta^rr)
  } else {
    p_x <- 0
  }
  if (log == TRUE) {
    return(log(p_x))
  } else {
    return(p_x)
  }
}
## get critical value for stratified_fisher_exact_test

## alpha: threshold of p-value, default 0.05
SFET_power <- function(strata, treatment, nsubjects, nresponders,
                       case, data = NULL, side = c("up"), alpha = .05,
                       theta = 1){

  nstrata <- strata %>% unique() %>% length()

  if(length(theta) < nstrata) {
    theta <- rep_len(theta, length.out = nstrata)
  }

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

  x_prob <- x_range <-
    x_prob_h1 <- vector(mode = "list", length = nstrata)

  for(i in 1:nstrata){
    x_range[[i]] <- mlower[i]:mupper[i]
  }

  ## calculate the hypergeometric density in each 2*2 table
  for(i in 1:nstrata){
    x_prob[[i]] <- dhyper(x_range[[i]], m = mm[i], n = nn[i] - mm[i], k = zz[i], log=TRUE)
    for(j in 1:mcounts[i]) {
      x_prob_h1[[i]][j] <- dhyper_h1(x_range[[i]][j], m = mm[i], n = nn[i] - mm[i], k = zz[i],
                                     theta = theta[i], log=TRUE)
    }
  }

  ## initialize pvalue and the iterating machine
  pmf_vec <- sum(mlower):sum(mupper) * 0
  pmf_vec_h1 <- sum(mlower):sum(mupper) * 0

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

    ## after the iteration, obtain the critical values and calculate the conditional power

    if (current[nstrata] > mcounts[nstrata]) {
      break();
    }

    x_samp_prob <- x_samp_prob_h1 <- rep(NA, nstrata)

    for (j in 1:nstrata) {
      x_samp_prob[j] <- x_prob[[j]][current[j]]
      x_samp_prob_h1[j] <- x_prob_h1[[j]][current[j]]
    }
    pmf_vec[sum(current) - nstrata + 1] %<>%
      magrittr::add(exp(sum(x_samp_prob)))
    pmf_vec_h1[sum(current) - nstrata + 1] %<>%
      magrittr::add(exp(sum(x_samp_prob_h1)))


    current[1] <- current[1] + 1
  }

  if (side == "up") {
    upper_tail_prob <- pmf_vec %>% rev %>% cumsum %>% rev
    c_alpha <- (upper_tail_prob <= alpha) %>% which %>% min %>%
      magrittr::add(sum(mlower) - 1)
    power <- pmf_vec_h1[(c_alpha - sum(mlower) + 1):(sum(mupper) + 1)] %>% sum
    return(power)
  } else if (side == "lower") {
    lower_tail_prob <- pmf_vec %>% cumsum
    c_alpha <- (lower_tail_prob <= alpha) %>% which %>% max %>%
      magrittr::add(sum(mlower) + 1)
    power <- pmf_vec_h1[1:(c_alpha - sum(mlower) + 1)] %>% sum
    return(power)
  } else {
    upper_tail_prob <- pmf_vec %>% rev %>% cumsum %>% rev
    c_alpha_up <- (upper_tail_prob <= alpha/2) %>% which %>% min %>%
      magrittr::add(sum(mlower) - 1)
    lower_tail_prob <- pmf_vec %>% cumsum
    c_alpha_lower <- (lower_tail_prob <= alpha/2) %>% which %>% max %>%
      magrittr::add(sum(mlower) + 1)
    power <- pmf_vec_h1[-((c_alpha_lower - sum(mlower) + 2):
                          (c_alpha_up - sum(mlower)))] %>% sum
    return(power)
  }
}


strata <- rep(1:3, each = 2)
treatment <- rep(1:2, 3)
nsubjects <- c(670, 530, 500, 710, 490, 610)
nresponders <- c(21, 22, 19, 16, 20, 24)
case <- 1

power_up <- SFET_power(strata, treatment, nsubjects, nresponders,
                       case, data = NULL, side = "up", alpha = .05,
                       theta = 1.5)
## 0.6657976
power_lower <- SFET_power(strata, treatment, nsubjects, nresponders,
                          case, data = NULL, side = "lower", alpha = .05,
                          theta = 1.5)
## 0.0001635896
power_both <- SFET_power(strata, treatment, nsubjects, nresponders,
                         case, data = NULL, side = "both", alpha = .05,
                         theta = 1.5)
## 0.5216413

power_up_h0 <- SFET_power(strata, treatment, nsubjects, nresponders,
                       case, data = NULL, side = "up", alpha = .05,
                       theta = 1)
## 0.03980743
