
library(dplyr)

library(magrittr)
## get critical value for stratified_fisher_exact_test

## alpha: threshold of p-value, default 0.05
SFET_c_alpha <- function(strata, treatment, nsubjects, nresponders,
                         case, data = NULL, side = c("up"), alpha = .05){

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

  ## initialize pvalue and the iterating machine
  pmf_vec <- sum(mlower):sum(mupper) * 0

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
      if (side == "up") {
        upper_tail_prob <- pmf_vec %>% rev %>% cumsum %>% rev
        c_alpha <- (upper_tail_prob <= alpha) %>% which %>% min %>%
          magrittr::add(sum(mlower) - 1)
        return(c_alpha)
      } else if (side == "lower") {
        lower_tail_prob <- pmf_vec %>% cumsum
        c_alpha <- (lower_tail_prob <= alpha) %>% which %>% max %>%
          magrittr::add(sum(mlower) + 1)
        return(c_alpha)
      } else {
        upper_tail_prob <- pmf_vec %>% rev %>% cumsum %>% rev
        c_alpha_up <- (upper_tail_prob <= alpha/2) %>% which %>% min %>%
          magrittr::add(sum(mlower) - 1)
        lower_tail_prob <- pmf_vec %>% cumsum
        c_alpha_lower <- (lower_tail_prob <= alpha/2) %>% which %>% max %>%
          magrittr::add(sum(mlower) + 1)
        return(c(c_alpha_lower, c_alpha_up))
      }
    }

    x_samp_prob <- rep(NA, nstrata)

    for (j in 1:nstrata) {
      x_samp_prob[j] <- x_prob[[j]][current[j]]
    }
    pmf_vec[sum(current) - nstrata + 1] %<>%
      magrittr::add(exp(sum(x_samp_prob)))


    current[1] <- current[1] + 1
  }
}


strata <- rep(1:3, each = 2)
treatment <- rep(1:2, 3)
nsubjects <- c(670, 530, 500, 710, 490, 610)
nresponders <- c(21, 22, 19, 16, 20, 24)
case <- 1

c_up <- SFET_c_alpha(strata, treatment, nsubjects, nresponders, case, data = NULL,
                        side = "up", alpha = 0.05)
## 68
c_lower <- SFET_c_alpha(strata, treatment, nsubjects, nresponders, case, data = NULL,
                     side = "lower", alpha = 0.05)
## 50

c_twoside <- SFET_c_alpha(strata, treatment, nsubjects, nresponders, case, data = NULL,
                     side = "both", alpha = 0.05)
## 49 70
