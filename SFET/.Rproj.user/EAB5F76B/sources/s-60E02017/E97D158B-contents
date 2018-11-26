#' Stratified Fisher Exact Test for Count Data
#'
#' Performs Fisher's exact test for testing the null of independence of
#' treatments and outcomes in contingency tables stratified by a categorical variable
#'
#' Stratified Fisher Exact Test provides a 'global test' for the association
#' between the treatment and response rate across strata.
#' In some cases sfet can be slow. Approximation by network algorithm or by the
#' analogous of Sodhi 2017 may be applied to sacrifice some accuracy for efficiency
#'
#'
#'
#'
#' @import dplyr
#'
#'
#' @param strata a vector specifying the strata of each treatment arm
#' @param treatment label of each arm, for now only support two-level factor vector
#' @param nsubjects an integer vector specifying the number of subjects in each arm
#' @param nresponders an integer vector specifying the number of responders in each arm
#' @param case a vector of length one specifying which label corresponds to the case group (as opposed to the control group)
#' @param data specify the parent data frame
#' @param side character specifying whether the function gives upper-tail, lower-tail or both-tail p-value
#'
#' @return The p-value of an upper-/lower-tailed or both-sided stratified fisher exact test.
#'
#' @export
#'
#' @examples
#' strata = c(1,1,2,2,3,3)
#' treatment = c(1,2,1,2,1,2)
#' nsubjects = c(11, 13, 9, 12, 8, 10)
#' nresponders = c(10, 12, 9, 11, 8, 7)
#' case = 1
#'
#' system.time(pv.up <- SFET(strata, treatment, nsubjects, nresponders, case = 2, side = "up"))
#' system.time(pv.lower <- SFET(strata, treatment, nsubjects, nresponders, case = 2, side = "lower"))
#' system.time(pv.both <- SFET(strata, treatment, nsubjects, nresponders, case = 2, side = "both"))

## case: specify which group is case

SFET <- function(strata, treatment, nsubjects, nresponders,
                 case, data = NULL, side = "up"){

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
  x_prob <- x_range <- vector(mode = "list", length = nstrata)
  for(i in 1:nstrata){
    x_range[[i]] <- m_lower[i]:m_upper[i]
  }

  for(i in 1:nstrata){
    x_prob[[i]] <- dhyper(x_range[[i]], m = mm[i], n = nn[i] - mm[i], k = zz[i], log=TRUE)
  }

  #  pv <- 0
  x_grid <- x_range %>% expand.grid() %>% rowSums()
  x_prob_grid <- x_prob %>% expand.grid() %>% rowSums() %>% exp()

  if (side == "up") {
    pv <- sum( x_prob_grid[x_grid >= ss] )
  } else if (side == "lower") {
    pv <- sum( x_prob_grid[x_grid <= ss] )
  } else if (side == "both") {
    pv.up <- sum( x_prob_grid[x_grid >= s] )
    pv.lower <- sum( x_prob_grid [x_grid <= ss] )
    pv <- min(1, 2 * min(pv.up, pv.lower))
  } else {
    stop("side should be one of \"up\", \"lower\" or \"both\".")
  }
  return(pv)
}


sfet_bayes <- function(strata, treatment, nsubjects, nresponders,
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


  event_grid <- matrix(0, nrow = sample_size, ncol = nstrata)
  for(i in 1:nstrata){
    event_grid[, i] <- rhyper(nn = sample_size,
                              m = zz[i],
                              n = nn[i] - zz[i],
                              k = nn[i] - mm[i])
  }

  event_incidence <- event_grid %>%
    apply(1, . %>% {
      any(. <= zz - xx) %>%
        as.numeric()
      })

  return(event_incidence %>% mean())
}

sfet_naive_mc <- function(strata, treatment, nsubjects, nresponders,
                          case, data = NULL, side = "up",
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
      return(mean(ss_null >= ss))
    } else if (side == "lower") {
      return(mean(ss_null < ss))
    } else if (side == "both") {
      return(min(1, 2 * min(mean(ss_null >= ss), mean(ss_null <= ss))))
    }
}

rm(list = ls())
strata <- rep(1:3, each = 2)
treatment <- rep(1:2, 3)
nsubjects <- c(670, 530, 500, 710, 490, 610)
nresponders <- c(21, 22, 19, 16, 20, 24)
case <- 1
tt1 <- system.time(test.p.bayes <- sfet_bayes(strata, treatment, nsubjects, nresponders, case))
## 0.93915
## 4.91

tt2 <- system.time(test.p.exact <- SFET(strata, treatment, nsubjects, nresponders, case))
## 0.395
## 0.17

tt3 <- system.time(test.p.mc <- sfet_naive_mc(strata, treatment, nsubjects, nresponders,
                                                          case, data = NULL, side = c("up"),
                                                          sample_size = 20000))
## 7.11
## 0.736

save.image()
