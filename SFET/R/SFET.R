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
#'
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
                 case, data = NULL, side = c("up", "lower", "both")){

  nstrata <- length(unique(strata))

  data <- with(data, data.frame(strata = strata, treatment = treatment,
                     nsubjects = nsubjects, nresponders = nresponders))
  s <- sum(nresponders[treatment == case])
  data.strata <- split(data, data$strata)
  x <- sapply(data.strata, function(tab) with(tab, nresponders[treatment == case])) ## responders in the treatment group
  z <- sapply(data.strata, function(tab) with(tab, sum(nresponders))) ## total number of responders
  n <- sapply(data.strata, function(tab) with(tab, sum(nsubjects))) ## total number of subjects in each stratum
  m <- sapply(data.strata, function(tab) with(tab, nsubjects[treatment == case])) ## subjects in the treatment group

  m_lower <- pmax(0, z-n+m)
  m_upper <- pmin(z, m)
  x_prob <- x_range <- vector(mode = "list", length = nstrata)
  for(i in 1:nstrata){
    x_range[[i]] <- m_lower[i]:m_upper[i]
  }

  for(i in 1:nstrata){
    x_prob[[i]] <- dhyper(x_range[[i]], m = m[i], n = n[i] - m[i], k = z[i], log=TRUE)
  }

  #  pv <- 0
  x_grid <- rowSums(expand.grid(x_range))
  x_prob_grid <- exp(rowSums(expand.grid(x_prob)))

  if (side == "up") {
    pv <- sum( x_prob_grid [ x_grid>= s ] )
  } else if (side == "lower") {
    pv <- sum( x_prob_grid [ x_grid<= s ] )
  } else if (side == "both") {
    pv.up <- sum( x_prob_grid [ x_grid>= s ] )
    pv.lower <- sum( x_prob_grid [ x_grid<= s ] )
    pv <- min(1, 2 * min(pv.up, pv.lower))
  } else {
    stop("side should be one of \"up\", \"lower\" or \"both\".")
  }
  return(pv)
}


sfet.bayes <- function(strata, treatment, nsubjects, nresponders,
                       case, data = NULL, side = c("up", "lower", "both")){

  nstrata <- length(unique(strata))


  data <- with(data, data.frame(strata = strata, treatment = treatment,
                                nsubjects = nsubjects, nresponders = nresponders))
  s <- sum(nresponders[treatment == case])
  data.strata <- split(data, data$strata)
  x <- sapply(data.strata, function(tab) with(tab, nresponders[treatment == case])) ## responders in the treatment group
  z <- sapply(data.strata, function(tab) with(tab, sum(nresponders))) ## total number of responders
  n <- sapply(data.strata, function(tab) with(tab, sum(nsubjects))) ## total number of subjects in each stratum
  m <- sapply(data.strata, function(tab) with(tab, nsubjects[treatment == case])) ## subjects in the treatment group

  m_lower <- pmax(0, z-n+m)
  m_upper <- pmin(z, m)
  x_prob <- x_range <- vector(mode = "list", length = nstrata)
  for(i in 1:nstrata){
    x_range[[i]] <- m_lower[i]:m_upper[i]
  }

  for(i in 1:nstrata){
    x_prob[[i]] <- dhyper(x_range[[i]], m = m[i], n = n[i] - m[i], k = z[i], log=TRUE)
  }

  #  pv <- 0
  x_grid <- rowSums(expand.grid(x_range))
  x_prob_grid <- exp(rowSums(expand.grid(x_prob)))

  if (side == "up") {
    pv <- sum( x_prob_grid [ x_grid>= s ] )
  } else if (side == "lower") {
    pv <- sum( x_prob_grid [ x_grid<= s ] )
  } else if (side == "both") {
    pv.up <- sum( x_prob_grid [ x_grid>= s ] )
    pv.lower <- sum( x_prob_grid [ x_grid<= s ] )
    pv <- min(1, 2 * min(pv.up, pv.lower))
  } else {
    stop("side should be one of \"up\", \"lower\" or \"both\".")
  }
  return(pv)
}
