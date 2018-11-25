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

library(dplyr)

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
