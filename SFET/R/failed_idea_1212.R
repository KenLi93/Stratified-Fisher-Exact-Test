library(dplyr)
library(magrittr)
load("./R/small_tab.RData")

## 0-0 tables have no contribution to the p-value
## 21 0-1 tables and one big table
nstrata <- small_tab$strata %>% unique %>% length

ss <- small_tab$nresponders %>% subset(treatment == case) %>% sum()

zz <- small_tab %>%
  group_by(strata) %>%
  summarise(tot_responders = sum(nresponders)) %>%
  .$tot_responders

mm <- small_tab %>%
  filter(treatment == case) %>%
  .$nsubjects

nn <- small_tab %>%
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
  x_prob[[i]] <- dhyper(x_range[[i]], m = mm[i], n = nn[i] - mm[i], k = zz[i], log = F)
}


pvalue <- 0

q0 <- x_prob %>% sapply(. %>% .[1]) %>% max(.[1:21])
q1 <- x_prob %>% sapply(. %>% .[2]) %>% max(.[1:21])


for(tt in max(ss - nstrata + 1, mlower[nstrata]):mupper[nstrata]) {
  big_prob <- dhyper(tt, m = mm[nstrata], n = nn[nstrata] - mm[nstrata], k = zz[nstrata])
  small_prob <- 0
  for(rr in max(ss - tt, 0):(nstrata - 1)) {
    small_prob <- small_prob + choose(nstrata - 1, rr) * q1^rr *
      q0^(nstrata - 1 - rr)
  }
  print(small_prob)
  pvalue <- pvalue + big_prob * small_prob
}

## pvalue: 19.74971
for(tt in )
