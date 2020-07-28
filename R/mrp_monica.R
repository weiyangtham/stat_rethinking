library(tidyverse)
library(magrittr)
library(haven)
library(ipumsr)

# 1. Read in and tidy up data ---------------------------------------------


# read in MNCS
# this file can be downloaded here: https://osf.io/uzqdn/
d <- read_dta("data/MNCS-PV2.dta")

# keep variables of interest and only those responses with a year and age

d <- d %>%
  select(yrmar,
         agemar,
         agemarc,
         genmar,
         spgenmar,
         namechg,
         ednow,
         state) %>%
  filter(!is.na(agemar), !is.na(yrmar))

# keep only hetero women, make age group, education and decade married variables

d <- d %>%
  filter(genmar==2&spgenmar==1) %>%
  mutate(kept_name = as.numeric(namechg==1),
         state_name = str_to_lower(as.character(factor(state,
                                                       levels = attributes(d$state)$labels,
                                                       labels = names(attributes(d$state)$labels)))),
         age = agemar + (2019-yrmar),
         age_group = (as.character(cut(age,
                                       breaks = c(seq(15, 80, by = 5), Inf),
                                       labels = seq(15, 80, by = 5), right = FALSE
         ))),
         decade_married = (as.character(cut(yrmar,
                                            breaks = c(seq(1969, 2019, by = 10), Inf),
                                            labels = seq(1969, 2019, by = 10), right = FALSE
         ))),
         educ_group = case_when(
           ednow<5 ~ "<BA",
           ednow==5 ~ "BA",
           ednow>5 ~ ">BA",
           TRUE ~ "NA"
         ))

d_mncs <- d %>%
  select(kept_name, state_name, age_group, decade_married, educ_group) %>%
  filter(age_group>20, age_group<80, decade_married>1969)

saveRDS(d_mncs, "data/d_mncs.RDS")

# ACS data

ddi <- read_ipums_ddi("data/usa_00002.xml")
dc <- read_ipums_micro(ddi)

dc %<>%
  janitor::clean_names()

d_acs <- dc %>%
  mutate(age_group = (as.character(cut(age,
                                       breaks = c(seq(15, 80, by = 5), Inf),
                                       labels = seq(15, 80, by = 5), right = FALSE
  ))),
  decade_married = (as.character(cut(yrmarr,
                                     breaks = c(seq(1969, 2019, by = 10), Inf),
                                     labels = seq(1969, 2019, by = 10), right = FALSE
  ))),
  educ_group = case_when(
    educ<10 ~ "<BA",
    educ==10~"BA",
    educ>10 ~">BA",
    TRUE~ "NA"
  ))

d_acs
