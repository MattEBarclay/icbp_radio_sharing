## Meta-analysis - dotplot of chemo/radio by country and cancer

library(metafor)
library(tidyverse)
library(broom)
library(gt)
library(lme4)

# Useful vectors for sorting things out
load("Rdata/useful_vectors.Rdata")

# Relevant data file
trt_unadj_n <- readRDS("Rdata/trt_unadj_n.RDS")

# Useful functions
source("radio_analysis0_repeated_functions.R")

# Create dataset for meta-analysis ----------------------------------------

# Pull out age split for each cancer site
# We do not have an age/sex split
dat <- trt_unadj_n |>
  filter(
    stage == "All stages",
    ((variable == "age") & (cancer != "All 8 cancers")),
    jurisdiction != "",
    trt == "Radio"
  ) |>
  mutate(age = case_when(
    variable == "age" & variable_value == "64" ~ "15-64",
    variable == "age" & variable_value == "74" ~ "65-74",
    variable == "age" & variable_value == "84" ~ "75-84",
    variable == "age" & variable_value == "99" ~ "85-99",
    TRUE ~ ""
  )) |>
  mutate(age = factor(age, levels = c("15-64","65-74","75-84","85-99"))) |>
  select(trt, jurisdiction, ordering, cancer, age, n, n_trt, prop, lower, upper) 

# apply common cleaning and exclusion rules
dat <- dat |>
  mutate(site_order  = factor(cancer, levels = site_rad)) |>
  radio_dat_clean() |>
  radio_exclusions_all() |>
  arrange(site_order, ordering)



# IPD meta-analysis -------------------------------------------------------
# Well,  "IPD" - adjusted for cancer and age

all_unadj   <- glmer(cbind(n_trt,n-n_trt) ~ (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(jurisdiction != "Prince Edward Island") |> filter(jurisdiction != "Saskatchewan"))

all_cas     <- glmer(cbind(n_trt,n-n_trt) ~ site_order + (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(jurisdiction != "Prince Edward Island") |> filter(jurisdiction != "Saskatchewan"))

all_cas_age <- glmer(cbind(n_trt,n-n_trt) ~ site_order*age + (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(jurisdiction != "Prince Edward Island") |> filter(jurisdiction != "Saskatchewan"))

summary(all_unadj)
summary(all_cas)
summary(all_cas_age)

oes_all <- glmer(cbind(n_trt,n-n_trt) ~       (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(site_order == "Oesophageal") |> filter(!exclude))
oes_age <- glmer(cbind(n_trt,n-n_trt) ~ age + (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(site_order == "Oesophageal") |> filter(!exclude))

sto_all <- glmer(cbind(n_trt,n-n_trt) ~       (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(site_order == "Stomach") |> filter(!exclude))
sto_age <- glmer(cbind(n_trt,n-n_trt) ~ age + (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(site_order == "Stomach") |> filter(!exclude))

col_all <- glmer(cbind(n_trt,n-n_trt) ~       (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(site_order == "Colon") |> filter(!exclude))
col_age <- glmer(cbind(n_trt,n-n_trt) ~ age + (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(site_order == "Colon") |> filter(!exclude))

rec_all <- glmer(cbind(n_trt,n-n_trt) ~       (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(site_order == "Rectal") |> filter(!exclude))
rec_age <- glmer(cbind(n_trt,n-n_trt) ~ age + (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(site_order == "Rectal") |> filter(!exclude))

liv_all <- glmer(cbind(n_trt,n-n_trt) ~       (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(site_order == "Liver") |> filter(!exclude))
liv_age <- glmer(cbind(n_trt,n-n_trt) ~ age + (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(site_order == "Liver") |> filter(!exclude))

pan_all <- glmer(cbind(n_trt,n-n_trt) ~       (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(site_order == "Pancreatic") |> filter(!exclude))
pan_age <- glmer(cbind(n_trt,n-n_trt) ~ age + (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(site_order == "Pancreatic") |> filter(!exclude))

lun_all <- glmer(cbind(n_trt,n-n_trt) ~       (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(site_order == "Lung") |> filter(!exclude))
lun_age <- glmer(cbind(n_trt,n-n_trt) ~ age + (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(site_order == "Lung") |> filter(!exclude))

ova_all <- glmer(cbind(n_trt,n-n_trt) ~       (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(site_order == "Ovarian") |> filter(!exclude))
ova_age <- glmer(cbind(n_trt,n-n_trt) ~ age + (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(site_order == "Ovarian") |> filter(!exclude))

# list RE standard deviations
VarCorr(all_unadj)
VarCorr(all_cas)
VarCorr(all_cas_age)
VarCorr(oes_all)
VarCorr(oes_age)
VarCorr(sto_all)
VarCorr(sto_age)
VarCorr(col_all)
VarCorr(col_age)
VarCorr(rec_all)
VarCorr(rec_age)
VarCorr(liv_all)
VarCorr(liv_age)
VarCorr(pan_all)
VarCorr(pan_age)
VarCorr(lun_all)
VarCorr(lun_age)
VarCorr(ova_all)
VarCorr(ova_age)



# Clean up ----------------------------------------------------------------
rm(list = ls())

