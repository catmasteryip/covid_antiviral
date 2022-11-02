# Introduction

This repo contains the R code on investigation of COVID antivirals treatment effect on organ dysfunction and mortality

# Notebook explanations

abbreviation: pp3 means paper3, pending publication on bacterial infection subjects

data_import.Rmd: raw data import and minor manipulation, libraries imports <br />
data_import_pp3.Rmd: data import notebook for paper 3 <br />
IP_data.Rmd: non-meds data wrangling <br />
IP_data_pp3.Rmd: non-meds data wrangling for paper 3 <br />
meds.Rmd: meds data wrangling <br />
culture_considered.Rmd: organ dysfunction and mortality dataframes considering bacterial culture status, paper 3 <br />
culture_not_considered.Rmd: organ dysfunction and mortality dataframes without considering bacterial culture status <br />
iptw_HR.Rmd: iptw event/risk, ir, ird and hr generation <br />
./scripts/iptwhelper.R: iptw helper functions <br />
./scripts/helper.R: helper functions of meds and comorbs imputation
