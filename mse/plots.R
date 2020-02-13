library(purrr)
library(dplyr)
library(ggplot2)
library(gfdlm)

sc <- tibble::tribble(
  ~scenario,             ~scenario_human,            ~scenario_type,
  "base",                 "Base",                     "Reference",
  "upweight_dogfish",     "Upweight dogfish survey",  "Reference",
  "low_catch",            "Low catch",                "Reference",
  "sporadic_recruitment", "Episodic recruitment",     "Reference",
  "lowM",                 "Low M (M=0.02)",           "Robustness",
  "pinniped",             "Pinniped mortality",       "Robustness"
)

mse <- map(sc$scenario,
  ~readRDS(paste0("data-generated/mse/MSE_", ., ".rds")))
names(mse) <- sc$scenario

reference_mp <- c("FMSYref", "NFref")
mps <- c(
  "FMSYref", "NF", # rename NF to NFref; FMSYref75?
  "CC_5t", "CC_10t", "CC_15t", # keep
  "ICI_YE", "ICI2_YE", # cut?
  "Iratio_YE", # keep
  "GB_slope_YE", # keep
  "IT5_YE", "IT10_YE", # what is the target?
  "Islope_YE", # 4 versions?
  "IDX_YE", "IDX_smooth_YE", # add floor
  "SP_4080_5f", "SP_4080_10f",
    "SP_2060_5f", "SP_2060_10f",
  "SP_interim" # explain priors
  # add:
  #   - Itargets
  #   - GB_slopes?
)

# Set up PMs ------------------------------------------------------------------

`1GT LRP` <- gfdlm::pm_factory("SBMSY", 0.4, c(38, 38))
`1.5GT LRP` <- gfdlm::pm_factory("SBMSY", 0.4, c(56, 56))
`1.5GT USR` <- gfdlm::pm_factory("SBMSY", 0.8, c(56, 56))

# `MT LRP` <- gfdlm::pm_factory("SBMSY", 0.4, c(38, 38))
# `LT LRP` <- gfdlm::pm_factory("SBMSY", 0.4, c(56, 56))
# `LT USR` <- gfdlm::pm_factory("SBMSY", 0.8, c(56, 56))
#
`LRP 1GT` <- gfdlm::pm_factory("SBMSY", 0.4, c(38, 38))
`LRP 1.5GT` <- gfdlm::pm_factory("SBMSY", 0.4, c(56, 56))
`USR 1.5GT` <- gfdlm::pm_factory("SBMSY", 0.8, c(56, 56))

FMSY <- DLMtool::PNOF
AAVC <- DLMtool::AAVY
STC <- gfdlm::pm_factory("LTY", 0.5, c(1, 10))
LTC <- gfdlm::pm_factory("LTY", 0.5, c(38, 38))
PM <- c("LRP 1.5GT", "USR 1.5GT", "LRP 1GT", "FMSY", "STC", "LTC", "AAVC")

custom_pal <- c(RColorBrewer::brewer.pal(3, "Set2"), "grey60", "grey30")
names(custom_pal) <- c(c("CC_10t", "CC_15t", "CC_5t"), "FMSYref", "NF")

plots <- plot_factory(
  mse,
  pm = PM,
  scenario_df = sc,
  this_year = 2019,
  mp_sat = c("CC_10t", "CC_15t", "CC_5t", "NF", "FMSYref"),
  mp_not_sat = c("CC_15t"),
  mp_not_sat_highlight = c("CC_15t"),
  mp_ref = c("FMSYref", "NF"),
  custom_pal = custom_pal,
  tradeoff = c("LRP 1.5GT", "STC"),
  eg_scenario = "upweight_dogfish",
  satisficed_criteria = c("LRP 1.5GT" = 0.8, "STC" = 0.7)
)
