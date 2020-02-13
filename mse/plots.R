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


pm_angle <- theme(axis.text.x.top = element_text(angle = 45, hjust = 0))
cobweb <- geom_line(alpha = 0.15, position = position_dodge(width = 0.6))

plots$tigure_all_scenarios_avg
plots$tigure_minimum + pm_angle
plots$tigure_refset + pm_angle
plots$tigure_robset + pm_angle
plots$kobe
plots$convergence + coord_cartesian(ylim = c(-0.02, 1.02), expand = FALSE)
plots$worms_proj
plots$worms_hist_proj
plots$tradeoff_refset
plots$tradeoff_robset
plots$radar_refset_avg
plots$radar_refset
plots$lollipops_refset + cobweb
plots$lollipops_robset + cobweb
plots$lollipops_refset_avg + cobweb
plots$lollipops_refset
plots$lollipops_robset
plots$lollipops_refset_avg
plots$parallel_refset_avg
plots$parallel_refset
plots$parallel_refset_avg
plots$projections$base
plots$projections$upweight_dogfish
plots$projections$low_catch
plots$projections$sporadic_recruitment

# -----------------------------------------------------------------------------
# FIXME: pull into package:

mse2 <- mse
names(mse2) <- sc$scenario_human

bbmsy_ffmsy <- purrr::map_dfr(names(mse2), ~ {
  ts_data <- gfdlm:::get_ts(object = mse2[[.]], type = c("SSB", "FM"), this_year = 2019)
  gfdlm:::get_ts_quantiles(ts_data, probs = c(0.1, 0.5)) %>%
    mutate(scenario = .x)
})

catch <- purrr::map_dfr(names(mse2), ~ {
  ts_data <- gfdlm:::get_ts(object = mse2[[.]], type = c("C"), this_year = 2019)
  gfdlm:::get_ts_quantiles(ts_data, probs = c(0.1, 0.5)) %>%
    mutate(scenario = .x)
})

this_year <- 2019
rel_widths <- c(2, 1.2)
bbmsy_ffmsy$Type <- gsub("_", "/", bbmsy_ffmsy$Type)
bbmsy_ffmsy$Type <- gsub("MSY", "[MSY]", bbmsy_ffmsy$Type)

g1 <- bbmsy_ffmsy %>%
  filter(scenario != "Pinniped mortality") %>%
  filter(mp_name %in% c(c("CC_10t", "CC_15t", "CC_5t"), "FMSYref", "NF")) %>%
  ggplot(aes_string("real_year", "m", colour = "scenario", fill = "scenario")) +
  geom_line(na.rm = TRUE) +
  facet_grid(mp_name ~ Type, labeller = ggplot2::label_parsed) +
  geom_ribbon(aes_string(x = "real_year", ymin = "l", ymax = "u"),
    colour = NA, alpha = 0.07, na.rm = TRUE
  ) +
  theme_pbs() +
  coord_cartesian(expand = FALSE, ylim = c(0, 5)) +
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + ylab("Value") + xlab("Year") +
  geom_vline(xintercept = this_year, lty = 2, alpha = 0.3) +
  ggplot2::theme(panel.spacing = grid::unit(-0.1, "lines"))

g2 <- catch %>%
  filter(scenario != "Pinniped mortality") %>%
  filter(mp_name %in% c(c("CC_10t", "CC_15t", "CC_5t"), "FMSYref", "NF")) %>%
  ggplot(aes_string("real_year", "m", colour = "scenario", fill = "scenario")) +
  geom_line(na.rm = TRUE) +
  facet_grid(mp_name ~ Type) +
  geom_ribbon(aes_string(x = "real_year", ymin = "l", ymax = "u"),
    colour = NA, alpha = 0.07, na.rm = TRUE
  ) +
  theme_pbs() +
  # coord_cartesian(expand = FALSE, ylim = c(0, 5)) +
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + ylab("Catch") + xlab("Year") +
  geom_vline(xintercept = this_year, lty = 2, alpha = 0.3) +
  ggplot2::theme(panel.spacing = grid::unit(-0.1, "lines"))

g3 <- cowplot::plot_grid(
  g1 + theme(legend.position = "none"),
  g2 + theme(legend.position = "none"),
  rel_widths = rel_widths, align = "h"
)

legend <- cowplot::get_legend(
  # create some space for the legend
  g1 + theme(legend.box.margin = margin(0.2, 0.2, 12, .2), legend.position = "bottom")
)

g <- cowplot::plot_grid(g3, legend, rel_heights = c(4, .2), nrow = 2)

ggsave("~/Desktop/sc.pdf", width = 8, height = 8)
