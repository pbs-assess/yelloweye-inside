library("DLMtool")
library("MSEtool")
library("dplyr")
library("purrr")
library("ggplot2")
library("gfdlm")
library("here")
# Settings --------------------------------------------------------------------
reference_mp <- c("FMSYref", "FMSYref75", "NFref")
sp <- "ye"

# Need functions to calculate:

# Set up PMs ------------------------------------------------------------------
`LRP 1.5GT` <- gfdlm::pm_factory("SBMSY", 0.4, c(56, 56))
`LRP 1GT` <- gfdlm::pm_factory("SBMSY", 0.4, c(38, 38))
`USR 1.5GT` <- gfdlm::pm_factory("SBMSY", 0.8, c(56, 56))
FMSY <- DLMtool::PNOF
`AAVC` <- DLMtool::AAVY
`ST C10` <- gfdlm::pm_factory("LTY", 0.5 - 1e-4, c(1, 10))
`ST C15` <- gfdlm::pm_factory("LTY", 0.75 - 1e-4, c(1, 10))
`LT C20` <- gfdlm::pm_factory("LTY", 1 - 1e-4, c(56, 56))

mse_temp <- readRDS("mse/om/MSE_upweight_dogfish.rds")
catch <- apply(mse_temp@CB_hist[1, , (102 - 8 + 1):102, ], 2, sum) # Catch since 2012
rm(mse_temp)
ref_aadc <- gfdlm:::get_aadc(catch)
`ST AADC` <- gfdlm::pm_factory("AADC", ref_aadc, c(1, 10))
# `AADC 1GT` <- gfdlm::pm_factory("AADC", ref_aadc, c(1, 38))

PM <- c("LRP 1.5GT", "USR 1.5GT", "LRP 1GT", "ST C10", "ST C15", "LT C20", "ST AADC")

# Set up and checks -----------------------------------------------------------
sc <- readRDS("mse/om/ye-scenarios.rds")
sc <- rename(sc, scenario_human = scenarios_human)
fig_dir <- "mse/figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir)

# # FIXME: temp.
# sc <- tibble::tribble(
#   ~scenario,     ~scenario_human,        ~scenario_type,
#   "upweight_dogfish_nr",    "(4) Estimate HBLL selectivity","Reference",
#   "upweight_dogfish",       "(4b) Estimate HBLL selectivity\n(Reconstructed catch until 2006)","Reference"
# )

.ggsave <- function(filename, width, height, ...) {
  ggsave(file.path(fig_dir, paste0(sp, "-", filename, ".png")),
    width = width, height = height, ...
  )
}

get_filtered_scenario <- function(type, column) {
  filter(sc, scenario_type == type) %>%
    pull(!!column) %>%
    set_names()
}
scenarios <- sc$scenario %>% set_names()
scenarios_human <- sc$scenario_human %>% set_names()
scenarios_ref <- get_filtered_scenario("Reference", "scenario")
scenarios_ref_human <- get_filtered_scenario("Reference", "scenario_human")
scenarios_rob <- get_filtered_scenario("Robustness", "scenario")
scenarios_rob_human <- get_filtered_scenario("Robustness", "scenario_human")

# Read OMs --------------------------------------------------------------------
# om <- lapply(scenarios, function(x) readRDS(paste0("mse/om/", x, ".rds"))@OM)
# names(om) <- scenarios

mse <- map(scenarios, ~ readRDS(paste0("mse/om/MSE_", .x, ".rds")))
for (i in seq_along(mse)) {
  mse[[i]]@OM$RefY <- 20
  mse[[i]]@MPs[mse[[i]]@MPs %in% c("IDX_", "IDX_smooth_")] <- c("IDX", "IDX_smooth")
}

# Report range in reference points for each OM
ref_pt <- lapply(mse, function(x) {
  do.call(rbind, lapply(x@Misc$MSYRefs$ByYear[1:2], function(xx) {
    mean_ref <- mean(xx[, 1, x@nyears])
    sd_ref <- sd(xx[, 1, x@nyears])
    structure(c(mean_ref, sd_ref/mean_ref), names = c("Mean", "CV"))
  }))
})

ref_pt_text <- do.call(rbind, lapply(ref_pt, function(x) {
  round(x, 2) %>% format(trim = TRUE) %>% apply(1, function(xx) paste0(xx[1], " (", xx[2], ")"))
  })) %>%
  structure(dimnames = list(scenarios_human, c("MSY", "FMSY")))
write.csv(ref_pt_text, file = "mse/figures/ref_pt.csv")
saveRDS(ref_pt_text, "mse/figures/ref_pt.rds")


# COSEWIC metric E, Probability that the biomass is at least 2, 5 % B0 within the projection period
COSEWIC_E <- function(MSEobj, Ref = 0.02, Yrs = c(1, 100)) {
  Yrs <- ChkYrs(Yrs, MSEobj)
  PMobj <- new("PMobj")
  PMobj@Name <- PMobj@Caption <- paste0("Probability that biomass is above ", Ref * 100, "% B0")
  PMobj@Stat <- MSEobj@SSB[, , Yrs[1]:Yrs[2]] / MSEobj@OM$SSB0
  PMobj@Ref <- Ref
  PMobj@Prob <- calcProb(PMobj@Stat >= PMobj@Ref, MSEobj) #
  PMobj@Mean <- calcMean(PMobj@Prob) # calculate mean probability by MP
  PMobj@MPs <- MSEobj@MPs
  PMobj
}
`2% B0` <- `5% B0` <- COSEWIC_E
formals(`5% B0`)$Ref <- 0.05

e_df_list <- map(mse, ~ gfdlm::get_probs(.x, c("2% B0", "5% B0"))) # List with all OMs and PMs

# For all scenarios
# walk(names(mse), ~ {
#   g <- plot_tigure(e_df_list[[.x]]) + ggtitle(sc$scenario_human[match(.x, sc$scenario)])
#   ggsave(paste0("mse/figures/tigure_cosewic_", .x, ".png"), width = 3.5, height = 6.5)
# })

# Averaged across reference set
e_avg <- e_df_list[sc$scenario_type == "Reference"] %>%
  bind_rows(.id = "scenario") %>%
  group_by(MP) %>%
  summarise_if(is.numeric, mean)
pm_df_list <- map(mse[scenarios_ref], ~ gfdlm::get_probs(.x, PM)) # List with all OMs and PMs

# g <- plot_tigure(e_avg) + ggtitle("Averaged across reference set")
# ggsave("mse/figures/tigure_cosewic_average_ref.png", width = 3.5, height = 6.5)

# Satisficing -----------------------------------------------------------------
pm_df_list <- map(mse[scenarios_ref], ~ gfdlm::get_probs(.x, PM)) # List with all OMs and PMs
pm_df_list_rob <- map(mse[scenarios_rob], ~ gfdlm::get_probs(.x, PM)) # Robustness only
pm_df <- bind_rows(pm_df_list, .id = "scenario") # All as a data.frame

saveRDS(pm_df, file = "mse/om/ye-pm-all.rds")

# # FIXME:
# pm_angle <- theme(
#   axis.text.x.top = element_text(angle = 60, hjust = 0)
# )
#
# x <- pm_df_list[[1]] %>% dplyr::filter(MP %in% pm_df_list[[2]]$MP) %>%
#   arrange(`LRP 1.5GT`, `USR 1.5GT`, `LRP 1GT`, `ST C10`)
# plot_tigure(x, mp_order = x$MP) + pm_angle
# ggsave("mse/figures/tigure_upweight_dogfish.png", width = 5, height = 6.5)
# plot_tigure(pm_df_list[[2]], mp_order = x$MP) + pm_angle
# ggsave("mse/figures/tigure_upweight_dogfish_recon_catch.png", width = 5, height = 6.5)
# plot_scenario_projections(mse)
# ggsave("mse/figures/proj.pdf", width = 10, height = 42, limitsize = FALSE)
# .pm_df_list %>% set_names(sc$scenario_human) %>% plot_tigure_facet() + pm_angle
# ggsave("mse/figures/tigure_upweight_dogfish_recon_catch_comparison.png", width = 7.5, height = 6.5)

# Average across OMs
pm_avg <- group_by(pm_df, MP) %>% summarise_if(is.numeric, mean)
pm_min <- group_by(pm_df, MP) %>% summarise_if(is.numeric, min)

satisficed_criteria <- c("LRP 1.5GT" = 0.9, "ST C10" = 0.5)

# Plot tigure averaged across OMs

# plot_tigure(pm_avg, satisficed = satisficed_criteria) + ggtitle("Averaged across reference OMs")
# ggsave("mse/figures/tigure_average_ref.png", width = 6.5, height = 6.5)

# plot_tigure(pm_min, satisficed = satisficed_criteria)

# Plot all tigures in reference set
# walk(names(mse), ~ {
#   g <- plot_tigure(pm_df_list[[.x]]) + ggtitle(sc$scenario_human[match(.x, sc$scenario)])
#   ggsave(paste0("mse/figures/tigure_", .x, ".png"), width = 6.5, height = 6.5)
# })

# All tigures in robustness set
# plot_tigure(pm_df_list_rob[[1]]) + ggtitle("(A) Low M")
# ggsave(paste0("mse/figures/tigure_lowM.png"), width = 6.5, height = 6.5)
#
# plot_tigure(pm_df_list_rob[[2]]) + ggtitle("(B) High CV HBLL")
# ggsave(paste0("mse/figures/tigure_high_CV_HBLL.png"), width = 6.5, height = 6.5)

# # For satisficed MPs
mp_sat <- dplyr::filter(
  pm_avg,
  `LRP 1.5GT` > satisficed_criteria[1],
  `ST C10` > satisficed_criteria[2]) %>%
  pull(MP)
mp_sat <- mp_sat[!mp_sat %in% reference_mp]
mp_sat

mp_not_sat <- pm_avg$MP[!pm_avg$MP %in% mp_sat & !pm_avg$MP %in% reference_mp]

stopifnot(any(!mp_not_sat %in% mp_sat))
stopifnot(any(!mp_sat %in% mp_not_sat))

# # Projections of non-satisficed Index MPs
# walk(names(mse), ~ {
#   g <- plot_main_projections(Sub(mse[[.x]], MPs = c(
#     "Iratio_23", "Iratio_55",
#     "GB_slope_lambda1", "GB_slope_lambda05", "GB_slope_yrs10",
#     "IDX", "IDX_smooth"
#   )),
#   catch_breaks = c(0, 10, 20, 30),
#   catch_ylim = c(0, 40)
#   )
#   ggsave(paste0("mse/figures/projections/projections_index_", .x, ".png"), width = 6.5, height = 6.5)
# })
#
# # Projections of SP MPs (non were satisficed)
# walk(names(mse), ~ {
#   g <- plot_main_projections(Sub(mse[[.x]], MPs = c("SP_8040_10u", "SP_8040_5u", "SP_4010_10u", "SP_4010_5u")),
#     catch_breaks = c(0, 10, 20, 30),
#     catch_ylim = c(0, 40)
#   )
#   ggsave(paste0("mse/figures/projections/projections_SP_", .x, ".png"), width = 6.5, height = 6.5)
# })


# Projections of reference MPs
#walk(names(mse), ~ {
#  g <- plot_main_projections(Sub(mse[[.x]], MPs = reference_mp),
#                             catch_breaks = c(0, 10, 20, 30),
#                             catch_ylim = c(0, 40))
#  ggsave(paste0("mse/figures/projections/projections_refmp_", .x, ".png"), width = 6.5, height = 3.5)
#}
#)
#
## Projections of satisficed MPs
#walk(names(mse), ~ {
#  g <- plot_main_projections(Sub(mse[[.x]], MPs = mp_sat),
#                             catch_breaks = c(0, 10, 20, 30),
#                             catch_ylim = c(0, 40))
#  ggsave(paste0("mse/figures/projections/projections_satisficed_", .x, ".png"), width = 6.5, height = 6.5)
#}

# # Dot plots
ref_mp_cols <- c("grey45", "grey10", "grey75") %>% set_names(reference_mp)
custom_pal <- c(RColorBrewer::brewer.pal(8, "Dark2")[seq_along(mp_sat)], ref_mp_cols) %>%
  set_names(c(mp_sat, reference_mp))
# plot_dots(filter(pm_avg, MP %in% c(mp_sat, reference_mp)), type = "facet", dodge = 0.8)
# ggsave("mse/figures/dot-refset-avg.png", width = 8, height = 3)
#
# # Convergence plot
# scenarios %>%
#   purrr::map(~ DLMtool::Sub(mse[[.x]], MPs = c(mp_sat, reference_mp))) %>%
#   set_names(scenarios_human) %>%
#   gfdlm::plot_convergence(pm = c("LRP 1.5GT", "ST C10"), ylim = c(0.3, 1.05), custom_pal = custom_pal)
# ggsave("mse/figures/convergence.png", height = 4, width = 8)
#

mp_eg_not_sat <- c("Iratio_23", "Iratio_55",
  "GB_slope_lambda1", "GB_slope_lambda05", "GB_slope_yrs10",
  "IDX", "IDX_smooth", "SP_8040_10u", "SP_8040_5u", "SP_4010_10u", "SP_4010_5u")

sc$scenario_human[sc$scenario == "updog_fixsel"] <- "(1) Upweight\ndogfish survey"
sc$scenario_human[sc$scenario == "episodic_recruitment"] <- "(3) Episodic\nrecruitment"
sc$scenario_human[sc$scenario == "upweight_dogfish"] <- "(4) Estimate\nHBLL selectivity"
sc$scenario_human[sc$scenario == "high_index_cv"] <- "(B) High CV HBLL\n(projected)"

g <- map(e_df_list, ~ dplyr::filter(.x, MP %in% union(mp_sat, "NFref"))) %>%
  set_names(sc$scenario_human) %>%
  plot_tigure_facet() +
  theme(
    plot.margin = margin(t = 11/2 - 5, r = 11/2 + 15, b = 11/2, l = 11/2 - 5)
  )
ggsave("mse/figures/ye-tigure-cosewic-all.png", width = 5.7, height = 5)

g <- dplyr::filter(e_avg, MP %in% union(mp_sat, "NFref")) %>%
  plot_tigure() +
  scale_fill_viridis_c(limits = c(0, 1), begin = 0.15, end = 1, alpha = 0.6,
    option = "C", direction = 1)
g
ggsave("mse/figures/ye-tigure-cosewic-avg.png", width = 2.5, height = 3)

plots <- gfdlm::plot_factory(
  mse_list = mse,
  pm = PM,
  scenario_df = sc,
  mp_sat = mp_sat,
  mp_not_sat = mp_not_sat,
  mp_not_sat2 = mp_eg_not_sat,
  mp_ref = reference_mp,
  custom_pal = custom_pal,
  eg_scenario = "updog_fixsel",
  tradeoff = c("LRP 1.5GT", "ST C10"),
  satisficed_criteria = satisficed_criteria,
  skip_projections = FALSE, # TRUE for speed!
  catch_breaks = seq(0, 30, 10),
  catch_ylim = c(0, 40),
  survey_type = "AddInd",
  skip_worms = TRUE # memory problems
)

# rm(mse) # memory problems
g <- plots$projections_index +
 scale_x_continuous(breaks = seq(1975, 2120, 25)) +
 scale_y_continuous(labels = function(x) x / 1e6)
.ggsave("projections-index", width = 12, height = 10, plot = g)

g <- purrr::map(scenarios, ~ DLMtool::Sub(mse[[.x]], MPs = mp_sat)) %>%
 set_names(sc$scenario_human) %>%
 gfdlm::plot_convergence(pm_list = names(satisficed_criteria)) +
 scale_colour_manual(values = custom_pal) +
  theme(legend.position = "bottom")
.ggsave("convergence", width = 9, height = 4.5, plot = g)

.ggsave("dot-refset-avg", width = 8, height = 4,
  plot = plots$dot_refset_avg)
.ggsave("tradeoff-refset-avg", width = 4.5, height = 4,
  plot = plots$tradeoff_avg + coord_equal(xlim = c(0.5, 1), ylim = c(0.5, 1), expand = FALSE))
.ggsave("tradeoff-robset", width = 6, height = 3,
  plot = plots$tradeoff_robset + coord_equal(xlim = c(0.5, 1), ylim = c(0.5, 1), expand = FALSE))

pm_angle <- theme(
  axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0),
  plot.margin = margin(t = 11/2 - 5, r = 11/2 + 10, b = 11/2, l = 11/2 - 5)
)
.ggsave("tigure-refset", width = 6.75, height = 5.5,
  plot = plots$tigure_refset + pm_angle)
.ggsave("tigure-robset", width = 6.5, height = 3,
  plot = plots$tigure_robset + pm_angle)
.ggsave("tigure-refset-min", width = 5, height = 7,
  plot = plots$tigure_refset_min + pm_angle)
.ggsave("tigure-refset-avg", width = 5, height = 7,
  plot = plots$tigure_refset_avg + pm_angle)
.ggsave("radar-refset", width = 10, height = 10,
  plot = plots$radar_refset)
.ggsave("radar-robset", width = 10, height = 5,
  plot = plots$radar_robset)

g <- plots$projections_index +
 scale_x_continuous(breaks = seq(1980, 2090, 20)) +
 coord_cartesian(ylim = c(0, 1e5), expand = FALSE) +
 scale_y_continuous(breaks = seq(0, 1e5, 5e4), labels = function(x) x / 1e6)
.ggsave("projections-index", width = 12, height = 8, plot = g)

walk(names(plots$projections), ~ {
 .ggsave(paste0("projections-", .x),
   width = 8.5, height = 10,
   plot = plots$projections[[.x]]
 )
})
.ggsave("projections-not-sat2", width = 8.5, height = 13,
  plot = plots$projections_not_sat2)
.ggsave("projections-not-sat", width = 6.5, height = 20,
  plot = plots$projections_not_sat)
.ggsave("projections-scenarios", width = 8, height = 11,
  plot = plots$projections_scenarios + scale_color_brewer())

optimize_png <- FALSE
if (optimize_png && !identical(.Platform$OS.type, "windows")) {
  files_per_core <- 2
  setwd("mse/figures")
  system(paste0(
    "find -X . -name '*.png' -print0 | xargs -0 -n ",
    files_per_core, " -P ", parallel::detectCores()/2, " optipng -strip all"
  ))
  setwd(here())
}
