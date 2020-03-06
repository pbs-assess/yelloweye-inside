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

# Satisficing rules:

# 80% above LRP within 56 years (1.5 MGT), 95% within 1.5 MGT
# 50% above USR
LRP_thresh <- 0.95

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

catch <- apply(mse[[1]]@CB_hist[1, , (102-8+1):102, ], 2, sum) # Catch since 2012
ref_aadc <- gfdlm:::get_aadc(catch)
`ST AADC` <- gfdlm::pm_factory("AADC", ref_aadc, c(1, 10))
#`AADC 1GT` <- gfdlm::pm_factory("AADC", ref_aadc, c(1, 38))


PM <- c("LRP 1.5GT", "USR 1.5GT", "LRP 1GT", "ST C10", "ST C15", "LT C20", "ST AADC")



# Set up and checks -----------------------------------------------------------
sc <- readRDS("mse/om/ye-scenarios.rds")
sc <- rename(sc, scenario_human = scenarios_human)
#stopifnot(all(reference_mp %in% mp$mp))
fig_dir <- "mse/figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir)

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
#om <- lapply(scenarios, function(x) readRDS(paste0("mse/om/", x, ".rds"))@OM)
#names(om) <- scenarios

mse <- map(scenarios, ~readRDS(paste0("mse/om/MSE_", .x, ".rds")))
names(mse) <- scenarios
for (i in seq_along(mse)) {
  mse[[i]]@OM$RefY <- 20
  mse[[i]]@MPs[c(23, 25)] <- c("IDX", "IDX_smooth")
}


# COSEWIC metric E, Probability that the biomass is at least 2, 5 % B0 within the projection period
COSEWIC_E <- function(MSEobj, Ref = 0.02, Yrs = c(1, 100)) {
  Yrs <- ChkYrs(Yrs, MSEobj)
  PMobj <- new("PMobj")
  PMobj@Name <- PMobj@Caption <- paste0("Probability that biomass is above ", Ref * 100, "% B0")
  PMobj@Stat <- MSEobj@SSB[, , Yrs[1]:Yrs[2]]/MSEobj@OM$SSB0
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
walk(names(mse), ~ {
  g <- plot_tigure(e_df_list[[.x]]) + ggtitle(sc$scenario_human[match(.x, sc$scenario)])
  ggsave(paste0("mse/figures/tigure_cosewic_", .x, ".png"), width = 3.5, height = 6.5)
}
)

# Averaged across reference set
e_avg <- e_df_list[1:4] %>% bind_rows(.id = "scenario") %>% group_by(MP) %>% summarise_if(is.numeric, mean)
pm_df_list <- map(mse[scenarios_ref], ~ gfdlm::get_probs(.x, PM)) # List with all OMs and PMs

g <- plot_tigure(e_avg) + ggtitle("Averaged across reference set")
ggsave("mse/figures/tigure_cosewic_average_ref.png", width = 3.5, height = 6.5)

# Satisficing -----------------------------------------------------------------
pm_df_list <- map(mse[scenarios_ref], ~ gfdlm::get_probs(.x, PM)) # List with all OMs and PMs
pm_df_list_rob <- map(mse[scenarios_rob], ~ gfdlm::get_probs(.x, PM)) # Robustness only
pm_df <- bind_rows(pm_df_list, .id = "scenario") # All as a data.frame

saveRDS(pm_df, file = "mse/om/ye-pm-all.rds")


# Average across OMs
pm_avg <- group_by(pm_df, MP) %>% summarise_if(is.numeric, mean)
pm_min <- group_by(pm_df, MP) %>% summarise_if(is.numeric, min)

satisficed_criteria <- c("LRP 1.5GT" = 0.9, "ST C10" = 0.5)

# Plot tigure averaged across OMs

plot_tigure(pm_avg, satisficed = satisficed_criteria) + ggtitle("Averaged across reference OMs")
ggsave("mse/figures/tigure_average_ref.png", width = 6.5, height = 6.5)

#plot_tigure(pm_min, satisficed = satisficed_criteria)

# Plot all tigures in reference set
walk(names(mse), ~ {
  g <- plot_tigure(pm_df_list[[.x]]) + ggtitle(sc$scenario_human[match(.x, sc$scenario)])
  ggsave(paste0("mse/figures/tigure_", .x, ".png"), width = 6.5, height = 6.5)
}
)

# All tigures in robustness set
plot_tigure(pm_df_list_rob[[1]]) + ggtitle("(A) Low M")
ggsave(paste0("mse/figures/tigure_lowM.png"), width = 6.5, height = 6.5)

plot_tigure(pm_df_list_rob[[2]]) + ggtitle("(B) High CV HBLL")
ggsave(paste0("mse/figures/tigure_high_CV_HBLL.png"), width = 6.5, height = 6.5)


# Subset
#mp_sat <- pm_df_list[[1]]$MP
mp_index <- pm_df_list[[1]]$MP[c(7:22)]
mp_other <- pm_df_list[[1]]$MP[-c(7:22)]

mse_index <- purrr::map(scenarios, ~ DLMtool::Sub(mse[[.x]], MPs = mp_index))
mse_other <- purrr::map(scenarios, ~ DLMtool::Sub(mse[[.x]], MPs = mp_other))

pm_df_list_index <- map(pm_df_list, ~filter(.x, MP %in% mp_index))
pm_df_list_other <- map(pm_df_list, ~filter(.x, MP %in% mp_other))

pm_df_list_index_rob <- map(pm_df_list_rob, ~filter(.x, MP %in% mp_index))
pm_df_list_other_rob <- map(pm_df_list_rob, ~filter(.x, MP %in% mp_other))

#
# For satisficed MPs
mp_sat <- dplyr::filter(pm_avg, `LRP 1.5GT` > satisficed_criteria[1], `ST C10` > satisficed_criteria[2]) %>%
  pull(MP)
mp_sat <- mp_sat[!mp_sat %in% reference_mp]
mp_sat

# Tigure of satisficed MPs
pm_angle <- theme(
  axis.text.x.top = element_text(angle = 45, hjust = 0)
)

plot_tigure(filter(pm_avg, MP %in% mp_sat)) + ggtitle("Averaged acros reference OMs") + pm_angle
ggsave("mse/figures/tigure_satisficed_average_ref.png", width = 4.5, height = 3)

r <- Map(function(x, y, z) plot_tigure(filter(x, MP %in% z)) + ggtitle(y) + pm_angle, x = c(pm_df_list, pm_df_list_rob), y = sc$scenario_human,
         MoreArgs = list(z = mp_sat))

do.call(plot_grid, c(r, list(nrow = 3)))
ggsave("mse/figures/tigure_satisficed.png", width = 8, height = 10)

# Projections of satisficed MPs
walk(names(mse), ~ {
  g <- plot_main_projections(Sub(mse[[.x]], MPs = mp_sat),
                             catch_breaks = c(0, 10, 20, 30),
                             catch_ylim = c(0, 40))
  ggsave(paste0("mse/figures/projections/projections_satisficed_", .x, ".png"), width = 6.5, height = 6.5)
}
)


# Projections of non-satisficed Index MPs
walk(names(mse), ~ {
  g <- plot_main_projections(Sub(mse[[.x]], MPs = c("Iratio_23", "Iratio_55",
                                                    "GB_slope_lambda1", "GB_slope_lambda05", "GB_slope_yrs10",
                                                    "IDX", "IDX_smooth")),
                             catch_breaks = c(0, 10, 20, 30),
                             catch_ylim = c(0, 40))
  ggsave(paste0("mse/figures/projections/projections_index_", .x, ".png"), width = 6.5, height = 6.5)
}
)

# Projections of SP MPs (non were satisficed)
walk(names(mse), ~ {
  g <- plot_main_projections(Sub(mse[[.x]], MPs = c("SP_8040_10u", "SP_8040_5u", "SP_4010_10u", "SP_4010_5u")),
                             catch_breaks = c(0, 10, 20, 30),
                             catch_ylim = c(0, 40))
  ggsave(paste0("mse/figures/projections/projections_SP_", .x, ".png"), width = 6.5, height = 6.5)
}
)

# Dot plots
ref_mp_cols <- c("grey45", "grey10", "grey75") %>% set_names(reference_mp)
custom_pal <- c(RColorBrewer::brewer.pal(8, "Set2")[seq_along(mp_sat)], ref_mp_cols) %>%
  set_names(c(mp_sat, reference_mp))
plot_dots(filter(pm_avg, MP %in% c(mp_sat, reference_mp)), type = "facet", dodge = 0.8)
ggsave("mse/figures/dot-refset-avg.png", width = 8, height = 3)

# Convergence plot
scenarios %>%
  purrr::map(~ DLMtool::Sub(mse[[.x]], MPs = c(mp_sat, reference_mp))) %>%
  set_names(scenarios_human) %>%
  gfdlm::plot_convergence(pm = c("LRP 1.5GT", "ST C10"), ylim = c(0.3, 1.05), custom_pal = custom_pal)
ggsave("mse/figures/convergence.png", height = 4, width = 8)



#plots <- gfdlm::plot_factory(
#  mse_list = mse,
#  pm = PM,
#  scenario_df = sc,
#  mp_sat = mp_sat,
#  mp_not_sat = mp_not_sat,
#  mp_not_sat2 = mp_eg_not_sat,
#  mp_ref = reference_mp,
#  custom_pal = custom_pal,
#  eg_scenario = "updog_fixsel",
#  tradeoff = c("LRP 1.5GT", "ST C10"),
#  satisficed_criteria = satisficed_criteria,
#  skip_projections = FALSE, # TRUE for speed!
#  catch_breaks = seq(0, 75, 25), catch_ylim = c(0, 75),
#  survey_type = "AddInd"
#)
#
#g <- plots$projections_index +
#  scale_x_continuous(breaks = seq(1975, 2120, 25)) +
#  scale_y_continuous(labels = function(x) x / 1e6)
#.ggsave("projections-index", width = 12, height = 10, plot = g)
#
#mp_sat_conv <- dplyr::filter(pm_avg, `LRP 1.5GT` > 0.7, `ST C10` > 0.05) %>%
#  pull(MP)
#mp_sat_conv <- mp_sat_conv[!mp_sat_conv %in% reference_mp]
#mp_sat_conv
#
#g <- purrr::map(scenarios, ~ DLMtool::Sub(mse[[.x]], MPs = mp_sat)) %>%
#  set_names(scenarios_human) %>%
#  gfdlm::plot_convergence(pm_list = names(satisficed_criteria)) +
#  scale_colour_manual(values = custom_pal)
#.ggsave("convergence", width = 8, height = 11, plot = g)
#
#.ggsave("dot-refset-avg", width = 8, height = 4.5, plot = plots$dot_refset_avg)
#
#g <- plots$tradeoff_refset + facet_wrap(~scenario, ncol = 4)
#.ggsave("tradeoff-refset", width = 7.5, height = 5, plot = g)
#.ggsave("tradeoff-robset", width = 6, height = 3, plot = plots$tradeoff_robset)
#
#pm_angle <- theme(
#  axis.text.x.top = element_text(angle = 45, hjust = 0)
#)
#
#.ggsave("tigure-refset", width = 6.75, height = 5.5, plot = plots$tigure_refset + pm_angle)
#.ggsave("tigure-robset", width = 6.5, height = 3, plot = plots$tigure_robset + pm_angle)
#.ggsave("tigure-refset-min", width = 4.5, height = 6.75, plot = plots$tigure_refset_min + pm_angle)
#.ggsave("tigure-refset-avg", width = 4.5, height = 6.75, plot = plots$tigure_refset_avg + pm_angle)
#
#.ggsave("radar-refset", width = 10, height = 10, plot = plots$radar_refset)
#.ggsave("radar-refset-avg", width = 6, height = 6, plot = plots$radar_refset_avg)
#
#g <- plots$projections_index +
#  scale_x_continuous(breaks = seq(1980, 2090, 20)) +
#  coord_cartesian(ylim = c(0, 16e6)) +
#  scale_y_continuous(labels = function(x) x / 1e6)
#.ggsave("projections-index", width = 12, height = 8, plot = g)
#
#walk(names(plots$projections), ~ {
#  .ggsave(paste0("projections-", .x),
#    width = 8.5, height = 10,
#    plot = plots$projections[[.x]]
#  )
#})
#.ggsave("projections-not-sat2", width = 6.5, height = 9, plot = plots$projections_not_sat2)
#.ggsave("projections-not-sat", width = 6.5, height = 20, plot = plots$projections_not_sat)
#.ggsave("projections-scenarios-ref", width = 8, height = 10, plot = plots$projections_scenarios)
#
#
# -----------------------------------------------------------------------------
#
# # Tigure plots ----------------------------------------------------------------
# # All MPs in reference set
# for(i in 1:length(scenarios_ref_human)) {
#   g <- gfdlm::plot_tigure(pm_df_list_index[[i]]) + ggtitle(scenarios_ref_human[i])
#   ggsave(paste0("mse/figures/pm_ref/", scenarios_ref[i], "_index.png"), width = 5.5, height = 6)
#
#   g <- gfdlm::plot_tigure(pm_df_list_other[[i]]) + ggtitle(scenarios_ref_human[i])
#   ggsave(paste0("mse/figures/pm_ref/", scenarios_ref[i], "_other.png"), width = 5.5, height = 5)
# }
#
# # All MPs in robustness set
# for(i in 1:length(scenarios_rob_human)) {
#   g <- gfdlm::plot_tigure(pm_df_list_index_rob[[i]]) + ggtitle(scenarios_rob_human[i])
#   ggsave(paste0("mse/figures/pm_rob/", scenarios_rob[i], "_index.png"), width = 5.5, height = 6)
#
#   g <- gfdlm::plot_tigure(pm_df_list_other_rob[[i]]) + ggtitle(scenarios_rob_human[i])
#   ggsave(paste0("mse/figures/pm_rob/", scenarios_rob[i], "_other.png"), width = 5.5, height = 5)
# }
#
#
#
# # Convergence -----------------------------------------------------------------
#
# walk(names(mse), ~ {
#   g <- gfdlm::plot_convergence(mse[[.x]], "LRP 1.5GT") +
#     #scale_color_brewer(palette = "Set2") +
#     geom_hline(yintercept = 0.95, linetype = 3)
#   ggsave(paste0("mse/figures/convergence/converge-", .x, ".png"), width = 6.5, height = 6.5)
# })
#
#
# # Projections -----------------------------------------------------------------
#
# walk(names(mse_index), ~ {
#   g <- plot_main_projections(mse[[.x]],
#     catch_breaks = c(0, 25, 50, 75, 100),
#     catch_ylim = c(0, 100))
#   ggsave(paste0("mse/figures/projections/projections-index-", .x, ".png"), width = 6.5, height = 6.5)
# }
# )
#
#
# walk(names(mse_other), ~ {
#   g <- plot_main_projections(mse_index[[.x]],
#                              catch_breaks = c(0, 100, 200, 300))
#   ggsave(paste0("mse/figures/projections/projections-other-", .x, ".png"), width = 6.5, height = 6.5)
# }
# )
#
#
#
#
# # Plot future HBLL
# future_hbll <- function(mse, MPs) {
#   MP_ind <- match(MPs, mse[[1]]@MPs)
#   hbll <- lapply(mse, function(x) {
#     lapply(x@Misc$Data[MP_ind], function(y) {
#       apply(y@AddInd[, 1, ], 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
#       })
#   })
#   par(mfrow = c(length(mse), length(MPs)), mar = c(0, 0, 0, 0), oma = c(5, 4, 3, 3))
#   for(i in 1:length(mse)) {
#     for(j in 1:length(MPs)) {
#       matplot(c(1917 + 1:ncol(hbll[[i]][[j]]))[-c(1:82)], t(hbll[[i]][[j]][, -c(1:82)]), typ = 'l',
#               lty = c(2, 1, 1, 2, 2), col = "black", lwd = c(1, 1, 2, 1, 1),
#               yaxt = ifelse(j == 1, 's', 'n'), xaxt = ifelse(i == length(mse), "s", "n"))
#       #if(j == 1) legend("topleft", names(mse)[i], bty = 'n')
#       if(i == 1) {
#         text(2060, 1.2 * max(hbll[[i]][[j]], na.rm = TRUE), MPs[j], xpd = NA, font = 2)
#       }
#       if(j == length(MPs)) {
#         text(x = 2150, y = mean(range(hbll[[i]][[j]], na.rm=TRUE)), names(mse)[i], xpd = NA, srt = -90, font = 2)
#       }
#     }
#   }
#   mtext("Year", side = 1, outer = TRUE, line = 3)
#   mtext("Index", side = 2, outer = TRUE, line = 3)
# }
#
#
# png("mse/figures/future_hbll.png", height = 8, width = 10, units = "in", res = 500)
# future_hbll(mse, mp_index)
# dev.off()
#
#
#
#
#
# names(sc)[2] <- "scenario_human"
# mp_sat_with_ref <- mp_sat
# custom_pal <- structure(gplots::rich.colors(length(mp_sat)), names = mp_sat)
# satisficed_criteria <- structure(rep(0, length(PM)), names = PM)
# plots <- gfdlm::make_typical_plots(mse_list = mse, pm = PM, scenario_df = sc, this_year = this_year,
#                                    mp_sat = mp_sat, mp_not_sat = mp_sat, mp_not_sat_highlight = mp_sat,
#                                    eg_scenario = "base", custom_pal = custom_pal, satisficed_criteria = satisficed_criteria)
#
#
#
#
