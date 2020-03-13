library("dplyr")
library("DLMtool")
library("MSEtool")
library("here")
library("purrr")
library("cowplot")
library("ggplot2")

species_name <- "Inside Yelloweye Rockfish"
starting_year <- 1918
ending_year <- 2019
all_years <- seq(starting_year, ending_year)
nyear <- length(all_years)

fig_dir <- "mse/figures"
#if (!dir.exists(fig_dir)) dir.create(fig_dir)

# Set up the scenario names ---------------------------------------------------

sc <- tibble::tribble(
  ~scenario,     ~scenario_human,        ~scenario_type,
  "updog_fixsel",           "(1) Upweight\ndogfish survey",  "Reference",
  "lowcatch_fixsel",        "(2) Low catch",                "Reference",
  "episodic_recruitment",   "(3) Episodic\nrecruitment",     "Reference",
  "upweight_dogfish",       "(4) Estimate\nHBLL selectivity","Reference",
  "lowM_fixsel",            "(A) Low M",                    "Robustness",
  "high_index_cv",          "(B) High CV HBLL\nprojected", "Robustness"
)
sc <- mutate(sc, order = seq_len(n()))
saveRDS(sc, file = "mse/om/ye-scenarios.rds")

sra_ye <- lapply(sc$scenario, function(x) readRDS(paste0("mse/om/", x, ".rds")))
names(sra_ye) <- sc$scenario

#get converged replicates

#Check this is working ... copied from Rex. All scenarios have the same proportion 0.38

scenarios <- sc$scenario %>% purrr::set_names(sc$scenario_human)
oms <- map(scenarios, ~ {
  readRDS(paste0("mse/om/", .x, ".rds"))@OM
})
yelloweye_converged <- map_dfr(oms, ~tibble(nsim = .x@nsim), .id = "scenario")
saveRDS(yelloweye_converged, file = here("mse/om/ye-converged.rds"))


# Some plots ------------------------------------------------------------------

# Depletion, SSB, B/BMSY
get_SSB <- function(x, scenario, mse = NULL, type = c("SSB", "depletion", "MSY")) {
  type <- match.arg(type)

  if(type == "depletion") depletion <- x@SSB / sapply(x@Misc, getElement, "E0_SR")
  if(type == "SSB") depletion <- x@SSB
  if(type == "MSY") depletion <- x@SSB / mse@Misc$MSYRefs$ByYear$SSBMSY[, 1, mse@nyears]

  #if(class(mse) == "Hist") depletion <- x@SSB / mse@Ref$SSBMSY

  d1 <- t(apply(depletion[, -nyear], 2,
    FUN = quantile,
    probs = c(0.025, 0.5, 0.975)
  )) %>%
    as.data.frame() %>%
    cbind(all_years) %>%
    mutate(scenario = scenario) %>%
    rename(lwr = 1, med = 2, upr = 3, year = all_years)
  d2 <- t(apply(depletion[, -nyear], 2, FUN = quantile, probs = c(0.25, 0.75))) %>%
    as.data.frame() %>%
    cbind(all_years) %>%
    rename(lwr50 = 1, upr50 = 2, year = all_years)

  left_join(d1, d2, by = "year")
}

get_F <- function(x, scenario) {

  .F1 <- map(x@Misc, "F_at_age")
  .F <- map_dfc(.F1, ~tibble(.F = apply(.x, 1, max)))
  .F <- t(as.matrix(.F))

  last_year <- dim(.F)[2]
  all_years <- seq(x@OM@CurrentYr - x@OM@nyears + 1, x@OM@CurrentYr)
  all_years <- all_years #[-length(all_years)]

  d1 <- t(apply(.F[, ], 2,
                FUN = quantile,
                probs = c(0.025, 0.5, 0.975)
  )) %>%
    as.data.frame() %>%
    cbind(all_years) %>%
    mutate(scenario = scenario) %>%
    rename(lwr = 1, med = 2, upr = 3, year = all_years)
  d2 <- t(apply(.F[, ], 2,
                FUN = quantile,
                probs = c(0.25, 0.75))) %>%
    as.data.frame() %>%
    cbind(all_years) %>%
    rename(lwr50 = 1, upr50 = 2, year = all_years)

  left_join(d1, d2, by = "year")
}

get_Perr_y <- function(x, scenario) {
  max_age <- x@OM@maxage
  nyears <- x@OM@nyears
  perr_y <- x@OM@cpars$Perr_y[,max_age:(max_age+nyears-1), drop=FALSE]
  all_years <- seq(x@OM@CurrentYr - x@OM@nyears + 1, x@OM@CurrentYr)
  reshape2::melt(perr_y) %>%
    rename(iteration = Var1) %>%
    mutate(year = rep(all_years, each = max(iteration))) %>%
    mutate(scenario = scenario)
}


# Depletion -------------------------------------------------------------------
g <- purrr::map2_df(sra_ye, sc$scenario_human, get_SSB, type = "depletion") %>%
  mutate(scenario = factor(scenario, levels = sc$scenario_human)) %>%
  ggplot(aes(year, med, ymin = lwr, ymax = upr)) +
  geom_ribbon(fill = "grey90") +
  geom_ribbon(fill = "grey70", mapping = aes(ymin = lwr50, ymax = upr50)) +
  geom_line() +
  facet_wrap(vars(scenario)) +
  gfplot::theme_pbs() +
  labs(x = "Year", y = "Depletion") + coord_cartesian(expand = FALSE, ylim = c(0, 1.1))
ggsave(file.path(fig_dir, paste0("ye-compare-SRA-depletion-panel.png")),
  width = 8, height = 4
)

# SSB -------------------------------------------------------------------------
g <- purrr::map2_df(sra_ye, sc$scenario_human, get_SSB, type = "SSB") %>%
  mutate(scenario = factor(scenario, levels = sc$scenario_human)) %>%
  ggplot(aes(year, med, ymin = lwr, ymax = upr)) +
  geom_ribbon(fill = "grey90") +
  geom_ribbon(fill = "grey70", mapping = aes(ymin = lwr50, ymax = upr50)) +
  geom_line() +
  facet_wrap(vars(scenario)) +
  gfplot::theme_pbs() +
  labs(x = "Year", y = "Spawning biomass") + coord_cartesian(expand = FALSE, ylim = c(0, 8e3))
ggsave(file.path(fig_dir, paste0("ye-compare-SRA-SSB-panel.png")),
       width = 8, height = 4
)

#F
g <- purrr::map2_df(sra_ye, sc$scenario_human, get_F) %>%
  mutate(scenario = factor(scenario, levels = sc$scenario_human)) %>%
  ggplot(aes(year, med, ymin = lwr, ymax = upr)) +
  geom_ribbon(fill = "grey90") +
  geom_ribbon(fill = "grey70", mapping = aes(ymin = lwr50, ymax = upr50)) +
  geom_line() +
  facet_wrap(vars(scenario)) +
  gfplot::theme_pbs() +
  labs(x = "Year", y = "F") +
  coord_cartesian(ylim = c(0, 0.4), expand = FALSE)
ggsave(here::here("mse/figures/ye-compare-SRA-F-panel.png"),
       width = 8, height = 6.75
)

#Recdevs
g <- purrr::map2_df(sra_ye, sc$scenario_human, get_Perr_y) %>%
  mutate(scenario = factor(scenario, levels = sc$scenario_human)) %>%
  dplyr::filter(iteration %in% 1:100) %>%
  ggplot(aes(year, y = log(value), group = iteration)) +
  geom_line(alpha = 0.1) +
  facet_wrap(vars(scenario)) +
  gfplot::theme_pbs() +
  labs(x = "Year", y = "Recruitment deviation in log space") +
  coord_cartesian(ylim = c(-1.5, 1.7), expand = FALSE) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.6)
# g
ggsave(here::here("mse/figures/ye-compare-SRA-recdev-panel.png"),
       width = 8, height = 6.75
)

# SSB/SSBMSY  -----------------------------------------------------------------
#sra <- readRDS("mse/om/upweight_dogfish.rds")
#Hist <- runMSE(sra@OM, parallel = TRUE, Hist = TRUE)
#saveRDS(Hist, file = 'mse/om/Hist_upweight_dogfish.rds')
mse_ye <- lapply(sc$scenario, function(x) readRDS(paste0("mse/om/MSE_", x, ".rds")))
names(mse_ye) <- sc$scenario

g <- do.call(rbind, Map(get_SSB, x = sra_ye, scenario = sc$scenario_human, mse = mse_ye, type = "MSY")) %>%
  mutate(scenario = factor(scenario, levels = sc$scenario_human)) %>%
  ggplot(aes(year, med, ymin = lwr, ymax = upr)) +
  geom_ribbon(fill = "grey90") +
  geom_ribbon(fill = "grey70", mapping = aes(ymin = lwr50, ymax = upr50)) +
  geom_line() +
  facet_wrap(vars(scenario)) +
  gfplot::theme_pbs() +
  labs(x = "Year", y = expression(SSB/SSB[MSY])) + coord_cartesian(expand = FALSE, ylim = c(0, 5)) +
  geom_hline(yintercept = c(0.4, 0.8), linetype = 3)
ggsave(file.path(fig_dir, paste0("ye-compare-SRA-MSY-panel.png")),
       width = 8, height = 4
)

# Survey plots  ---------------------------------------------------------------
survey_names = c("HBLL", "Dogfish", "CPUE 86-90", "CPUE 95-01", "CPUE 03-05")
get_sra_survey <- function(sra, sc_name, survey_names) {
  n_surv <- dim(sra@Misc[[1]]$Ipred)[2]
  out2 <- purrr::map(seq_len(n_surv), function(i) {
    surveys <- do.call(cbind, purrr::map(sra@Misc, ~ .$Ipred[,i,drop=FALSE]))
    out <- reshape2::melt(surveys) %>%
      rename(year = Var1, iter = Var2)
    out$year <- out$year
    out$scenario <- sc_name
    out$survey <- survey_names[i]
    out
  })
  bind_rows(out2)
}
#surv <- get_sra_survey(sra, "pinniped")
#surv <- left_join(surv, sc, by = "scenario")

surv <- purrr::map2_dfr(sra_ye, sc$scenario, get_sra_survey, survey_names = survey_names)
surv <- left_join(surv, sc, by = "scenario")
surv$scenario_human <- factor(surv$scenario_human, levels = sc$scenario_human)
surv$year <- surv$year + min(all_years) - 1
surv$survey <- factor(surv$survey, levels = survey_names)

I_sd <- sra_ye[[1]]@data$I_sd %>% structure(dimnames = list(all_years, survey_names)) %>%
  as.data.frame() %>% cbind(data.frame(Year = all_years)) %>%
  reshape2::melt(id.vars = c("Year"), variable.name = "survey", value.name = "SD")

Index <- sra_ye[[1]]@data$Index %>% structure(dimnames = list(all_years, survey_names)) %>%
  as.data.frame() %>% cbind(data.frame(Year = all_years)) %>%
  reshape2::melt(id.vars = c("Year"), variable.name = "survey") %>% left_join(I_sd, by = c("Year", "survey"))

Index$lower <- exp(log(Index$value) - 2 * Index$SD)
Index$upper <- exp(log(Index$value) + 2 * Index$SD)

# Plot all surveys and OMs  ---------------------------------------------------
g <- ggplot(filter(surv, year >= 1980), aes(year, value, group = paste(iter, survey),
                                            colour = as.character(survey))) +
  geom_line(alpha = 0.05) +
  geom_pointrange(data = Index, mapping = aes(x = Year, y = value, ymin = lower, ymax = upper,
                                              fill = as.character(survey)), inherit.aes = FALSE, pch = 21, colour = "grey40") +
  facet_grid(survey~scenario_human, scales = "free_y") +
  gfplot::theme_pbs() +
  scale_color_brewer(palette = "Set2", direction = -1) +
  scale_fill_brewer(palette = "Set2", direction = -1) +
  ylab("Index") + xlab("Year") + labs(colour = "Survey", fill = "Survey") + coord_cartesian(xlim = c(1980, 2020))
ggsave("mse/figures/ye-index-fits.png", width = 15, height = 10)

# Plot HBLL  ------------------------------------------------------------------
#plot(1:5, col = RColorBrewer::brewer.pal(5, "Set2"), pch = 16)
g <- ggplot(filter(surv, year >= 1980 & survey == "HBLL"), aes(year, value, group = paste(iter))) +
  geom_line(alpha = 0.05, colour = "#66C2A5") +
  geom_pointrange(data = filter(Index, survey == "HBLL"), mapping = aes(x = Year, y = value, ymin = lower, ymax = upper),
                  inherit.aes = FALSE, pch = 21, colour = "grey40", fill = "#66C2A5") +
  facet_wrap(~scenario_human) +
  gfplot::theme_pbs() +
  scale_color_brewer(palette = "Set2", direction = -1) +
  scale_fill_brewer(palette = "Set2", direction = -1) +
  ylab("Index") + xlab("Year") + coord_cartesian(xlim = c(1980, 2020))
ggsave("mse/figures/ye-index-HBLL.png", width = 8, height = 5)

g <- ggplot(filter(surv, year >= 2000 & survey == "HBLL"), aes(year, value, group = paste(iter))) +
  geom_line(alpha = 0.05, colour = "#66C2A5") +
  geom_pointrange(data = filter(Index, survey == "HBLL"), mapping = aes(x = Year, y = value, ymin = lower, ymax = upper),
                  inherit.aes = FALSE, pch = 21, colour = "grey40", fill = "#66C2A5") +
  facet_wrap(~scenario_human) +
  gfplot::theme_pbs() +
  scale_color_brewer(palette = "Set2", direction = -1) +
  scale_fill_brewer(palette = "Set2", direction = -1) +
  ylab("Index") + xlab("Year") + coord_cartesian(xlim = c(2000, 2020))
ggsave("mse/figures/ye-index-HBLL2.png", width = 8, height = 5)

# Plot dogfish  ---------------------------------------------------------------
g <- ggplot(filter(surv, year >= 1980 & survey == "Dogfish"), aes(year, value, group = paste(iter))) +
  geom_line(alpha = 0.05, colour = "#FC8D62") +
  geom_pointrange(data = filter(Index, survey == "Dogfish"), mapping = aes(x = Year, y = value, ymin = lower, ymax = upper),
                  inherit.aes = FALSE, pch = 21, colour = "grey40", fill = "#FC8D62") +
  facet_wrap(~scenario_human) +
  gfplot::theme_pbs() +
  scale_color_brewer(palette = "Set2", direction = -1) +
  scale_fill_brewer(palette = "Set2", direction = -1) +
  ylab("Index") + xlab("Year") + coord_cartesian(xlim = c(1980, 2020))
ggsave("mse/figures/ye-index-dogfish.png", width = 8, height = 5)

# Plot HBLL age comps, report N = sampling trips for the multinomial likelihood in the SRA -----------
SRA_data <- readRDS("mse/scoping/SRA_data.rds")
Ntrips <- rowSums(SRA_data$s_CAA[,,1], na.rm = TRUE)
yr_ind <- Ntrips > 0
yr_plot <- SRA_data$Year[yr_ind]
Ntrips <- data.frame(Year = yr_plot, N = paste("N =", Ntrips[Ntrips > 0]))

HBLL_obs <- structure(SRA_data$s_CAA[yr_ind,,1]/rowSums(SRA_data$s_CAA[yr_ind,,1]), dimnames = list(Year = yr_plot, Age = 1:80)) %>%
  reshape2::melt(value.name = "Frequency")

nsim <- 250
HBLL_pred <- lapply(sc$scenario, function(xx) {
  do.call(rbind, lapply(1:nsim, function(x) {
    y <- sra_ye[[xx]]@Misc[[x]]$s_CAApred[match(yr_plot, SRA_data$Year),,1]
    y <- structure(y/rowSums(y), dimnames = list(Year = yr_plot, Age = 1:80))
    y <- reshape2::melt(y, value.name = "Frequency")
    y$Iter <- x
    y
  }))
})

dir.create("mse/figures/conditioning", showWarnings = FALSE)
walk(1:length(HBLL_pred), ~{
  g <- ggplot(HBLL_pred[[.x]], aes(Age, Frequency, group = Iter)) + facet_wrap(~Year, scales = "free_y") +
    geom_line(alpha = 0.05, colour = "#66C2A5") +
    geom_line(data = HBLL_obs, mapping = aes(x = Age, y = Frequency), inherit.aes = FALSE) +
    geom_point(data = HBLL_obs, mapping = aes(x = Age, y = Frequency), inherit.aes = FALSE, pch = 21, colour = "grey40", fill = "#66C2A5") +
    geom_label(data = Ntrips, mapping = aes(label = N), x = Inf, y = Inf, hjust = "right", vjust = "top", inherit.aes = FALSE) +
    gfplot::theme_pbs() + ggtitle(sc$scenario_human[.x])
  ggsave(paste0("mse/figures/conditioning/HBLL_age_comp_", sc$scenario[.x], ".png"), width = 10, height = 8, dpi = 500)
})

# Selectivity HBLL and dogfish
get_sra_selectivity <- function(sc_name) {
  sra <- sra_ye[[sc_name]]
  x <- do.call(cbind, purrr::map(sra@Misc, ~ .$s_vul[101,,1]))
  out <- reshape2::melt(x) %>%
    rename(Length = Var1, iter = Var2)
  out$scenario <- sc_name
  out
}
sel <- map_dfr(sc$scenario, get_sra_selectivity) # pick one
sel <- left_join(sel, sc, by = "scenario")
sel$scenario_human <- factor(sel$scenario_human, levels = sc$scenario_human)

g <- sel %>%
  ggplot(aes(Length, value, group = paste(iter))) +
  geom_line(alpha = 0.15) +
  gfplot::theme_pbs() + facet_wrap(~scenario_human) +
  ylab("HBLL selectivity") + xlab("Age") +
  coord_cartesian(expand = FALSE, ylim = c(-0.01, 1.1))
ggsave("mse/figures/HBLL-selectivity.png", width = 8, height = 5)

# Fishery selectivity  --------------------------------------------------------
sel2 <- data.frame(Age = 1:80, Commercial = sra_ye[[1]]@Misc[[1]]$vul[1,,1],
                   Recreational = sra_ye[[1]]@Misc[[1]]$vul[1,,2],
                   HBLL = sra_ye[[1]]@Misc[[1]]$s_vul[1,,1]) %>%
  reshape2::melt(id.var = "Age", variable.name = "Fleet", value.name = "Selectivity")

ggplot(sel2, aes(Age, Selectivity, colour = Fleet)) + geom_line() + gfplot::theme_pbs() +
  coord_cartesian(xlim = c(0, 40), expand = FALSE, ylim = c(0, 1.1))
ggsave("mse/figures/fishery-selectivity.png", width = 7, height = 4)

ggplot(filter(sel2, Fleet != "HBLL"), aes(Age, Selectivity, colour = Fleet)) + geom_line() + gfplot::theme_pbs() +
  coord_cartesian(xlim = c(0, 40), expand = FALSE, ylim = c(0, 1.1))
ggsave("mse/figures/fishery-selectivity2.png", width = 7, height = 4)

# Histograms of M and h  ------------------------------------------------------

#histogram of M and h
samps <- data.frame(M = sra_ye[[1]]@OM@cpars$M, h = sra_ye[[1]]@OM@cpars$h, lowM = sra_ye$lowM_fixsel@OM@cpars$M)
#png("mse/figures/M.png", res = 400, units = "in", height = 4, width = 5)
#hist(samps$M, xlab = "M", main = "")
#abline(v = 0.045, lty = 3, lwd = 2)
#dev.off()

#png("mse/figures/h.png", res = 400, units = "in", height = 4, width = 5)
#hist(samps$h, xlab = "steepness", main = "")
#abline(v = 0.71, lty = 3, lwd = 2)
#dev.off()

ggplot(samps, aes(M)) + geom_histogram(bins = 20) + geom_vline(xintercept = 0.045, linetype = 2) +
  gfplot::theme_pbs() + labs(x = "Natural mortality", ylab = "Frequency")
ggsave("mse/figures/M.png", height = 4, width = 5)

ggplot(samps, aes(lowM)) + geom_histogram(bins = 20) + geom_vline(xintercept = 0.045, linetype = 2) +
  gfplot::theme_pbs() + labs(x = "Natural mortality", ylab = "Frequency")
ggsave("mse/figures/lowM.png", height = 4, width = 5)

ggplot(samps, aes(h)) + geom_histogram(bins = 20) + geom_vline(xintercept = 0.71, linetype = 2) +
  gfplot::theme_pbs() + labs(x = "Steepness", ylab = "Frequency")
ggsave("mse/figures/steepness.png", height = 4, width = 5)

samps <- data.frame(M = c(samps$M, samps$lowM), Scenario = rep(c("All others", "Low M"), each = 250))
ggplot(samps, aes(M, colour = Scenario)) + geom_freqpoly(bins = 20) +
  gfplot::theme_pbs() + labs(x = "Natural mortality", ylab = "Frequency")
ggsave("mse/figures/lowM.png", height = 4, width = 5)

#### Low/high catch  ----------------------------------------------------------
cat <- data.frame(Year = rep(1918:2019, 3),
                  Catch = c(sra_ye[[1]]@data$Chist[, 2], sra_ye[[1]]@data$Chist[, 1], sra_ye$lowcatch_fixsel@data$Chist[, 1]),
                  Fleet = c(rep("Recreational", 102), rep("Catch", 2 * 102)),
                  Scenario = c(rep("All others", 2 * 102), rep("Low catch", 102)))

ggplot(cat, aes(Year, Catch, colour = Fleet, linetype = Scenario)) + geom_line() + gfplot::theme_pbs()
ggsave("mse/figures/catch.png", width = 5.5, height = 3.5, dpi = 500)

ggplot(filter(cat, Scenario != "Low catch"), aes(Year, Catch, colour = Fleet)) + geom_line() + gfplot::theme_pbs()
ggsave("mse/figures/catch2.png", width = 5.5, height = 3.5, dpi = 500)

### COSEWIC indicators and probability below LRP/USR in 2019
COSEWIC_Bdecline_hist <- function(MSEobj, Ref = 0.7, Yr = NULL) {

  # Historic
  stopifnot(!is.null(Yr))
  if(length(Yr) > 1) stop("Length of Yr is one.")

  year_start <- MSEobj@nyears - abs(Yr) + 1
  if(year_start < 1) year_start <- 1
  SSB <- apply(MSEobj@SSB_hist, c(1, 3), sum)
  metric <- 1 - SSB[, MSEobj@nyears]/SSB[, year_start]
  out <- sum(metric < Ref)/length(metric)
  return(out)
}

P_LRP <- function(MSEobj, LRP = 0.4, Yr = NULL) {
  if(is.null(Yr)) Yr <- MSEobj@nyears
  SSB <- apply(MSEobj@SSB_hist, c(1, 3), sum)
  metric <- SSB[, Yr]/MSEobj@Misc$MSYRefs$ByYear$SSBMSY[, 1, Yr]
  out <- sum(metric >= LRP)/length(metric)
  return(out)
}

P_USR <- P_LRP
formals(P_USR)$LRP <- 0.8

COSEWIC_P70 <- COSEWIC_P50 <- COSEWIC_P30 <- COSEWIC_Bdecline_hist
formals(COSEWIC_P50)$Ref <- 0.5
formals(COSEWIC_P30)$Ref <- 0.3

# Historical indicators  ------------------------------------------------------
# P70 - probability that SSB has not declined at least 70% within the past 3 GT
# P50 - probability that SSB has not declined at least 50% within the past 3 GT
# P30 - probability that SSB has not declined at least 30% within the past 3 GT

# LRP - probability that SSB > 0.4 BMSY in 2019
# USR - probability that SSB > 0.8 BMSY in 2019
mse_ye <- map(sc$scenario, ~readRDS(paste0("mse/om/MSE_", .x, ".rds")))
YE_historical <- rbind(vapply(mse_ye, COSEWIC_P70, numeric(1), Yr = 3 * 38),
                       vapply(mse_ye, COSEWIC_P50, numeric(1), Yr = 3 * 38),
                       vapply(mse_ye, COSEWIC_P30, numeric(1), Yr = 3 * 38),
                       vapply(mse_ye, P_LRP, numeric(1)),
                       vapply(mse_ye, P_USR, numeric(1))) %>%
  structure(dimnames = list(c("P70", "P50", "P30", "LRP", "USR"), sc$scenario_human)) %>%
  t() %>% as.data.frame()

YE_historical$MP <- factor(sc$scenario_human, levels = sc$scenario_human)

g <- gfdlm::plot_tigure(YE_historical, mp_order = rev(sc$scenario_human))
ggsave("mse/figures/historical_indicators.png", height = 4, width = 5)

# Plot episodic recruitment ---------------------------------------------------
png("mse/figures/rec_dev.png", height = 5, width = 5, units = "in", res = 400)
par(mfrow = c(2, 1), mar = c(2, 3, 1, 1), oma = c(3, 2, 0, 0))
matplot(1918:2119, t(sra_ye[[1]]@OM@cpars$Perr_y[26:28, -c(1:79)]), typ = 'l', lty = 1, xlab = "", ylab = "")
abline(v = 2019, lty = 2)
legend("topleft", "All other scenarios", bty = "n")

matplot(1918:2119, t(sra_ye$episodic_recruitment@OM@cpars$Perr_y[26:28, -c(1:79)]), typ = 'l', lty = 1, xlab = "", ylab = "")
abline(v = 2019, lty = 2)
mtext("Year", side = 1, outer = TRUE)
mtext("Recruitment deviations (normal space)", side = 2, outer = TRUE)
legend("topleft", "Episodic recruitment", bty = "n")
dev.off()

# Plot higher Isd and AC of HBLL in (B) High index CV  ------------------------
SRA <- readRDS("mse/om/high_index_cv.rds")

Iobs <- SRA@OM@cpars$Data@AddInd[1, 1, ]
Ipred <- lapply(SRA@Misc, getElement, "Ipred")

Isd <- vapply(Ipred, function(x, y) sd(log(y/x[, 1]), na.rm = TRUE), numeric(1), y = Iobs)
IAC <- vapply(Ipred, function(x, y) {xx <- log(y/x[, 1]); acf(xx[!is.na(xx)], lag.max = 1, plot = FALSE)$acf[2, 1, 1]},
              numeric(1), y = Iobs)

high_error <- data.frame(Isd = Isd, IAC = IAC)

g1 <- ggplot(high_error, aes(Isd)) + geom_histogram(bins = 15) + gfplot::theme_pbs() + labs(x = "HBLL standard deviation")
g2 <- ggplot(high_error, aes(IAC)) + geom_histogram(bins = 15) + gfplot::theme_pbs() + labs(x = "HBLL autocorrelation")

cowplot::plot_grid(g1, g2)
ggsave("mse/figures/HBLL_high_CV.png", width = 5, height = 3)

sc2 <- readRDS(here("mse", "om", "ye-scenarios.rds"))
# sc2$scenario_human <- paste0(sc2$order, " - ", sc2$scenario_human)
x <- oms %>% set_names(sc2$scenario_human) %>%
  map_dfr(~tibble(
    D = .x@cpars$D,
    h = .x@cpars$h,
    R0 = .x@cpars$R0,
    sigma_R = .x@cpars$Perr,
    AC = .x@cpars$AC,
    #L50 = .x@cpars$L50,
    #L50_95 = .x@cpars$L50_95,
    t0 = .x@cpars$t0,
    k = .x@cpars$K,
    Linf = .x@cpars$Linf,
    M = .x@cpars$M_ageArray[,1,1],
  ), .id = "Scenario") %>%
  reshape2::melt(id.vars = "Scenario") %>%
  dplyr::filter(!(variable == "R0" & value > 1e7))

x %>% dplyr::filter(variable %in% c("R0", "AC", "D")) %>%
  ggplot(aes(value)) +
  geom_histogram(bins = 30, colour = "grey40", fill = "white", lwd = 0.4) +
  # geom_freqpoly(aes(colour = Scenario), bins = 30) +
  facet_grid(Scenario~variable, scales = "free_x")+
  gfdlm::theme_pbs() +
  coord_cartesian(ylim = c(0, 200), expand = FALSE) +
  xlab("Parameter value") + ylab("Count") +
  # theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()) +
  scale_colour_brewer(palette = "Dark2")

ggsave(here::here("mse/figures/ye-sra-estimated.png"),
       width = 6.5, height = 8.5)

x %>% dplyr::filter(variable %in% c("sigma_R", "h", "L50", "L50_95", "t0", "k", "Linf", "M")) %>%
  ggplot(aes(value)) +
  # geom_histogram(bins = 40, colour = "grey60") +
  geom_freqpoly(aes(colour = Scenario), bins = 30) +
  facet_wrap(~variable, scales = "free_x")+
  gfdlm::theme_pbs() +
  coord_cartesian(ylim = c(0, 300), expand = FALSE) +
  xlab("Parameter value") + ylab("Count") +
  # theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()) +
  scale_colour_brewer(palette = "Dark2")

ggsave(here::here("mse/figures/ye-sra-filtered.png"),
       width = 6.5, height = 5.5)
