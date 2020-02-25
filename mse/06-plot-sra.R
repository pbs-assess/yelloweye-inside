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
if (!dir.exists(fig_dir)) dir.create(fig_dir)

# Set up the scenario names ---------------------------------------------------

sc <- tibble::tribble(
  ~scenario,     ~scenarios_human,        ~scenario_type,
  "updog_fixsel",           "(1) Upweight dogfish survey",  "Reference",
  "lowcatch_fixsel",        "(2) Low catch",                "Reference",
  "episodic_recruitment",   "(3) Episodic recruitment",     "Reference",
  "upweight_dogfish",       "(4) Estimate HBLL selectivity","Reference",
  "lowM_fixsel",            "(A) Low M",                    "Robustness",
  "high_index_cv",          "(B) High CV HBLL (projected)", "Robustness"
)
sc <- mutate(sc, order = seq_len(n()))
saveRDS(sc, file = "mse/om/ye-scenarios.rds")

sra_ye <- lapply(sc$scenario, function(x) readRDS(paste0("mse/om/", x, ".rds")))
names(sra_ye) <- sc$scenario


# Some plots ------------------------------------------------------------------

# Depletion, SSB, B/BMSY
get_SSB <- function(x, scenario, mse = NULL, type = c("SSB", "depletion", "MSY")) {
  type <- match.arg(type)

  if(type == "depletion") depletion <- x@SSB / sapply(x@Misc, getElement, "E0_SR")
  if(type == "SSB") depletion <- x@SSB
  if(type == "MSY") {
    if(class(mse) == "MSE") depletion <- x@SSB / mse@Misc$MSYRefs[[1]]$Refs$SSBMSY
    if(class(mse) == "Hist") depletion <- x@SSB / mse@Ref$SSBMSY
  }

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

# Depletion
g <- purrr::map2_df(sra_ye, sc$scenarios_human, get_SSB, type = "depletion") %>%
  mutate(scenario = factor(scenario, levels = sc$scenarios_human)) %>%
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

# SSB
g <- purrr::map2_df(sra_ye, sc$scenarios_human, get_SSB, type = "SSB") %>%
  mutate(scenario = factor(scenario, levels = sc$scenarios_human)) %>%
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

# SSB/SSBMSY
#sra <- readRDS("mse/om/upweight_dogfish.rds")
#Hist <- runMSE(sra@OM, parallel = TRUE, Hist = TRUE)
#saveRDS(Hist, file = 'mse/om/Hist_upweight_dogfish.rds')
mse_ye <- lapply(sc$scenario, function(x) {
  if(x == "upweight_dogfish") {
    readRDS(paste0("mse/om/Hist_", x, ".rds"))
  } else {
    readRDS(paste0("mse/om/MSE_", x, ".rds"))
  }
})
names(mse_ye) <- sc$scenario

g <- do.call(rbind, Map(get_SSB, x = sra_ye, scenario = sc$scenarios_human, mse = mse_ye, type = "MSY")) %>%
  mutate(scenario = factor(scenario, levels = sc$scenarios_human)) %>%
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


# Survey plots
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
surv$scenarios_human <- factor(surv$scenarios_human, levels = sc$scenarios_human)
surv$year <- surv$year + min(all_years) - 1
surv$survey <- factor(surv$survey, levels = survey_names)

I_sd <- sra_ye[[1]]@data$I_sd %>% structure(dimnames = list(all_years, survey_names)) %>%
  as.data.frame() %>% cbind(data.frame(Year = all_years)) %>%
  reshape2::melt(id.vars = c("Year"), variable.name = "survey", value.name = "SD")

Index <- sra_ye[[1]]@data$Index %>% structure(dimnames = list(all_years, survey_names)) %>%
  as.data.frame() %>% cbind(data.frame(Year = all_years)) %>%
  reshape2::melt(id.vars = c("Year"), variable.name = "survey") %>% left_join(I_sd, by = c("Year", "survey"))

Index$lower <- exp(log(Index$value) - 2 * Index$SD * 1.5)
Index$upper <- exp(log(Index$value) + 2 * Index$SD * 1.5)

# Plot all surveys and OMs
g <- ggplot(filter(surv, year >= 1980), aes(year, value, group = paste(iter, survey),
                                            colour = as.character(survey))) +
  geom_line(alpha = 0.05) +
  geom_pointrange(data = Index, mapping = aes(x = Year, y = value, ymin = lower, ymax = upper,
                                              fill = as.character(survey)), inherit.aes = FALSE, pch = 21, colour = "grey40") +
  facet_grid(survey~scenarios_human, scales = "free_y") +
  gfplot::theme_pbs() +
  scale_color_brewer(palette = "Set2", direction = -1) +
  scale_fill_brewer(palette = "Set2", direction = -1) +
  ylab("Index") + xlab("Year") + labs(colour = "Survey", fill = "Survey") + coord_cartesian(xlim = c(1980, 2020))
ggsave("mse/figures/ye-index-fits.png", width = 15, height = 10)

# Plot HBLL
#plot(1:5, col = RColorBrewer::brewer.pal(5, "Set2"), pch = 16)
g <- ggplot(filter(surv, year >= 1980 & survey == "HBLL"), aes(year, value, group = paste(iter))) +
  geom_line(alpha = 0.05, colour = "#66C2A5") +
  geom_pointrange(data = filter(Index, survey == "HBLL"), mapping = aes(x = Year, y = value, ymin = lower, ymax = upper),
                  inherit.aes = FALSE, pch = 21, colour = "grey40", fill = "#66C2A5") +
  facet_wrap(~scenarios_human) +
  gfplot::theme_pbs() +
  scale_color_brewer(palette = "Set2", direction = -1) +
  scale_fill_brewer(palette = "Set2", direction = -1) +
  ylab("Index") + xlab("Year") + coord_cartesian(xlim = c(1980, 2020))
ggsave("mse/figures/ye-index-HBLL.png", width = 8, height = 5)

g <- ggplot(filter(surv, year >= 2000 & survey == "HBLL"), aes(year, value, group = paste(iter))) +
  geom_line(alpha = 0.05, colour = "#66C2A5") +
  geom_pointrange(data = filter(Index, survey == "HBLL"), mapping = aes(x = Year, y = value, ymin = lower, ymax = upper),
                  inherit.aes = FALSE, pch = 21, colour = "grey40", fill = "#66C2A5") +
  facet_wrap(~scenarios_human) +
  gfplot::theme_pbs() +
  scale_color_brewer(palette = "Set2", direction = -1) +
  scale_fill_brewer(palette = "Set2", direction = -1) +
  ylab("Index") + xlab("Year") + coord_cartesian(xlim = c(2000, 2020))
ggsave("mse/figures/ye-index-HBLL2.png", width = 8, height = 5)

# Plot dogfish
g <- ggplot(filter(surv, year >= 1980 & survey == "Dogfish"), aes(year, value, group = paste(iter))) +
  geom_line(alpha = 0.05, colour = "#FC8D62") +
  geom_pointrange(data = filter(Index, survey == "Dogfish"), mapping = aes(x = Year, y = value, ymin = lower, ymax = upper),
                  inherit.aes = FALSE, pch = 21, colour = "grey40", fill = "#FC8D62") +
  facet_wrap(~scenarios_human) +
  gfplot::theme_pbs() +
  scale_color_brewer(palette = "Set2", direction = -1) +
  scale_fill_brewer(palette = "Set2", direction = -1) +
  ylab("Index") + xlab("Year") + coord_cartesian(xlim = c(1980, 2020))
ggsave("mse/figures/ye-index-dogfish.png", width = 8, height = 5)



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
sel$scenarios_human <- factor(sel$scenarios_human, levels = sc$scenarios_human)

g <- sel %>%
  ggplot(aes(Length, value, group = paste(iter))) +
  geom_line(alpha = 0.15) +
  gfplot::theme_pbs() + facet_wrap(~scenarios_human) +
  ylab("HBLL selectivity") + xlab("Age") +
  coord_cartesian(expand = FALSE, ylim = c(-0.01, 1.1))
ggsave("mse/figures/HBLL-selectivity.png", width = 8, height = 5)

# Fishery selectivity
sel2 <- data.frame(Age = 1:80, Commercial = sra_ye[[1]]@Misc[[1]]$vul[1,,1],
                   Recreational = sra_ye[[1]]@Misc[[1]]$vul[1,,2],
                   HBLL = sra_ye[[1]]@Misc[[1]]$s_vul[1,,1]) %>%
  reshape2::melt(id.var = "Age", variable.name = "Fleet", value.name = "Selectivity")

ggplot(sel2, aes(Age, Selectivity, colour = Fleet)) + geom_line() + gfplot::theme_pbs() +
  coord_cartesian(xlim = c(0, 40), expand = FALSE, ylim = c(0, 1.1))
ggsave("mse/figures/fishery-selectivity.png", width = 7, height = 4)


