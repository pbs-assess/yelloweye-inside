setwd(here::here()) # for RStudio Jobs
library(dplyr)
library(ggplot2)
library(sdmTMB)
# devtools::install_github("seananderson/ggsidekick")
library(ggsidekick) # for fourth_root_power_trans
theme_set(ggsidekick::theme_sleek())
source("analysis-survey/utils.R")
dir.create("data-generated", showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)

f <- "data-generated/yelloweye-rockfish-inside.rds"
if (file.exists(f)) {
  d <- readRDS(f)
} else {
  d <- gfdata::get_survey_sets("yelloweye rockfish", ssid = c(39, 40))
  block_ids <- gfdata::get_survey_blocks(ssid = c(39, 40))
  d <- dplyr::left_join(d, block_ids)
  d$block_designation <- as.numeric(d$block_designation)
  saveRDS(d, file = f)
}

hook <- readRDS("data-generated/yelloweye_hook_adjusted_data.rds") %>%
  select(
    survey_series_id, year, fishing_event_id, total_hooks, prop_bait_hooks,
    hook_adjust_factor, expected_catch
  )

d <- left_join(d, hook, by = c("survey_series_id", "year", "fishing_event_id"))
d <- d %>%
  select(
    survey_abbrev, year, longitude, latitude, catch_count, hook_count,
    grouping_code, depth_m, block_designation, hook_adjust_factor
  ) %>%
  rename(survey = survey_abbrev, block = block_designation)

d_utm <- convert2utm(d, coords = c("longitude", "latitude"))
d_utm <- filter(d_utm, grouping_code > 270 & grouping_code < 330) # not all part of survey

joint_grid <- readRDS("data-generated/hbll-inside-grid.rds")

d_utm$depth_log <- log(d_utm$depth_m)
d_utm$depth_centred <- d_utm$depth_log - mean(d_utm$depth_log)
d_utm$depth_scaled <- d_utm$depth_centred / sd(d_utm$depth_centred)
d_utm$Y_cent <- d_utm$Y - mean(d_utm$Y)
d_utm$X_cent <- d_utm$X - mean(d_utm$X)
d_utm$area_swept <- d_utm$hook_count * 0.0024384 * 0.009144 * 1000
d_utm$offset_area_swept <- log(d_utm$area_swept)
d_utm$offset <- log(d_utm$area_swept / d_utm$hook_adjust_factor)

joint_grid_utm <- convert2utm(joint_grid, coords = c("longitude", "latitude"))
joint_grid_utm$offset <- mean(d_utm$offset)

years <- sort(unique(d_utm$year))
joint_grid_utm <- expand_prediction_grid(joint_grid_utm, years = years) %>%
  mutate(depth_centred = log(depth) - mean(d_utm$depth_log)) %>%
  mutate(depth_scaled = depth_centred / sd(d_utm$depth_centred))
joint_grid_utm <- mutate(joint_grid_utm, Y_cent = Y - mean(d_utm$Y))
joint_grid_utm <- mutate(joint_grid_utm, X_cent = X - mean(d_utm$X))
north_grid_utm <- filter(joint_grid_utm, survey %in% "HBLL INS N")
south_grid_utm <- filter(joint_grid_utm, survey %in% "HBLL INS S")

g <- ggplot(d_utm, aes(X, Y)) +
  geom_tile(
    data = north_grid_utm, aes(x = X, y = Y), size = 0.5,
    colour = "grey50", fill = "white"
  ) +
  geom_tile(
    data = south_grid_utm, aes(x = X, y = Y), size = 0.5,
    colour = "grey80", fill = "white"
  ) +
  facet_wrap(~year) +
  geom_point(pch = 21, mapping = aes(
    size = catch_count / area_swept,
    colour = catch_count / area_swept
  ), alpha = 1) +
  coord_fixed() +
  scale_color_viridis_c() +
  scale_fill_viridis_c() +
  scale_size_area(max_size = 8) +
  labs(
    colour = "Count density\n(units TODO)", size = "Count density\n(units TODO)",
    fill = "Count density\n(units TODO)"
  )
ggsave("figs/hbll-joint-raw-data.pdf", width = 10, height = 7)

g <- ggplot(filter(joint_grid_utm, year == 2019), aes(X, Y)) +
  geom_tile(aes(x = X, y = Y, fill = area), width = 0.02, height = 0.02) +
  scale_fill_distiller(palette = "BrBG", direction = 1) +
  coord_fixed() +
  labs(fill = expression(Area ~ "in" ~ water ~ (km^2))) +
  xlab("UTMs East (100km)") + ylab("UTMs West (100km)")
ggsave("figs/hbll-area-in-water.pdf", width = 7, height = 5)

sp <- make_spde(d_utm$X, d_utm$Y, n_knots = 500)
pdf("figs/hbll-joint-spde.pdf", width = 7, height = 7)
plot_spde(sp)
dev.off()

# Fit model -----------------------------------------------------

model_file <- "data-generated/hbll-inside-joint-hook-eps-depth.rds"
model_file_depth <- "data-generated/hbll-inside-joint-hook-eps.rds"
if (!file.exists(model_file) || !file.exists(model_file_depth)) {
  tictoc::tic()
  m_nodepth <- sdmTMB(
    formula = catch_count ~ 0 +
      as.factor(year) +
      offset,
    data = d_utm,
    spde = sp,
    time = "year",
    silent = FALSE,
    anisotropy = FALSE,
    cores = 4L,
    reml = TRUE,
    family = nbinom2(link = "log")
  )
  tictoc::toc()
  saveRDS(m_nodepth, file = model_file)

  tictoc::tic()
  m <- sdmTMB(
    formula = catch_count ~ 0 +
      as.factor(year) +
      depth_scaled + I(depth_scaled^2) +
      offset,
    data = d_utm,
    spde = sp,
    time = "year",
    silent = FALSE,
    anisotropy = FALSE,
    cores = 4L,
    reml = TRUE,
    family = nbinom2(link = "log")
  )
  tictoc::toc()
  saveRDS(m, file = model_file_depth)
} else {
  m_nodepth <- readRDS(model_file)
  m <- readRDS(model_file_depth)
  m_nodepth$tmb_obj$retape()
  m$tmb_obj$retape()
}
m_nodepth
m

# AIC(m_nodepth)
# AIC(m)

# Project density ------------------------------

s_years <- filter(d_utm, survey == "HBLL INS S") %>%
  pull(year) %>%
  unique()

predictions_nodepth <- predict(m_nodepth,
  newdata = joint_grid_utm,
  return_tmb_object = TRUE, xy_cols = c("X", "Y"), area = joint_grid_utm$area
)
ind_nodepth <- get_index(predictions_nodepth, bias_correct = FALSE)
saveRDS(ind_nodepth, file = "data-generated/hbll-joint-index-no-depth.rds")

predictions <- predict(m,
  newdata = joint_grid_utm,
  return_tmb_object = TRUE, xy_cols = c("X", "Y"), area = joint_grid_utm$area
)
ind <- get_index(predictions, bias_correct = FALSE)
saveRDS(ind, file = "data-generated/hbll-joint-index.rds")

# Diagnostics -----------------------------------

set.seed(93817)
d_utm$resids <- residuals(m) # randomized quantile residuals
hist(d_utm$resids)
pdf("figs/hbll-joint-residuals-qq.pdf", width = 5, height = 5)
par(cex = 0.75)
qqnorm(d_utm$resids)
abline(a = 0, b = 1)
dev.off()

ggplot(d_utm, aes(X, Y, col = resids)) +
  scale_colour_gradient2() +
  geom_point(size = 0.9) + facet_wrap(~year) + coord_fixed() +
  labs(colour = "Residual")
ggsave("figs/hbll-joint-residual-map.pdf", width = 10, height = 10)

plot_map <- function(dat, column, wrap = TRUE) {
  gg <- ggplot(data = dat) +
    geom_tile(mapping = aes(X, Y, fill = {{ column }}), width = 0.025, height = 0.025) +
    coord_fixed() +
    scale_fill_viridis_c(option = "D")
  if (wrap) gg + facet_wrap(~year) else gg
}

g <- plot_map(predictions$data, exp(est)) +
  scale_fill_viridis_c(trans = "sqrt", option = "D") +
  labs(
    fill = "Estimated\nrelative\nabundance",
    size = "Observed\nrelative\nabundance"
  ) +
  geom_point(
    data = d_utm, pch = 21, mapping = aes(
      x = X, y = Y,
      size = catch_count / area_swept
    ),
    inherit.aes = FALSE, colour = "grey20", alpha = 0.5
  ) +
  scale_size_area(max_size = 7)
ggsave("figs/hbll-joint-prediction-sqrt.pdf", width = 10, height = 10)

g <- g + scale_fill_viridis_c(trans = "log10", option = "D")
ggsave("figs/hbll-joint-prediction-log.pdf", width = 10, height = 10)

plot_map(predictions$data, exp(est_non_rf)) +
  scale_fill_viridis_c(trans = "sqrt", option = "D") +
  labs(fill = "Fixed effect\nestimate")
ggsave("figs/hbll-joint-non-rf.pdf", width = 10, height = 10)

# plot_map(predictions$data, est_rf) +
#   ggtitle("All spatial and spatiotemporal random effects") +
#   scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "grey90")
# ggsave("figs/hbll-joint-rf.pdf", width = 10, height = 10)

plot_map(filter(predictions$data, year == 2018), omega_s, wrap = FALSE) +
  scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "grey90")
ggsave("figs/hbll-joint-omega.pdf", width = 5, height = 5)

plot_map(predictions$data, epsilon_st) +
  scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "grey90")
ggsave("figs/hbll-joint-epsilon.pdf", width = 10, height = 10)

# What about the individual surveys? -----------------------

pred_north <- predict(m,
  newdata = north_grid_utm,
  return_tmb_object = TRUE, xy_cols = c("X", "Y"), area = north_grid_utm$area
)
ind_north <- get_index(pred_north, bias_correct = FALSE)

pred_south <- predict(m,
  newdata = south_grid_utm,
  return_tmb_object = TRUE, xy_cols = c("X", "Y"), area = south_grid_utm$area
)
ind_south <- get_index(pred_south, bias_correct = FALSE)

ind_north$type <- "HBLL INS N"
ind_south$type <- "HBLL INS S"
ind$type <- "HBLL INS"
ind_nodepth$type <- "HBLL INS no depth"

n_years <- filter(d_utm, survey %in% "HBLL INS N") %>%
  pull(year) %>%
  unique()
s_years <- filter(d_utm, survey %in% "HBLL INS S") %>%
  pull(year) %>%
  unique()

ind_north_plot <- filter(ind_north, year %in% n_years)
ind_south_plot <- filter(ind_south, year %in% s_years)

.geomean <- exp(mean(log(ind$est)))
.ratio <- (exp(mean(log(ind_nodepth$est))) / .geomean)
.ind_nodepth <- ind_nodepth %>%
  mutate(est = est / .ratio) %>%
  mutate(upr = upr / .ratio) %>%
  mutate(lwr = lwr / .ratio)

all_plot <- bind_rows(ind_north_plot, ind_south_plot) %>%
  bind_rows(ind) %>%
  mutate(type2 = type) %>%
  bind_rows(mutate(.ind_nodepth, type2 = "HBLL INS")) %>%
  filter(year > 2003)

scale <- 1
g <- bind_rows(ind_north, ind_south) %>%
  bind_rows(ind) %>%
  mutate(type2 = type) %>%
  bind_rows(mutate(.ind_nodepth, type2 = "HBLL INS")) %>%
  filter(year > 2003) %>%
  ggplot(aes(year, est * scale)) +
  geom_line(aes(colour = type)) +
  geom_point(aes(colour = type), pch = 21) +
  geom_point(data = all_plot, aes(colour = type)) +
  geom_ribbon(aes(ymin = lwr * scale, ymax = upr * scale, fill = type),
    alpha = 0.2, colour = NA
  ) +
  facet_wrap(~type2, ncol = 1) +
  labs(colour = "Type", fill = "Type") +
  xlab("Year") + ylab("Estimated relative abundance") +
  geom_vline(xintercept = s_years, lty = 2, alpha = 0.2, lwd = 0.2) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  scale_x_continuous(breaks = seq(2004, 2018, 2)) +
  coord_cartesian(ylim = c(0, max(ind$upr) * 0.8)) +
  theme(legend.position = c(0.22, 0.56))
# g
ggsave("figs/hbll-index-components-eps-depth.pdf", width = 4.1, height = 7.5)

# Design based comparison: -----------------------------------
# out <- d_utm %>%
#   boot_biomass(reps = 1000L)
# all_modelled <- bind_rows(ind_north, ind_south) %>%
#   bind_rows(ind)
# design_based <- all_modelled %>%
#   group_by(type) %>%
#   mutate(est_scaled = est * scale) %>%
#   summarise(geomean = exp(mean(log(est_scaled)))) %>%
#   rename(survey = type) %>%
#   right_join(out, by = "survey") %>%
#   group_by(survey) %>%
#   mutate(scaled_biomass = biomass / (exp(mean(log(biomass))) / geomean)) %>%
#   mutate(scaled_max = upr / (exp(mean(log(biomass))) / geomean)) %>%
#   mutate(scaled_min = lwr / (exp(mean(log(biomass))) / geomean)) %>%
#   rename(type = survey) %>%
#   mutate(type2 = type)
#
# g + geom_line(data = design_based, aes(x = year, y = scaled_biomass),
#   inherit.aes = FALSE) +
#   geom_point(data = design_based, aes(x = year, y = scaled_biomass),
#     inherit.aes = FALSE, pch = 4) +
#   geom_ribbon(data = design_based,
#     aes(x = year, ymin = scaled_min, ymax = scaled_max), inherit.aes = FALSE, alpha = 0.1) +
#   labs(colour = "Type", fill = "Type") + theme(legend.position = "none")
# ggsave("figs/hbll-index-components-with-design-hook-eps-depth.pdf", width = 5.3, height = 8.5)
