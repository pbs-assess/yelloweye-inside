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

d <- d %>%
  filter(survey_series_id %in% c(39, 40)) %>% # inside only
  select(
    survey_abbrev, year, longitude, latitude, density_ppkm2,
    grouping_code, depth_m, block_designation
  ) %>%
  rename(survey = survey_abbrev, block = block_designation) %>%
  mutate(density_1000ppkm2 = density_ppkm2 / 1000)

d_utm <- convert2utm(d, coords = c("longitude", "latitude"))
d_utm <- filter(d_utm, !(X < 9.2 & Y < 54)) # not part of the survey... too far south and west

g <- ggplot(d_utm, aes(X, Y,
  size = density_1000ppkm2,
  colour = survey)) +
  facet_wrap(~year) +
  geom_point(pch = 21) +
  scale_size_area()
ggsave("figs/hbll-joint-raw-data.pdf", width = 10, height = 10)

joint_grid <- readRDS("data-generated/hbll-inside-grid.rds")
# d_utm <- left_join(d_utm, select(joint_grid, block, rock100, rock20), by = "block")
# sum(is.na(d_utm$rock100))
# sum(is.na(d_utm$rock20))
# sum(!is.na(d_utm$rock20))

# g <- ggplot(d_utm, aes(X, Y, colour = is.na(rock100))) +i
#   facet_wrap(~year) +
#   geom_point(pch = 21)
# ggsave("figs/hbll-missing-blocks.pdf", width = 10, height = 10)

# d_utm <- filter(d_utm, !is.na(rock100))

d_utm$depth_log <- log(d_utm$depth_m)
d_utm$depth_centred <- d_utm$depth_log - mean(d_utm$depth_log)
d_utm$depth_scaled <- d_utm$depth_centred / sd(d_utm$depth_centred)
d_utm$Y_cent <- d_utm$Y - mean(d_utm$Y)
d_utm$X_cent <- d_utm$X - mean(d_utm$X)

# d_utm$rock20_scaled <- sqrt(d_utm$rock20) / sd(sqrt(d_utm$rock20))

joint_grid_utm <- convert2utm(joint_grid, coords = c("longitude", "latitude"))
years <- sort(unique(d_utm$year))
joint_grid_utm <- expand_prediction_grid(joint_grid_utm, years = years) %>%
  mutate(depth_centred = log(depth) - mean(d_utm$depth_log)) %>%
  mutate(depth_scaled = depth_centred / sd(d_utm$depth_centred))
  # mutate(rock20_scaled = sqrt(rock20) / sd(sqrt(d_utm$rock20)))

joint_grid_utm <- mutate(joint_grid_utm, Y_cent = Y - mean(d_utm$Y))
joint_grid_utm <- mutate(joint_grid_utm, X_cent = X - mean(d_utm$X))
north_grid_utm <- filter(joint_grid_utm, survey %in% "HBLL INS N")
south_grid_utm <- filter(joint_grid_utm, survey %in% "HBLL INS S")

sp <- make_spde(d_utm$X, d_utm$Y, n_knots = 250)
pdf("figs/hbll-joint-spde.pdf", width = 7, height = 7)
plot_spde(sp)
dev.off()

# experiment with offsetting catchability:
# d_utm$north <- ifelse(d_utm$survey == "HBLL INS N", 1, 0)
# joint_grid_utm$north <- ifelse(joint_grid_utm$survey == "HBLL INS N", 1, 0)
# d_utm$real_year <- d_utm$year
# d_utm$year[d_utm$year == 2005] <- 2004
# joint_grid_utm$real_year <- joint_grid_utm$year
# joint_grid_utm <- filter(joint_grid_utm, year != 2005)

# Fit model -----------------------------------------------------

model_file <- "data-generated/hbll-inside-joint.rds"
# if (!file.exists(model_file)) {
  tictoc::tic()
  m <- sdmTMB(
    formula = density_1000ppkm2 ~ 0 +
      Y_cent + I(Y_cent^2) + X_cent +
      # north +
      as.factor(year) + depth_scaled + I(depth_scaled^2),
    data = d_utm,
    spde = sp,
    time = "year",
    silent = FALSE,
    anisotropy = TRUE,
    ar1_fields = FALSE,
    include_spatial = TRUE,
    family = tweedie(link = "log")
  )
  tictoc::toc()
  saveRDS(m, file = model_file)
# } else {
#   m <- readRDS(model_file)
#   m$tmb_obj$retape()
# }
m

# Project density ------------------------------

predictions <- predict(m,
  newdata = joint_grid_utm,
  return_tmb_object = TRUE, xy_cols = c("X", "Y")
)
saveRDS(predictions, file = "data-generated/hbll-inside-predictions.rds")
ind <- get_index(predictions, bias_correct = FALSE)
save(ind, file = "data-generated/hbll-joint-index.rds")

d_utm$resids <- residuals(m) # randomized quantile residuals
hist(d_utm$resids)
pdf("figs/hbll-joint-residuals-qq.pdf", width = 5, height = 5)
par(cex = 0.75)
qqnorm(d_utm$resids)
abline(a = 0, b = 1)
dev.off()

ggplot(d_utm, aes(X, Y, col = resids)) +
  scale_colour_gradient2() +
  geom_point(size = 0.9) + facet_wrap(~year) + coord_fixed()
ggsave("figs/hbll-joint-residual-map.pdf", width = 10, height = 10)

plot_map <- function(dat, column) {
  ggplot() +
    geom_point(data = dat, aes_string("X", "Y", colour = column), size = 0.5) +
    facet_wrap(~year) +
    coord_fixed()
}

plot_map(predictions$data, "exp(est)") +
  scale_colour_viridis_c(option = "C") +
  ggtitle("Prediction (fixed effects + all random effects)")
ggsave("figs/hbll-joint-prediction.pdf", width = 10, height = 10)

plot_map(predictions$data, "exp(est)") +
  scale_colour_viridis_c(trans = "fourth_root_power", option = "C") +
  ggtitle(paste0(
    "Prediction (fourth root transformed colour; ",
    "fixed effects + all random effects)"
  ))
ggsave("figs/hbll-joint-prediction-sqrt.pdf", width = 10, height = 10)

plot_map(predictions$data, "exp(est_non_rf)") +
  ggtitle("Fixed effects only") +
  scale_colour_viridis_c(trans = "fourth_root_power", option = "C")
ggsave("figs/hbll-joint-non-rf.pdf", width = 10, height = 10)

plot_map(predictions$data, "est_rf") +
  ggtitle("All spatial and spatiotemporal random effects") +
  scale_colour_gradient2()
ggsave("figs/hbll-joint-rf.pdf", width = 10, height = 10)

plot_map(filter(predictions$data, year == 2018), "omega_s") +
  ggtitle("Spatial random effects only") +
  scale_colour_gradient2()
ggsave("figs/hbll-joint-omega.pdf", width = 5, height = 5)

plot_map(predictions$data, "epsilon_st") +
  ggtitle("Spatiotemporal random effects only") +
  scale_colour_gradient2()
ggsave("figs/hbll-joint-epsilon.pdf", width = 10, height = 10)

scale <- 2 * 2 # 2 x 2 km grid
ggplot(ind, aes(year, est * scale)) + geom_line() +
  geom_ribbon(aes(ymin = lwr * scale, ymax = upr * scale), alpha = 0.4) +
  xlab("Year") + ylab(expression(Estimated~density~(1000*s~of~fish/km^2))) +
  geom_vline(xintercept = seq(2006, 2018, 2), lty = 2, alpha = 0.2)
ggsave("figs/hbll-index.pdf", width = 8, height = 5)

# What about the individual surveys? -----------------------

pred_north <- predict(m,
  newdata = north_grid_utm,
  return_tmb_object = TRUE, xy_cols = c("X", "Y")
)
ind_north <- get_index(pred_north, bias_correct = FALSE)

pred_south <- predict(m,
  newdata = south_grid_utm,
  return_tmb_object = TRUE, xy_cols = c("X", "Y")
)
ind_south <- get_index(pred_south, bias_correct = FALSE)

ind_north$type <- "HBLL INS N"
ind_south$type <- "HBLL INS S"
ind$type <- "HBLL INS all"

n_years <- filter(d_utm, survey %in% "HBLL INS N") %>%
  pull(year) %>%
  unique()
s_years <- filter(d_utm, survey %in% "HBLL INS S") %>%
  pull(year) %>%
  unique()

ind_north_plot <- filter(ind_north, year %in% n_years)
ind_south_plot <- filter(ind_south, year %in% s_years)
all_plot <- bind_rows(ind_north_plot, ind_south_plot) %>%
  bind_rows(ind)

g <- bind_rows(ind_north, ind_south) %>%
  bind_rows(ind) %>%
  ggplot(aes(year, est * scale)) +
  geom_line(aes(colour = type)) +
  geom_point(aes(colour = type), pch = 21) +
  geom_point(data = all_plot, aes(colour = type)) +
  geom_ribbon(aes(ymin = lwr * scale, ymax = upr * scale, fill = type),
    alpha = 0.2, colour = NA
  ) +
  facet_wrap(~type, ncol = 1) +
  xlab("Year") + ylab(expression(Estimated~density~(1000*s~of~fish/km^2))) +
  geom_vline(xintercept = seq(2004, 2018, 2), lty = 2, alpha = 0.2, lwd = 0.2) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  scale_x_continuous(breaks = seq(2004, 2018, 2))
g
ggsave("figs/hbll-index-components.pdf", width = 5.5, height = 8.5)

# Design based comparison: -----------------------------------
out <- boot_biomass(d, reps = 500L)
all_modelled <- bind_rows(ind_north, ind_south) %>%
  bind_rows(ind)
design_based <- all_modelled %>%
  group_by(type) %>%
  mutate(est_scaled = est * scale) %>%
  summarise(geomean = exp(mean(log(est_scaled)))) %>%
  rename(survey = type) %>%
  right_join(out, by = "survey") %>%
  group_by(survey) %>%
  mutate(scaled_biomass = biomass / (exp(mean(log(biomass))) / geomean)) %>%
  mutate(scaled_max = upr / (exp(mean(log(biomass))) / geomean)) %>%
  mutate(scaled_min = lwr / (exp(mean(log(biomass))) / geomean)) %>%
  rename(type = survey)

g + geom_line(data = design_based, aes(x = year, y = scaled_biomass), inherit.aes = FALSE) +
  geom_point(data = design_based, aes(x = year, y = scaled_biomass), inherit.aes = FALSE, pch = 4) +
  geom_ribbon(data = design_based, aes(x = year, ymin = scaled_min, ymax = scaled_max), inherit.aes = FALSE, alpha = 0.1)
ggsave("figs/hbll-index-components-with-design.pdf", width = 5.5, height = 8.5)
