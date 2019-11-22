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
cores <- floor(parallel::detectCores() / 2 - 2)
TMB::openmp(cores)

# Load data and prep. -----------------------------------------------------

d <- readr::read_csv("analysis-survey/dogfish/all_dogfish_sets_data.csv")
joint_grid <- gfplot::dogfish_grid$grid

group_by(d, year) %>% summarise(n_na = sum(!is.na(ye_count)))
group_by(d, year) %>% summarise(n_na = sum(is.na(ye_count)))

# FIXME: For now assume that the NAs should be zeros
d$ye_count[is.na(d$ye_count)] <- 0
d$dogfish_count[is.na(d$dogfish_count)] <- 0

# coast <- gfplot:::load_coastline(range(d$longitude) + c(-.2, .2),
#   range(d$latitude) + c(-0, 0),
#   utm_zone = 9
# )
# coast <- gfplot:::utm2ll(coast)

# ggplot() + geom_polygon(
#   data = coast, aes_string(x = "X", y = "Y", group = "PID"),
#   fill = "grey87", col = "grey70", lwd = 0.2
# ) + coord_fixed() +
#   geom_text(data = filter(d, year == 2004),
#     mapping = aes(longitude, latitude, label = fe_end_deployment_time)) +
#   geom_point(data = joint_grid,
#     mapping = aes(X, Y), col = "red")

# Weird locations, e.g. on land or far away!
d <- filter(d, fe_end_deployment_time != "11/22/2004 8:04")
d <- filter(d, fe_end_deployment_time != "11/20/2004 7:38")
d <- filter(d, fe_end_deployment_time != "11/20/2004 9:02")
d <- filter(d, fe_end_deployment_time != "11/20/2004  14:32:00 PM")

# ggplot() + geom_polygon(
#   data = coast, aes_string(x = "X", y = "Y", group = "PID"),
#   fill = "grey87", col = "grey70", lwd = 0.2
# ) + coord_fixed() +
#   geom_point(data = d,
#     mapping = aes(longitude, latitude, colour = ye_count)) +
#   geom_point(data = joint_grid,
#     mapping = aes(X, Y), col = "red", pch = 21)

d_utm <- convert2utm(d, coords = c("longitude", "latitude"))

d_utm$hook_code <- as.factor(d_utm$hook_code)
d_utm$dogfish_count[d_utm$dogfish_count == 0] <- 1

joint_grid_utm <- joint_grid %>%
  mutate(X = X / 100, Y = Y / 100)

joint_grid_utm$area <- 1
joint_grid_utm$hook_code <- filter(d_utm, year == 2019)$hook_code[1]
joint_grid_utm$dogfish_count <- mean(d_utm$dogfish_count)

nrow(d)
sp <- make_spde(d_utm$X, d_utm$Y, n_knots = 250)
pdf("figs/dogfish-joint-spde.pdf", width = 7, height = 7)
plot_spde(sp)
dev.off()

# Fit models -----------------------------------------------------

d_utm2004 <- filter(d_utm, year == 2004)
nrow(d_utm2004)
sp2004 <- make_spde(d_utm2004$X, d_utm2004$Y, n_knots = 15)
d_utm2004$offset <- log(d_utm2004$area_swept_km2 * 100)

plot_spde(sp2004)
m2004 <- sdmTMB(
  formula = ye_count ~ 1 + offset + as.factor(hook_code),
  data = d_utm2004,
  spde = sp2004,
  silent = FALSE,
  anisotropy = FALSE,
  reml = TRUE,
  spatial_only = T,
  family = nbinom2(link = "log")
)
m2004

m2004 <- MASS::glm.nb(
  formula = ye_count ~ 1 + offset(offset) + as.factor(hook_code) + log(dogfish_count),
  data = d_utm2004
)
m2004
exp(coef(m2004)[["as.factor(hook_code)3"]])
exp(confint(m2004)["as.factor(hook_code)3", ])

ratio <- filter(d_utm, year == 2004) %>% group_by(hook_code) %>% summarize(m = mean(ye_count))
ratio
.ratio <- ratio$m[2] / ratio$m[1]
.ratio

d_utm$area_swept_km2[d_utm$hook_code == 3] <- d_utm$area_swept_km2[d_utm$hook_code == 3] * exp(coef(m2004)[["as.factor(hook_code)3"]])
d_utm$offset <- log(d_utm$area_swept_km2 * 100)
joint_grid_utm$offset <- mean(d_utm$offset)
years <- sort(unique(d_utm$year))
joint_grid_utm <- expand_prediction_grid(joint_grid_utm, years = years)

m <- sdmTMB(
  # formula = ye_count ~ 0 + as.factor(year) + as.factor(hook_code) + log(dogfish_count) + offset,
  formula = ye_count ~ 0 + as.factor(year) + log(dogfish_count) + offset,
  data = d_utm,
  spde = sp,
  time = "year",
  silent = FALSE,
  anisotropy = FALSE,
  reml = TRUE,
  spatial_only = T,
  family = nbinom2(link = "log")
)
m
sink("figs/dogfish-model.txt")
print(m)
sink()
# saveRDS(m, file = "data-generated/dogfish-model.rds")

set.seed(82302)
d_utm$resids <- residuals(m) # randomized quantile residuals

# Project density ------------------------------

get_predictions <- function(model) {
  predict(model,
    newdata = joint_grid_utm,
    return_tmb_object = TRUE, xy_cols = c("X", "Y"), area = joint_grid_utm$area
  )
}

bias_correct <- FALSE
predictions <- get_predictions(m)
ind <- get_index(predictions, bias_correct = bias_correct)
saveRDS(ind, file = "data-generated/dogfish-index.rds")

# Raw data plots -----------------------------------------------------

ggplot(d_utm, aes(X, Y)) +
  geom_tile(
    data = joint_grid_utm, aes(x = X, y = Y), size = 0.5,
    colour = "grey80", fill = "white"
  ) +
  facet_wrap(~year) +
  geom_point(pch = 21, mapping = aes(
    size = ye_count / area_swept_km2
  ), alpha = 1) +
  coord_fixed() +
  scale_color_viridis_c() +
  scale_fill_viridis_c() +
  scale_size_area(max_size = 8) +
  labs(
    colour = "Count density\n(units TODO)", size = "Count density\n(units TODO)",
    fill = "Count density\n(units TODO)"
  )
ggsave("figs/dogfish-ye-raw-data.pdf", width = 10, height = 7)

ggplot(d_utm, aes(X, Y)) +
  geom_tile(
    data = joint_grid_utm, aes(x = X, y = Y), size = 0.5,
    colour = "grey80", fill = "white"
  ) +
  facet_wrap(~year) +
  geom_point(pch = 21, mapping = aes(
    size = dogfish_count / area_swept_km2
  ), alpha = 1) +
  coord_fixed() +
  scale_color_viridis_c() +
  scale_fill_viridis_c() +
  scale_size_area(max_size = 8) +
  labs(
    colour = "Count density\n(units TODO)", size = "Count density\n(units TODO)",
    fill = "Count density\n(units TODO)"
  )
ggsave("figs/dogfish-dogfish-raw-data.pdf", width = 10, height = 7)

# Diagnostics and plots -----------------------------------

diverging_scale <- scale_fill_gradient2(high = scales::muted("red"),
  low = scales::muted("blue"), mid = "grey90")

ggplot(d_utm, aes(X, Y, col = resids)) +
  scale_colour_gradient2(high = scales::muted("red"),
    low = scales::muted("blue"), mid = "grey90") +
  geom_jitter(size = 0.9, width = 0.03, height = 0.03) + facet_wrap(~year) + coord_fixed() +
  labs(colour = "Residual")
ggsave("figs/dogfish-residual-map.pdf", width = 10, height = 10)

qqnorm(d_utm$resids)
qqline(d_utm$resids)

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
      size = ye_count / exp(offset)
    ),
    inherit.aes = FALSE, colour = "grey20", alpha = 0.5
  ) +
  scale_size_area(max_size = 7)
# ggsave("figs/dogfish-prediction-sqrt.pdf", width = 10, height = 10)

g + scale_fill_viridis_c(trans = "log10", option = "D")
ggsave("figs/dogfish-prediction-log.pdf", width = 10, height = 10)

plot_map(predictions$data, exp(est_non_rf)) +
  scale_fill_viridis_c(trans = "sqrt", option = "D") +
  labs(fill = "Fixed effect\nestimate")
ggsave("figs/dogfish-non-rf.pdf", width = 10, height = 10)

plot_map(filter(predictions$data, year == 2019), omega_s, wrap = FALSE) +
  diverging_scale
ggsave("figs/dogfish-omega.pdf", width = 5, height = 5)

# plot_map(predictions$data, epsilon_st) +
#   diverging_scale
# ggsave("figs/hbll-joint-epsilon.pdf", width = 10, height = 10)

ind %>%
  ggplot(aes(year, est)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  xlab("Year") + ylab("Estimated relative abundance")

  # scale_x_continuous(breaks = seq(2004, 2018, 2)) +
  # coord_cartesian(ylim = c(0, max(ind$upr) * 0.65), expand = FALSE,
  #   xlim = range(ind$year) + c(-0.3, 0.3))
ggsave("figs/dogfish-index.pdf", width = 5, height = 4)

