library(dplyr)
library(ggplot2)
library(sdmTMB)
# devtools::install_github("seananderson/ggsidekick")
library(ggsidekick) # for fourth_root_power_trans
theme_set(ggsidekick::theme_sleek())
source("analysis-survey/utils.R")
dir.create("data-generated", showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)
cores_star <- parallel::detectCores()
cores <- if (cores_star == 16L) 4 else 2
TMB::openmp(cores)

# Load data and prep. -----------------------------------------------------

d <- readr::read_csv("dogfish/all_dogfish_sets_data.csv")
joint_grid <- gfplot::dogfish_grid$grid

group_by(d, year) %>% summarise(n_na = sum(!is.na(ye_count)))
group_by(d, year) %>% summarise(n_na = sum(is.na(ye_count)))

# FIXME: Assume that the NAs should be zeros
d$ye_count[is.na(d$ye_count)] <- 0
d$dogfish_count[is.na(d$dogfish_count)] <- 0
d <- d %>% arrange(year, fishing_event_id)
d$circle_hook <- ifelse(d$hook_code == 1, 1, 0)
d$circle_hook_centered <- d$circle_hook - mean(d$circle_hook)
d$hook_code <- as.factor(d$hook_code)
d$dogfish_count[d$dogfish_count == 0] <- 1

d_utm <- convert2utm(d, coords = c("longitude", "latitude"))
joint_grid_utm <- joint_grid %>%
  mutate(X = X / (1000*100), Y = Y / (1000*100))

joint_grid_utm$area <- gfplot::dogfish_grid$cell_area # .5x.5 km grid
# joint_grid_utm$hook_code <- filter(d_utm, year == 1986)$hook_code[1]
# joint_grid_utm$hook_code <- filter(d_utm, year == 2019)$hook_code[1]
joint_grid_utm$circle_hook_centered <- 0 # mean(d_utm$circle_hook_centered)
joint_grid_utm$dogfish_count <- mean(d_utm$dogfish_count)

d_utm$offset <- log(d_utm$area_swept_km2 * 100)
# d_utm$offset <- log(d_utm$lglsp_hook_count / 100)
joint_grid_utm$offset <- mean(d_utm$offset)
years <- sort(unique(d_utm$year))
joint_grid_utm <- expand_prediction_grid(joint_grid_utm, years = years)

# Fit models -----------------------------------------------------
length(unique(paste(d_utm$X, d_utm$Y)))
sp <- make_spde(d_utm$X, d_utm$Y, n_knots = 300)
png("figs/dogfish-joint-spde.png", width = 7, height = 7,
  units = "in", res = 200)
plot_spde(sp)
dev.off()

m <- sdmTMB(
  formula =
    ye_count ~ 0 + circle_hook_centered +
    as.factor(year) +
    log(dogfish_count) + offset,
  data = d_utm,
  spde = sp,
  time = "year",
  silent = FALSE,
  anisotropy = TRUE,
  reml = TRUE,
  family = nbinom2(link = "log")
)
print(m)
sink("figs/dogfish-model.txt")
print(m)
sink()
saveRDS(m, file = "data-generated/dogfish-model.rds")

set.seed(82302)
d_utm$resids <- residuals(m) # randomized quantile residuals

# Project density ------------------------------

get_predictions <- function(model) {
  predict(model,
    newdata = joint_grid_utm,
    return_tmb_object = TRUE, xy_cols = c("X", "Y"),
    area = joint_grid_utm$area
  )
}
predictions <- get_predictions(m)
ind <- get_index(predictions, bias_correct = FALSE)
saveRDS(ind, file = "data-generated/dogfish-index.rds")

ind %>%
  ggplot(aes(year, est)) +
  geom_line(alpha = 0.5) +
  # geom_point() +
  geom_pointrange(aes(ymin = lwr, ymax = upr)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.08) +
  xlab("Year") + ylab("Estimated relative abundance") +
  scale_x_continuous(breaks = seq(1985, 2020, 5)) +
  coord_cartesian(
    ylim = c(0, max(ind$upr) * 0.85), expand = FALSE,
    xlim = range(ind$year) + c(-0.3, 0.3)
  ) +
  geom_vline(xintercept = 2004, lty = 2, alpha = 0.3)
ggsave("figs/dogfish-index-estimated-hook.png", width = 6.5, height = 4)

# Raw data plots -----------------------------------------------------

coast <- gfplot:::load_coastline(range(d$longitude),
  range(d$latitude),
  utm_zone = 9
)
coast$X <- coast$X/100
coast$Y <- coast$Y/100
ggcoast <-  geom_polygon(
  data = coast, aes_string(x = "X", y = "Y", group = "PID"),
  fill = "grey92", col = "grey80", lwd = 0.2
)
# ggplot() + ggcoast + coord_fixed()
coords <- coord_fixed(xlim = range(joint_grid_utm$X),
  ylim = range(joint_grid_utm$Y))

plot_raw_data <- function(column) {
  lab <- "Count density\n(units TODO)"
  ggplot(d_utm, aes(X, Y)) +
    coords +
    ggcoast +
    geom_tile(
      data = joint_grid_utm, aes(x = X, y = Y), width = 1/100, height = 1/100,
      colour = "grey70", fill = "grey70"
    ) +
    facet_wrap(~year) +
    geom_point(pch = 21, mapping = aes_string(size = column)) +
    scale_color_viridis_c() +
    scale_fill_viridis_c() +
    scale_size_area(max_size = 8) +
    labs(colour = lab, size = lab, fill = lab)
}

g <- plot_raw_data("ye_count / area_swept_km2")
ggsave("figs/dogfish-yelloweye-per-area-data.png", width = 10, height = 7)
g <- plot_raw_data("ye_count")
ggsave("figs/dogfish-yelloweye-raw-data.png", width = 10, height = 7)

g <- plot_raw_data("dogfish_count / area_swept_km2")
ggsave("figs/dogfish-dogfish-per-area-data.png", width = 10, height = 7)
g <- plot_raw_data("dogfish_count")
ggsave("figs/dogfish-dogfish-raw-data.png", width = 10, height = 7)

g <- plot_raw_data("area_swept_km2")
ggsave("figs/dogfish-area_swept_km2-raw-data.png", width = 10, height = 7)

g <- plot_raw_data("lglsp_hook_count")
ggsave("figs/dogfish-lglsp_hook_count-raw-data.png", width = 10, height = 7)

# Diagnostics and plots -----------------------------------

diverging_scale <- scale_fill_gradient2(
  high = scales::muted("red"),
  low = scales::muted("blue"), mid = "grey90"
)

g <- ggplot(d_utm, aes(X, Y, col = resids)) +
  scale_colour_gradient2(
    high = scales::muted("red"),
    low = scales::muted("blue"), mid = "grey90"
  ) +
  geom_jitter(size = 0.9, width = 0.03, height = 0.03) + facet_wrap(~year) +
  coords +
  ggcoast +
  labs(colour = "Residual")
ggsave("figs/dogfish-residual-map.png", width = 10, height = 10)

qqnorm(d_utm$resids)
qqline(d_utm$resids)

plot_map <- function(dat, column, wrap = TRUE) {
  gg <- ggplot(data = dat) +
    geom_tile(mapping = aes(X, Y, fill = {{ column }}),
      width = 0.01, height = 0.01) +
    coords +
    ggcoast +
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
g <- g + scale_fill_viridis_c(trans = "log10", option = "D")
ggsave("figs/dogfish-prediction-log.png", width = 10, height = 10)

g <- plot_map(predictions$data, exp(est_non_rf)) +
  scale_fill_viridis_c(trans = "sqrt", option = "D") +
  labs(fill = "Fixed effect\nestimate")
ggsave("figs/dogfish-non-rf.png", width = 10, height = 10)

g <- plot_map(filter(predictions$data, year == 2019), omega_s, wrap = FALSE) +
  diverging_scale
ggsave("figs/dogfish-omega.png", width = 5, height = 5)

g <- plot_map(predictions$data, epsilon_st) +
  diverging_scale
ggsave("figs/dogfish-epsilon.png", width = 10, height = 10)
