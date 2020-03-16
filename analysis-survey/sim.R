# Some simulation testing of stitching a survey together
# with missing areas each year

library(sdmTMB)
library(ggplot2)
library(dplyr)

set.seed(29210)

d <- tibble::tibble(X = runif(1000, 0, 1), Y = runif(1000, 0, 1))
d$year <- rep(1:10, each = 100)
d$type <- "a_obs"

nd <- expand.grid(X = seq(0, 1, length.out = 25), Y = seq(0, 1, length.out = 25), year = 1:10)
nd$type <- "b_grid"

.d <- rbind(d, nd)

.rf_sim <- function(model, x, y) {
  out <- sdmTMB:::rf_sim(model, x, y)
  out - mean(out)
}

sigma_O <- 2.2
sigma_E <- 0.3
kappa <- 0.1
x = .d$X
y = .d$Y
rf_omega <- RandomFields::RMmatern(nu = 1, var = sigma_O^2, scale = 1/kappa)
omega_s <- .rf_sim(model = rf_omega, x, y)
.d$omega_s <- omega_s

ggplot(filter(.d, type == "b_grid", year == 1), aes_string("X", "Y", fill = "omega_s")) +
  geom_raster() +
  coord_equal(expand = FALSE) +
  scale_fill_gradient2() +
  gfdlm::theme_pbs()

rf_epsilon <- RandomFields::RMmatern(nu = 1, var = sigma_E^2, scale = 1/kappa)

.d <- .d %>% group_by(year) %>%
  group_split() %>%
  purrr::map_dfr(
    ~tibble(X = .x$X, Y = .x$Y, type = .x$type, omega_s = .x$omega_s, year = .x$year,
      eps_st = .rf_sim(model = rf_epsilon, .x$X, .x$Y))) %>%
  ungroup() %>%
  arrange(type)

ggplot(dplyr::filter(.d, year %in% 1:5),
  aes_string("X", "Y", colour = "eps_st")) +
  geom_point() +
  coord_equal(expand = FALSE) +
  scale_color_gradient2() +
  facet_grid(type~year) +
  gfdlm::theme_pbs()

f <- ~ 0 + as.factor(year)
X_ij <- model.matrix(f, .d)

set.seed(2939)
b_j <- rlnorm(10, meanlog = 0.1, sdlog = 0.2)
b_j
log(b_j)
plot(b_j, type = "o")

.d$eta <- as.numeric(.d$omega_s + .d$eps_st + X_ij %*% b_j)
.d$y <- exp(.d$eta)
# .d$y[1:1000] <- MASS::rnegbin(1000, mu = exp(.d$eta[1:1000]), theta = 1.2) # observations
.d$y[1:1000] <- rpois(1000, lambda = exp(.d$eta[1:1000])) # observations

ggplot(dplyr::filter(.d, year %in% 1:5, type == "a_obs"), aes_string("X", "Y", colour = "y")) +
  geom_point() +
  scale_color_viridis_c(trans = "sqrt") +
  facet_grid(~year) +
  coord_equal(expand = FALSE) +
  gfdlm::theme_pbs()

ggplot(dplyr::filter(.d, year %in% 1:5, type != "a_obs"), aes_string("X", "Y", colour = "y")) +
  geom_point() +
  scale_color_viridis_c(trans = "sqrt") +
  facet_grid(~year) +
  coord_equal(expand = FALSE) +
  gfdlm::theme_pbs()

d <- filter(.d, type == "a_obs")
nd <- filter(.d, type != "a_obs")

d$omega_s_true <- d$omega_s
d$omega_s <- NULL

# throw out data
d <- d[!(d$year %in% seq(1, 9, 2) & d$Y >= 0.5), ]
d <- d[!(d$year %in% seq(2, 10, 2) & d$Y < 0.5), ]
# d %>% group_by(year) %>% count()
# d <- d[-sample(which(d$year %in% seq(2, 10, 2)), 180), ]
# d %>% group_by(year) %>% count()

actual <- group_by(nd, year) %>%
  summarise(total = sum(y))

naive <- group_by(d, year) %>%
  summarize(total = sum(y)) %>%
  mutate(total = total / exp(mean(log(total))) * exp(mean(log(actual$total))))

plot(actual$year, actual$total, type = "o")
lines(naive$year, naive$total)

ggplot(d, aes_string("X", "Y", colour = "y")) +
  geom_point() +
  scale_color_viridis_c(trans = "sqrt") +
  facet_wrap(~year, nrow = 2) +
  coord_equal(expand = FALSE) +
  gfdlm::theme_pbs() +
  labs(colour = "Observed count")
# ggsave(here::here("figs/geostatistical-sim-count.png"), width = 8, height = 3)

sp <- make_spde(d$X, d$Y, n_knots = 200)
plot_spde(sp)

m <- sdmTMB(y ~ 0 + as.factor(year), data = d, spde = sp,
  family = poisson(link = "log"), include_spatial = TRUE,
  silent = FALSE, time = "year", reml = TRUE)
m

p <- predict(m, newdata = nd, return_tmb_object = TRUE, xy_cols = c("X", "Y"))

g1 <- ggplot(p$data, aes_string("X", "Y", fill = "exp(est)")) +
  geom_raster() +
  scale_fill_viridis_c(limits = c(0, exp(max(p$data$est)))) +
  facet_wrap(~year, ncol = 5) +
  geom_point(data = d, inherit.aes = FALSE, aes(X, Y, size = y), pch = 21, alpha = 0.38, colour = "black") +
  scale_size_area() + coord_equal(expand = FALSE) +
  labs(fill = "Expected\ncount", size = "Observed\ncount")  +
scale_x_continuous(breaks = seq(0, 1, .5)) +
  scale_y_continuous(breaks = seq(0, 1, .5)) +
  ggtitle("(B) Observed (dots) and estimated (colour)")

g2 <- ggplot(filter(.d, type == "b_grid"), aes_string("X", "Y", fill = "exp(eta)")) +
  geom_raster() +
  scale_fill_viridis_c(limits = c(0, exp(max(p$data$est)))) +
  facet_wrap(~year, ncol = 5) + coord_equal(expand = FALSE) +
  labs(fill = "Count") +
  scale_x_continuous(breaks = seq(0, 1, .5)) +
  scale_y_continuous(breaks = seq(0, 1, .5)) +
  ggtitle("(A) Actual (simulated)")
g <- cowplot::plot_grid(g2, g1, nrow = 2)
ggsave(here::here("figs/geostatistical-sim-predicted.png"), width = 8.5, height = 8)

ggplot(filter(p$data, year == 1), aes_string("X", "Y", fill = "omega_s")) +
  geom_raster() +
  scale_fill_gradient2()

ggplot(filter(.d, type == "b_grid", year == 1), aes_string("X", "Y", fill = "omega_s")) +
  geom_raster() +
  scale_fill_gradient2()

ggplot(p$data, aes_string("X", "Y", fill = "est_rf")) +
  geom_raster() +
  scale_fill_gradient2() +
  facet_wrap(~year)

ggplot(filter(.d, type == "b_grid"), aes_string("X", "Y", fill = "omega_s + eps_st")) +
  geom_raster() +
  scale_fill_gradient2() +
  facet_wrap(~year)

ggplot(p$data, aes_string("X", "Y", fill = "est_rf")) +
  geom_raster() +
  scale_fill_gradient2() +
  facet_wrap(~year)

.i <- get_index(p, bias_correct = FALSE)

actual <- group_by(nd, year) %>%
  summarise(total = sum(y)) %>%
  mutate(type = "Actual")

naive <- group_by(d, year) %>%
  summarize(total = sum(y)) %>%
  mutate(total = total / exp(mean(log(total))) * exp(mean(log(actual$total)))) %>%
  mutate(type = "Naive design-based")

.i. <- .i %>%
  mutate(total = est) %>%
  mutate(total = total / exp(mean(log(total))) * exp(mean(log(actual$total)))) %>%
  mutate(type = "Geostatistical")

d_ <- bind_rows(actual, naive) %>%
  bind_rows(.i.)

ggplot(d_, aes(year, colour = type, fill = type, lty = type)) +
  geom_vline(xintercept = seq(2, 10, 2), lty = 2, col = "grey60", lwd = 0.3) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), colour = NA, alpha = 0.2) +
  geom_line(aes(y = total)) +
  # geom_line(data = actual, mapping = aes(y = total), col = "red", lwd = 1.2, lty = 1) +
  # geom_line(data = naive, colour = "blue", aes(y = total)) +
  ggsidekick::theme_sleek() +
  ylab("Total abundance") +  scale_x_continuous(breaks = seq(1, 10, 1)) +
  xlab("Year") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(fill = "Type", colour = "Type", lty = "Type") +
  scale_linetype_manual(values = c(1, 2, 6))

ggsave(here::here("figs/geostatistical-sim-stitched-index.png"), width = 6, height = 3.4)

optimize_png <- TRUE
if (optimize_png && !identical(.Platform$OS.type, "windows")) {
  files_per_core <- 2
  setwd("figs")
  system(paste0(
    "find -X . -name '*.png' -print0 | xargs -0 -n ",
    files_per_core, " -P ", parallel::detectCores() / 2, " optipng -strip all"
  ))
  setwd("../")
}
