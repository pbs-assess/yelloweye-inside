library(reshape2)
library(dplyr)
library(gfplot)
library(ggplot2)

dat <- readRDS("data-generated/yelloweye-rockfish-ins-privacy.rds")

# GROWTH
check_convergence_tmb <- TRUE

tmb_init <- list(k = 0.9, linf = 55, log_sigma = log(0.1), t0 = -1)

# Linf = 65.9, K = 0.04, t0 = -8.59, sigma = 0.14
vb_m <- fit_vb(dat$survey_samples, sex = "male", method = "tmb",
               too_high_quantile = 1, check_convergence_tmb = check_convergence_tmb,
               tmb_init = tmb_init)
# Linf = 66.8, K = 0.03, t0 = -10.3, sigma = 0.13
vb_f <- fit_vb(dat$survey_samples, sex = "female", method = "tmb",
               too_high_quantile = 1, check_convergence_tmb = check_convergence_tmb,
               tmb_init = tmb_init)

# Linf = 65.2, K = 0.04, t0 = -9.04, sigma = 0.14
vb_a <- gfplot::fit_vb(dat$survey_samples, sex = "all", method = "tmb",
                       too_high_quantile = 1, check_convergence_tmb = check_convergence_tmb,
                       tmb_init = tmb_init)
plot_growth(object_all = vb_a) + geom_line(data = vb_a$predictions, aes(age, length)) +
  guides(col = FALSE, lty = FALSE)  +
  ggtitle("Growth")
ggsave("mse/figures/vb.png", width = 4, height = 3)


# Plot figure and residual
datpl <- data.frame(age = dat$survey_samples$age, length = dat$survey_samples$length) %>%
  mutate(pred = 65.2 * (1 - exp(-0.04 * (age + 9.04))) * exp(-0.5 * 0.14^2)) %>% mutate(resid = length - pred) %>%
  mutate(resid_log = log(length/pred))
plot(length ~ jitter(age), datpl)
lines(length ~ age, vb_a[[1]], col = "red", lwd = 3)

plot(resid ~ jitter(age), datpl)
abline(h = 0)

plot(resid_log ~ jitter(age), datpl)
abline(h = 0)



lw_m <- fit_length_weight(dat$survey_samples, sex = "male", method = "tmb",
                          too_high_quantile = 1)
lw_f <- fit_length_weight(dat$survey_samples, sex = "female", method = "tmb",
                          too_high_quantile = 1)

# log_a = -11.24, b = 3.09
lw_a <- gfplot::fit_length_weight(dat$survey_samples, sex = "all", method = "tmb", too_high_quantile = 1)
plot_length_weight(object_all = lw_a) + geom_line(data = lw_a$predictions, aes(length, weight)) +
  guides(col = FALSE, lty = FALSE)
ggsave("mse/figures/length-weight.png", width = 4, height = 3)





#### MATURITY
if (sum(!is.na(dat$survey_samples$maturity_code)) > 10) {
  mat_age <- dat$survey_samples %>%
    fit_mat_ogive(
      type = "age",
      months = seq(1, 12))
} else {
  mat_age <- NA
}

ogive_age <- plot_mat_ogive(mat_age)
ggsave("mse/figures/mat-ogive-age.png", width = 5, height = 3)

if (sum(!is.na(dat$survey_samples$maturity_code)) > 10) {
  mat_length <- dat$survey_samples %>%
    fit_mat_ogive(
      type = "length",
      months = seq(1, 12))
} else {
  mat_length <- NA
}

ogive_length <- plot_mat_ogive(mat_length)
ggsave("mse/figures/mat-ogive-length.png", width = 5, height = 3)

cowplot::plot_grid(ogive_age, ogive_length, nrow = 2, ncol = 1)
ggsave("mse/figures/mat-ogives.png", width = 5, height = 6)


# Show predicted vs. observed proportions
prop <- mat_age[[3]]$data %>% group_by(female, age_or_length) %>% summarise(prop_mature = sum(mature)/n()) %>% reshape2::acast(list("female", "age_or_length"))

prop_female <- data.frame(age = as.numeric(colnames(prop)), prop = prop[2, ])

ggplot(prop_female, aes(age, prop)) + geom_point() +
  geom_line(data = filter(mat_age[[2]], female == 1), aes(age_or_length, glmm_fe)) + gfplot::theme_pbs() +
  coord_cartesian(xlim = c(0, 60)) + labs(x = "Age (years)", y = "Probability mature")
ggsave("mse/figures/mat-prop.png", width = 5, height = 3)

# A50 = 14.4, A95 = 27.4, Mat = 0 if age <=7

# Plot maturity by month: ----------------------------------------------------
gfplot::tidy_maturity_months(dat$survey_samples) %>% gfplot::plot_maturity_months()
ggsave("mse/figures/mat-months.png", width = 5, height = 3)

# Age frequencies: ----------------------------------------------------------

ages <- gfplot::tidy_ages_raw(dat$survey_samples,
  survey = c("HBLL INS N", "HBLL INS S"),
  sample_type = "survey")

survey_col_names = c("HBLL INS N", "HBLL INS S")
survey_cols = c(RColorBrewer::brewer.pal(length(survey_col_names), "Set2"))
survey_cols <- stats::setNames(survey_cols, survey_col_names)

g_ages <- gfplot::plot_ages(ages, survey_cols = survey_cols) +
  guides(fill = FALSE, colour = FALSE) +
  ggtitle("Age frequencies") +
  labs(y = "Age (years)")
ggsave("mse/figures/age-freq.png", width = 3, height = 5)

# Length frequencies: ----------------------------------------------------------

len <- gfplot::tidy_lengths_raw(dat$survey_samples,
  sample_type = "survey",
  survey = c("HBLL INS N", "HBLL INS S"))

len$survey_abbrev <- factor(len$survey_abbrev,
  levels = c("HBLL INS N", "HBLL INS S"))

g_lengths <- gfplot::plot_lengths(len, survey_cols = survey_cols,
  bin_size = 2) +
  guides(colour = FALSE, fill = FALSE) +
  ggtitle("Length frequencies") +
  ggplot2::xlab(paste("Length", "(cm)")) +
  ggplot2::ylab("Relative length frequency")
ggsave("mse/figures/length-freq.png", width = 3, height = 5)

optimize_png <- TRUE
if (optimize_png && !identical(.Platform$OS.type, "windows")) {
  files_per_core <- 4
  setwd("mse/figures")
  system(paste0(
    "find -X . -name 'ye*.png' -print0 | xargs -0 -n ",
    files_per_core, " -P ", parallel::detectCores() / 2, " optipng -strip all"
  ))
  setwd("../../")
}
