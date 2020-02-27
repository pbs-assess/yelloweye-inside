library(reshape2)
library(dplyr)
library(gfplot)

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
plot_growth(object_all = vb_a) + geom_line(data = vb_a$predictions, aes(age, length))
ggsave("mse/figures/ye-length.png", width = 4, height = 4)


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
plot_length_weight(object_all = lw_a) + geom_line(data = lw_a$predictions, aes(length, weight))
ggsave("mse/figures/ye-weight.png", width = 4, height = 4)





#### MATURITY
if (sum(!is.na(dat$survey_samples$maturity_code)) > 10) {
  mat_age <- dat$survey_samples %>%
    fit_mat_ogive(
      type = "age",
      months = seq(1, 12))
} else {
  mat_age <- NA
}

plot_mat_ogive(mat_age)
ggsave("mse/figures/maturity2.png", width = 5, height = 3)


# Show predicted vs. observed proportions
prop <- mat_age[[3]]$data %>% group_by(female, age_or_length) %>% summarise(prop_mature = sum(mature)/n()) %>% acast(list("female", "age_or_length"))

prop_female <- data.frame(age = as.numeric(colnames(prop)), prop = prop[2, ])

ggplot(prop_female, aes(age, prop)) + geom_point() +
  geom_line(data = filter(mat_age[[2]], female == 1), aes(age_or_length, glmm_fe)) + gfplot::theme_pbs() +
  coord_cartesian(xlim = c(0, 60)) + labs(x = "Age (years)", y = "Probability mature")
ggsave("mse/figures/ye-maturity3.png", width = 5, height = 3)

# A50 = 14.4, A95 = 27.4, Mat = 0 if age <=7
