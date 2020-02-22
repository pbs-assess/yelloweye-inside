
library(MSEtool)
library(dplyr)
setup(12)

############ Condition operating models with SRA_scope and data
SRA_data <- readRDS("mse/scoping/SRA_data.rds")
data_names <- c("Chist", "Index", "I_sd", "I_type", "length_bin", "s_CAA", "CAA", "CAL", "I_basis")
data_ind <- match(data_names, names(SRA_data))

OM_condition <- readRDS("mse/scoping/OM_250sim.rds")

# Generate error in projected HBLL indices (CV = 0.25)
AddIbeta <- matrix(1, OM_condition@nsim, 5)

set.seed(24)
AddIerr <- rnorm(SRA@OM@nsim * 5 * (SRA@OM@nyears + SRA@OM@proyears), -0.5 * 0.25^2, 0.25) %>% exp() %>%
  array(dim = c(SRA@OM@nsim, 5, SRA@OM@nyears + SRA@OM@proyears))

add_Ierr <- function(SRA) {
  SRA@OM@cpars$AddIbeta <- AddIbeta
  SRA@OM@cpars$AddIerr <- AddIerr
  return(SRA)
}

## SRA for HBLL sel
SRA_for_selectivity <- readRDS("mse/scoping/scoping_base.rds")[[1]]
base_s_vul_par <- c(SRA_for_selectivity@mean_fit$report$s_LFS[1], SRA_for_selectivity@mean_fit$report$s_L5[1])
s_vul_par <- matrix(c(base_s_vul_par, 0.5), 3, 5)
map_s_vul_par <- matrix(NA, 3, 5)

# Upweight dogfish
SRA <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("free", 2),
                 s_selectivity = rep("logistic", 5), cores = 8,
                 vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                 map_s_vul_par = SRA_data$map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                 LWT = list(CAL = 0, CAA = 0, Index = c(1, 4, 1, 1, 1)))
SRA <- add_Ierr(SRA)

saveRDS(SRA, file = "mse/om/upweight_dogfish.rds")
SRA <- readRDS(file = "mse/om/upweight_dogfish.rds")

# Upweight dogfish survey, fix HBLL sel from base
SRA <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("free", 2),
                 s_selectivity = rep("logistic", 5), cores = 8,
                 vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                 s_vul_par = s_vul_par, map_s_vul_par = map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                 LWT = list(CAL = 0, CAA = 0, Index = c(1, 4, 1, 1, 1)))
SRA <- add_Ierr(SRA)

saveRDS(SRA, file = "mse/om/updog_fixsel.rds")
SRA <- readRDS(file = "mse/om/updog_fixsel.rds")

# Low catch - use gfdatabase estimates of commerical catch in 1986-2005
SRA_data2 <- SRA_data
SRA_data2$Chist[match(1986:2005, SRA_data2$Year), 1] <- 0.5 * SRA_data2$Chist[match(1986:2005, SRA_data2$Year), 1]
SRA <- SRA_scope(OM_condition, data = SRA_data2[data_ind], condition = "catch2", selectivity = rep("free", 2),
                 s_selectivity = rep("logistic", 5), cores = 1,
                 vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                 map_s_vul_par = SRA_data$map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                 LWT = list(CAL = 0, CAA = 0, Index = c(1, 4, 1, 1, 1)))
SRA <- add_Ierr(SRA)

saveRDS(SRA, file = "mse/om/lowcatch.rds")
SRA <- readRDS(file = "mse/om/lowcatch.rds")

# Low catch - fix HBLL sel from base
SRA <- SRA_scope(OM_condition, data = SRA_data2[data_ind], condition = "catch2", selectivity = rep("free", 2),
                 s_selectivity = rep("logistic", 5), cores = 8,
                 vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                 s_vul_par = s_vul_par, map_s_vul_par = map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                 LWT = list(CAL = 0, CAA = 0, Index = c(1, 4, 1, 1, 1)))
SRA <- add_Ierr(SRA)

saveRDS(SRA, file = "mse/om/lowcatch_fixsel.rds")
SRA <- readRDS(file = "mse/om/lowcatch_fixsel.rds")


#### Episodic recruitment
SRA <- readRDS("mse/om/updog_fixsel.rds")

set.seed(324)
sporadic_recruitment <- function(x, years = length(x), rate = 1/38, high_sigmaR = 2, plus_one = FALSE) {
  require(dplyr)
  nhigh <- rbinom(1, years, rate)
  high_ind <- sample(1:years, nhigh)

  new_samp <- rnorm(nhigh, -0.5 * high_sigmaR^2, high_sigmaR) %>% exp()
  if(plus_one) new_samp <- new_samp + 1

  x[high_ind] <- new_samp

  return(x)
}

new_Perr_y <- apply(SRA@OM@cpars$Perr_y[, 182:281], 1, sporadic_recruitment) %>% t()

#matplot(t(SRA@OM@cpars$Perr_y[1:5, 182:281]), typ = 'l', ylim = c(0, 5))
#matlines(t(new_Perr_y[1:5, ]), typ = 'o')
#
#matplot(t(SRA@OM@cpars$Perr_y[, 182:281]), typ = 'l')
#matlines(t(new_Perr_y), typ = 'o')
#
#diff <- new_Perr_y/SRA@OM@cpars$Perr_y[, 182:281]
#matplot(t(diff[, 1:56]), typ = 'l')

#plot(apply(new_Perr_y, 2, mean), typ = 'o')
#lines(apply(new_Perr_y, 2, mean), typ = 'o', col = 'red')
#lines(apply(SRA@OM@cpars$Perr_y[, 182:281], 2, mean), typ = 'o', col = 'blue')
#abline(h = 1)

SRA@OM@cpars$Perr_y[, 182:281] <- new_Perr_y
saveRDS(SRA, file = "mse/om/episodic_recruitment.rds")



# Low M robustness scenario
set.seed(91283)
M_samps <- rlnorm(OM_condition@nsim, log(0.025) - 0.5 * 0.2^2, 0.2)
OM_condition@cpars$M <- M_samps

SRA <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("free", 2),
                 s_selectivity = rep("logistic", 5), cores = 1,
                 vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                 map_s_vul_par = SRA_data$map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                 LWT = list(CAL = 0, CAA = 0, Index = c(1, 4, 1, 1, 1)))
SRA <- add_Ierr(SRA)

saveRDS(SRA, file = "mse/om/lowM.rds")
SRA <- readRDS("mse/om/lowM.rds")


# Low M fix sel robustness scenario
set.seed(91283)
M_samps <- rlnorm(OM_condition@nsim, log(0.025) - 0.5 * 0.2^2, 0.2)
OM_condition@cpars$M <- M_samps

SRA <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("free", 2),
                 s_selectivity = rep("logistic", 5), cores = 12,
                 vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                 s_vul_par = s_vul_par, map_s_vul_par = map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                 LWT = list(CAL = 0, CAA = 0, Index = c(1, 4, 1, 1, 1)))
SRA <- add_Ierr(SRA)

saveRDS(SRA, file = "mse/om/lowM_fixsel.rds")
SRA <- readRDS("mse/om/lowM_fixsel.rds")

# Higher index CV
SRA <- readRDS("mse/om/updog_fixsel.rds")

generate_high_Ierr <- function(SRA) {

  Iobs <- SRA@OM@cpars$Data@AddInd[1, , ] %>% t()
  Ipred <- lapply(SRA@Misc, getElement, "Ipred")

  Isd <- rbind %>% do.call(lapply(Ipred, function(x, y) apply(log(y/x), 2, function(xx) sd(xx, na.rm = TRUE)), y = Iobs))
  IAC <- rbind %>% do.call(lapply(Ipred, function(x, y) apply(log(y/x), 2, function(xx) acf(xx[!is.na(xx)], lag.max = 1, plot = FALSE)$acf[2, 1, 1]),
                                  y = Iobs))

  plot(Iobs[, 1], ylim = c(0, 80000), pch = 16, typ = 'o')
  matlines(do.call(cbind, lapply(Ipred, function(x) x[, 1])))

  #hist(Isd[, 1])
  #hist(IAC[, 1])

  set.seed(24)
  I_dev_mu <- -0.5 * Isd^2 * (1 - IAC)/sqrt(1 - IAC^2)
  I_devs <- rnorm(SRA@OM@nsim * 5 * (SRA@OM@proyears + SRA@OM@nyears), I_dev_mu, Isd) %>%
    array(c(SRA@OM@nsim, 5, SRA@OM@proyears + SRA@OM@nyears))
  for(i in 2:dim(I_devs)[3]) I_devs[,,i] <- IAC * I_devs[,,i-1] + I_devs[,,i] * sqrt(1 - IAC^2)

  res <- exp(I_devs)
  #matplot(t(res[,1,]), type = 'l')
  #hist(log(Iobs[, 1]/Ipred[[6]][, 1]))
  #hist(log(res[,1,]))
  SRA@OM@cpars$AddIerr <- res
  return(SRA)
}

SRA <- generate_high_Ierr(SRA)

saveRDS(SRA, file = "mse/om/high_index_cv.rds")



