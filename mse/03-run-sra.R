
library(MSEtool)

############ Condition operating models with SRA_scope and data
SRA_data <- readRDS("mse/scoping/SRA_data.rds")
data_names <- c("Chist", "Index", "I_sd", "I_type", "length_bin", "s_CAA", "CAA", "CAL", "I_basis")
data_ind <- match(data_names, names(SRA_data))

OM_condition <- readRDS("mse/scoping/OM2.rds")


# Base
SRA <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("free", 2),
                 s_selectivity = rep("logistic", 5), cores = 1,
                 vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                 map_s_vul_par = SRA_data$map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                 LWT = list(CAL = 0, CAA = 0))
saveRDS(SRA, file = "mse/OM/base.rds")
#SRA <- readRDS("mse/OM/base.rds")

plot(SRA, file = "mse/OM/OM_base", dir = getwd(), open_file = FALSE, f_name = SRA_data$f_name, s_name = SRA_data$s_name,
     MSY_ref = c(0.4, 0.8))

# Upweight dogfish
SRA2 <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("free", 2),
                 s_selectivity = rep("logistic", 5), cores = 1,
                 vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                 map_s_vul_par = SRA_data$map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                 LWT = list(CAL = 0, CAA = 0, Index = c(1, 4, 1, 1, 1)))

saveRDS(SRA, file = "mse/OM/upweight_dogfish.rds")
#SRA <- readRDS("mse/OM/upweight_dogfish.rds")
plot(SRA, file = "mse/OM/OM_upweight_dogfish", dir = getwd(), open_file = FALSE, f_name = SRA_data$f_name, s_name = SRA_data$s_name,
     MSY_ref = c(0.4, 0.8))

#UPweight dogfish survey, fix HBLL sel

s_vul_par <- matrix(c(46.27726, 30.55939, 0.5), 3, 5)
map_s_vul_par <- matrix(NA, 3, 5)

SRA <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("free", 2),
                 s_selectivity = rep("logistic", 5), cores = 1,
                 vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                 s_vul_par = s_vul_par, map_s_vul_par = map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                 LWT = list(CAL = 0, CAA = 0, Index = c(1, 4, 1, 1, 1)))

compare_SRA(SRA, SRA2, scenario = list(names = c("fix HBLL sel", "est HBLL sel")))

# Low catch - use gfdatabase estimates of commerical catch in 1986-2005
SRA_data2 <- SRA_data
SRA_data2$Chist[match(1986:2005, SRA_data2$Year), 1] <- 0.5 * SRA_data2$Chist[match(1986:2005, SRA_data2$Year), 1]
SRA <- SRA_scope(OM_condition, data = SRA_data2[data_ind], condition = "catch2", selectivity = rep("free", 2),
                 s_selectivity = rep("logistic", 5), cores = 1,
                 vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                 map_s_vul_par = SRA_data$map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                 LWT = list(CAL = 0, CAA = 0, Index = c(1, 4, 1, 1, 1)))

# Low catch - fix HBLL sel
SRA <- SRA_scope(OM_condition, data = SRA_data2[data_ind], condition = "catch2", selectivity = rep("free", 2),
                 s_selectivity = rep("logistic", 5), cores = 1,
                 vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                 s_vul_par = s_vul_par, map_s_vul_par = map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                 LWT = list(CAL = 0, CAA = 0, Index = c(1, 4, 1, 1, 1)))


saveRDS(SRA, file = "mse/OM/low_catch.rds")
#SRA <- readRDS("mse/OM/OM_low_catch.rds")
plot(SRA, file = "mse/OM/OM_low_catch", dir = getwd(), open_file = FALSE, f_name = SRA_data$f_name, s_name = SRA_data$s_name,
     MSY_ref = c(0.4, 0.8))

# No HBLL age comps

s_vul_par <- matrix(c(46.27726, 30.55939, 0.5), 3, 5)
map_s_vul_par <- matrix(NA, 3, 5)

SRA <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("free", 2),
                 s_selectivity = rep("logistic", 5), cores = 1,
                 vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                 s_vul_par = s_vul_par, map_s_vul_par = map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                 LWT = list(CAL = 0, CAA = 0, s_CAA = 0, Index = c(1, 4, 1, 1, 1)))



# Est rec sel
SRA_data$vul_par[1:2, 2] <- c(30, 25)
map_vul_par <- matrix(NA, 80, 2)
map_vul_par[1:2, 2] <- 1:2
SRA <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = c("free", "logistic"),
                 s_selectivity = rep("logistic", 5), cores = 1,
                 vul_par = SRA_data$vul_par, map_vul_par = map_vul_par,
                 map_s_vul_par = SRA_data$map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                 LWT = list(CAL = 1, CAA = 0, Index = c(1, 4, 1, 1, 1)))

# Est both rec and com sel
SRA_data$vul_par[1:2, ] <- c(50, 40, 30, 25)
map_vul_par <- matrix(NA, 80, 2)
map_vul_par[1:2, ] <- 1:4
SRA <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("logistic", 2),
                 s_selectivity = rep("logistic", 5), cores = 1,
                 vul_par = SRA_data$vul_par, map_vul_par = map_vul_par,
                 map_s_vul_par = SRA_data$map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                 LWT = list(CAL = 1, CAA = 20, Index = c(1, 4, 1, 1, 1)))






#### Episodic recruitment
SRA <- readRDS("mse/OM/upweight_dogfish.rds")
set.seed(324)

sporadic_recruitment2 <- function(x, years = length(x), low_sigmaR = 0.4, high_sigmaR = 0.8) {
  require(dplyr)
  nhigh <- 25

  high_ind <- sample(1:years, nhigh)
  new_samp <- rnorm(nhigh, -0.5 * high_sigmaR^2, high_sigmaR) %>% exp()

  x[high_ind] <- new_samp

  return(x)
}

new_Perr_y <- apply(SRA@OM@cpars$Perr_y[, 182:281], 1, sporadic_recruitment2)
SRA@OM@cpars$Perr_y[, 182:281] <- t(new_Perr_y)
saveRDS(SRA, file = "mse/OM/sporadic_recruitment.rds")

# M = 0.02
OM_condition@cpars$M <- rep(0.02, OM@nsim)
SRA <- SRA_scope(OM_condition, condition = "catch2", Chist = SRA_data$Chist, Index = SRA_data$Index, I_sd = SRA_data$I_sd, I_type = SRA_data$I_type,
                 selectivity = rep("logistic", 2), s_selectivity = rep("logistic", 5), length_bin = 0.1 * SRA_data$length_bin, cores = 12,
                 s_CAA = SRA_data$s_CAA, vul_par = SRA_data$vul_par, map_s_vul_par = SRA_data$map_s_vul_par,
                 map_log_rec_dev = SRA_data$map_log_rec_dev)
saveRDS(SRA, file = "mse/OM/lowM.rds")
SRA <- readRDS("mse/OM/lowM.rds")
plot(SRA, file = "mse/OM/OM_lowM", dir = getwd(), open_file = FALSE, f_name = SRA_data$f_name, s_name = SRA_data$s_name,
     MSY_ref = c(0.4, 0.8))


# Pinniped mortality
Mp <- readxl::read_excel("mse/scoping/pinniped_M.xlsx") + 0.02 # spreadsheet includes estimate of pinniped M from 2011 assessment
Mp$mu <- log(Mp$Median)

find_sigma <- function(Mp) {
  uniroot_fn <- function(sd, mu, upper_q) {
    mean_ <- exp(mu + 0.5 * sd *sd)
    mean_log <- log(mean_)
    plnorm(upper_q, mean_log, sd) - 0.975
  }
  out <- mapply(function(x, y) uniroot(uniroot_fn, interval = c(1e-8, 10), mu = x, upper_q = y)$root,
                x = Mp$mu, y = Mp$Upper)
  return(out)
}
Mp$sigma <- find_sigma(Mp)

set.seed(205)
Msamps_p <- runif(OM@nsim)
out_fn <- function(x, y) qlnorm(x, -0.5 * y^2, y)
M_dev <- outer(Msamps_p, Mp$sigma, out_fn) %>% t() %>% "*"(Mp$Median)
matplot(M_dev, type = 'l')
apply(M_dev, 1, mean)

Msamps_extra <- M_dev[nrow(M_dev), ] %>% matrix(OM@proyears + 10, OM@nsim, byrow = TRUE)

OM@cpars$M <- NULL
OM@cpars$M_ageArray <- rbind(M_dev, Msamps_extra) %>% array(c(OM@nyears + OM@proyears, OM@nsim, OM@maxage)) %>%
  aperm(c(2, 3, 1))
OM_condition <- OM

SRA <- SRA_scope(OM_condition, condition = "catch2", Chist = SRA_data$Chist, Index = SRA_data$Index, I_sd = SRA_data$I_sd, I_type = SRA_data$I_type,
                 selectivity = rep("logistic", 2), s_selectivity = rep("logistic", 5), length_bin = 0.1 * SRA_data$length_bin, cores = 12, #mean_fit = TRUE,
                 s_CAA = SRA_data$s_CAA, vul_par = SRA_data$vul_par, map_s_vul_par = SRA_data$map_s_vul_par,
                 map_log_rec_dev = SRA_data$map_log_rec_dev)
saveRDS(SRA, file = "mse/OM/pinniped.rds")
#SRA <- readRDS("mse/OM/pinniped.rds")
plot(SRA, file = "mse/OM/OM_pinniped", dir = getwd(), open_file = FALSE, f_name = SRA_data$f_name, s_name = SRA_data$s_name,
     MSY_ref = c(0.4, 0.8))




