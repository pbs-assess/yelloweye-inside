library(DLMtool) # version 5.4.1 to be released on CRAN for CSAS
library(MSEtool) # version 1.5 to be released on CRAN for CSAS
library(dplyr)
library(reshape2)

SRA_data <- readRDS("mse/scoping/SRA_data.rds")

################## OM setup
#### Stock
stock_yelloweye <- new("Stock")
stock_yelloweye@Name <- "Inside Yelloweye Rockfish"
stock_yelloweye@Common_Name <- "Yelloweye Rockfish"
stock_yelloweye@Species <- "Sebastes ruberrimus"
stock_yelloweye@maxage <- SRA_data$maxage # 80 as the plus group

stock_yelloweye@R0 <- 1200 # Start value of R0 for SRA_scope, will be estimated by conditioning

## Will be updated in cpars section below
stock_yelloweye@M <- c(0.045, 0.045)
stock_yelloweye@Msd <- stock_yelloweye@Linfsd <- stock_yelloweye@Ksd <-
  stock_yelloweye@L50 <- stock_yelloweye@L50_95 <- c(0, 0)
stock_yelloweye@h <- c(0.71, 0.71)

stock_yelloweye@SRrel <- 1
stock_yelloweye@Perr <- c(0.4, 0.4)

# Will be updated by SRA_scope, need a placeholder for now
#stock_yelloweye@AC <- stock_yelloweye@D <- c(0, 0)

## gfsynopsis
stock_yelloweye@Linf <- c(65.1, 65.1)
stock_yelloweye@K <- c(0.04, 0.04)
stock_yelloweye@t0 <- c(-9.04, -9.04)
stock_yelloweye@LenCV <- c(0.14, 0.14) #only used if we use length data
stock_yelloweye@a <- 1.31e-5
stock_yelloweye@b <- 3.09

stock_yelloweye@Fdisc <- c(1, 1) # 100% discard mortality

# Effectively parameterize a single area model
stock_yelloweye@Frac_area_1 <- stock_yelloweye@Size_area_1 <- stock_yelloweye@Prob_staying <- c(0.5, 0.5)


#### Fleet
fleet_yelloweye <- new("Fleet")
fleet_yelloweye@nyears <- length(SRA_data$Year) #101
fleet_yelloweye@Spat_targ <- c(1, 1) # No spatial dynamics
fleet_yelloweye@DR <- c(0, 0) # All removals will be accounted for in the catch
fleet_yelloweye@CurrentYr <- 2019

# Placeholder, not used. Effort/F/selectivity will be updated by SRA_scope, qinc/qcv only affect effort MPs
#fleet_yelloweye@EffYears <- c(1, fleet_yelloweye@nyears)
fleet_yelloweye@qinc <- fleet_yelloweye@qcv <- c(0, 0)
#fleet_yelloweye@Esd <- fleet_yelloweye@qinc <- fleet_yelloweye@qcv <- c(0, 0)

# Selectivity parameters will be starting values for HBLL survey selectivity (fixed for all other fleets/surveys)
fleet_yelloweye@L5 <- c(33, 33)
fleet_yelloweye@LFS <- c(45, 45)
fleet_yelloweye@Vmaxlen <- c(1, 1)
fleet_yelloweye@isRel <- FALSE

#### Observation and Implementation
# We will only care about Cobs and Cbiascv. Index error is handled internally based on SRA_scope inputs and model fit.
# We don't use age/length data, or any life history information for MPs.
obs_yelloweye <- DLMtool::Precise_Unbiased
imp_yelloweye <- DLMtool::Perfect_Imp


OM <- new("OM", stock_yelloweye, fleet_yelloweye, obs_yelloweye, imp_yelloweye)
OM@Name <- "Inside Yelloweye Rockfish"
OM@proyears <- 100

################# OM with 2 sims
OM2 <- OM
OM2@nsim <- 2

Mat_age <- ifelse(1:OM2@maxage <= 7, 0, 1/(1 + exp(-log(19) * (1:OM2@maxage - 14.4)/(27.4-14.4))))
OM2@cpars$Mat_age <- array(Mat_age, c(OM2@maxage, OM2@nyears + OM2@proyears, OM2@nsim)) %>% aperm(perm = c(3, 1, 2))

saveRDS(OM2, file = "mse/scoping/OM2.rds")


################## OM with 250 sims
OM@nsim <- 250

################# cpars (custom parameters)
# Maturity: from gfsynopsis 50% maturity at 14.4 years, 95% maturity at 27.4, immature for all ages <= 7
Mat_age <- ifelse(1:OM@maxage <= 7, 0, 1/(1 + exp(-log(19) * (1:OM@maxage - 14.4)/(27.4-14.4))))
OM@cpars$Mat_age <- array(Mat_age, c(OM@maxage, OM@nyears + OM@proyears, OM@nsim)) %>% aperm(perm = c(3, 1, 2))

# Natural mortality M ~ lognormal mean 0.045, sd = 0.2
set.seed(91283)
M_samps <- rlnorm(OM@nsim, log(0.045) - 0.5 * 0.2^2, 0.2)
OM@cpars$M <- M_samps

# Sample steepness h ~ transformed beta with mean = 0.71, sd = 0.15.
# x = (h - 0.2)/0.8 where x ~ Beta(mean = 0.6375, sd = 0.12)
h_alpha <- alphaconv(0.6375, 0.12)
h_beta <- betaconv(0.6375, 0.12)
set.seed(65423)

h_samps <- rbeta(OM@nsim, h_alpha, h_beta)
h_samps <- 0.8 * h_samps + 0.2
OM@cpars$h <- h_samps

saveRDS(OM, file = "mse/scoping/OM.rds")

### Plot M and h samples
#hist(M_samps, xlab = "Natural mortality", main = ""); abline(v = 0.045, lwd = 2, lty = 3)
#hist(h_samps, xlab = "Steepness", main = ""); abline(v = 0.71, lwd = 2, lty = 3)


############ Condition operating models with SRA_scope and data
#setup(12)

# Base
OM_condition <- OM
SRA <- SRA_scope(OM_condition, condition = "catch2", Chist = SRA_data$Chist, Index = SRA_data$Index, I_sd = SRA_data$I_sd, I_type = SRA_data$I_type,
                 selectivity = rep("logistic", 2), s_selectivity = rep("logistic", 5), length_bin = 0.1 * SRA_data$length_bin, cores = 1,
                 s_CAA = SRA_data$s_CAA, vul_par = SRA_data$vul_par, map_s_vul_par = SRA_data$map_s_vul_par,
                 map_log_rec_dev = SRA_data$map_log_rec_dev)
saveRDS(SRA, file = "mse/OM/base.rds")
#SRA <- readRDS("mse/OM/base.rds")

plot(SRA, file = "mse/OM/OM_base", dir = getwd(), open_file = FALSE, f_name = SRA_data$f_name, s_name = SRA_data$s_name,
     MSY_ref = c(0.4, 0.8))

# Upweight dogfish
SRA <- SRA_scope(OM_condition, condition = "catch2", Chist = SRA_data$Chist, Index = SRA_data$Index, I_sd = SRA_data$I_sd, I_type = SRA_data$I_type,
                 selectivity = rep("logistic", 2), s_selectivity = rep("logistic", 5), length_bin = 0.1 * SRA_data$length_bin, cores = 12,
                 s_CAA = SRA_data$s_CAA, LWT = list(Index = c(1, 4, 1, 1, 1)),
                 vul_par = SRA_data$vul_par, map_s_vul_par = SRA_data$map_s_vul_par,
                 map_log_rec_dev = SRA_data$map_log_rec_dev)

saveRDS(SRA, file = "mse/OM/upweight_dogfish.rds")
#SRA <- readRDS("mse/OM/upweight_dogfish.rds")
plot(SRA, file = "mse/OM/OM_upweight_dogfish", dir = getwd(), open_file = FALSE, f_name = SRA_data$f_name, s_name = SRA_data$s_name,
     MSY_ref = c(0.4, 0.8))

# Low catch - use gfdatabase estimates of commerical catch in 1986-2005
SRA_data$Chist[match(1986:2005, SRA_data$Year), 1] <- 0.5 * SRA_data$Chist[match(1986:2005, SRA_data$Year), 1]
sfExportAll()
SRA <- SRA_scope(OM_condition, condition = "catch2", Chist = SRA_data$Chist, Index = SRA_data$Index, I_sd = SRA_data$I_sd, I_type = SRA_data$I_type,
                 selectivity = rep("logistic", 2), s_selectivity = rep("logistic", 5), length_bin = 0.1 * SRA_data$length_bin, cores = 12,
                 s_CAA = SRA_data$s_CAA, LWT = list(Index = c(1, 4, 1, 1, 1)),
                 vul_par = SRA_data$vul_par, map_s_vul_par = SRA_data$map_s_vul_par,
                 map_log_rec_dev = SRA_data$map_log_rec_dev)

saveRDS(SRA, file = "mse/OM/low_catch.rds")
#SRA <- readRDS("mse/OM/OM_low_catch.rds")
plot(SRA, file = "mse/OM/OM_low_catch", dir = getwd(), open_file = FALSE, f_name = SRA_data$f_name, s_name = SRA_data$s_name,
     MSY_ref = c(0.4, 0.8))


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




