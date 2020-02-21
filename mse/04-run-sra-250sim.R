
library(MSEtool)
setup(8)

############ Condition operating models with SRA_scope and data
SRA_data <- readRDS("mse/scoping/SRA_data.rds")
data_names <- c("Chist", "Index", "I_sd", "I_type", "length_bin", "s_CAA", "CAA", "CAL", "I_basis")
data_ind <- match(data_names, names(SRA_data))

OM_condition <- readRDS("mse/scoping/OM_250sim.rds")

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

saveRDS(SRA, file = "mse/om/upweight_dogfish.rds")
#SRA_list <- readRDS("mse/scoping/scoping_upweight_dogfish.rds")
#SRA2 <- SRA_list[[1]]; ret2 <- SRA_list[[2]]
plot(SRA2, retro = ret2, file = "mse/scoping/scoping_upweight_dogfish", dir = getwd(), open_file = FALSE, f_name = SRA_data$f_name, s_name = SRA_data$s_name,
     MSY_ref = c(0.4, 0.8), render_args = list(output_format = "word_document"))

# Upweight dogfish survey, fix HBLL sel from base
SRA <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("free", 2),
                 s_selectivity = rep("logistic", 5), cores = 8,
                 vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                 s_vul_par = s_vul_par, map_s_vul_par = map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                 LWT = list(CAL = 0, CAA = 0, Index = c(1, 4, 1, 1, 1)))

saveRDS(SRA, file = "mse/om/updog_fixsel.rds")
#SRA_list <- readRDS("mse/scoping/scoping_updog_fixsel.rds")
#SRA3 <- SRA_list[[1]]; ret3 <- SRA_list[[2]]

plot(SRA3, retro = ret3, file = "mse/scoping/scoping_updog_fixsel", dir = getwd(), open_file = FALSE, f_name = SRA_data$f_name, s_name = SRA_data$s_name,
     MSY_ref = c(0.4, 0.8), render_args = list(output_format = "word_document"))

# Low catch - use gfdatabase estimates of commerical catch in 1986-2005
SRA_data2 <- SRA_data
SRA_data2$Chist[match(1986:2005, SRA_data2$Year), 1] <- 0.5 * SRA_data2$Chist[match(1986:2005, SRA_data2$Year), 1]
SRA <- SRA_scope(OM_condition, data = SRA_data2[data_ind], condition = "catch2", selectivity = rep("free", 2),
                 s_selectivity = rep("logistic", 5), cores = 1,
                 vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                 map_s_vul_par = SRA_data$map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                 LWT = list(CAL = 0, CAA = 0, Index = c(1, 4, 1, 1, 1)))

saveRDS(SRA, file = "mse/om/lowcatch.rds")


# Low catch - fix HBLL sel from base
SRA <- SRA_scope(OM_condition, data = SRA_data2[data_ind], condition = "catch2", selectivity = rep("free", 2),
                 s_selectivity = rep("logistic", 5), cores = 8,
                 vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                 s_vul_par = s_vul_par, map_s_vul_par = map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                 LWT = list(CAL = 0, CAA = 0, Index = c(1, 4, 1, 1, 1)))

saveRDS(SRA, file = "mse/om/lowcatch_fixsel.rds")


## Try to estimate fishery selectivity
SRA_data$vul_par[1:2, ] <- c(50, 40, 30, 25)
map_vul_par <- matrix(NA, 80, 2)
map_vul_par[1:2, ] <- 1:4
SRA6 <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("logistic", 2),
                 s_selectivity = rep("logistic", 5), cores = 1,
                 vul_par = SRA_data$vul_par, map_vul_par = map_vul_par,
                 map_s_vul_par = SRA_data$map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                 LWT = list(CAL = 1, CAA = 20, Index = c(1, 4, 1, 1, 1)))

ret6 <- retrospective(SRA6, 11)
saveRDS(list(SRA6, ret6), file = "mse/scoping/scoping_estfisherysel_esthbllsel.rds")
SRA_list <- readRDS("mse/scoping/scoping_estfisherysel_esthbllsel.rds")
SRA6 <- SRA_list[[1]]; ret6 <- SRA_list[[2]]

plot(SRA6, retro = ret6, file = "mse/scoping/scoping_estfisherysel_esthbll_sel", dir = getwd(), open_file = FALSE, f_name = SRA_data$f_name, s_name = SRA_data$s_name,
     MSY_ref = c(0.4, 0.8), render_args = list(output_format = "word_document"))


## Compare plots
compare_SRA(SRA, SRA2, SRA3, SRA4, SRA5, SRA6,
            scenario = list(names = c("base", "upweight dogfish", "up.dog. fix HBLL sel",
                                      "low catch", "low catch fix HBLL sel", "est fishery/HBLL sel")),
            filename = "mse/scoping/compare_scoping", dir = getwd(), open_file = FALSE, f_name = SRA_data$f_name, s_name = SRA_data$s_name,
            MSY_ref = c(0.4, 0.8), render_args = list(output_format = "word_document"))



# Grid M and steepness
library(dplyr)
DLMtool::setup(8)

LH_grid <- expand.grid(M = seq(0.02, 0.07, 0.01), h = seq(0.65, 0.75, 0.01))

OM_condition <- readRDS("mse/scoping/OM_2sim.rds")
OM_condition@nsim <- nrow(LH_grid)
OM_condition@cpars$M <- LH_grid$M
OM_condition@cpars$h <- LH_grid$h

Mat_age <- OM_condition@cpars$Mat_age[1,,1]
OM_condition@cpars$Mat_age <- array(Mat_age,
                                    c(OM_condition@maxage, OM_condition@nyears + OM_condition@proyears, OM_condition@nsim)) %>%
  aperm(perm = c(3, 1, 2))


# Upweight dogfish
SRA7 <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("free", 2),
                 s_selectivity = rep("logistic", 5), cores = 8,
                 vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                 map_s_vul_par = SRA_data$map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                 LWT = list(CAL = 0, CAA = 0, Index = c(1, 4, 1, 1, 1)))
saveRDS(SRA7, file = "mse/scoping/profile_M_and_h.rds")
SRA7 <- readRDS("mse/scoping/profile_M_and_h.rds")
plot(SRA7, sims = LH_grid$h == 0.71, file = "mse/scoping/profile_M", dir = getwd(), open_file = FALSE, f_name = SRA_data$f_name, s_name = SRA_data$s_name,
     MSY_ref = c(0.4, 0.8), scenarios = list(names = paste0("M = 0.0", 2:7), col = 1:6))

# Upweight dog. fix HBLL sel
base_s_vul_par <- c(SRA@mean_fit$report$s_LFS[1], SRA@mean_fit$report$s_L5[1])
s_vul_par <- matrix(c(base_s_vul_par, 0.5), 3, 5)
map_s_vul_par <- matrix(NA, 3, 5)





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


