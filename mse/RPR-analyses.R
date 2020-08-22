cores <- parallel::detectCores() / 2
library(MSEtool)

############ Condition operating models with SRA_scope and data
SRA_data <- readRDS("mse/scoping/SRA_data.rds")
data_names <- c("Chist", "Index", "I_sd", "I_type", "length_bin", "s_CAA", "CAA", "CAL", "I_units")
data_ind <- match(data_names, names(SRA_data))

OM_condition <- readRDS("mse/scoping/OM_2sim.rds")


SRA <- list()

# initial fix 1x weight dogfish
SRA[[1]] <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("free", 2),
                 s_selectivity = rep("logistic", 5), cores = 1,
                 vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                 map_s_vul_par = SRA_data$map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                 LWT = list(CAL = 0, CAA = 0))

# 4x upweight dogfish
SRA[[2]] <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("free", 2),
                      s_selectivity = rep("logistic", 5), cores = 1,
                      vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                      map_s_vul_par = SRA_data$map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                      LWT = list(CAL = 0, CAA = 0, Index = c(1, 4, 1, 1, 1)))

# 0.1x downweight dogfish
SRA[[3]] <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("free", 2),
                      s_selectivity = rep("logistic", 5), cores = 1,
                      vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                      map_s_vul_par = SRA_data$map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                      LWT = list(CAL = 0, CAA = 0, Index = c(1, 0.1, 1, 1, 1)))

# 4x weight dogfish - ex 2019
SRA_data2 <- SRA_data
SRA_data2$Index[102, 2] <- NA
SRA[[4]] <- SRA_scope(OM_condition, data = SRA_data2[data_ind], condition = "catch2", selectivity = rep("free", 2),
                      s_selectivity = rep("logistic", 5), cores = 1,
                      vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                      map_s_vul_par = SRA_data$map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                      LWT = list(CAL = 0, CAA = 0, Index = c(1, 4, 1, 1, 1)))

# 0.1x downweight survey age comps
SRA[[5]] <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("free", 2),
                      s_selectivity = rep("logistic", 5), cores = 1,
                      vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                      map_s_vul_par = SRA_data$map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                      LWT = list(CAL = 0, CAA = 0, s_CAA = 0.1, Index = c(1, 1, 1, 1, 1)))

# No survey age comps - fix HBLL selectivity
base_s_vul_par <- c(SRA[[1]]@mean_fit$report$s_LFS[1], SRA[[1]]@mean_fit$report$s_L5[1])
s_vul_par <- matrix(c(base_s_vul_par, 0.5), 3, 5)
map_s_vul_par <- matrix(NA, 3, 5)
SRA[[6]] <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("free", 2),
                      s_selectivity = rep("logistic", 5), cores = 1,
                      vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                      s_vul_par = s_vul_par, map_s_vul_par = map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                      LWT = list(CAL = 0, CAA = 0, s_CAA = 0, Index = c(1, 1, 1, 1, 1)))

# No survey age comps - fix HBLL selectivity and don't estimate rec devs
base_s_vul_par <- c(SRA[[1]]@mean_fit$report$s_LFS[1], SRA[[1]]@mean_fit$report$s_L5[1])
s_vul_par <- matrix(c(base_s_vul_par, 0.5), 3, 5)
map_s_vul_par <- matrix(NA, 3, 5)
SRA[[7]] <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("free", 2),
                      s_selectivity = rep("logistic", 5), cores = 1,
                      vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                      s_vul_par = s_vul_par, map_s_vul_par = map_s_vul_par,
                      map_log_rec_dev = rep(NA, 102),
                      LWT = list(CAL = 0, CAA = 0, s_CAA = 0, Index = c(1, 1, 1, 1, 1)))

# Initial fit with com sel = rec sel
vul_par2 <- SRA_data$vul_par
vul_par2[, 1] <- vul_par2[, 2]
SRA[[8]] <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("free", 2),
                      s_selectivity = rep("logistic", 5), cores = 1,
                      vul_par = vul_par2, map_vul_par = matrix(NA, 80, 2),
                      map_s_vul_par = SRA_data$map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                      LWT = list(CAL = 0, CAA = 0))

saveRDS(SRA, file = "mse/om/RPR_SRA.rds")
SRA <- readRDS("mse/om/RPR_SRA.rds")

##### Generate reports and plots
SRA <- readRDS("mse/om/RPR_SRA.rds")

library(latex2exp)
SRA_names <- c("Base OM", TeX("Dogfish survey $\\lambda = 0.1$"), "Exclude 2019 dogfish",
               TeX("HBLL ages $\\lambda = 0.1$"), "Exclude HBLL ages")
SRA <- SRA[c(2, 3, 4, 5, 6)]
Hist <- lapply(SRA, function(x) runMSE(x@OM, Hist = TRUE, silent = TRUE))

B_BMSY <- Map(function(SRA, Hist) {
  SSB <- SRA@mean_fit$report$E
  SSBMSY <- Hist@Ref$SSBMSY[1]
  SSB/SSBMSY
}, SRA = SRA, Hist = Hist)

png("mse/figures/alt_SRA_fit.png", height = 4, width = 6, units = "in", res = 400)
par(mar = c(5, 4, 1, 1))
matplot(1918:2020, do.call(cbind, B_BMSY), type = "n", lty = 1, xlab = "Year", ylab = expression(B/B[MSY]))
abline(h = 0, col = "grey")
abline(h = c(0.4, 0.8), lty = 2)
matlines(1918:2020, do.call(cbind, B_BMSY), lty = 1, lwd = 2)
legend("bottomleft", SRA_names, col = 1:5, lty = 1, lwd = 2)
dev.off()

#SRA_names <- c("1xDog", "4xDog", "0.1xDog", "4xDog_Ex2019", "0.1xCAA", "0xCAA", "0xCAA_NoDev", "LowComSel")
#
#f_name <- c("Com", "Rec")
#s_name <- c("HBLL", "Dogfish", "CPUE1", "CPUE2", "CPUE3")
#compare_SRA(SRA[[1]], SRA[[2]], SRA[[3]], SRA[[4]], filename = "sensitivity_dogfish", dir = getwd(),
#            f_name = f_name, s_name = s_name, MSY_ref = c(0.4, 0.8),
#            scenario = list(names = SRA_names[1:4], col = gplots::rich.colors(5)[-4], lty = rep(1, 4)))
#
#
#compare_SRA(SRA[[1]], SRA[[5]], SRA[[6]], SRA[[7]], filename = "sensitivity_CAA", dir = getwd(),
#            f_name = f_name, s_name = s_name, MSY_ref = c(0.4, 0.8),
#            scenario = list(names = SRA_names[c(1, 5:7)], col = gplots::rich.colors(5)[-4], lty = rep(1, 4)))
#
#compare_SRA(SRA[[1]], SRA[[8]], filename = "sensitivity_sel", dir = getwd(),
#            f_name = f_name, s_name = s_name, MSY_ref = c(0.4, 0.8),
#            scenario = list(names = SRA_names[c(1, 8)]))


