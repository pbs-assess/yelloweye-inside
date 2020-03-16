
# Miscellaneous plots generally associated with conditioning OM
library(MSEtool)

# Plot B/BMSY estimated from state-space surplus production model
SRA_list <- readRDS("mse/scoping/scoping_base.rds")
SRA <- SRA_list[[1]]; ret <- SRA_list[[2]]

SP_mod <- SP_SS(Data = SRA@OM@cpars$Data, AddInd = 1:5, use_r_prior = TRUE, start = list(r_prior = c(0.068, 0.03)))
#plot(SP_mod, dir = getwd(), filename = "mse/scoping/scoping_SP", open_file = FALSE)

# Little to no difference with reconstructed catch (no doubling 1986-2005)
#Data_lowcatch <- SRA@OM@cpars$Data
#Data_lowcatch@Cat[, match(1986:2005, Data_lowcatch@Year)] <- 0.5 * Data_lowcatch@Cat[, match(1986:2005, Data_lowcatch@Year)]
#SP_mod2 <- SP_SS(Data = Data_lowcatch, AddInd = 1:5, use_r_prior = TRUE, start = list(r_prior = c(0.068, 0.03)))

png("mse/figures/SP_fit.png", height = 3.5, width = 4.5, units = "in", res = 220)
par(mar = c(5, 4, 1, 1))
plot(as.numeric(names(SP_mod@B_BMSY)), SP_mod@B_BMSY, typ = "l", lwd = 3, xlab = "Year", ylab = expression(B/B[MSY]),
     ylim = c(0, 2.1))
abline(h = 0, col = "grey")
abline(h = c(0.4, 0.8), lty = 3)
dev.off()


# Plot retrospective bias with respect to M from 3 operating models:
#- Initial fit not used in reference/robustness set
#- Upweight dogfish (OM #1)
#- Upweight dogfish, estimate HBLL sel (OM #4)

############ Condition operating models with SRA_scope and data
SRA_data <- readRDS("mse/scoping/SRA_data.rds")
data_names <- c("Chist", "Index", "I_sd", "I_type", "length_bin", "s_CAA", "CAA", "CAL", "I_units")
data_ind <- match(data_names, names(SRA_data))

OM_condition <- readRDS("mse/scoping/OM_2sim.rds")


# Initial fit
M_vec <- seq(0.02, 0.06, 0.01)
setup(length(M_vec))
sfExportAll()
ret_init_fit <- sfLapply(M_vec, function(x) {
  OM_condition@M <- rep(x, 2)
  SRA <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("free", 2),
                   s_selectivity = rep("logistic", 5), cores = 1,
                   vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                   map_s_vul_par = SRA_data$map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                   LWT = list(CAL = 0, CAA = 0))
  retrospective(SRA, nyr = 11, figure = FALSE)
})


ret_updog_fixsel <- sfLapply(M_vec, function(x) {
  OM_condition@M <- rep(x, 2)

  SRA <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("free", 2),
                   s_selectivity = rep("logistic", 5), cores = 1,
                   vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                   map_s_vul_par = SRA_data$map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                   LWT = list(CAL = 0, CAA = 0))

  base_s_vul_par <- c(SRA@mean_fit$report$s_LFS[1], SRA@mean_fit$report$s_L5[1])
  s_vul_par <- matrix(c(base_s_vul_par, 0.5), 3, 5)
  map_s_vul_par <- matrix(NA, 3, 5)

  SRA3 <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("free", 2),
                    s_selectivity = rep("logistic", 5), cores = 1,
                    vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                    s_vul_par = s_vul_par, map_s_vul_par = map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                    LWT = list(CAL = 0, CAA = 0, Index = c(1, 4, 1, 1, 1)))
  retrospective(SRA3, nyr = 11, figure = FALSE)
})


ret_updog_estsel <- sfLapply(M_vec, function(x) {
  OM_condition@M <- rep(x, 2)
  SRA <- SRA_scope(OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("free", 2),
                   s_selectivity = rep("logistic", 5), cores = 1,
                   vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 80, 2),
                   map_s_vul_par = SRA_data$map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
                   LWT = list(CAL = 0, CAA = 0, Index = c(1, 4, 1, 1, 1)))
  retrospective(SRA, nyr = 11, figure = FALSE)
})

rho <- cbind %>% do.call(lapply(list(ret_init_fit, ret_updog_fixsel, ret_updog_estsel),
                                function(x) vapply(x, function(xx) summary(xx)[3, 1], numeric(1))))


png("mse/figures/conditioning/retrospective_M.png", height = 4, width = 4, units = "in", res = 220)
par(mar = c(5, 4, 1, 1))
matplot(M_vec, rho, typ = "o", pch = 16, lty = 1, xlab = "Natural mortality", ylab = "SSB Mohn's rho")
abline(h = 0, lty = 3)
legend("bottomleft", c("Initial fit", "(1) Upweight dogfish", "(4) Estimate HBLL selectivity"), col = 1:3, pch = 16, cex=0.6, bty="n")
dev.off()
