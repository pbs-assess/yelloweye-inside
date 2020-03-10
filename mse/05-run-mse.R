
library(MSEtool)
library(gfdlm)
source("mse/YE_MPs.R")
setup(parallel::detectCores()/2)
sfExportAll()
sfLibrary(gfdlm)


scenario <- c("updog_fixsel", "lowcatch_fixsel", "episodic_recruitment", "lowM_fixsel", "high_index_cv", "upweight_dogfish")

SRA <- readRDS("mse/scoping/scoping_base.rds")[[1]]
nsim <- 250
set.seed(24)
AddIerr <- rnorm(nsim * 5 * (SRA@OM@nyears + SRA@OM@proyears), -0.5 * 0.25^2, 0.25) %>% exp() %>%
  array(dim = c(nsim, 5, SRA@OM@nyears + SRA@OM@proyears))
rm(SRA)

for(i in 1:length(scenario)) {
  SRA <- readRDS(paste0("mse/om/", scenario[i], ".rds"))
  myOM <- SRA@OM

  message("Running ", scenario[i], " scenario...")

  # Constant catch scenarios and reference MPs
  myOM@interval <- c(1, 1, rep(200, 4))
  myOM@cpars$AddIerr <- AddIerr
  message("Constant catch and ref MPs...")
  myMSE1 <- runMSE(myOM, MPs = c("FMSYref", "FMSYref75", "NFref", "CC_5t", "CC_10t", "CC_15t"), parallel = TRUE)
  saveRDS(myMSE1, file = paste0("mse/om/MSE_", scenario[i], "_fixedTAC.rds"))
  rm(myMSE1)

  # Index MPs annual update
  myOM@interval <- 1
  message("Index MPs...")
  myMSE2 <- runMSE(myOM, MPs = c("GB_slope_lambda1", "GB_slope_lambda05", "GB_slope_yrs10",
                                 "Islope_5_lambda04", "Islope_10_lambda04", "Islope_10_lambda08",
                                 "Iratio_23", "Iratio_55"), parallel = TRUE)
  saveRDS(myMSE2, file = paste0("mse/om/MSE_", scenario[i], "_index.rds"))
  rm(myMSE2)

  # Index MPs 5 year interval
  myOM@interval <- 5
  message("Index MPs (5 year interval)...")
  myMSE3 <- runMSE(myOM, MPs = c("GB_slope_lambda1_5u", "GB_slope_lambda05_5u", "GB_slope_yrs10_5u",
                                 "Islope_5_lambda04_5u", "Islope_10_lambda04_5u", "Islope_10_lambda08_5u",
                                 "Iratio_23_5u", "Iratio_55_5u"), parallel = TRUE)
  saveRDS(myMSE3, file = paste0("mse/om/MSE_", scenario[i], "_index_5y.rds"))
  rm(myMSE3)

  # IDX
  myOM@interval <- c(1, 1, 1, 1, 5, 5, 5, 5)
  message("IDX...")
  myMSE4 <- runMSE(myOM, MPs = c("IDX_", "IDX_yrs5", "IDX_smooth_", "IDX_smooth_yrs5",
                                 "IDX_5u", "IDX_yrs5_5u", "IDX_smooth_5u", "IDX_smooth_yrs5_5u"), parallel = TRUE)
  saveRDS(myMSE4, file = paste0("mse/om/MSE_", scenario[i], "_IDX.rds"))
  rm(myMSE4)

  # SP with 5 and 10 yr intervals
  myOM@interval <- c(5, 5, 10, 10)
  message("SP...")
  myMSE5 <- runMSE(myOM, MPs = c("SP_8040_5u", "SP_4010_5u", "SP_8040_10u", "SP_4010_10u"), parallel = TRUE)
  saveRDS(myMSE5, file = paste0("mse/OM/MSE_", scenario[i], "_SP.rds"))
  rm(myMSE5)

  #out <- merge_MSE(myMSE1, myMSE2, myMSE3, myMSE4)
  #saveRDS(out, file = paste0("mse/om/MSE_", scenario[i], ".rds"))

  message("Done for ", scenario[i], ".\n")
}

sfStop()

#for(i in 1:length(scenario)) {
#
#  myMSE1 <- readRDS(paste0("mse/om/MSE_", scenario[i], "_fixedTAC.rds"))
#  myMSE2 <- readRDS(paste0("mse/om/MSE_", scenario[i], "_index_slope.rds"))
#  myMSE3 <- readRDS(paste0("mse/om/MSE_", scenario[i], "_index_ratio.rds"))
#  myMSE4 <- readRDS(paste0("mse/om/MSE_", scenario[i], "_IDX_SP.rds"))
#
#  out <- merge_MSE(myMSE1, myMSE2, myMSE3, myMSE4)
#  saveRDS(out, file = paste0("mse/om/MSE_", scenario[i], ".rds"))
#
#  message("Done for ", scenario[i], ".\n")
#}
#
