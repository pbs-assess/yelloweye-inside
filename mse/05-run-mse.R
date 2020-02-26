
library(MSEtool)
library(gfdlm)
source("mse/YE_MPs.R")
setup(12)
sfExportAll()
sfLibrary(gfdlm)


scenario <- c("updog_fixsel", "lowcatch_fixsel", "episodic_recruitment", "lowM_fixsel", "high_index_cv", "upweight_dogfish")
#scenario <- "upweight_dogfish"
subsets <- c("fixedTAC", "index", "IDXSP", "SP")

for(i in 1:length(scenario)) {
  SRA <- readRDS(paste0("mse/om/", scenario[i], ".rds"))
  myOM <- SRA@OM

  message("Running ", scenario[i], " scenario...")

  # Constant catch scenarios
  myOM@interval <- c(1, 1, rep(200, 4))
  message("Constant catch and ref MPs...")
  myMSE1 <- runMSE(myOM, MPs = c("FMSYref", "FMSYref75", "NFref", "CC_5t", "CC_10t", "CC_15t"), parallel = TRUE)
  saveRDS(myMSE1, file = paste0("mse/om/MSE_", scenario[i], "_fixedTAC.rds"))

  # Index slope MPs
  myOM@interval <- 1
  message("Index slope MPs...")
  myMSE2 <- runMSE(myOM, MPs = c("GB_slope_lambda1", "GB_slope_lambda05", "GB_slope_yrs10",
                                 "Islope_5_lambda04", "Islope_10_lambda04", "Islope_10_lambda08"), parallel = TRUE)
  saveRDS(myMSE2, file = paste0("mse/om/MSE_", scenario[i], "_index_slope.rds"))

  # Index ratio MPs
  myOM@interval <- 1
  message("Index ratio MPs...")
  myMSE3 <- runMSE(myOM, MPs = c("Iratio_23", "Iratio_510", "IT5_mc05", "IT5_mc025", "IT10_mc05", "IT10_mc025",
                                 "Itarget_5", "Itarget_10"), parallel = TRUE)
  saveRDS(myMSE3, file = paste0("mse/om/MSE_", scenario[i], "_index_ratio.rds"))

  # IDX and SP_4080 with 5 and 10 yr intervals
  myOM@interval <- c(1, 1, 5, 10, 5, 10)
  message("IDX and SP...")
  myMSE4 <- runMSE(myOM, MPs = c("IDX_YE", "IDX_smooth_YE", "SP_4080_5f", "SP_4080_10f", "SP_2060_5f", "SP_2060_10f"), parallel = TRUE)
  saveRDS(myMSE4, file = paste0("mse/om/MSE_", scenario[i], "_IDX_SP.rds"))

  # SP_2060 with 5 and 10 yr intervals and interim SP
  #myOM@interval <- c(5, 10, 1)
  #message("SP_2060 and interim SP_MSY...")
  #myMSE4 <- runMSE(myOM, MPs = c("SP_2060_5f", "SP_2060_10f", "SP_interim"), parallel = TRUE)
  #saveRDS(myMSE4, file = paste0("mse/OM/MSE_", scenario[i], "_SP.rds"))

  out <- merge_MSE(myMSE1, myMSE2, myMSE3, myMSE4)
  saveRDS(out, file = paste0("mse/om/MSE_", scenario[i], ".rds"))

  message("Done for ", scenario[i], ".\n")
}

sfStop()


for(i in 1:length(scenario)) {

  myMSE1 <- readRDS(paste0("mse/om/MSE_", scenario[i], "_fixedTAC.rds"))
  myMSE2 <- readRDS(paste0("mse/om/MSE_", scenario[i], "_index_slope.rds"))
  myMSE3 <- readRDS(paste0("mse/om/MSE_", scenario[i], "_index_ratio.rds"))
  myMSE4 <- readRDS(paste0("mse/om/MSE_", scenario[i], "_IDX_SP.rds"))

  out <- merge_MSE(myMSE1, myMSE2, myMSE3, myMSE4)
  saveRDS(out, file = paste0("mse/om/MSE_", scenario[i], ".rds"))

  message("Done for ", scenario[i], ".\n")
}


