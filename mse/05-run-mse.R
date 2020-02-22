
library(MSEtool)
source("mse/YE_MPs.R")
setup(12)
sfExportAll()
sfLibrary(gfdlm)


scenario <- c("updog_fixsel", "lowcatch_fixsel", "episodic_recruitment", "lowM_fixsel", "high_index_cv")
subsets <- c("fixedTAC", "index", "IDXSP", "SP")

for(i in 1:length(scenario)) {
  SRA <- readRDS(paste0("mse/om/", scenario[i], ".rds"))
  myOM <- SRA@OM

  message("Running ", scenario[i], " scenario...")

  # Constant catch scenarios
  myOM@interval <- c(1, rep(200, 4))
  message("Constant catch and ref MPs...")
  myMSE1 <- runMSE(myOM, MPs = c("FMSYref", "NF", "CC_5t", "CC_10t", "CC_15t"), parallel = TRUE)
  saveRDS(myMSE1, file = paste0("mse/om/MSE_", scenario[i], "_fixedTAC.rds"))

  ## All other index-based MPs
  #myOM@interval <- 1
  #message("Other index-based MPs...")
  #myMSE2 <- runMSE(myOM, MPs = c("ICI_YE", "ICI2_YE", "Iratio_YE", "GB_slope_YE", "IT5_YE", "IT10_YE", "Islope_YE"), parallel = TRUE)
  #saveRDS(myMSE2, file = paste0("mse/om/MSE_", scenario[i], "_index.rds"))

  # IDX and SP_4080 with 5 and 10 yr intervals
  myOM@interval <- c(1, 1, 5, 10, 5, 10)
  message("IDX and SP...")
  myMSE3 <- runMSE(myOM, MPs = c("IDX_YE", "IDX_smooth_YE", "SP_4080_5f", "SP_4080_10f", "SP_2060_5f", "SP_2060_10f"), parallel = TRUE)
  saveRDS(myMSE3, file = paste0("mse/om/MSE_", scenario[i], "_IDX_SP.rds"))

  # SP_2060 with 5 and 10 yr intervals and interim SP
  #myOM@interval <- c(5, 10, 1)
  #message("SP_2060 and interim SP_MSY...")
  #myMSE4 <- runMSE(myOM, MPs = c("SP_2060_5f", "SP_2060_10f", "SP_interim"), parallel = TRUE)
  #saveRDS(myMSE4, file = paste0("mse/OM/MSE_", scenario[i], "_SP.rds"))

  #out <- merge_MSE(myMSE1, myMSE2, myMSE3, myMSE4)
  #saveRDS(out, file = paste0("mse/om/MSE_", scenario[i], ".rds"))

  message("Done for ", scenario[i], ".\n")
}

sfStop()
