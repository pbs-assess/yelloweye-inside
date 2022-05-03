
# Run on DLMtool v.5.4.3 in R 3.6.2
# Can't run in parallel!
library(MSEtool)
#library(gfdlm)
#source("mse/YE_MPs.R")
#sfExportAll()
#sfLibrary(gfdlm)


make_effort_OM <- function(FMSYvec, F_terminal, relF = 1) {
  force(FMSYvec)
  force(F_terminal)

  MP <- function(x = 1, Data, reps = 1) {
    Ftarget <- relF * FMSYvec[x]
    Rec <- new("Rec")
    Rec@Effort <- rep(Ftarget/F_terminal[x], reps)
    return(Rec)
  }
  formals(MP)$relF <- relF
  structure(MP, class = "MP")
}


scenario <- c("updog_fixsel", "lowcatch_fixsel", "episodic_recruitment", "lowM_fixsel", "high_index_cv", "upweight_dogfish")

#SRA <- readRDS("mse/scoping/scoping_base.rds")[[1]]
#nsim <- 250
#set.seed(24)
#AddIerr <- rnorm(nsim * 5 * (SRA@OM@nyears + SRA@OM@proyears), -0.5 * 0.25^2, 0.25) %>% exp() %>%
#  array(dim = c(nsim, 5, SRA@OM@nyears + SRA@OM@proyears))
#rm(SRA)

for(i in 1:length(scenario)) {
  SRA <- readRDS(paste0("mse/om/", scenario[i], ".rds"))
  myOM <- SRA@OM
  myOM@interval <- 200

  myOM@Prob_staying <- myOM@Frac_area_1 <- myOM@Size_area_1 <- c(0.5, 0.5)
  myOM@TAEFrac <- c(1, 1)
  myOM@TAESD <- myOM@qinc <- myOM@qcv <- c(0, 0)

  message("Running ", scenario[i], " scenario...")
  message("Getting hist")
  Hist <- runMSE(myOM, silent = TRUE, Hist = TRUE)

  message("Calculating FMSY")
  FMSYvec <- sapply(1:myOM@nsim, function(x) {
    optMSY_eq(x = x, M_ageArray = Hist@SampPars$M_ageArray,
              Wt_age = Hist@SampPars$Wt_age, Mat_age = Hist@SampPars$Mat_age,
              V = Hist@SampPars$V, maxage = Hist@SampPars$maxage[1],
              R0 = Hist@SampPars$R0, SRrel = Hist@SampPars$SRrel, hs = Hist@SampPars$hs,
              yr.ind = myOM@nyears, plusgroup = 1)["F"]
  })

  F_terminal <- Hist@AtAge$FM[, , myOM@nyears, 1] %>% apply(1, max)

  MP_FMSY <- make_effort_OM(FMSYvec, F_terminal)
  MP_50FMSY <- make_effort_OM(FMSYvec, F_terminal, relF = 0.5)
  MP_75FMSY <- make_effort_OM(FMSYvec, F_terminal, relF = 0.75)
  MP_125FMSY <- make_effort_OM(FMSYvec, F_terminal, relF = 1.25)

  message("Projecting...")
  myMSE <- runMSE(myOM, MPs = c("MP_FMSY", "MP_50FMSY", "MP_75FMSY", "MP_125FMSY"), silent = TRUE)
  saveRDS(myMSE, file = paste0("mse/om/MSE_", scenario[i], "_FMSYproj.rds"))

  message("Done.\n")
}

