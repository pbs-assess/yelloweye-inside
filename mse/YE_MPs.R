

##### Constant catch MPs of 5, 10, or 15 tonnes
CC_5t <- function(x, Data, reps) {
  Rec <- new("Rec")
  Rec@TAC <- rep(5, reps)
  return(Rec)
}
class(CC_5t) <- "MP"

CC_10t <- function(x, Data, reps) {
  Rec <- new("Rec")
  Rec@TAC <- rep(10, reps)
  return(Rec)
}
class(CC_10t) <- "MP"

CC_15t <- function(x, Data, reps) {
  Rec <- new("Rec")
  Rec@TAC <- rep(15, reps)
  return(Rec)
}
class(CC_15t) <- "MP"


##### MP wrapper function so that it uses HBLL
YE_wrapper <- function(MP, add_min_TAC = 0.5, ...) {
  MP <- substitute(MP)
  MP_call <- as.call(c(MP, x = quote(x), Data = quote(Data), reps = quote(reps), list(...)))

  MP_body <- bquote({
    Data@Ind <- Data@AddInd[, 1, ]
    Data@CV_Ind <- Data@CV_AddInd[, 1, ]
    Data@Iref <- 2 * apply(Data@AddInd[, 1, 93:102], 1, mean, na.rm = TRUE)
    Rec <- .(MP_call)
    if(!all(is.na(Rec@TAC))) Rec@TAC[Rec@TAC < .(add_min_TAC)] <- .(add_min_TAC)
    return(Rec)
  })

  YE_MP <- eval(call("function", as.pairlist(alist(x = 1, Data = , reps = 1)), MP_body))
  class(YE_MP) <- "MP"
  return(YE_MP)
}

##### Index-based HCRs that uses HBLL
# IDX
IDX_ <- IDX_5u <- YE_wrapper(gfdlm::IDX, tac_floor = 5)
IDX_yrs5 <- IDX_yrs5_5u <- YE_wrapper(gfdlm::IDX, tac_floor = 5, year_ref = 5)
IDX_smooth_ <- IDX_smooth_5u <- YE_wrapper(gfdlm::IDX_smooth, tac_floor = 5)
IDX_smooth_yrs5 <- IDX_smooth_yrs5_5u <- YE_wrapper(gfdlm::IDX_smooth, tac_floor = 5, year_ref = 5)

# Iratio
Iratio_23 <- Iratio_23_5u <- YE_wrapper(Iratio)
Iratio_55 <- Iratio_55_5u <- YE_wrapper(Iratio, yrs = c(5, 10))

# I slope
GB_slope_lambda1 <- GB_slope_lambda1_5u <- YE_wrapper(GB_slope)
GB_slope_lambda05 <- GB_slope_lambda05_5u <- YE_wrapper(GB_slope, lambda = 0.5) # Slower changes in catch relative to index
GB_slope_yrs10 <- GB_slope_yrs10_5u <- YE_wrapper(GB_slope, yrsmth = 10) # Longer time window for smoothing

Islope_5_lambda04 <- Islope_5_lambda04_5u <- YE_wrapper(Islope1, xx = 0)
Islope_10_lambda04 <- Islope_10_lambda04_5u <- YE_wrapper(Islope1, xx = 0, yrsmth = 10)
Islope_10_lambda08 <- Islope_10_lambda08_5u <- YE_wrapper(Islope1, xx = 0, yrsmth = 10, lambda = 0.8)

# Ratio of mean index to Iref = 2 * mean 2010-2019 index
#IT_YE <- function(x, Data, reps = 100, plot=FALSE, yrsmth = 5, mc = 0.05) {
#  dependencies = "Data@Ind, Data@MPrec, Data@CV_Ind, Data@Iref"
#  ind <- max(1, (length(Data@Year) - yrsmth + 1)):length(Data@Year)
#  if(is.na(Data@Iref[x])) return(list(TAC=rep(as.numeric(NA), reps)))
#  deltaI <- mean(Data@Ind[x, ind], na.rm=TRUE)/Data@Iref[x]
#  if (deltaI < (1 - mc)) deltaI <- 1 - mc
#  if (deltaI > (1 + mc)) deltaI <- 1 + mc
#  TAC <- Data@MPrec[x] * deltaI * trlnorm(reps, 1, if(is.na(Data@CV_Ind[x, 1])) 0.1 else Data@CV_Ind[x, 1])
#  TAC <- TACfilter(TAC)
#  if (plot) {
#    op <- par(no.readonly = TRUE)
#    on.exit(op)
#    par(mfrow=c(1,2), oma=c(1,1,1,1), mar=c(5,4,1,4))
#    ylim <- range(c(Data@Ind[x,ind], Data@Iref[x]))
#    plot(Data@Year[ind], Data@Ind[x,ind], xlab="Year",
#         ylab= paste0("Index (previous ", yrsmth, "years)"), bty="l", type="l",
#         lwd=2, ylim=ylim)
#    lines(Data@Year[ind], rep(mean(Data@Ind[x, ind], na.rm=TRUE), length(ind)), lty=2)
#    text(quantile(Data@Year[ind],0.15), mean(Data@Ind[x, ind], na.rm=TRUE), "Mean Index", pos=3)
#    lines(Data@Year[ind], rep(mean(Data@Iref[x], na.rm=TRUE), length(ind)), lty=3)
#    text(quantile(Data@Year[ind],0.15), Data@Iref[x], "Reference Index", pos=3)
#
#    boxplot(TAC, ylab=paste0("TAC (", Data@Units, ")"))
#    points(1, Data@MPrec[x], cex=2, col="blue", pch=16)
#    text(1, Data@MPrec[x], cex=1, col="blue", "Last Catch", pos=1)
#  }
#  Rec <- new("Rec")
#  Rec@TAC <- TAC
#  return(Rec)
#}
#class(IT_YE) <- "MP"
#
#IT5_mc05 <- YE_wrapper(IT_YE) #
#IT5_mc025 <- YE_wrapper(IT_YE, mc = 0.025)
#
#IT10_mc05 <- YE_wrapper(IT_YE, mc = 0.05)
#IT10_mc025 <- YE_wrapper(IT_YE, mc = 0.025)

# Iref = 1.5 * mean 2010-2015 index
#Itarget_5 <- YE_wrapper(Itarget1)
#Itarget_10 <- YE_wrapper(Itarget1, yrsmth = 10)

##### SP model that only uses future HBLL and with r-prior
SP_YE <- function(x, Data, ...) {
  nyears <- ncol(Data@Cat)
  if(nyears > 102) Data@AddInd[, 2:5, 103:nyears] <- Data@CV_AddInd[, 2:5, 103:nyears] <- NA
  SP_SS(x = x, Data = Data, AddInd = 1:5, start = list(r_prior = c(0.068, 0.03)), ...)
}
class(SP_YE) <- "Assess"

##### SP with 80-40 BMSY HCR
SP_8040_ <- function(x, Data, reps = 1) {
  do_Assessment <- SP_YE(x = x, Data = Data)
  Rec <- HCR_ramp(Assessment = do_Assessment, reps = reps,
                  LRP = 0.4, TRP = 0.8, RP_type = "SSB_SSBMSY")
  Rec@Misc <- MSEtool:::Assess_diagnostic(x, Data, do_Assessment, include_assessment = FALSE)
  return(Rec)
}
class(SP_8040_) <- "MP"
SP_8040 <- YE_wrapper(SP_8040_)
SP_8040_5u <- SP_8040_10u <- SP_8040

##### SP with 40-10 B0 HCR
SP_4010_ <- function(x, Data, reps = 1) {
  do_Assessment <- SP_YE(x = x, Data = Data)
  Rec <- HCR_ramp(Assessment = do_Assessment, reps = reps,
                  LRP = 0.1, TRP = 0.4, RP_type = "SSB_SSBMSY")
  Rec@Misc <- MSEtool:::Assess_diagnostic(x, Data, do_Assessment, include_assessment = FALSE)
  return(Rec)
}
class(SP_4010_) <- "MP"
SP_4010 <- YE_wrapper(SP_4010_)
SP_4010_5u <- SP_4010_10u <- SP_4010



merge_MSE <- function(...) {
  dots <- list(...)
  if(length(dots) == 1) dots <- dots[[1]]

  slots_identical <- function(slotname, x = dots, is_logical = FALSE) {
    res <- lapply(x, getElement, slotname)
    is_identical <- all(vapply(res[-1], identical, logical(1), res[[1]]))
    if(is_logical) {
      return(is_identical)
    } else return(unique(do.call(c, res)))
  }

  slots_identical("Name")
  slots_identical("nyears")
  slots_identical("proyears")
  slots_identical("nsim")

  stopifnot(slots_identical("OM", is_logical = TRUE))
  stopifnot(slots_identical("Obs", is_logical = TRUE))
  stopifnot(slots_identical("SSB_hist", is_logical = TRUE))
  stopifnot(slots_identical("CB_hist", is_logical = TRUE))
  stopifnot(slots_identical("FM_hist", is_logical = TRUE))

  #nMPs <- vapply(dots, getElement, numeric(1), "nMPs")

  slotvec <- c("B_BMSY", "F_FMSY", "B", "SSB", "VB", "FM", "C", "TAC", "Effort", "PAA", "CAA", "CAL")
  res <- lapply(slotvec, function(x) do.call(abind::abind, c(lapply(dots, getElement, x), along = 2)))
  names(res) <- slotvec

  Misc <- lapply(dots, slot, "Misc")
  #names(Misc[[1]])

  Misc_identical <- function(x) all(vapply(x[-1], identical, logical(1), x[[1]]))

  Data <- do.call(c, lapply(Misc, getElement, "Data"))
  TryMP <- do.call(cbind, lapply(Misc, getElement, "TryMP"))

  Unfished <- lapply(Misc, getElement, "Unfished")
  Unfished_Refs <- lapply(Unfished, getElement, "Refs")
  stopifnot(Misc_identical(Unfished_Refs))

  Unfished_ByYear <- lapply(Unfished, getElement, "ByYear")
  stopifnot(Misc_identical(Unfished_ByYear))

  MSYRefs <- lapply(Misc, getElement, "MSYRefs")
  MSYRefs_Refs <- lapply(MSYRefs, getElement, "Refs")
  stopifnot(Misc_identical(MSYRefs_Refs))

  MSYRefs_ByYear <- lapply(MSYRefs, getElement, "ByYear")
  MSYRefs_ByYear2 <- lapply(names(MSYRefs_ByYear[[1]]),
                            function(x) do.call(abind::abind, c(lapply(MSYRefs_ByYear, getElement, x), along = 2)))
  names(MSYRefs_ByYear2) <- names(MSYRefs_ByYear[[1]])

  Misc_new <- list(Data = Data, TryMP = TryMP,
                   Unfished = list(Refs = Unfished_Refs[[1]], ByYear = Unfished_ByYear[[1]]),
                   MSYRefs = list(Refs = MSYRefs_Refs[[1]], ByYear = MSYRefs_ByYear2))

  slotvec_Misc <- c("LatEffort", "Revenue", "Cost", "TAE")
  Misc2 <- lapply(slotvec_Misc, function(x) do.call(abind::abind, c(lapply(Misc, getElement, x), along = 2)))
  names(Misc2) <- slotvec_Misc

  ## Create MSE Object ####
  MSEout <- new("MSE", Name = slots_identical("Name"), nyears = slots_identical("nyears"),
                proyears = slots_identical("proyears"), nMPs = length(slots_identical("MPs")),
                MPs = slots_identical("MPs"), nsim = slots_identical("nsim"),
                OM = dots[[1]]@OM, Obs = dots[[1]]@Obs, B_BMSY = res$B_BMSY, F_FMSY = res$F_FMSY, B = res$B, SSB = res$SSB,
                VB = res$VB, FM = res$FM, res$C, TAC = res$TAC, SSB_hist = dots[[1]]@SSB_hist, CB_hist = dots[[1]]@CB_hist,
                FM_hist = dots[[1]]@FM_hist, Effort = res$Effort, PAA = res$PAA, CAA = res$CAA, CAL = res$CAL,
                CALbins = slots_identical("CALbins"), Misc = c(Misc_new, Misc2))

  # Store MSE info
  attr(MSEout, "version") <- packageVersion("DLMtool")
  attr(MSEout, "date") <- date()
  attr(MSEout, "R.version") <- R.version

  MSEout
}
