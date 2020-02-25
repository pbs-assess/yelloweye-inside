

##### Constant catch MPs of 0, 5, 10, or 15 tonnes
NF <- function(x, Data, reps) {
  Rec <- new("Rec")
  Rec@TAC <- rep(1e-8, reps)
  return(Rec)
}
class(NF) <- "MP"

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
IDX_YE <- YE_wrapper(gfdlm::IDX, tac_floor = 5)
IDX_smooth_YE <- YE_wrapper(gfdlm::IDX_smooth, tac_floor = 5)

# Iratio
Iratio_23 <- YE_wrapper(Iratio)
Iratio_510 <- YE_wrapper(Iratio, yrs = c(5, 10))

# I slope
GB_slope_lambda1 <- YE_wrapper(GB_slope)
GB_slope_lambda05 <- YE_wrapper(GB_slope, lambda = 0.5) # Slower changes in catch relative to index
GB_slope_yrs10 <- YE_wrapper(GB_slope, yrsmth = 10) # Longer time window for smoothing

Islope_5_lambda04 <- YE_wrapper(Islope1, xx = 0)
Islope_10_lambda04 <- YE_wrapper(Islope1, xx = 0, yrsmth = 10)
Islope_10_lambda08 <- YE_wrapper(Islope1, xx = 0, yrsmth = 10, lambda = 0.8)

# Ratio of mean index to Iref = 2 * mean 2010-2019 index
IT_YE <- function(x, Data, reps = 100, plot=FALSE, yrsmth = 5, mc = 0.05) {
  dependencies = "Data@Ind, Data@MPrec, Data@CV_Ind, Data@Iref"
  ind <- max(1, (length(Data@Year) - yrsmth + 1)):length(Data@Year)
  if(is.na(Data@Iref[x])) return(list(TAC=rep(as.numeric(NA), reps)))
  deltaI <- mean(Data@Ind[x, ind], na.rm=TRUE)/Data@Iref[x]
  if (deltaI < (1 - mc)) deltaI <- 1 - mc
  if (deltaI > (1 + mc)) deltaI <- 1 + mc
  TAC <- Data@MPrec[x] * deltaI * trlnorm(reps, 1, if(is.na(Data@CV_Ind[x, 1])) 0.1 else Data@CV_Ind[x, 1])
  TAC <- TACfilter(TAC)
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(op)
    par(mfrow=c(1,2), oma=c(1,1,1,1), mar=c(5,4,1,4))
    ylim <- range(c(Data@Ind[x,ind], Data@Iref[x]))
    plot(Data@Year[ind], Data@Ind[x,ind], xlab="Year",
         ylab= paste0("Index (previous ", yrsmth, "years)"), bty="l", type="l",
         lwd=2, ylim=ylim)
    lines(Data@Year[ind], rep(mean(Data@Ind[x, ind], na.rm=TRUE), length(ind)), lty=2)
    text(quantile(Data@Year[ind],0.15), mean(Data@Ind[x, ind], na.rm=TRUE), "Mean Index", pos=3)
    lines(Data@Year[ind], rep(mean(Data@Iref[x], na.rm=TRUE), length(ind)), lty=3)
    text(quantile(Data@Year[ind],0.15), Data@Iref[x], "Reference Index", pos=3)

    boxplot(TAC, ylab=paste0("TAC (", Data@Units, ")"))
    points(1, Data@MPrec[x], cex=2, col="blue", pch=16)
    text(1, Data@MPrec[x], cex=1, col="blue", "Last Catch", pos=1)
  }
  Rec <- new("Rec")
  Rec@TAC <- TAC
  return(Rec)
}
class(IT_YE) <- "MP"

IT5_mc05 <- YE_wrapper(IT_YE) #
IT5_mc025 <- YE_wrapper(IT_YE, mc = 0.025)

IT10_mc05 <- YE_wrapper(IT_YE, mc = 0.05)
IT10_mc025 <- YE_wrapper(IT_YE, mc = 0.025)

# Iref = 1.5 * mean 2010-2015 index
Itarget_5 <- YE_wrapper(Itarget1)
Itarget_10 <- YE_wrapper(Itarget1, yrsmth = 10)

##### SP model that only uses future HBLL and with r-rprio
SP_YE <- function(x, Data, ...) {
  nyears <- ncol(Data@Cat)
  if(nyears > 102) Data@AddInd[, 2:5, 103:nyears] <- Data@CV_AddInd[, 2:5, 103:nyears] <- NA
  SP_SS(x = x, Data = Data, AddInd = 1:5, start = list(r_prior = c(0.068, 0.03)), ...)
}
class(SP_YE) <- "Assess"

##### SP with 40-80 HCR
SP_4080_ <- function(x, Data, reps = 1) {
  do_Assessment <- SP_YE(x = x, Data = Data)
  Rec <- HCR_ramp(Assessment = do_Assessment, reps = reps,
                  LRP = 0.4, TRP = 0.8, RP_type = "SSB_SSBMSY")
  Rec@Misc <- MSEtool:::Assess_diagnostic(x, Data, do_Assessment, include_assessment = FALSE)
  return(Rec)
}
class(SP_4080_) <- "MP"
SP_4080 <- YE_wrapper(SP_4080_)
SP_4080_5f <- SP_4080_10f <- SP_4080

##### SP with 20-60 HCR
SP_2060_ <- function(x, Data, reps = 1) {
  do_Assessment <- SP_YE(x = x, Data = Data)
  Rec <- HCR_ramp(Assessment = do_Assessment, reps = reps,
                  LRP = 0.2, TRP = 0.6, RP_type = "SSB_SSBMSY")
  Rec@Misc <- MSEtool:::Assess_diagnostic(x, Data, do_Assessment, include_assessment = FALSE)
  return(Rec)
}
class(SP_2060_) <- "MP"
SP_2060 <- YE_wrapper(SP_2060_)
SP_2060_5f <- SP_2060_10f <- SP_2060

##### Interim MP with SP
interim_MP <- eval(bquote(function(x, Data, reps = 1, assess, assessment_interval, assess_args = list(), HCR = HCR_MSY, HCR_args = list(),
                                   I_smooth = c("none", "loess", "mean", "buffer"), smooth_par = NULL, use_ramp = FALSE, cap_TAC = FALSE) {

  #dependencies <- .(MSEtool:::get_dependencies("SP"))
  I_smooth <- match.arg(I_smooth)

  current_yr <- Data@Year[length(Data@Year)]
  run_assess <- current_yr == Data@LHYear
  if (current_yr > Data@LHYear) {
    run_assess <- current_yr == Data@Misc[[x]]$next_assess_yr
  }
  if (!run_assess) run_interim_MP <- TRUE else run_interim_MP <- FALSE

  if (run_assess) {
    # Return TAC = UMSY * VB_current when run_SCA = TRUE
    assess_formals <- list(x = x, Data = Data)
    do_Assessment <- do.call(assess, c(assess_formals, assess_args))
    Assess_output <- MSEtool:::Assess_diagnostic(x, Data, do_Assessment, include_assessment = FALSE)

    if (do_Assessment@conv) {
      Rec <- do.call(match.fun(HCR), c(Assessment = do_Assessment, reps = reps, HCR_args))

      # Set-up references for interim MP
      q <- as.numeric(do_Assessment@SD$value[names(do_Assessment@SD$value) == "q"][1])
      sigma_buffer <- sd(log(do_Assessment@Obs_Index[, 1]/do_Assessment@Index[, 1]), na.rm = TRUE)
      buffer_quantities <- structure(c(Rec@TAC, do_Assessment@Index[nrow(do_Assessment@Index), 1], sigma_buffer),
                                     names = c("C_ref", "I_ref", "sigma_buffer"))
      Rec@Misc <- c(list(I = Data@Ind[x, ], q = q, buffer_quantities = buffer_quantities, Last_Assess = current_yr,
                         UMSY = 1 - exp(-do_Assessment@FMSY), MSY = do_Assessment@MSY, SSB_dep = do_Assessment@SSB_SSB0[length(do_Assessment@SSB_SSB0)],
                         next_assess_yr = current_yr + assessment_interval), Assess_output)

      if (cap_TAC) {
        Rec@TAC <- max(Rec@TAC, Rec@Misc$MSY)
      }
      run_interim_MP <- FALSE

    } else {

      if (current_yr == Data@LHYear || length(Data@Misc[[x]]) == 2) {
        Rec <- new("Rec")
        Rec@TAC <- TACfilter(rep(NA, reps))
        Rec@Misc <- c(list(next_assess_yr = current_yr + 1), Assess_output)

        run_interim_MP <- FALSE
      } else {
        run_interim_MP <- TRUE
        next_assess_yr <- current_yr + 1
      }
    }
  }

  if (run_interim_MP) {
    # Estimate new_VB as new_I/q, then new TAC = UMSY * new_VB - also equivalent: TAC = MSY * I_y / I_MSY
    q_ratio <- Data@Ind[x, 1]/Data@Misc[[x]]$I[1]
    q_update <- Data@Misc[[x]]$q * q_ratio

    if (I_smooth == "none") new_Index <- Data@Ind[x, length(Data@Ind[x, ])]
    if (I_smooth == "loess") {
      I_df <- data.frame(Year = Data@Year, Ind = Data@Ind[x, ])
      fit <- loess(Ind ~ Year, I_df)
      new_Index <- fit$fitted[length(fit$fitted)]
    }

    if (I_smooth == "mean") {
      nyr <- smooth_par[1]
      new_Index <- mean(Data@Ind[x, (length(Data@Ind[x, ])-nyr+1):length(Data@Ind[x, ])], na.rm = TRUE)
    }

    if (I_smooth == "buffer") {
      new_Index <- Data@Ind[x, length(Data@Ind[x, ])]
      new_Iref <- q_ratio * Data@Misc[[x]]$buffer_quantities["I_ref"]
      sd_buffer <- Data@Misc[[x]]$buffer_quantities["sigma_buffer"]
      new_TAC <- Data@Misc[[x]]$buffer_quantities["C_ref"] * (new_Index + smooth_par * sd_buffer)/
        (new_Iref + smooth_par * sd_buffer)
    } else {
      VB_curr <- as.numeric(new_Index/q_update)
      new_TAC <- Data@Misc[[x]]$UMSY * VB_curr
    }

    if (cap_TAC) {
      TAC_used <- max(new_TAC, Data@Misc[[x]]$MSY)
    } else TAC_used <- new_TAC

    if (use_ramp) { #40-10 control rule based on previous assessment
      TAC_used <- HCRlin(Data@Misc[[x]]$SSB_dep, 0.1, 0.4) * TAC_used
    }

    Rec <- new("Rec")
    if(is.infinite(TAC_used) || is.na(TAC_used)) stop("Error in TAC during interim")
    Rec@TAC <- TACfilter(TAC_used)
    if (run_assess) {
      Rec@Misc <- c(list(I = Data@Ind[x, ], q = q_update, buffer_quantities = Data@Misc[[x]]$buffer_quantities,
                         Last_Assess = Data@Misc[[x]]$Last_Assess,
                         UMSY = Data@Misc[[x]]$UMSY, MSY = Data@Misc[[x]]$MSY,
                         SSB_dep = Data@Misc[[x]]$SSB_dep,
                         next_assess_yr = ifelse(exists("next_assess_yr"), next_assess_yr, Data@Misc[[x]]$next_assess_yr)),
                    Assess_output)

    } else {

      Rec@Misc <- list(I = Data@Ind[x, ], q = q_update, buffer_quantities = Data@Misc[[x]]$buffer_quantities,
                       Last_Assess = Data@Misc[[x]]$Last_Assess,
                       UMSY = Data@Misc[[x]]$UMSY, MSY = Data@Misc[[x]]$MSY,
                       SSB_dep = Data@Misc[[x]]$SSB_dep,
                       next_assess_yr = ifelse(exists("next_assess_yr"), next_assess_yr, Data@Misc[[x]]$next_assess_yr),
                       diagnostic = Data@Misc[[x]]$diagnostic)
    }
  }
  return(Rec)
}))
class(interim_MP) <- "MP"
environment(interim_MP) <- asNamespace("MSEtool")


# Function to make MP with variable options
make_iMP <- function(...) {
  fn <- interim_MP
  dots <- list(...)
  arg_ind <- pmatch(names(dots), names(formals(fn)))
  formals(fn)[arg_ind] <- dots
  class(fn) <- "MP"
  return(fn)
}

SP_interim <- make_iMP(assess = SP_YE, assessment_interval = 10, I_smooth = "buffer", smooth_par = 1)



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
  MSYRefs_ByYear2 <- lapply(c("MSY", "FMSY", "SSBMSY", "BMSY", "VBMSY"),
                            function(x) do.call(abind::abind, c(lapply(MSYRefs_ByYear, getElement, x), along = 2)))

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
