library(dplyr)
library(reshape2)

#### Setup SRA data inputs
#SRA_data <- list(Year = 1918:2019)
#SRA_data$Len_bins <- seq(100, 850, 50)
#SRA_data$length_bin <- SRA_data$Len_bins[-length(SRA_data$Len_bins)] + 25
#SRA_data$maxage <- 80
#
## Data list
#MO <- readRDS("mse/data/yelloweye-rockfish-ins-privacy.rds")
#
####### Indices: HBLL, dogfish, commercial CPUE
#generate_indices <- function() {
#
#  Index <- I_sd <- matrix(NA, length(SRA_data$Year), 5)
#
#  #HBLL
#  HBLL <- readRDS("mse/data/hbll-joint-index.rds")
#  Index[, 1] <- vapply(SRA_data$Year, function(x) ifelse(!is.na(match(x, HBLL$year)), HBLL$est[match(x, HBLL$year)], NA), numeric(1))
#  I_sd[, 1] <- vapply(SRA_data$Year, function(x) ifelse(!is.na(match(x, HBLL$year)), HBLL$se[match(x, HBLL$year)], NA), numeric(1))
#
#  # Dogfish
#  dogfish <- readRDS("mse/data/dogfish-index.rds")
#  Index[, 2] <- vapply(SRA_data$Year, function(x) ifelse(!is.na(match(x, dogfish$year)), dogfish$est[match(x, dogfish$year)], NA), numeric(1))
#  I_sd[, 2] <- vapply(SRA_data$Year, function(x) ifelse(!is.na(match(x, dogfish$year)), dogfish$se[match(x, dogfish$year)], NA), numeric(1))
#
#  # Preliminary dogfish from survey (use selectivity from?)
#  Ind_2011 <- readxl::read_excel("mse/data/commercial_cpue.xlsx")
#
#  # Commerical CPUE with 2 breakpoints (3 series)
#  Index[, 3] <- vapply(SRA_data$Year, function(x) ifelse(!is.na(match(x, Ind_2011$Year)), Ind_2011$CC1[match(x, Ind_2011$Year)], NA), numeric(1))
#  I_sd[, 3] <- vapply(SRA_data$Year, function(x) ifelse(!is.na(match(x, Ind_2011$Year)), Ind_2011$CC1_SE[match(x, Ind_2011$Year)], NA), numeric(1))
#
#  Index[, 4] <- vapply(SRA_data$Year, function(x) ifelse(!is.na(match(x, Ind_2011$Year)), Ind_2011$CC2[match(x, Ind_2011$Year)], NA), numeric(1))
#  I_sd[, 4] <- vapply(SRA_data$Year, function(x) ifelse(!is.na(match(x, Ind_2011$Year)), Ind_2011$CC2_SE[match(x, Ind_2011$Year)], NA), numeric(1))
#
#  Index[, 5] <- vapply(SRA_data$Year, function(x) ifelse(!is.na(match(x, Ind_2011$Year)), Ind_2011$CC3[match(x, Ind_2011$Year)], NA), numeric(1))
#  I_sd[, 5] <- vapply(SRA_data$Year, function(x) ifelse(!is.na(match(x, Ind_2011$Year)), Ind_2011$CC3_SE[match(x, Ind_2011$Year)], NA), numeric(1))
#
#
#  SRA_data$I_type <<- c("est", "est", 1, 1, 1)
#  SRA_data$Index <<- Index
#  SRA_data$I_sd <<- I_sd
#
#  ##### Estimate HBLL selectivity and match it to the dogfish survey
#  map_s_vul_par <- matrix(NA, 3, 5)
#  map_s_vul_par[1:2, 1:2] <- c(1, 2, 1, 2)
#
#  SRA_data$map_s_vul_par <<- map_s_vul_par
#}
#generate_indices()
#
#
#################### Data
###### Rec lengths
#generate_CAL <- function() {
#
#  res <- readxl::read_excel("mse/data/SC Rec YellowEye Biodata.xlsx", sheet = 2) %>% filter(!is.na(`LENGTH(MM)`))
#  res$Len_bin <- cut(res$`LENGTH(MM)`, breaks = SRA_data$Len_bins, labels = SRA_data$Len_bins[-length(SRA_data$Len_bins)])
#  rec_len_matrix <- summarise(group_by(res, YEAR, Len_bin), n = n()) %>% acast(list("YEAR", "Len_bin"), value.var = "n", fill = 0)
#
#  Len_bins_avail <- colnames(rec_len_matrix) %>% as.numeric()
#  Nareas <- summarise(group_by(res, YEAR), n = length(unique(SUBAREA)))
#
#  #plot_composition(as.numeric(rownames(rec_len_matrix)), rec_len_matrix, CAL_bins = Len_bins_avail)
#  #plot(Len_bins_avail, colSums(rec_len_matrix), typ = 'o')
#
#  rec_CAL_fn <- function(x) {
#    ind <- match(x, as.numeric(rownames(rec_len_matrix)))
#    out <- rep(NA, length(SRA_data$Len_bins)-1) %>% as.numeric()
#    if(!is.na(ind)) {
#      L_ind <- match(Len_bins_avail, SRA_data$Len_bins)
#      out[L_ind] <- rec_len_matrix[ind, ] * Nareas$n[ind]/sum(rec_len_matrix[ind, ])
#    }
#    return(out)
#  }
#
#  CAL <- array(dim = c(length(SRA_data$Year), length(SRA_data$Len_bins)-1, 2))
#  CAL[,,2] <- do.call(rbind, lapply(SRA_data$Year, rec_CAL_fn))
#  #plot_composition(1918:2018, SRA_data$CAL[,,2], CAL_bins = SRA_data$Len_bins[-length(SRA_data$Len_bins)])
#
#  SRA_data$CAL <<- CAL
#}
#generate_CAL()
#
####### Commerical data: handline, trawl, comps
##HL <- readRDS("mse/data/ye-ins-catch-hl.rds")
##TR <- readRDS("mse/data/ye-ins-catch-tr.rds")
##MO$catch$total <- MO$catch$sum_landed_kg + MO$catch$sum_discarded_kg
##comm <- summarise(group_by(MO$catch, year, gear), total = sum(total)) %>% right_join(data.frame(year = 1969:2019)) %>% acast(list("year", "gear"))
##write.csv(1e-3 * rowSums(cbind(HL, TR), na.rm = TRUE), file = "mse/data/comm.csv")
#
#
####### Catch data. Fleet 1 = commerical (Hook and line plus trawl). Fleet 2 = Rec (FSC should be included in both fleets)
#generate_catch <- function() {
#  cat <- read.csv("mse/data/catch.csv")
#
#  Chist <- matrix(NA, length(SRA_data$Year), 2)
#  Chist[, 1] <- vapply(SRA_data$Year, function(x) ifelse(!is.na(match(x, cat$Year)), cat$Comm[match(x, cat$Year)], NA), numeric(1))
#  Chist[, 2] <- vapply(SRA_data$Year, function(x) ifelse(!is.na(match(x, cat$Year)), cat$Rec[match(x, cat$Year)], NA), numeric(1))
#
#  SRA_data$Chist <<- Chist
#}
#generate_catch()
#
### Double the commercial catch from 1986-2005
#double_commercial_catch <- function() {
#  ind <- match(1986:2005, SRA_data$Year)
#  SRA_data$Chist[ind, 1] <<- 2 * SRA_data$Chist[ind, 1]
#}
#double_commercial_catch()
#
### Comm age comps - 1989 HL
#generate_CAA <- function() {
#  comm_age <- filter(MO[[1]], !is.na(age) & gear_desc == "LONGLINE") %>% mutate(age2 = ifelse(age >= SRA_data$maxage, SRA_data$maxage, age)) %>%
#    group_by(year, age2) %>% summarise(n = n())
#
#  n_obs <- comm_age %>% group_by(year) %>% summarise(n = sum(n))
#
#  comm_age <- comm_age %>% full_join(data.frame(age2 = 1:SRA_data$maxage)) %>% acast(list("year", "age2"), value.var = "n", fill = 0)
#
#  n_fishing_events <- filter(MO[[1]], !is.na(age) & gear_desc == "LONGLINE")  %>%
#    group_by(year) %>% summarise(ntrips = length(unique(fishing_event_id))) %>% full_join(n_obs)
#
#  CAA <- array(dim = c(length(SRA_data$Year), SRA_data$maxage, 2))
#
#  CAA_fn <- function(x, mat, ntrips) {
#    ind <- match(x, as.numeric(rownames(mat)))
#    out <- rep(NA, SRA_data$maxage) %>% as.numeric()
#    if(!is.na(ind)) {
#      out[as.numeric(colnames(mat))] <- mat[ind, ]
#
#      if(sum(out, na.rm = TRUE) > 0) {
#        Nt <- min(ntrips$ntrips[match(x, ntrips$year)], ntrips$n[match(x, ntrips$year)])
#        out <- out * Nt/sum(out, na.rm = TRUE)
#      }
#    }
#    return(out)
#  }
#  CAA[,,1] <- do.call(rbind, lapply(SRA_data$Year, CAA_fn, mat = comm_age, ntrips = n_fishing_events))
#  SRA_data$CAA <<- CAA
#}
#generate_CAA()
#
### Survey age comps
#generate_survey_CAA <- function() {
#  s_age <- filter(MO[[2]], !is.na(age) & (survey_abbrev == "HBLL INS N" | survey_abbrev == "HBLL INS S")) %>%
#    mutate(age2 = ifelse(age >= SRA_data$maxage, SRA_data$maxage, age)) %>% group_by(year, age2) %>%
#    summarise(n = n())
#
#  n_obs <- s_age %>% group_by(year) %>% summarise(n = sum(n))
#  s_age <- s_age %>% full_join(data.frame(age2 = 1:SRA_data$maxage)) %>% acast(list("year", "age2"), value.var = "n", fill = 0)
#
#  n_fishing_events <- filter(MO[[2]], !is.na(age) & (survey_abbrev == "HBLL INS N" | survey_abbrev == "HBLL INS S"))  %>%
#    group_by(year) %>% summarise(ntrips = length(unique(fishing_event_id))) %>% full_join(n_obs)
#
#  s_CAA <- array(dim = c(length(SRA_data$Year), SRA_data$maxage, 5))
#
#  CAA_fn <- function(x, mat, ntrips) {
#    ind <- match(x, as.numeric(rownames(mat)))
#    out <- rep(NA, SRA_data$maxage) %>% as.numeric()
#    if(!is.na(ind)) {
#      out[as.numeric(colnames(mat))] <- mat[ind, ]
#
#      if(sum(out, na.rm = TRUE) > 0) {
#        Nt <- ntrips$ntrips[match(x, ntrips$year)]
#        out <- out * Nt/sum(out, na.rm = TRUE)
#      }
#    }
#    return(out)
#  }
#
#  s_CAA[,,1] <- do.call(rbind, lapply(SRA_data$Year, CAA_fn, mat = s_age, ntrips = n_fishing_events))
#
#  SRA_data$s_CAA <<- s_CAA
#}
#generate_survey_CAA()
#
## Selectivity parameters for fishery
#Linf <- 65.1
#K <- 0.04
#t0 <- -9.04
#ages <- c(10, 20, 5, 10) # ages of 5% and full selectivity for commercial, then for rec.
#Linf * (1 - exp(-K * (ages - t0)))
#SRA_data$vul_par <- matrix(c(45, 35, 0.99, 35, 28, 0.99), 3, 2)
#
## Fix rec devs - NA = fixed to zero for that year, number = estimate
#SRA_data$map_log_rec_dev <- c(rep(NA, 30), 1:67, rep(NA, 5))
#
## For growth and maturity, see growth_maturity.R (gfsynopsis code)
#SRA_data$f_name <- c("Commercial", "Recreational")
#SRA_data$s_name <- c("HBLL", "Dogfish", "Com_CPUE", "Com_CPUE", "Com_CPUE")
#
#saveRDS(SRA_data, file = "mse/scoping/SRA_data.rds")
