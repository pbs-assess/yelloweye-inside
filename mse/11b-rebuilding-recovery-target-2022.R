
FRENCH <- FALSE

library(tidyverse)
library(rosettafish)
library(gfdlm)

if(FRENCH) {
  options(french = TRUE)
}

sc <- readRDS("mse/om/ye-scenarios2.rds")

sc_french <- c(
  "(1) Base",
  "(2) Faibles prises",
  "(3) Recrutement\népisodique",
  "(4) Estimation de la\nsélectivité du RPFD",
  "(A) Faible M ",
  "(B) CV élevé\ndu RPFD")

if (FRENCH) sc$scenario_human <- sc_french

get_filtered_scenario <- function(type, column) {
  filter(sc, scenario_type == type) %>%
    pull(!!column) %>%
    set_names()
}
scenarios <- sc$scenario %>% set_names()
scenario_human <- sc$scenario_human %>% set_names()
scenarios_ref <- get_filtered_scenario("Reference", "scenario")
scenarios_ref_human <- get_filtered_scenario("Reference", "scenario_human")
scenarios_rob <- get_filtered_scenario("Robustness", "scenario")
scenarios_rob_human <- get_filtered_scenario("Robustness", "scenario_human")

# Read MSE objects --------------------------------------------------------------------
targ_proj <- function(x, y, BMSYtarget = seq(0.4, 1.6, 0.2), MP, ymax = 100) {

  MPind <- match(MP, x@MPs)
  lapply(BMSYtarget, function(xx) {
    x@B_BMSY[, MPind, 1:ymax, drop = FALSE] %>%
      apply(c(2, 3), function(x) mean(x >= xx)) %>%
      structure(dimnames = list(MP = x@MPs[MPind], Year = x@OM$CurrentYr[1] + 1:ymax)) %>%
      reshape2::melt() %>%
      mutate(OM = y,
             target_num = xx,
             target = factor(paste0(100 * xx, "~\"%\"~B[", en2fr("MSY", FRENCH), "]"),
                             levels = paste0(100 * BMSYtarget, "~\"%\"~B[", en2fr("MSY", FRENCH), "]")),
             target2 = factor(paste0(100 * xx, "%\n", en2fr("BMSY", FRENCH)),
                              levels = paste0(100 * BMSYtarget, "%\n", en2fr("BMSY", FRENCH))))
  }) %>% bind_rows()
}
parser <- function(x) parse(text = x)

mse_fmsyproj <- map(scenarios, ~ readRDS(paste0("mse/om/MSE_", .x, "_FMSYproj.rds")))
mse <- map(scenarios, ~ readRDS(paste0("mse/om/MSE_", .x, ".rds")))

# Rebuilding target --------------------------------------------------------------------
.ggsave <- function(filename, plot, ...) {
  if(FRENCH) {
    filename <- paste0("figs_2022/fr/", filename)
  } else {
    filename <- paste0("figs_2022/", filename)
  }
  ggsave(filename, plot, ...)
}

# Annual probability of being above X% MSY with no fishing
pNFref <- Map(targ_proj, x = mse, y = scenario_human, MoreArgs = list(MP = "NFref")) %>% bind_rows()

# Annual probability of being above X% MSY with 15 t TAC
p15t <- Map(targ_proj, x = mse, y = scenario_human, MoreArgs = list(MP = "CC_15t")) %>% bind_rows()

# Annual probability of being above X% MSY with FMSY fishing
pFMSY <- Map(targ_proj, x = mse_fmsyproj, y = scenario_human, MoreArgs = list(MP = "MP_FMSY")) %>% bind_rows()

# Summary table of probabilities after 1.5 generations
MPout <- c("NFref" = en2fr("No fishing", FRENCH, custom_terms = data.frame(english = "No fishing", french = "Aucun pêche")),
           "CC_15t" = en2fr("CC_15t", FRENCH, custom_terms = data.frame(english = "CC_15t", french = "PC_15t")),
           "MP_FMSY" = en2fr("FMSY", FRENCH))

tig <- local({
  tt <- rbind(pNFref, p15t, pFMSY) %>% filter(Year == 2019 + 56) %>%
    mutate(MP = MPout[match(MP, names(MPout))]) %>%
    select(!"target" & !"Year")

  OM_names <- unique(tt$OM)
  tt_list <- lapply(OM_names, function(x) {
    dplyr::filter(tt, OM == x) %>% reshape2::dcast(MP ~ target2)
  }) %>% structure(names = OM_names)

  g <- plot_tigure_facet(tt_list, mp_order = MPout %>% rev()) +
    theme(axis.text.x = element_text(size = 8))
  g$facet$params$ncol <- 2
  g
})
.ggsave("rebuilding_table.png", tig, height = 5, width = 6)


tt_avg <- local({
  tt <- rbind(pNFref, p15t, pFMSY) %>% filter(Year == 2019 + 56) %>%
    mutate(MP = MPout[match(MP, names(MPout))]) %>%
    select(!"target" & !"Year")

  OM_names <- unique(tt$OM)[1:4]

  tt_avg <- tt %>% filter(OM %in% OM_names) %>% group_by(MP, target2) %>% summarise(value = mean(value)) %>%
    reshape2::dcast(MP ~ target2)

  g <- plot_tigure(tt_avg, mp_order = MPout %>% rev()) +
    theme(axis.text.x = element_text(size = 8))
  g
})
.ggsave("rebuilding_table_avg.png", tt_avg, height = 2, width = 4)




# Figures
y_trans <- en2fr("Probability above target", FRENCH,
                 custom_terms = data.frame(english = "Probability above target",
                                           french = "Probabilité en haut de cible"))
col_trans <- en2fr("Candidate\nrebuilding targets", FRENCH,
                   custom_terms = data.frame(english = "Candidate\nrebuilding targets",
                                             french = "Cibles de\nreconstruction\npotentielles"))
g <- ggplot(pNFref, aes(Year, value, colour = target)) +
  geom_path() +
  gfdlm::theme_pbs() +
  facet_wrap(~ OM) +
  expand_limits(y = 0) + coord_cartesian(xlim = 2019 + c(1, 56)) +
  labs(x = en2fr("Year", FRENCH), y = y_trans, colour = col_trans) +
  scale_color_viridis_d(labels = parser)
.ggsave("projection_rebuilding_target_NFref.png", g, height = 4, width = 8)

g <- ggplot(p15t, aes(Year, value, colour = target)) +
  geom_path() +
  gfdlm::theme_pbs() +
  facet_wrap(~ OM) +
  expand_limits(y = 0) + coord_cartesian(xlim = 2019 + c(1, 56)) +
  labs(x = en2fr("Year", FRENCH), y = y_trans, colour = col_trans) +
  scale_color_viridis_d(labels = parser)
.ggsave("/projection_rebuilding_target_15t.png", g, height = 4, width = 8)

g <- ggplot(pFMSY, aes(Year, value, colour = target)) +
  geom_path() +
  gfdlm::theme_pbs() +
  facet_wrap(~ OM) +
  expand_limits(y = 0) + coord_cartesian(xlim = 2019 + c(1, 56)) +
  labs(x = en2fr("Year", FRENCH), y = y_trans, colour = col_trans) +
  scale_color_viridis_d(labels = parser)
.ggsave("projection_rebuilding_target_FMSY.png", g, height = 4, width = 8)


# Recovery target --------------------------------------------------------------------
Hist_SSB <- sapply(sc$scenario, function(x) readRDS(paste0("mse/om/", x, ".rds"))@SSB, simplify = "array")

cosewic_proj <- function(i, x, y, threshold = c(0.3, 0.5, 0.7), MP, MGT = 38, ymax = 100, Hist) {

  MPind <- match(MP, x[[i]]@MPs)
  B_hist <- Hist[, 1:102, i]
  B_proj <- x[[i]]@SSB[, MPind, ]
  p <- sapply(1:ymax, function(yy) {
    arr <- cbind(B_hist, B_proj[, 1:yy, drop = FALSE])
    yend <- ncol(arr)
    ystart <- max(1, ncol(arr) - 3 * MGT + 1)

    metric <- 1 - arr[, yend]/arr[, ystart]
    sapply(threshold, function(Ref) sum(metric > Ref)/length(metric))

  }) %>% structure(dimnames = list(threshold = threshold, Year = x[[i]]@OM$CurrentYr[1] + 1:ymax)) %>%
    reshape2::melt(value.name = "tvalue") %>% mutate(OM = y[i], MP = x[[i]]@MPs[MPind])

  return(p)
}

# Probability of decline during projection
rNFref <- lapply(1:length(mse), cosewic_proj, x = mse, y = scenario_human, MP = "NFref", Hist = Hist_SSB) %>% bind_rows()

r15t <- lapply(1:length(mse), cosewic_proj, x = mse, y = scenario_human, MP = "CC_15t", Hist = Hist_SSB) %>% bind_rows()

rFMSY <- lapply(1:length(mse_fmsyproj), cosewic_proj, x = mse_fmsyproj,
                y = scenario_human, MP = "MP_FMSY", Hist = Hist_SSB) %>% bind_rows()

col_trans <- en2fr("COSEWIC\nthreshold", FRENCH,
                   custom_terms = data.frame(english = "COSEWIC\nthreshold",
                                             french = "Seuil du\nCOSEPAC"))
y_trans <- en2fr("Probability of X % decline", FRENCH,
                 custom_terms = data.frame(english = "Probability of X % decline",
                                           french = "Probabilité de déclin (X%)"))
g <- ggplot(rNFref, aes(Year, tvalue, colour = paste0(100 * threshold, "%"))) +
  geom_path() +
  #geom_point() +
  gfdlm::theme_pbs() +
  facet_wrap(~ OM) +
  expand_limits(y = 0) +
  geom_hline(yintercept = 0.5, linetype = 3) + geom_vline(xintercept = 2019 + 56, linetype = 3) +
  labs(x = en2fr("Year", FRENCH), y = y_trans, colour = paste(col_trans, "(X%)"))
.ggsave("projection_recovery_target_NFref.png", g, height = 4, width = 8)

g <- ggplot(r15t, aes(Year, tvalue, colour = paste0(100 * threshold, "%"))) +
  geom_path() +
  #geom_point() +
  gfdlm::theme_pbs() +
  facet_wrap(~ OM) +
  expand_limits(y = 0) +
  geom_hline(yintercept = 0.5, linetype = 3) + geom_vline(xintercept = 2019 + 56, linetype = 3) +
  labs(x = en2fr("Year", FRENCH), y = y_trans, colour = paste(col_trans, "(X%)"))
.ggsave("projection_recovery_target_15t.png", g, height = 4, width = 8)


g <- ggplot(rFMSY, aes(Year, tvalue, colour = paste0(100 * threshold, "%"))) +
  geom_path() +
  #geom_point() +
  gfdlm::theme_pbs() +
  facet_wrap(~ OM) +
  expand_limits(y = 0) +
  geom_hline(yintercept = 0.5, linetype = 3) + geom_vline(xintercept = 2019 + 56, linetype = 3) +
  labs(x = en2fr("Year", FRENCH), y = y_trans, colour = paste(col_trans, "(X%)"))
.ggsave("projection_recovery_target_FMSY.png", g, height = 4, width = 8)

# Summary table of probabilities after 1.5 generations
tig <- local({
  tt <- rbind(rNFref, r15t, rFMSY) %>% filter(Year == 2019 + 56) %>%
    mutate(MP = MPout[match(MP, names(MPout))],
           threshold = paste0(threshold * 100, "%\n", en2fr("decline", FRENCH, case = "lower"))) %>%
    select(!"Year")

  OM_names <- unique(tt$OM)
  tt_list <- lapply(OM_names, function(x) {
    dplyr::filter(tt, OM == x) %>% reshape2::dcast(MP ~ threshold, value.var = "tvalue")
  }) %>% structure(names = OM_names)
  tt_list
  g <- plot_tigure_facet(tt_list, mp_order = MPout %>% rev()) +
    theme(axis.text.x = element_text(size = 8)) +
    scale_fill_viridis_c(limits = c(0, 1), begin = 0.15, end = 1, alpha = 0.6, option = "D", direction = -1)
  g
})
.ggsave("recovery_table.png", tig, height = 3.5, width = 6)

tt_avg <- local({
  tt <- rbind(rNFref, r15t, rFMSY) %>% filter(Year == 2019 + 56) %>%
    mutate(MP = MPout[match(MP, names(MPout))],
           threshold = paste0(threshold * 100, "%\n", en2fr("decline", FRENCH, case = "lower"))) %>%
    select(!"Year")

  OM_names <- unique(tt$OM)[1:4]

  tt_avg <- tt %>% filter(OM %in% OM_names) %>% group_by(MP, threshold) %>% summarise(tvalue = mean(tvalue)) %>%
    reshape2::dcast(MP ~ threshold, value.var = "tvalue")

  g <- plot_tigure(tt_avg, mp_order = MPout %>% rev()) +
    theme(axis.text.x = element_text(size = 8)) +
    scale_fill_viridis_c(limits = c(0, 1), begin = 0.15, end = 1, alpha = 0.6, option = "D", direction = -1)
  g
})
.ggsave("recovery_table_avg.png", tt_avg, height = 2, width = 3)

# Recovery table after 100 years
tig_100 <- local({
  tt <- rbind(rNFref, r15t, rFMSY) %>% filter(Year == max(Year)) %>%
    mutate(MP = MPout[match(MP, names(MPout))],
           threshold = paste0(threshold * 100, "%\n", en2fr("decline", FRENCH, case = "lower"))) %>%
    select(!"Year")

  OM_names <- unique(tt$OM)
  tt_list <- lapply(OM_names, function(x) {
    dplyr::filter(tt, OM == x) %>% reshape2::dcast(MP ~ threshold, value.var = "tvalue")
  }) %>% structure(names = OM_names)
  tt_list
  g <- plot_tigure_facet(tt_list, mp_order = MPout %>% rev()) +
    theme(axis.text.x = element_text(size = 8))
  g
})
.ggsave("recovery_table_100.png", tig_100, height = 3.5, width = 6)

tt_avg_100 <- local({
  tt <- rbind(rNFref, r15t, rFMSY) %>% filter(Year == max(Year)) %>%
    mutate(MP = MPout[match(MP, names(MPout))],
           threshold = paste0(threshold * 100, "%\n", en2fr("decline", FRENCH, case = "lower"))) %>%
    select(!"Year")

  OM_names <- unique(tt$OM)[1:4]

  tt_avg <- tt %>% filter(OM %in% OM_names) %>% group_by(MP, threshold) %>% summarise(tvalue = mean(tvalue)) %>%
    reshape2::dcast(MP ~ threshold, value.var = "tvalue")

  g <- plot_tigure(tt_avg, mp_order = MPout %>% rev()) +
    theme(axis.text.x = element_text(size = 8))
  g
})
.ggsave("recovery_table_avg_100.png", tt_avg_100, height = 2, width = 3)


# Probability of decline vs. biomass
cNFref <- left_join(pNFref, rNFref, by = c("Year", "OM", "MP")) %>%
  mutate(target3 = paste0("Y = ", 100 * target_num, "% BMSY"), threshold3 = paste0(100 * threshold, "%"))

c15t <- left_join(p15t, r15t, by = c("Year", "OM", "MP")) %>%
  mutate(target3 = paste0("Y = ", 100 * target_num, "% BMSY"), threshold3 = paste0(100 * threshold, "%"))

cFMSY <- left_join(pFMSY, rFMSY, by = c("Year", "OM", "MP")) %>%
  mutate(target3 = paste0("Y = ", 100 * target_num, "% BMSY"), threshold3 = paste0(100 * threshold, "%"))


#dat_end <- cNFref %>% filter(Year %in% c(range(Year), 2019 + 56))
#g <- ggplot(cNFref, aes(tvalue, value, colour = threshold3, group = threshold3)) +
#  facet_grid(OM ~ target3, scales = "free_y") +
#  geom_vline(xintercept = 0.5, linetype = 3) +
#  geom_hline(yintercept = 0.5, linetype = 3) +
#  geom_path() +
#  geom_point(data = dat_end, aes(shape = as.factor(Year))) +
#  labs(y = "Probability above Y % BMSY", x = "Probability of X % decline",
#       shape = "Year", colour = "COSEWIC\nthreshold (X%)") +
#  theme_bw() +
#  scale_shape_manual(values = c(16, 8, 17)) +
#  theme(panel.spacing = unit(0, "lines"), axis.text.x = element_text(angle = 45, vjust = 0.5)) +
#  expand_limits(y = 0.4) +
#  ggtitle(expression("Projections with no fishing"))
#ggsave("figs_2022/projection_recovery_phase_NFref.png", g, height = 7.25, width = 7.5)
#
#
#dat_end <- c15t %>% filter(Year %in% c(range(Year), 2019 + 56))
#g <- ggplot(c15t, aes(tvalue, value, colour = threshold3, group = threshold3)) +
#  facet_grid(OM ~ target3, scales = "free_y") +
#  geom_vline(xintercept = 0.5, linetype = 3) +
#  geom_hline(yintercept = 0.5, linetype = 3) +
#  geom_path() +
#  geom_point(data = dat_end, aes(shape = as.factor(Year))) +
#  labs(y = "Probability above Y % BMSY", x = "Probability of X % decline",
#       shape = "Year", colour = "COSEWIC\nthreshold (X%)") +
#  theme_bw() +
#  scale_shape_manual(values = c(16, 8, 17)) +
#  theme(panel.spacing = unit(0, "lines"), axis.text.x = element_text(angle = 45, vjust = 0.5)) +
#  expand_limits(y = 0.4) +
#  ggtitle(expression("Projections with 15 t catch"))
#ggsave("figs_2022/projection_recovery_phase_15t.png", g, height = 7.25, width = 7.5)
#
#dat_end <- cFMSY %>% filter(Year %in% c(range(Year), 2019 + 56))
#g <- ggplot(cFMSY, aes(tvalue, value, colour = threshold3, group = threshold3)) +
#  facet_grid(OM ~ target3, scales = "free_y") +
#  geom_vline(xintercept = 0.5, linetype = 3) +
#  geom_hline(yintercept = 0.5, linetype = 3) +
#  geom_path() +
#  geom_point(data = dat_end, aes(shape = as.factor(Year))) +
#  labs(y = "Probability above Y % BMSY", x = "Probability of X % decline",
#       shape = "Year", colour = "COSEWIC\nthreshold (X%)") +
#  theme_bw() +
#  scale_shape_manual(values = c(16, 8, 17)) +
#  theme(panel.spacing = unit(0, "lines"), axis.text.x = element_text(angle = 45, vjust = 0.5)) +
#  expand_limits(y = 0.4) +
#  ggtitle(expression("Projections with F ="~F[MSY]))
#ggsave("figs_2022/projection_recovery_phase_FMSY.png", g, height = 7.25, width = 7.5)


# Scatter plot
scat_plot <- rbind(cNFref, c15t, cFMSY) %>% filter(Year == 2019 + 56, MP == "CC_15t") %>%
  mutate(target = paste0(target_num * 100, "% B", en2fr("MSY", FRENCH)) %>%
           factor(levels = paste0(seq(40, 160, 20), "% B", en2fr("MSY", FRENCH))))

col_trans <- en2fr("Candidate\ntargets", FRENCH,
                   custom_terms = data.frame(english = "Candidate\ntargets",
                                             french = "Cibles\npotentielles"))
shape_trans <- en2fr("COSEWIC\nthreshold", FRENCH,
                   custom_terms = data.frame(english = "COSEWIC\nthreshold",
                                             french = "Seuil du\nCOSEPAC"))
x_trans <- en2fr("Probability of decline", FRENCH,
                 custom_terms = data.frame(english = "Probability of decline",
                                               french = "Probabilité de déclin"))
y_trans <- en2fr("Probability above target", FRENCH,
                 custom_terms = data.frame(english = "Probability above target",
                                           french = "Probabilité en haut de cible"))

g <- ggplot(scat_plot, aes(tvalue, value, group = target, colour = target, shape = threshold3)) +
  geom_path() +
  geom_point() +
  theme_bw() + facet_wrap(~ OM) +
  labs(colour = col_trans,
       shape = shape_trans,
       x = x_trans,
       y = y_trans) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1))
.ggsave("projection_recovery_scatter.png", g, height = 4, width = 7.5)

