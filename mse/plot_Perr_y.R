library("dplyr")
library("DLMtool")
library("MSEtool")
library("here")
library("purrr")
library("ggplot2")
library("reshape2")
library("cowplot")
library("future")
cores <- floor(future::availableCores()/2)
plan(multisession, workers = cores)

sp <- "ye"

sc <- readRDS("mse/om/ye-scenarios2.rds")
fig_dir <- "mse/figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir)

.ggsave <- function(filename, width, height, ...) {
  ggsave(file.path(fig_dir, paste0(sp, "-", filename, ".png")),
         width = width, height = height, ...
  )
}

sra_ye <- lapply(sc$scenario, function(x) readRDS(paste0("mse/om/", x, ".rds")))
names(sra_ye) <- sc$scenario

scenarios <- sc$scenario %>% purrr::set_names(sc$scenario_human)

get_Perr_y <- function(x, scenario) {
  max_age <- x@OM@maxage
  nyears <- x@OM@nyears+x@OM@proyears
  perr_y <- x@OM@cpars$Perr_y[,max_age:(max_age+nyears-1), drop=FALSE]
  all_years <- seq(x@OM@CurrentYr - x@OM@nyears + 1, x@OM@CurrentYr+x@OM@proyears)
  reshape2::melt(perr_y) %>%
    rename(iteration = Var1) %>%
    mutate(year = rep(all_years, each = max(iteration))) %>%
    mutate(scenario = scenario)
}

g <- purrr::map2_df(sra_ye, sc$scenario_human, get_Perr_y) %>%
  mutate(scenario = factor(scenario, levels = sc$scenario_human)) %>%
  dplyr::filter(iteration %in% 1:50) %>%
  ggplot(aes(year, y = log(value), group = iteration)) +
  geom_line(alpha = 0.15) +
  facet_wrap(vars(scenario)) +
  gfplot::theme_pbs() +
  labs(x = "Year", y = "Recruitment deviations in log space") +
  coord_cartesian(ylim = c(-3.5, 3.5), expand = FALSE) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.6)
 #print(g)
.ggsave("ye-projections-recdev-panel-rep_1-50.png",
       width = 8, height = 5)

g <- purrr::map2_df(sra_ye, sc$scenario_human, get_Perr_y) %>%
  mutate(scenario = factor(scenario, levels = sc$scenario_human)) %>%
  dplyr::filter(iteration %in% 1:50) %>%
  ggplot(aes(year, y = value, group = iteration)) +
  geom_line(alpha = 0.15) +
  facet_wrap(vars(scenario)) +
  gfplot::theme_pbs() +
  labs(x = "Year", y = "Recruitment deviations") +
  coord_cartesian(ylim = c(0, exp(3.)), expand = FALSE) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.6)
#print(g)
.ggsave("ye-projections-recdev-panel-exp-rep_1-50.png",
        width = 8, height = 5)
