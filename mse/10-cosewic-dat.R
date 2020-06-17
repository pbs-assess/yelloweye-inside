library(purrr)
library(dplyr)
library(ggplot2)

starting_year <- 1918
ending_year <- 2019
all_years <- seq(starting_year, ending_year)
nyear <- length(all_years)
sc <- readRDS("mse/om/ye-scenarios2.rds")
sra_ye <- map(sc$scenario, function(x) readRDS(paste0("mse/om/", x, ".rds")))
names(sra_ye) <- sc$scenario
scenarios <- sc$scenario %>% set_names(sc$scenario_human)
names(sra_ye) <- sc$scenario_human

get_ts <- function(x, type = c("B", "SSB", "depletion")) {
  type <- match.arg(type)
  if (type == "B") dat <- t(sapply(x@Misc, getElement, "B"))
  if (type == "depletion") dat <- x@SSB / sapply(x@Misc, getElement, "E0_SR")
  if (type == "SSB") dat <- x@SSB
  reshape2::melt(dat[, -nyear]) %>%
    rename(iteration = Var1) %>%
    select(-Var2) %>%
    mutate(year = rep(all_years, each = max(iteration)))
}

get_apical_f <- function(x) {
  .F1 <- map(x@Misc, "F_at_age")
  .F <- map_dfc(.F1, ~ tibble(.F = apply(.x, 1, max)))
  .F <- t(as.matrix(.F))
  row.names(.F) <- NULL
  last_year <- dim(.F)[2]
  all_years <- seq(x@OM@CurrentYr - x@OM@nyears + 1, x@OM@CurrentYr)
  reshape2::melt(.F) %>%
    rename(iteration = Var1) %>%
    select(-Var2) %>%
    mutate(year = rep(all_years, each = max(iteration)))
}

ssb <- map_dfr(sra_ye, get_ts, type = "SSB", .id = "scenario") %>%
  rename(ssb = value)
biomass <- map_dfr(sra_ye, get_ts, type = "B", .id = "scenario") %>%
  rename(biomass = value)
apical_f <- map_dfr(sra_ye, get_apical_f, .id = "scenario") %>%
  rename(apical_f = value)

# all iterations of the 3 values:
out <- left_join(ssb, biomass) %>%
  left_join(apical_f) %>%
  as_tibble()

# get median and 95% intervals for ssb, b, f
p <- c(0.025, 0.5, 0.975)
p_names <- c("lwr", "med", "upr")

p_funs <- map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>%
  set_names(nm = p_names)

summaries <- out %>%
  group_by(scenario, year) %>%
  summarize_at(vars(biomass, ssb, apical_f), p_funs)
