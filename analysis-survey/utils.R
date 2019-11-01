# sf::st_crs(3156)
convert2utm <- function(df, coords = c("X", "Y"), out_crs = 3156) {
  x <- sf::st_as_sf(df,
    coords = coords, crs = 4326
  ) %>%
    sf::st_transform(crs = out_crs) %>%
    sf::st_coordinates() %>%
    as.data.frame()
  x$X <- x$X / 100000
  x$Y <- x$Y / 100000
  dplyr::bind_cols(x,
    df[, which(!names(df) %in% coords), drop = FALSE])
}

expand_prediction_grid <- function(grid, years) {
  nd <- do.call("rbind",
    replicate(length(years), grid, simplify = FALSE))
  nd[["year"]] <- rep(years, each = nrow(grid))
  nd
}

boot_biomass <- function(dat, reps = 100) {
  out <- dat %>%
    group_by(year, survey) %>%
    do({
      b <- boot::boot(., statistic = calc_bio, strata = .$grouping_code, R = reps)
      suppressWarnings(bci <- boot::boot.ci(b, type = "perc"))
      tibble::tibble(
        mean_boot = mean(b$t),
        median_boot = median(b$t),
        lwr = bci$percent[[4]],
        upr = bci$percent[[5]],
        cv = sd(b$t)/mean(b$t),
        biomass = calc_bio(.))
    })
}

calc_bio <- function(dat, i = seq_len(nrow(dat))) {
  dat[i, ] %>% group_by(year, survey, grouping_code) %>%
    summarise(density = mean(density_ppkm2)) %>%
    group_by(year) %>%
    summarise(biomass = sum(density * 2 * 2)) %>%
    pull(biomass)
}

