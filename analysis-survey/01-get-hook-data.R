library(dplyr)

data <- gfdata::get_ll_hook_data(species = 442, ssid = c(39, 40))

d2 <- gfdata::get_survey_sets("Yelloweye Rockfish", ssid = c(39, 40))

# the following code should deal with the issue of 0 baited hooks being observed.
adjust <- data %>%
  group_by(year, fishing_event_id) %>%
  mutate(total_hooks = count_target_species + count_non_target_species +
    count_bait_only + count_empty_hooks - count_bent_broken) %>%
  mutate(count_bait_only = replace(count_bait_only, which(count_bait_only == 0), 1)) %>%
  mutate(prop_bait_hooks = count_bait_only / total_hooks) %>%
  mutate(hook_adjust_factor = -log(prop_bait_hooks) / (1 - prop_bait_hooks)) %>%
  mutate(expected_catch = round(count_target_species * hook_adjust_factor))

# find out which events are in hook data but not set data:
setdiff(data$fishing_event_id, d2$fishing_event_id)

hook_adjusted_data <- left_join(d2, adjust, by = c("fishing_event_id", "year"))

saveRDS(hook_adjusted_data, file = "data-generated/yelloweye_hook_adjusted_data.rds")
