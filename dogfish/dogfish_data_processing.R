
# "smashing" dogfish data together
# need to add 2004 and 2019 data to the csv file provided by Elise that was extracted using "dogfish survey sets data.R"
# remember that 2004 is the comparison year wherein both circle and J hooks were used.
# 2019 has hook by hoook data.
# written by Midoli Bresch
# last updated 2019-11-21

# clear workspace
rm(list=ls())

library(dplyr)
library(sp)

# read in the data
dat2019 <- read.csv(file = "hbll_dogfish_2019_data.csv", header=TRUE, stringsAsFactors = FALSE)
dat     <- readr::read_csv(file = "dogfish_survey_sets_data.csv")
dat2004 <- readr::read_csv(file = "dogfish_2004.csv")

# extract just the dogfish survey data from dat2019
dat2019 <- dat2019 %>%
  filter(survey == "Dogfish Longline")

# header names that Elise's extraction has
headers <- names(dat)

# set some constants (these are from Elise's email)
hook_spacing <- 0.024384
gangion_length <- 0.009144
year <- 2019
hook_code <- 1 # for circle hooks

# get what headers are missing from 2004 and 2019 data sets

miss2019 <- setdiff(headers, names(dat2019) )
miss2004 <- setdiff(headers, names(dat2004) )

## -------------2019 Data --

ye2019 <- dat2019 %>%
  filter(species_code == 442) %>%
  rename(ye_count = catch_count) %>%
  select(fishing_event_id, ye_count)

dog2019 <- dat2019 %>%
  filter(species_code == 44) %>%
  rename(dogfish_count = catch_count)

# join dogfish_counts to ye2019 data frame, replace NAs with 0s
all2019 <- left_join(dog2019, ye2019, by="fishing_event_id")
all2019 <- replace(all2019, is.na(all2019),0)

# rename some headers, add in year and hook code (1) and area swept km2
data2019 <- all2019 %>%
  rename(latitude = start_latitude,
         longitude = start_longitude,
         depth_m = start_depth_m,
         fe_end_deployment_time = end_deploy_time,
         fe_begin_retrieval_time = start_deploy_time,
         lglsp_hook_count = num_hooks) %>%
  mutate(year = year) %>%
  mutate(hook_code = hook_code) %>%
  mutate(area_swept_km2 = hook_spacing * gangion_length * lglsp_hook_count)

# what is still missing?
miss2019 <- setdiff(headers, names(data2019) )

# what doesn't need to be there?
extra2019 <- setdiff(names(data2019), headers)

# remove extra headings from 2019 data (listed these out so they can easily be changed)
data2019 <- data2019 %>%
  select(-survey, -set_number, -species_code, -species_desc, -survey_abbrev, -fe_major_level_id, -trip_id,
         -block_designation, -end_depth_m, -mean_bottom_depth_m, -end_latitude, -end_longitude,
         -num_empty_hooks, -num_bait_only_hooks)


## ---------------- 2004 Data ------------------------------------------------------

# change lat and long from hours mins secs to decimal degrees (found the conversion in Elise's SQL code)
# add in area_swept_km2
# remove columns with previous location format
data2004 <- dat2004 %>%
    mutate_at("lat_mins", as.numeric) %>% 
    mutate_at("long_mins", as.numeric) %>% 
    mutate(latitude = lat_hours + lat_mins/60,
              longitude = -(long_hours + long_mins/60) ) %>% 
    mutate(area_swept_km2 = hook_spacing * gangion_length * lglsp_hook_count) %>% 
    select(-lat_mins, -lat_hours, -long_mins, -long_hours)
    

# what is missing from 2004?
missing2004 <- setdiff(headers, names(data2004))

## --------------Combine the datasets ---------------------

# remove some headers from dat
dat <- dat %>%
  select(-grouping_code, -grouping_desc, -grouping_depth_id, -fe_bottom_water_temperature)

# remove some headers from data2019
data2019 <- data2019 %>%
  select(-grouping_code)

# remove some headers from data2004
data2004 <- data2004 %>%
  select(-set_number, -location, -depth_stratum)

alldata <- bind_rows(dat, data2019, data2004)

readr::write_csv(alldata, "all_dogfish_sets_data.csv")
