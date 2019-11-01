# email from Dana:
d2 <- readxl::read_xls("data-raw/complete_research_blocks_table.xls")

library(dplyr)

d2 <- d2 %>%
  filter(exported_research_blocks.SS_ID %in% c(39, 40)) %>%
  select(exported_research_blocks.SS_ALT_DES, ZonalStats_bathy100.MEAN,
    exported_research_blocks.LATITUDE, exported_research_blocks.LONGITUDE, exported_research_blocks.percent100_rock, exported_research_blocks.Percent20m_rock, exported_research_blocks.BLOCK_DESI) %>%
  rename(survey = exported_research_blocks.SS_ALT_DES,
    depth = ZonalStats_bathy100.MEAN, longitude = exported_research_blocks.LONGITUDE,
    latitude = exported_research_blocks.LATITUDE,
    rock20 = exported_research_blocks.Percent20m_rock,
    rock100 = exported_research_blocks.percent100_rock,
    block = exported_research_blocks.BLOCK_DESI) %>%
  mutate(depth = -1 * depth) %>%
  filter(!is.na(depth)) %>% # fixme (removing 8)
  mutate(survey = paste0("HBLL ", survey))

# plot(d2$exported_research_blocks.percent100_rock, d2$exported_research_blocks.Percent20m_rock)

saveRDS(d2, file = "data-generated/hbll-inside-grid.rds")

