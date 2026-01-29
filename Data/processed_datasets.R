# Load First -------------------------------------------------------------------

rm(list=ls())

require(librarian)
librarian::shelf(tidyverse,here, janitor, googlesheets4, lubridate, splitstackshape,
                 googledrive,httpuv,dplyr,ggplot2,pwr2,tidyr,broom,ggpubr, paletteer)

googlesheets4::gs4_auth(cache = TRUE)

# Data Sets --------------------------------------------------------------------

#dissection data
gonad_raw <- read_sheet("https://docs.google.com/spreadsheets/d/18R1F5KkILws3e-8CYz83BdGYMG8zeuvOIeKcaSXeqC4/edit?gid=1621016945#gid=1621016945")

#kelp_data, quad_data, urchin_sizefq
load("Data/kelp_recovery_data.rda")

#correct patch types 2024-2025
load("Data/lda_patch_transitionsv2.rda")

#benthic survey data averaged to zone level
load(file.path("/Volumes/enhydra/data/students/sofia/zone_level_data.rda"))

# Data Wrangling ----------------------------------------------------------

#update patch types
patch_types <- transitions_tbl_constrained %>%
  rename_with(~ gsub("patch_", "", .x)) %>% 
  pivot_longer(cols = c(`2024`, `2025`), names_to = "year", values_to = "new_type") %>%
  mutate(year = as.numeric(year)) %>% 
  mutate(site_id = paste (site,site_type,zone,year), #used to join 
         site_id_final = paste (site,new_type,zone,year)) #updated patch types

#filter and clean up
urchin_sizefq1 <- urchin_sizefq %>%
  filter(species=="Purple") %>% 
  mutate(year = year(survey_date),
         patch_type = site_type) %>% 
  select(-c(survey_date, site_type, survey_type, region, latitude, longitude, depth_m)) %>% 
  mutate(site_id = paste(site, patch_type,zone,year)) 

#join patch_types with urchin_sizefq1
urchin_sizefq_joined <- left_join(urchin_sizefq1, patch_types, by = "site_id") %>%
  mutate(site_type.x = if_else(year.y == 2024, "2024", "2025")) %>%
  mutate(site_id = paste(site_id, year.y)) %>% 
  rename(year = year.y) %>% 
  mutate(site_id = word(site_id,1,4)) %>% 
  select(-site_id) 

urchin_sizefq_joined <- full_join(patch_types, urchin_sizefq1, by = "site_id") %>% 
  mutate(site_id_final = coalesce(site_id_final, site_id)) %>% #If site_id_final is NA, replace it with the value in site_id on the same row (for REC 13 and 14 since they were only done in 2025)
  select(-site_id) #site_id_final has correct patch types 

#saved dataframe to project directory --> Data folder
saveRDS(urchin_sizefq_joined, file = "Data/urchin_sizefq_joined.rds")

#select focal variables from dissection data 
gonad_working <- gonad_raw %>% 
  select(
    date_collected, survey_type, site_number, site_type, zone, species, sex, 
    test_diameter_mm, gonad_mass_g, animal_24hr_mass_g) %>%
  filter(species %in% c("purple_urchin", "purple_urchins")) %>%
  mutate(species = "purple_urchin") %>%
  filter(!is.na(gonad_mass_g)) %>% 
  mutate(
    year = as.numeric(year(date_collected)),
    site_id = paste(site_number, site_type, zone, year),
    site_id = str_squish(site_id)) %>% 
  filter(gonad_mass_g < 30, survey_type == "Recovery") 

#join patch_types to gonad_working
gonad_joined <- left_join(gonad_working, patch_types, by = "site_id") %>% 
  mutate(site_id_final = coalesce(site_id_final, site_id)) %>% #If site_id_final is NA, replace it with the value in site_id on the same row (for REC 13 and 14 since they were only done for dissections in 2024)
  select(-site_id) #site_id_final has correct patch types 

#saved dataframe to project directory --> Data folder
saveRDS(gonad_joined, file = "Data/gonad_joined.rds")

#add site_id to quad data
quad_working <- quad_data %>% 
  mutate(year = year(survey_date),
         site_id = paste(site, site_type, zone, year))

quad_joined <- left_join(quad_working, patch_types, by = "site_id") %>% 
  mutate(site_id_final = coalesce(site_id_final, site_id)) %>% #If site_id_final is NA, replace it with the value in site_id on the same row (for REC 13 and 14 since they were only done in 2025)
  select(-site_id) #site_id_final has correct patch types 

#saved dataframe to project directory --> Data folder
saveRDS(quad_joined, file = "Data/quad_joined.rds")

#saved dataframe to project directory --> Data folder ** has not been wrangled 
saveRDS(quad_build3, file = "Data/quad_build3.rds")




