
# load first --------------------------------------------------------------
rm(list=ls())

require(librarian)
librarian::shelf(tidyverse,here, janitor, googlesheets4, lubridate, splitstackshape,
                 googledrive,googlesheets4,httpuv,dplyr,ggplot2,pwr2, tidyr, broom,
                 ggpubr)

gonad_raw <- read_sheet("https://docs.google.com/spreadsheets/d/18R1F5KkILws3e-8CYz83BdGYMG8zeuvOIeKcaSXeqC4/edit?gid=1621016945#gid=1621016945")

spawn_raw <- read_sheet("https://docs.google.com/spreadsheets/d/1_uwJqlLdNA3VaHjxfXEMuBcUgZVhKT99hzrV2IXosT4/edit?usp=sharing") 

load("/Users/sofiarivas/Downloads/kelp_recovery_data (1).rda") 

load("/Users/sofiarivas/Downloads/lda_patch_transitionsv2.rda")

# reworking data frames -------------------------------------------------------------
patch_types <- transitions_tbl_constrained %>%
  rename_with(~ gsub("patch_", "", .x)) %>%
  pivot_longer(cols = c(`2024`, `2025`), names_to = "year", values_to = "new_type") %>%
  mutate(year = as.numeric(year)) %>% 
  mutate(site_id = paste (site,new_type,zone, year))

spawn_joined <- spawn_raw %>% 
  select(Date_Collected, Site_Number, State, Transect, Sex, Test_Height_mm, Test_Diameter_mm,
         Animal_24hr_Mass_g, Gonad_Mass_g, Spawn_Mass_g, Spawn_Mass_false) %>% 
  rename_with(tolower) %>% 
  mutate(
    year = year(date_collected),
    site_type = state,
    zone = str_to_title(str_to_lower(transect)),
    site = site_number %>%
      str_remove("-.*") %>%                      
      str_replace("([A-Za-z]+)([0-9]+)", "\\1_\\2"),  
    site_id = paste(site, site_type, zone, year)) %>% 
  left_join(patch_types %>% 
          select(site_id,new_type),
  by = "site_id") %>% 
  mutate(
    site_type = if_else(year == 2025 & !is.na(new_type), new_type, site_type)) %>%
  select(-new_type)
    