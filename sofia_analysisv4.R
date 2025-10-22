#just gonad data, no spawning data
#do i need to convert to cm? 
rm(list=ls())

# Load First --------------------------------------------------------------
#install.packages('librarian')
require(librarian)
librarian::shelf(tidyverse,here, janitor, googlesheets4, lubridate, splitstackshape,
                 googledrive,googlesheets4,httpuv,dplyr,ggplot2,pwr2, tidyr, broom)


# Load Data Sets ---------------------------------------------------------------

gonad_raw <- read_sheet("https://docs.google.com/spreadsheets/d/18R1F5KkILws3e-8CYz83BdGYMG8zeuvOIeKcaSXeqC4/edit?gid=1621016945#gid=1621016945")

#load("/Users/sofiarivas/Downloads/kelp_recovery_data (1).rda") #Sofia
load("/Volumes/enhydra/data/kelp_recovery/MBA_kelp_forest_database/processed/recovery/kelp_recovery_data.rda") #Josh

#load("/Users/sofiarivas/Downloads/lda_patch_transitionsv2.rda") #Sofia
load("/Users/jossmith/code_respositories/kelp_recovery/output/lda_patch_transitionsv2.rda") #Josh


# Reworking Data ----------------------------------------------------------

#clean up patch types and make wider
patch_types <- transitions_tbl_constrained %>%
  rename_with(~ gsub("patch_", "", .x)) %>%
  pivot_longer(cols = c(`2024`, `2025`), names_to = "year", values_to = "new_type") %>%
  mutate(year = as.numeric(year)) %>% 
  mutate(site_id = paste (site,new_type,zone, year))

#filter and clean sizes
urchin_sizefq_1 <- urchin_sizefq %>%
  filter(species=="Purple") %>% 
  mutate(year = year(survey_date)) %>% 
  mutate(site_id = paste(site, site_type,zone, year)) 

#select focal vars from dissection data and join patch types
gonad_build1 <- gonad_raw %>% 
  select("date_collected","survey_type","site_number","site_type","zone",
         "transect","species","sex","test_diameter_mm","gonad_mass_g","animal_24hr_mass_g") %>%
  filter(species %in% c("purple_urchin","purple_urchins")) %>%
  mutate(species = "purple_urchin")%>%
  filter(!is.na(gonad_mass_g)) %>% 
  mutate(year = year(date_collected),
         #year2 = is.character(year),
         site_id = paste(site_number,site_type,zone,year)) %>% 
  filter(gonad_mass_g<30) %>% 
  filter(survey_type == "Recovery") %>%
  #join patch types
  left_join(., patch_types, by = c("site_id") %>% 
  rename(year = year.y) %>% 
  mutate(site_type.x = if_else(year == 2024, "2024", "2025")) %>% 
  mutate(site_id = paste(site_id, year))
 

#Join patch types with urchin size fq
urchin_sizefq_joined <- left_join(urchin_sizefq_1, patch_types, by = "site_id") %>%
  mutate(site_type.x = if_else(year.y == 2024, "2024", "2025")) %>%
  mutate(site_id = paste(site_id, year.y)) %>% 
  rename (year = year.y)

#add year and unique identifier to quad data
quad_working <- quad_data %>% 
  mutate(year = year(survey_date)) %>% 
  mutate(site_id = paste(site, site_type,zone, year))

quad_joined <- left_join(quad_working, patch_types, by = "site_id") %>% 
  mutate(site_type.x = if_else(year.y == 2024, "2024", "2025")) %>% 
  mutate(site_id = paste(site_id, year.y))


# Stats and Calculations-------------------------------------------------------------------

#avg urchin density per site/type/zone 
avg_urchin_density <- quad_joined %>% 
  group_by(site_id) %>%
  # mutate(site_id = paste(site, site_type, zone)) %>% 
  summarize(avg_density = mean(purple_urchin_densitym2, na.rm = TRUE)) %>% 
  ungroup() %>% 
  #unite(col = site_id, site, site_type, zone, sep=" ", remove = FALSE) %>%
  mutate(density80m2 = avg_density*80) #%>% 
#mutate(site_id = toupper(site_id))


#model for urchin size 
coeff_table <- gonad_build1 %>%
  group_by(site_id) %>%
  nest() %>%
  mutate(model = map(data, ~ lm(gonad_mass_g ~ test_diameter_mm, data = .)),
         tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  select(site_id, term, estimate, std.error, statistic, p.value)

coeff_wide <- coeff_table %>%
  select(site_id, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>% 
  rename(b = "(Intercept)", a = "test_diameter_mm")


sampled_urchins <- avg_urchin_density %>%
  left_join(urchin_sizefq_joined, by = "site_id") %>%
  group_by(site_id) %>%
  summarise(
    sampled_sizes = list({
      valid <- !is.na(count) & count > 0
      clean_counts <- count[valid]
      clean_sizes  <- size_cm[valid]
      n_to_sample <- round(first(avg_density))
      if (length(clean_counts) == 0 || sum(clean_counts) == 0 || n_to_sample == 0) {rep(NA, n_to_sample)} 
      else {sample(
        clean_sizes,
        size = n_to_sample,
        replace = TRUE,
        prob = clean_counts / sum(clean_counts))}}),
    .groups = "drop") %>%
  unnest(cols = sampled_sizes)

sampled_urchins_80 <- avg_urchin_density %>%
  left_join(urchin_sizefq_joined, by = "site_id") %>%
  group_by(site_id) %>%
  summarise(
    sampled_sizes = list({
      valid <- !is.na(count) & count > 0
      clean_counts <- count[valid]
      clean_sizes  <- size_cm[valid]
      n_to_sample <- round(first(density80m2))
      if (length(clean_counts) == 0 || sum(clean_counts) == 0 || n_to_sample == 0) {rep(NA, n_to_sample)} 
      else {
        sample(
          clean_sizes,
          size = n_to_sample,
          replace = TRUE,
          prob = clean_counts / sum(clean_counts))}}),
    .groups = "drop") %>%
  unnest(cols = sampled_sizes)

#convert diameter to mass
converted_measurements <- sampled_urchins_80 %>%    
  left_join(coeff_wide, by = "site_id") %>%
  mutate(
    size_mm = sampled_sizes*10,
    biomass_g = -14.2 + 7.44 * exp(0.04 * size_mm)
  ) %>%
  #tidy
  mutate(biomass_g = ifelse(biomass_g < 0,1,biomass_g))

#calculate site level mean gonad mass
urchin_gsi <- gonad_build1 %>%
  mutate(
    GSI = (gonad_mass_g / animal_24hr_mass_g) * 100,
    #GSI_somatic = (gonad_mass_g / (animal_24hr_mass_g - gonad_mass_g)) * 100
  ) %>%
  group_by(site_id) %>%
  summarize(u_GSI = mean(GSI),
            sd_GSI = sd(GSI),
            n_GSI = n(GSI))

#calcualte gonad mass
gonad_mass_site_zone <- converted_measurement %>%
                        left_join(., urchin_gsi, by = "site_id") %>% #warning many-to-many is ok
                        mutate(gonad_mass_g = biomass_g*u_GSI) 

#calculate total gonad mass per site
gonad_mass_site <- goand_mass_site_zone %>%
                  group_by(site_id) %>%
                  summarize(t_gonad_mass = sum(gonad_mass_g))








