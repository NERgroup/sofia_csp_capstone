rm(list=ls())

# Load First --------------------------------------------------------------
install.packages('librarian')
require(librarian)
librarian::shelf(tidyverse,here, janitor, googlesheets4, lubridate, splitstackshape,
                 googledrive,googlesheets4,httpuv,dplyr,ggplot2,pwr2, tidyr, broom)


# Load Data Sets ---------------------------------------------------------------

gonad_raw <- read_sheet("https://docs.google.com/spreadsheets/d/18R1F5KkILws3e-8CYz83BdGYMG8zeuvOIeKcaSXeqC4/edit?gid=1621016945#gid=1621016945")

load("/Users/sofiarivas/Downloads/kelp_recovery_data (1).rda") #Sofia
#load("/Volumes/enhydra/data/kelp_recovery/MBA_kelp_forest_database/processed/recovery/kelp_recovery_data.rda") #Josh

load("/Users/sofiarivas/Downloads/lda_patch_transitionsv2.rda") #Sofia
#load("/Users/jossmith/code_respositories/kelp_recovery/output/lda_patch_transitionsv2.rda") #Josh


# Reworking Data ----------------------------------------------------------

#includes year in site_id, not used for joining patch types
gonad_working <- gonad_raw %>% 
  select("date_collected","survey_type","site_number","site_type","zone",
         "transect","species","sex","test_diameter_mm","gonad_mass_g") %>%
  filter(species %in% c("purple_urchin","purple_urchins")) %>%
  mutate(species = "purple_urchin")%>%
  filter(!is.na(gonad_mass_g)) %>% 
  mutate(year = year(date_collected),
         #year2 = is.character(year),
         site_id = paste(site_number,site_type,zone,year)) %>% 
  filter(gonad_mass_g<30) %>% 
  filter(survey_type == "Recovery")

#used for joining patch types
gonad_working2 <- gonad_raw %>% 
  select("date_collected","survey_type","site_number","site_type","zone",
         "transect","species","sex","test_diameter_mm","gonad_mass_g","animal_24hr_mass_g") %>%
  filter(species %in% c("purple_urchin","purple_urchins")) %>%
  mutate(species = "purple_urchin")%>%
  filter(!is.na(gonad_mass_g)) %>% 
  mutate(year = year(date_collected),
         #year2 = is.character(year),
         site_id = paste(site_number,site_type,zone)) %>% 
  filter(gonad_mass_g<30) %>% 
  filter(survey_type == "Recovery")

#patch types
patch_types <- transitions_tbl_constrained %>%
  rename_with(~ gsub("patch_", "", .x)) %>%
  pivot_longer(cols = c(`2024`, `2025`), names_to = "year", values_to = "new_type") %>%
  mutate(year = as.numeric(year)) %>% 
  mutate(site_id = paste (site,new_type,zone))

#joined gonad data with patch types 
gonad_patch_joined <- left_join(gonad_working2, patch_types, by = "site_id") %>% 
  rename(year = year.y) %>% 
  mutate(site_type.x = if_else(year == 2024, "2024", "2025")) %>% 
  mutate(site_id = paste(site_id, year)) %>% 
  mutate(test_diameter_cm = test_diameter_mm / 10)

#with rounded diameters 
gonad_patch_joined2 <- left_join(gonad_working2, patch_types, by = "site_id") %>% 
  rename(year = year.y) %>% 
  mutate(site_type.x = if_else(year == 2024, "2024", "2025")) %>% 
  mutate(site_id = paste(site_id, year)) %>% 
  mutate(test_diameter_cm = round(test_diameter_mm / 10))


urchin_sizefq_1 <- urchin_sizefq %>%
  filter(species=="Purple") %>% 
  mutate(year = year(survey_date)) %>% 
  mutate(site_id = paste(site, site_type,zone)) 

#joined urchin size frequency data with patch types
urchin_sizefq_joined <- left_join(urchin_sizefq_1, patch_types, by = "site_id") %>%
  mutate(site_type.x = if_else(year.y == 2024, "2024", "2025")) %>%
  mutate(site_id = paste(site_id, year.y)) %>% 
  rename (year = year.y)

quad_working <- quad_data %>% 
  mutate(year = year(survey_date)) %>% 
  mutate(site_id = paste(site, site_type,zone))


quad_joined <- left_join(quad_working, patch_types, by = "site_id") %>% 
  mutate(site_type.x = if_else(year.y == 2024, "2024", "2025")) %>% 
  mutate(site_id = paste(site_id, year.y))

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

#changed sampled_sizes to test_diameter_cm
sampled_urchins2 <- avg_urchin_density %>%
  left_join(urchin_sizefq_joined, by = "site_id") %>%
  group_by(site_id) %>%
  summarise(
    test_diameter_cm = list({
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
  unnest(cols = test_diameter_cm)


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

sampled_urchins_802 <- avg_urchin_density %>%
  left_join(urchin_sizefq_joined, by = "site_id") %>%
  group_by(site_id) %>%
  summarise(
    test_diameter_cm = list({
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
  unnest(cols = test_diameter_cm)

# Gonad Index -----------------------------------------------------------------

urchin_gsi <- gonad_patch_joined2 %>%
  mutate(
    GSI = (gonad_mass_g / animal_24hr_mass_g) * 100,
    GSI_somatic = (gonad_mass_g / (animal_24hr_mass_g - gonad_mass_g)) * 100
  )
#1st try 
#regression for GSI-diameter relationship
model <- lm(GSI ~ test_diameter_cm, data = urchin_gsi)
sampled_with_gsi <- sampled_urchins2 %>%
  mutate(GSI_pred = predict(model, newdata = sampled_urchins2))
#w/ random biological variation?
sampled_with_gsi2 <- sampled_with_gsi %>%
  mutate(GSI_sampled = GSI_pred + rnorm(n(), mean = 0, sd = sd(urchin_gsi$GSI, na.rm = TRUE)))

#2nd try 
#closest test diameter within the same site (diameters were previously rounded)
sampled_with_gsi3 <- sampled_urchins_802 %>%
  rowwise() %>%
  mutate(
    matched_GSI = {site_data <- urchin_gsi %>% filter(site_id == site_id)
      idx <- which.min(abs(site_data$test_diameter_cm - test_diameter_cm))
      if (length(idx) == 0) NA_real_ else site_data$GSI[idx]}) %>%
  ungroup()

