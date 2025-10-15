#just gonad data, no spawning data
rm(list=ls())

# Load First --------------------------------------------------------------
install.packages('librarian')
require(librarian)
librarian::shelf(tidyverse,here, janitor, googlesheets4, lubridate, splitstackshape,
                 googledrive,googlesheets4,httpuv,dplyr,ggplot2,pwr2, tidyr, broom)


# Load Data Sets ---------------------------------------------------------------

gonad_raw <- read_sheet("https://docs.google.com/spreadsheets/d/18R1F5KkILws3e-8CYz83BdGYMG8zeuvOIeKcaSXeqC4/edit?gid=1621016945#gid=1621016945")

load("/Users/sofiarivas/Downloads/kelp_recovery_data (1).rda")

load("/Users/sofiarivas/Downloads/lda_patch_transitionsv2.rda")

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
         "transect","species","sex","test_diameter_mm","gonad_mass_g") %>%
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
  mutate(site_id = paste (site,site_type,zone)) %>% 
  rename_with(~ gsub("patch_", "", .x))
  pivot_longer(cols = c(`2024`, `2025`), names_to = "year", values_to = "new_zone") %>%
  mutate(year = as.numeric(year))
  
#joined gonad data with patch types 
gonad_patch_joined <- left_join(gonad_working2, patch_types, by = "site_id") %>% 
  mutate(site_type.x = if_else(year == 2024, `2024`, `2025`)) %>% 
  mutate(site_id = paste(site_id, year))
  
urchin_sizefq_1 <- urchin_sizefq %>%
  filter(species=="Purple") %>% 
  mutate(year = year(survey_date)) %>% 
  mutate(site_id = paste(site, site_type,zone)) 

#joined urchin size frequency data with patch types
urchin_sizefq_joined <- left_join(urchin_sizefq_1, patch_types, by = "site_id") %>% 
  mutate(site_type.x = if_else(year == 2024, `2024`, `2025`)) %>% 
  mutate(site_id = paste(site_id, year))

quad_working <- quad_data %>% 
  mutate(year = year(survey_date)) %>% 
  mutate(site_id = paste(site, site_type,zone))

quad_joined <- left_join(quad_working, patch_types, by = "site_id") %>% 
  mutate(site_type.x = if_else(year == 2024, `2024`, `2025`)) %>% 
  mutate(site_id = paste(site_id, year))

# Figures -----------------------------------------------------------------

#urchin sizefq for site/type/zone 
ggplot(urchin_sizefq_joined,
       aes(x=size_cm, y=count))+
  facet_wrap(~site_id)+
  geom_col(fill="skyblue3")+
  theme_classic()

#urchin diameter x gonad mass for site/type/zone (y=mx+b)
ggplot(gonad_patch_joined,
       aes(x=test_diameter_mm, y=gonad_mass_g))+
  facet_wrap(~site_id)+
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "aquamarine4")+
  theme_classic()

#sampling 
ggplot(sampled_urchins_80, aes(x = sampled_sizes)) +
  geom_histogram(binwidth = 1, color = "white") +
  facet_wrap(~ site_id, scales = "free_y") +
  labs(
    title = "Simulated Urchin Size Distributions by Site",
    x = "Urchin size (cm)",
    y = "Count of sampled individuals"
  ) +
  theme_classic(base_size = 13)

ggplot(sampled_urchins, aes(x = sampled_sizes, fill = site_id)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ site_id, scales = "free_y") +
  labs(x = "Urchin size (cm)", y = "Density") +
  theme_classic(base_size = 13)


# Stats -------------------------------------------------------------------

#avg urchin density per site/type/zone 
avg_urchin_density <- quad_data %>% 
  group_by(site, site_type, zone) %>%
  # mutate(site_id = paste(site, site_type, zone)) %>% 
  summarize(avg_density = mean(purple_urchin_densitym2, na.rm = TRUE)) %>% 
  ungroup() %>% 
  unite(col = site_id, site, site_type, zone, sep=" ", remove = FALSE) %>%
  mutate(density80m2 = avg_density*80) #%>% 
#mutate(site_id = toupper(site_id))


#model for urchin size 
coeff_table <- gonad_working %>%
  group_by(site_id) %>%
  nest() %>%
  mutate(model = map(data, ~ lm(Gonad_Mass_g ~ Test_Diameter_mm, data = .)),
         tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  select(site_id, term, estimate, std.error, statistic, p.value)

coeff_wide <- coeff_table %>%
  select(site_id, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>% 
  rename(b = "(Intercept)", a = "Test_Diameter_mm")

#sample from size distribution
library(purrr)

sampled_urchins <- avg_urchin_density %>%
  left_join(urchin_sizefq_1, by = "site_id") %>%
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
  left_join(urchin_sizefq_1, by = "site_id") %>%
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

coeff_wide2 <- coeff_wide %>%
  mutate(site_id = site_id %>%
           tolower() %>%                 # make lowercase
           str_replace_all(" ", "_") %>% # replace spaces with _
           str_replace_all("-", "_")     # replace hyphens with _
  )

sampled_urchins_802 <- sampled_urchins_80 %>%
  mutate(site_id = site_id %>%
           # Remove any underscore immediately after letters at start (e.g., REC_)
           str_replace("^([A-Za-z]+)_", "\\1") %>%
           # Make everything lowercase
           tolower() %>%
           # Replace spaces and hyphens with underscores
           str_replace_all("[ -]", "_")
  )

#convert diameter to mass
converted_measurements <- sampled_urchins_802 %>%    
  left_join(coeff_wide2, by = "site_id") %>%
  mutate(
    predicted_mass = a * sampled_sizes + b) #%>%
#drop sites without dissected urchins
#filter(!(is.na(predicted_mass)))

setdiff(sampled_urchins_802$site_id, coeff_wide2$site_id)

