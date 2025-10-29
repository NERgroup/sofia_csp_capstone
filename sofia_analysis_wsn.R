
# load first --------------------------------------------------------------
rm(list=ls())

require(librarian)
librarian::shelf(tidyverse,here, janitor, googlesheets4, lubridate, splitstackshape,
                 googledrive,googlesheets4,httpuv,dplyr,ggplot2,pwr2, tidyr, broom,
                 ggpubr,scales)

gonad_raw <- read_sheet("https://docs.google.com/spreadsheets/d/18R1F5KkILws3e-8CYz83BdGYMG8zeuvOIeKcaSXeqC4/edit?gid=1621016945#gid=1621016945")

spawn_raw <- read_sheet("https://docs.google.com/spreadsheets/d/1_uwJqlLdNA3VaHjxfXEMuBcUgZVhKT99hzrV2IXosT4/edit?usp=sharing") 

load("/Users/sofiarivas/Downloads/kelp_recovery_data (1).rda") 

load("/Users/sofiarivas/Downloads/lda_patch_transitionsv2.rda")

# reworking data frames -------------------------------------------------------------
patch_types <- transitions_tbl_constrained %>%
  rename_with(~ gsub("patch_", "", .x)) %>%
  pivot_longer(cols = c(`2024`, `2025`), names_to = "year", values_to = "new_type") %>%
  mutate(year = as.numeric(year),
         patch_type = site_type) %>% 
  mutate(site_id = paste (site,new_type,zone, year))

spawn_joined <- spawn_raw %>% 
  select(
    Date_Collected, Site_Number, State, Transect, Sex, Test_Height_mm, Test_Diameter_mm,
    Animal_24hr_Mass_g, Gonad_Mass_g, Spawn_Mass_g, Spawn_Mass_false) %>% 
  rename_with(tolower) %>% 
  mutate(
    year = year(date_collected),
    patch_type = state,
    zone = str_to_title(str_to_lower(transect)),
    spawn_mass_g = as.numeric(unlist(spawn_mass_g)),
    spawn_mass_false = as.numeric(unlist(spawn_mass_false)),
    site = site_number %>%
      str_remove("-.*") %>%
      str_replace("([A-Za-z]+)([0-9]+)", "\\1_\\2"),
    site_id = paste(site, patch_type, zone, year)) %>% 
  left_join(
    patch_types %>% select(site_id, new_type),
    by = "site_id") %>% 
  mutate(
    site_type = if_else(year == 2025 & !is.na(new_type), new_type, patch_type)) %>% 
  select(-new_type)

counts <- spawn_joined %>%
  filter(spawn_mass_g < 4) %>%
  group_by(patch_type) %>%
  summarise(n = n())

spawn_sex <- spawn_joined %>%
  mutate(sex = as.character(sex)) %>%
  mutate(sex = ifelse(is.na(sex) | sex == "" | sex == "NA", "Female", sex),
         sex = str_to_title(str_trim(sex))) %>%
  group_by(patch_type, sex) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(patch_type) %>%
  mutate(prop = n / sum(n),
         total_n = sum(n)) %>%
  ungroup()

# plots -------------------------------------------------------------------

#spawn mass x patch type (raw spawn mass)
ggplot(spawn_joined %>% filter(spawn_mass_g < 4),
       aes(x = patch_type, y = spawn_mass_g, fill = patch_type)) +
  geom_point(aes(color = patch_type),
             position = position_jitter(width = 0.1),
             alpha = 0.3, size = 2) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(
    name = "Patch Type",
    values = c("BAR" = "slateblue", "FOR" = "darkolivegreen", "INCIP" = "coral")) +
  scale_color_manual(
    name = "Patch Type",
    values = c("BAR" = "slateblue", "FOR" = "darkolivegreen", "INCIP" = "coral")) +
  scale_x_discrete(
    labels = c(
      "BAR" = paste0("Barren\n(n=", spawn_joined %>% filter(patch_type == "BAR", spawn_mass_g < 4) %>% nrow(), ")"),
      "FOR" = paste0("Forest\n(n=", spawn_joined %>% filter(patch_type == "FOR", spawn_mass_g < 4) %>% nrow(), ")"),
      "INCIP" = paste0("Incipient\n(n=", spawn_joined %>% filter(patch_type == "INCIP", spawn_mass_g < 4) %>% nrow(), ")"))) +
  labs(
    x = "Patch Type",
    y = "Spawned Mass (g)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10))+
  stat_compare_means(
    comparisons = list(
      c("BAR", "FOR"),
      c("BAR", "INCIP"),
      c("FOR", "INCIP")),
    method = "t.test",       
    label = "p.format",     
    size = 3.5,
    y.position = c(10,10,10))
ggsave("patchtype_spawnmass.png")

#sex ratio from spawn_joined
ggplot(spawn_sex, aes(x = patch_type, y = prop, fill = sex)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  geom_text(aes(label = paste0(round(prop*100), "%")), 
            position = position_fill(vjust = 0.5), size = 3) +  # percentages in middle of segments
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(
    name = "Sex",
    values = c("Male" = "steelblue", "Female" = "indianred")) +
  scale_x_discrete(labels = c("BAR" = "Barren", "FOR" = "Forest", "INCIP" = "Incipient")) +
  labs(
    x = "Patch Type",
    y = "Proportion of Individuals") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10))
ggsave("patchtype_sexratio.png")


ggplot(spawn_sex, aes(x = patch_type, y = n, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("Male" = "skyblue", "Female" = "pink"))


