# Load First -------------------------------------------------------------------

rm(list=ls())

require(librarian)
librarian::shelf(tidyverse,here, janitor, googlesheets4, lubridate, splitstackshape,
                 googledrive,httpuv,dplyr,ggplot2,pwr2,tidyr,broom,ggpubr, paletteer)

# Data Sets --------------------------------------------------------------------

quad_joined <- readRDS("Data/quad_joined.rds")
gonad_joined <- readRDS("Data/gonad_joined.rds")
urchin_sizefq_joined <- readRDS("Data/urchin_sizefq_joined.rds")

# Stats and Calculations -------------------------------------------------------

#avg urchin density per m2 and 80m2 for each site_id 
avg_urchin_density <- quad_joined %>% 
  group_by(site_id_final) %>%
  summarize(avg_density = mean(purple_urchin_densitym2, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(density80m2 = avg_density*80) 

#for each site, estimate how many urchins are present in an 80 m² area, then simulate individual urchin sizes by sampling from observed size-frequency distributions, weighted by observed counts
sampled_urchins_80 <- left_join(avg_urchin_density, 
                                urchin_sizefq_joined,
                                by = "site_id_final") %>% 
  select(-c(site.x, zone.x, site_type, year.x, site.y, zone.y, year.y, patch_type)) %>%  
  group_by(site_id_final) %>% 
  summarize(sampled_sizes = list ({
                                    valid <- !is.na(count) & count > 0
                                    clean_counts <- count[valid]
                                    clean_sizes  <- size_cm[valid]
                                    n_to_sample <- round(first(density80m2))
                                 if (length(clean_counts) == 0 || sum(clean_counts) == 0 || n_to_sample == 0) 
                                    {rep(NA, n_to_sample)} 
                               else {sample(clean_sizes,
                                             size = n_to_sample,
                                             replace = TRUE,
                                             prob = clean_counts / sum(clean_counts))}}),
                                   .groups = "drop") %>%
  unnest(cols = sampled_sizes)







#model for urchin size 
coeff_table <- gonad_joined %>%
  group_by(site_id_final) %>%
  nest() %>%
  mutate(model = map(data, ~ lm(gonad_mass_g ~ test_diameter_mm, data = .)),
         tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  select(site_id_final, term, estimate, std.error, statistic, p.value)

coeff_wide <- coeff_table %>%
  select(site_id_final, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>% 
  rename(b = "(Intercept)", a = "test_diameter_mm")

# sampled_urchins <- avg_urchin_density %>%
#   left_join(urchin_sizefq_joined, by = "site_id") %>%
#   group_by(site_id) %>%
#   summarise(
#     sampled_sizes = list({
#       valid <- !is.na(count) & count > 0
#       clean_counts <- count[valid]
#       clean_sizes  <- size_cm[valid]
#       n_to_sample <- round(first(avg_density))
#       if (length(clean_counts) == 0 || sum(clean_counts) == 0 || n_to_sample == 0) {rep(NA, n_to_sample)} 
#       else {sample(
#         clean_sizes,
#         size = n_to_sample,
#         replace = TRUE,
#         prob = clean_counts / sum(clean_counts))}}),
#     .groups = "drop") %>%
#   unnest(cols = sampled_sizes)



#convert diameter to mass
converted_measurements <- sampled_urchins_80 %>%    
  left_join(coeff_wide, by = "site_id_final") %>%
  mutate(
    size_mm = sampled_sizes*10,
    biomass_g = -14.2 + 7.44 * exp(0.04 * size_mm)) %>%
  mutate(biomass_g = ifelse(biomass_g < 0,1,biomass_g))

#calculate site level mean gonad mass
urchin_gsi <- gonad_joined %>%
  mutate(GSI = gonad_mass_g / animal_24hr_mass_g) %>%
  group_by(site_id_final) %>%
  summarize(
    u_GSI = mean(GSI, na.rm = TRUE),    # mean GSI per site
    sd_GSI = sd(GSI, na.rm = TRUE),     # SD per site
    n_GSI = n())                        # number of urchins per site

set.seed(12)
#calcualte gonad mass
gonad_mass_site_zone <- converted_measurements %>%
  left_join(., urchin_gsi, by = "site_id_final") %>% #warning many-to-many is ok
  mutate(gonad_mass_g = biomass_g*u_GSI) %>% 
  group_by(site_id_final) %>% 
  mutate(GSI_sim = rnorm(n(), mean = u_GSI, sd = sd_GSI))%>% #draw from normal dist 
  ungroup() %>% 
  mutate(GSI_sim = ifelse(GSI_sim<0,0,GSI_sim),
         gonad_mass_sim = biomass_g*GSI_sim) %>% 
  drop_na()

#calculate total gonad mass per site
gonad_mass_site_total <- gonad_mass_site_zone %>%
  group_by(site_id_final) %>%
  summarize(t_gonad_mass = sum(gonad_mass_g),
            t_gonad_mass_sim = sum(gonad_mass_sim),
            n_urch = n(),
            t_biomass_g = sum(biomass_g),
            t_biomass_kg = t_biomass_g/1000) %>%
  mutate(site_type = word(site_id_final, 2),
         t_gonad_mass_kg = t_gonad_mass/1000,
         t_gonad_mass_sim_kg = t_gonad_mass_sim/1000,
         t_gonad_mass_sim_g_m2 = t_gonad_mass_sim/80,
         t_biomass_sim_g_m2 = t_biomass_g/80) %>% 
  filter(site_type %in% c("BAR", "FOR","INCIP"))
#filter(t_gonad_mass_kg<400)

#total urchins x total gonad mass (natural variability)
plot(t_gonad_mass_sim_kg~n_urch, data = gonad_mass_site_total)#per 80m2 

#total urchins x total gonad mass (natural variability)
ggplot(gonad_mass_site_total, aes(x = n_urch, y = t_gonad_mass_sim_kg, fill = site_type))+
  geom_point()+
  geom_smooth(method = "loess", span = 1)



#total biomass x total gonad mass (with incipient)
ggplot(
  gonad_mass_site_total %>% filter(t_biomass_kg < 60),
  aes(x = t_biomass_kg, y = t_gonad_mass_sim_kg, color = site_type, fill = site_type)) +
  geom_point(alpha = 0.3, size = 2) +      
  geom_smooth(aes(color = site_type), method = "loess", span = 1, se = TRUE, size = 0.8) +  
  scale_fill_manual(
    name = "Patch Type",
    values = c(
      "BAR" = "slateblue",
      "FOR" = "darkolivegreen",
      "INCIP" = "coral"),
    labels = c("BAR" = "Barren", "FOR" = "Forest", "INCIP" = "Incipient")) +
  scale_color_manual(
    name = "Patch Type",
    values = c(
      "BAR" = "slateblue",
      "FOR" = "darkolivegreen",
      "INCIP" = "coral"),
    labels = c("BAR" = "Barren", "FOR" = "Forest", "INCIP" = "Incipient")) +
  labs(
    x = "Total Biomass (kg)",
    y = "Total Gonad Mass (kg)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10))
ggsave("biomass_gonadmass_incipient.png")

#total biomass x total gonad mass (without incipient)
ggplot(
  gonad_mass_site_total %>% filter(site_type != "INCIP", t_gonad_mass_sim_g_m2 <58),
  aes(x = t_biomass_sim_g_m2, y = t_gonad_mass_sim_g_m2, color = site_type, fill = site_type)
) +
  geom_point(alpha = 0.3, size = 2) +
  geom_smooth(aes(color = site_type), method = "loess", span = 1, se = TRUE, size = 1, alpha = 0.4) +
  scale_fill_manual(
    name = "Patch Type",
    values = c("BAR" = "#7570B3FF", "FOR" = "#1B9E77FF"),
    labels = c("BAR" = "Barren", "FOR" = "Forest")
  ) +
  scale_color_manual(
    name = "Patch Type",
    values = c("BAR" = "#7570B3FF", "FOR" = "#1B9E77FF"),
    labels = c("BAR" = "Barren", "FOR" = "Forest")
  ) +
  labs(
    x = "Total Urchin Biomass (g/m²)",
    y = "Total Gonad Mass (g/m²)"
  ) +
  #coord_cartesian(xlim = c(0, 624)) +  
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
ggsave("biomass_gonadmass.png")

#site type x gonad mass per 80m2 boxplot (natural variability)
ggplot(gonad_mass_site_total, aes(x = site_type, y = t_gonad_mass_sim_kg, fill = site_type)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.7) +
  scale_fill_manual(values = c(
    "BAR" = "mediumpurple3",
    "FOR" = "seagreen",
    "INCIP" = "lightblue")) +
  labs(x = "Site Type", y = "Total Gonad Mass per 80m² (kg)") +
  theme_classic() +
  stat_compare_means(
    comparisons = list(
      c("BAR", "FOR"),
      c("BAR", "INCIP"),
      c("FOR", "INCIP")),
    method = "t.test",
    label = "p.signif")

#site type x gonad mass per m2 boxplot (natural variability)
ggplot(gonad_mass_site_total %>% filter(site_type != "INCIP"), aes(x = site_type, y = t_gonad_mass_sim/80, fill = site_type)) +
  geom_point(aes(color = site_type), position = position_jitter(width = 0.1), alpha = 0.3, size = 2)+
  geom_boxplot(width = 0.6) +
  scale_fill_manual(
    name = "Patch Type",   
    values = c(
      "BAR" = "#7570B3FF",
      "FOR" = "#1B9E77FF"),
    labels = c("BAR" = "Barren", "FOR" = "Forest", "INCIP" = "Incipient")) +
  scale_color_manual(
    name = "Patch Type",   
    values = c(
      "BAR" = "#7570B3FF",
      "FOR" = "#1B9E77FF"),
    labels = c("BAR" = "Barren", "FOR" = "Forest")) +
  scale_x_discrete(labels = c("BAR" = "Barren", "FOR" = "Forest")) +
  labs(
    x = "Patch Type",
    y = "Total Gonad Mass (g/m²)") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))+
  stat_compare_means(
    comparisons = list(
      c("BAR", "FOR")),
    method = "t.test",
    label = "p.format",
    size = 3.5)
ggsave("gonadmass_m2_patchtype_boxplot.png", width = 6, height = 4)


#site type x gonad mass per m2 boxplot (not natural variability)
ggplot(gonad_mass_site_total, aes(x = site_type, y = t_gonad_mass/80, fill = site_type)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.7) +
  scale_fill_manual(values = c(
    "BAR" = "mediumpurple3",
    "FOR" = "seagreen",
    "INCIP" = "lightblue")) +
  labs(x = "Site Type", y = "Total Gonad Mass per m² (g)") +
  theme_classic() +
  stat_compare_means(
    comparisons = list(
      c("BAR", "FOR"),
      c("BAR", "INCIP"),
      c("FOR", "INCIP")),
    method = "t.test",
    label = "p.signif")


#stats
anova_model <- aov(t_gonad_mass ~ site_type, data = gonad_mass_site_total)
summary(anova_model)

TukeyHSD(anova_model)

gonad_mass_site_total %>%
  group_by(site_type) %>%
  summarise(shapiro_p = shapiro.test(t_gonad_mass)$p.value)

kruskal.test(t_gonad_mass ~ site_type, data = gonad_mass_site_total)





