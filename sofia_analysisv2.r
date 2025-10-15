rm(list=ls())

# Load First --------------------------------------------------------------
install.packages('librarian')
require(librarian)
librarian::shelf(tidyverse,here, janitor, googlesheets4, lubridate, splitstackshape,
                 googledrive,googlesheets4,httpuv,dplyr,ggplot2,pwr2, tidyr, broom)


# Data Sets ---------------------------------------------------------------

spawn_raw <- read_sheet("https://docs.google.com/spreadsheets/d/11DLr38iVRDcvWiDoBY1hl_Cn5M9oAxDwBRQsx_eMbHk/edit?usp=sharing") 

gonad_raw <- read_sheet("https://docs.google.com/spreadsheets/d/18R1F5KkILws3e-8CYz83BdGYMG8zeuvOIeKcaSXeqC4/edit?gid=1621016945#gid=1621016945")


load("/Users/sofiarivas/Downloads/kelp_recovery_data (1).rda")

load("/Users/sofiarivas/Downloads/lda_patch_transitionsv2.rda")


#delete unnecessary rows
spawn_working <- subset(spawn_raw, select = -c(Data_Enterer,Gonad_Mass_total,Spawn_Mass_total,Treatment))

gonad_working <- subset(gonad_raw, select = -c(Institution,Name_of_Data_Enterer,Date_Fixed,Treatment,Soft_Tissue_Mass_g,Notes,Date_Entered)) %>%
  filter(Species %in% c("PUR", "purple_urchin")) %>%
  mutate(Species = "purple_urchin")%>%
  filter(!is.na(Gonad_Mass_g)) %>% 
  mutate(site_id = paste(Site_Number, Transect)) %>% 
  filter(year(Date_Collected) == 2025,
         Gonad_Mass_g<30) %>% 
  filter(!str_starts(Site_Number, "MAR")) %>% 
  mutate(site_id = gsub("DEEP$", "Deep", site_id)) %>%
  mutate(site_id = gsub("SHALLOW$", "Shallow", site_id))


urchin_sizefq_1 <- urchin_sizefq %>%
  filter(species=="Purple") %>% 
  mutate(site_id = paste(site, site_type,zone)) %>% 
  filter(year(survey_date) == 2025) 

#make values numeric
spawn_working$Spawn_Mass_g <- as.numeric(as.character(spawn_working$Spawn_Mass_g))
spawn_working$Spawn_Mass_false <- as.numeric(as.character(spawn_working$Spawn_Mass_false))
spawn_working$Animal_24hr_Mass_g <- as.numeric(spawn_working$Animal_24hr_Mass_g)

#total count of urchins spawned for each ecosystem state
total_counts<- spawn_working%>%
  data.frame()%>%
  filter(!is.na(State),!is.na(Sex))%>%
  filter(Spawn_Mass_false != 0)%>%
  filter(!(Sex == "NA"))%>%
  mutate(State=factor(State),
         Sex=factor(Sex))%>%
  group_by(State)%>%
  summarize(total_n=n())

#total count of each sex spawned for each ecosystem state
sex_counts<- spawn_working%>%
  data.frame()%>%
  filter(!is.na(State),!is.na(Sex))%>%
  filter(Spawn_Mass_false != 0)%>%
  filter(!(Sex == "NA"))%>%
  mutate(State=factor(State),
         Sex=factor(Sex))%>%
  group_by(State,Sex)%>%
  summarize(sex_counts=n(),.groups="drop")%>%
  left_join(total_counts,by="State")%>%
  mutate(prop= (sex_counts/total_n))


# Figures -----------------------------------------------------------------

#boxplot overall spawn mass per ecosystem state
ggplot(data = spawn_working, aes(x = State, y = Spawn_Mass_false)) +
  geom_boxplot() +
  #facet_wrap(~ State) +
  labs(
    x = "Ecosystem State",
    y = "Mass (g) ",
    fill = "Sex") +
  theme_classic()

#boxplot spawn mass per ecosystem state by sex
ggplot(data = spawn_working %>% 
         filter(Spawn_Mass_false != 0)
       , aes(x = State, y = Spawn_Mass_g)) +
  geom_boxplot() +
  facet_wrap(~ Sex) +
  labs(
    x = "Ecosystem State",
    y = "Mass (g) ",
    fill = "Sex") +
  theme_classic()

#boxplot 24hr urchin mass per ecosystem state by sex
ggplot(data = spawn_working, aes(x = Sex, y = Animal_24hr_Mass_g, fill = Sex)) +
  geom_boxplot() +
  facet_wrap(~ State) +
  labs(
    x = "Sex",
    y = "Mass (g) ",
    fill = "Sex") +
  theme_classic()

#sex distribution for each state 
ggplot(data = sex_counts
       , aes(x = State, y = prop, fill=Sex)) +
  geom_col() +
  #facet_wrap(~ Sex) +
  labs(
    x = "Ecosystem State",
    y = "Percent (%) ",
    fill = "Sex") +
  theme_classic()

#urchin sizefq for site/type/zone (just first zone in dataset)
ggplot(urchin_sizefq [1:18,],
       aes(x=size_cm, y=count))+
  geom_col(fill="skyblue3")+
  theme_classic()

#urchin sizefq for site/type/zone 
ggplot(urchin_sizefq_1,
       aes(x=size_cm, y=count))+
  facet_wrap(~site_id)+
  geom_col(fill="skyblue3")+
  theme_classic()

#urchin diameter x gonad mass for site/type/zone (y=mx+b)
ggplot(gonad_working,
       aes(x=Test_Diameter_mm, y=Gonad_Mass_g))+
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

#y=mx+b for the first set of site/type/zone 
regression <- lm(Gonad_Mass_g ~ Test_Diameter_mm, data = gonad_working [1:27,])
summary(regression)

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

