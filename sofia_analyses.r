
install.packages('librarian')
require(librarian)
librarian::shelf(tidyverse,here, janitor, googlesheets4, lubridate, splitstackshape,
                 googledrive,googlesheets4,httpuv,dplyr,ggplot2,pwr2)



# Data Sets ---------------------------------------------------------------

spawn_raw <- read_sheet("https://docs.google.com/spreadsheets/d/11DLr38iVRDcvWiDoBY1hl_Cn5M9oAxDwBRQsx_eMbHk/edit?usp=sharing") 

gonad_raw <- read_sheet("https://docs.google.com/spreadsheets/d/1Ih-hBXRtfXVMdxw5ibZnXy_dZErdcx5FfeKMSc0HEc4/edit?usp=sharing")

load("kelp_recovery_data.rda")

#delete unnecessary rows
spawn_working <- subset(spawn_raw, select = -c(Data_Enterer,Gonad_Mass_total,Spawn_Mass_total,Treatment))

gonad_working <- subset(gonad_raw, select = -c(Institution,Name_of_Data_Enterer,Date_Fixed,Treatment,Soft_Tissue_Mass_g,Notes,Date_Entered)) %>%
  filter(Species %in% c("PUR", "purple_urchin")) %>%
  mutate(Species = "purple_urchin")%>%
  filter(!is.na(Gonad_Mass_g)) %>% 
  mutate(site_id = paste(Site_Number, Transect)) %>% 
  filter(year(Date_Collected) == 2025,
         Gonad_Mass_g<20)

urchin_sizefq_1 <- urchin_sizefq %>%
  filter(species=="purple_urchin") %>% 
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

ggplot(urchin_sizefq_1,
       aes(x=size_cm, y=count))+
  facet_wrap(~site_id)+
  geom_col(fill="skyblue3")+
  theme_classic()

#urchin diameter x gonad mass for site/type/zone (just first zone in dataset)
ggplot(gonad_working,
       aes(x=Test_Diameter_mm, y=Gonad_Mass_g))+
  facet_wrap(~site_id)+
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "aquamarine4")+
  theme_classic()

# Stats -------------------------------------------------------------------

#y=mx+b for the first set of site/type/zone
regression <- lm(Gonad_Mass_g ~ Test_Diameter_mm, data = gonad_working [1:27,])
summary(regression)

#avg urchin density per site/type/zone 
avg_urchin_density <- quad_data %>% 
  mutate(site_id = paste(site, site_type, zone)) %>% 
  group_by(site, site_type, zone) %>% 
  summarize(avg_density = mean(purple_urchin_densitym2, na.rm = TRUE)) %>% 
  mutate(density80m2 = avg_density*80)


library(dplyr)
library(broom)

coeff_table <- gonad_working %>%
  group_by(site_id) %>%
  do(model = lm(Gonad_Mass_g ~ Test_Diameter_mm, data = .)) %>%
  tidy(model) %>%   # extract coefficients
  ungroup() %>%
  select(site_id, term, estimate, std.error, statistic, p.value)

coeff_table

library(dplyr)
library(tidyr)
library(broom)

coeff_table <- gonad_working %>%
  group_by(site_id) %>%
  nest() %>%
  mutate(model = map(data, ~ lm(Gonad_Mass_g ~ Test_Diameter_mm, data = .)),
         tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  select(site_id, term, estimate, std.error, statistic, p.value)

coeff_table
view(coeff_table)

coeff_wide <- coeff_table %>%
  select(site_id, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>% 
  rename(b = "(Intercept)", a = "Test_Diameter_mm")

coeff_wide
view(coeff_wide)

#sample from size distribution

