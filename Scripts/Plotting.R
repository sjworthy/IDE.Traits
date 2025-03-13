# Plotting and basic stats

# load libraries
library(tidyverse)

#### Load imputed data without woody species ####

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)

# imputed traits with woody species
imputed.NW = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.outliersRM.csv", row.names = 1) %>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.NW.annual = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.annuals.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.NW.perennial = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennials.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.NW.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.forbs.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.NW.grass = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.grass.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.NW.annual.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.annual.forbs.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.NW.perennial.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennial.forb.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.NW.perennial.grass = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennial.grass.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)

#### Drop NAs ####
# drop NAs so data matches exactly that used in the models

imputed.NW.2 = imputed.NW %>%
  select(cover.change,site_code,Taxon, leafN.final:mean.MAP) %>%
  drop_na()

# number of sites

unique(imputed.NW.2$site_code) # 63 sites

#### Woody species excluded ####

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)

woody = imputed.traits %>%
  filter(functional_group %in% c("WOODY") | local_lifeform %in% c("SHRUB","TREE","WOODY"))
# 112 populations removed, 81 species removed prior to outlier removal

# outlier removal
imputed.traits.2 = imputed.traits %>%
  filter(round(leafN.final, 5) < 46.00560) %>% # 4
  filter(round(height.final, 5) < 7.787850) %>% # 11
  filter(round(rootN.final, 5) < 25.55616) %>% # 18
  filter(round(SLA.final, 5) < 49.05439) %>% # 8
  filter(round(root.depth.final, 5) < 2.915000) %>% # 15
  filter(round(RTD.final, 5) < 0.7460162) %>% # 11
  filter(round(SRL.final, 5) < 437.3523) %>% # 21
  filter(round(rootDiam.final, 5) < 1.200000) %>% # 6
  filter(round(RMF.final, 5) < 0.7670084) # 10

woody.2 = imputed.traits.2 %>%
  filter(functional_group %in% c("WOODY") | local_lifeform %in% c("SHRUB","TREE","WOODY"))
# 81 populations remained, 62 species remained after outlier removal.

# 28% of populations removed, 23% of species removed as outliers
# out of all outliers removed (94 populations) 50 were woody (53%)



imputed.traits.NW= read.csv("./Formatted.Data/Revisions/imputed.traits.NW.final.csv", row.names = 1)



