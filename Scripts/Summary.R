# Get Summary stats

# load libraries
library(tidyverse)
library(corrplot)

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
imputed.NW.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.forbs.legume.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.NW.graminoid = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.graminoid.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.NW.annual.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.annual.forb.legume.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.NW.perennial.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennial.forb.legume.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.NW.perennial.graminoid = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennial.graminoid.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)

#### Drop NAs ####
# drop NAs so data matches exactly that used in the models

imputed.NW.2 = imputed.NW %>%
  select(cover.change,site_code,Taxon, leafN.final:mean.MAP) %>%
  drop_na()

imputed.NW.2 = imputed.NW %>%
  select(cover.change,site_code,Taxon, leafN.final:mean.MAP, local_lifeform,
         functional_group) %>%
  drop_na()

# number of sites

unique(imputed.NW.2$site_code) # 63 sites

# number of species

unique(imputed.NW.2$Taxon) # 421 species

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

#### Load imputed data without woody species not scaled ####

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)

# imputed traits with woody species
imputed.NW = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.outliersRM.csv", row.names = 1) %>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.NW.annual = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.annuals.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.NW.perennial = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennials.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.NW.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.forbs.legume.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.NW.graminoid = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.graminoid.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.NW.annual.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.annual.forb.legume.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.NW.perennial.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennial.forb.legume.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.NW.perennial.graminoid = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennial.graminoid.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)

# Remove NAs

imputed.NW.2 = imputed.NW %>%
  dplyr::select(cover.change,site_code,Taxon, leafN.final:mean.MAP) %>%
  drop_na()
imputed.NW.annual.2 = imputed.NW.annual %>%
  dplyr::select(cover.change,site_code,Taxon, functional_group,leafN.final:mean.MAP) %>%
  drop_na()
imputed.NW.perennial.2 = imputed.NW.perennial %>%
  dplyr::select(cover.change,site_code,Taxon,functional_group,leafN.final:mean.MAP) %>%
  drop_na()
imputed.NW.graminoid.2 = imputed.NW.graminoid %>%
  dplyr::select(cover.change,site_code,Taxon, leafN.final:mean.MAP) %>%
  drop_na()
imputed.NW.forb.2 = imputed.NW.forb %>%
  dplyr::select(cover.change,site_code,Taxon, leafN.final:mean.MAP) %>%
  drop_na()
imputed.NW.annual.forb.2 = imputed.NW.annual.forb %>%
  dplyr::select(cover.change,site_code,Taxon, leafN.final:mean.MAP) %>%
  drop_na()
imputed.NW.perennial.forb.2 = imputed.NW.perennial.forb%>%
  dplyr::select(cover.change,site_code,Taxon, leafN.final:mean.MAP) %>%
  drop_na()
imputed.NW.perennial.graminoid.2 = imputed.NW.perennial.graminoid %>%
  dplyr::select(cover.change,site_code,Taxon, leafN.final:mean.MAP) %>%
  drop_na()

#### Table S2-S6 ####
mean(imputed.NW.2$leafN.final, na.rm = TRUE)
mean(imputed.NW.2$height.final, na.rm = TRUE)
mean(imputed.NW.2$rootN.final, na.rm = TRUE)
mean(imputed.NW.2$SLA.final, na.rm = TRUE)
mean(imputed.NW.2$root.depth.final, na.rm = TRUE)
mean(imputed.NW.2$RTD.final, na.rm = TRUE)
mean(imputed.NW.2$SRL.final, na.rm = TRUE)
mean(imputed.NW.2$rootDiam.final, na.rm = TRUE)
mean(imputed.NW.2$RMF.final, na.rm = TRUE)
mean(imputed.NW.2$mean.MAP, na.rm = TRUE)
mean(imputed.NW.2$mean.DSI, na.rm = TRUE)
mean(imputed.NW.2$cover.change)

sd(imputed.NW.2$leafN.final, na.rm = TRUE)
sd(imputed.NW.2$height.final, na.rm = TRUE)
sd(imputed.NW.2$rootN.final, na.rm = TRUE)
sd(imputed.NW.2$SLA.final, na.rm = TRUE)
sd(imputed.NW.2$RMF.final, na.rm = TRUE)
sd(imputed.NW.2$root.depth.final, na.rm = TRUE)
sd(imputed.NW.2$RTD.final, na.rm = TRUE)
sd(imputed.NW.2$SRL.final, na.rm = TRUE)
sd(imputed.NW.2$rootDiam.final, na.rm = TRUE)
sd(imputed.NW.2$mean.MAP, na.rm = TRUE)
sd(imputed.NW.2$mean.DSI, na.rm = TRUE)
sd(imputed.NW.2$cover.change)

range(imputed.NW.2$leafN.final, na.rm = TRUE)
range(imputed.NW.2$height.final, na.rm = TRUE)
range(imputed.NW.2$rootN.final, na.rm = TRUE)
range(imputed.NW.2$SLA.final, na.rm = TRUE)
range(imputed.NW.2$RMF.final, na.rm = TRUE)
range(imputed.NW.2$root.depth.final, na.rm = TRUE)
range(imputed.NW.2$RTD.final, na.rm = TRUE)
range(imputed.NW.2$SRL.final, na.rm = TRUE)
range(imputed.NW.2$rootDiam.final, na.rm = TRUE)
range(imputed.NW.2$mean.MAP, na.rm = TRUE)
range(imputed.NW.2$mean.DSI, na.rm = TRUE)
range(imputed.NW.2$cover.change)

sd(imputed.NW.2$leafN.final, na.rm = TRUE)/mean(imputed.NW.2$leafN.final, na.rm = TRUE)
sd(imputed.NW.2$height.final, na.rm = TRUE)/mean(imputed.NW.2$height.final, na.rm = TRUE)
sd(imputed.NW.2$rootN.final, na.rm = TRUE)/mean(imputed.NW.2$rootN.final, na.rm = TRUE)
sd(imputed.NW.2$SLA.final, na.rm = TRUE)/mean(imputed.NW.2$SLA.final, na.rm = TRUE)
sd(imputed.NW.2$RMF.final, na.rm = TRUE)/mean(imputed.NW.2$RMF.final, na.rm = TRUE)
sd(imputed.NW.2$root.depth.final, na.rm = TRUE)/mean(imputed.NW.2$root.depth.final, na.rm = TRUE)
sd(imputed.NW.2$RTD.final, na.rm = TRUE)/mean(imputed.NW.2$RTD.final, na.rm = TRUE)
sd(imputed.NW.2$SRL.final, na.rm = TRUE)/mean(imputed.NW.2$SRL.final, na.rm = TRUE)
sd(imputed.NW.2$rootDiam.final, na.rm = TRUE)/mean(imputed.NW.2$rootDiam.final, na.rm = TRUE)
sd(imputed.NW.2$mean.MAP, na.rm = TRUE)/mean(imputed.NW.2$mean.MAP, na.rm = TRUE)
sd(imputed.NW.2$mean.DSI, na.rm = TRUE)/mean(imputed.NW.2$mean.DSI, na.rm = TRUE)
sd(imputed.NW.2$cover.change)/mean(imputed.NW.2$cover.change)

mean(imputed.NW.annual.2$leafN.final, na.rm = TRUE)
mean(imputed.NW.annual.2$height.final, na.rm = TRUE)
mean(imputed.NW.annual.2$rootN.final, na.rm = TRUE)
mean(imputed.NW.annual.2$SLA.final, na.rm = TRUE)
mean(imputed.NW.annual.2$RMF.final, na.rm = TRUE)
mean(imputed.NW.annual.2$root.depth.final, na.rm = TRUE)
mean(imputed.NW.annual.2$RTD.final, na.rm = TRUE)
mean(imputed.NW.annual.2$SRL.final, na.rm = TRUE)
mean(imputed.NW.annual.2$rootDiam.final, na.rm = TRUE)
mean(imputed.NW.annual.2$mean.MAP, na.rm = TRUE)
mean(imputed.NW.annual.2$mean.DSI, na.rm = TRUE)
mean(imputed.NW.annual.2$cover.change)

sd(imputed.NW.annual.2$leafN.final, na.rm = TRUE)
sd(imputed.NW.annual.2$height.final, na.rm = TRUE)
sd(imputed.NW.annual.2$rootN.final, na.rm = TRUE)
sd(imputed.NW.annual.2$SLA.final, na.rm = TRUE)
sd(imputed.NW.annual.2$RMF.final, na.rm = TRUE)
sd(imputed.NW.annual.2$root.depth.final, na.rm = TRUE)
sd(imputed.NW.annual.2$RTD.final, na.rm = TRUE)
sd(imputed.NW.annual.2$SRL.final, na.rm = TRUE)
sd(imputed.NW.annual.2$rootDiam.final, na.rm = TRUE)
sd(imputed.NW.annual.2$mean.MAP, na.rm = TRUE)
sd(imputed.NW.annual.2$mean.DSI, na.rm = TRUE)
sd(imputed.NW.annual.2$cover.change)

range(imputed.NW.annual.2$leafN.final, na.rm = TRUE)
range(imputed.NW.annual.2$height.final, na.rm = TRUE)
range(imputed.NW.annual.2$rootN.final, na.rm = TRUE)
range(imputed.NW.annual.2$SLA.final, na.rm = TRUE)
range(imputed.NW.annual.2$RMF.final, na.rm = TRUE)
range(imputed.NW.annual.2$root.depth.final, na.rm = TRUE)
range(imputed.NW.annual.2$RTD.final, na.rm = TRUE)
range(imputed.NW.annual.2$SRL.final, na.rm = TRUE)
range(imputed.NW.annual.2$rootDiam.final, na.rm = TRUE)
range(imputed.NW.annual.2$mean.MAP, na.rm = TRUE)
range(imputed.NW.annual.2$mean.DSI, na.rm = TRUE)
range(imputed.NW.annual.2$cover.change)

mean(imputed.NW.perennial.2$leafN.final, na.rm = TRUE)
mean(imputed.NW.perennial.2$height.final, na.rm = TRUE)
mean(imputed.NW.perennial.2$rootN.final, na.rm = TRUE)
mean(imputed.NW.perennial.2$SLA.final, na.rm = TRUE)
mean(imputed.NW.perennial.2$RMF.final, na.rm = TRUE)
mean(imputed.NW.perennial.2$root.depth.final, na.rm = TRUE)
mean(imputed.NW.perennial.2$RTD.final, na.rm = TRUE)
mean(imputed.NW.perennial.2$SRL.final, na.rm = TRUE)
mean(imputed.NW.perennial.2$rootDiam.final, na.rm = TRUE)
mean(imputed.NW.perennial.2$mean.MAP, na.rm = TRUE)
mean(imputed.NW.perennial.2$mean.DSI, na.rm = TRUE)
mean(imputed.NW.perennial.2$cover.change)

sd(imputed.NW.perennial.2$leafN.final, na.rm = TRUE)
sd(imputed.NW.perennial.2$height.final, na.rm = TRUE)
sd(imputed.NW.perennial.2$rootN.final, na.rm = TRUE)
sd(imputed.NW.perennial.2$SLA.final, na.rm = TRUE)
sd(imputed.NW.perennial.2$RMF.final, na.rm = TRUE)
sd(imputed.NW.perennial.2$root.depth.final, na.rm = TRUE)
sd(imputed.NW.perennial.2$RTD.final, na.rm = TRUE)
sd(imputed.NW.perennial.2$SRL.final, na.rm = TRUE)
sd(imputed.NW.perennial.2$rootDiam.final, na.rm = TRUE)
sd(imputed.NW.perennial.2$mean.MAP, na.rm = TRUE)
sd(imputed.NW.perennial.2$mean.DSI, na.rm = TRUE)
sd(imputed.NW.perennial.2$cover.change)

range(imputed.NW.perennial.2$leafN.final, na.rm = TRUE)
range(imputed.NW.perennial.2$height.final, na.rm = TRUE)
range(imputed.NW.perennial.2$rootN.final, na.rm = TRUE)
range(imputed.NW.perennial.2$SLA.final, na.rm = TRUE)
range(imputed.NW.perennial.2$RMF.final, na.rm = TRUE)
range(imputed.NW.perennial.2$root.depth.final, na.rm = TRUE)
range(imputed.NW.perennial.2$RTD.final, na.rm = TRUE)
range(imputed.NW.perennial.2$SRL.final, na.rm = TRUE)
range(imputed.NW.perennial.2$rootDiam.final, na.rm = TRUE)
range(imputed.NW.perennial.2$mean.MAP, na.rm = TRUE)
range(imputed.NW.perennial.2$mean.DSI, na.rm = TRUE)
range(imputed.NW.perennial.2$cover.change)

mean(imputed.NW.graminoid.2$leafN.final, na.rm = TRUE)
mean(imputed.NW.graminoid.2$height.final, na.rm = TRUE)
mean(imputed.NW.graminoid.2$rootN.final, na.rm = TRUE)
mean(imputed.NW.graminoid.2$SLA.final, na.rm = TRUE)
mean(imputed.NW.graminoid.2$RMF.final, na.rm = TRUE)
mean(imputed.NW.graminoid.2$root.depth.final, na.rm = TRUE)
mean(imputed.NW.graminoid.2$RTD.final, na.rm = TRUE)
mean(imputed.NW.graminoid.2$SRL.final, na.rm = TRUE)
mean(imputed.NW.graminoid.2$rootDiam.final, na.rm = TRUE)
mean(imputed.NW.graminoid.2$mean.MAP, na.rm = TRUE)
mean(imputed.NW.graminoid.2$mean.DSI, na.rm = TRUE)
mean(imputed.NW.graminoid.2$cover.change)

sd(imputed.NW.graminoid.2$leafN.final, na.rm = TRUE)
sd(imputed.NW.graminoid.2$height.final, na.rm = TRUE)
sd(imputed.NW.graminoid.2$rootN.final, na.rm = TRUE)
sd(imputed.NW.graminoid.2$SLA.final, na.rm = TRUE)
sd(imputed.NW.graminoid.2$RMF.final, na.rm = TRUE)
sd(imputed.NW.graminoid.2$root.depth.final, na.rm = TRUE)
sd(imputed.NW.graminoid.2$RTD.final, na.rm = TRUE)
sd(imputed.NW.graminoid.2$SRL.final, na.rm = TRUE)
sd(imputed.NW.graminoid.2$rootDiam.final, na.rm = TRUE)
sd(imputed.NW.graminoid.2$mean.MAP, na.rm = TRUE)
sd(imputed.NW.graminoid.2$mean.DSI, na.rm = TRUE)
sd(imputed.NW.graminoid.2$cover.change)

range(imputed.NW.graminoid.2$leafN.final, na.rm = TRUE)
range(imputed.NW.graminoid.2$height.final, na.rm = TRUE)
range(imputed.NW.graminoid.2$rootN.final, na.rm = TRUE)
range(imputed.NW.graminoid.2$SLA.final, na.rm = TRUE)
range(imputed.NW.graminoid.2$RMF.final, na.rm = TRUE)
range(imputed.NW.graminoid.2$root.depth.final, na.rm = TRUE)
range(imputed.NW.graminoid.2$RTD.final, na.rm = TRUE)
range(imputed.NW.graminoid.2$SRL.final, na.rm = TRUE)
range(imputed.NW.graminoid.2$rootDiam.final, na.rm = TRUE)
range(imputed.NW.graminoid.2$mean.MAP, na.rm = TRUE)
range(imputed.NW.graminoid.2$mean.DSI, na.rm = TRUE)
range(imputed.NW.graminoid.2$cover.change)

mean(imputed.NW.forb.2$leafN.final, na.rm = TRUE)
mean(imputed.NW.forb.2$height.final, na.rm = TRUE)
mean(imputed.NW.forb.2$rootN.final, na.rm = TRUE)
mean(imputed.NW.forb.2$SLA.final, na.rm = TRUE)
mean(imputed.NW.forb.2$RMF.final, na.rm = TRUE)
mean(imputed.NW.forb.2$root.depth.final, na.rm = TRUE)
mean(imputed.NW.forb.2$RTD.final, na.rm = TRUE)
mean(imputed.NW.forb.2$SRL.final, na.rm = TRUE)
mean(imputed.NW.forb.2$rootDiam.final, na.rm = TRUE)
mean(imputed.NW.forb.2$mean.MAP, na.rm = TRUE)
mean(imputed.NW.forb.2$mean.DSI, na.rm = TRUE)
mean(imputed.NW.forb.2$cover.change)

sd(imputed.NW.forb.2$leafN.final, na.rm = TRUE)
sd(imputed.NW.forb.2$height.final, na.rm = TRUE)
sd(imputed.NW.forb.2$rootN.final, na.rm = TRUE)
sd(imputed.NW.forb.2$SLA.final, na.rm = TRUE)
sd(imputed.NW.forb.2$RMF.final, na.rm = TRUE)
sd(imputed.NW.forb.2$root.depth.final, na.rm = TRUE)
sd(imputed.NW.forb.2$RTD.final, na.rm = TRUE)
sd(imputed.NW.forb.2$SRL.final, na.rm = TRUE)
sd(imputed.NW.forb.2$rootDiam.final, na.rm = TRUE)
sd(imputed.NW.forb.2$mean.MAP, na.rm = TRUE)
sd(imputed.NW.forb.2$mean.DSI, na.rm = TRUE)
sd(imputed.NW.forb.2$cover.change)

range(imputed.NW.forb.2$leafN.final, na.rm = TRUE)
range(imputed.NW.forb.2$height.final, na.rm = TRUE)
range(imputed.NW.forb.2$rootN.final, na.rm = TRUE)
range(imputed.NW.forb.2$SLA.final, na.rm = TRUE)
range(imputed.NW.forb.2$RMF.final, na.rm = TRUE)
range(imputed.NW.forb.2$root.depth.final, na.rm = TRUE)
range(imputed.NW.forb.2$RTD.final, na.rm = TRUE)
range(imputed.NW.forb.2$SRL.final, na.rm = TRUE)
range(imputed.NW.forb.2$rootDiam.final, na.rm = TRUE)
range(imputed.NW.forb.2$mean.MAP, na.rm = TRUE)
range(imputed.NW.forb.2$mean.DSI, na.rm = TRUE)
range(imputed.NW.forb.2$cover.change)

mean(imputed.NW.annual.forb.2$leafN.final, na.rm = TRUE)
mean(imputed.NW.annual.forb.2$height.final, na.rm = TRUE)
mean(imputed.NW.annual.forb.2$rootN.final, na.rm = TRUE)
mean(imputed.NW.annual.forb.2$SLA.final, na.rm = TRUE)
mean(imputed.NW.annual.forb.2$RMF.final, na.rm = TRUE)
mean(imputed.NW.annual.forb.2$root.depth.final, na.rm = TRUE)
mean(imputed.NW.annual.forb.2$RTD.final, na.rm = TRUE)
mean(imputed.NW.annual.forb.2$SRL.final, na.rm = TRUE)
mean(imputed.NW.annual.forb.2$rootDiam.final, na.rm = TRUE)
mean(imputed.NW.annual.forb.2$mean.MAP, na.rm = TRUE)
mean(imputed.NW.annual.forb.2$mean.DSI, na.rm = TRUE)
mean(imputed.NW.annual.forb.2$cover.change)

sd(imputed.NW.annual.forb.2$leafN.final, na.rm = TRUE)
sd(imputed.NW.annual.forb.2$height.final, na.rm = TRUE)
sd(imputed.NW.annual.forb.2$rootN.final, na.rm = TRUE)
sd(imputed.NW.annual.forb.2$SLA.final, na.rm = TRUE)
sd(imputed.NW.annual.forb.2$RMF.final, na.rm = TRUE)
sd(imputed.NW.annual.forb.2$root.depth.final, na.rm = TRUE)
sd(imputed.NW.annual.forb.2$RTD.final, na.rm = TRUE)
sd(imputed.NW.annual.forb.2$SRL.final, na.rm = TRUE)
sd(imputed.NW.annual.forb.2$rootDiam.final, na.rm = TRUE)
sd(imputed.NW.annual.forb.2$mean.MAP, na.rm = TRUE)
sd(imputed.NW.annual.forb.2$mean.DSI, na.rm = TRUE)
sd(imputed.NW.annual.forb.2$cover.change)

range(imputed.NW.annual.forb.2$leafN.final, na.rm = TRUE)
range(imputed.NW.annual.forb.2$height.final, na.rm = TRUE)
range(imputed.NW.annual.forb.2$rootN.final, na.rm = TRUE)
range(imputed.NW.annual.forb.2$SLA.final, na.rm = TRUE)
range(imputed.NW.annual.forb.2$RMF.final, na.rm = TRUE)
range(imputed.NW.annual.forb.2$root.depth.final, na.rm = TRUE)
range(imputed.NW.annual.forb.2$RTD.final, na.rm = TRUE)
range(imputed.NW.annual.forb.2$SRL.final, na.rm = TRUE)
range(imputed.NW.annual.forb.2$rootDiam.final, na.rm = TRUE)
range(imputed.NW.annual.forb.2$mean.MAP, na.rm = TRUE)
range(imputed.NW.annual.forb.2$mean.DSI, na.rm = TRUE)
range(imputed.NW.annual.forb.2$cover.change)

mean(imputed.NW.perennial.graminoid.2$leafN.final, na.rm = TRUE)
mean(imputed.NW.perennial.graminoid.2$height.final, na.rm = TRUE)
mean(imputed.NW.perennial.graminoid.2$rootN.final, na.rm = TRUE)
mean(imputed.NW.perennial.graminoid.2$SLA.final, na.rm = TRUE)
mean(imputed.NW.perennial.graminoid.2$RMF.final, na.rm = TRUE)
mean(imputed.NW.perennial.graminoid.2$root.depth.final, na.rm = TRUE)
mean(imputed.NW.perennial.graminoid.2$RTD.final, na.rm = TRUE)
mean(imputed.NW.perennial.graminoid.2$SRL.final, na.rm = TRUE)
mean(imputed.NW.perennial.graminoid.2$rootDiam.final, na.rm = TRUE)
mean(imputed.NW.perennial.graminoid.2$mean.MAP, na.rm = TRUE)
mean(imputed.NW.perennial.graminoid.2$mean.DSI, na.rm = TRUE)
mean(imputed.NW.perennial.graminoid.2$cover.change)

sd(imputed.NW.perennial.graminoid.2$leafN.final, na.rm = TRUE)
sd(imputed.NW.perennial.graminoid.2$height.final, na.rm = TRUE)
sd(imputed.NW.perennial.graminoid.2$rootN.final, na.rm = TRUE)
sd(imputed.NW.perennial.graminoid.2$SLA.final, na.rm = TRUE)
sd(imputed.NW.perennial.graminoid.2$RMF.final, na.rm = TRUE)
sd(imputed.NW.perennial.graminoid.2$root.depth.final, na.rm = TRUE)
sd(imputed.NW.perennial.graminoid.2$RTD.final, na.rm = TRUE)
sd(imputed.NW.perennial.graminoid.2$SRL.final, na.rm = TRUE)
sd(imputed.NW.perennial.graminoid.2$rootDiam.final, na.rm = TRUE)
sd(imputed.NW.perennial.graminoid.2$mean.MAP, na.rm = TRUE)
sd(imputed.NW.perennial.graminoid.2$mean.DSI, na.rm = TRUE)
sd(imputed.NW.perennial.graminoid.2$cover.change)

range(imputed.NW.perennial.graminoid.2$leafN.final, na.rm = TRUE)
range(imputed.NW.perennial.graminoid.2$height.final, na.rm = TRUE)
range(imputed.NW.perennial.graminoid.2$rootN.final, na.rm = TRUE)
range(imputed.NW.perennial.graminoid.2$SLA.final, na.rm = TRUE)
range(imputed.NW.perennial.graminoid.2$RMF.final, na.rm = TRUE)
range(imputed.NW.perennial.graminoid.2$root.depth.final, na.rm = TRUE)
range(imputed.NW.perennial.graminoid.2$RTD.final, na.rm = TRUE)
range(imputed.NW.perennial.graminoid.2$SRL.final, na.rm = TRUE)
range(imputed.NW.perennial.graminoid.2$rootDiam.final, na.rm = TRUE)
range(imputed.NW.perennial.graminoid.2$mean.MAP, na.rm = TRUE)
range(imputed.NW.perennial.graminoid.2$mean.DSI, na.rm = TRUE)
range(imputed.NW.perennial.graminoid.2$cover.change)

mean(imputed.NW.perennial.forb.2$leafN.final, na.rm = TRUE)
mean(imputed.NW.perennial.forb.2$height.final, na.rm = TRUE)
mean(imputed.NW.perennial.forb.2$rootN.final, na.rm = TRUE)
mean(imputed.NW.perennial.forb.2$SLA.final, na.rm = TRUE)
mean(imputed.NW.perennial.forb.2$RMF.final, na.rm = TRUE)
mean(imputed.NW.perennial.forb.2$root.depth.final, na.rm = TRUE)
mean(imputed.NW.perennial.forb.2$RTD.final, na.rm = TRUE)
mean(imputed.NW.perennial.forb.2$SRL.final, na.rm = TRUE)
mean(imputed.NW.perennial.forb.2$rootDiam.final, na.rm = TRUE)
mean(imputed.NW.perennial.forb.2$mean.MAP, na.rm = TRUE)
mean(imputed.NW.perennial.forb.2$mean.DSI, na.rm = TRUE)
mean(imputed.NW.perennial.forb.2$cover.change)

sd(imputed.NW.perennial.forb.2$leafN.final, na.rm = TRUE)
sd(imputed.NW.perennial.forb.2$height.final, na.rm = TRUE)
sd(imputed.NW.perennial.forb.2$rootN.final, na.rm = TRUE)
sd(imputed.NW.perennial.forb.2$SLA.final, na.rm = TRUE)
sd(imputed.NW.perennial.forb.2$RMF.final, na.rm = TRUE)
sd(imputed.NW.perennial.forb.2$root.depth.final, na.rm = TRUE)
sd(imputed.NW.perennial.forb.2$RTD.final, na.rm = TRUE)
sd(imputed.NW.perennial.forb.2$SRL.final, na.rm = TRUE)
sd(imputed.NW.perennial.forb.2$rootDiam.final, na.rm = TRUE)
sd(imputed.NW.perennial.forb.2$mean.MAP, na.rm = TRUE)
sd(imputed.NW.perennial.forb.2$mean.DSI, na.rm = TRUE)
sd(imputed.NW.perennial.forb.2$cover.change)

range(imputed.NW.perennial.forb.2$leafN.final, na.rm = TRUE)
range(imputed.NW.perennial.forb.2$height.final, na.rm = TRUE)
range(imputed.NW.perennial.forb.2$rootN.final, na.rm = TRUE)
range(imputed.NW.perennial.forb.2$SLA.final, na.rm = TRUE)
range(imputed.NW.perennial.forb.2$RMF.final, na.rm = TRUE)
range(imputed.NW.perennial.forb.2$root.depth.final, na.rm = TRUE)
range(imputed.NW.perennial.forb.2$RTD.final, na.rm = TRUE)
range(imputed.NW.perennial.forb.2$SRL.final, na.rm = TRUE)
range(imputed.NW.perennial.forb.2$rootDiam.final, na.rm = TRUE)
range(imputed.NW.perennial.forb.2$mean.MAP, na.rm = TRUE)
range(imputed.NW.perennial.forb.2$mean.DSI, na.rm = TRUE)
range(imputed.NW.perennial.forb.2$cover.change)

# get functional groups in annuals and perennials

table(imputed.NW.annual.2$functional_group)
# forb 124, graminoid 44 legume 10
table(imputed.NW.perennial.2$functional_group)
# forb 224, graminoid 31, graminoid 174, legume 33

#### T-test of cover differences ####
t.test(imputed.NW.annual.2$cover.change,imputed.NW.perennial.2$cover.change)
# t = 0.59583, df = 324.25, p-value = 0.5517
t.test(imputed.NW.graminoid.2$cover.change,imputed.NW.forb.2$cover.change)
# t = -0.89215, df = 420.74, p-value = 0.3728

### Correlation Plots ####

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)

# imputed traits with woody species
imputed.NW = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.outliersRM.csv", row.names = 1) %>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.NW.annual = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.annuals.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.NW.perennial = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennials.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.NW.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.forbs.legume.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.NW.graminoid = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.graminoid.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.NW.annual.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.annual.forb.legume.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.NW.perennial.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennial.forb.legume.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.NW.perennial.graminoid = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennial.graminoid.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)

# Remove NAs

imputed.NW.2 = imputed.NW %>%
  select(cover.change,site_code,Taxon, leafN.final:mean.MAP) %>%
  drop_na()
imputed.NW.annual.2 = imputed.NW.annual %>%
  select(cover.change,site_code,Taxon, leafN.final:mean.MAP) %>%
  drop_na()
imputed.NW.perennial.2 = imputed.NW.perennial %>%
  select(cover.change,site_code,Taxon,leafN.final:mean.MAP) %>%
  drop_na()
imputed.NW.graminoid.2 = imputed.NW.graminoid %>%
  select(cover.change,site_code,Taxon, leafN.final:mean.MAP) %>%
  drop_na()
imputed.NW.forb.2 = imputed.NW.forb %>%
  select(cover.change,site_code,Taxon, leafN.final:mean.MAP) %>%
  drop_na()
imputed.NW.annual.forb.2 = imputed.NW.annual.forb %>%
  select(cover.change,site_code,Taxon, leafN.final:mean.MAP) %>%
  drop_na()
imputed.NW.perennial.forb.2 = imputed.NW.perennial.forb%>%
  select(cover.change,site_code,Taxon, leafN.final:mean.MAP) %>%
  drop_na()
imputed.NW.perennial.graminoid.2 = imputed.NW.perennial.graminoid %>%
  select(cover.change,site_code,Taxon, leafN.final:mean.MAP) %>%
  drop_na()

colnames(imputed.NW.2)[4:14] = c("Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter","SRL","RTD","RMF",
                                 "DSI","Precipitation")
colnames(imputed.NW.annual.2)[4:14] = c("Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter","SRL","RTD","RMF",
                                 "DSI","Precipitation")
colnames(imputed.NW.perennial.2)[4:14] = c("Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter","SRL","RTD","RMF",
                                 "DSI","Precipitation")
colnames(imputed.NW.forb.2)[4:14] = c("Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter","SRL","RTD","RMF",
                                 "DSI","Precipitation")
colnames(imputed.NW.graminoid.2)[4:14] = c("Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter","SRL","RTD","RMF",
                                 "DSI","Precipitation")
colnames(imputed.NW.annual.forb.2)[4:14] = c("Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter","SRL","RTD","RMF",
                                 "DSI","Precipitation")
colnames(imputed.NW.perennial.graminoid.2)[4:14] = c("Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter","SRL","RTD","RMF",
                                             "DSI","Precipitation")
colnames(imputed.NW.perennial.forb.2)[4:14] = c("Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter","SRL","RTD","RMF",
                                             "DSI","Precipitation")
# output 8 x 8
all.traits.cor = cor(imputed.NW.2[,c(4:14)],use = "pairwise") 
corrplot(all.traits.cor, method="number",tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")
# highest correlation is -0.28 and 0.34

write.csv(all.traits.cor, file = "./Results/all.traits.cor.csv")

annual.traits.cor = cor(imputed.NW.annual.2[,c(4:14)],use = "pairwise") 
corrplot(annual.traits.cor, method="number",tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")
# highest correlation is -0.31 and 0.35

write.csv(annual.traits.cor, file = "./Results/annual.traits.cor.csv")

perennial.traits.cor = cor(imputed.NW.perennial.2[,c(4:14)],use = "pairwise") 
corrplot(perennial.traits.cor, method="number",tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")
# highest correlation is -0.29 and 0.39

write.csv(perennial.traits.cor, file = "./Results/perennial.traits.cor.csv")

forb.traits.cor = cor(imputed.NW.forb.2[,c(4:14)],use = "pairwise") 
corrplot(forb.traits.cor, method="number",tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")
# highest correlation is -0.24 and 0.33

write.csv(forb.traits.cor, file = "./Results/forb.traits.cor.csv")

graminoid.traits.cor = cor(imputed.NW.graminoid.2[,c(4:14)],use = "pairwise") 
corrplot(graminoid.traits.cor, method="number",tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")
# highest correlation is -0.47 and 0.25

write.csv(graminoid.traits.cor, file = "./Results/graminoid.traits.cor.csv")

annual.forb.traits.cor = cor(imputed.NW.annual.forb.2[,c(4:14)],use = "pairwise") 
corrplot(annual.forb.traits.cor, method="number",tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")
# highest correlation is -0.43 and 0.35

write.csv(annual.forb.traits.cor, file = "./Results/annual.forb.traits.cor.csv")

perennial.forb.traits.cor = cor(imputed.NW.perennial.forb.2[,c(4:14)],use = "pairwise") 
corrplot(perennial.forb.traits.cor, method="number",tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")
# highest correlation is -0.27 and 0.38

write.csv(perennial.forb.traits.cor, file = "./Results/perennial.forb.traits.cor.csv")

perennial.graminoid.traits.cor = cor(imputed.NW.perennial.graminoid.2[,c(4:14)],use = "pairwise") 
corrplot(perennial.graminoid.traits.cor, method="number",tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")
# highest correlation is -0.44 and 0.25

write.csv(perennial.graminoid.traits.cor, file = "./Results/perennial.graminoid.traits.cor.csv")


#### Height of covers ####

site.data = read.csv("./Raw.Data/Site_Elev-Disturb.csv")

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)

# imputed traits with woody species
imputed.NW = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.outliersRM.csv", row.names = 1) %>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.NW.2 = imputed.NW %>%
  select(cover.change,site_code,Taxon, leafN.final:mean.MAP) %>%
  drop_na()

unq.sites = unique(imputed.NW.2$site_code)

site.data.2 = site.data %>%
  filter(site_code %in% unq.sites)

table(site.data.2$habitat.type)
# graminoidland Shrubland 
# 44        19 

