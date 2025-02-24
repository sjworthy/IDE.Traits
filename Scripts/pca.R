library(funspace)
library(tidyverse)
library(ggbiplot)

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)

# imputed traits with woody species
imputed = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.outliersRM.csv", row.names = 1) %>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)

### PCA Ellipses ###

lifspan = imputed %>%
  filter(local_lifespan %in% c("ANNUAL","PERENNIAL"))

pc.traits = prcomp(lifspan[,c(26:34)], scale  = TRUE, center = TRUE)
summary(pc.traits)

ggbiplot(pc.traits, groups = lifspan$local_lifespan, ellipse = TRUE, labels.size = 4, varname.adjust = 1.2, varname.size = 5) +
  labs(fill = "local_lifespan", color = "local_lifespan")+
  theme_classic(base_size = 15)

functional_group = imputed %>%
  filter(functional_group %in% c("FORB","GRASS"))

pc.traits.2 = prcomp(functional_group[,c(26:34)], scale  = TRUE, center = TRUE)
summary(pc.traits.2)

ggbiplot(pc.traits.2, groups = functional_group$functional_group, ellipse = TRUE, labels.size = 4, varname.adjust = 1.2, varname.size = 5) +
  labs(fill = "functional_group", color = "functional_group")+
  theme_classic(base_size = 15)

#### PCA plots ####

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)

# imputed traits with woody species
imputed = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.outliersRM.csv", row.names = 1) %>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.annual = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.annuals.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.perennial = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.perennials.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.forbs.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.grass = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.grass.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.annual.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.annual.forbs.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.perennial.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.perennial.forb.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.perennial.grass = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.perennial.grass.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) 

pc.traits.all = prcomp(imputed[,c(26:34)], scale  = TRUE, center = TRUE)
summary(pc.traits.all)

ggbiplot(pc.traits.all, labels.size = 4,varname.adjust = 1.2, varname.size = 5, varname.color = "red") +
  theme_classic(base_size = 15)

pc.traits.annual = prcomp(imputed.annual[,c(26:34)], scale  = TRUE, center = TRUE)
summary(pc.traits.annual)

ggbiplot(pc.traits.annual, labels.size = 4,varname.adjust = 1.2, varname.size = 5, varname.color = "red") +
  theme_classic(base_size = 15)

pc.traits.perennial = prcomp(imputed.perennial[,c(26:34)], scale  = TRUE, center = TRUE)
summary(pc.traits.perennial)

ggbiplot(pc.traits.perennial, labels.size = 4,varname.adjust = 1.2, varname.size = 5, varname.color = "red") +
  theme_classic(base_size = 15)


  