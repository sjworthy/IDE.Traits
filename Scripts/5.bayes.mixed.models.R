## Script to evaluate linear mixed effects models of cover change and traits
# cover outliers removed

# load libraries 
library(tidyverse)
library(brms)

#### Load imputed data with woody species ####

imputed.traits = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.outliersRM.csv", row.names = 1)

hist(imputed.traits$cover.change)
boxplot(imputed.traits$cover.change)
mean = mean(imputed.traits$cover.change, na.rm = TRUE)
std = sd(imputed.traits$cover.change, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$cover.change[which(imputed.traits$cover.change <Tmin | imputed.traits$cover.change > Tmax)])
# -79.47619 - -24.50000 and 23.52000 - 55.00000
# percent removed
19/813*100 # 2.34%

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)

# imputed traits with woody species
imputed = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.outliersRM.csv", row.names = 1) %>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.annual = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.annuals.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.perennial = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.perennials.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.forbs.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.grass = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.grass.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.annual.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.annual.forbs.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.perennial.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.perennial.forb.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.perennial.grass = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.perennial.grass.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)

#### imputed traits model ####

priors <- c(prior(normal(0, 10), class = b))

imputed.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                             SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                           family = gaussian(),
                           prior = priors,
                           data = imputed)

summary(imputed.traits.model)
# nothing significant

imputed.traits.interact = brm(cover.change ~ RTD.final*SRL.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                           family = gaussian(),
                           prior = priors,
                           data = imputed)

summary(imputed.traits.interact)
# nothing significant

#### imputed traits lifespan model ####

piors <- c(prior(normal(0, 10), class = b))

annual.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                            SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                          family = gaussian(),
                          prior = priors,
                          data = imputed.annual)

summary(annual.traits.model)
# MAP


annual.traits.interact = brm(cover.change ~ SRL.final*RTD.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                          family = gaussian(),
                          prior = priors,
                          data = imputed.annual)

summary(annual.traits.interact)
# nothing significant

perennial.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                               SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                             family = gaussian(),
                             prior = priors,
                             data = imputed.perennial)

summary(perennial.traits.model)
# nothing significant

perennial.traits.interact = brm(cover.change ~ SRL.final*RTD.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                             family = gaussian(),
                             prior = priors,
                             data = imputed.perennial)

summary(perennial.traits.interact)
# nothing significant

#### impute traits functional group ####

priors <- c(prior(normal(0, 10), class = b))

forb.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                          SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                        family = gaussian(),
                        prior = priors,
                        data = imputed.forb)

summary(forb.traits.model)
# leafN

forb.traits.interact = brm(cover.change ~ SRL.final*RTD.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                        family = gaussian(),
                        prior = priors,
                        data = imputed.forb)

summary(forb.traits.interact)
# nothing significant

grass.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                           SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                         family = gaussian(),
                         prior = priors,
                         data = imputed.grass)

summary(grass.traits.model)
# nothing significant

grass.traits.interact = brm(cover.change ~ SRL.final*RTD.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                         family = gaussian(),
                         prior = priors,
                         data = imputed.grass)

summary(grass.traits.interact)
# nothing significant

#### imputed traits lifespan x functional group ####

priors <- c(prior(normal(0, 10), class = b))

annual.forb.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                 SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                               family = gaussian(),
                               prior = priors,
                               data = imputed.annual.forb)

summary(annual.forb.traits.model)
# nothing significant

annual.forb.traits.interact = brm(cover.change ~ SRL.final*RTD.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                               family = gaussian(),
                               prior = priors,
                               data = imputed.annual.forb)

summary(annual.forb.traits.interact)
# nothing significant

perennial.forb.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                    SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                  family = gaussian(),
                                  prior = priors,
                                  data = imputed.perennial.forb)

summary(perennial.forb.traits.model)
# nothing significant

perennial.forb.traits.interact = brm(cover.change ~ SRL.final*RTD.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                  family = gaussian(),
                                  prior = priors,
                                  data = imputed.perennial.forb)

summary(perennial.forb.traits.interact)
# nothing significant

perennial.grass.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                     SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                   family = gaussian(),
                                   prior = priors,
                                   data = imputed.perennial.grass)

summary(perennial.grass.traits.model)
# nothing significant

perennial.grass.traits.interact = brm(cover.change ~ SRL.final*RTD.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                   family = gaussian(),
                                   prior = priors,
                                   data = imputed.perennial.grass)

summary(perennial.grass.traits.interact)
# nothing significant


#### Load imputed data without woody species ####

imputed.traits = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.outliersRM.csv", row.names = 1)

hist(imputed.traits$cover.change)
boxplot(imputed.traits$cover.change)
mean = mean(imputed.traits$cover.change, na.rm = TRUE)
std = sd(imputed.traits$cover.change, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$cover.change[which(imputed.traits$cover.change <Tmin | imputed.traits$cover.change > Tmax)])
# -79.47619 - -24.50000 and 23.52000 - 55.00000
# percent removed
18/719*100 # 2.50%

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

#### imputed traits model NW ####

priors <- c(prior(normal(0, 10), class = b))

imputed.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                              family = gaussian(),
                              prior = priors,
                              data = imputed.NW)

summary(imputed.traits.NW.model)
# leafN

imputed.traits.NW.interact = brm(cover.change ~ SRL.final*RTD.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                              family = gaussian(),
                              prior = priors,
                              data = imputed.NW)

summary(imputed.traits.NW.interact)
# nothing significant

#### imputed traits lifespan model NW ####

piors <- c(prior(normal(0, 10), class = b))

annual.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                               SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                             family = gaussian(),
                             prior = priors,
                             data = imputed.NW.annual)

summary(annual.traits.NW.model)
# MAP

annual.traits.NW.interact = brm(cover.change ~ SRL.final*RTD.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                             family = gaussian(),
                             prior = priors,
                             data = imputed.NW.annual)

summary(annual.traits.NW.interact)
# MAP

perennial.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                  SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                family = gaussian(),
                                prior = priors,
                                data = imputed.NW.perennial)

summary(perennial.traits.NW.model)
# nothing significant

perennial.traits.NW.interact = brm(cover.change ~ SRL.final*RTD.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                family = gaussian(),
                                prior = priors,
                                data = imputed.NW.perennial)

summary(perennial.traits.NW.interact)
# nothing significant

#### impute traits functional group NW ####

priors <- c(prior(normal(0, 10), class = b))

forb.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                             SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                           family = gaussian(),
                           prior = priors,
                           data = imputed.NW.forb)

summary(forb.traits.NW.model)
# leafN significant

forb.traits.NW.interact = brm(cover.change ~ SRL.final*RTD.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                           family = gaussian(),
                           prior = priors,
                           data = imputed.NW.forb)

summary(forb.traits.NW.interact)
# RTD*SRL significant

conditional_effects(forb.traits.NW.interact)

grass.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                              SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                            family = gaussian(),
                            prior = priors,
                            data = imputed.NW.grass)

summary(grass.traits.NW.model)
# nothing significant

grass.traits.NW.interact = brm(cover.change ~ SRL.final*RTD.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                            family = gaussian(),
                            prior = priors,
                            data = imputed.NW.grass)

summary(grass.traits.NW.interact)
# nothing significant

#### imputed traits lifespan x functional group NW ####

priors <- c(prior(normal(0, 10), class = b))

annual.forb.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                    SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                  family = gaussian(),
                                  prior = priors,
                                  data = imputed.NW.annual.forb)

summary(annual.forb.traits.NW.model)
# MAP significant

annual.forb.traits.NW.interact = brm(cover.change ~ SRL.final*RTD.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                  family = gaussian(),
                                  prior = priors,
                                  data = imputed.NW.annual.forb)

summary(annual.forb.traits.NW.interact)
# MAP significant

perennial.forb.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                       SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                     family = gaussian(),
                                     prior = priors,
                                     data = imputed.NW.perennial.forb)

summary(perennial.forb.traits.NW.model)
# nothing significant

perennial.forb.traits.NW.interact= brm(cover.change ~ SRL.final*RTD.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                     family = gaussian(),
                                     prior = priors,
                                     data = imputed.NW.perennial.forb)

summary(perennial.forb.traits.NW.interact)
# nothing significant

perennial.grass.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                        SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                      family = gaussian(),
                                      prior = priors,
                                      data = imputed.NW.perennial.grass)

summary(perennial.grass.traits.NW.model)
# nothing significant

perennial.grass.traits.NW.interact = brm(cover.change ~ SRL.final*RTD.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                      family = gaussian(),
                                      prior = priors,
                                      data = imputed.NW.perennial.grass)

summary(perennial.grass.traits.NW.interact)
# nothing significant

#### Load CC data with woody species ####

all.traits.CC = read.csv("./Formatted.Data/Revisions/Final.Data/final.data.CC.outliersRM.csv", row.names = 1)

hist(all.traits.CC$cover.change)
boxplot(all.traits.CC$cover.change)
mean = mean(all.traits.CC$cover.change, na.rm = TRUE)
std = sd(all.traits.CC$cover.change, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(all.traits.CC$cover.change[which(all.traits.CC$cover.change <Tmin | all.traits.CC$cover.change > Tmax)])
# -79.47619 -54.15358 -38.00000  28.91667
# percent removed
4/211*100 #1.89%

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)

# imputed traits with woody species
all.CC = read.csv("./Formatted.Data/Revisions/Final.Data/final.data.CC.outliersRM.csv", row.names = 1) %>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) < 28.91667) %>%
  mutate_at(vars(4:12,16,17),scale)
perennial.CC = read.csv("./Formatted.Data/Revisions/Final.Data/final.data.CC.perennials.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) < 28.91667) %>%
  mutate_at(vars(4:12,16,17),scale)
forb.CC = read.csv("./Formatted.Data/Revisions/Final.Data/final.data.CC.forbs.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) < 28.91667) %>%
  mutate_at(vars(4:12,16,17),scale)
grass.CC = read.csv("./Formatted.Data/Revisions/Final.Data/final.data.CC.grass.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) < 28.91667) %>%
  mutate_at(vars(4:12,16,17),scale)

#### CC traits ####
priors <- c(prior(normal(0, 10), class = b))

CC.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                 SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
               family = gaussian(),
               prior = priors,
               data = all.CC)

summary(CC.model)
# nothing significant

CC.model.interact = brm(cover.change ~ SRL.groot.cahill.merge*RTD.groot.cahill.merge + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
               family = gaussian(),
               prior = priors,
               data = all.CC)

summary(CC.model.interact)
# nothing significant

#### CC traits lifespan ####

priors <- c(prior(normal(0, 10), class = b))

perennial.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                        SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                      family = gaussian(),
                      prior = priors,
                      data = perennial.CC)

summary(perennial.model)
# RTD

perennial.model.interact = brm(cover.change ~ SRL.groot.cahill.merge*RTD.groot.cahill.merge + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                      family = gaussian(),
                      prior = priors,
                      data = perennial.CC)

summary(perennial.model.interact)
# RTD

#### CC traits functional group ####

priors <- c(prior(normal(0, 10), class = b))

forb.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                   SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                 family = gaussian(),
                 prior = priors,
                 data = forb.CC)

summary(forb.model)
# nothing significant

forb.model.interact = brm(cover.change ~ SRL.groot.cahill.merge*RTD.groot.cahill.merge + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                 family = gaussian(),
                 prior = priors,
                 data = forb.CC)

summary(forb.model.interact)
# nothing significant

grass.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                    SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                  family = gaussian(),
                  prior = priors,
                  data = grass.CC)

summary(grass.model)
# nothing significant

grass.model.interact = brm(cover.change ~ SRL.groot.cahill.merge*RTD.groot.cahill.merge + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                  family = gaussian(),
                  prior = priors,
                  data = grass.CC)

summary(grass.model.interact)
# nothing significant

#### Load CC data without woody species ####

all.traits.CC = read.csv("./Formatted.Data/Revisions/Final.Data/final.data.CC.NW.outliersRM.csv", row.names = 1)

hist(all.traits.CC$cover.change)
boxplot(all.traits.CC$cover.change)
mean = mean(all.traits.CC$cover.change, na.rm = TRUE)
std = sd(all.traits.CC$cover.change, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(all.traits.CC$cover.change[which(all.traits.CC$cover.change <Tmin | all.traits.CC$cover.change > Tmax)])
# -79.47619 -54.15358 -38.00000  28.91667
# percent removed
4/204*100 #1.96%

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)

# imputed traits with woody species
all.CC.NW = read.csv("./Formatted.Data/Revisions/Final.Data/final.data.CC.NW.outliersRM.csv", row.names = 1) %>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) < 28.91667) %>%
  mutate_at(vars(4:12,16,17),scale)
perennial.CC.NW = read.csv("./Formatted.Data/Revisions/Final.Data/final.data.CC.NW.perennials.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) < 28.91667) %>%
  mutate_at(vars(4:12,16,17),scale)
forb.CC.NW = read.csv("./Formatted.Data/Revisions/Final.Data/final.data.CC.NW.forbs.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) < 28.91667) %>%
  mutate_at(vars(4:12,16,17),scale)
grass.CC.NW = read.csv("./Formatted.Data/Revisions/Final.Data/final.data.CC.NW.grass.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) < 28.91667) %>%
  mutate_at(vars(4:12,16,17),scale)

#### CC traits NW ####

priors <- c(prior(normal(0, 10), class = b))

CC.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                 SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
               family = gaussian(),
               prior = priors,
               data = all.CC.NW)

summary(CC.model)
# nothing significant

CC.model.interact = brm(cover.change ~ SRL.groot.cahill.merge*RTD.groot.cahill.merge + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
               family = gaussian(),
               prior = priors,
               data = all.CC.NW)

summary(CC.model.interact)
# nothing significant

#### CC traits lifespan NW ####

priors <- c(prior(normal(0, 10), class = b))

perennial.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                        SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                      family = gaussian(),
                      prior = priors,
                      data = perennial.CC.NW)

summary(perennial.model)
# nothing significant

perennial.model.interact = brm(cover.change ~ SRL.groot.cahill.merge*RTD.groot.cahill.merge + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                      family = gaussian(),
                      prior = priors,
                      data = perennial.CC.NW)

summary(perennial.model.interact)
# nothing significant

#### CC traits NW functional group ####

priors <- c(prior(normal(0, 10), class = b))

forb.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                   SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                 family = gaussian(),
                 prior = priors,
                 data = forb.CC.NW)

summary(forb.model)
# nothing significant

forb.model.interact = brm(cover.change ~ SRL.groot.cahill.merge*RTD.groot.cahill.merge + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                 family = gaussian(),
                 prior = priors,
                 data = forb.CC.NW)

summary(forb.model.interact)
# nothing significant

grass.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                    SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                  family = gaussian(),
                  prior = priors,
                  data = grass.CC.NW)

summary(grass.model)
# nothing significant

grass.model.interact = brm(cover.change ~ SRL.groot.cahill.merge*RTD.groot.cahill.merge  + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                  family = gaussian(),
                  prior = priors,
                  data = grass.CC.NW)

summary(grass.model.interact)
# nothing significant


#### Load original CC data with woody species ####

all.traits.CC = read.csv("./Formatted.Data/Revisions/final.data.CC.csv", row.names = 1)

hist(all.traits.CC$cover.change)
boxplot(all.traits.CC$cover.change)
mean = mean(all.traits.CC$cover.change, na.rm = TRUE)
std = sd(all.traits.CC$cover.change, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(all.traits.CC$cover.change[which(all.traits.CC$cover.change <Tmin | all.traits.CC$cover.change > Tmax)])
# -79.47619 -54.15358 -38.00000  28.91667
# percent removed
4/236*100 #1.69%

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)

# imputed traits with woody species
all.CC.OG = read.csv("./Formatted.Data/Revisions/final.data.CC.csv", row.names = 1) %>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) < 28.91667) %>%
  mutate_at(vars(4:12,16,17),scale)
perennial.CC.OG = all.traits.CC %>%
  filter(local_lifespan == "PERENNIAL") %>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) < 28.91667) %>%
  mutate_at(vars(4:12,16,17),scale)
forb.CC.OG = all.traits.CC %>%
  filter(functional_group == "FORB") %>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) < 28.91667) %>%
  mutate_at(vars(4:12,16,17),scale)
grass.CC.OG = all.traits.CC %>%
  filter(functional_group == "GRASS") %>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) < 28.91667) %>%
  mutate_at(vars(4:12,16,17),scale)

#### CC traits OG ####

priors <- c(prior(normal(0, 10), class = b))

CC.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                 SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
               family = gaussian(),
               prior = priors,
               data = all.CC.OG)

summary(CC.model)
# MAP, RTD

#### CC traits lifespan OG ####

priors <- c(prior(normal(0, 10), class = b))

perennial.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                        SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                      family = gaussian(),
                      prior = priors,
                      data = perennial.CC.OG)

summary(perennial.model)
# RTD

#### CC traits OG functional group ####

priors <- c(prior(normal(0, 10), class = b))

forb.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                   SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                 family = gaussian(),
                 prior = priors,
                 data = forb.CC.OG)

summary(forb.model)
# nothing significant

grass.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                    SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                  family = gaussian(),
                  prior = priors,
                  data = grass.CC.OG)

summary(grass.model)
# nothing significant




#### Load original CC data without woody species ####

all.traits.CC = read.csv("./Formatted.Data/Revisions/final.data.CC.NW.csv", row.names = 1)

hist(all.traits.CC$cover.change)
boxplot(all.traits.CC$cover.change)
mean = mean(all.traits.CC$cover.change, na.rm = TRUE)
std = sd(all.traits.CC$cover.change, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(all.traits.CC$cover.change[which(all.traits.CC$cover.change <Tmin | all.traits.CC$cover.change > Tmax)])
# -79.47619 -54.15358 -38.00000  28.91667
# percent removed
4/224*100 #1.79%

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)

# imputed traits with woody species
all.CC.OG = read.csv("./Formatted.Data/Revisions/final.data.CC.NW.csv", row.names = 1) %>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) < 28.91667) %>%
  mutate_at(vars(4:12,16,17),scale)
perennial.CC.OG = all.traits.CC %>%
  filter(local_lifespan == "PERENNIAL") %>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) < 28.91667) %>%
  mutate_at(vars(4:12,16,17),scale)
forb.CC.OG = all.traits.CC %>%
  filter(functional_group == "FORB") %>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) < 28.91667) %>%
  mutate_at(vars(4:12,16,17),scale)
grass.CC.OG = all.traits.CC %>%
  filter(functional_group == "GRASS") %>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) < 28.91667) %>%
  mutate_at(vars(4:12,16,17),scale)

#### CC traits OG ####

priors <- c(prior(normal(0, 10), class = b))

CC.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                 SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
               family = gaussian(),
               prior = priors,
               data = all.CC.OG)

summary(CC.model)
# nothing significant

#### CC traits lifespan OG ####

priors <- c(prior(normal(0, 10), class = b))

perennial.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                        SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                      family = gaussian(),
                      prior = priors,
                      data = perennial.CC.OG)

summary(perennial.model)
# nothing significant

#### CC traits OG functional group ####

priors <- c(prior(normal(0, 10), class = b))

forb.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                   SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                 family = gaussian(),
                 prior = priors,
                 data = forb.CC.OG)

summary(forb.model)
# nothing significant

grass.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                    SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                  family = gaussian(),
                  prior = priors,
                  data = grass.CC.OG)

summary(grass.model)
# nothing significant