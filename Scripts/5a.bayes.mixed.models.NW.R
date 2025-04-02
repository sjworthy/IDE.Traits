## Script to evaluate linear mixed effects models of cover change and traits
# cover outliers removed

# load libraries 
library(tidyverse)
library(brms)

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
#saveRDS(imputed.traits.NW.model, file = "./Results/all.imputed.traits.no_woody.rds")
bayes_R2(imputed.traits.NW.model)
# R2 0.06008998 0.02492269 0.02400752 0.1181871

imputed.NW.traits.height.leafN = brm(cover.change ~ height.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                     family = gaussian(),
                                     prior = priors,
                                     data = imputed.NW)

summary(imputed.NW.traits.height.leafN)
# height x leafN

saveRDS(imputed.NW.traits.height.leafN, file = "./Results/all.imputed.traits.NW.height.leafN.rds")
bayes_R2(imputed.NW.traits.height.leafN)
# R2 0.04979303 0.02318114 0.01594568 0.1060524

conditional_effects(imputed.NW.traits.height.leafN)
# positive cover change with taller height and higher leafN

imputed.NW.traits.depth.SLA= brm(cover.change ~ root.depth.final*SLA.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                 family = gaussian(),
                                 prior = priors,
                                 data = imputed.NW)

summary(imputed.NW.traits.depth.SLA)
# nothing significant

imputed.NW.traits.depth.leafN = brm(cover.change ~ root.depth.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                    family = gaussian(),
                                    prior = priors,
                                    data = imputed.NW)

summary(imputed.NW.traits.depth.leafN)
# leafN

imputed.NW.traits.RTD.SRL = brm(cover.change ~ RTD.final*SRL.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                family = gaussian(),
                                prior = priors,
                                data = imputed.NW)

summary(imputed.NW.traits.RTD.SRL)
# nothing significant

imputed.NW.traits.SLA.RMF = brm(cover.change ~ SLA.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                family = gaussian(),
                                prior = priors,
                                data = imputed.NW)

summary(imputed.NW.traits.SLA.RMF)
# nothing significant

imputed.NW.traits.leafN.RMF = brm(cover.change ~ leafN.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                  family = gaussian(),
                                  prior = priors,
                                  data = imputed.NW)

summary(imputed.NW.traits.leafN.RMF)
# leaf significant positive
# MAP significant positive
# leafN:RMF significant positive

saveRDS(imputed.NW.traits.leafN.RMF, file = "./Results/all.imputed.NW.traits.leafN.RMF.rds")
bayes_R2(imputed.NW.traits.leafN.RMF)
# R2 0.05056307 0.02389004 0.01587269 0.1095646

conditional_effects(imputed.NW.traits.leafN.RMF)
# higher cover change with high leafN and high RMF or low leafN and low RMF


#### imputed ANNUAL traits lifespan model NW ####

piors <- c(prior(normal(0, 10), class = b))

annual.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                               SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                             family = gaussian(),
                             prior = priors,
                             data = imputed.NW.annual)

summary(annual.traits.NW.model)
# MAP
saveRDS(annual.traits.NW.model, file = "./Results/annual.imputed.traits.no_woody.rds")
bayes_R2(annual.traits.NW.model)
# R2 0.1542627 0.04448325 0.07546093 0.2476329

imputed.NW.annual.traits.height.leafN = brm(cover.change ~ height.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                            family = gaussian(),
                                            prior = priors,
                                            data = imputed.NW.annual)

summary(imputed.NW.annual.traits.height.leafN)
# MAP

imputed.NW.annual.traits.depth.SLA= brm(cover.change ~ root.depth.final*SLA.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                        family = gaussian(),
                                        prior = priors,
                                        data = imputed.NW.annual)

summary(imputed.NW.annual.traits.depth.SLA)
# MAP

imputed.NW.annual.traits.depth.leafN = brm(cover.change ~ root.depth.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                           family = gaussian(),
                                           prior = priors,
                                           data = imputed.NW.annual)

summary(imputed.NW.annual.traits.depth.leafN)
# MAP is positive, significant

imputed.NW.annual.traits.RTD.SRL = brm(cover.change ~ RTD.final*SRL.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                       family = gaussian(),
                                       prior = priors,
                                       data = imputed.NW.annual)

summary(imputed.NW.annual.traits.RTD.SRL)
# nothing significant

imputed.NW.annual.traits.SLA.RMF = brm(cover.change ~ SLA.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                       family = gaussian(),
                                       prior = priors,
                                       data = imputed.NW.annual)

summary(imputed.NW.annual.traits.SLA.RMF)
# MAP

imputed.NW.annual.traits.leafN.RMF = brm(cover.change ~ leafN.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                         family = gaussian(),
                                         prior = priors,
                                         data = imputed.NW.annual)

summary(imputed.NW.annual.traits.leafN.RMF)
# MAP significant positive
# leafN:RMF significant positive

saveRDS(imputed.NW.annual.traits.leafN.RMF, file = "./Results/annual.imputed.NW.traits.leafN.RMF.rds")
bayes_R2(imputed.NW.annual.traits.leafN.RMF)
# R2 0.1310284 0.04710103 0.0522819 0.2401864

conditional_effects(imputed.NW.annual.traits.leafN.RMF)
# higher cover change with high leafN and high RMF or low leafN and low RMF

#### imputed PERENNIAL traits lifespan model NW ####
perennial.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                  SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                family = gaussian(),
                                prior = priors,
                                data = imputed.NW.perennial)

summary(perennial.traits.NW.model)
# nothing significant

saveRDS(perennial.traits.NW.model, file = "./Results/perennial.imputed.traits.no_woody.rds")
bayes_R2(perennial.traits.NW.model)
# R2 0.081794 0.03372474 0.03279104 0.1602133

imputed.NW.perennial.traits.height.leafN = brm(cover.change ~ height.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                               family = gaussian(),
                                               prior = priors,
                                               data = imputed.NW.perennial)

summary(imputed.NW.perennial.traits.height.leafN)
# height x leafN significant positive

saveRDS(imputed.NW.perennial.traits.height.leafN, file = "./Results/perennial.imputed.traits.NW.height.leafN.rds")
bayes_R2(imputed.NW.perennial.traits.height.leafN)
#R2 0.06295824 0.03318674 0.01677983 0.1430209

conditional_effects(imputed.NW.perennial.traits.height.leafN)
# higher cover for taller plants with higher leafN

imputed.NW.perennial.traits.depth.SLA= brm(cover.change ~ root.depth.final*SLA.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                           family = gaussian(),
                                           prior = priors,
                                           data = imputed.NW.perennial)

summary(imputed.NW.perennial.traits.depth.SLA)
# nothing significant

imputed.NW.perennial.traits.depth.leafN = brm(cover.change ~ root.depth.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                              family = gaussian(),
                                              prior = priors,
                                              data = imputed.NW.perennial)

summary(imputed.NW.perennial.traits.depth.leafN)
# nothing significant

imputed.NW.perennial.traits.RTD.SRL = brm(cover.change ~ RTD.final*SRL.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                          family = gaussian(),
                                          prior = priors,
                                          data = imputed.NW.perennial)

summary(imputed.NW.perennial.traits.RTD.SRL)
# nothing significant

imputed.NW.perennial.traits.SLA.RMF = brm(cover.change ~ SLA.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                          family = gaussian(),
                                          prior = priors,
                                          data = imputed.NW.perennial)

summary(imputed.NW.perennial.traits.SLA.RMF)
# nothing significant

imputed.NW.perennial.traits.leafN.RMF = brm(cover.change ~ leafN.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                            family = gaussian(),
                                            prior = priors,
                                            data = imputed.NW.perennial)

summary(imputed.NW.perennial.traits.leafN.RMF)
# nothing significant

#### impute FORB traits functional group NW ####

priors <- c(prior(normal(0, 10), class = b))

forb.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                             SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                           family = gaussian(),
                           prior = priors,
                           data = imputed.NW.forb)

summary(forb.traits.NW.model)
# leafN significant

saveRDS(forb.traits.NW.model, file = "./Results/forb.imputed.traits.no_woody.rds")
bayes_R2(forb.traits.NW.model)
# R2 0.1148149 0.04304049 0.04828436 0.2129581

imputed.NW.forb.traits.height.leafN = brm(cover.change ~ height.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                          family = gaussian(),
                                          prior = priors,
                                          data = imputed.NW.forb)

summary(imputed.NW.forb.traits.height.leafN)
# leafN

imputed.NW.forb.traits.depth.SLA= brm(cover.change ~ root.depth.final*SLA.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                      family = gaussian(),
                                      prior = priors,
                                      data = imputed.NW.forb)

summary(imputed.NW.forb.traits.depth.SLA)
# nothing significant

imputed.NW.forb.traits.depth.leafN = brm(cover.change ~ root.depth.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                         family = gaussian(),
                                         prior = priors,
                                         data = imputed.NW.forb)

summary(imputed.NW.forb.traits.depth.leafN)
# leafN

imputed.NW.forb.traits.RTD.SRL = brm(cover.change ~ RTD.final*SRL.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                     family = gaussian(),
                                     prior = priors,
                                     data = imputed.NW.forb)

summary(imputed.NW.forb.traits.RTD.SRL)
# RTD x SRL significant, negative

saveRDS(imputed.NW.forb.traits.RTD.SRL, file = "./Results/forb.imputed.traits.NW.RTD.SRL.rds")
bayes_R2(imputed.NW.forb.traits.RTD.SRL)
# R2 0.0804646 0.04240653 0.02177166 0.1872374

conditional_effects(imputed.NW.forb.traits.RTD.SRL)
# higher cover with higher RTD and lower SRL or lower RTD with higher SRL

imputed.NW.forb.traits.SLA.RMF = brm(cover.change ~ SLA.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                     family = gaussian(),
                                     prior = priors,
                                     data = imputed.NW.forb)

summary(imputed.NW.forb.traits.SLA.RMF)
# nothing significant

imputed.NW.forb.traits.leafN.RMF = brm(cover.change ~ leafN.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                       family = gaussian(),
                                       prior = priors,
                                       data = imputed.NW.forb)

summary(imputed.NW.forb.traits.leafN.RMF)
# leafN significant positive

#### impute GRASS traits functional group NW ####
grass.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                              SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                            family = gaussian(),
                            prior = priors,
                            data = imputed.NW.grass)

summary(grass.traits.NW.model)
# nothing significant

saveRDS(grass.traits.NW.model, file = "./Results/grass.imputed.traits.no_woody.rds")
bayes_R2(grass.traits.NW.model)
# R2 0.1131223 0.04113157 0.04815195 0.2051698

imputed.NW.grass.traits.height.leafN = brm(cover.change ~ height.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                           family = gaussian(),
                                           prior = priors,
                                           data = imputed.NW.grass)

summary(imputed.NW.grass.traits.height.leafN)
# nothing significant

imputed.NW.grass.traits.depth.SLA= brm(cover.change ~ root.depth.final*SLA.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                       family = gaussian(),
                                       prior = priors,
                                       data = imputed.NW.grass)

summary(imputed.NW.grass.traits.depth.SLA)
# nothing significant

imputed.NW.grass.traits.depth.leafN = brm(cover.change ~ root.depth.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                          family = gaussian(),
                                          prior = priors,
                                          data = imputed.NW.grass)

summary(imputed.NW.grass.traits.depth.leafN)
# MAP is positive, significant

imputed.NW.grass.traits.RTD.SRL = brm(cover.change ~ RTD.final*SRL.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                      family = gaussian(),
                                      prior = priors,
                                      data = imputed.NW.grass)

summary(imputed.NW.grass.traits.RTD.SRL)
# nothing significant

imputed.NW.grass.traits.SLA.RMF = brm(cover.change ~ SLA.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                      family = gaussian(),
                                      prior = priors,
                                      data = imputed.NW.grass)

summary(imputed.NW.grass.traits.SLA.RMF)
# nothing significant

imputed.NW.grass.traits.leafN.RMF = brm(cover.change ~ leafN.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                        family = gaussian(),
                                        prior = priors,
                                        data = imputed.NW.grass)

summary(imputed.NW.grass.traits.leafN.RMF)
# nothing significant

#### imputed traits ANNUAL FORB NW ####

priors <- c(prior(normal(0, 10), class = b))

annual.forb.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                    SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                  family = gaussian(),
                                  prior = priors,
                                  data = imputed.NW.annual.forb)

summary(annual.forb.traits.NW.model)
# MAP significant

saveRDS(annual.forb.traits.NW.model, file = "./Results/annual.forb.imputed.traits.no_woody.rds")
bayes_R2(annual.forb.traits.NW.model)
# R2 0.2252887 0.0659826 0.1116487 0.3650369

imputed.NW.annual.forb.traits.height.leafN = brm(cover.change ~ height.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                 family = gaussian(),
                                                 prior = priors,
                                                 data = imputed.NW.annual.forb)

summary(imputed.NW.annual.forb.traits.height.leafN)
# nothing significant

imputed.NW.annual.forb.traits.depth.SLA= brm(cover.change ~ root.depth.final*SLA.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                             family = gaussian(),
                                             prior = priors,
                                             data = imputed.NW.annual.forb)

summary(imputed.NW.annual.forb.traits.depth.SLA)
# MAP

imputed.NW.annual.forb.traits.depth.leafN = brm(cover.change ~ root.depth.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                family = gaussian(),
                                                prior = priors,
                                                data = imputed.NW.annual.forb)

summary(imputed.NW.annual.forb.traits.depth.leafN)
# MAP is positive, significant

imputed.NW.annual.forb.traits.RTD.SRL = brm(cover.change ~ RTD.final*SRL.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                            family = gaussian(),
                                            prior = priors,
                                            data = imputed.NW.annual.forb)

summary(imputed.NW.annual.forb.traits.RTD.SRL)
# nothing significant

imputed.NW.annual.forb.traits.SLA.RMF = brm(cover.change ~ SLA.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                            family = gaussian(),
                                            prior = priors,
                                            data = imputed.NW.annual.forb)

summary(imputed.NW.annual.forb.traits.SLA.RMF)
# MAP

imputed.NW.annual.forb.traits.leafN.RMF = brm(cover.change ~ leafN.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                              family = gaussian(),
                                              prior = priors,
                                              data = imputed.NW.annual.forb)

summary(imputed.NW.annual.forb.traits.leafN.RMF)
# MAP significant positive

#### imputed traits PERENNIAL FORB NW ####

perennial.forb.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                       SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                     family = gaussian(),
                                     prior = priors,
                                     data = imputed.NW.perennial.forb)

summary(perennial.forb.traits.NW.model)
# nothing significant

saveRDS(perennial.forb.traits.NW.model, file = "./Results/perennial.forb.imputed.traits.no_woody.rds")
bayes_R2(perennial.forb.traits.NW.model)
#R2 0.1475479 0.05794628 0.06079698 0.2874462

imputed.NW.perennial.forb.traits.height.leafN = brm(cover.change ~ height.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                    family = gaussian(),
                                                    prior = priors,
                                                    data = imputed.NW.perennial.forb)

summary(imputed.NW.perennial.forb.traits.height.leafN)
# leafN

imputed.NW.perennial.forb.traits.depth.SLA= brm(cover.change ~ root.depth.final*SLA.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                family = gaussian(),
                                                prior = priors,
                                                data = imputed.NW.perennial.forb)

summary(imputed.NW.perennial.forb.traits.depth.SLA)
# nothing significant

imputed.NW.perennial.forb.traits.depth.leafN = brm(cover.change ~ root.depth.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                   family = gaussian(),
                                                   prior = priors,
                                                   data = imputed.NW.perennial.forb)

summary(imputed.NW.perennial.forb.traits.depth.leafN)
# leafN

imputed.NW.perennial.forb.traits.RTD.SRL = brm(cover.change ~ RTD.final*SRL.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                               family = gaussian(),
                                               prior = priors,
                                               data = imputed.NW.perennial.forb)

summary(imputed.NW.perennial.forb.traits.RTD.SRL)
# nothing significant

imputed.NW.perennial.forb.traits.SLA.RMF = brm(cover.change ~ SLA.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                               family = gaussian(),
                                               prior = priors,
                                               data = imputed.NW.perennial.forb)

summary(imputed.NW.perennial.forb.traits.SLA.RMF)
# nothing significant

imputed.NW.perennial.forb.traits.leafN.RMF = brm(cover.change ~ leafN.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                 family = gaussian(),
                                                 prior = priors,
                                                 data = imputed.NW.perennial.forb)

summary(imputed.NW.perennial.forb.traits.leafN.RMF)
# leafN:RMF significant positive





#### imputed traits PERENNIAL GRASS NW #####

perennial.grass.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                        SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                      family = gaussian(),
                                      prior = priors,
                                      data = imputed.NW.perennial.grass)

summary(perennial.grass.traits.NW.model)
# nothing significant

saveRDS(perennial.grass.traits.NW.model, file = "./Results/perennial.grass.imputed.traits.no_woody.rds")
bayes_R2(perennial.grass.traits.NW.model)
#R2 0.1499307 0.06219888 0.05914333 0.2792502

imputed.NW.perennial.grass.traits.height.leafN = brm(cover.change ~ height.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                     family = gaussian(),
                                                     prior = priors,
                                                     data = imputed.NW.perennial.grass)

summary(imputed.NW.perennial.grass.traits.height.leafN)
# nothing significant

imputed.NW.perennial.grass.traits.depth.SLA= brm(cover.change ~ root.depth.final*SLA.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                 family = gaussian(),
                                                 prior = priors,
                                                 data = imputed.NW.perennial.grass)

summary(imputed.NW.perennial.grass.traits.depth.SLA)
# nothing significant

imputed.NW.perennial.grass.traits.depth.leafN = brm(cover.change ~ root.depth.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                    family = gaussian(),
                                                    prior = priors,
                                                    data = imputed.NW.perennial.grass)

summary(imputed.NW.perennial.grass.traits.depth.leafN)
# nothing significant

imputed.NW.perennial.grass.traits.RTD.SRL = brm(cover.change ~ RTD.final*SRL.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                family = gaussian(),
                                                prior = priors,
                                                data = imputed.NW.perennial.grass)

summary(imputed.NW.perennial.grass.traits.RTD.SRL)
# nothing significant

imputed.NW.perennial.grass.traits.SLA.RMF = brm(cover.change ~ SLA.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                family = gaussian(),
                                                prior = priors,
                                                data = imputed.NW.perennial.grass)

summary(imputed.NW.perennial.grass.traits.SLA.RMF)
# nothing significant

imputed.NW.perennial.grass.traits.leafN.RMF = brm(cover.change ~ leafN.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                  family = gaussian(),
                                                  prior = priors,
                                                  data = imputed.NW.perennial.grass)

summary(imputed.NW.perennial.grass.traits.leafN.RMF)

#### ALL Backtransform ####

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)

# imputed traits with woody species
imputed.NW = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.outliersRM.csv", row.names = 1) %>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)

imputed.traits.NW.model = readRDS("./Results/all.imputed.traits.no_woody.rds")

imputed.NW$leafN.bt =(imputed.NW$leafN.final*7.114868) + 21.16957
imputed.NW$height.bt = (imputed.NW$height.final*0.2353906) + 0.3476824
imputed.NW$rootN.bt = (imputed.NW$rootN.final*4.488773) + 9.855027
imputed.NW$SLA.bt = (imputed.NW$SLA.final*8.581554) + 19.93883
imputed.NW$Depth.bt = (imputed.NW$root.depth.final*0.5333377) + 0.5784653
imputed.NW$Diam.bt = (imputed.NW$rootDiam.final*0.1598629) + 0.3501538
imputed.NW$SRL.bt = (imputed.NW$SRL.final*74.26099) + 106.4204
imputed.NW$RTD.bt = (imputed.NW$RTD.final*0.1159369) + 0.2349392
imputed.NW$RMF.bt = (imputed.NW$RMF.final*0.1058164) + 0.4007788
imputed.NW$DSI.bt = (imputed.NW$mean.DSI*0.1950438) + -0.4604709
imputed.NW$MAP.bt = (imputed.NW$mean.MAP*487.1451) + 691.6107

imputed.NW.2 = imputed.NW %>%
  dplyr::select(cover.change,leafN.bt:MAP.bt) %>%
  drop_na()

range(imputed.NW.2$cover.change) # -20.66667  21.24500
range(imputed.NW.2$leafN.bt) # 2.68952 44.19000
range(imputed.NW.2$height.bt) # 0.009000065 1.322333168
range(imputed.NW.2$rootN.bt) # 0.3800001 25.6283267
range(imputed.NW.2$SLA.bt) # 2.617597 46.320000
range(imputed.NW.2$Depth.bt) # 0.048000 3.709334
range(imputed.NW.2$Diam.bt) # 0.07840003 1.12601234
range(imputed.NW.2$SRL.bt) # 2.497204 417.749987
range(imputed.NW.2$RTD.bt) # 0.02325659 0.61999992
range(imputed.NW.2$RMF.bt) # 0.1211194 0.7328658
range(imputed.NW.2$DSI.bt) # -0.8722616  0.1242845
range(imputed.NW.2$MAP.bt) # 132.8 2366.0

all.effects = conditional_effects(imputed.traits.NW.model)

#### ALL leafN ####

leafN.effects = all.effects$leafN.final
leafN.effects$leafN.bt = (leafN.effects$leafN.final*7.114868) + 21.16957

x.value = c(-2,-1,0,1,2,3)
(x.value*7.114868) + 21.16957

all.leafN.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = leafN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.bt, y = estimate__), color = "#769370", size = 1.5) +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Leaf N", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.68951,44.19001)+
  theme_classic()
all.leafN.plot

ggsave("./Plots/all.traits.NW.leafN.pdf", height = 3, width = 3)

#scale_x_continuous(limits = c(-2.6,3.3),
                   #breaks = seq(-2, 3, by = 1),
                   #labels = c(6.9, 14.1, 21.2, 28.3, 35.4, 42.5))+

#### All height ####

height.effects = all.effects$height.final
height.effects$height.bt = (height.effects$height.final*0.2353906) + 0.3476824

x.value = c(-1,0,1,2,3,4)
(x.value*0.2353906) + 0.3476824

all.height.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = height.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.bt, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = height.effects, aes(x = height.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Height", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.009000064,1.322333169)+
  theme_classic()
all.height.plot

ggsave("./Plots/all.traits.NW.height.pdf", height = 3, width = 3)

#### All rootN ####

rootN.effects = all.effects$rootN.final
rootN.effects$rootN.bt = (rootN.effects$rootN.final*4.488773) + 9.855027

x.value = c(-2,-1,0,1,2,3)
(x.value*4.488773) + 9.855027

all.rootN.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = rootN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.bt, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Root N", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.3800000,25.6283268)+
  theme_classic()
all.rootN.plot

ggsave("./Plots/all.traits.NW.rootN.pdf", height = 3, width = 3)

#### All SLA ####

SLA.effects = all.effects$SLA.final
SLA.effects$SLA.bt = (SLA.effects$SLA.final*8.581554) + 19.93883

x.value = c(-2,-1,0,1,2,3)
(x.value*8.581554) + 19.93883

all.SLA.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = SLA.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.bt, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Specific Leaf Area", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.617596, 46.320001)+
  theme_classic()
all.SLA.plot

ggsave("./Plots/all.traits.NW.SLA.pdf", height = 3, width = 3)

#### All Depth ####

root.depth.effects = all.effects$root.depth.final
root.depth.effects$Depth.bt = (root.depth.effects$root.depth.final*0.5333377) + 0.5784653

x.value = c(0,1,2,3,4,5)
(x.value*0.5333377) + 0.5784653

all.root.depth.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = Depth.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = Depth.bt, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = Depth.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Rooting Depth", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.0479999,3.709335)+
  theme_classic()
all.root.depth.plot

ggsave("./Plots/all.traits.NW.root.depth.pdf", height = 3, width = 3)

#### All Diameter ####

rootDiam.effects = all.effects$rootDiam.final
rootDiam.effects$Diam.bt = (rootDiam.effects$rootDiam.final*0.1598629) + 0.3501538

x.value = c(-1,0,1,2,3,4)
(x.value*0.1598629) + 0.3501538

all.rootDiam.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = Diam.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = Diam.bt, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = Diam.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Root Diameter", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.07840002,1.12601235)+
  theme_classic()
all.rootDiam.plot

ggsave("./Plots/all.traits.NW.rootDiam.pdf", height = 3, width = 3)

#### All SRL ####

SRL.effects = all.effects$SRL.final
SRL.effects$SRL.bt = (SRL.effects$SRL.final*74.26099) + 106.4204

x.value = c(-1,0,1,2,3,4)
(x.value*74.26099) + 106.4204

all.SRL.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = SRL.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.bt, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Specific Root Length", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.497203,417.749988)+
  theme_classic()
all.SRL.plot

ggsave("./Plots/all.traits.NW.SRL.pdf", height = 3, width = 3)

#### All RTD ####

RTD.effects = all.effects$RTD.final
RTD.effects$RTD.bt = (RTD.effects$RTD.final*0.1159369) + 0.2349392

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1159369) + 0.2349392

all.RTD.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = RTD.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.bt, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Root Tissue Density", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.02325658, 0.61999993)+
  theme_classic()
all.RTD.plot

ggsave("./Plots/all.traits.NW.RTD.pdf", height = 3, width = 3)

#### All RMF ####

RMF.effects = all.effects$RMF.final
RMF.effects$RMF.bt = (RMF.effects$RMF.final*0.1058164) + 0.4007788

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1058164) + 0.4007788

all.RMF.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = RMF.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RMF.effects, aes(x = RMF.bt, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RMF.effects, aes(x = RMF.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Root Mass Fraction", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.1211193,0.7328659)+
  theme_classic()
all.RMF.plot

ggsave("./Plots/all.traits.NW.RMF.pdf", height = 3, width = 3)

#### All DSI ####

DSI.effects = all.effects$mean.DSI
DSI.effects$DSI.bt = (DSI.effects$mean.DSI*0.1950438) + -0.4604709

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1950438) + -0.4604709

all.DSI.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = DSI.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = DSI.effects, aes(x = DSI.bt, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = DSI.effects, aes(x = DSI.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Drought Severity Index", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(-0.8722616,0.1242846)+
   theme_classic()
all.DSI.plot

ggsave("./Plots/all.traits.NW.DSI.pdf", height = 3, width = 3)

#### All MAP ####

MAP.effects = all.effects$mean.MAP
MAP.effects$MAP.bt= (MAP.effects$mean.MAP*487.1451) + 691.6107

x.value = c(-1,0,1,2,3)
(x.value*487.1451) + 691.6107

all.MAP.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = MAP.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = MAP.effects, aes(x = MAP.bt, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = MAP.effects, aes(x = MAP.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Mean Annual Precipitation", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(132.7,2366.1)+
  theme_classic()
all.MAP.plot

ggsave("./Plots/all.traits.NW.MAP.pdf", height = 3, width = 3)

#### ALL interactions ####

imputed.NW.traits.height.leafN = readRDS("./Results/all.imputed.traits.NW.height.leafN.rds")

height.leafN.effect = conditional_effects(imputed.NW.traits.height.leafN, effects = "height.final:leafN.final")$`height.final:leafN.final`
height.leafN.effect$height.bt = (height.leafN.effect$height.final*0.2353906) + 0.3476824

# only plot high and low values

height.leafN.effect.2 = height.leafN.effect %>%
  filter(effect2__ %in% c(-0.98,1.03))

# height
x.value = c(-1,0,1,2,3,4)
(x.value*0.2353906) + 0.3476824

# leafN
x.value = c(-0.98,1.03)
(x.value*7.114868) + 21.16957


all.height.x.leafN.plot = ggplot(data = height.leafN.effect.2, aes(x = height.bt, y = estimate__, group = as.factor(effect2__))) +
  geom_line(aes(color = as.factor(effect2__)), size = 1.5,show.legend = FALSE) +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = as.factor(effect2__)), alpha = 0.5,show.legend = FALSE) +  
  geom_line(aes(y = upper__, color = factor(effect2__)), size = 1.2,show.legend = FALSE) +
  geom_line(aes(y = lower__, color = factor(effect2__)), size = 1.2,show.legend = FALSE) + 
  labs(x = "Height", y = "Percent Cover Change", color = "Leaf N") +
  scale_colour_manual(values = c("black", "#769370"), labels = c("28.5","14.2"))+
  scale_fill_manual(values = c("black", "#769370"))+
  theme_classic()
all.height.x.leafN.plot

ggsave("./Plots/all.traits.NW.height.leafN.pdf", height = 3, width = 3)
ggsave("./Plots/all.traits.NW.height.leafN.legend.pdf", height = 3, width = 3)

imputed.NW.traits.leafN.RMF = readRDS("./Results/all.imputed.NW.traits.leafN.RMF.rds")

leafN.RMF.effect = conditional_effects(imputed.NW.traits.leafN.RMF, effects = "leafN.final:RMF.final")$`leafN.final:RMF.final`
leafN.RMF.effect$leafN.bt = (leafN.RMF.effect$leafN.final*7.114868) + 21.16957

# only plot high and low values

leafN.RMF.effect.2 = leafN.RMF.effect %>%
  filter(effect2__ %in% c(-1.02,0.99))
# leafN
x.value = c(-2,-1,0,1,2,3)
(x.value*7.114868) + 21.16957

# RMF
x.value = c(-1.02,0.99)
(x.value*0.1058164) + 0.4007788

all.leafN.x.RMF.plot = ggplot(data = leafN.RMF.effect.2, aes(x = leafN.bt, y = estimate__, group = as.factor(effect2__))) +
  geom_line(aes(color = as.factor(effect2__)), size = 1.5,show.legend = FALSE) +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = as.factor(effect2__)), alpha = 0.5,show.legend = FALSE) +  
  geom_line(aes(y = upper__, color = factor(effect2__)), size = 1.2,show.legend = FALSE) +
  geom_line(aes(y = lower__, color = factor(effect2__)), size = 1.2,show.legend = FALSE) + 
  labs(x = "Leaf N", y = "Percent Cover Change", color = "RMF") +
  scale_colour_manual(values = c("black", "#769370"), labels = c("0.5","0.3"))+
  scale_fill_manual(values = c("black", "#769370"))+
  xlim(2.68951,44.19001)+
  theme_classic()
all.leafN.x.RMF.plot

ggsave("./Plots/all.traits.NW.leafN.RMF.pdf", height = 3, width = 3)
ggsave("./Plots/all.traits.NW.leafN.RMF.legend.pdf", height = 3, width = 3)

#### ANNUAL Backtransform ####

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)

# imputed traits with woody species
imputed.NW.annual = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.annuals.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)

annual.traits.NW.model = readRDS("./Results/annual.imputed.traits.no_woody.rds")

imputed.NW.annual$leafN.bt =(imputed.NW.annual$leafN.final*7.384253) + 22.54886
imputed.NW.annual$height.bt = (imputed.NW.annual$height.final*0.2308751) + 0.3156657
imputed.NW.annual$rootN.bt = (imputed.NW.annual$rootN.final*4.33202) + 9.997142
imputed.NW.annual$SLA.bt = (imputed.NW.annual$SLA.final*9.312314) + 22.81754
imputed.NW.annual$Depth.bt = (imputed.NW.annual$root.depth.final*0.4346265) + 0.5402869
imputed.NW.annual$Diam.bt = (imputed.NW.annual$rootDiam.final*0.1367222) + 0.3407223
imputed.NW.annual$SRL.bt = (imputed.NW.annual$SRL.final*73.49115) + 125.697
imputed.NW.annual$RTD.bt = (imputed.NW.annual$RTD.final*0.1111574) + 0.2026466
imputed.NW.annual$RMF.bt = (imputed.NW.annual$RMF.final*0.1037188) + 0.3934184
imputed.NW.annual$DSI.bt = (imputed.NW.annual$mean.DSI*0.2300417) + -0.4690337
imputed.NW.annual$MAP.bt = (imputed.NW.annual$mean.MAP*293.1676) + 477.7734

annual.imputed.NW.2 = imputed.NW.annual %>%
  dplyr::select(cover.change,leafN.bt:MAP.bt) %>%
  drop_na()

annual.effects = conditional_effects(annual.traits.NW.model)

#### Annual leafN ####
leafN.effects = annual.effects$leafN.final
leafN.effects$leafN.bt= (leafN.effects$leafN.final*7.384253) + 22.54886

attr(annual.imputed.NW.2$leafN.final, "scaled:scale")
attr(annual.imputed.NW.2$leafN.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*7.384253) + 22.54886

annual.leafN.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = leafN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Leaf N", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.68951,44.19001)+
  theme_classic()
annual.leafN.plot

ggsave("./Plots/annual.traits.NW.leafN.pdf", height = 3, width = 3)

#### Annual height ####

height.effects = annual.effects$height.final
height.effects$height.bt= (height.effects$height.final*0.2308751) + 0.3156657

attr(annual.imputed.NW.2$height.final, "scaled:scale")
attr(annual.imputed.NW.2$height.final, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*0.2308751) + 0.3156657

annual.height.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = height.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = height.effects, aes(x = height.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Height", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.009000064,1.322333169)+
  theme_classic()
annual.height.plot

ggsave("./Plots/annual.traits.NW.height.pdf", height = 3, width = 3)


#### Annual RootN ####
rootN.effects = annual.effects$rootN.final
rootN.effects$rootN.bt= (rootN.effects$rootN.final*4.33202) + 9.997142

attr(annual.imputed.NW.2$rootN.final, "scaled:scale")
attr(annual.imputed.NW.2$rootN.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*4.33202) + 9.997142

annual.rootN.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = rootN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Root N", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.3800000,25.6283268)+
  theme_classic()
annual.rootN.plot

ggsave("./Plots/annual.traits.NW.rootN.pdf", height = 3, width = 3)

#### Annual SLA ####

SLA.effects = annual.effects$SLA.final
SLA.effects$SLA.bt= (SLA.effects$SLA.final*9.312314) + 22.81754

attr(annual.imputed.NW.2$SLA.final, "scaled:scale")
attr(annual.imputed.NW.2$SLA.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*9.312314) + 22.81754

annual.SLA.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = SLA.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Specific Leaf Area", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.617596,46.320001)+
  theme_classic()
annual.SLA.plot

ggsave("./Plots/annual.traits.NW.SLA.pdf", height = 3, width = 3)

#### Annual depth ####

root.depth.effects = annual.effects$root.depth.final
root.depth.effects$Depth.bt= (root.depth.effects$root.depth.final*0.4346265) + 0.5402869

attr(annual.imputed.NW.2$root.depth.final, "scaled:scale")
attr(annual.imputed.NW.2$root.depth.final, "scaled:center")

x.value = c(0,1,2,3,4,5)
(x.value*0.4346265) + 0.5402869

annual.root.depth.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = Depth.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = Depth.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = Depth.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Rooting Depth", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.0479999,3.709335)+
  theme_classic()
annual.root.depth.plot

ggsave("./Plots/annual.traits.NW.root.depth.pdf", height = 3, width = 3)

#### Annual Diam ####

rootDiam.effects = annual.effects$rootDiam.final
rootDiam.effects$Diam.bt= (rootDiam.effects$rootDiam.final*0.1367222) + 0.3407223

attr(annual.imputed.NW.2$rootDiam.final, "scaled:scale")
attr(annual.imputed.NW.2$rootDiam.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1367222) + 0.3407223

annual.rootDiam.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = Diam.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = Diam.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = Diam.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Root Diameter", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.07840002,1.12601235)+
  theme_classic()
annual.rootDiam.plot

ggsave("./Plots/annual.traits.NW.rootDiam.pdf", height = 3, width = 3)

#### Annual SRL ####

SRL.effects = annual.effects$SRL.final
SRL.effects$SRL.bt= (SRL.effects$SRL.final*73.49115) + 125.697

attr(annual.imputed.NW.2$SRL.final, "scaled:scale")
attr(annual.imputed.NW.2$SRL.final, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*73.49115) + 125.697

annual.SRL.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = SRL.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Specific Root Length", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.497203,417.749988)+
  theme_classic()
annual.SRL.plot

ggsave("./Plots/annual.traits.NW.SRL.pdf", height = 3, width = 3)

#### Annual RTD ####

RTD.effects = annual.effects$RTD.final
RTD.effects$RTD.bt= (RTD.effects$RTD.final*0.1111574) + 0.2026466

attr(annual.imputed.NW.2$RTD.final, "scaled:scale")
attr(annual.imputed.NW.2$RTD.final, "scaled:center")

x.value = c(-1,0,1,2,3)
(x.value*0.1111574) + 0.2026466

annual.RTD.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = RTD.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Root Tissue Density", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.02325658,0.61999993)+
  theme_classic()
annual.RTD.plot

ggsave("./Plots/annual.traits.NW.RTD.pdf", height = 3, width = 3)

#### Annual RMF ####

RMF.effects = annual.effects$RMF.final
RMF.effects$RMF.bt= (RMF.effects$RMF.final*0.1037188) + 0.3934184

attr(annual.imputed.NW.2$RMF.final, "scaled:scale")
attr(annual.imputed.NW.2$RMF.final, "scaled:center")

x.value = c(-3,-2,-1,0,1,2,3)
(x.value*0.1037188) + 0.3934184

annual.RMF.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = RMF.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RMF.effects, aes(x = RMF.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RMF.effects, aes(x = RMF.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Root Mass Fraction", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.1211193,0.7328659)+
  theme_classic()
annual.RMF.plot

ggsave("./Plots/annual.traits.NW.RMF.pdf", height = 3, width = 3)

#### Annual DSI ####

DSI.effects = annual.effects$mean.DSI
DSI.effects$DSI.bt= (DSI.effects$mean.DSI*0.2300417) + -0.4690337

attr(annual.imputed.NW.2$mean.DSI, "scaled:scale")
attr(annual.imputed.NW.2$mean.DSI, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.2300417) + -0.4690337

annual.DSI.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = DSI.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = DSI.effects, aes(x = DSI.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = DSI.effects, aes(x = DSI.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Drought Severity Index", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(-0.8722616,0.1242846)+
  theme_classic()
annual.DSI.plot

ggsave("./Plots/annual.traits.NW.DSI.pdf", height = 3, width = 3)

#### Annual MAP ####

MAP.effects = annual.effects$mean.MAP
MAP.effects$MAP.bt= (MAP.effects$mean.MAP*293.1676) + 477.7734

attr(annual.imputed.NW.2$mean.MAP, "scaled:scale")
attr(annual.imputed.NW.2$mean.MAP, "scaled:center")

x.value = c(-1,0,1,2,3)
(x.value*293.1676) + 477.7734

annual.MAP.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = MAP.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = MAP.effects, aes(x = MAP.bt, y = estimate__), color = "#979461", size = 1.5) +  
  geom_ribbon(data = MAP.effects, aes(x = MAP.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Mean Annual Precipitation", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(132.7,2366.1)+
  theme_classic()
annual.MAP.plot

ggsave("./Plots/annual.traits.NW.MAP.pdf", height = 3, width = 3)

#### ANNUAL interactions ####

annual.imputed.NW.traits.leafN.RMF = readRDS("./Results/annual.imputed.NW.traits.leafN.RMF.rds")

leafN.RMF.effect = conditional_effects(annual.imputed.NW.traits.leafN.RMF, effects = "leafN.final:RMF.final")$`leafN.final:RMF.final`
leafN.RMF.effect$leafN.bit= (leafN.RMF.effect$leafN.final*7.384253) + 22.54886

# only plot high and low values

leafN.RMF.effect.2 = leafN.RMF.effect %>%
  filter(effect2__ %in% c(-1.05,0.97))

# leafN
x.value = c(-2,-1,0,1,2,3)
(x.value*7.384253) + 22.54886

# RMF
x.value = c(-1.05,0.97)
(x.value*0.1037188) + 0.3934184

annual.leafN.x.RMF.plot = ggplot(data = leafN.RMF.effect.2, aes(x = leafN.bit, y = estimate__, group = as.factor(effect2__))) +
  geom_line(aes(color = as.factor(effect2__)), size = 1.5,show.legend = FALSE) +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = as.factor(effect2__)), alpha = 0.5,show.legend = FALSE) +  
  geom_line(aes(y = upper__, color = factor(effect2__)), size = 1.2,show.legend = FALSE) +
  geom_line(aes(y = lower__, color = factor(effect2__)), size = 1.2,show.legend = FALSE) + 
  labs(x = "Leaf N", y = "Percent Cover Change", color = "RMF") +
  scale_colour_manual(values = c("black", "#979461"), labels = c("0.49","0.28"))+
  scale_fill_manual(values = c("black", "#979461"))+
  xlim(2.68951,44.19001)+
  theme_classic()
annual.leafN.x.RMF.plot

ggsave("./Plots/annual.traits.NW.leafN.RMF.pdf", height = 3, width = 3)
ggsave("./Plots/annual.traits.NW.leafN.RMF.legend.pdf", height = 3, width = 3)

#### PERENNIAL Backtransform ####

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)
imputed.NW.perennial = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennials.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)

perennial.traits.NW.model = readRDS("./Results/perennial.imputed.traits.no_woody.rds")

imputed.NW.perennial$leafN.bt =(imputed.NW.perennial$leafN.final*6.812995) + 20.47022
imputed.NW.perennial$height.bt = (imputed.NW.perennial$height.final*0.2365653) + 0.3627378
imputed.NW.perennial$rootN.bt = (imputed.NW.perennial$rootN.final*4.532848) + 9.765547
imputed.NW.perennial$SLA.bt = (imputed.NW.perennial$SLA.final*7.966249) + 18.75888
imputed.NW.perennial$Depth.bt = (imputed.NW.perennial$root.depth.final*0.5706892) + 0.5903092
imputed.NW.perennial$Diam.bt = (imputed.NW.perennial$rootDiam.final*0.166311) + 0.3535401
imputed.NW.perennial$SRL.bt = (imputed.NW.perennial$SRL.final*73.08878) + 98.7131
imputed.NW.perennial$RTD.bt = (imputed.NW.perennial$RTD.final*0.1151264) + 0.249861
imputed.NW.perennial$RMF.bt = (imputed.NW.perennial$RMF.final*0.1053099) + 0.4042604
imputed.NW.perennial$DSI.bt = (imputed.NW.perennial$mean.DSI*0.1779261) + -0.4547629
imputed.NW.perennial$MAP.bt = (imputed.NW.perennial$mean.MAP*522.5359) + 772.0729

perennial.imputed.NW.2 = imputed.NW.perennial %>%
  dplyr::select(cover.change,leafN.bt:MAP.bt) %>%
  drop_na()

perennial.effects = conditional_effects(perennial.traits.NW.model)

#### Perennial LeafN ####

leafN.effects = perennial.effects$leafN.final
leafN.effects$leafN.bt= (leafN.effects$leafN.final*6.812995) + 20.47022

attr(perennial.imputed.NW.2$leafN.final, "scaled:scale")
attr(perennial.imputed.NW.2$leafN.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*6.812995) + 20.47022

perennial.leafN.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = leafN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Leaf N", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.68951,44.19001)+
  theme_classic()
perennial.leafN.plot

ggsave("./Plots/perennial.traits.NW.leafN.pdf", height = 3, width = 3)

#### Perennial Height ####

height.effects = perennial.effects$height.final
height.effects$height.bt= (height.effects$height.final*0.2365653) + 0.3627378

attr(perennial.imputed.NW.2$height.final, "scaled:scale")
attr(perennial.imputed.NW.2$height.final, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*0.2365653) + 0.3627378

perennial.height.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = height.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = height.effects, aes(x = height.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Height", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.009000064,1.322333169)+
  theme_classic()
perennial.height.plot

ggsave("./Plots/perennial.traits.NW.height.pdf", height = 3, width = 3)

#### Perennial RootN ####

rootN.effects = perennial.effects$rootN.final
rootN.effects$rootN.bt= (rootN.effects$rootN.final*4.532848) + 9.765547

attr(perennial.imputed.NW.2$rootN.final, "scaled:scale")
attr(perennial.imputed.NW.2$rootN.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*4.532848) + 9.765547

perennial.rootN.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = rootN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Root N", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.3800000,25.6283268)+
  theme_classic()
perennial.rootN.plot

ggsave("./Plots/perennial.traits.NW.rootN.pdf", height = 3, width = 3)

#### Perennial SLA ####

SLA.effects = perennial.effects$SLA.final
SLA.effects$SLA.bt= (SLA.effects$SLA.final*7.966249) + 18.75888

attr(perennial.imputed.NW.2$SLA.final, "scaled:scale")
attr(perennial.imputed.NW.2$SLA.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*7.966249) + 18.75888

perennial.SLA.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = SLA.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Specific Leaf Area", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.617596,46.320001)+
  theme_classic()
perennial.SLA.plot

ggsave("./Plots/perennial.traits.NW.SLA.pdf", height = 3, width = 3)

#### Perennial Depth ####

root.depth.effects = perennial.effects$root.depth.final
root.depth.effects$Depth.bt= (root.depth.effects$root.depth.final*0.5706892) + 0.5903092

attr(perennial.imputed.NW.2$root.depth.final, "scaled:scale")
attr(perennial.imputed.NW.2$root.depth.final, "scaled:center")

x.value = c(-1,0,1,2,3,4,5)
(x.value*0.5706892) + 0.5903092

perennial.root.depth.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = Depth.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = Depth.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = Depth.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Rooting Depth", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.0479999,3.709335)+
  theme_classic()
perennial.root.depth.plot

ggsave("./Plots/perennial.traits.NW.root.depth.pdf", height = 3, width = 3)

#### Perennial Diam ####

rootDiam.effects = perennial.effects$rootDiam.final
rootDiam.effects$Diam.bt= (rootDiam.effects$rootDiam.final*0.166311) + 0.3535401

attr(perennial.imputed.NW.2$rootDiam.final, "scaled:scale")
attr(perennial.imputed.NW.2$rootDiam.final, "scaled:center")

x.value = c(-1,0,1,2,3,4,5)
(x.value*0.166311) + 0.3535401

perennial.rootDiam.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = Diam.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = Diam.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = Diam.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Root Diameter", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.07840002,1.12601235)+
  theme_classic()
perennial.rootDiam.plot

ggsave("./Plots/perennial.traits.NW.rootDiam.pdf", height = 3, width = 3)

#### Perennial SRL ####

SRL.effects = perennial.effects$SRL.final
SRL.effects$SRL.bt= (SRL.effects$SRL.final*73.08878) + 98.7131

attr(perennial.imputed.NW.2$SRL.final, "scaled:scale")
attr(perennial.imputed.NW.2$SRL.final, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*73.08878) + 98.7131

perennial.SRL.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = SRL.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Specific Root Length", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.497203,417.749988)+
  theme_classic()
perennial.SRL.plot

ggsave("./Plots/perennial.traits.NW.SRL.pdf", height = 3, width = 3)

#### Perennial RTD ####

RTD.effects = perennial.effects$RTD.final
RTD.effects$RTD.bt= (RTD.effects$RTD.final*0.1151264) + 0.249861

attr(perennial.imputed.NW.2$RTD.final, "scaled:scale")
attr(perennial.imputed.NW.2$RTD.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1151264) + 0.249861

perennial.RTD.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = RTD.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Root Tissue Density", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.02325658,0.61999993)+
  theme_classic()
perennial.RTD.plot

ggsave("./Plots/perennial.traits.NW.RTD.pdf", height = 3, width = 3)

#### Perennial RMF ####

RMF.effects = perennial.effects$RMF.final
RMF.effects$RMF.bt= (RMF.effects$RMF.final*0.1053099) + 0.4042604

attr(perennial.imputed.NW.2$RMF.final, "scaled:scale")
attr(perennial.imputed.NW.2$RMF.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1053099) + 0.4042604

perennial.RMF.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = RMF.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RMF.effects, aes(x = RMF.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RMF.effects, aes(x = RMF.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Root Mass Fraction", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.1211193,0.7328659)+
  theme_classic()
perennial.RMF.plot

ggsave("./Plots/perennial.traits.NW.RMF.pdf", height = 3, width = 3)

#### Perennial DSI ####

DSI.effects = perennial.effects$mean.DSI
DSI.effects$DSI.bt= (DSI.effects$mean.DSI*0.1779261) + -0.4547629

attr(perennial.imputed.NW.2$mean.DSI, "scaled:scale")
attr(perennial.imputed.NW.2$mean.DSI, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1779261) + -0.4547629

perennial.DSI.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = DSI.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = DSI.effects, aes(x = DSI.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = DSI.effects, aes(x = DSI.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Drought Severity Index", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(-0.8722616,0.1242846)+
  theme_classic()
perennial.DSI.plot

ggsave("./Plots/perennial.traits.NW.DSI.pdf", height = 3, width = 3)

#### Perennial MAP ####

MAP.effects = perennial.effects$mean.MAP
MAP.effects$MAP.bt= (MAP.effects$mean.MAP*522.5359) + 772.0729

attr(perennial.imputed.NW.2$mean.MAP, "scaled:scale")
attr(perennial.imputed.NW.2$mean.MAP, "scaled:center")

x.value = c(-1,0,1,2,3)
(x.value*522.5359) + 772.0729

perennial.MAP.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = MAP.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = MAP.effects, aes(x = MAP.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = MAP.effects, aes(x = MAP.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Mean Annual Precipitation", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(132.7,2366.1)+
  theme_classic()
perennial.MAP.plot

ggsave("./Plots/perennial.traits.NW.MAP.pdf", height = 3, width = 3)

#### PERENNIAL interactions ####

perennial.imputed.NW.traits.height.leafN = readRDS("./Results/perennial.imputed.traits.NW.height.leafN.rds")

height.leafN.effect = conditional_effects(perennial.imputed.NW.traits.height.leafN, effects = "height.final:leafN.final")$`height.final:leafN.final`
height.leafN.effect$height.bt= (height.leafN.effect$height.final*0.2365653) + 0.3627378

# only plot high and low values

height.leafN.effect.2 = height.leafN.effect %>%
  filter(effect2__ %in% c(-0.97,1.02))

# height
x.value = c(-1,0,1,2,3,4)
(x.value*0.2365653) + 0.3627378

# leafN
x.value = c(-0.97,1.02)
(x.value*6.812995) + 20.47022

perennial.height.x.leafN.plot = ggplot(data = height.leafN.effect.2, aes(x = height.bt, y = estimate__, group = as.factor(effect2__))) +
  geom_line(aes(color = as.factor(effect2__)), size = 1.5,show.legend = FALSE) +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = as.factor(effect2__)), alpha = 0.5,show.legend = FALSE) +  
  geom_line(aes(y = upper__, color = factor(effect2__)), size = 1.2,show.legend = FALSE) +
  geom_line(aes(y = lower__, color = factor(effect2__)), size = 1.2,show.legend = FALSE) + 
  labs(x = "Height", y = "Percent Cover Change", color = "Leaf N") +
  scale_colour_manual(values = c("black", "#F1C646"), labels = c("27.41","13.86"))+
  scale_fill_manual(values = c("black", "#F1C646"))+
  xlim(0.009000064,1.322333169)+
  theme_classic()
perennial.height.x.leafN.plot

ggsave("./Plots/perennial.traits.NW.height.leafN.pdf", height = 3, width = 3)
ggsave("./Plots/perennial.traits.NW.height.leafN.legend.pdf", height = 3, width = 3)

#### FORB Backtransform ####

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)
imputed.NW.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.forbs.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)

forb.traits.NW.model = readRDS("./Results/forb.imputed.traits.no_woody.rds")

imputed.NW.forb$leafN.bt =(imputed.NW.forb$leafN.final*6.403163) + 21.66854
imputed.NW.forb$height.bt = (imputed.NW.forb$height.final*0.2396416) + 0.3214792
imputed.NW.forb$rootN.bt = (imputed.NW.forb$rootN.final*4.173927) + 10.5821
imputed.NW.forb$SLA.bt = (imputed.NW.forb$SLA.final*8.362926) + 20.14641
imputed.NW.forb$Depth.bt = (imputed.NW.forb$root.depth.final*0.6095345) + 0.6503508
imputed.NW.forb$Diam.bt = (imputed.NW.forb$rootDiam.final*0.1590925) + 0.3579696
imputed.NW.forb$SRL.bt = (imputed.NW.forb$SRL.final*71.25811) + 105.5323
imputed.NW.forb$RTD.bt = (imputed.NW.forb$RTD.final*0.1184977) + 0.227384
imputed.NW.forb$RMF.bt = (imputed.NW.forb$RMF.final*0.1176732) + 0.3992223
imputed.NW.forb$DSI.bt = (imputed.NW.forb$mean.DSI*0.1949549) + -0.4644806
imputed.NW.forb$MAP.bt = (imputed.NW.forb$mean.MAP*448.2603) + 657.7821

forb.imputed.NW.2 = imputed.NW.forb %>%
  dplyr::select(cover.change,leafN.bt:MAP.bt) %>%
  drop_na()

forb.effects = conditional_effects(forb.traits.NW.model)

#### Forb leafN ####

leafN.effects = forb.effects$leafN.final
leafN.effects$leafN.bt= (leafN.effects$leafN.final*6.403163) + 21.66854

attr(forb.imputed.NW.2$leafN.final, "scaled:scale")
attr(forb.imputed.NW.2$leafN.final, "scaled:center")

x.value = c(-3,-2,-1,0,1,2,3)
(x.value*6.403163) + 21.66854

forb.leafN.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = leafN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.bt, y = estimate__), color = "#F17236", size = 1.5) +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Leaf N", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.68951,44.19001)+
  theme_classic()
forb.leafN.plot

ggsave("./Plots/forb.traits.NW.leafN.pdf", height = 3, width = 3)

#### Forb height ####

height.effects = forb.effects$height.final
height.effects$height.bt= (height.effects$height.final*0.2396416) + 0.3214792

attr(forb.imputed.NW.2$height.final, "scaled:scale")
attr(forb.imputed.NW.2$height.final, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*0.2396416) + 0.3214792

forb.height.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = height.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = height.effects, aes(x = height.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Height", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.009000064,1.322333169)+
  theme_classic()
forb.height.plot

ggsave("./Plots/forb.traits.NW.height.pdf", height = 3, width = 3)

#### Forb rootN ####

rootN.effects = forb.effects$rootN.final
rootN.effects$rootN.bt= (rootN.effects$rootN.final*4.173927) + 10.5821

attr(forb.imputed.NW.2$rootN.final, "scaled:scale")
attr(forb.imputed.NW.2$rootN.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*4.173927) + 10.5821

forb.rootN.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = rootN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Root N", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.3800000,25.6283268) +
  theme_classic()
forb.rootN.plot

ggsave("./Plots/forb.traits.NW.rootN.pdf", height = 3, width = 3)

#### Forb SLA ####

SLA.effects = forb.effects$SLA.final
SLA.effects$SLA.bt= (SLA.effects$SLA.final*8.362926) + 20.14641

attr(forb.imputed.NW.2$SLA.final, "scaled:scale")
attr(forb.imputed.NW.2$SLA.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*8.362926) + 20.14641

forb.SLA.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = SLA.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Specific Leaf Area", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.617596,46.320001)+
  theme_classic()
forb.SLA.plot

ggsave("./Plots/forb.traits.NW.SLA.pdf", height = 3, width = 3)

#### Forb Depth ####

root.depth.effects = forb.effects$root.depth.final
root.depth.effects$Depth.bt= (root.depth.effects$root.depth.final*0.6095345) + 0.6503508

attr(forb.imputed.NW.2$root.depth.final, "scaled:scale")
attr(forb.imputed.NW.2$root.depth.final, "scaled:center")

x.value = c(-1,0,1,2,3,4,5)
(x.value*0.6095345) + 0.6503508

forb.root.depth.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = Depth.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = Depth.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = Depth.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Rooting Depth", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.0479999,3.709335)+
  theme_classic()
forb.root.depth.plot

ggsave("./Plots/forb.traits.NW.root.depth.pdf", height = 3, width = 3)

#### Forb Diam ####

rootDiam.effects = forb.effects$rootDiam.final
rootDiam.effects$Diam.bt= (rootDiam.effects$rootDiam.final*0.1590925) + 0.3579696

attr(forb.imputed.NW.2$rootDiam.final, "scaled:scale")
attr(forb.imputed.NW.2$rootDiam.final, "scaled:center")

x.value = c(-1,0,1,2,3,4,5)
(x.value*0.1590925) + 0.3579696

forb.rootDiam.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = Diam.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = Diam.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = Diam.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Root Diameter", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.07840002,1.12601235)+
  theme_classic()
forb.rootDiam.plot

ggsave("./Plots/forb.traits.NW.rootDiam.pdf", height = 3, width = 3)

#### Forb SRL ####

SRL.effects = forb.effects$SRL.final
SRL.effects$SRL.bt= (SRL.effects$SRL.final*71.25811) + 105.5323

attr(forb.imputed.NW.2$SRL.final, "scaled:scale")
attr(forb.imputed.NW.2$SRL.final, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*71.25811) + 105.5323

forb.SRL.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = SRL.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Specific Root Length", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.497203,417.749988)+
  theme_classic()
forb.SRL.plot

ggsave("./Plots/forb.traits.NW.SRL.pdf", height = 3, width = 3)

#### Forb RTD ####

RTD.effects = forb.effects$RTD.final
RTD.effects$RTD.bt= (RTD.effects$RTD.final*0.1184977) + 0.227384

attr(forb.imputed.NW.2$RTD.final, "scaled:scale")
attr(forb.imputed.NW.2$RTD.final, "scaled:center")

x.value = c(-1,0,1,2,3)
(x.value*0.1184977) + 0.227384

forb.RTD.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = RTD.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Root Tissue Density", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.02325658,0.61999993)+
  theme_classic()
forb.RTD.plot

ggsave("./Plots/forb.traits.NW.RTD.pdf", height = 3, width = 3)

#### Forb RMF ####

RMF.effects = forb.effects$RMF.final
RMF.effects$RMF.bt= (RMF.effects$RMF.final*0.1176732) + 0.3992223

attr(forb.imputed.NW.2$RMF.final, "scaled:scale")
attr(forb.imputed.NW.2$RMF.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1176732) + 0.3992223

forb.RMF.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = RMF.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RMF.effects, aes(x = RMF.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RMF.effects, aes(x = RMF.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Root Mass Fraction", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.1211193,0.7328659)+
  theme_classic()
forb.RMF.plot

ggsave("./Plots/forb.traits.NW.RMF.pdf", height = 3, width = 3)

#### Forb DSI ####

DSI.effects = forb.effects$mean.DSI
DSI.effects$DSI.bt= (DSI.effects$mean.DSI*0.1949549) + -0.4644806

attr(forb.imputed.NW.2$mean.DSI, "scaled:scale")
attr(forb.imputed.NW.2$mean.DSI, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1949549) + -0.4644806

forb.DSI.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = DSI.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = DSI.effects, aes(x = DSI.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = DSI.effects, aes(x = DSI.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Drought Severity Index", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(-0.8722616,0.1242846)+
  theme_classic()
forb.DSI.plot

ggsave("./Plots/forb.traits.NW.DSI.pdf", height = 3, width = 3)

#### Forb MAP ####

MAP.effects = forb.effects$mean.MAP
MAP.effects$MAP.bt= (MAP.effects$mean.MAP*448.2603) + 657.7821

attr(forb.imputed.NW.2$mean.MAP, "scaled:scale")
attr(forb.imputed.NW.2$mean.MAP, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*448.2603) + 657.7821

forb.MAP.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = MAP.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = MAP.effects, aes(x = MAP.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = MAP.effects, aes(x = MAP.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Mean Annual Precipitation", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(132.7,2366.1)+
  theme_classic()
forb.MAP.plot

ggsave("./Plots/forb.traits.NW.MAP.pdf", height = 3, width = 3)

#### FORB interactions ####

forb.imputed.NW.traits.RTD.SRL = readRDS("./Results/forb.imputed.traits.NW.RTD.SRL.rds")

RTD.SRL.effect = conditional_effects(forb.imputed.NW.traits.RTD.SRL, effects = "RTD.final:SRL.final")$`RTD.final:SRL.final`
RTD.SRL.effect$RTD.bt= (RTD.SRL.effect$RTD.final*0.1184977) + 0.227384

# only plot high and low values

RTD.SRL.effect.2 = RTD.SRL.effect %>%
  filter(effect2__ %in% c(-1,1.05))

# RTD
x.value = c(-1,0,1,2,3)
(x.value*0.1184977) + 0.227384

# SRL
x.value = c(-1,1.05)
(x.value*71.25811) + 105.5323

forb.RTD.x.SRL.plot = ggplot(data = RTD.SRL.effect.2, aes(x = RTD.bt, y = estimate__, group = as.factor(effect2__))) +
  geom_line(aes(color = as.factor(effect2__)), size = 1.5) +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = as.factor(effect2__)), alpha = 0.5,show.legend = FALSE) +  
  geom_line(aes(y = upper__, color = factor(effect2__)), size = 1.2) +
  geom_line(aes(y = lower__, color = factor(effect2__)), size = 1.2) + 
  labs(x = "Root Tissue Density", y = "Percent Cover Change", color = "SRL") +
  scale_colour_manual(values = c("black", "#F17236"), labels = c("180.35","34.27"))+
  scale_fill_manual(values = c("black", "#F17236"))+
  xlim(0.02325658,0.61999993)+
  theme_classic()
forb.RTD.x.SRL.plot

ggsave("./Plots/forb.traits.NW.RTD.SRL.pdf", height = 3, width = 3)
ggsave("./Plots/forb.traits.NW.RTD.SRL.legend.pdf", height = 3, width = 3)

#### GRASS Backtransform ####

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)
imputed.NW.grass = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.grass.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)

grass.traits.NW.model = readRDS("./Results/grass.imputed.traits.no_woody.rds")

imputed.NW.grass$leafN.bt =(imputed.NW.grass$leafN.final*6.281442) + 18.41922
imputed.NW.grass$height.bt = (imputed.NW.grass$height.final*0.2266876) + 0.3988779
imputed.NW.grass$rootN.bt = (imputed.NW.grass$rootN.final*2.712129) + 7.59045
imputed.NW.grass$SLA.bt = (imputed.NW.grass$SLA.final*9.570923) + 20.3176
imputed.NW.grass$Depth.bt = (imputed.NW.grass$root.depth.final*0.3671355) + 0.4529895
imputed.NW.grass$Diam.bt = (imputed.NW.grass$rootDiam.final*0.1707781) + 0.3486828
imputed.NW.grass$SRL.bt = (imputed.NW.grass$SRL.final*82.36672) + 115.2106
imputed.NW.grass$RTD.bt = (imputed.NW.grass$RTD.final*0.1085128) + 0.2621843
imputed.NW.grass$RMF.bt = (imputed.NW.grass$RMF.final*0.07914561) + 0.4233237
imputed.NW.grass$DSI.bt = (imputed.NW.grass$mean.DSI*0.2007843) + -0.4496194
imputed.NW.grass$MAP.bt = (imputed.NW.grass$mean.MAP*441.5584) + 637.2953

grass.imputed.NW.2 = imputed.NW.grass %>%
  dplyr::select(cover.change,leafN.bt:MAP.bt) %>%
  drop_na()

grass.effects = conditional_effects(grass.traits.NW.model)

#### Grass leafN ####

leafN.effects = grass.effects$leafN.final
leafN.effects$leafN.bt= (leafN.effects$leafN.final*6.281442) + 18.41922

attr(grass.imputed.NW.2$leafN.final, "scaled:scale")
attr(grass.imputed.NW.2$leafN.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*6.281442) + 18.41922

grass.leafN.plot = ggplot() +
  geom_point(data = grass.imputed.NW.2, aes(x = leafN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.bt, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Leaf N", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.68951,44.19001)+
  theme_classic()
grass.leafN.plot

ggsave("./Plots/grass.traits.NW.leafN.pdf", height = 3, width = 3)

#### Grass height ####

height.effects = grass.effects$height.final
height.effects$height.bt= (height.effects$height.final*0.2266876) + 0.3988779

attr(grass.imputed.NW.2$height.final, "scaled:scale")
attr(grass.imputed.NW.2$height.final, "scaled:center")

x.value = c(-1,0,1,2,3)
(x.value*0.2266876) + 0.3988779

grass.height.plot = ggplot() +
  geom_point(data = grass.imputed.NW.2, aes(x = height.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.bt, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = height.effects, aes(x = height.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Height", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.009000064,1.322333169)+
  theme_classic()
grass.height.plot

ggsave("./Plots/grass.traits.NW.height.pdf", height = 3, width = 3)

#### Grass rootN ####

rootN.effects = grass.effects$rootN.final
rootN.effects$rootN.bt= (rootN.effects$rootN.final*2.712129) + 7.59045

attr(grass.imputed.NW.2$rootN.final, "scaled:scale")
attr(grass.imputed.NW.2$rootN.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3,4)
(x.value*2.712129) + 7.59045

grass.rootN.plot = ggplot() +
  geom_point(data = grass.imputed.NW.2, aes(x = rootN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.bt, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Root N", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.3800000,25.6283268)+
  theme_classic()
grass.rootN.plot

ggsave("./Plots/grass.traits.NW.rootN.pdf", height = 3, width = 3)

#### Grass SLA ####

SLA.effects = grass.effects$SLA.final
SLA.effects$SLA.bt= (SLA.effects$SLA.final*9.570923) + 20.3176

attr(grass.imputed.NW.2$SLA.final, "scaled:scale")
attr(grass.imputed.NW.2$SLA.final, "scaled:center")

x.value = c(-1,0,1,2,3)
(x.value*9.570923) + 20.31763

grass.SLA.plot = ggplot() +
  geom_point(data = grass.imputed.NW.2, aes(x = SLA.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.bt, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Specific Leaf Area", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.617596,46.320001)+
  theme_classic()
grass.SLA.plot

ggsave("./Plots/grass.traits.NW.SLA.pdf", height = 3, width = 3)

#### Grass Depth ####

root.depth.effects = grass.effects$root.depth.final
root.depth.effects$Depth.bt= (root.depth.effects$root.depth.final*0.3671355) + 0.4529895

attr(grass.imputed.NW.2$root.depth.final, "scaled:scale")
attr(grass.imputed.NW.2$root.depth.final, "scaled:center")

x.value = c(-1,0,1,2,3,4,5,6)
(x.value*0.3671355) + 0.4529895

grass.root.depth.plot = ggplot() +
  geom_point(data = grass.imputed.NW.2, aes(x = Depth.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = Depth.bt, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = Depth.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Rooting Depth", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.0479999,3.709335)+
  theme_classic()
grass.root.depth.plot

ggsave("./Plots/grass.traits.NW.root.depth.pdf", height = 3, width = 3)

#### Grass Diam ####

rootDiam.effects = grass.effects$rootDiam.final
rootDiam.effects$Diam.bt= (rootDiam.effects$rootDiam.final*0.1707781) + 0.3486828

attr(grass.imputed.NW.2$rootDiam.final, "scaled:scale")
attr(grass.imputed.NW.2$rootDiam.final, "scaled:center")

x.value = c(-1,0,1,2,3)
(x.value*0.1707781) + 0.3486828

grass.rootDiam.plot = ggplot() +
  geom_point(data = grass.imputed.NW.2, aes(x = Diam.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = Diam.bt, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = Diam.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Root Diameter", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.07840002,1.12601235)+
  theme_classic()
grass.rootDiam.plot

ggsave("./Plots/grass.traits.NW.rootDiam.pdf", height = 3, width = 3)

#### Grass SRL ####

SRL.effects = grass.effects$SRL.final
SRL.effects$SRL.bt= (SRL.effects$SRL.final*82.36672) + 115.2106

attr(grass.imputed.NW.2$SRL.final, "scaled:scale")
attr(grass.imputed.NW.2$SRL.final, "scaled:center")

x.value = c(-1,0,1,2,3)
(x.value*82.36672) + 115.2106

grass.SRL.plot = ggplot() +
  geom_point(data = grass.imputed.NW.2, aes(x = SRL.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.bt, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Specific Root Length", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.497203,417.749988)+
  theme_classic()
grass.SRL.plot

ggsave("./Plots/grass.traits.NW.SRL.pdf", height = 3, width = 3)

#### Grass RTD ####

RTD.effects = grass.effects$RTD.final
RTD.effects$RTD.bt= (RTD.effects$RTD.final*0.1085128) + 0.2621843

attr(grass.imputed.NW.2$RTD.final, "scaled:scale")
attr(grass.imputed.NW.2$RTD.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1085128) + 0.2621843

grass.RTD.plot = ggplot() +
  geom_point(data = grass.imputed.NW.2, aes(x = RTD.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.bt, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Root Tissue Density", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.02325658,0.61999993)+
  theme_classic()
grass.RTD.plot

ggsave("./Plots/grass.traits.NW.RTD.pdf", height = 3, width = 3)

#### Grass RMF ####

RMF.effects = grass.effects$RMF.final
RMF.effects$RMF.bt= (RMF.effects$RMF.final*0.07914561) + 0.4233237

attr(grass.imputed.NW.2$RMF.final, "scaled:scale")
attr(grass.imputed.NW.2$RMF.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.07914561) + 0.4233237

grass.RMF.plot = ggplot() +
  geom_point(data = grass.imputed.NW.2, aes(x = RMF.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RMF.effects, aes(x = RMF.bt, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RMF.effects, aes(x = RMF.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Root Mass Fraction", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.1211193,0.7328659)+
  theme_classic()
grass.RMF.plot

ggsave("./Plots/grass.traits.NW.RMF.pdf", height = 3, width = 3)

#### Grass DSI ####

DSI.effects = grass.effects$mean.DSI
DSI.effects$DSI.bt= (DSI.effects$mean.DSI*0.2007843) + -0.4496194

attr(grass.imputed.NW.2$mean.DSI, "scaled:scale")
attr(grass.imputed.NW.2$mean.DSI, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.2007843) + -0.4496194

grass.DSI.plot = ggplot() +
  geom_point(data = grass.imputed.NW.2, aes(x = DSI.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = DSI.effects, aes(x = DSI.bt, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = DSI.effects, aes(x = DSI.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Drought Severity Index", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(-0.8722616,0.1242846)+
  theme_classic()
grass.DSI.plot

ggsave("./Plots/grass.traits.NW.DSI.pdf", height = 3, width = 3)

#### Grass MAP ####

MAP.effects = grass.effects$mean.MAP
MAP.effects$MAP.bt= (MAP.effects$mean.MAP*441.5584) + 637.2953

attr(grass.imputed.NW.2$mean.MAP, "scaled:scale")
attr(grass.imputed.NW.2$mean.MAP, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*441.5584) + 637.2953

grass.MAP.plot = ggplot() +
  geom_point(data = grass.imputed.NW.2, aes(x = MAP.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = MAP.effects, aes(x = MAP.bt, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = MAP.effects, aes(x = MAP.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Mean Annual Precipitation", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(132.7,2366.1)+
  theme_classic()
grass.MAP.plot

ggsave("./Plots/grass.traits.NW.MAP.pdf", height = 3, width = 3)

#### PERENNIAL GRASS Backtransform ####

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)
imputed.NW.perennial.grass = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennial.grass.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)

perennial.grass.traits.NW.model = readRDS("./Results/perennial.grass.imputed.traits.no_woody.rds")

imputed.NW.perennial.grass$leafN.bt =(imputed.NW.perennial.grass$leafN.final*5.50471) + 17.71951
imputed.NW.perennial.grass$height.bt = (imputed.NW.perennial.grass$height.final*0.23344) + 0.4066261
imputed.NW.perennial.grass$rootN.bt = (imputed.NW.perennial.grass$rootN.final*2.940311) + 7.794231
imputed.NW.perennial.grass$SLA.bt = (imputed.NW.perennial.grass$SLA.final*9.175084) + 19.32517
imputed.NW.perennial.grass$Depth.bt = (imputed.NW.perennial.grass$root.depth.final*0.3665233) + 0.4460136
imputed.NW.perennial.grass$Diam.bt = (imputed.NW.perennial.grass$rootDiam.final*0.1693868) + 0.3418054
imputed.NW.perennial.grass$SRL.bt = (imputed.NW.perennial.grass$SRL.final*85.10344) + 113.0663
imputed.NW.perennial.grass$RTD.bt = (imputed.NW.perennial.grass$RTD.final*0.1033287) + 0.2760841
imputed.NW.perennial.grass$RMF.bt = (imputed.NW.perennial.grass$RMF.final*0.08185729) + 0.4278072
imputed.NW.perennial.grass$DSI.bt = (imputed.NW.perennial.grass$mean.DSI*0.1851188) + -0.4469141
imputed.NW.perennial.grass$MAP.bt = (imputed.NW.perennial.grass$mean.MAP*466.5336) + 686.9575

perennial.grass.imputed.NW.2 = imputed.NW.perennial.grass %>%
  dplyr::select(cover.change,leafN.bt:MAP.bt) %>%
  drop_na()

perennial.grass.effects = conditional_effects(perennial.grass.traits.NW.model)

#### Perennial Grass leafN ####

leafN.effects = perennial.grass.effects$leafN.final
leafN.effects$leafN.bt= (leafN.effects$leafN.final*5.50471) + 17.71951

x.value = c(-2,-1,0,1,2,3)
(x.value*5.50471) + 17.71951

perennial.grass.leafN.plot = ggplot() +
  geom_point(data = perennial.grass.imputed.NW.2, aes(x = leafN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.bt, y = estimate__), color = "#6089B5", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6089B5") +  
  labs(x = "Leaf N", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.68951,44.19001)+
  theme_classic()
perennial.grass.leafN.plot

ggsave("./Plots/perennial.grass.traits.NW.leafN.pdf", height = 3, width = 3)

#### Perennial Grass height ####

height.effects = perennial.grass.effects$height.final
height.effects$height.bt= (height.effects$height.final*0.23344) + 0.4066261

x.value = c(-1,0,1,2,3)
(x.value*0.23344) + 0.4066261

perennial.grass.height.plot = ggplot() +
  geom_point(data = perennial.grass.imputed.NW.2, aes(x = height.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.bt, y = estimate__), color = "#6089B5", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = height.effects, aes(x = height.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6089B5") +  
  labs(x = "Height", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.009000064,1.322333169)+
  theme_classic()
perennial.grass.height.plot

ggsave("./Plots/perennial.grass.traits.NW.height.pdf", height = 3, width = 3)

#### Perennial Grass rootN ####

rootN.effects = perennial.grass.effects$rootN.final
rootN.effects$rootN.bt= (rootN.effects$rootN.final*2.940311) + 7.794231

x.value = c(-2,-1,0,1,2,3,4)
(x.value*2.940311) + 7.794231

perennial.grass.rootN.plot = ggplot() +
  geom_point(data = perennial.grass.imputed.NW.2, aes(x = rootN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.bt, y = estimate__), color = "#6089B5", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6089B5") +  
  labs(x = "Root N", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.3800000,25.6283268)+
  theme_classic()
perennial.grass.rootN.plot

ggsave("./Plots/perennial.grass.traits.NW.rootN.pdf", height = 3, width = 3)

#### Perennial Grass SLA ####

SLA.effects = perennial.grass.effects$SLA.final
SLA.effects$SLA.bt= (SLA.effects$SLA.final*9.175084) + 19.32517

x.value = c(-1,0,1,2,3)
(x.value*9.175084) + 19.32517

perennial.grass.SLA.plot = ggplot() +
  geom_point(data = perennial.grass.imputed.NW.2, aes(x = SLA.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.bt, y = estimate__), color = "#6089B5", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6089B5") +  
  labs(x = "Specific Leaf Area", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.617596,46.320001)+
  theme_classic()
perennial.grass.SLA.plot

ggsave("./Plots/perennial.grass.traits.NW.SLA.pdf", height = 3, width = 3)

#### Perennial Grass depth ####

root.depth.effects = perennial.grass.effects$root.depth.final
root.depth.effects$Depth.bt= (root.depth.effects$root.depth.final*0.3665233) + 0.4460136

x.value = c(-1,0,1,2,3,4,5,6)
(x.value*0.3665233) + 0.4460136

perennial.grass.root.depth.plot = ggplot() +
  geom_point(data = perennial.grass.imputed.NW.2, aes(x = Depth.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = Depth.bt, y = estimate__), color = "#6089B5", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = Depth.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6089B5") +  
  labs(x = "Rooting Depth", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.0479999,3.709335)+
  theme_classic()
perennial.grass.root.depth.plot

ggsave("./Plots/perennial.grass.traits.NW.root.depth.pdf", height = 3, width = 3)

#### Perennial Grass Diam ####

rootDiam.effects = perennial.grass.effects$rootDiam.final
rootDiam.effects$Diam.bt= (rootDiam.effects$rootDiam.final*0.1693868) + 0.3418054

x.value = c(-1,0,1,2,3)
(x.value*0.1693868) + 0.3418054

perennial.grass.rootDiam.plot = ggplot() +
  geom_point(data = perennial.grass.imputed.NW.2, aes(x = Diam.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = Diam.bt, y = estimate__), color = "#6089B5", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = Diam.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6089B5") +  
  labs(x = "Root Diameter", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.07840002,1.12601235)+
  theme_classic()
perennial.grass.rootDiam.plot

ggsave("./Plots/perennial.grass.traits.NW.rootDiam.pdf", height = 3, width = 3)

#### Perennial Grass SRL ####

SRL.effects = perennial.grass.effects$SRL.final
SRL.effects$SRL.bt= (SRL.effects$SRL.final*85.10344) + 113.0663

x.value = c(-1,0,1,2,3)
(x.value*85.10344) + 113.0663

perennial.grass.SRL.plot = ggplot() +
  geom_point(data = perennial.grass.imputed.NW.2, aes(x = SRL.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.bt, y = estimate__), color = "#6089B5", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6089B5") +  
  labs(x = "Specific Root Length", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.497203,417.749988)+
  theme_classic()
perennial.grass.SRL.plot

ggsave("./Plots/perennial.grass.traits.NW.SRL.pdf", height = 3, width = 3)

#### Perennial Grass RTD ####

RTD.effects = perennial.grass.effects$RTD.final
RTD.effects$RTD.bt= (RTD.effects$RTD.final*0.1033287) + 0.2760841

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1033287) + 0.2760841

perennial.grass.RTD.plot = ggplot() +
  geom_point(data = perennial.grass.imputed.NW.2, aes(x = RTD.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.bt, y = estimate__), color = "#6089B5", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6089B5") +  
  labs(x = "Root Tissue Density", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.02325658,0.61999993)+
  theme_classic()
perennial.grass.RTD.plot

ggsave("./Plots/perennial.grass.traits.NW.RTD.pdf", height = 3, width = 3)

#### Perennial Grass RMF ####

RMF.effects = perennial.grass.effects$RMF.final
RMF.effects$RMF.bt= (RMF.effects$RMF.final*0.08185729) + 0.4278072

x.value = c(-2,-1,0,1,2,3)
(x.value*0.08185729) + 0.4278072

perennial.grass.RMF.plot = ggplot() +
  geom_point(data = perennial.grass.imputed.NW.2, aes(x = RMF.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RMF.effects, aes(x = RMF.bt, y = estimate__), color = "#6089B5", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RMF.effects, aes(x = RMF.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6089B5") +  
  labs(x = "Root Mass Fraction", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.1211193,0.7328659)+
  theme_classic()
perennial.grass.RMF.plot

ggsave("./Plots/perennial.grass.traits.NW.RMF.pdf", height = 3, width = 3)

#### Perennial Grass DSI ####

DSI.effects = perennial.grass.effects$mean.DSI
DSI.effects$DSI.bt= (DSI.effects$mean.DSI*0.1851188) + -0.4469141

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1851188) + -0.4469141

perennial.grass.DSI.plot = ggplot() +
  geom_point(data = perennial.grass.imputed.NW.2, aes(x = DSI.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = DSI.effects, aes(x = DSI.bt, y = estimate__), color = "#6089B5", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = DSI.effects, aes(x = DSI.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6089B5") +  
  labs(x = "Drought Severity Index", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(-0.8722616,0.1242846)+
  theme_classic()
perennial.grass.DSI.plot

ggsave("./Plots/perennial.grass.traits.NW.DSI.pdf", height = 3, width = 3)

#### Perennial Grass MAP ####

MAP.effects = perennial.grass.effects$mean.MAP
MAP.effects$MAP.bt= (MAP.effects$mean.MAP*466.5336) + 686.9575

x.value = c(-1,0,1,2,3,4)
(x.value*466.5336) + 686.9575

perennial.grass.MAP.plot = ggplot() +
  geom_point(data = perennial.grass.imputed.NW.2, aes(x = MAP.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = MAP.effects, aes(x = MAP.bt, y = estimate__), color = "#6089B5", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = MAP.effects, aes(x = MAP.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6089B5") +  
  labs(x = "Mean Annual Precipitation", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(132.7,2366.1)+
  theme_classic()
perennial.grass.MAP.plot

ggsave("./Plots/perennial.grass.traits.NW.MAP.pdf", height = 3, width = 3)

#### PERENNIAL forb Backtransform ####

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)
imputed.NW.perennial.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennial.forb.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)

perennial.forb.traits.NW.model = readRDS("./Results/perennial.forb.imputed.traits.no_woody.rds")

imputed.NW.perennial.forb$leafN.bt =(imputed.NW.perennial.forb$leafN.final*6.163125) + 21.2703
imputed.NW.perennial.forb$height.bt = (imputed.NW.perennial.forb$height.final*0.2398901) + 0.3335025
imputed.NW.perennial.forb$rootN.bt = (imputed.NW.perennial.forb$rootN.final*4.060239) + 10.55429
imputed.NW.perennial.forb$SLA.bt = (imputed.NW.perennial.forb$SLA.final*7.330182) + 18.8354
imputed.NW.perennial.forb$Depth.bt = (imputed.NW.perennial.forb$root.depth.final*0.6839351) + 0.6887514
imputed.NW.perennial.forb$Diam.bt = (imputed.NW.perennial.forb$rootDiam.final*0.1702321) + 0.3728557
imputed.NW.perennial.forb$SRL.bt = (imputed.NW.perennial.forb$SRL.final*66.81515) + 92.88134
imputed.NW.perennial.forb$RTD.bt = (imputed.NW.perennial.forb$RTD.final*0.1226583) + 0.241658
imputed.NW.perennial.forb$RMF.bt = (imputed.NW.perennial.forb$RMF.final*0.1195843) + 0.4031475
imputed.NW.perennial.forb$DSI.bt = (imputed.NW.perennial.forb$mean.DSI*0.1767211) + -0.4540636
imputed.NW.perennial.forb$MAP.bt = (imputed.NW.perennial.forb$mean.MAP*483.3507) + 751.5249

perennial.forb.imputed.NW.2 = imputed.NW.perennial.forb %>%
  dplyr::select(cover.change,leafN.bt:MAP.bt) %>%
  drop_na()

perennial.forb.effects = conditional_effects(perennial.forb.traits.NW.model)

#### Perennial forb leafN ####

leafN.effects = perennial.forb.effects$leafN.final
leafN.effects$leafN.bt= (leafN.effects$leafN.final*6.163125) + 21.2703

x.value = c(-2,-1,0,1,2,3)
(x.value*6.163125) + 21.2703

perennial.forb.leafN.plot = ggplot() +
  geom_point(data = perennial.forb.imputed.NW.2, aes(x = leafN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.bt, y = estimate__), color = "black", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "black") +  
  labs(x = "Leaf N", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.68951,44.19001)+
  theme_classic()
perennial.forb.leafN.plot

ggsave("./Plots/perennial.forb.traits.NW.leafN.pdf", height = 3, width = 3)

#### Perennial forb height ####

height.effects = perennial.forb.effects$height.final
height.effects$height.bt= (height.effects$height.final*0.2398901) + 0.3335025

x.value = c(-1,0,1,2,3)
(x.value*0.2398901) + 0.3335025

perennial.forb.height.plot = ggplot() +
  geom_point(data = perennial.forb.imputed.NW.2, aes(x = height.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.bt, y = estimate__), color = "black", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = height.effects, aes(x = height.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "black") +  
  labs(x = "Height", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.009000064,1.322333169)+
  theme_classic()
perennial.forb.height.plot

ggsave("./Plots/perennial.forb.traits.NW.height.pdf", height = 3, width = 3)

#### Perennial forb rootN ####

rootN.effects = perennial.forb.effects$rootN.final
rootN.effects$rootN.bt= (rootN.effects$rootN.final*4.060239) + 10.55429

x.value = c(-2,-1,0,1,2,3,4)
(x.value*4.060239) + 10.55429

perennial.forb.rootN.plot = ggplot() +
  geom_point(data = perennial.forb.imputed.NW.2, aes(x = rootN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.bt, y = estimate__), color = "black", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "black") +  
  labs(x = "Root N", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.3800000,25.6283268)+
  theme_classic()
perennial.forb.rootN.plot

ggsave("./Plots/perennial.forb.traits.NW.rootN.pdf", height = 3, width = 3)

#### Perennial forb SLA ####

SLA.effects = perennial.forb.effects$SLA.final
SLA.effects$SLA.bt= (SLA.effects$SLA.final*7.330182) + 18.8354

x.value = c(-1,0,1,2,3)
(x.value*7.330182) + 18.8354

perennial.forb.SLA.plot = ggplot() +
  geom_point(data = perennial.forb.imputed.NW.2, aes(x = SLA.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.bt, y = estimate__), color = "black", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "black") +  
  labs(x = "Specific Leaf Area", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.617596,46.320001)+
  theme_classic()
perennial.forb.SLA.plot

ggsave("./Plots/perennial.forb.traits.NW.SLA.pdf", height = 3, width = 3)

#### Perennial forb depth ####

root.depth.effects = perennial.forb.effects$root.depth.final
root.depth.effects$Depth.bt= (root.depth.effects$root.depth.final*0.6839351) + 0.6887514

x.value = c(-1,0,1,2,3,4,5,6)
(x.value*0.6839351) + 0.6887514

perennial.forb.root.depth.plot = ggplot() +
  geom_point(data = perennial.forb.imputed.NW.2, aes(x = Depth.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = Depth.bt, y = estimate__), color = "black", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = Depth.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "black") +  
  labs(x = "Rooting Depth", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.0479999,3.709335)+
  theme_classic()
perennial.forb.root.depth.plot

ggsave("./Plots/perennial.forb.traits.NW.root.depth.pdf", height = 3, width = 3)

#### Perennial forb Diam ####

rootDiam.effects = perennial.forb.effects$rootDiam.final
rootDiam.effects$Diam.bt= (rootDiam.effects$rootDiam.final*0.1702321) + 0.3728557

x.value = c(-1,0,1,2,3)
(x.value*0.1702321) + 0.3728557

perennial.forb.rootDiam.plot = ggplot() +
  geom_point(data = perennial.forb.imputed.NW.2, aes(x = Diam.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = Diam.bt, y = estimate__), color = "black", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = Diam.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "black") +  
  labs(x = "Root Diameter", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.07840002,1.12601235)+
  theme_classic()
perennial.forb.rootDiam.plot

ggsave("./Plots/perennial.forb.traits.NW.rootDiam.pdf", height = 3, width = 3)

#### Perennial forb SRL ####

SRL.effects = perennial.forb.effects$SRL.final
SRL.effects$SRL.bt= (SRL.effects$SRL.final*66.81515) + 92.88134

x.value = c(-1,0,1,2,3)
(x.value*66.81515) + 92.88134

perennial.forb.SRL.plot = ggplot() +
  geom_point(data = perennial.forb.imputed.NW.2, aes(x = SRL.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.bt, y = estimate__), color = "black", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "black") +  
  labs(x = "Specific Root Length", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.497203,417.749988)+
  theme_classic()
perennial.forb.SRL.plot

ggsave("./Plots/perennial.forb.traits.NW.SRL.pdf", height = 3, width = 3)

#### Perennial forb RTD ####

RTD.effects = perennial.forb.effects$RTD.final
RTD.effects$RTD.bt= (RTD.effects$RTD.final*0.1226583) + 0.241658

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1226583) + 0.241658

perennial.forb.RTD.plot = ggplot() +
  geom_point(data = perennial.forb.imputed.NW.2, aes(x = RTD.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.bt, y = estimate__), color = "black", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "black") +  
  labs(x = "Root Tissue Density", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.02325658,0.61999993)+
  theme_classic()
perennial.forb.RTD.plot

ggsave("./Plots/perennial.forb.traits.NW.RTD.pdf", height = 3, width = 3)

#### Perennial forb RMF ####

RMF.effects = perennial.forb.effects$RMF.final
RMF.effects$RMF.bt= (RMF.effects$RMF.final*0.1195843) + 0.4031475

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1195843) + 0.4031475

perennial.forb.RMF.plot = ggplot() +
  geom_point(data = perennial.forb.imputed.NW.2, aes(x = RMF.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RMF.effects, aes(x = RMF.bt, y = estimate__), color = "black", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RMF.effects, aes(x = RMF.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "black") +  
  labs(x = "Root Mass Fraction", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.1211193,0.7328659)+
  theme_classic()
perennial.forb.RMF.plot

ggsave("./Plots/perennial.forb.traits.NW.RMF.pdf", height = 3, width = 3)

#### Perennial forb DSI ####

DSI.effects = perennial.forb.effects$mean.DSI
DSI.effects$DSI.bt= (DSI.effects$mean.DSI*0.1767211) + -0.4540636

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1767211) + -0.4540636

perennial.forb.DSI.plot = ggplot() +
  geom_point(data = perennial.forb.imputed.NW.2, aes(x = DSI.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = DSI.effects, aes(x = DSI.bt, y = estimate__), color = "black", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = DSI.effects, aes(x = DSI.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "black") +  
  labs(x = "Drought Severity Index", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(-0.8722616,0.1242846)+
  theme_classic()
perennial.forb.DSI.plot

ggsave("./Plots/perennial.forb.traits.NW.DSI.pdf", height = 3, width = 3)

#### Perennial forb MAP ####

MAP.effects = perennial.forb.effects$mean.MAP
MAP.effects$MAP.bt= (MAP.effects$mean.MAP*483.3507) + 751.5249

x.value = c(-1,0,1,2,3,4)
(x.value*483.3507) + 751.5249

perennial.forb.MAP.plot = ggplot() +
  geom_point(data = perennial.forb.imputed.NW.2, aes(x = MAP.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = MAP.effects, aes(x = MAP.bt, y = estimate__), color = "black", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = MAP.effects, aes(x = MAP.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "black") +  
  labs(x = "Mean Annual Precipitation", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(132.7,2366.1)+
  theme_classic()
perennial.forb.MAP.plot

ggsave("./Plots/perennial.forb.traits.NW.MAP.pdf", height = 3, width = 3)

#### ANNUAL FORB Backtransform ####

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)
imputed.NW.annual.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.annual.forbs.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)

annual.forb.traits.NW.model = readRDS("./Results/annual.forb.imputed.traits.no_woody.rds")

imputed.NW.annual.forb$leafN.bt =(imputed.NW.annual.forb$leafN.final*6.789956) + 22.1907
imputed.NW.annual.forb$height.bt = (imputed.NW.annual.forb$height.final*0.2390744) + 0.3002752
imputed.NW.annual.forb$rootN.bt = (imputed.NW.annual.forb$rootN.final*4.242741) + 10.49856
imputed.NW.annual.forb$SLA.bt = (imputed.NW.annual.forb$SLA.final*9.43772) + 22.54378
imputed.NW.annual.forb$Depth.bt = (imputed.NW.annual.forb$root.depth.final*0.4594328) + 0.576239
imputed.NW.annual.forb$Diam.bt = (imputed.NW.annual.forb$rootDiam.final*0.1309004) + 0.33456
imputed.NW.annual.forb$SRL.bt = (imputed.NW.annual.forb$SRL.final*73.2686) + 124.8044
imputed.NW.annual.forb$RTD.bt = (imputed.NW.annual.forb$RTD.final*0.110191) + 0.2068199
imputed.NW.annual.forb$RMF.bt = (imputed.NW.annual.forb$RMF.final*0.1125971) + 0.3882978
imputed.NW.annual.forb$DSI.bt = (imputed.NW.annual.forb$mean.DSI*0.2195586) + -0.4745653
imputed.NW.annual.forb$MAP.bt = (imputed.NW.annual.forb$mean.MAP*297.4522) + 478.9622

annual.forb.imputed.NW.2 = imputed.NW.annual.forb %>%
  dplyr::select(cover.change,leafN.bt:MAP.bt) %>%
  drop_na()

annual.forb.effects = conditional_effects(annual.forb.traits.NW.model)

#### annual forb leafN ####

leafN.effects = annual.forb.effects$leafN.final
leafN.effects$leafN.bt= (leafN.effects$leafN.final*6.789956) + 22.1907

x.value = c(-2,-1,0,1,2,3)
(x.value*6.789956) + 22.1907

annual.forb.leafN.plot = ggplot() +
  geom_point(data = annual.forb.imputed.NW.2, aes(x = leafN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.bt, y = estimate__), color = "#B50200", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#B50200") +  
  labs(x = "Leaf N", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.68951,44.19001)+
  theme_classic()
annual.forb.leafN.plot

ggsave("./Plots/annual.forb.traits.NW.leafN.pdf", height = 3, width = 3)

#### annual forb height ####

height.effects = annual.forb.effects$height.final
height.effects$height.bt= (height.effects$height.final*0.2390744) + 0.3002752

x.value = c(-1,0,1,2,3)
(x.value*0.2390744) + 0.3002752

annual.forb.height.plot = ggplot() +
  geom_point(data = annual.forb.imputed.NW.2, aes(x = height.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.bt, y = estimate__), color = "#B50200", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = height.effects, aes(x = height.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#B50200") +  
  labs(x = "Height", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.009000064,1.322333169)+
  theme_classic()
annual.forb.height.plot

ggsave("./Plots/annual.forb.traits.NW.height.pdf", height = 3, width = 3)

#### annual forb rootN ####

rootN.effects = annual.forb.effects$rootN.final
rootN.effects$rootN.bt= (rootN.effects$rootN.final*4.242741) + 10.49856

x.value = c(-2,-1,0,1,2,3,4)
(x.value*4.242741) + 10.49856

annual.forb.rootN.plot = ggplot() +
  geom_point(data = annual.forb.imputed.NW.2, aes(x = rootN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.bt, y = estimate__), color = "#B50200", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#B50200") +  
  labs(x = "Root N", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.3800000,25.6283268)+
  theme_classic()
annual.forb.rootN.plot

ggsave("./Plots/annual.forb.traits.NW.rootN.pdf", height = 3, width = 3)

#### annual forb SLA ####

SLA.effects = annual.forb.effects$SLA.final
SLA.effects$SLA.bt= (SLA.effects$SLA.final*9.43772) + 22.54378

x.value = c(-1,0,1,2,3)
(x.value*9.43772) + 22.54378

annual.forb.SLA.plot = ggplot() +
  geom_point(data = annual.forb.imputed.NW.2, aes(x = SLA.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.bt, y = estimate__), color = "#B50200", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#B50200") +  
  labs(x = "Specific Leaf Area", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.617596,46.320001)+
  theme_classic()
annual.forb.SLA.plot

ggsave("./Plots/annual.forb.traits.NW.SLA.pdf", height = 3, width = 3)

#### annual forb depth ####

root.depth.effects = annual.forb.effects$root.depth.final
root.depth.effects$Depth.bt= (root.depth.effects$root.depth.final*0.4594328) + 0.576239

x.value = c(-1,0,1,2,3,4,5,6)
(x.value*0.4594328) + 0.576239

annual.forb.root.depth.plot = ggplot() +
  geom_point(data = annual.forb.imputed.NW.2, aes(x = Depth.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = Depth.bt, y = estimate__), color = "#B50200", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = Depth.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#B50200") +  
  labs(x = "Rooting Depth", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.0479999,3.709335)+
  theme_classic()
annual.forb.root.depth.plot

ggsave("./Plots/annual.forb.traits.NW.root.depth.pdf", height = 3, width = 3)

#### annual forb Diam ####

rootDiam.effects = annual.forb.effects$rootDiam.final
rootDiam.effects$Diam.bt= (rootDiam.effects$rootDiam.final*0.1309004) + 0.33456

x.value = c(-1,0,1,2,3)
(x.value*0.1309004) + 0.33456

annual.forb.rootDiam.plot = ggplot() +
  geom_point(data = annual.forb.imputed.NW.2, aes(x = Diam.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = Diam.bt, y = estimate__), color = "#B50200", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = Diam.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#B50200") +  
  labs(x = "Root Diameter", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.07840002,1.12601235)+
  theme_classic()
annual.forb.rootDiam.plot

ggsave("./Plots/annual.forb.traits.NW.rootDiam.pdf", height = 3, width = 3)

#### annual forb SRL ####

SRL.effects = annual.forb.effects$SRL.final
SRL.effects$SRL.bt= (SRL.effects$SRL.final*73.2686) + 124.8044

x.value = c(-1,0,1,2,3)
(x.value*73.2686) + 124.8044

annual.forb.SRL.plot = ggplot() +
  geom_point(data = annual.forb.imputed.NW.2, aes(x = SRL.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.bt, y = estimate__), color = "#B50200", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#B50200") +  
  labs(x = "Specific Root Length", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(2.497203,417.749988)+
  theme_classic()
annual.forb.SRL.plot

ggsave("./Plots/annual.forb.traits.NW.SRL.pdf", height = 3, width = 3)

#### annual forb RTD ####

RTD.effects = annual.forb.effects$RTD.final
RTD.effects$RTD.bt= (RTD.effects$RTD.final*0.110191) + 0.2068199

x.value = c(-2,-1,0,1,2,3)
(x.value*0.110191) + 0.2068199

annual.forb.RTD.plot = ggplot() +
  geom_point(data = annual.forb.imputed.NW.2, aes(x = RTD.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.bt, y = estimate__), color = "#B50200", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#B50200") +  
  labs(x = "Root Tissue Density", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.02325658,0.61999993)+
  theme_classic()
annual.forb.RTD.plot

ggsave("./Plots/annual.forb.traits.NW.RTD.pdf", height = 3, width = 3)

#### annual forb RMF ####

RMF.effects = annual.forb.effects$RMF.final
RMF.effects$RMF.bt= (RMF.effects$RMF.final*0.1125971) + 0.3882978

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1125971) + 0.3882978

annual.forb.RMF.plot = ggplot() +
  geom_point(data = annual.forb.imputed.NW.2, aes(x = RMF.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RMF.effects, aes(x = RMF.bt, y = estimate__), color = "#B50200", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RMF.effects, aes(x = RMF.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#B50200") +  
  labs(x = "Root Mass Fraction", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(0.1211193,0.7328659)+
  theme_classic()
annual.forb.RMF.plot

ggsave("./Plots/annual.forb.traits.NW.RMF.pdf", height = 3, width = 3)

#### annual forb DSI ####

DSI.effects = annual.forb.effects$mean.DSI
DSI.effects$DSI.bt= (DSI.effects$mean.DSI*0.2195586) + -0.4745653

x.value = c(-2,-1,0,1,2,3)
(x.value*0.2195586) + -0.4745653

annual.forb.DSI.plot = ggplot() +
  geom_point(data = annual.forb.imputed.NW.2, aes(x = DSI.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = DSI.effects, aes(x = DSI.bt, y = estimate__), color = "#B50200", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = DSI.effects, aes(x = DSI.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#B50200") +  
  labs(x = "Drought Severity Index", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(-0.8722616,0.1242846)+
  theme_classic()
annual.forb.DSI.plot

ggsave("./Plots/annual.forb.traits.NW.DSI.pdf", height = 3, width = 3)

#### annual forb MAP ####

MAP.effects = annual.forb.effects$mean.MAP
MAP.effects$MAP.bt= (MAP.effects$mean.MAP*297.4522) + 478.9622

x.value = c(-1,0,1,2,3,4)
(x.value*297.4522) + 478.9622

annual.forb.MAP.plot = ggplot() +
  geom_point(data = annual.forb.imputed.NW.2, aes(x = MAP.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = MAP.effects, aes(x = MAP.bt, y = estimate__), color = "#B50200", size = 1.5) +  
  geom_ribbon(data = MAP.effects, aes(x = MAP.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#B50200") +  
  labs(x = "Mean Annual Precipitation", y = "Percent Cover Change") +
  ylim(-20.66667,21.24500)+
  xlim(132.7,2366.1)+
  theme_classic()
annual.forb.MAP.plot

ggsave("./Plots/annual.forb.traits.NW.MAP.pdf", height = 3, width = 3)


#### Dot Plots ####
imputed.traits.NW.model = readRDS("./Results/all.imputed.traits.no_woody.rds")
sum = summary(imputed.traits.NW.model)
coef = sum$fixed
row.names(coef) = c("Intercept","Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter",
                    "SRL", "RTD","RMF","DSI","Precipitation")
coef.2 = coef %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#769370"))
  
All.species.plot = ggplot(data = coef.2,
       aes(y = factor(row.names(coef.2), levels = rev(row.names(coef.2))),
           x = Estimate,
           xmin = `l-95% CI`,
           xmax = `u-95% CI`,
           fill = resp.fill)) +
  geom_pointrange(size = 1, shape = 21, color = "#769370") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "All-Species")
All.species.plot

annual.traits.NW.model = readRDS("./Results/annual.imputed.traits.no_woody.rds")

sum = summary(annual.traits.NW.model)
coef = sum$fixed
row.names(coef) = c("Intercept","Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter",
                    "SRL", "RTD","RMF","DSI","Precipitation")
coef.annual = coef %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#979461"))

annual.plot = ggplot(data = coef.annual,
                          aes(y = factor(row.names(coef.annual), levels = rev(row.names(coef.annual))),
                              x = Estimate,
                              xmin = `l-95% CI`,
                              xmax = `u-95% CI`,
                              fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#979461") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Annuals")
annual.plot

perennial.traits.NW.model = readRDS("./Results/perennial.imputed.traits.no_woody.rds")

sum = summary(perennial.traits.NW.model)
coef = sum$fixed
row.names(coef) = c("Intercept","Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter",
                    "SRL", "RTD","RMF","DSI","Precipitation")
coef.perennial = coef %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#F1C646"))

perennial.plot = ggplot(data = coef.perennial,
                     aes(y = factor(row.names(coef.perennial), levels = rev(row.names(coef.perennial))),
                         x = Estimate,
                         xmin = `l-95% CI`,
                         xmax = `u-95% CI`,
                         fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#F1C646") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Perennials")
perennial.plot

forb.traits.NW.model = readRDS("./Results/forb.imputed.traits.no_woody.rds")

sum = summary(forb.traits.NW.model)
coef = sum$fixed
row.names(coef) = c("Intercept","Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter",
                    "SRL", "RTD","RMF","DSI","Precipitation")
coef.forb = coef %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#F17236"))

forb.plot = ggplot(data = coef.forb,
                        aes(y = factor(row.names(coef.forb), levels = rev(row.names(coef.forb))),
                            x = Estimate,
                            xmin = `l-95% CI`,
                            xmax = `u-95% CI`,
                            fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#F17236") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Forbs")
forb.plot

grass.traits.NW.model = readRDS("./Results/grass.imputed.traits.no_woody.rds")

sum = summary(grass.traits.NW.model)
coef = sum$fixed
row.names(coef) = c("Intercept","Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter",
                    "SRL", "RTD","RMF","DSI","Precipitation")
coef.grass = coef %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#6E687E"))

grass.plot = ggplot(data = coef.grass,
                   aes(y = factor(row.names(coef.grass), levels = rev(row.names(coef.grass))),
                       x = Estimate,
                       xmin = `l-95% CI`,
                       xmax = `u-95% CI`,
                       fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#6E687E") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Grasses")
grass.plot

annual.forb.traits.NW.model = readRDS("./Results/annual.forb.imputed.traits.no_woody.rds")

sum = summary(annual.forb.traits.NW.model)
coef = sum$fixed
row.names(coef) = c("Intercept","Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter",
                    "SRL", "RTD","RMF","DSI","Precipitation")
coef.annual.forb = coef %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#B50200"))

annual.forb.plot = ggplot(data = coef.annual.forb,
                    aes(y = factor(row.names(coef.annual.forb), levels = rev(row.names(coef.annual.forb))),
                        x = Estimate,
                        xmin = `l-95% CI`,
                        xmax = `u-95% CI`,
                        fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#B50200") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Annual Forbs")
annual.forb.plot

perennial.grass.traits.NW.model = readRDS("./Results/perennial.grass.imputed.traits.no_woody.rds")

sum = summary(perennial.grass.traits.NW.model)
coef = sum$fixed
row.names(coef) = c("Intercept","Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter",
                    "SRL", "RTD","RMF","DSI","Precipitation")
coef.perennial.grass = coef %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#6089B5"))

perennial.grass.plot = ggplot(data = coef.perennial.grass,
                    aes(y = factor(row.names(coef.perennial.grass), levels = rev(row.names(coef.perennial.grass))),
                        x = Estimate,
                        xmin = `l-95% CI`,
                        xmax = `u-95% CI`,
                        fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#6089B5") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Perennial Grasses")
perennial.grass.plot

perennial.forb.traits.NW.model = readRDS("./Results/perennial.forb.imputed.traits.no_woody.rds")

sum = summary(perennial.forb.traits.NW.model)
coef = sum$fixed
row.names(coef) = c("Intercept","Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter",
                    "SRL", "RTD","RMF","DSI","Precipitation")
coef.perennial.forb = coef %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "black"))

perennial.forb.plot = ggplot(data = coef.perennial.forb,
                              aes(y = factor(row.names(coef.perennial.forb), levels = rev(row.names(coef.perennial.forb))),
                                  x = Estimate,
                                  xmin = `l-95% CI`,
                                  xmax = `u-95% CI`,
                                  fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "black") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Perennial Forbs")
perennial.forb.plot

library(cowplot)

png(filename = "./Plots/coef.plot.png",  width = 11, 
    height = 11, units = "in", res = 600)

plot_grid(All.species.plot,annual.plot,perennial.plot,grass.plot,forb.plot,
                     annual.forb.plot,perennial.forb.plot,perennial.grass.plot)

dev.off()

#### Plotting random effects ####

imputed.traits.NW.model = readRDS("./Results/all.imputed.traits.no_woody.rds")

rf = ranef(imputed.traits.NW.model)
species.rf = as.data.frame(rf$Taxon)
row.names(species.rf) = stringr::str_to_sentence(row.names(species.rf))
site.rf = as.data.frame(rf$site_code)


species.rf.2 = species.rf %>%
  mutate(coef.zero = Q2.5.Intercept/Q97.5.Intercept  < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#769370"))

species.rf.plot = ggplot(data = species.rf.2,
                          aes(y = row.names(species.rf.2),
                              x = Estimate.Intercept,
                              xmin = Q2.5.Intercept,
                              xmax = Q97.5.Intercept,
                              fill = resp.fill)) +
  geom_pointrange(size = 1, shape = 21, color = "#769370") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic()+
  labs(y = "",x = "Mean with 95% CI", title = "All-Species")
species.rf.plot
