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

#### ALL Imputed plots ####

imputed.traits.NW.model = readRDS("./Results/all.imputed.traits.no_woody.rds")

imputed.NW.2 = imputed.NW %>%
  dplyr::select(cover.change,leafN.final:mean.MAP) %>%
  drop_na()

range(imputed.NW.2$leafN.final) # -2.597385  3.235539
range(imputed.NW.2$height.final) # -1.438810  4.140568
range(imputed.NW.2$rootN.final) # -2.110828  3.513945
range(imputed.NW.2$SLA.final) # -2.018426  3.074172
range(imputed.NW.2$root.depth.final) # -0.9946143  5.8703299
range(imputed.NW.2$rootDiam.final) # -1.699918  4.853275
range(imputed.NW.2$SRL.final) # -1.399432  4.192371
range(imputed.NW.2$RTD.final) # -1.825843  3.321296
range(imputed.NW.2$RMF.final) # -2.642873  3.138332
range(imputed.NW.2$mean.DSI) # -2.111273  2.998072
range(imputed.NW.2$mean.MAP) # -1.147113  3.437147

all.effects = conditional_effects(imputed.traits.NW.model)

leafN.effects = all.effects$leafN.final

x.value = c(-2,-1,0,1,2,3)
(x.value*7.114868) + 21.16957

all.leafN.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = leafN.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.final, y = estimate__), color = "#769370", size = 1.5) +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Leaf N", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2.6,3.3),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(6.93, 14.05, 21.17, 28.28, 35.40, 42.51))+
  theme_classic(base_size = 22)
all.leafN.plot

ggsave("./Plots/all.traits.NW.leafN.pdf", height = 7, width = 7)

height.effects = all.effects$height.final

x.value = c(-1,0,1,2,3,4)
(x.value*0.2353906) + 0.3476824

all.height.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = height.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.final, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = height.effects, aes(x = height.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Height", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.5, 4.2),
                     breaks = seq(-1, 4, by = 1),
                     labels = c(0.11, 0.35, 0.58, 0.82, 1.05, 1.29))+
  theme_classic(base_size = 22)
all.height.plot

ggsave("./Plots/all.traits.NW.height.pdf", height = 7, width = 7)

rootN.effects = all.effects$rootN.final

x.value = c(-2,-1,0,1,2,3)
(x.value*4.488773) + 9.855027

all.rootN.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = rootN.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.final, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Root N", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2.2, 3.6),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(0.88, 5.36, 9.86, 14.34, 18.83, 23.32))+
  theme_classic(base_size = 22)
all.rootN.plot

ggsave("./Plots/all.traits.NW.rootN.pdf", height = 7, width = 7)

SLA.effects = all.effects$SLA.final

x.value = c(-2,-1,0,1,2,3)
(x.value*8.581554) + 19.93883

all.SLA.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = SLA.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.final, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "SLA", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2.2, 3.1),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(2.78, 11.36, 19.94, 28.52, 37.11, 45.68))+
  theme_classic(base_size = 22)
all.SLA.plot

ggsave("./Plots/all.traits.NW.SLA.pdf", height = 7, width = 7)

root.depth.effects = all.effects$root.depth.final

x.value = c(-1,0,1,2,3,4,5,6)
(x.value*0.5333377) + 0.5784653

all.root.depth.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = root.depth.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = root.depth.final, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = root.depth.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Rooting Depth", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1, 5.9),
                     breaks = seq(-1,6, by = 1),
                     labels = c(0.05,0.58, 1.11, 1.65, 2.18, 2.71, 3.24,3.79))+
  theme_classic(base_size = 22)
all.root.depth.plot

ggsave("./Plots/all.traits.NW.root.depth.pdf", height = 7, width = 7)

rootDiam.effects = all.effects$rootDiam.final

x.value = c(-1,0,1,2,3,4,5)
(x.value*0.1598629) + 0.3501538

all.rootDiam.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = rootDiam.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = rootDiam.final, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = rootDiam.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Root Diameter", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.7, 4.9),
                     breaks = seq(-1, 5, by = 1),
                     labels = c(0.19, 0.35, 0.51, 0.67, 0.83, 0.99,1.15))+
  theme_classic(base_size = 22)
all.rootDiam.plot

ggsave("./Plots/all.traits.NW.rootDiam.pdf", height = 7, width = 7)

SRL.effects = all.effects$SRL.final

x.value = c(-1,0,1,2,3,4)
(x.value*74.26099) + 106.4204

all.SRL.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = SRL.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.final, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Specific Root Length", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.4, 4.2),
                     breaks = seq(-1, 4, by = 1),
                     labels = c(32.16, 106.42, 180.68, 254.94, 329.20, 403.46))+
  theme_classic(base_size = 22)
all.SRL.plot

ggsave("./Plots/all.traits.NW.SRL.pdf", height = 7, width = 7)

RTD.effects = all.effects$RTD.final

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1159369) + 0.2349392

all.RTD.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = RTD.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.final, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Root Tissue Density", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.9, 3.4),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(0.12, 0.23, 0.35, 0.47, 0.58, 0.70))+
  theme_classic(base_size = 22)
all.RTD.plot

ggsave("./Plots/all.traits.NW.RTD.pdf", height = 7, width = 7)

RMF.effects = all.effects$RMF.final

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1058164) + 0.4007788

all.RMF.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = RMF.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RMF.effects, aes(x = RMF.final, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RMF.effects, aes(x = RMF.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Root Mass Fraction", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2.7, 3.2),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(0.19, 0.29, 0.40, 0.51, 0.61, 0.72))+
  theme_classic(base_size = 22)
all.RMF.plot

ggsave("./Plots/all.traits.NW.RMF.pdf", height = 7, width = 7)

DSI.effects = all.effects$mean.DSI

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1950438) + -0.4604709

all.DSI.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = mean.DSI, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = DSI.effects, aes(x = mean.DSI, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = DSI.effects, aes(x = mean.DSI, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Drought Severity Index", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2.2, 3),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(-0.85, -0.66, -0.46, -0.27, -0.07, 0.12))+
  theme_classic(base_size = 22)
all.DSI.plot

ggsave("./Plots/all.traits.NW.DSI.pdf", height = 7, width = 7)

MAP.effects = all.effects$mean.MAP

x.value = c(-1,0,1,2,3)
(x.value*487.1451) + 691.6107

all.MAP.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = mean.MAP, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = MAP.effects, aes(x = mean.MAP, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = MAP.effects, aes(x = mean.MAP, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Mean Annual Precipitation", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.2, 3.5),
                     breaks = seq(-1, 3, by = 1),
                     labels = c(204.47, 691.61, 1178.17, 1665.90, 2153.05))+
  theme_classic(base_size = 22)
all.MAP.plot

ggsave("./Plots/all.traits.NW.MAP.pdf", height = 7, width = 7)

#### ALL interactions ####

imputed.NW.traits.height.leafN = readRDS("./Results/all.imputed.traits.NW.height.leafN.rds")

height.leafN.effect = conditional_effects(imputed.NW.traits.height.leafN, effects = "height.final:leafN.final")$`height.final:leafN.final`

# only plot high and low values

height.leafN.effect.2 = height.leafN.effect %>%
  filter(effect2__ %in% c(-0.98,1.03))
# height
x.value = c(-1,0,1,2,3,4)
(x.value*0.2353906) + 0.3476824

# leafN
x.value = c(-0.98,1.03)
(x.value*7.114868) + 21.16957


all.height.x.leafN.plot = ggplot(data = height.leafN.effect.2, aes(x = height.final, y = estimate__, group = as.factor(effect2__))) +
  geom_line(aes(color = as.factor(effect2__)), size = 1.5) +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = as.factor(effect2__)), alpha = 0.5,show.legend = FALSE) +  
  geom_line(aes(y = upper__, color = factor(effect2__)), size = 1.2) +
  geom_line(aes(y = lower__, color = factor(effect2__)), size = 1.2) + 
  labs(x = "Height", y = "Percent Cover Change", color = "Leaf N") +
  scale_colour_manual(values = c("black", "#769370"), labels = c("28.50","14.20"))+
  scale_fill_manual(values = c("black", "#769370"))+
  scale_x_continuous(breaks = seq(-1, 4, by = 1),
                     labels = c(0.11, 0.35, 0.58, 0.82, 1.05, 1.29))+
  theme_classic(base_size = 22)
all.height.x.leafN.plot

ggsave("./Plots/all.traits.NW.height.leafN.pdf", height = 7, width = 7)
ggsave("./Plots/all.traits.NW.height.leafN.legend.pdf", height = 7, width = 7)

imputed.NW.traits.leafN.RMF = readRDS("./Results/all.imputed.NW.traits.leafN.RMF.rds")

leafN.RMF.effect = conditional_effects(imputed.NW.traits.leafN.RMF, effects = "leafN.final:RMF.final")$`leafN.final:RMF.final`

# only plot high and low values

leafN.RMF.effect.2 = leafN.RMF.effect %>%
  filter(effect2__ %in% c(-1.02,0.99))
# leafN
x.value = c(-2,-1,0,1,2,3)
(x.value*7.114868) + 21.16957

# RMF
x.value = c(-1.02,0.99)
(x.value*0.1058164) + 0.4007788

all.leafN.x.RMF.plot = ggplot(data = leafN.RMF.effect.2, aes(x = leafN.final, y = estimate__, group = as.factor(effect2__))) +
  geom_line(aes(color = as.factor(effect2__)), size = 1.5) +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = as.factor(effect2__)), alpha = 0.5,show.legend = FALSE) +  
  geom_line(aes(y = upper__, color = factor(effect2__)), size = 1.2) +
  geom_line(aes(y = lower__, color = factor(effect2__)), size = 1.2) + 
  labs(x = "Leaf N", y = "Percent Cover Change", color = "RMF") +
  scale_colour_manual(values = c("black", "#769370"), labels = c("0.51","0.29"))+
  scale_fill_manual(values = c("black", "#769370"))+
  scale_x_continuous(breaks = seq(-2, 3, by = 1),
                     labels = c(6.93, 14.05, 21.17, 28.28, 35.40, 42.51))+
  theme_classic(base_size = 22)
all.leafN.x.RMF.plot

ggsave("./Plots/all.traits.NW.leafN.RMF.pdf", height = 7, width = 7)
ggsave("./Plots/all.traits.NW.leafN.RMF.legend.pdf", height = 7, width = 7)

#### ANNUAL Imputed plots ####

annual.traits.NW.model = readRDS("./Results/annual.imputed.traits.no_woody.rds")

annual.imputed.NW.2 = imputed.NW.annual %>%
  dplyr::select(cover.change,leafN.final:mean.MAP) %>%
  drop_na()

annual.effects = conditional_effects(annual.traits.NW.model)

range(annual.imputed.NW.2$leafN.final) # -2.220848  2.930715
range(annual.imputed.NW.2$height.final) # -1.328275  3.539346
range(annual.imputed.NW.2$rootN.final) # -2.099977  3.608290
range(annual.imputed.NW.2$SLA.final) # -2.169165  2.342324
range(annual.imputed.NW.2$root.depth.final) # -1.132666  3.075555
range(annual.imputed.NW.2$rootDiam.final) # -1.760667  2.969117
range(annual.imputed.NW.2$SRL.final) # -1.405576  3.973989
range(annual.imputed.NW.2$RTD.final) # -1.613838  3.181998
range(annual.imputed.NW.2$RMF.final) # -2.625359  2.581537
range(annual.imputed.NW.2$mean.DSI) # -1.752847  2.579177
range(annual.imputed.NW.2$mean.MAP) # -1.176711  2.739820

leafN.effects = annual.effects$leafN.final

attr(annual.imputed.NW.2$leafN.final, "scaled:scale")
attr(annual.imputed.NW.2$leafN.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*7.384253) + 22.54886

annual.leafN.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = leafN.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.final, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Leaf N", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2.3,3.0),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(7.78, 15.16, 22.55, 29.93, 37.31, 44.70))+
  theme_classic(base_size = 22)
annual.leafN.plot

ggsave("./Plots/annual.traits.NW.leafN.pdf", height = 7, width = 7)

height.effects = annual.effects$height.final

attr(annual.imputed.NW.2$height.final, "scaled:scale")
attr(annual.imputed.NW.2$height.final, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*0.2308751) + 0.3156657

annual.height.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = height.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.final, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = height.effects, aes(x = height.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Height", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.4, 3.6),
                     breaks = seq(-1, 4, by = 1),
                     labels = c(0.08, 0.32, 0.55, 0.78, 1.01, 1.24))+
  theme_classic(base_size = 22)
annual.height.plot

ggsave("./Plots/annual.traits.NW.height.pdf", height = 7, width = 7)

rootN.effects = annual.effects$rootN.final

attr(annual.imputed.NW.2$rootN.final, "scaled:scale")
attr(annual.imputed.NW.2$rootN.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*4.33202) + 9.997142

annual.rootN.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = rootN.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.final, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Root N", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2.1, 3.7),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(1.33, 5.67, 9.99, 14.33, 18.66, 22.99))+
  theme_classic(base_size = 22)
annual.rootN.plot

ggsave("./Plots/annual.traits.NW.rootN.pdf", height = 7, width = 7)

SLA.effects = annual.effects$SLA.final

attr(annual.imputed.NW.2$SLA.final, "scaled:scale")
attr(annual.imputed.NW.2$SLA.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*9.312314) + 22.81754

annual.SLA.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = SLA.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.final, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "SLA", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2.2, 2.4),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(4.19, 13.51, 22.82, 32.13, 41.44, 50.75))+
  theme_classic(base_size = 22)
annual.SLA.plot

ggsave("./Plots/annual.traits.NW.SLA.pdf", height = 7, width = 7)

root.depth.effects = annual.effects$root.depth.final

attr(annual.imputed.NW.2$root.depth.final, "scaled:scale")
attr(annual.imputed.NW.2$root.depth.final, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*0.4346265) + 0.5402869

annual.root.depth.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = root.depth.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = root.depth.final, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = root.depth.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Rooting Depth", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.2, 3.1),
                     breaks = seq(-1,4, by = 1),
                     labels = c(0.11,0.54, 0.97, 1.41, 1.84, 2.28))+
  theme_classic(base_size = 22)
annual.root.depth.plot

ggsave("./Plots/annual.traits.NW.root.depth.pdf", height = 7, width = 7)

rootDiam.effects = annual.effects$rootDiam.final

attr(annual.imputed.NW.2$rootDiam.final, "scaled:scale")
attr(annual.imputed.NW.2$rootDiam.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1367222) + 0.3407223

annual.rootDiam.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = rootDiam.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = rootDiam.final, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = rootDiam.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Root Diameter", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.8, 3.0),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(0.07, 0.20, 0.34, 0.48, 0.61, 0.75))+
  theme_classic(base_size = 22)
annual.rootDiam.plot

ggsave("./Plots/annual.traits.NW.rootDiam.pdf", height = 7, width = 7)

SRL.effects = annual.effects$SRL.final

attr(annual.imputed.NW.2$SRL.final, "scaled:scale")
attr(annual.imputed.NW.2$SRL.final, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*73.49115) + 125.697

annual.SRL.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = SRL.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.final, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Specific Root Length", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.5, 4.0),
                     breaks = seq(-1, 4, by = 1),
                     labels = c(52.21, 125.69, 199.19, 272.68, 346.17, 419.66))+
  theme_classic(base_size = 22)
annual.SRL.plot

ggsave("./Plots/annual.traits.NW.SRL.pdf", height = 7, width = 7)

RTD.effects = annual.effects$RTD.final

attr(annual.imputed.NW.2$RTD.final, "scaled:scale")
attr(annual.imputed.NW.2$RTD.final, "scaled:center")

x.value = c(-1,0,1,2,3)
(x.value*0.1111574) + 0.2026466

annual.RTD.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = RTD.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.final, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Root Tissue Density", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.7, 3.2),
                     breaks = seq(-1, 3, by = 1),
                     labels = c(0.09, 0.21, 0.32, 0.42, 0.54))+
  theme_classic(base_size = 22)
annual.RTD.plot

ggsave("./Plots/annual.traits.NW.RTD.pdf", height = 7, width = 7)

RMF.effects = annual.effects$RMF.final

attr(annual.imputed.NW.2$RMF.final, "scaled:scale")
attr(annual.imputed.NW.2$RMF.final, "scaled:center")

x.value = c(-3,-2,-1,0,1,2,3)
(x.value*0.1037188) + 0.3934184

annual.RMF.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = RMF.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RMF.effects, aes(x = RMF.final, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RMF.effects, aes(x = RMF.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Root Mass Fraction", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2.7, 2.6),
                     breaks = seq(-3, 3, by = 1),
                     labels = c(0.08, 0.19, 0.29, 0.39, 0.49, 0.61, 0.71))+
  theme_classic(base_size = 22)
annual.RMF.plot

ggsave("./Plots/annual.traits.NW.RMF.pdf", height = 7, width = 7)

DSI.effects = annual.effects$mean.DSI

attr(annual.imputed.NW.2$mean.DSI, "scaled:scale")
attr(annual.imputed.NW.2$mean.DSI, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.2300417) + -0.4690337

annual.DSI.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = mean.DSI, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = DSI.effects, aes(x = mean.DSI, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = DSI.effects, aes(x = mean.DSI, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Drought Severity Index", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.8,2.6),
                     breaks = seq(-1, 3, by = 1),
                     labels = c(-0.69, -0.46, -0.24, -0.01, 0.22))+
  theme_classic(base_size = 22)
annual.DSI.plot

ggsave("./Plots/annual.traits.NW.DSI.pdf", height = 7, width = 7)

MAP.effects = annual.effects$mean.MAP

attr(annual.imputed.NW.2$mean.MAP, "scaled:scale")
attr(annual.imputed.NW.2$mean.MAP, "scaled:center")

x.value = c(-1,0,1,2,3)
(x.value*293.1676) + 477.7734

annual.MAP.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = mean.MAP, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = MAP.effects, aes(x = mean.MAP, y = estimate__), color = "#979461", size = 1.5) +  
  geom_ribbon(data = MAP.effects, aes(x = mean.MAP, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Mean Annual Precipitation", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.2, 2.8),
                     breaks = seq(-1, 2, by = 1),
                     labels = c(184.61, 477.77, 770.94, 1064.11))+
  theme_classic(base_size = 22)
annual.MAP.plot

ggsave("./Plots/annual.traits.NW.MAP.pdf", height = 7, width = 7)

#### ANNUAL interactions ####

annual.imputed.NW.traits.leafN.RMF = readRDS("./Results/annual.imputed.NW.traits.leafN.RMF.rds")

leafN.RMF.effect = conditional_effects(annual.imputed.NW.traits.leafN.RMF, effects = "leafN.final:RMF.final")$`leafN.final:RMF.final`

# only plot high and low values

leafN.RMF.effect.2 = leafN.RMF.effect %>%
  filter(effect2__ %in% c(-1.05,0.97))

# leafN
x.value = c(-2,-1,0,1,2,3)
(x.value*7.384253) + 22.54886

# RMF
x.value = c(-1.05,0.97)
(x.value*0.1037188) + 0.3934184

annual.leafN.x.RMF.plot = ggplot(data = leafN.RMF.effect.2, aes(x = leafN.final, y = estimate__, group = as.factor(effect2__))) +
  geom_line(aes(color = as.factor(effect2__)), size = 1.5,show.legend = FALSE) +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = as.factor(effect2__)), alpha = 0.5,show.legend = FALSE) +  
  geom_line(aes(y = upper__, color = factor(effect2__)), size = 1.2,show.legend = FALSE) +
  geom_line(aes(y = lower__, color = factor(effect2__)), size = 1.2,show.legend = FALSE) + 
  labs(x = "Leaf N", y = "Percent Cover Change", color = "RMF") +
  scale_colour_manual(values = c("black", "#979461"), labels = c("0.49","0.28"))+
  scale_fill_manual(values = c("black", "#979461"))+
  scale_x_continuous(breaks = seq(-2, 3, by = 1),
                     labels = c(7.78, 15.16, 22.55, 29.93, 37.32, 44.70))+
  theme_classic(base_size = 22)
annual.leafN.x.RMF.plot

ggsave("./Plots/annual.traits.NW.leafN.RMF.pdf", height = 7, width = 7)
ggsave("./Plots/annual.traits.NW.leafN.RMF.legend.pdf", height = 7, width = 7)

#### PERENNIAL Imputed plots ####

perennial.traits.NW.model = readRDS("./Results/perennial.imputed.traits.no_woody.rds")

perennial.imputed.NW.2 = imputed.NW.perennial %>%
  dplyr::select(cover.change,leafN.final:mean.MAP) %>%
  drop_na()

perennial.effects = conditional_effects(perennial.traits.NW.model)

range(perennial.imputed.NW.2$leafN.final) # -2.609822  3.130349
range(perennial.imputed.NW.2$height.final) # -1.448808  4.056365
range(perennial.imputed.NW.2$rootN.final) # -2.070563  3.386362
range(perennial.imputed.NW.2$SLA.final) # -1.933604  3.459736
range(perennial.imputed.NW.2$root.depth.final) # -0.9029594  5.4653640
range(perennial.imputed.NW.2$rootDiam.final) # -1.654370  4.644744
range(perennial.imputed.NW.2$SRL.final) # -1.316425  3.643188
range(perennial.imputed.NW.2$RTD.final) # -1.954990  3.215066
range(perennial.imputed.NW.2$RMF.final) # -2.352661  3.120367
range(perennial.imputed.NW.2$mean.DSI) # -2.346473  3.254426
range(perennial.imputed.NW.2$mean.MAP) # -1.110111  3.050369

leafN.effects = perennial.effects$leafN.final

attr(perennial.imputed.NW.2$leafN.final, "scaled:scale")
attr(perennial.imputed.NW.2$leafN.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*6.812995) + 20.47022

perennial.leafN.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = leafN.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.final, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Leaf N", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2.7,3.2),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(6.84, 13.66, 20.47, 27.28, 34.10, 40.91))+
  theme_classic(base_size = 22)
perennial.leafN.plot

ggsave("./Plots/perennial.traits.NW.leafN.pdf", height = 7, width = 7)

height.effects = perennial.effects$height.final

attr(perennial.imputed.NW.2$height.final, "scaled:scale")
attr(perennial.imputed.NW.2$height.final, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*0.2365653) + 0.3627378

perennial.height.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = height.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.final, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = height.effects, aes(x = height.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Height", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.5, 4.1),
                     breaks = seq(-1, 4, by = 1),
                     labels = c(0.13, 0.36, 0.59, 0.84, 1.07, 1.31))+
  theme_classic(base_size = 22)
perennial.height.plot

ggsave("./Plots/perennial.traits.NW.height.pdf", height = 7, width = 7)

rootN.effects = perennial.effects$rootN.final

attr(perennial.imputed.NW.2$rootN.final, "scaled:scale")
attr(perennial.imputed.NW.2$rootN.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*4.532848) + 9.765547

perennial.rootN.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = rootN.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.final, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Root N", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2.1, 3.4),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(0.70, 5.23, 9.77, 14.29, 18.83, 23.36))+
  theme_classic(base_size = 22)
perennial.rootN.plot

ggsave("./Plots/perennial.traits.NW.rootN.pdf", height = 7, width = 7)

SLA.effects = perennial.effects$SLA.final

attr(perennial.imputed.NW.2$SLA.final, "scaled:scale")
attr(perennial.imputed.NW.2$SLA.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*7.966249) + 18.75888

perennial.SLA.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = SLA.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.final, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "SLA", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2.0, 3.5),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(2.83, 10.79, 18.75, 26.73, 34.69, 42.66))+
  theme_classic(base_size = 22)
perennial.SLA.plot

ggsave("./Plots/perennial.traits.NW.SLA.pdf", height = 7, width = 7)

root.depth.effects = perennial.effects$root.depth.final

attr(perennial.imputed.NW.2$root.depth.final, "scaled:scale")
attr(perennial.imputed.NW.2$root.depth.final, "scaled:center")

x.value = c(-1,0,1,2,3,4,5)
(x.value*0.5706892) + 0.5903092

perennial.root.depth.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = root.depth.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = root.depth.final, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = root.depth.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Rooting Depth", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1, 5.5),
                     breaks = seq(-1,5, by = 1),
                     labels = c(0.02,0.59, 1.16, 1.73, 2.30, 2.87, 3.44))+
  theme_classic(base_size = 22)
perennial.root.depth.plot

ggsave("./Plots/perennial.traits.NW.root.depth.pdf", height = 7, width = 7)

rootDiam.effects = perennial.effects$rootDiam.final

attr(perennial.imputed.NW.2$rootDiam.final, "scaled:scale")
attr(perennial.imputed.NW.2$rootDiam.final, "scaled:center")

x.value = c(-1,0,1,2,3,4,5)
(x.value*0.166311) + 0.3535401

perennial.rootDiam.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = rootDiam.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = rootDiam.final, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = rootDiam.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Root Diameter", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.7, 4.7),
                     breaks = seq(-1, 4, by = 1),
                     labels = c(0.19, 0.35, 0.52, 0.69, 0.85, 1.02))+
  theme_classic(base_size = 22)
perennial.rootDiam.plot

ggsave("./Plots/perennial.traits.NW.rootDiam.pdf", height = 7, width = 7)

SRL.effects = perennial.effects$SRL.final

attr(perennial.imputed.NW.2$SRL.final, "scaled:scale")
attr(perennial.imputed.NW.2$SRL.final, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*73.08878) + 98.7131

perennial.SRL.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = SRL.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.final, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Specific Root Length", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.4, 3.7),
                     breaks = seq(-1, 4, by = 1),
                     labels = c(25.62, 98.71, 171.80, 244.89, 317.98, 391.07))+
  theme_classic(base_size = 22)
perennial.SRL.plot

ggsave("./Plots/perennial.traits.NW.SRL.pdf", height = 7, width = 7)

RTD.effects = perennial.effects$RTD.final

attr(perennial.imputed.NW.2$RTD.final, "scaled:scale")
attr(perennial.imputed.NW.2$RTD.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1151264) + 0.249861

perennial.RTD.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = RTD.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.final, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Root Tissue Density", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2.0, 3.3),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(0.02, 0.13, 0.25, 0.36, 0.48, 0.59))+
  theme_classic(base_size = 22)
perennial.RTD.plot

ggsave("./Plots/perennial.traits.NW.RTD.pdf", height = 7, width = 7)

RMF.effects = perennial.effects$RMF.final

attr(perennial.imputed.NW.2$RMF.final, "scaled:scale")
attr(perennial.imputed.NW.2$RMF.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1053099) + 0.4042604

perennial.RMF.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = RMF.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RMF.effects, aes(x = RMF.final, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RMF.effects, aes(x = RMF.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Root Mass Fraction", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2.4, 3.2),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(0.19, 0.29, 0.41, 0.51, 0.61, 0.72))+
  theme_classic(base_size = 22)
perennial.RMF.plot

ggsave("./Plots/perennial.traits.NW.RMF.pdf", height = 7, width = 7)

DSI.effects = perennial.effects$mean.DSI

attr(perennial.imputed.NW.2$mean.DSI, "scaled:scale")
attr(perennial.imputed.NW.2$mean.DSI, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1779261) + -0.4547629

perennial.DSI.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = mean.DSI, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = DSI.effects, aes(x = mean.DSI, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = DSI.effects, aes(x = mean.DSI, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Drought Severity Index", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2.4, 3.3),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(-0.81, -0.63, -0.45, -0.28, -0.09, 0.08))+
  theme_classic(base_size = 22)
perennial.DSI.plot

ggsave("./Plots/perennial.traits.NW.DSI.pdf", height = 7, width = 7)

MAP.effects = perennial.effects$mean.MAP

attr(perennial.imputed.NW.2$mean.MAP, "scaled:scale")
attr(perennial.imputed.NW.2$mean.MAP, "scaled:center")

x.value = c(-1,0,1,2,3)
(x.value*522.5359) + 772.0729

perennial.MAP.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = mean.MAP, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = MAP.effects, aes(x = mean.MAP, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = MAP.effects, aes(x = mean.MAP, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Mean perennial Precipitation", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.2, 3.1),
                     breaks = seq(-1, 3, by = 1),
                     labels = c(249.54, 772.07, 1294.61, 1817.14, 2339.68))+
  theme_classic(base_size = 22)
perennial.MAP.plot

ggsave("./Plots/perennial.traits.NW.MAP.pdf", height = 7, width = 7)


#### PERENNIAL interactions ####

perennial.imputed.NW.traits.height.leafN = readRDS("./Results/perennial.imputed.traits.NW.height.leafN.rds")

height.leafN.effect = conditional_effects(perennial.imputed.NW.traits.height.leafN, effects = "height.final:leafN.final")$`height.final:leafN.final`

# only plot high and low values

height.leafN.effect.2 = height.leafN.effect %>%
  filter(effect2__ %in% c(-0.97,1.02))

# height
x.value = c(-1,0,1,2,3,4)
(x.value*0.2365653) + 0.3627378

# leafN
x.value = c(-0.97,1.02)
(x.value*6.812995) + 20.47022

perennial.height.x.leafN.plot = ggplot(data = height.leafN.effect.2, aes(x = height.final, y = estimate__, group = as.factor(effect2__))) +
  geom_line(aes(color = as.factor(effect2__)), size = 1.5,show.legend = FALSE) +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = as.factor(effect2__)), alpha = 0.5,show.legend = FALSE) +  
  geom_line(aes(y = upper__, color = factor(effect2__)), size = 1.2,show.legend = FALSE) +
  geom_line(aes(y = lower__, color = factor(effect2__)), size = 1.2,show.legend = FALSE) + 
  labs(x = "Height", y = "Percent Cover Change", color = "Leaf N") +
  scale_colour_manual(values = c("black", "#F1C646"), labels = c("27.41","13.86"))+
  scale_fill_manual(values = c("black", "#F1C646"))+
  scale_x_continuous(breaks = seq(-1, 4, by = 1),
                     labels = c(0.13, 0.36, 0.59, 0.84, 1.07, 1.31))+
  theme_classic(base_size = 22)
perennial.height.x.leafN.plot

ggsave("./Plots/perennial.traits.NW.height.leafN.pdf", height = 7, width = 7)
ggsave("./Plots/perennial.traits.NW.height.leafN.legend.pdf", height = 7, width = 7)

#### FORB Imputed plots ####

forb.traits.NW.model = readRDS("./Results/forb.imputed.traits.no_woody.rds")

forb.imputed.NW.2 = imputed.NW.forb %>%
  dplyr::select(cover.change,leafN.final:mean.MAP) %>%
  drop_na()

forb.effects = conditional_effects(forb.traits.NW.model)

range(forb.imputed.NW.2$leafN.final) # -2.964007  3.517241
range(forb.imputed.NW.2$height.final) # -1.303944  4.176462
range(forb.imputed.NW.2$rootN.final) # -2.319661  3.604814
range(forb.imputed.NW.2$SLA.final) # -2.096015  3.002170
range(forb.imputed.NW.2$root.depth.final) # -0.9521214  5.0185553
range(forb.imputed.NW.2$rootDiam.final) # -1.551108  4.827648
range(forb.imputed.NW.2$SRL.final) # -1.429975  4.218829
range(forb.imputed.NW.2$RTD.final) # -1.722628  2.986269
range(forb.imputed.NW.2$RMF.final) # -2.363349  2.835340
range(forb.imputed.NW.2$mean.DSI) # -2.091669  3.020008
range(forb.imputed.NW.2$mean.MAP) # -1.171155  3.810773

leafN.effects = forb.effects$leafN.final

attr(forb.imputed.NW.2$leafN.final, "scaled:scale")
attr(forb.imputed.NW.2$leafN.final, "scaled:center")

x.value = c(-3,-2,-1,0,1,2,3)
(x.value*6.403163) + 21.66854

forb.leafN.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = leafN.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.final, y = estimate__), color = "#F17236", size = 1.5) +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Leaf N", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-3.0,3.6),
                     breaks = seq(-3, 3, by = 1),
                     labels = c(2.46, 8.86,15.27, 21.67, 28.07, 34.47, 40.88))+
  theme_classic(base_size = 22)
forb.leafN.plot

ggsave("./Plots/forb.traits.NW.leafN.pdf", height = 7, width = 7)

height.effects = forb.effects$height.final

attr(forb.imputed.NW.2$height.final, "scaled:scale")
attr(forb.imputed.NW.2$height.final, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*0.2396416) + 0.3214792

forb.height.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = height.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.final, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = height.effects, aes(x = height.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Height", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.4, 4.2),
                     breaks = seq(-1, 4, by = 1),
                     labels = c(0.08, 0.32, 0.56, 0.81, 1.04, 1.28))+
  theme_classic(base_size = 22)
forb.height.plot

ggsave("./Plots/forb.traits.NW.height.pdf", height = 7, width = 7)

rootN.effects = forb.effects$rootN.final

attr(forb.imputed.NW.2$rootN.final, "scaled:scale")
attr(forb.imputed.NW.2$rootN.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*4.173927) + 10.5821

forb.rootN.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = rootN.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.final, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Root N", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2.4, 3.7),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(2.23, 6.41, 10.58, 14.76, 18.93, 23.11))+
  theme_classic(base_size = 22)
forb.rootN.plot

ggsave("./Plots/forb.traits.NW.rootN.pdf", height = 7, width = 7)

SLA.effects = forb.effects$SLA.final

attr(forb.imputed.NW.2$SLA.final, "scaled:scale")
attr(forb.imputed.NW.2$SLA.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*8.362926) + 20.14641

forb.SLA.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = SLA.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.final, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "SLA", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2.1, 3.1),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(3.42, 11.78, 20.15, 28.51, 36.87, 45.23))+
  theme_classic(base_size = 22)
forb.SLA.plot

ggsave("./Plots/forb.traits.NW.SLA.pdf", height = 7, width = 7)

root.depth.effects = forb.effects$root.depth.final

attr(forb.imputed.NW.2$root.depth.final, "scaled:scale")
attr(forb.imputed.NW.2$root.depth.final, "scaled:center")

x.value = c(-1,0,1,2,3,4,5)
(x.value*0.6095345) + 0.6503508

forb.root.depth.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = root.depth.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = root.depth.final, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = root.depth.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Rooting Depth", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1, 5.1),
                     breaks = seq(-1,5, by = 1),
                     labels = c(0.04,0.65, 1.26, 1.87, 2.48, 3.09, 3.69))+
  theme_classic(base_size = 22)
forb.root.depth.plot

ggsave("./Plots/forb.traits.NW.root.depth.pdf", height = 7, width = 7)

rootDiam.effects = forb.effects$rootDiam.final

attr(forb.imputed.NW.2$rootDiam.final, "scaled:scale")
attr(forb.imputed.NW.2$rootDiam.final, "scaled:center")

x.value = c(-1,0,1,2,3,4,5)
(x.value*0.1590925) + 0.3579696

forb.rootDiam.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = rootDiam.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = rootDiam.final, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = rootDiam.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Root Diameter", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.6, 4.9),
                     breaks = seq(-1, 5, by = 1),
                     labels = c(0.19, 0.36, 0.52, 0.68, 0.84, 0.99,1.15))+
  theme_classic(base_size = 22)
forb.rootDiam.plot

ggsave("./Plots/forb.traits.NW.rootDiam.pdf", height = 7, width = 7)

SRL.effects = forb.effects$SRL.final

attr(forb.imputed.NW.2$SRL.final, "scaled:scale")
attr(forb.imputed.NW.2$SRL.final, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*71.25811) + 105.5323

forb.SRL.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = SRL.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.final, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Specific Root Length", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.5, 4.3),
                     breaks = seq(-1, 4, by = 1),
                     labels = c(34.27, 105.53, 176.79, 248.05, 319.31, 390.56))+
  theme_classic(base_size = 22)
forb.SRL.plot

ggsave("./Plots/forb.traits.NW.SRL.pdf", height = 7, width = 7)

RTD.effects = forb.effects$RTD.final

attr(forb.imputed.NW.2$RTD.final, "scaled:scale")
attr(forb.imputed.NW.2$RTD.final, "scaled:center")

x.value = c(-1,0,1,2,3)
(x.value*0.1184977) + 0.227384

forb.RTD.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = RTD.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.final, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Root Tissue Density", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.8, 3.0),
                     breaks = seq(-1, 3, by = 1),
                     labels = c(0.11, 0.23, 0.35, 0.46, 0.58))+
  theme_classic(base_size = 22)
forb.RTD.plot

ggsave("./Plots/forb.traits.NW.RTD.pdf", height = 7, width = 7)

RMF.effects = forb.effects$RMF.final

attr(forb.imputed.NW.2$RMF.final, "scaled:scale")
attr(forb.imputed.NW.2$RMF.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1176732) + 0.3992223

forb.RMF.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = RMF.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RMF.effects, aes(x = RMF.final, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RMF.effects, aes(x = RMF.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Root Mass Fraction", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2.4, 2.9),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(0.16, 0.28, 0.39, 0.52, 0.31, 0.75))+
  theme_classic(base_size = 22)
forb.RMF.plot

ggsave("./Plots/forb.traits.NW.RMF.pdf", height = 7, width = 7)

DSI.effects = forb.effects$mean.DSI

attr(forb.imputed.NW.2$mean.DSI, "scaled:scale")
attr(forb.imputed.NW.2$mean.DSI, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1949549) + -0.4644806

forb.DSI.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = mean.DSI, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = DSI.effects, aes(x = mean.DSI, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = DSI.effects, aes(x = mean.DSI, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Drought Severity Index", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2.1, 3.1),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(-0.85, -0.66, -0.46, -0.27, -0.07, 0.12))+
  theme_classic(base_size = 22)
forb.DSI.plot

ggsave("./Plots/forb.traits.NW.DSI.pdf", height = 7, width = 7)

MAP.effects = forb.effects$mean.MAP

attr(forb.imputed.NW.2$mean.MAP, "scaled:scale")
attr(forb.imputed.NW.2$mean.MAP, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*448.2603) + 657.7821

forb.MAP.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = mean.MAP, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = MAP.effects, aes(x = mean.MAP, y = estimate__), color = "#F17236", size = 1.5) +  
  geom_ribbon(data = MAP.effects, aes(x = mean.MAP, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Mean forb Precipitation", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.2, 3.9),
                     breaks = seq(-1, 3, by = 1),
                     labels = c(209.52, 657.78, 1106.04, 1554.31, 2022.56))+
  theme_classic(base_size = 22)
forb.MAP.plot

ggsave("./Plots/forb.traits.NW.MAP.pdf", height = 7, width = 7)

#### FORB interactions ####

forb.imputed.NW.traits.RTD.SRL = readRDS("./Results/forb.imputed.traits.NW.RTD.SRL.rds")

RTD.SRL.effect = conditional_effects(forb.imputed.NW.traits.RTD.SRL, effects = "RTD.final:SRL.final")$`RTD.final:SRL.final`

# only plot high and low values

RTD.SRL.effect.2 = RTD.SRL.effect %>%
  filter(effect2__ %in% c(-1,1.05))

# RTD
x.value = c(-1,0,1,2,3)
(x.value*0.1184977) + 0.227384

# SRL
x.value = c(-1,1.05)
(x.value*71.25811) + 105.5323

forb.RTD.x.SRL.plot = ggplot(data = RTD.SRL.effect.2, aes(x = RTD.final, y = estimate__, group = as.factor(effect2__))) +
  geom_line(aes(color = as.factor(effect2__)), size = 1.5) +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = as.factor(effect2__)), alpha = 0.5,show.legend = FALSE) +  
  geom_line(aes(y = upper__, color = factor(effect2__)), size = 1.2) +
  geom_line(aes(y = lower__, color = factor(effect2__)), size = 1.2) + 
  labs(x = "RTD", y = "Percent Cover Change", color = "SRL") +
  scale_colour_manual(values = c("black", "#F17236"), labels = c("180.35","34.27"))+
  scale_fill_manual(values = c("black", "#F17236"))+
  scale_x_continuous(breaks = seq(-1, 3, by = 1),
                     labels = c(0.11, 0.23, 0.35, 0.47, 0.58))+
  theme_classic(base_size = 22)
forb.RTD.x.SRL.plot

ggsave("./Plots/forb.traits.NW.RTD.SRL.pdf", height = 7, width = 7)
ggsave("./Plots/forb.traits.NW.RTD.SRL.legend.pdf", height = 7, width = 7)
#### GRASS Imputed plots ####

grass.traits.NW.model = readRDS("./Results/grass.imputed.traits.no_woody.rds")

grass.imputed.NW.2 = imputed.NW.grass %>%
  dplyr::select(cover.change,leafN.final:mean.MAP) %>%
  drop_na()

grass.effects = conditional_effects(grass.traits.NW.model)

range(grass.imputed.NW.2$leafN.final) # -1.937329  3.576683
range(grass.imputed.NW.2$height.final) # -1.605196  3.534036
range(grass.imputed.NW.2$rootN.final) # -2.658593  4.697596
range(grass.imputed.NW.2$SLA.final) # -1.772276  2.716809
range(grass.imputed.NW.2$root.depth.final) # -1.029564  6.072718
range(grass.imputed.NW.2$rootDiam.final) # -1.456175  3.228266
range(grass.imputed.NW.2$SRL.final) # -1.261202  3.673079
range(grass.imputed.NW.2$RTD.final) # -2.187709  3.297452
range(grass.imputed.NW.2$RMF.final) # -2.216468  2.901228
range(grass.imputed.NW.2$mean.DSI) # -2.104957  2.858311
range(grass.imputed.NW.2$mean.MAP) # -1.142534  3.915008

leafN.effects = grass.effects$leafN.final

attr(grass.imputed.NW.2$leafN.final, "scaled:scale")
attr(grass.imputed.NW.2$leafN.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*6.281442) + 18.41922

grass.leafN.plot = ggplot() +
  geom_point(data = grass.imputed.NW.2, aes(x = leafN.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.final, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Leaf N", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2,3.6),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(5.86, 12.14,18.42, 24.71, 30.98, 37.26))+
  theme_classic(base_size = 22)
grass.leafN.plot

ggsave("./Plots/grass.traits.NW.leafN.pdf", height = 7, width = 7)

height.effects = grass.effects$height.final

attr(grass.imputed.NW.2$height.final, "scaled:scale")
attr(grass.imputed.NW.2$height.final, "scaled:center")

x.value = c(-1,0,1,2,3)
(x.value*0.2266876) + 0.3988779

grass.height.plot = ggplot() +
  geom_point(data = grass.imputed.NW.2, aes(x = height.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.final, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = height.effects, aes(x = height.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Height", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.7, 3.6),
                     breaks = seq(-1, 3, by = 1),
                     labels = c(0.17, 0.39, 0.63, 0.85, 1.08))+
  theme_classic(base_size = 22)
grass.height.plot

ggsave("./Plots/grass.traits.NW.height.pdf", height = 7, width = 7)

rootN.effects = grass.effects$rootN.final

attr(grass.imputed.NW.2$rootN.final, "scaled:scale")
attr(grass.imputed.NW.2$rootN.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3,4)
(x.value*2.712129) + 7.59045

grass.rootN.plot = ggplot() +
  geom_point(data = grass.imputed.NW.2, aes(x = rootN.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.final, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Root N", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2.7, 4.7),
                     breaks = seq(-2, 4, by = 1),
                     labels = c(2.17, 4.88, 7.59, 10.31, 13.01, 15.73, 18.44))+
  theme_classic(base_size = 22)
grass.rootN.plot

ggsave("./Plots/grass.traits.NW.rootN.pdf", height = 7, width = 7)

SLA.effects = grass.effects$SLA.final

attr(grass.imputed.NW.2$SLA.final, "scaled:scale")
attr(grass.imputed.NW.2$SLA.final, "scaled:center")

x.value = c(-1,0,1,2,3)
(x.value*9.570923) + 20.31763

grass.SLA.plot = ggplot() +
  geom_point(data = grass.imputed.NW.2, aes(x = SLA.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.final, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "SLA", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.8, 2.8),
                     breaks = seq(-1, 2, by = 1),
                     labels = c(10.75, 20.32, 29.89, 39.46))+
  theme_classic(base_size = 22)
grass.SLA.plot

ggsave("./Plots/grass.traits.NW.SLA.pdf", height = 7, width = 7)

root.depth.effects = grass.effects$root.depth.final

attr(grass.imputed.NW.2$root.depth.final, "scaled:scale")
attr(grass.imputed.NW.2$root.depth.final, "scaled:center")

x.value = c(-1,0,1,2,3,4,5,6)
(x.value*0.3671355) + 0.4529895

grass.root.depth.plot = ggplot() +
  geom_point(data = grass.imputed.NW.2, aes(x = root.depth.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = root.depth.final, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = root.depth.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Rooting Depth", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.1, 6.1),
                     breaks = seq(-1,6, by = 1),
                     labels = c(0.08,0.45, 0.82, 1.19, 1.55, 1.92, 2.29, 2.66))+
  theme_classic(base_size = 22)
grass.root.depth.plot

ggsave("./Plots/grass.traits.NW.root.depth.pdf", height = 7, width = 7)

rootDiam.effects = grass.effects$rootDiam.final

attr(grass.imputed.NW.2$rootDiam.final, "scaled:scale")
attr(grass.imputed.NW.2$rootDiam.final, "scaled:center")

x.value = c(-1,0,1,2,3)
(x.value*0.1707781) + 0.3486828

grass.rootDiam.plot = ggplot() +
  geom_point(data = grass.imputed.NW.2, aes(x = rootDiam.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = rootDiam.final, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = rootDiam.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Root Diameter", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.5, 3.3),
                     breaks = seq(-1, 3, by = 1),
                     labels = c(0.18, 0.35, 0.52, 0.69, 0.86))+
  theme_classic(base_size = 22)
grass.rootDiam.plot

ggsave("./Plots/grass.traits.NW.rootDiam.pdf", height = 7, width = 7)

SRL.effects = grass.effects$SRL.final

attr(grass.imputed.NW.2$SRL.final, "scaled:scale")
attr(grass.imputed.NW.2$SRL.final, "scaled:center")

x.value = c(-1,0,1,2,3)
(x.value*82.36672) + 115.2106

grass.SRL.plot = ggplot() +
  geom_point(data = grass.imputed.NW.2, aes(x = SRL.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.final, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Specific Root Length", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.3, 3.7),
                     breaks = seq(-1, 3, by = 1),
                     labels = c(32.84, 115.21, 197.58, 279.94, 362.31))+
  theme_classic(base_size = 22)
grass.SRL.plot

ggsave("./Plots/grass.traits.NW.SRL.pdf", height = 7, width = 7)

RTD.effects = grass.effects$RTD.final

attr(grass.imputed.NW.2$RTD.final, "scaled:scale")
attr(grass.imputed.NW.2$RTD.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1085128) + 0.2621843

grass.RTD.plot = ggplot() +
  geom_point(data = grass.imputed.NW.2, aes(x = RTD.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.final, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Root Tissue Density", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2.2, 3.3),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(0.05, 0.15, 0.26, 0.37, 0.48,0.59))+
  theme_classic(base_size = 22)
grass.RTD.plot

ggsave("./Plots/grass.traits.NW.RTD.pdf", height = 7, width = 7)

RMF.effects = grass.effects$RMF.final

attr(grass.imputed.NW.2$RMF.final, "scaled:scale")
attr(grass.imputed.NW.2$RMF.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.07914561) + 0.4233237

grass.RMF.plot = ggplot() +
  geom_point(data = grass.imputed.NW.2, aes(x = RMF.final, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RMF.effects, aes(x = RMF.final, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RMF.effects, aes(x = RMF.final, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Root Mass Fraction", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2.3, 3.0),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(0.26, 0.34, 0.42, 0.51, 0.58, 0.66))+
  theme_classic(base_size = 22)
grass.RMF.plot

ggsave("./Plots/grass.traits.NW.RMF.pdf", height = 7, width = 7)

DSI.effects = grass.effects$mean.DSI

attr(grass.imputed.NW.2$mean.DSI, "scaled:scale")
attr(grass.imputed.NW.2$mean.DSI, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.2007843) + -0.4496194

grass.DSI.plot = ggplot() +
  geom_point(data = grass.imputed.NW.2, aes(x = mean.DSI, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = DSI.effects, aes(x = mean.DSI, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = DSI.effects, aes(x = mean.DSI, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Drought Severity Index", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-2.2, 2.9),
                     breaks = seq(-2, 3, by = 1),
                     labels = c(-0.85, -0.66, -0.45, -0.25, -0.05, 0.15))+
  theme_classic(base_size = 22)
grass.DSI.plot

ggsave("./Plots/grass.traits.NW.DSI.pdf", height = 7, width = 7)

MAP.effects = grass.effects$mean.MAP

attr(grass.imputed.NW.2$mean.MAP, "scaled:scale")
attr(grass.imputed.NW.2$mean.MAP, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*441.5584) + 637.2953

grass.MAP.plot = ggplot() +
  geom_point(data = grass.imputed.NW.2, aes(x = mean.MAP, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = MAP.effects, aes(x = mean.MAP, y = estimate__), color = "#6E687E", size = 1.5) +  
  geom_ribbon(data = MAP.effects, aes(x = mean.MAP, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Mean grass Precipitation", y = "Percent Cover Change") +
  ylim(-20.7,21.3)+
  scale_x_continuous(limits = c(-1.2, 4.0),
                     breaks = seq(-1, 4, by = 1),
                     labels = c(195.74, 637.29, 1078.85, 1520.41, 1961.97,2403.52))+
  theme_classic(base_size = 22)
grass.MAP.plot

ggsave("./Plots/grass.traits.NW.MAP.pdf", height = 7, width = 7)
