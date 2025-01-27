## Script to evaluate linear mixed effects models of cover change and traits

# load libraries 
library(tidyverse)
library(brms)

#### Load imputed data with woody species ####

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)

# imputed traits with woody species
imputed = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.outliersRM.csv", row.names = 1) %>%
  left_join(., enviro, by = "site_code") %>%
  mutate_at(vars(26:36),scale)
imputed.annual = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.annuals.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  mutate_at(vars(26:36),scale)
imputed.perennial = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.perennials.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  mutate_at(vars(26:36),scale)
imputed.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.forbs.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  mutate_at(vars(26:36),scale)
imputed.grass = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.grass.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  mutate_at(vars(26:36),scale)
imputed.annual.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.annual.forbs.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  mutate_at(vars(26:36),scale)
imputed.perennial.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.perennial.forb.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  mutate_at(vars(26:36),scale)
imputed.perennial.grass = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.perennial.grass.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
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

#### imputed traits lifespan model ####

piors <- c(prior(normal(0, 10), class = b))

annual.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                            SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                          family = gaussian(),
                          prior = priors,
                          data = imputed.annual)

summary(annual.traits.model)
# nothing significant

perennial.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                               SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                             family = gaussian(),
                             prior = priors,
                             data = imputed.perennial)

summary(perennial.traits.model)
# nothing significant

#### impute traits functional group ####

priors <- c(prior(normal(0, 10), class = b))

forb.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                          SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                        family = gaussian(),
                        prior = priors,
                        data = imputed.forb)

summary(forb.traits.model)
# nothing significant

grass.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                           SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                         family = gaussian(),
                         prior = priors,
                         data = imputed.grass)

summary(grass.traits.model)
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

perennial.forb.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                    SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                  family = gaussian(),
                                  prior = priors,
                                  data = imputed.perennial.forb)

summary(perennial.forb.traits.model)
# nothing significant

perennial.grass.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                     SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                   family = gaussian(),
                                   prior = priors,
                                   data = imputed.perennial.grass)

summary(perennial.grass.traits.model)
# nothing significant




#### Load imputed data without woody species ####

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)

# imputed traits with woody species
imputed.NW = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.outliersRM.csv", row.names = 1) %>%
  left_join(., enviro, by = "site_code") %>%
  mutate_at(vars(26:36),scale)
imputed.NW.annual = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.annuals.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  mutate_at(vars(26:36),scale)
imputed.NW.perennial = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennials.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  mutate_at(vars(26:36),scale)
imputed.NW.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.forbs.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  mutate_at(vars(26:36),scale)
imputed.NW.grass = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.grass.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  mutate_at(vars(26:36),scale)
imputed.NW.annual.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.annual.forbs.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  mutate_at(vars(26:36),scale)
imputed.NW.perennial.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennial.forb.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  mutate_at(vars(26:36),scale)
imputed.NW.perennial.grass = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennial.grass.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  mutate_at(vars(26:36),scale)

#### imputed traits model NW ####

priors <- c(prior(normal(0, 10), class = b))

imputed.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                             SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                           family = gaussian(),
                           prior = priors,
                           data = imputed.NW)

summary(imputed.traits.NW.model)
# nothing significant

#### imputed traits lifespan model NW ####

piors <- c(prior(normal(0, 10), class = b))

annual.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                            SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                          family = gaussian(),
                          prior = priors,
                          data = imputed.NW.annual)

summary(annual.traits.NW.model)
# nothing significant

perennial.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                               SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                             family = gaussian(),
                             prior = priors,
                             data = imputed.NW.perennial)

summary(perennial.traits.NW.model)
# nothing significant

#### impute traits functional group NW ####

priors <- c(prior(normal(0, 10), class = b))

forb.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                          SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                        family = gaussian(),
                        prior = priors,
                        data = imputed.NW.forb)

summary(forb.traits.NW.model)
# height significant

grass.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                           SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                         family = gaussian(),
                         prior = priors,
                         data = imputed.NW.grass)

summary(grass.traits.NW.model)
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

perennial.forb.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                    SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                  family = gaussian(),
                                  prior = priors,
                                  data = imputed.NW.perennial.forb)

summary(perennial.forb.traits.NW.model)
# nothing significant

perennial.grass.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                     SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                   family = gaussian(),
                                   prior = priors,
                                   data = imputed.NW.perennial.grass)

summary(perennial.grass.traits.NW.model)
# nothing significant

#### Load CC data with woody species ####

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)

# imputed traits with woody species
all.CC = read.csv("./Formatted.Data/Revisions/Final.Data/final.data.CC.outliersRM.csv", row.names = 1) %>%
  left_join(., enviro, by = "site_code") %>%
  mutate_at(vars(4:12,16,17),scale)
perennial.CC = read.csv("./Formatted.Data/Revisions/Final.Data/final.data.CC.perennials.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  mutate_at(vars(4:12,16,17),scale)
forb.CC = read.csv("./Formatted.Data/Revisions/Final.Data/final.data.CC.forbs.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  mutate_at(vars(4:12,16,17),scale)
grass.CC = read.csv("./Formatted.Data/Revisions/Final.Data/final.data.CC.grass.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  mutate_at(vars(4:12,16,17),scale)

#### CC traits ####
priors <- c(prior(normal(0, 10), class = b))

CC.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                 SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
               family = gaussian(),
               prior = priors,
               data = all.CC)

summary(CC.model)
# RTD significant

#### CC traits lifespan ####

priors <- c(prior(normal(0, 10), class = b))

perennial.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                        SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                      family = gaussian(),
                      prior = priors,
                      data = perennial.CC)

summary(perennial.model)
# SLA

#### CC traits functional group ####

priors <- c(prior(normal(0, 10), class = b))

forb.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                   SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                 family = gaussian(),
                 prior = priors,
                 data = forb.CC)

summary(forb.model)
# nothing significant

grass.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                   SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                 family = gaussian(),
                 prior = priors,
                 data = grass.CC)

summary(grass.model)
# MAP

#### Load CC data without woody species ####

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)

# imputed traits with woody species
all.CC.NW = read.csv("./Formatted.Data/Revisions/Final.Data/final.data.CC.NW.outliersRM.csv", row.names = 1) %>%
  left_join(., enviro, by = "site_code") %>%
  mutate_at(vars(4:12,16,17),scale)
perennial.CC.NW = read.csv("./Formatted.Data/Revisions/Final.Data/final.data.CC.NW.perennials.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  mutate_at(vars(4:12,16,17),scale)
forb.CC.NW = read.csv("./Formatted.Data/Revisions/Final.Data/final.data.CC.NW.forbs.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  mutate_at(vars(4:12,16,17),scale)
grass.CC.NW = read.csv("./Formatted.Data/Revisions/Final.Data/final.data.CC.NW.grass.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  mutate_at(vars(4:12,16,17),scale)

#### CC traits NW ####

priors <- c(prior(normal(0, 10), class = b))

CC.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                 SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
               family = gaussian(),
               prior = priors,
               data = all.CC.NW)

summary(CC.model)
# RTD

#### CC traits lifespan NW ####

priors <- c(prior(normal(0, 10), class = b))

perennial.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                        SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                      family = gaussian(),
                      prior = priors,
                      data = perennial.CC.NW)

summary(perennial.model)
# SLA

#### CC traits NW functional group ####

priors <- c(prior(normal(0, 10), class = b))

forb.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                   SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                 family = gaussian(),
                 prior = priors,
                 data = forb.CC.NW)

summary(forb.model)

grass.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                   SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                 family = gaussian(),
                 prior = priors,
                 data = grass.CC.NW)

summary(grass.model)
# RTD, MAP

