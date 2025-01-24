
# enviro data

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)

#### imputed traits ####

imputed = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)
imputed.traits = left_join(imputed, enviro, by="site_code")

imputed.traits.2 = scale(imputed.traits[,c(26:36)])
imputed.traits.3 = cbind(imputed.traits[,c(1:3)],imputed.traits.2)

priors <- c(prior(normal(0, 10), class = b))

imputed.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                             SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                           family = gaussian(),
                           prior = priors,
                           data = imputed.traits.3)

summary(imputed.traits.model)
# nothing significant

#### imputed traits outlier #####

imputed = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)
imputed.traits = left_join(imputed, enviro, by="site_code")

# 94 total removed because some removed with multiple traits
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

imputed.traits.3 = scale(imputed.traits.2[,c(26:36)])
imputed.traits.4 = cbind(imputed.traits.2[,c(1:3)],imputed.traits.3)

priors <- c(prior(normal(0, 10), class = b))

imputed.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                             SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                           family = gaussian(),
                           prior = priors,
                           data = imputed.traits.4)

summary(imputed.traits.model)
# nothing significant

#### imputed traits lifespan ####

imputed = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)
imputed.traits = left_join(imputed, enviro, by="site_code")

annual.imputed = imputed.traits %>%
  filter(local_lifespan == "ANNUAL")
perennial.imputed = imputed.traits %>%
  filter(local_lifespan == "PERENNIAL")

annual.imputed.2 = scale(annual.imputed[,c(26:36)])
annual.imputed.3 = cbind(annual.imputed[,c(1:3)],annual.imputed.2)

perennial.imputed.2 = scale(perennial.imputed[,c(26:36)])
perennial.imputed.3 = cbind(perennial.imputed[,c(1:3)],perennial.imputed.2)

priors <- c(prior(normal(0, 10), class = b))

annual.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                            SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                          family = gaussian(),
                          prior = priors,
                          data = annual.imputed.3)

summary(annual.traits.model)
# nothing significant

perennial.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                               SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                             family = gaussian(),
                             prior = priors,
                             data = perennial.imputed.3)

summary(perennial.traits.model)
# RTD and SLA significant, positive

#### impute traits functional group ####

imputed = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)
imputed.traits = left_join(imputed, enviro, by="site_code")

forb.imputed = imputed.traits %>%
  filter(functional_group == "FORB")
grass.imputed = imputed.traits %>%
  filter(functional_group == "GRASS")

forb.imputed.2 = scale(forb.imputed[,c(26:36)])
forb.imputed.3 = cbind(forb.imputed[,c(1:3)],forb.imputed.2)

grass.imputed.2 = scale(grass.imputed[,c(26:36)])
grass.imputed.3 = cbind(grass.imputed[,c(1:3)],grass.imputed.2)

priors <- c(prior(normal(0, 10), class = b))

forb.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                          SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                        family = gaussian(),
                        prior = priors,
                        data = forb.imputed.3)

summary(forb.traits.model)

grass.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                           SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                         family = gaussian(),
                         prior = priors,
                         data = grass.imputed.3)

summary(grass.traits.model)
# MAP positive significant

#### imputed traits lifespan x functional group ####

imputed = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)
imputed.traits = left_join(imputed, enviro, by="site_code")

annual.forb.imputed = imputed.traits %>%
  filter(functional_group == "FORB" & local_lifespan == "ANNUAL")
perennial.forb.imputed = imputed.traits %>%
  filter(functional_group == "FORB" & local_lifespan == "PERENNIAL")
perennial.grass.imputed = imputed.traits %>%
  filter(functional_group == "GRASS" & local_lifespan == "PERENNIAL")

annual.forb.imputed.2 = scale(annual.forb.imputed[,c(26:36)])
annual.forb.imputed.3 = cbind(annual.forb.imputed[,c(1:3)],annual.forb.imputed.2)

perennial.forb.imputed.2 = scale(perennial.forb.imputed[,c(26:36)])
perennial.forb.imputed.3 = cbind(perennial.forb.imputed[,c(1:3)],perennial.forb.imputed.2)

perennial.grass.imputed.2 = scale(perennial.grass.imputed[,c(26:36)])
perennial.grass.imputed.3 = cbind(perennial.grass.imputed[,c(1:3)],perennial.grass.imputed.2)

priors <- c(prior(normal(0, 10), class = b))

annual.forb.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                 SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                               family = gaussian(),
                               prior = priors,
                               data = annual.forb.imputed.3)

summary(annual.forb.traits.model)
# MAP significant

perennial.forb.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                    SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                  family = gaussian(),
                                  prior = priors,
                                  data = perennial.forb.imputed.3)

summary(perennial.forb.traits.model)

perennial.grass.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                     SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                   family = gaussian(),
                                   prior = priors,
                                   data = perennial.grass.imputed.3)

summary(perennial.grass.traits.model)

#### imputed traits lifespan outlier ####

imputed = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)
imputed.traits = left_join(imputed, enviro, by="site_code")

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

annual.imputed = imputed.traits.2 %>%
  filter(local_lifespan == "ANNUAL")
perennial.imputed = imputed.traits.2 %>%
  filter(local_lifespan == "PERENNIAL")

annual.imputed.2 = scale(annual.imputed[,c(26:36)])
annual.imputed.3 = cbind(annual.imputed[,c(1:3)],annual.imputed.2)

perennial.imputed.2 = scale(perennial.imputed[,c(26:36)])
perennial.imputed.3 = cbind(perennial.imputed[,c(1:3)],perennial.imputed.2)

priors <- c(prior(normal(0, 10), class = b))


annual.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                            SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                          family = gaussian(),
                          prior = priors,
                          data = annual.imputed.3)

summary(annual.traits.model)
# nothing significant

perennial.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                               SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                             family = gaussian(),
                             prior = priors,
                             data = perennial.imputed.3)

summary(perennial.traits.model)

#### impute traits functional group outlier ####

imputed = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)
imputed.traits = left_join(imputed, enviro, by="site_code")

imputed.traits.2 = imputed.traits %>%
  filter(round(leafN.final, 5) < 46.00560) %>% # 4
  filter(round(height.final, 5) < 7.787850) %>% # 11
  filter(round(rootN.final, 5) < 25.55616) %>% # 18
  filter(round(SLA.final, 5) < 49.05439) %>% # 8
  filter(round(root.depth.final, 5) < 2.915000) %>% # 15
  filter(round(RTD.final, 5) < 0.7460162) %>% # 11
  filter(round(SRL.final, 5) < 437.3523) %>% # 21
  filter(round(rootDiam.final, 5) < 1.200000) %>% # 6
  filter(round(RMF.final, 5) < 0.7670084)

forb.imputed = imputed.traits.2 %>%
  filter(functional_group == "FORB")
grass.imputed = imputed.traits.2 %>%
  filter(functional_group == "GRASS")

forb.imputed.2 = scale(forb.imputed[,c(26:36)])
forb.imputed.3 = cbind(forb.imputed[,c(1:3)],forb.imputed.2)

grass.imputed.2 = scale(grass.imputed[,c(26:36)])
grass.imputed.3 = cbind(grass.imputed[,c(1:3)],grass.imputed.2)


priors <- c(prior(normal(0, 10), class = b))

forb.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                          SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                        family = gaussian(),
                        prior = priors,
                        data = forb.imputed.3)

summary(forb.traits.model)

grass.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                           SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                         family = gaussian(),
                         prior = priors,
                         data = grass.imputed.3)

summary(grass.traits.model)

#### imputed traits lifespan x functional group outlier ####

imputed = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)
imputed.traits = left_join(imputed, enviro, by="site_code")

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

annual.forb.imputed = imputed.traits.2 %>%
  filter(functional_group == "FORB" & local_lifespan == "ANNUAL")
perennial.forb.imputed = imputed.traits.2 %>%
  filter(functional_group == "FORB" & local_lifespan == "PERENNIAL")
perennial.grass.imputed = imputed.traits.2 %>%
  filter(functional_group == "GRASS" & local_lifespan == "PERENNIAL")


annual.forb.imputed.2 = scale(annual.forb.imputed[,c(26:36)])
annual.forb.imputed.3 = cbind(annual.forb.imputed[,c(1:3)],annual.forb.imputed.2)

perennial.forb.imputed.2 = scale(perennial.forb.imputed[,c(26:36)])
perennial.forb.imputed.3 = cbind(perennial.forb.imputed[,c(1:3)],perennial.forb.imputed.2)

perennial.grass.imputed.2 = scale(perennial.grass.imputed[,c(26:36)])
perennial.grass.imputed.3 = cbind(perennial.grass.imputed[,c(1:3)],perennial.grass.imputed.2)

priors <- c(prior(normal(0, 10), class = b))

annual.forb.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                 SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                               family = gaussian(),
                               prior = priors,
                               data = annual.forb.imputed.3)

summary(annual.forb.traits.model)

perennial.forb.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                    SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                  family = gaussian(),
                                  prior = priors,
                                  data = perennial.forb.imputed.3)

summary(perennial.forb.traits.model)

perennial.grass.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                     SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                   family = gaussian(),
                                   prior = priors,
                                   data = perennial.grass.imputed.3)

summary(perennial.grass.traits.model)

#### CC traits ####
final.data = read.csv("./Formatted.Data/Revisions/final.data.CC.csv",row.names = 1)
final.data.CC = left_join(final.data, enviro, by="site_code")

final.data.CC.2 = scale(final.data.CC[,c(4:12,16,17)])
final.data.CC.3 = cbind(final.data.CC[,c(1:3)],final.data.CC.2)

# 236 data points estimating 11 variables

priors <- c(prior(normal(0, 10), class = b))

CC.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                 SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
               family = gaussian(),
               prior = priors,
               data = final.data.CC.3)

summary(CC.model)
# RTD significant

#### CC traits outlier ####

final.data = read.csv("./Formatted.Data/Revisions/final.data.CC.csv",row.names = 1)
final.data.CC = left_join(final.data, enviro, by="site_code")

# 25 total removed because some removed with multiple traits
final.data.CC.2 = final.data.CC %>%
  filter(round(leafN.mg.g, 5) != 54.71597) %>% # 2
  filter(round(height.m, 5) < 12.84847) %>% # 3
  filter(round(rootN.mg.g, 5) < 31.17241) %>% # 5
  filter(round(root.depth_m, 5) < 2.426850) %>% # 7
  filter(round(RTD.groot.cahill.merge, 5) < 1.19455) %>% # 4
  filter(round(SRL.groot.cahill.merge, 5) < 471.2364) %>% # 5
  filter(round(rootDiam.mm, 5) != 1.3965) # 1

final.data.CC.3 = scale(final.data.CC.2[,c(4:12,16,17)])
final.data.CC.4 = cbind(final.data.CC.2[,c(1:3)],final.data.CC.3)

priors <- c(prior(normal(0, 10), class = b))


CC.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                 SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
               family = gaussian(),
               prior = priors,
               data = final.data.CC.4)

summary(CC.model)
# RTD significant

#### CC traits lifespan ####

final.data = read.csv("./Formatted.Data/Revisions/final.data.CC.csv",row.names = 1)
final.data.CC = left_join(final.data, enviro, by="site_code")

table(final.data.CC$local_lifespan)

perennial = final.data.CC %>%
  filter(local_lifespan == "PERENNIAL")

perennial.2 = scale(perennial[,c(4:12,16,17)])
perennial.3 = cbind(perennial[,c(1:3)],perennial.2)

priors <- c(prior(normal(0, 10), class = b))

perennial.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                        SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                      family = gaussian(),
                      prior = priors,
                      data = perennial.3)

summary(perennial.model)

#### CC traits functional group ####

final.data = read.csv("./Formatted.Data/Revisions/final.data.CC.csv",row.names = 1)
final.data.CC = left_join(final.data, enviro, by="site_code")

table(final.data.CC$functional_group)

forb = final.data.CC %>%
  filter(functional_group == "FORB")

grass = final.data.CC %>%
  filter(functional_group == "GRASS")

forb.2 = scale(forb[,c(4:12,16,17)])
forb.3 = cbind(forb[,c(1:3)],forb.2)

grass.2 = scale(grass[,c(4:12,16,17)])
grass.3 = cbind(grass[,c(1:3)],grass.2)

priors <- c(prior(normal(0, 10), class = b))

forb.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                   SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                 family = gaussian(),
                 prior = priors,
                 data = forb.3)

summary(forb.model)

grass.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                   SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                 family = gaussian(),
                 prior = priors,
                 data = grass.3)

summary(grass.model)

#### CC traits lifespan outlier ####

final.data = read.csv("./Formatted.Data/Revisions/final.data.CC.csv",row.names = 1)
final.data.CC = left_join(final.data, enviro, by="site_code")

# 25 total removed because some removed with multiple traits
final.data.CC.2 = final.data.CC %>%
  filter(round(leafN.mg.g, 5) != 54.71597) %>% # 2
  filter(round(height.m, 5) < 12.84847) %>% # 3
  filter(round(rootN.mg.g, 5) < 31.17241) %>% # 5
  filter(round(root.depth_m, 5) < 2.426850) %>% # 7
  filter(round(RTD.groot.cahill.merge, 5) < 1.19455) %>% # 4
  filter(round(SRL.groot.cahill.merge, 5) < 471.2364) %>% # 5
  filter(round(rootDiam.mm, 5) != 1.3965) # 1

perennial = final.data.CC.2 %>%
  filter(local_lifespan == "PERENNIAL")

perennial.2 = scale(perennial[,c(4:12,16,17)])
perennial.3 = cbind(perennial[,c(1:3)],perennial.2)

priors <- c(prior(normal(0, 10), class = b))

perennial.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                        SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                      family = gaussian(),
                      prior = priors,
                      data = perennial.3)

summary(perennial.model)
# SLA significant

#### CC traits functional group outlier ####

final.data = read.csv("./Formatted.Data/Revisions/final.data.CC.csv",row.names = 1)
final.data.CC = left_join(final.data, enviro, by="site_code")

final.data.CC.2 = final.data.CC %>%
  filter(round(leafN.mg.g, 5) != 54.71597) %>% # 2
  filter(round(height.m, 5) < 12.84847) %>% # 3
  filter(round(rootN.mg.g, 5) < 31.17241) %>% # 5
  filter(round(root.depth_m, 5) < 2.426850) %>% # 7
  filter(round(RTD.groot.cahill.merge, 5) < 1.19455) %>% # 4
  filter(round(SRL.groot.cahill.merge, 5) < 471.2364) %>% # 5
  filter(round(rootDiam.mm, 5) != 1.3965) # 1

forb = final.data.CC.2 %>%
  filter(functional_group == "FORB")

grass = final.data.CC.2 %>%
  filter(functional_group == "GRASS")

forb.2 = scale(forb[,c(4:12,16,17)])
forb.3 = cbind(forb[,c(1:3)],forb.2)

grass.2 = scale(grass[,c(4:12,16,17)])
grass.3 = cbind(grass[,c(1:3)],grass.2)

priors <- c(prior(normal(0, 10), class = b))

forb.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                   SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                 family = gaussian(),
                 prior = priors,
                 data = forb.3)

summary(forb.model)

grass.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                    SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                  family = gaussian(),
                  prior = priors,
                  data = grass.3)

summary(grass.model)


#### imputed traits NW ####

imputed = read.csv("./Formatted.Data/Revisions/imputed.traits.NW.final.csv",row.names = 1)
imputed.traits = left_join(imputed, enviro, by="site_code")

imputed.traits.2 = scale(imputed.traits[,c(26:34)])
imputed.traits.3 = cbind(imputed.traits[,c(1:3)],imputed.traits.2)

priors <- c(prior(normal(0, 10), class = b))

imputed.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                             SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                           family = gaussian(),
                           prior = priors,
                           data = imputed.traits.3)

summary(imputed.traits.model)
# nothing significant

#### imputed traits NW outlier ####

imputed = read.csv("./Formatted.Data/Revisions/imputed.traits.NW.final.csv",row.names = 1)
imputed.traits = left_join(imputed, enviro, by="site_code")

imputed.traits.2 = imputed.traits %>%
  filter(round(leafN.final, 5) < 45.18532) %>% # 4
  filter(round(height.final, 5) < 1.371600) %>% # 11
  filter(round(rootN.final, 5) < 27.18327) %>% # 18
  filter(round(SLA.final, 5) < 49.05439) %>% # 8
  filter(round(root.depth.final, 5) < 6.451600) %>% # 15
  filter(round(RTD.final, 5) < 0.6525000) %>% # 11
  filter(round(SRL.final, 5) < 437.3523) %>% # 21
  filter(round(rootDiam.final, 5) < 1.126012) %>% # 6
  filter(round(RMF.final, 5) > 0.04850791 & round(RMF.final, 5) < 0.76700839)

imputed.traits.3 = scale(imputed.traits.2[,c(26:34)])
imputed.traits.4 = cbind(imputed.traits.2[,c(1:3)],imputed.traits.3)

priors <- c(prior(normal(0, 10), class = b))

imputed.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                             SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                           family = gaussian(),
                           prior = priors,
                           data = imputed.traits.4)

summary(imputed.traits.model)
# nothing significant

#### imputed traits lifespan NW ####

imputed = read.csv("./Formatted.Data/Revisions/imputed.traits.NW.final.csv",row.names = 1)
imputed.traits = left_join(imputed, enviro, by="site_code")

annual.imputed = imputed.traits %>%
  filter(local_lifespan == "ANNUAL")
perennial.imputed = imputed.traits %>%
  filter(local_lifespan == "PERENNIAL")

annual.imputed.2 = scale(annual.imputed[,c(26:34)])
annual.imputed.3 = cbind(annual.imputed[,c(1:3)],annual.imputed.2)

perennial.imputed.2 = scale(perennial.imputed[,c(26:34)])
perennial.imputed.3 = cbind(perennial.imputed[,c(1:3)],perennial.imputed.2)

priors <- c(prior(normal(0, 10), class = b))


annual.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                            SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                          family = gaussian(),
                          prior = priors,
                          data = annual.imputed.3)

summary(annual.traits.model)
# nothing significant

perennial.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                               SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                             family = gaussian(),
                             prior = priors,
                             data = perennial.imputed.3)

summary(perennial.traits.model)

#### impute traits functional group NW ####

imputed = read.csv("./Formatted.Data/Revisions/imputed.traits.NW.final.csv",row.names = 1)
imputed.traits = left_join(imputed, enviro, by="site_code")

forb.imputed = imputed.traits %>%
  filter(functional_group == "FORB")
grass.imputed = imputed.traits %>%
  filter(functional_group == "GRASS")

forb.imputed.2 = scale(forb.imputed[,c(26:34)])
forb.imputed.3 = cbind(forb.imputed[,c(1:3)],forb.imputed.2)

grass.imputed.2 = scale(grass.imputed[,c(26:34)])
grass.imputed.3 = cbind(grass.imputed[,c(1:3)],grass.imputed.2)

priors <- c(prior(normal(0, 10), class = b))


forb.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                          SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                        family = gaussian(),
                        prior = priors,
                        data = forb.imputed.3)

summary(forb.traits.model)

grass.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                           SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                         family = gaussian(),
                         prior = priors,
                         data = grass.imputed.3)

summary(grass.traits.model)
# nothing significant

#### imputed traits lifespan x functional group NW ####

imputed = read.csv("./Formatted.Data/Revisions/imputed.traits.NW.final.csv",row.names = 1)
imputed.traits = left_join(imputed, enviro, by="site_code")

annual.forb.imputed = imputed.traits %>%
  filter(functional_group == "FORB" & local_lifespan == "ANNUAL")
perennial.forb.imputed = imputed.traits %>%
  filter(functional_group == "FORB" & local_lifespan == "PERENNIAL")
perennial.grass.imputed = imputed.traits %>%
  filter(functional_group == "GRASS" & local_lifespan == "PERENNIAL")

annual.forb.imputed.2 = scale(annual.forb.imputed[,c(26:34)])
annual.forb.imputed.3 = cbind(annual.forb.imputed[,c(1:3)],annual.forb.imputed.2)

perennial.forb.imputed.2 = scale(perennial.forb.imputed[,c(26:34)])
perennial.forb.imputed.3 = cbind(perennial.forb.imputed[,c(1:3)],perennial.forb.imputed.2)

perennial.grass.imputed.2 = scale(perennial.grass.imputed[,c(26:34)])
perennial.grass.imputed.3 = cbind(perennial.grass.imputed[,c(1:3)],perennial.grass.imputed.2)

priors <- c(prior(normal(0, 10), class = b))

annual.forb.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                 SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                               family = gaussian(),
                               prior = priors,
                               data = annual.forb.imputed.3)

summary(annual.forb.traits.model)
# nothing significant

perennial.forb.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                    SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                                  family = gaussian(),
                                  prior = priors,
                                  data = perennial.forb.imputed.3)

summary(perennial.forb.traits.model)

perennial.grass.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                     SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                                   family = gaussian(),
                                   prior = priors,
                                   data = perennial.grass.imputed.3)

summary(perennial.grass.traits.model)


#### imputed traits lifespan NW outlier ####

imputed = read.csv("./Formatted.Data/Revisions/imputed.traits.NW.final.csv",row.names = 1)
imputed.traits = left_join(imputed, enviro, by="site_code")

imputed.traits.2 = imputed.traits %>%
  filter(round(leafN.final, 5) < 45.18532) %>% # 4
  filter(round(height.final, 5) < 1.371600) %>% # 11
  filter(round(rootN.final, 5) < 27.18327) %>% # 18
  filter(round(SLA.final, 5) < 49.05439) %>% # 8
  filter(round(root.depth.final, 5) < 6.451600) %>% # 15
  filter(round(RTD.final, 5) < 0.6525000) %>% # 11
  filter(round(SRL.final, 5) < 437.3523) %>% # 21
  filter(round(rootDiam.final, 5) < 1.126012) %>% # 6
  filter(round(RMF.final, 5) > 0.04850791 & round(RMF.final, 5) < 0.76700839) # 10

annual.imputed = imputed.traits.2 %>%
  filter(local_lifespan == "ANNUAL")
perennial.imputed = imputed.traits.2 %>%
  filter(local_lifespan == "PERENNIAL")

annual.imputed.2 = scale(annual.imputed[,c(26:34)])
annual.imputed.3 = cbind(annual.imputed[,c(1:3)],annual.imputed.2)

perennial.imputed.2 = scale(perennial.imputed[,c(26:34)])
perennial.imputed.3 = cbind(perennial.imputed[,c(1:3)],perennial.imputed.2)

priors <- c(prior(normal(0, 10), class = b))


annual.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                            SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                          family = gaussian(),
                          prior = priors,
                          data = annual.imputed.3)

summary(annual.traits.model)
# rootN significant

perennial.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                               SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                             family = gaussian(),
                             prior = priors,
                             data = perennial.imputed.3)

summary(perennial.traits.model)

#### impute traits functional group NW outlier ####

imputed = read.csv("./Formatted.Data/Revisions/imputed.traits.NW.final.csv",row.names = 1)
imputed.traits = left_join(imputed, enviro, by="site_code")

imputed.traits.2 = imputed.traits %>%
  filter(round(leafN.final, 5) < 45.18532) %>% # 4
  filter(round(height.final, 5) < 1.371600) %>% # 11
  filter(round(rootN.final, 5) < 27.18327) %>% # 18
  filter(round(SLA.final, 5) < 49.05439) %>% # 8
  filter(round(root.depth.final, 5) < 6.451600) %>% # 15
  filter(round(RTD.final, 5) < 0.6525000) %>% # 11
  filter(round(SRL.final, 5) < 437.3523) %>% # 21
  filter(round(rootDiam.final, 5) < 1.126012) %>% # 6
  filter(round(RMF.final, 5) > 0.04850791 & round(RMF.final, 5) < 0.76700839) # 10

forb.imputed = imputed.traits.2 %>%
  filter(functional_group == "FORB")
grass.imputed = imputed.traits.2 %>%
  filter(functional_group == "GRASS")

forb.imputed.2 = scale(forb.imputed[,c(26:34)])
forb.imputed.3 = cbind(forb.imputed[,c(1:3)],forb.imputed.2)

grass.imputed.2 = scale(grass.imputed[,c(26:34)])
grass.imputed.3 = cbind(grass.imputed[,c(1:3)],grass.imputed.2)

priors <- c(prior(normal(0, 10), class = b))


forb.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                          SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                        family = gaussian(),
                        prior = priors,
                        data = forb.imputed.3)

summary(forb.traits.model)

grass.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                           SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                         family = gaussian(),
                         prior = priors,
                         data = grass.imputed.3)

summary(grass.traits.model)
# nothing significant

#### imputed traits lifespan x functional group NW outlier ####

imputed = read.csv("./Formatted.Data/Revisions/imputed.traits.NW.final.csv",row.names = 1)
imputed.traits = left_join(imputed, enviro, by="site_code")

imputed.traits.2 = imputed.traits %>%
  filter(round(leafN.final, 5) < 45.18532) %>% # 4
  filter(round(height.final, 5) < 1.371600) %>% # 11
  filter(round(rootN.final, 5) < 27.18327) %>% # 18
  filter(round(SLA.final, 5) < 49.05439) %>% # 8
  filter(round(root.depth.final, 5) < 6.451600) %>% # 15
  filter(round(RTD.final, 5) < 0.6525000) %>% # 11
  filter(round(SRL.final, 5) < 437.3523) %>% # 21
  filter(round(rootDiam.final, 5) < 1.126012) %>% # 6
  filter(round(RMF.final, 5) > 0.04850791 & round(RMF.final, 5) < 0.76700839) # 10

annual.forb.imputed = imputed.traits.2 %>%
  filter(functional_group == "FORB" & local_lifespan == "ANNUAL")
perennial.forb.imputed = imputed.traits.2 %>%
  filter(functional_group == "FORB" & local_lifespan == "PERENNIAL")
perennial.grass.imputed = imputed.traits.2 %>%
  filter(functional_group == "GRASS" & local_lifespan == "PERENNIAL")

annual.forb.imputed.2 = scale(annual.forb.imputed[,c(26:34)])
annual.forb.imputed.3 = cbind(annual.forb.imputed[,c(1:3)],annual.forb.imputed.2)

perennial.forb.imputed.2 = scale(perennial.forb.imputed[,c(26:34)])
perennial.forb.imputed.3 = cbind(perennial.forb.imputed[,c(1:3)],perennial.forb.imputed.2)

perennial.grass.imputed.2 = scale(perennial.grass.imputed[,c(26:34)])
perennial.grass.imputed.3 = cbind(perennial.grass.imputed[,c(1:3)],perennial.grass.imputed.2)

priors <- c(prior(normal(0, 10), class = b))

annual.forb.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                 SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                               family = gaussian(),
                               prior = priors,
                               data = annual.forb.imputed.3)

summary(annual.forb.traits.model)
# nothing significant

perennial.forb.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                    SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                                  family = gaussian(),
                                  prior = priors,
                                  data = perennial.forb.imputed.3)

summary(perennial.forb.traits.model)

perennial.grass.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                     SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                                   family = gaussian(),
                                   prior = priors,
                                   data = perennial.grass.imputed.3)

summary(perennial.grass.traits.model)


#### CC traits NW ####
final.data = read.csv("./Formatted.Data/Revisions/final.data.CC.NW.csv",row.names = 1)
final.data.CC = left_join(final.data, enviro, by="site_code")

final.data.CC.2 = scale(final.data.CC[,c(4:12)])
final.data.CC.3 = cbind(final.data.CC[,c(1:3)],final.data.CC.2)

# 236 data points estimating 11 variables

priors <- c(prior(normal(0, 10), class = b))


CC.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                 SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + (1|site_code) + (1|Taxon), 
               family = gaussian(),
               prior = priors,
               data = final.data.CC.3)

summary(CC.model)

#### CC traits lifespan NW ####

final.data = read.csv("./Formatted.Data/Revisions/final.data.CC.NW.csv",row.names = 1)
final.data.CC = left_join(final.data, enviro, by="site_code")

table(final.data.CC$local_lifespan)

perennial = final.data.CC %>%
  filter(local_lifespan == "PERENNIAL")

perennial.2 = scale(perennial[,c(4:12)])
perennial.3 = cbind(perennial[,c(1:3)],perennial.2)

priors <- c(prior(normal(0, 10), class = b))

perennial.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                        SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + (1|site_code) + (1|Taxon), 
                      family = gaussian(),
                      prior = priors,
                      data = perennial.3)

summary(perennial.model)

#### CC traits NW functional group ####

final.data = read.csv("./Formatted.Data/Revisions/final.data.CC.NW.csv",row.names = 1)
final.data.CC = left_join(final.data, enviro, by="site_code")

table(final.data.CC$functional_group)

forb = final.data.CC %>%
  filter(functional_group == "FORB")

forb.2 = scale(forb[,c(4:12)])
forb.3 = cbind(forb[,c(1:3)],forb.2)

priors <- c(prior(normal(0, 10), class = b))

forb.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                   SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + (1|site_code) + (1|Taxon), 
                 family = gaussian(),
                 prior = priors,
                 data = forb.3)

summary(forb.model)


#### CC traits NW outlier ####
final.data = read.csv("./Formatted.Data/Revisions/final.data.CC.NW.csv",row.names = 1)
final.data.CC = left_join(final.data, enviro, by="site_code")

final.data.CC.2 = final.data.CC %>%
  filter(round(leafN.mg.g, 5) != 54.71597) %>% # 2
  filter(round(height.m, 5) < 1.800048) %>% # 3
  filter(round(rootN.mg.g, 5) < 31.17241) %>% # 5
  filter(round(root.depth_m, 5) < 2.426850) %>% # 7
  filter(round(SRL.groot.cahill.merge, 5) < 471.2364) # 5

final.data.CC.3 = scale(final.data.CC.2[,c(4:12)])
final.data.CC.4 = cbind(final.data.CC.2[,c(1:3)],final.data.CC.3)

# 236 data points estimating 11 variables

priors <- c(prior(normal(0, 10), class = b))


CC.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                 SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + (1|site_code) + (1|Taxon), 
               family = gaussian(),
               prior = priors,
               data = final.data.CC.4)

summary(CC.model)

#### CC traits lifespan NW outlier ####

final.data = read.csv("./Formatted.Data/Revisions/final.data.CC.NW.csv",row.names = 1)
final.data.CC = left_join(final.data, enviro, by="site_code")

table(final.data.CC$local_lifespan)

final.data.CC.2 = final.data.CC %>%
  filter(round(leafN.mg.g, 5) != 54.71597) %>% # 2
  filter(round(height.m, 5) < 1.800048) %>% # 3
  filter(round(rootN.mg.g, 5) < 31.17241) %>% # 5
  filter(round(root.depth_m, 5) < 2.426850) %>% # 7
  filter(round(SRL.groot.cahill.merge, 5) < 471.2364) # 5


perennial = final.data.CC.2 %>%
  filter(local_lifespan == "PERENNIAL")

perennial.2 = scale(perennial[,c(4:12)])
perennial.3 = cbind(perennial[,c(1:3)],perennial.2)

priors <- c(prior(normal(0, 10), class = b))

perennial.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                        SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + (1|site_code) + (1|Taxon), 
                      family = gaussian(),
                      prior = priors,
                      data = perennial.3)

summary(perennial.model)

#### CC traits NW functional group outlier ####

final.data = read.csv("./Formatted.Data/Revisions/final.data.CC.NW.csv",row.names = 1)
final.data.CC = left_join(final.data, enviro, by="site_code")

table(final.data.CC$functional_group)

final.data.CC.2 = final.data.CC %>%
  filter(round(leafN.mg.g, 5) != 54.71597) %>% # 2
  filter(round(height.m, 5) < 1.800048) %>% # 3
  filter(round(rootN.mg.g, 5) < 31.17241) %>% # 5
  filter(round(root.depth_m, 5) < 2.426850) %>% # 7
  filter(round(SRL.groot.cahill.merge, 5) < 471.2364) # 5

forb = final.data.CC.2 %>%
  filter(functional_group == "FORB")

forb.2 = scale(forb[,c(4:12)])
forb.3 = cbind(forb[,c(1:3)],forb.2)

priors <- c(prior(normal(0, 10), class = b))

forb.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                   SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + (1|site_code) + (1|Taxon), 
                 family = gaussian(),
                 prior = priors,
                 data = forb.3)

summary(forb.model)

