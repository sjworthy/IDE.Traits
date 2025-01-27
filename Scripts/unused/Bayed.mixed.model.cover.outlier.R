#### imputed traits ####

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)

hist(imputed.traits$cover.change)
boxplot(imputed.traits$cover.change)
mean = mean(imputed.traits$cover.change, na.rm = TRUE)
std = sd(imputed.traits$cover.change, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$cover.change[which(imputed.traits$cover.change <Tmin | imputed.traits$cover.change > Tmax)])
# -79.47619 - -23.66667 and 23.52000 - 55.00000

imputed.traits.2 = imputed.traits %>%
  filter(round(cover.change, 5) > -23.66667 & round(cover.change, 5) < 23.52000)

imputed.traits.3 = scale(imputed.traits.2[,c(26:34)])
imputed.traits.4 = cbind(imputed.traits.2[,c(1:3)],imputed.traits.3)

priors <- c(prior(normal(0, 10), class = b))

imputed.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                             SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                           family = gaussian(),
                           prior = priors,
                           data = imputed.traits.4)

summary(imputed.traits.model)
# nothing significant

#### imputed traits outlier #####

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)

imputed.traits.2 = imputed.traits %>%
  filter(round(cover.change, 5) > -23.66667 & round(cover.change, 5) < 23.52000)

# 94 total removed because some removed with multiple traits
imputed.traits.3 = imputed.traits.2 %>%
  filter(round(leafN.final, 5) < 46.00560) %>% # 4
  filter(round(height.final, 5) < 7.787850) %>% # 11
  filter(round(rootN.final, 5) < 25.55616) %>% # 18
  filter(round(SLA.final, 5) < 49.05439) %>% # 8
  filter(round(root.depth.final, 5) < 2.915000) %>% # 15
  filter(round(RTD.final, 5) < 0.7460162) %>% # 11
  filter(round(SRL.final, 5) < 437.3523) %>% # 21
  filter(round(rootDiam.final, 5) < 1.200000) %>% # 6
  filter(round(RMF.final, 5) < 0.7670084) # 10

imputed.traits.4 = scale(imputed.traits.3[,c(26:34)])
imputed.traits.5 = cbind(imputed.traits.3[,c(1:3)],imputed.traits.4)

priors <- c(prior(normal(0, 10), class = b))

imputed.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                             SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                           family = gaussian(),
                           prior = priors,
                           data = imputed.traits.5)

summary(imputed.traits.model)
# nothing significant

#### imputed traits lifespan ####

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)

imputed.traits.2 = imputed.traits %>%
  filter(round(cover.change, 5) > -23.66667 & round(cover.change, 5) < 23.52000)

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
# nothing significant

perennial.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                               SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                             family = gaussian(),
                             prior = priors,
                             data = perennial.imputed.3)

summary(perennial.traits.model)
# RTD significant, positive

#### impute traits functional group ####

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)

imputed.traits.2 = imputed.traits %>%
  filter(round(cover.change, 5) > -23.66667 & round(cover.change, 5) < 23.52000)

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
# leafN significant positive

grass.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                           SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                         family = gaussian(),
                         prior = priors,
                         data = grass.imputed.3)

summary(grass.traits.model)
# nothing significant

#### imputed traits lifespan x functional group ####

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)

imputed.traits.2 = imputed.traits %>%
  filter(round(cover.change, 5) > -23.66667 & round(cover.change, 5) < 23.52000)

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

#### imputed traits lifespan outlier ####

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)

imputed.traits.2 = imputed.traits %>%
  filter(round(cover.change, 5) > -23.66667 & round(cover.change, 5) < 23.52000)

# 19 total removed because some removed with multiple traits
imputed.traits.3 = imputed.traits.2 %>%
  filter(round(leafN.final, 5) < 46.00560) %>% # 4
  filter(round(height.final, 5) < 7.787850) %>% # 11
  filter(round(rootN.final, 5) < 25.55616) %>% # 18
  filter(round(SLA.final, 5) < 49.05439) %>% # 8
  filter(round(root.depth.final, 5) < 2.915000) %>% # 15
  filter(round(RTD.final, 5) < 0.7460162) %>% # 11
  filter(round(SRL.final, 5) < 437.3523) %>% # 21
  filter(round(rootDiam.final, 5) < 1.200000) %>% # 6
  filter(round(RMF.final, 5) < 0.7670084) # 10

annual.imputed = imputed.traits.3 %>%
  filter(local_lifespan == "ANNUAL")
perennial.imputed = imputed.traits.3 %>%
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

#### impute traits functional group outlier ####

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)

imputed.traits.2 = imputed.traits %>%
  filter(round(cover.change, 5) > -23.66667 & round(cover.change, 5) < 23.52000)

imputed.traits.3 = imputed.traits.2 %>%
  filter(round(leafN.final, 5) < 46.00560) %>% # 4
  filter(round(height.final, 5) < 7.787850) %>% # 11
  filter(round(rootN.final, 5) < 25.55616) %>% # 18
  filter(round(SLA.final, 5) < 49.05439) %>% # 8
  filter(round(root.depth.final, 5) < 2.915000) %>% # 15
  filter(round(RTD.final, 5) < 0.7460162) %>% # 11
  filter(round(SRL.final, 5) < 437.3523) %>% # 21
  filter(round(rootDiam.final, 5) < 1.200000) %>% # 6
  filter(round(RMF.final, 5) < 0.7670084) # 10

forb.imputed = imputed.traits.3 %>%
  filter(functional_group == "FORB")
grass.imputed = imputed.traits.3 %>%
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

#### imputed traits lifespan x functional group outlier ####

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)

imputed.traits.2 = imputed.traits %>%
  filter(round(cover.change, 5) > -23.66667 & round(cover.change, 5) < 23.52000)

imputed.traits.3 = imputed.traits.2 %>%
  filter(round(leafN.final, 5) < 46.00560) %>% # 4
  filter(round(height.final, 5) < 7.787850) %>% # 11
  filter(round(rootN.final, 5) < 25.55616) %>% # 18
  filter(round(SLA.final, 5) < 49.05439) %>% # 8
  filter(round(root.depth.final, 5) < 2.915000) %>% # 15
  filter(round(RTD.final, 5) < 0.7460162) %>% # 11
  filter(round(SRL.final, 5) < 437.3523) %>% # 21
  filter(round(rootDiam.final, 5) < 1.200000) %>% # 6
  filter(round(RMF.final, 5) < 0.7670084) # 10

annual.forb.imputed = imputed.traits.3 %>%
  filter(functional_group == "FORB" & local_lifespan == "ANNUAL")
perennial.forb.imputed = imputed.traits.3 %>%
  filter(functional_group == "FORB" & local_lifespan == "PERENNIAL")
perennial.grass.imputed = imputed.traits.3 %>%
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

#### CC traits ####
final.data.CC = read.csv("./Formatted.Data/Revisions/final.data.CC.csv",row.names = 1)

hist(final.data.CC$cover.change)
boxplot(final.data.CC$cover.change)
mean = mean(final.data.CC$cover.change, na.rm = TRUE)
std = sd(final.data.CC$cover.change, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$cover.change[which(final.data.CC$cover.change <Tmin | final.data.CC$cover.change > Tmax)])
# -79.47619 - -38.00000 and 28.91667

final.data.CC.2 = final.data.CC %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) != 28.91667)

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

#### CC traits outlier ####

final.data.CC = read.csv("./Formatted.Data/Revisions/final.data.CC.csv",row.names = 1)

final.data.CC.2 = final.data.CC %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) != 28.91667)

# 25 total removed because some removed with multiple traits
final.data.CC.3 = final.data.CC.2 %>%
  filter(round(leafN.mg.g, 5) != 54.71597) %>% # 2
  filter(round(height.m, 5) < 12.84847) %>% # 3
  filter(round(rootN.mg.g, 5) < 31.17241) %>% # 5
  filter(round(root.depth_m, 5) < 2.426850) %>% # 7
  filter(round(RTD.groot.cahill.merge, 5) < 1.19455) %>% # 4
  filter(round(SRL.groot.cahill.merge, 5) < 471.2364) %>% # 5
  filter(round(rootDiam.mm, 5) != 1.3965) # 1

final.data.CC.4 = scale(final.data.CC.3[,c(4:12)])
final.data.CC.5 = cbind(final.data.CC.3[,c(1:3)],final.data.CC.4)

priors <- c(prior(normal(0, 10), class = b))


CC.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                 SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + SRL.groot.cahill.merge*RTD.groot.cahill.merge + (1|site_code) + (1|Taxon), 
               family = gaussian(),
               prior = priors,
               data = final.data.CC.5)

summary(CC.model)

#### CC traits lifespan ####

final.data.CC = read.csv("./Formatted.Data/Revisions/final.data.CC.csv",row.names = 1)

final.data.CC.2 = final.data.CC %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) != 28.91667)

table(final.data.CC$local_lifespan)

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

#### CC traits functional group ####

final.data.CC = read.csv("./Formatted.Data/Revisions/final.data.CC.csv",row.names = 1)

final.data.CC.2 = final.data.CC %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) != 28.91667)

table(final.data.CC$functional_group)

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

#### CC traits lifespan outlier ####

final.data.CC = read.csv("./Formatted.Data/Revisions/final.data.CC.csv",row.names = 1)

final.data.CC.2 = final.data.CC %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) != 28.91667)

perennial = final.data.CC.2 %>%
  filter(local_lifespan == "PERENNIAL")

# 25 total removed because some removed with multiple traits
perennial.2 = perennial %>%
  filter(round(leafN.mg.g, 5) != 54.71597) %>% # 2
  filter(round(height.m, 5) < 12.84847) %>% # 3
  filter(round(rootN.mg.g, 5) < 31.17241) %>% # 5
  filter(round(root.depth_m, 5) < 2.426850) %>% # 7
  filter(round(RTD.groot.cahill.merge, 5) < 1.19455) %>% # 4
  filter(round(SRL.groot.cahill.merge, 5) < 471.2364) %>% # 5
  filter(round(rootDiam.mm, 5) != 1.3965) # 1


perennial.3 = scale(perennial.2 [,c(4:12)])
perennial.4 = cbind(perennial.2 [,c(1:3)],perennial.3)

priors <- c(prior(normal(0, 10), class = b))

perennial.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                        SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + (1|site_code) + (1|Taxon), 
                      family = gaussian(),
                      prior = priors,
                      data = perennial.4)

summary(perennial.model)


#### imputed traits NW ####

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.NW.final.csv",row.names = 1)

hist(imputed.traits$cover.change)
boxplot(imputed.traits$cover.change)
mean = mean(imputed.traits$cover.change, na.rm = TRUE)
std = sd(imputed.traits$cover.change, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$cover.change[which(imputed.traits$cover.change <Tmin | imputed.traits$cover.change > Tmax)])
# -79.47619 - -24.50000 and 23.52000 - 55.00000

imputed.traits.2 = imputed.traits %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)

imputed.traits.3 = scale(imputed.traits.2[,c(26:34)])
imputed.traits.4 = cbind(imputed.traits.2[,c(1:3)],imputed.traits.3)

priors <- c(prior(normal(0, 10), class = b))

imputed.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                             SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                           family = gaussian(),
                           prior = priors,
                           data = imputed.traits.4)

summary(imputed.traits.model)
# nothing significant

#### imputed traits NW outlier ####

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.NW.final.csv",row.names = 1)

imputed.traits.2 = imputed.traits %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)

imputed.traits.3 = imputed.traits.2 %>%
  filter(round(leafN.final, 5) < 45.18532) %>% # 4
  filter(round(height.final, 5) < 1.371600) %>% # 11
  filter(round(rootN.final, 5) < 27.18327) %>% # 18
  filter(round(SLA.final, 5) < 49.05439) %>% # 8
  filter(round(root.depth.final, 5) < 6.451600) %>% # 15
  filter(round(RTD.final, 5) < 0.6525000) %>% # 11
  filter(round(SRL.final, 5) < 437.3523) %>% # 21
  filter(round(rootDiam.final, 5) < 1.126012) %>% # 6
  filter(round(RMF.final, 5) > 0.04850791 & round(RMF.final, 5) < 0.76700839)

imputed.traits.4 = scale(imputed.traits.3[,c(26:34)])
imputed.traits.5 = cbind(imputed.traits.3[,c(1:3)],imputed.traits.4)

priors <- c(prior(normal(0, 10), class = b))

imputed.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                             SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                           family = gaussian(),
                           prior = priors,
                           data = imputed.traits.5)

summary(imputed.traits.model)
# nothing significant

#### imputed traits lifespan NW ####

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.NW.final.csv",row.names = 1)

imputed.traits.2 = imputed.traits %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)

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
# nothing significant

perennial.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                               SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                             family = gaussian(),
                             prior = priors,
                             data = perennial.imputed.3)

summary(perennial.traits.model)

#### impute traits functional group NW ####

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.NW.final.csv",row.names = 1)

imputed.traits.2 = imputed.traits %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)

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

#### imputed traits lifespan x functional group NW ####

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.NW.final.csv",row.names = 1)

imputed.traits.2 = imputed.traits %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)

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


#### imputed traits lifespan NW outlier ####

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.NW.final.csv",row.names = 1)

imputed.traits.2 = imputed.traits %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)

imputed.traits.3 = imputed.traits.2 %>%
  filter(round(leafN.final, 5) < 45.18532) %>% # 4
  filter(round(height.final, 5) < 1.371600) %>% # 11
  filter(round(rootN.final, 5) < 27.18327) %>% # 18
  filter(round(SLA.final, 5) < 49.05439) %>% # 8
  filter(round(root.depth.final, 5) < 6.451600) %>% # 15
  filter(round(RTD.final, 5) < 0.6525000) %>% # 11
  filter(round(SRL.final, 5) < 437.3523) %>% # 21
  filter(round(rootDiam.final, 5) < 1.126012) %>% # 6
  filter(round(RMF.final, 5) > 0.04850791 & round(RMF.final, 5) < 0.76700839) # 10

annual.imputed = imputed.traits.3 %>%
  filter(local_lifespan == "ANNUAL")
perennial.imputed = imputed.traits.3 %>%
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

perennial.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                               SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                             family = gaussian(),
                             prior = priors,
                             data = perennial.imputed.3)

summary(perennial.traits.model)

#### impute traits functional group NW outlier ####

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.NW.final.csv",row.names = 1)

imputed.traits.2 = imputed.traits %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)

imputed.traits.3 = imputed.traits.2 %>%
  filter(round(leafN.final, 5) < 45.18532) %>% # 4
  filter(round(height.final, 5) < 1.371600) %>% # 11
  filter(round(rootN.final, 5) < 27.18327) %>% # 18
  filter(round(SLA.final, 5) < 49.05439) %>% # 8
  filter(round(root.depth.final, 5) < 6.451600) %>% # 15
  filter(round(RTD.final, 5) < 0.6525000) %>% # 11
  filter(round(SRL.final, 5) < 437.3523) %>% # 21
  filter(round(rootDiam.final, 5) < 1.126012) %>% # 6
  filter(round(RMF.final, 5) > 0.04850791 & round(RMF.final, 5) < 0.76700839) # 10

forb.imputed = imputed.traits.3 %>%
  filter(functional_group == "FORB")
grass.imputed = imputed.traits.3 %>%
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

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.NW.final.csv",row.names = 1)

imputed.traits.2 = imputed.traits %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)

imputed.traits.3 = imputed.traits.2 %>%
  filter(round(leafN.final, 5) < 45.18532) %>% # 4
  filter(round(height.final, 5) < 1.371600) %>% # 11
  filter(round(rootN.final, 5) < 27.18327) %>% # 18
  filter(round(SLA.final, 5) < 49.05439) %>% # 8
  filter(round(root.depth.final, 5) < 6.451600) %>% # 15
  filter(round(RTD.final, 5) < 0.6525000) %>% # 11
  filter(round(SRL.final, 5) < 437.3523) %>% # 21
  filter(round(rootDiam.final, 5) < 1.126012) %>% # 6
  filter(round(RMF.final, 5) > 0.04850791 & round(RMF.final, 5) < 0.76700839) # 10

annual.forb.imputed = imputed.traits.3 %>%
  filter(functional_group == "FORB" & local_lifespan == "ANNUAL")
perennial.forb.imputed = imputed.traits.3 %>%
  filter(functional_group == "FORB" & local_lifespan == "PERENNIAL")
perennial.grass.imputed = imputed.traits.3 %>%
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
final.data.CC = read.csv("./Formatted.Data/Revisions/final.data.CC.NW.csv",row.names = 1)

hist(final.data.CC$cover.change)
boxplot(final.data.CC$cover.change)
mean = mean(final.data.CC$cover.change, na.rm = TRUE)
std = sd(final.data.CC$cover.change, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$cover.change[which(final.data.CC$cover.change <Tmin | final.data.CC$cover.change > Tmax)])
# -79.47619 - -38.00000 and 28.91667

final.data.CC.2 = final.data.CC %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) != 28.91667)

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

#### CC traits lifespan NW ####

final.data.CC = read.csv("./Formatted.Data/Revisions/final.data.CC.NW.csv",row.names = 1)

final.data.CC.2 = final.data.CC %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) != 28.91667)

table(final.data.CC$local_lifespan)

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

#### CC traits NW functional group ####

final.data.CC = read.csv("./Formatted.Data/Revisions/final.data.CC.NW.csv",row.names = 1)

final.data.CC.2 = final.data.CC %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) != 28.91667)

table(final.data.CC$functional_group)

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


#### CC traits NW outlier ####
final.data.CC = read.csv("./Formatted.Data/Revisions/final.data.CC.NW.csv",row.names = 1)

final.data.CC.2 = final.data.CC %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) != 28.91667)

final.data.CC.3 = final.data.CC.2 %>%
  filter(round(leafN.mg.g, 5) != 54.71597) %>% # 2
  filter(round(height.m, 5) < 1.800048) %>% # 3
  filter(round(rootN.mg.g, 5) < 31.17241) %>% # 5
  filter(round(root.depth_m, 5) < 2.426850) %>% # 7
  filter(round(SRL.groot.cahill.merge, 5) < 471.2364) # 5

final.data.CC.4 = scale(final.data.CC.3[,c(4:12)])
final.data.CC.5 = cbind(final.data.CC.3[,c(1:3)],final.data.CC.4)

# 236 data points estimating 11 variables

priors <- c(prior(normal(0, 10), class = b))


CC.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                 SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + (1|site_code) + (1|Taxon), 
               family = gaussian(),
               prior = priors,
               data = final.data.CC.5)

summary(CC.model)

#### CC traits lifespan NW outlier ####

final.data.CC = read.csv("./Formatted.Data/Revisions/final.data.CC.NW.csv",row.names = 1)

final.data.CC.2 = final.data.CC %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) != 28.91667)

table(final.data.CC$local_lifespan)

final.data.CC.3 = final.data.CC.2 %>%
  filter(round(leafN.mg.g, 5) != 54.71597) %>% # 2
  filter(round(height.m, 5) < 1.800048) %>% # 3
  filter(round(rootN.mg.g, 5) < 31.17241) %>% # 5
  filter(round(root.depth_m, 5) < 2.426850) %>% # 7
  filter(round(SRL.groot.cahill.merge, 5) < 471.2364) # 5


perennial = final.data.CC.3 %>%
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

final.data.CC = read.csv("./Formatted.Data/Revisions/final.data.CC.NW.csv",row.names = 1)

final.data.CC.2 = final.data.CC %>%
  filter(round(cover.change, 5) > -38.00000 & round(cover.change, 5) != 28.91667)

table(final.data.CC$functional_group)

final.data.CC.3 = final.data.CC.2 %>%
  filter(round(leafN.mg.g, 5) != 54.71597) %>% # 2
  filter(round(height.m, 5) < 1.800048) %>% # 3
  filter(round(rootN.mg.g, 5) < 31.17241) %>% # 5
  filter(round(root.depth_m, 5) < 2.426850) %>% # 7
  filter(round(SRL.groot.cahill.merge, 5) < 471.2364) # 5

forb = final.data.CC.3 %>%
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

