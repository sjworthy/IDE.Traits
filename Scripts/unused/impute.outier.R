#### imputed.traits.final ####
# model of trait data with complete cases for all species

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)

#### checking for trait outliers from trait data without woody species ####
# not going to remove outliers since need complete cases of traits for individuals to be included in analyses using lmer

hist(imputed.traits$leafN.final)
boxplot(imputed.traits$leafN.final)
mean = mean(imputed.traits$leafN.final, na.rm = TRUE)
std = sd(imputed.traits$leafN.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$leafN.final[which(imputed.traits$leafN.final <Tmin | imputed.traits$leafN.final > Tmax)])
# removed leafN 46.00560 54.71597 54.71597 81.93849
# percent removed 
table(is.na(imputed.traits$leafN.final))
958-368
2/590*100 #0.34%

hist(imputed.traits$height.final)
boxplot(imputed.traits$height.final)
mean = mean(imputed.traits$height.final, na.rm = TRUE)
std = sd(imputed.traits$height.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$height.final[which(imputed.traits$height.final <Tmin | imputed.traits$height.final > Tmax)])
# removed height 7.787850  9.457333 10.000000 10.000000 12.848467 19.050000 20.170667 20.501000 23.394667 27.101857 32.573656
# percent removed 
table(is.na(imputed.traits$height.final))
958-216
14/742*100 #1.89%

hist(imputed.traits$rootN.final)
boxplot(imputed.traits$rootN.final)
mean = mean(imputed.traits$rootN.final, na.rm = TRUE)
std = sd(imputed.traits$rootN.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$rootN.final[which(imputed.traits$rootN.final <Tmin | imputed.traits$rootN.final > Tmax)])
# remove rootN 25.55616 - 39.26274
# percent removed 
table(is.na(imputed.traits$rootN.final))
958-613
6/345*100 #1.74%

hist(imputed.traits$SLA.final)
boxplot(imputed.traits$SLA.final)
mean = mean(imputed.traits$SLA.final, na.rm = TRUE)
std = sd(imputed.traits$SLA.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$SLA.final[which(imputed.traits$SLA.final <Tmin | imputed.traits$SLA.final > Tmax)])
# remove SLA: 49.05439 53.50000 53.50000 57.70000 61.28000 63.13000 68.90000 68.90000
# percent removed 
table(is.na(imputed.traits$SLA.final))
958-244
8/714*100 #1.12%

hist(imputed.traits$root.depth.final)
boxplot(imputed.traits$root.depth.final)
mean = mean(imputed.traits$root.depth.final, na.rm = TRUE)
std = sd(imputed.traits$root.depth.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$root.depth.final[which(imputed.traits$root.depth.final <Tmin | imputed.traits$root.depth.final > Tmax)])
# remove depth 2.915000 - 6.451600
# percent removed 
table(is.na(imputed.traits$root.depth.final))
958-478
15/480*100 #3.13%

hist(imputed.traits$RTD.final)
boxplot(imputed.traits$RTD.final)
mean = mean(imputed.traits$RTD.final, na.rm = TRUE)
std = sd(imputed.traits$RTD.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$RTD.final[which(imputed.traits$RTD.final <Tmin | imputed.traits$RTD.final > Tmax)])
# remove RTD 0.7460162 - 1.4757011
# percent removed 
table(is.na(imputed.traits$RTD.final))
958-581
1/377*100 #0.27%

hist(imputed.traits$SRL.final)
boxplot(imputed.traits$SRL.final)
mean = mean(imputed.traits$SRL.final, na.rm = TRUE)
std = sd(imputed.traits$SRL.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$SRL.final[which(imputed.traits$SRL.final <Tmin | imputed.traits$SRL.final > Tmax)])
# remove SRL 437.3523-808.0298
# percent removed 
table(is.na(imputed.traits$SRL.final))
958-575
10/383*100 #2.61%

hist(imputed.traits$rootDiam.final)
boxplot(imputed.traits$rootDiam.final)
mean = mean(imputed.traits$rootDiam.final, na.rm = TRUE)
std = sd(imputed.traits$rootDiam.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$rootDiam.final[which(imputed.traits$rootDiam.final <Tmin | imputed.traits$rootDiam.final > Tmax)])
# remove diam 1.200000 1.396500 1.640000 1.850000 1.963333 4.830000
# percent removed 
table(is.na(imputed.traits$rootDiam.final))
958-556
4/402*100 #0.995%

hist(imputed.traits$RMF.final)
boxplot(imputed.traits$RMF.final)
mean = mean(imputed.traits$RMF.final, na.rm = TRUE)
std = sd(imputed.traits$RMF.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$RMF.final[which(imputed.traits$RMF.final <Tmin | imputed.traits$RMF.final > Tmax)])
# remove RMF: 0.7670084 0.7670084 0.7670084 0.7670084 0.7768355 0.7995028 0.7995028 0.8274671 0.8274671 0.8630137
# percent removed 
table(is.na(imputed.traits$RMF.final))
958-556
4/402*100 #0.995%


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


# look at correlations between traits
cor.traits.impute = cor(imputed.traits.2[,c(26:34)],use = "pairwise") 
corrplot(cor.traits.impute, method="number",tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")
# highest correlation is -0.31 and 0.35

imputed.traits.3 = scale(imputed.traits.2[,c(26:34)])
imputed.traits.4 = cbind(imputed.traits.2[,c(1:3)],imputed.traits.3)

# 907 data points estimating 11 variables

imputed.traits.model = lmer(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                              SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), data = imputed.traits.4)
summary(imputed.traits.model)
# nothing significant

plot(imputed.traits.model)
output.plot = allEffects(imputed.traits.model)
plot(output.plot)
hist(resid(imputed.traits.model))

# Add interaction between RTD and SRL

imputed.traits.model.interact = lmer(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                       SRL.final + RTD.final + RMF.final + SRL.final*RTD.final + (1|site_code) + (1|Taxon), data = imputed.traits.4)
summary(imputed.traits.model.interact)
# nothing is significant

#### imputed.traits.final lifespan ####

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)

table(imputed.traits$local_lifespan)

annual.imputed = imputed.traits %>%
  filter(local_lifespan == "ANNUAL")
perennial.imputed = imputed.traits %>%
  filter(local_lifespan == "PERENNIAL")


# 19 total removed because some removed with multiple traits
annual.imputed.2 = annual.imputed %>%
  filter(round(leafN.final, 5) < 46.00560) %>% # 4
  filter(round(height.final, 5) < 7.787850) %>% # 11
  filter(round(rootN.final, 5) < 25.55616) %>% # 18
  filter(round(SLA.final, 5) < 49.05439) %>% # 8
  filter(round(root.depth.final, 5) < 2.915000) %>% # 15
  filter(round(RTD.final, 5) < 0.7460162) %>% # 11
  filter(round(SRL.final, 5) < 437.3523) %>% # 21
  filter(round(rootDiam.final, 5) < 1.200000) %>% # 6
  filter(round(RMF.final, 5) < 0.7670084) # 10

annual.imputed.3 = scale(annual.imputed.2[,c(26:34)])
annual.imputed.4 = cbind(annual.imputed.2[,c(1:3)],annual.imputed.3)

annual.impute.model = lmer(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                             SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), data = annual.imputed.4)
summary(annual.impute.model)

# 74 total removed because some removed with multiple traits
perennial.imputed.2 = perennial.imputed %>%
  filter(round(leafN.final, 5) < 46.00560) %>% # 4
  filter(round(height.final, 5) < 7.787850) %>% # 11
  filter(round(rootN.final, 5) < 25.55616) %>% # 18
  filter(round(SLA.final, 5) < 49.05439) %>% # 8
  filter(round(root.depth.final, 5) < 2.915000) %>% # 15
  filter(round(RTD.final, 5) < 0.7460162) %>% # 11
  filter(round(SRL.final, 5) < 437.3523) %>% # 21
  filter(round(rootDiam.final, 5) < 1.200000) %>% # 6
  filter(round(RMF.final, 5) < 0.7670084) # 10

perennial.imputed.3 = scale(perennial.imputed.2[,c(26:34)])
perennial.imputed.4 = cbind(perennial.imputed.2[,c(1:3)],perennial.imputed.3)

perennial.impute.model = lmer(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), data = perennial.imputed.4)
summary(perennial.impute.model)
# RTD is significant

#### imputed.traits.final functional group ####

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)

table(imputed.traits$functional_group)

forb.imputed = imputed.traits %>%
  filter(functional_group == "FORB")
grass.imputed = imputed.traits %>%
  filter(functional_group == "GRASS")

# 29 total removed because some removed with multiple traits
forb.imputed.2 = forb.imputed %>%
  filter(round(leafN.final, 5) < 46.00560) %>% # 4
  filter(round(height.final, 5) < 7.787850) %>% # 11
  filter(round(rootN.final, 5) < 25.55616) %>% # 18
  filter(round(SLA.final, 5) < 49.05439) %>% # 8
  filter(round(root.depth.final, 5) < 2.915000) %>% # 15
  filter(round(RTD.final, 5) < 0.7460162) %>% # 11
  filter(round(SRL.final, 5) < 437.3523) %>% # 21
  filter(round(rootDiam.final, 5) < 1.200000) %>% # 6
  filter(round(RMF.final, 5) < 0.7670084) # 10

forb.imputed.3 = scale(forb.imputed.2[,c(26:34)])
forb.imputed.4 = cbind(forb.imputed.2[,c(1:3)],forb.imputed.3)

forb.impute.model = lmer(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                           SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), data = forb.imputed.4)
summary(forb.impute.model)

# 14 total removed because some removed with multiple traits
grass.imputed.2 = grass.imputed %>%
  filter(round(leafN.final, 5) < 46.00560) %>% # 4
  filter(round(height.final, 5) < 7.787850) %>% # 11
  filter(round(rootN.final, 5) < 25.55616) %>% # 18
  filter(round(SLA.final, 5) < 49.05439) %>% # 8
  filter(round(root.depth.final, 5) < 2.915000) %>% # 15
  filter(round(RTD.final, 5) < 0.7460162) %>% # 11
  filter(round(SRL.final, 5) < 437.3523) %>% # 21
  filter(round(rootDiam.final, 5) < 1.200000) %>% # 6
  filter(round(RMF.final, 5) < 0.7670084) # 10

grass.imputed.3 = scale(grass.imputed.2[,c(26:34)])
grass.imputed.4 = cbind(grass.imputed.2[,c(1:3)],grass.imputed.3)

grass.impute.model = lmer(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                            SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), data = grass.imputed.4)
summary(grass.impute.model)


#### imputed.traits.final lifespan x functional group ####

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)

table(imputed.traits$functional_group,imputed.traits$local_lifespan)
# annual forb
# perennial forb
# perennial grass

annual.forb.imputed = imputed.traits %>%
  filter(functional_group == "FORB" & local_lifespan == "ANNUAL")
perennial.forb.imputed = imputed.traits %>%
  filter(functional_group == "FORB" & local_lifespan == "PERENNIAL")
perennial.grass.imputed = imputed.traits %>%
  filter(functional_group == "GRASS" & local_lifespan == "PERENNIAL")


# 3 total removed because some removed with multiple traits
annual.forb.imputed.2 = annual.forb.imputed %>%
  filter(round(leafN.final, 5) < 46.00560) %>% # 4
  filter(round(height.final, 5) < 7.787850) %>% # 11
  filter(round(rootN.final, 5) < 25.55616) %>% # 18
  filter(round(SLA.final, 5) < 49.05439) %>% # 8
  filter(round(root.depth.final, 5) < 2.915000) %>% # 15
  filter(round(RTD.final, 5) < 0.7460162) %>% # 11
  filter(round(SRL.final, 5) < 437.3523) %>% # 21
  filter(round(rootDiam.final, 5) < 1.200000) %>% # 6
  filter(round(RMF.final, 5) < 0.7670084) # 10

annual.forb.imputed.3 = scale(annual.forb.imputed.2[,c(26:34)])
annual.forb.imputed.4 = cbind(annual.forb.imputed.2[,c(1:3)],annual.forb.imputed.3)

annual.forb.impute.model = lmer(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                  SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), data = annual.forb.imputed.4)
summary(annual.forb.impute.model)

# 25 removed
perennial.forb.imputed.2 = perennial.forb.imputed %>%
  filter(round(leafN.final, 5) < 46.00560) %>% # 4
  filter(round(height.final, 5) < 7.787850) %>% # 11
  filter(round(rootN.final, 5) < 25.55616) %>% # 18
  filter(round(SLA.final, 5) < 49.05439) %>% # 8
  filter(round(root.depth.final, 5) < 2.915000) %>% # 15
  filter(round(RTD.final, 5) < 0.7460162) %>% # 11
  filter(round(SRL.final, 5) < 437.3523) %>% # 21
  filter(round(rootDiam.final, 5) < 1.200000) %>% # 6
  filter(round(RMF.final, 5) < 0.7670084) # 10

perennial.forb.imputed.3 = scale(perennial.forb.imputed.2[,c(26:34)])
perennial.forb.imputed.4 = cbind(perennial.forb.imputed.2[,c(1:3)],perennial.forb.imputed.3)

perennial.forb.impute.model = lmer(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                     SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), data = perennial.forb.imputed.4)
summary(perennial.forb.impute.model)

# 4 removed
perennial.grass.imputed.2 = perennial.grass.imputed %>%
  filter(round(leafN.final, 5) < 46.00560) %>% # 4
  filter(round(height.final, 5) < 7.787850) %>% # 11
  filter(round(rootN.final, 5) < 25.55616) %>% # 18
  filter(round(SLA.final, 5) < 49.05439) %>% # 8
  filter(round(root.depth.final, 5) < 2.915000) %>% # 15
  filter(round(RTD.final, 5) < 0.7460162) %>% # 11
  filter(round(SRL.final, 5) < 437.3523) %>% # 21
  filter(round(rootDiam.final, 5) < 1.200000) %>% # 6
  filter(round(RMF.final, 5) < 0.7670084) # 10

perennial.grass.imputed.3 = scale(perennial.grass.imputed.2[,c(26:34)])
perennial.grass.imputed.4 = cbind(perennial.grass.imputed.2[,c(1:3)],perennial.grass.imputed.3)

perennial.grass.impute.model = lmer(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                      SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), data = perennial.grass.imputed.4)
summary(perennial.grass.impute.model)

### Bayes ####

priors <- c(prior(normal(0, 10), class = b))


imputed.traits.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                 SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
               family = gaussian(),
               prior = priors,
               data = imputed.traits.4)

summary(imputed.traits.model)
plot(imputed.traits.model, nvariables=6, ask=FALSE)
pp_check(imputed.traits.model, ndraws = 500)

annual.impute.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                            SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                          family = gaussian(),
                          prior = priors,
                          data = annual.imputed.4)

summary(annual.impute.model)

perennial.impute.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                            SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                          family = gaussian(),
                          prior = priors,
                          data = perennial.imputed.4)

summary(perennial.impute.model)

grass.impute.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                             SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                           family = gaussian(),
                           prior = priors,
                           data = grass.imputed.4)

forb.impute.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                           SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), 
                         family = gaussian(),
                         prior = priors,
                         data = forb.imputed.4)

#### checking for trait outliers from trait data without woody species ####
# not going to remove outliers since need complete cases of traits for individuals to be included in analyses using lmer

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.NW.final.csv", row.names = 1)

hist(imputed.traits$leafN.final)
boxplot(imputed.traits$leafN.final)
mean = mean(imputed.traits$leafN.final, na.rm = TRUE)
std = sd(imputed.traits$leafN.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$leafN.final[which(imputed.traits$leafN.final <Tmin | imputed.traits$leafN.final > Tmax)])
# removed leafN 45.18532 45.32000 46.00560 54.71597 54.71597
# percent removed 
table(is.na(imputed.traits$leafN.final))
958-368
2/590*100 #0.34%

hist(imputed.traits$height.final)
boxplot(imputed.traits$height.final)
mean = mean(imputed.traits$height.final, na.rm = TRUE)
std = sd(imputed.traits$height.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$height.final[which(imputed.traits$height.final <Tmin | imputed.traits$height.final > Tmax)])
# removed height 1.371600 1.383400 1.392886 1.439388 1.462500 1.493520 1.500000 1.750000 1.800000 1.800048 1.828800 1.996670 2.847733 3.800000
# percent removed 
table(is.na(imputed.traits$height.final))
958-216
14/742*100 #1.89%

hist(imputed.traits$rootN.final)
boxplot(imputed.traits$rootN.final)
mean = mean(imputed.traits$rootN.final, na.rm = TRUE)
std = sd(imputed.traits$rootN.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$rootN.final[which(imputed.traits$rootN.final <Tmin | imputed.traits$rootN.final > Tmax)])
# remove rootN 27.18327 - 39.39083
# percent removed 
table(is.na(imputed.traits$rootN.final))
958-613
6/345*100 #1.74%

hist(imputed.traits$SLA.final)
boxplot(imputed.traits$SLA.final)
mean = mean(imputed.traits$SLA.final, na.rm = TRUE)
std = sd(imputed.traits$SLA.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$SLA.final[which(imputed.traits$SLA.final <Tmin | imputed.traits$SLA.final > Tmax)])
# remove SLA: 49.05439 53.50000 53.50000 57.70000 61.28000 63.13000 68.90000 68.90000
# percent removed 
table(is.na(imputed.traits$SLA.final))
958-244
8/714*100 #1.12%

hist(imputed.traits$root.depth.final)
boxplot(imputed.traits$root.depth.final)
mean = mean(imputed.traits$root.depth.final, na.rm = TRUE)
std = sd(imputed.traits$root.depth.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$root.depth.final[which(imputed.traits$root.depth.final <Tmin | imputed.traits$root.depth.final > Tmax)])
# remove depth 6.451600 - 6.451600
# percent removed 
table(is.na(imputed.traits$root.depth.final))
958-478
15/480*100 #3.13%

hist(imputed.traits$RTD.final)
boxplot(imputed.traits$RTD.final)
mean = mean(imputed.traits$RTD.final, na.rm = TRUE)
std = sd(imputed.traits$RTD.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$RTD.final[which(imputed.traits$RTD.final <Tmin | imputed.traits$RTD.final > Tmax)])
# remove RTD 0.6525000 - 0.9708695
# percent removed 
table(is.na(imputed.traits$RTD.final))
958-581
1/377*100 #0.27%

hist(imputed.traits$SRL.final)
boxplot(imputed.traits$SRL.final)
mean = mean(imputed.traits$SRL.final, na.rm = TRUE)
std = sd(imputed.traits$SRL.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$SRL.final[which(imputed.traits$SRL.final <Tmin | imputed.traits$SRL.final > Tmax)])
# remove SRL 437.3523-760.8500
# percent removed 
table(is.na(imputed.traits$SRL.final))
958-575
10/383*100 #2.61%

hist(imputed.traits$rootDiam.final)
boxplot(imputed.traits$rootDiam.final)
mean = mean(imputed.traits$rootDiam.final, na.rm = TRUE)
std = sd(imputed.traits$rootDiam.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$rootDiam.final[which(imputed.traits$rootDiam.final <Tmin | imputed.traits$rootDiam.final > Tmax)])
# remove diam 1.126012 1.396500 1.640000 1.850000 1.963333 4.830000
# percent removed 
table(is.na(imputed.traits$rootDiam.final))
958-556
4/402*100 #0.995%

hist(imputed.traits$RMF.final)
boxplot(imputed.traits$RMF.final)
mean = mean(imputed.traits$RMF.final, na.rm = TRUE)
std = sd(imputed.traits$RMF.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$RMF.final[which(imputed.traits$RMF.final <Tmin | imputed.traits$RMF.final > Tmax)])
# remove RMF: 0.04850791- 0.8630137
# percent removed 
table(is.na(imputed.traits$RMF.final))
958-556
4/402*100 #0.995%


# 94 total removed because some removed with multiple traits
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






