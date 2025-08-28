## Script to evaluate linear mixed effects models of cover change and traits
# cover outliers removed

# load libraries 
library(tidyverse)
library(brms)
library(cowplot)

# https://imraugh.wordpress.com/2020/06/19/deconstructing-and-plotting-moderation-in-r/

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

# change forbs and legumes to just forbs
# change grass and graminoids to graminoids

imputed.NW$functional_group[imputed.NW$functional_group == "LEGUME"] <- "FORB"
imputed.NW$functional_group[imputed.NW$functional_group == "GRASS"] <- "GRAMINOID"

# remove all rows that are not annual or perennial
imputed.NW.2 = imputed.NW %>%
  filter(local_lifespan %in% c("ANNUAL","PERENNIAL"))
# make them factors
imputed.NW.2$local_lifespan = as.factor(imputed.NW.2$local_lifespan)
imputed.NW.2$functional_group = as.factor(imputed.NW.2$functional_group)

# subset for columns we need
imputed.NW.3 = imputed.NW.2 %>%
  select(cover.change,local_lifespan,functional_group,leafN.final,height.final,rootN.final,SLA.final,root.depth.final,
         rootDiam.final,SRL.final,RTD.final,RMF.final,mean.DSI,mean.MAP,site_code,Taxon) %>%
  drop_na()

#### Model with lifespan x functional group ####
# this model will included fewer data points since some populations
# were not annual or perennial
# 640 populations of 409 species compared to 661 populations of 421 species

priors <- c(prior(normal(0, 10), class = b))

imputed.traits.NW.all.cats.model = brm(cover.change ~ local_lifespan + functional_group + local_lifespan*functional_group + 
                                         leafN.final*local_lifespan*functional_group + height.final*local_lifespan*functional_group + 
                                         rootN.final*local_lifespan*functional_group + SLA.final*local_lifespan*functional_group +
                                         root.depth.final*local_lifespan*functional_group + rootDiam.final*local_lifespan*functional_group +
                                         SRL.final*local_lifespan*functional_group + RTD.final*local_lifespan*functional_group + 
                                         RMF.final*local_lifespan*functional_group + mean.DSI*local_lifespan*functional_group + 
                                         mean.MAP*local_lifespan*functional_group + (1|site_code) + (1|Taxon), 
                                       family = gaussian(),
                                       prior = priors,
                                       data = imputed.NW.3)

summary(imputed.traits.NW.all.cats.model)
bayes_R2(imputed.traits.NW.all.cats.model)
# 0.1451649 0.02516069 0.1009502 0.1995333
# saveRDS(imputed.traits.NW.all.cats.model, file = "./Results/all.cats.imputed.traits.no_woody.rds")

imputed.traits.NW.all.cats.model = readRDS("./Results/all.cats.imputed.traits.no_woody.rds")

# test for differences between combinations of lifespan and functional group for each trait
# also see if any particular group is significant for each trait
emtrends(imputed.traits.NW.all.cats.model, pairwise ~ local_lifespan*functional_group, var = "leafN.final")
# not significant

emtrends(imputed.traits.NW.all.cats.model, pairwise ~ local_lifespan*functional_group, var = "height.final")
# annual graminoid significant negative
# annual forbs and annual graminoids significantly differ
# perennial forbs and annual graminoids significantly differ

emtrends(imputed.traits.NW.all.cats.model, pairwise ~ local_lifespan*functional_group, var = "rootN.final")
# not significant

emtrends(imputed.traits.NW.all.cats.model, pairwise ~ local_lifespan*functional_group, var = "SLA.final")
# annual graminoid significant negative
# perennial forbs and annual graminoids significantly differ
# annual graminoids and perennial graminoids significantly differ

emtrends(imputed.traits.NW.all.cats.model, pairwise ~ local_lifespan*functional_group, var = "root.depth.final")
# not significant

emtrends(imputed.traits.NW.all.cats.model, pairwise ~ local_lifespan*functional_group, var = "rootDiam.final")
# not significant

emtrends(imputed.traits.NW.all.cats.model, pairwise ~ local_lifespan*functional_group, var = "SRL.final")
# not significant

emtrends(imputed.traits.NW.all.cats.model, pairwise ~ local_lifespan*functional_group, var = "RTD.final")
# perennial forb alone significant, positive
# contrasts not significant

emtrends(imputed.traits.NW.all.cats.model, pairwise ~ local_lifespan*functional_group, var = "RMF.final")
# not significant

emtrends(imputed.traits.NW.all.cats.model, pairwise ~ local_lifespan*functional_group, var = "mean.DSI")
# not significant

emtrends(imputed.traits.NW.all.cats.model, pairwise ~ local_lifespan*functional_group, var = "mean.MAP")
# not significant


### Model with lifespan groups ####

priors <- c(prior(normal(0, 10), class = b))

imputed.traits.NW.lifespan.model = brm(cover.change ~ local_lifespan + 
                                         leafN.final*local_lifespan + height.final*local_lifespan + 
                                         rootN.final*local_lifespan + SLA.final*local_lifespan +
                                         root.depth.final*local_lifespan + rootDiam.final*local_lifespan +
                                         SRL.final*local_lifespan + RTD.final*local_lifespan + 
                                         RMF.final*local_lifespan + mean.DSI*local_lifespan + 
                                         mean.MAP*local_lifespan + (1|site_code) + (1|Taxon), 
                                       family = gaussian(),
                                       prior = priors,
                                       data = imputed.NW.3)

summary(imputed.traits.NW.lifespan.model)
bayes_R2(imputed.traits.NW.lifespan.model)
# R2 0.09073373 0.02491132 0.04998194 0.1476335
# saveRDS(imputed.traits.NW.lifespan.model, file = "./Results/lifespan.cats.imputed.traits.no_woody.rds")

imputed.traits.NW.lifespan.model = readRDS("./Results/lifespan.cats.imputed.traits.no_woody.rds")

# test for differences between annuals and perennials for each trait
emtrends(imputed.traits.NW.lifespan.model, pairwise ~ local_lifespan, var = "leafN.final")
# not significant

emtrends(imputed.traits.NW.lifespan.model, pairwise ~ local_lifespan, var = "height.final")
# not significant

emtrends(imputed.traits.NW.lifespan.model, pairwise ~ local_lifespan, var = "rootN.final")
# not significant

emtrends(imputed.traits.NW.lifespan.model, pairwise ~ local_lifespan, var = "SLA.final")
# not significant

emtrends(imputed.traits.NW.lifespan.model, pairwise ~ local_lifespan, var = "root.depth.final")
# not significant

emtrends(imputed.traits.NW.lifespan.model, pairwise ~ local_lifespan, var = "rootDiam.final")
# not significant

emtrends(imputed.traits.NW.lifespan.model, pairwise ~ local_lifespan, var = "SRL.final")
# not significant

emtrends(imputed.traits.NW.lifespan.model, pairwise ~ local_lifespan, var = "RTD.final")
# not significant

emtrends(imputed.traits.NW.lifespan.model, pairwise ~ local_lifespan, var = "RMF.final")
# not significant

emtrends(imputed.traits.NW.lifespan.model, pairwise ~ local_lifespan, var = "mean.DSI")
# not significant

emtrends(imputed.traits.NW.lifespan.model, pairwise ~ local_lifespan, var = "mean.MAP")
# contrasts not significant
# annual significant, positive


#### Model with functional groups ####

imputed.NW.functional.group = imputed.NW %>%
  select(cover.change,functional_group,leafN.final,height.final,rootN.final,SLA.final,root.depth.final,
         rootDiam.final,SRL.final,RTD.final,RMF.final,mean.DSI,mean.MAP,site_code,Taxon) %>%
  drop_na()

priors <- c(prior(normal(0, 10), class = b))

imputed.traits.NW.functional.group.model = brm(cover.change ~ functional_group + 
                                                 leafN.final*functional_group + height.final*functional_group + 
                                                 rootN.final*functional_group + SLA.final*functional_group +
                                                 root.depth.final*functional_group + rootDiam.final*functional_group +
                                                 SRL.final*functional_group + RTD.final*functional_group + 
                                                 RMF.final*functional_group + mean.DSI*functional_group + 
                                                 mean.MAP*functional_group + (1|site_code) + (1|Taxon), 
                                               family = gaussian(),
                                               prior = priors,
                                               data = imputed.NW.functional.group)

summary(imputed.traits.NW.functional.group.model)
bayes_R2(imputed.traits.NW.functional.group.model)
#R2 0.09210356 0.02421696 0.052281 0.1472661
#saveRDS(imputed.traits.NW.functional.group.model, file = "./Results/functional.group.cats.imputed.traits.no_woody.rds")

imputed.traits.NW.functional.group.model = readRDS("./Results/functional.group.cats.imputed.traits.no_woody.rds")

# test for differences between annuals and perennials for each trait
emtrends(imputed.traits.NW.functional.group.model, pairwise ~ functional_group, var = "leafN.final")
# contrasts not significant
# forb significant, positive

emtrends(imputed.traits.NW.functional.group.model, pairwise ~ functional_group, var = "height.final")
# forb and graminoid significantly different

emtrends(imputed.traits.NW.functional.group.model, pairwise ~ functional_group, var = "rootN.final")
# not significant

emtrends(imputed.traits.NW.functional.group.model, pairwise ~ functional_group, var = "SLA.final")
# not significant

emtrends(imputed.traits.NW.functional.group.model, pairwise ~ functional_group, var = "root.depth.final")
# not significant

emtrends(imputed.traits.NW.functional.group.model, pairwise ~ functional_group, var = "rootDiam.final")
# not significant

emtrends(imputed.traits.NW.functional.group.model, pairwise ~ functional_group, var = "SRL.final")
# not significant

emtrends(imputed.traits.NW.functional.group.model, pairwise ~ functional_group, var = "RTD.final")
# not significant

emtrends(imputed.traits.NW.functional.group.model, pairwise ~ functional_group, var = "RMF.final")
# not significant

emtrends(imputed.traits.NW.functional.group.model, pairwise ~ functional_group, var = "mean.DSI")
# not significant

emtrends(imputed.traits.NW.functional.group.model, pairwise ~ functional_group, var = "mean.MAP")
# not significant

#### All data, no group models ####

imputed.NW = imputed.NW %>%
  select(cover.change,functional_group,leafN.final,height.final,rootN.final,SLA.final,root.depth.final,
         rootDiam.final,SRL.final,RTD.final,RMF.final,mean.DSI,mean.MAP,site_code,Taxon) %>%
  drop_na()

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

imputed.NW.traits.depth.leafN = brm(cover.change ~ root.depth.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                    family = gaussian(),
                                    prior = priors,
                                    data = imputed.NW)

summary(imputed.NW.traits.depth.leafN)
# leafN

saveRDS(imputed.NW.traits.depth.leafN, file = "./Results/all.imputed.NW.traits.depth.leafN.rds")
bayes_R2(imputed.NW.traits.depth.leafN)
# 0.04902346 0.02481235 0.01394253 0.1094119

imputed.NW.traits.RTD.SRL = brm(cover.change ~ RTD.final*SRL.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                family = gaussian(),
                                prior = priors,
                                data = imputed.NW)

summary(imputed.NW.traits.RTD.SRL)
# nothing significant

saveRDS(imputed.NW.traits.RTD.SRL, file = "./Results/all.imputed.NW.traits.RTD.SRL.rds")
bayes_R2(imputed.NW.traits.RTD.SRL)
# R2 0.04620861 0.02366104 0.01308091 0.1012018

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



#### Lifespan interaction models ####

priors <- c(prior(normal(0, 10), class = b))
imputed.NW.traits.height.leafN.lifespan = brm(cover.change ~ local_lifespan + height.final*leafN.final*local_lifespan + 
                                                mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                              family = gaussian(),
                                              prior = priors,
                                              data = imputed.NW.3)

summary(imputed.NW.traits.height.leafN.lifespan)
bayes_R2(imputed.NW.traits.height.leafN.lifespan)
# R2 0.05374738 0.02257641 0.01989944 0.1089797
# saveRDS(imputed.NW.traits.height.leafN.lifespan, file = "./Results/imputed.traits.NW.height.leafN.lifespan.rds")

imputed.NW.traits.height.leafN.lifespan = readRDS("./Results/imputed.traits.NW.height.leafN.lifespan.rds")

# test for differences between annuals and perennials for each trait
emtrends(imputed.NW.traits.height.leafN.lifespan, 
         pairwise ~ local_lifespan | leafN.final,
         var = "height.final",
         at = list(leafN.final = c(mean(imputed.NW.3$leafN.final),min(imputed.NW.3$leafN.final),
                                   max(imputed.NW.3$leafN.final))))
# annual and perennial not significantly different

conditional_effects(imputed.NW.traits.height.leafN.lifespan,
                    effects = "leafN.final:height.final",
                    conditions = data.frame(local_lifespan = c("ANNUAL", "PERENNIAL")))
# positive cover change with taller height and higher leafN

imputed.NW.traits.depth.leafN.lifespan = brm(cover.change ~ local_lifespan + root.depth.final*leafN.final*local_lifespan + 
                                               mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                             family = gaussian(),
                                             prior = priors,
                                             data = imputed.NW.3)

summary(imputed.NW.traits.depth.leafN.lifespan)
# saveRDS(imputed.NW.traits.depth.leafN.lifespan, file = "./Results/imputed.NW.traits.depth.leafN.lifespan.rds")
bayes_R2(imputed.NW.traits.depth.leafN.lifespan)
# R2 0.05340336 0.02401003 0.01838627 0.1112438

imputed.NW.traits.depth.leafN.lifespan = readRDS("./Results/imputed.NW.traits.depth.leafN.lifespan.rds")

# test for differences between annuals and perennials for each trait
emtrends(imputed.NW.traits.depth.leafN.lifespan, 
         pairwise ~ local_lifespan | leafN.final,
         var = "root.depth.final",
         at = list(leafN.final = c(mean(imputed.NW.3$leafN.final),min(imputed.NW.3$leafN.final),
                                   max(imputed.NW.3$leafN.final))))
# not significant

imputed.NW.traits.RTD.SRL.lifespan = brm(cover.change ~ local_lifespan + RTD.final*SRL.final*local_lifespan +
                                           mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                         family = gaussian(),
                                         prior = priors,
                                         data = imputed.NW.3)

summary(imputed.NW.traits.RTD.SRL.lifespan)
# saveRDS(imputed.NW.traits.RTD.SRL.lifespan, file = "./Results/imputed.NW.traits.RTD.SRL.lifespan.rds")
bayes_R2(imputed.NW.traits.RTD.SRL.lifespan)
# R2 0.05214755 0.02356031 0.0186025 0.1091153

imputed.NW.traits.RTD.SRL.lifespan = readRDS("./Results/imputed.NW.traits.RTD.SRL.lifespan.rds")

# test for differences between annuals and perennials for each trait
emtrends(imputed.NW.traits.RTD.SRL.lifespan, 
               pairwise ~ local_lifespan | SRL.final,
               var = "RTD.final",
               at = list(SRL.final = c(mean(imputed.NW.3$SRL.final),min(imputed.NW.3$SRL.final),
                                       max(imputed.NW.3$SRL.final))))
# lifespans not different

conditional_effects(imputed.NW.traits.RTD.SRL.lifespan,
                    effects = "RTD.final:SRL.final",
                    conditions = data.frame(local_lifespan = c("ANNUAL", "PERENNIAL")))

imputed.NW.traits.leafN.RMF.lifespan = brm(cover.change ~ local_lifespan + leafN.final*RMF.final*local_lifespan +
                                             mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                           family = gaussian(),
                                           prior = priors,
                                           data = imputed.NW.3)

summary(imputed.NW.traits.leafN.RMF.lifespan)
# saveRDS(imputed.NW.traits.leafN.RMF.lifespan, file = "./Results/imputed.NW.traits.leafN.RMF.lifespan.rds")
bayes_R2(imputed.NW.traits.leafN.RMF.lifespan)
# R2 0.05891678 0.02582159 0.02123988 0.1184489

imputed.NW.traits.leafN.RMF.lifespan = readRDS("./Results/imputed.NW.traits.leafN.RMF.lifespan.rds")

# test for differences between annuals and perennials for each trait
emtrends(imputed.NW.traits.leafN.RMF.lifespan, 
               pairwise ~ local_lifespan | leafN.final,
               var = "RMF.final",
               at = list(leafN.final = c(mean(imputed.NW.3$leafN.final),min(imputed.NW.3$leafN.final),
                                         max(imputed.NW.3$leafN.final))))
# lifespans don't differ

## Environment Interactions

imputed.NW.traits.height.MAP.DSI.lifespan = brm(cover.change ~ local_lifespan + height.final*mean.MAP*local_lifespan +
                                                  height.final*mean.DSI*local_lifespan + (1|site_code) + (1|Taxon), 
                                                family = gaussian(),
                                                prior = priors,
                                                data = imputed.NW.3)

summary(imputed.NW.traits.height.MAP.DSI.lifespan)
bayes_R2(imputed.NW.traits.height.MAP.DSI.lifespan)
# R2 0.0612661 0.02346778 0.02552157 0.1157723
#saveRDS(imputed.NW.traits.height.MAP.DSI.lifespan, file = "./Results/imputed.NW.traits.height.MAP.DSI.lifespan.rds")

imputed.NW.traits.height.MAP.DSI.lifespan = readRDS("./Results/imputed.NW.traits.height.MAP.DSI.lifespan.rds")

# test for differences between annuals and perennials for each trait
emtrends(imputed.NW.traits.height.MAP.DSI.lifespan, 
               pairwise ~ local_lifespan | height.final,
               var = "mean.MAP",
               at = list(height.final = c(mean(imputed.NW.3$height.final),min(imputed.NW.3$height.final),
                                          max(imputed.NW.3$height.final))))
# significant different in MAP x height relationship between annuals and perennials
# annuals are taller with higher MAP, significant

conditional_effects(imputed.NW.traits.height.MAP.DSI.lifespan,
                    effects = "height.final:mean.MAP",
                    conditions = data.frame(local_lifespan = c("ANNUAL", "PERENNIAL")))
# positive cover change with taller height and higher leafN

emtrends(imputed.NW.traits.height.MAP.DSI.lifespan, 
               pairwise ~ local_lifespan | height.final,
               var = "mean.DSI",
               at = list(height.final = c(mean(imputed.NW.3$height.final),min(imputed.NW.3$height.final),
                                          max(imputed.NW.3$height.final))))
# not significant

imputed.NW.traits.leafN.MAP.DSI.lifespan = brm(cover.change ~ local_lifespan + leafN.final*mean.MAP*local_lifespan +
                                                 leafN.final*mean.DSI*local_lifespan + (1|site_code) + (1|Taxon), 
                                               family = gaussian(),
                                               prior = priors,
                                               data = imputed.NW.3)

summary(imputed.NW.traits.leafN.MAP.DSI.lifespan)
# leafN x DSI

# saveRDS(imputed.NW.traits.leafN.MAP.DSI.lifespan, file = "./Results/imputed.NW.traits.leafN.MAP.DSI.lifespan.rds")
bayes_R2(imputed.NW.traits.leafN.MAP.DSI.lifespan)
# R2 0.07813172 0.02717293 0.03685992 0.1424638

imputed.NW.traits.leafN.MAP.DSI.lifespan = readRDS("./Results/imputed.NW.traits.leafN.MAP.DSI.lifespan.rds")

emtrends(imputed.NW.traits.leafN.MAP.DSI.lifespan, 
               pairwise ~ local_lifespan | leafN.final,
               var = "mean.MAP",
               at = list(leafN.final = c(mean(imputed.NW.3$leafN.final),min(imputed.NW.3$leafN.final),
                                          max(imputed.NW.3$leafN.final))))
# not different

emtrends(imputed.NW.traits.leafN.MAP.DSI.lifespan, 
               pairwise ~ local_lifespan | leafN.final,
               var = "mean.DSI",
               at = list(leafN.final = c(mean(imputed.NW.3$leafN.final),min(imputed.NW.3$leafN.final),
                                         max(imputed.NW.3$leafN.final))))
# not different
# annuals are significant

conditional_effects(imputed.NW.traits.leafN.MAP.DSI.lifespan,
                    effects = "leafN.final:mean.DSI",
                    conditions = data.frame(local_lifespan = c("ANNUAL", "PERENNIAL")))
# positive cover change with higher leafN in less severe drought

imputed.NW.traits.root.depth.MAP.DSI.lifespan = brm(cover.change ~ local_lifespan + root.depth.final*mean.MAP*local_lifespan +
                                                      root.depth.final*mean.DSI*local_lifespan + (1|site_code) + (1|Taxon), 
                                           family = gaussian(),
                                           prior = priors,
                                           data = imputed.NW.3)

summary(imputed.NW.traits.root.depth.MAP.DSI.lifespan)
# MAP
#saveRDS(imputed.NW.traits.root.depth.MAP.DSI.lifespan, file = "./Results/imputed.NW.traits.depth.MAP.DSI.lifespan.rds")

emtrends(imputed.NW.traits.root.depth.MAP.DSI.lifespan, 
               pairwise ~ local_lifespan | root.depth.final,
               var = "mean.MAP",
               at = list(root.depth.final = c(mean(imputed.NW.3$root.depth.final),min(imputed.NW.3$root.depth.final),
                                         max(imputed.NW.3$root.depth.final))))
# annuals are significant only at intermediate root depth

conditional_effects(imputed.NW.traits.root.depth.MAP.DSI.lifespan,
                    effects = "root.depth.final:mean.MAP",
                    conditions = data.frame(local_lifespan = c("ANNUAL", "PERENNIAL")))

emtrends(imputed.NW.traits.root.depth.MAP.DSI.lifespan, 
               pairwise ~ local_lifespan | root.depth.final,
               var = "mean.DSI",
               at = list(root.depth.final = c(mean(imputed.NW.3$root.depth.final),min(imputed.NW.3$root.depth.final),
                                         max(imputed.NW.3$root.depth.final))))
# not different

imputed.NW.traits.RTD.MAP.DSI.lifespan = brm(cover.change ~ local_lifespan + RTD.final*mean.MAP*local_lifespan + 
                                               RTD.final*mean.DSI*local_lifespan + (1|site_code) + (1|Taxon), 
                                    family = gaussian(),
                                    prior = priors,
                                    data = imputed.NW.3)

summary(imputed.NW.traits.RTD.MAP.DSI.lifespan)
#saveRDS(imputed.NW.traits.RTD.MAP.DSI.lifespan, file = "./Results/imputed.NW.traits.RTD.MAP.DSI.lifespan.rds")

emtrends(imputed.NW.traits.RTD.MAP.DSI.lifespan, 
               pairwise ~ local_lifespan | RTD.final,
               var = "mean.MAP",
               at = list(RTD.final = c(mean(imputed.NW.3$RTD.final),min(imputed.NW.3$RTD.final),
                                              max(imputed.NW.3$RTD.final))))
# annual significant at low values

conditional_effects(imputed.NW.traits.RTD.MAP.DSI.lifespan,
                    effects = "RTD.final:mean.MAP",
                    conditions = data.frame(local_lifespan = c("ANNUAL", "PERENNIAL")))

emtrends(imputed.NW.traits.RTD.MAP.DSI.lifespan, 
               pairwise ~ local_lifespan | RTD.final,
               var = "mean.DSI",
               at = list(RTD.final = c(mean(imputed.NW.3$RTD.final),min(imputed.NW.3$RTD.final),
                                       max(imputed.NW.3$RTD.final))))
# not different

imputed.NW.traits.SRL.MAP.DSI.lifespan = brm(cover.change ~ local_lifespan + SRL.final*mean.MAP*local_lifespan + 
                                               SRL.final*mean.DSI*local_lifespan + (1|site_code) + (1|Taxon), 
                                    family = gaussian(),
                                    prior = priors,
                                    data = imputed.NW.3)

summary(imputed.NW.traits.SRL.MAP.DSI.lifespan)
#saveRDS(imputed.NW.traits.SRL.MAP.DSI.lifespan, file = "./Results/imputed.NW.traits.SRL.MAP.DSI.lifespan.rds")

emtrends(imputed.NW.traits.SRL.MAP.DSI.lifespan, 
         pairwise ~ local_lifespan | SRL.final,
         var = "mean.MAP",
         at = list(SRL.final = c(mean(imputed.NW.3$SRL.final),min(imputed.NW.3$SRL.final),
                                 max(imputed.NW.3$SRL.final))))
# annual significant at intermediate values

conditional_effects(imputed.NW.traits.SRL.MAP.DSI.lifespan,
                    effects = "SRL.final:mean.MAP",
                    conditions = data.frame(local_lifespan = c("ANNUAL", "PERENNIAL")))

emtrends(imputed.NW.traits.SRL.MAP.DSI.lifespan, 
         pairwise ~ local_lifespan | SRL.final,
         var = "mean.DSI",
         at = list(SRL.final = c(mean(imputed.NW.3$SRL.final),min(imputed.NW.3$SRL.final),
                                 max(imputed.NW.3$SRL.final))))
# not significnat

imputed.NW.traits.RMF.MAP.DSI.lifespan = brm(cover.change ~ local_lifespan + RMF.final*mean.MAP*local_lifespan + 
                                               RMF.final*mean.DSI*local_lifespan + (1|site_code) + (1|Taxon), 
                                    family = gaussian(),
                                    prior = priors,
                                    data = imputed.NW.3)

summary(imputed.NW.traits.RMF.MAP.DSI.lifespan)
saveRDS(imputed.NW.traits.RMF.MAP.DSI.lifespan, file = "./Results/imputed.NW.traits.RMF.MAP.DSI.lifespan.rds")

emtrends(imputed.NW.traits.RMF.MAP.DSI.lifespan, 
         pairwise ~ local_lifespan | RMF.final,
         var = "mean.MAP",
         at = list(RMF.final = c(mean(imputed.NW.3$RMF.final),min(imputed.NW.3$RMF.final),
                                 max(imputed.NW.3$RMF.final))))
# annual significant at intermediate values

conditional_effects(imputed.NW.traits.RTD.MAP.DSI.lifespan,
                    effects = "RTD.final:mean.MAP",
                    conditions = data.frame(local_lifespan = c("ANNUAL", "PERENNIAL")))

emtrends(imputed.NW.traits.RMF.MAP.DSI.lifespan, 
         pairwise ~ local_lifespan | RMF.final,
         var = "mean.DSI",
         at = list(RMF.final = c(mean(imputed.NW.3$RMF.final),min(imputed.NW.3$RMF.final),
                                 max(imputed.NW.3$RMF.final))))
# annual and perennials significantly differ at hight RMF values

