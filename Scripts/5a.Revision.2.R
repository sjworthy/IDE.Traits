## Script to evaluate linear mixed effects models of cover change and traits
# cover outliers removed

# load libraries 
library(tidyverse)
library(brms)
library(emmeans)
library(cowplot)
library(stringr)


# https://imraugh.wordpress.com/2020/06/19/deconstructing-and-plotting-moderation-in-r/

#### Load imputed data without woody species ####

# determine outlier cover change 
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

# read in environmental data
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

# add relative BACI column to other data - for revisions round 2
all.cover = read.csv("./Formatted.Data/Revisions/BACI.data.final.csv") %>%
  select(site_code,Taxon,relative.cover.change)

# remove whitespace to help merging
imputed.NW$site_code <- str_trim(imputed.NW$site_code)
all.cover$site_code <- str_trim(all.cover$site_code)

# join together
imputed.NW = left_join(imputed.NW,all.cover)

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

#### All species imputed traits model NW ####

priors <- c(prior(normal(0, 10), class = b))

imputed.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                              family = gaussian(),
                              prior = priors,
                              data = imputed.NW)

summary(imputed.traits.NW.model)
# leafN significant positive
#saveRDS(imputed.traits.NW.model, file = "./Results/all.imputed.traits.no_woody.rds")
bayes_R2(imputed.traits.NW.model)
# R2 0.06008998 0.02492269 0.02400752 0.1181871

imputed.traits.NW.model = readRDS("./Results/all.imputed.traits.no_woody.rds")


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

#### Model with lifespan x functional group ####
# this model will included fewer data points since some populations
# were not annual or perennial
# 640 populations of 409 species compared to 661 populations of 421 species

# remove annual graminoids
imputed.NW.4 = imputed.NW.3 %>%
  filter(!(local_lifespan %in% c("ANNUAL") & functional_group %in% c("GRAMINOID")))
# group names together
imputed.NW.4$all.cats = paste(imputed.NW.4$local_lifespan, imputed.NW.4$functional_group, sep = "_")

priors <- c(prior(normal(0, 10), class = b))

imputed.traits.NW.all.cats.model = brm(cover.change ~ all.cats + 
                                         leafN.final*all.cats + height.final*all.cats + 
                                         rootN.final*all.cats+ SLA.final*all.cats +
                                         root.depth.final*all.cats + rootDiam.final*all.cats +
                                         SRL.final*all.cats + RTD.final*all.cats + 
                                         RMF.final*all.cats + mean.DSI*all.cats + 
                                         mean.MAP*all.cats + (1|site_code) + (1|Taxon), 
                                       family = gaussian(),
                                       prior = priors,
                                       data = imputed.NW.4)

summary(imputed.traits.NW.all.cats.model)
bayes_R2(imputed.traits.NW.all.cats.model)
# 0.1317751 0.02920505 0.08324678 0.1985978
#saveRDS(imputed.traits.NW.all.cats.model, file = "./Results/all.cats.imputed.traits.no_woody.rds")

imputed.traits.NW.all.cats.model = readRDS("./Results/all.cats.imputed.traits.no_woody.rds")

# checking model assumptions
pp_check(imputed.traits.NW.all.cats.model, type = "error_hist") # histogram of residuals
pp_check(imputed.traits.NW.all.cats.model, type = "error_scatter_avg") # # fitted vs residuals
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.traits.NW.all.cats.model))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
fitted_vals <- as.data.frame(fitted(imputed.traits.NW.all.cats.model, summary = TRUE))
plot(fitted_vals$Estimate, model_residuals$Estimate,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)


# test for differences between combinations of lifespan and functional group for each trait
# also see if any particular group is significant for each trait
emtrends(imputed.traits.NW.all.cats.model, pairwise ~ all.cats, var = "leafN.final")
# not significant

emtrends(imputed.traits.NW.all.cats.model, pairwise ~ all.cats, var = "height.final")
# not significant

emtrends(imputed.traits.NW.all.cats.model, pairwise ~ all.cats, var = "rootN.final")
# not significant

emtrends(imputed.traits.NW.all.cats.model, pairwise ~ all.cats, var = "SLA.final")
# not significant

emtrends(imputed.traits.NW.all.cats.model, pairwise ~ all.cats, var = "root.depth.final")
# not significant

emtrends(imputed.traits.NW.all.cats.model, pairwise ~ all.cats, var = "rootDiam.final")
# not significant

emtrends(imputed.traits.NW.all.cats.model, pairwise ~ all.cats, var = "SRL.final")
# not significant

emtrends(imputed.traits.NW.all.cats.model, pairwise ~ all.cats, var = "RTD.final")
# perennial forb alone significant, positive
# contrasts not significant

emtrends(imputed.traits.NW.all.cats.model, pairwise ~ all.cats, var = "RMF.final")
# not significant

emtrends(imputed.traits.NW.all.cats.model, pairwise ~ all.cats, var = "mean.DSI")
# not significant

emtrends(imputed.traits.NW.all.cats.model, pairwise ~ all.cats, var = "mean.MAP")
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

# checking model assumptions
pp_check(imputed.traits.NW.lifespan.model)
pp_check(imputed.traits.NW.lifespan.model, type = "error_scatter_avg")


instasummary(imputed.traits.NW.lifespan.model)
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

imputed.NW.functional.group = imputed.NW %>%
  select(cover.change,functional_group,leafN.final,height.final,rootN.final,SLA.final,root.depth.final,
         rootDiam.final,SRL.final,RTD.final,RMF.final,mean.DSI,mean.MAP,site_code,Taxon,relative.cover.change
         ) %>%
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



test = brm(cover.change ~ functional_group + 
             leafN.final*functional_group + height.final*functional_group + 
             RTD.final*functional_group + (1|site_code) + (1|Taxon), 
           family = gaussian(),
           prior = priors,
           data = imputed.NW.functional.group)

summary(imputed.traits.NW.functional.group.model)
bayesa_R2(imputed.traits.NW.functional.group.model)
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

#### Plotting Functional Group Interactions ####
imputed.NW.functional.group = imputed.NW %>%
  select(cover.change,functional_group,leafN.final,height.final,rootN.final,SLA.final,root.depth.final,
         rootDiam.final,SRL.final,RTD.final,RMF.final,mean.DSI,mean.MAP,site_code,Taxon) %>%
  drop_na()

imputed.traits.NW.functional.group.model = readRDS("./Results/functional.group.cats.imputed.traits.no_woody.rds")

height.functional.group.effect = conditional_effects(imputed.traits.NW.functional.group.model, 
                                                     effects = "height.final:functional_group")$`height.final:functional_group`

attr(imputed.NW.functional.group$height.final, "scaled:scale")
attr(imputed.NW.functional.group$height.final, "scaled:center")

height.functional.group.effect$height.bt = (height.functional.group.effect$height.final*0.2353906) + 0.3476824

height.x.functional.group.plot = ggplot(data = height.functional.group.effect, aes(x = height.bt, y = estimate__, group = as.factor(effect2__))) +
  geom_line(aes(color = as.factor(effect2__)), size = 1.5,show.legend = FALSE) +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = as.factor(effect2__)), alpha = 0.5,show.legend = FALSE) +  
  geom_line(aes(y = upper__, color = factor(effect2__)), size = 1.2,show.legend = FALSE) +
  geom_line(aes(y = lower__, color = factor(effect2__)), size = 1.2,show.legend = FALSE) + 
  labs(x = "Height (m)", y = "Cover Change (%)", color = "Functional Group") +
  scale_colour_manual(values = c("#F17236", "#6E687E"), labels = c("Forbs","Graminoids"))+
  scale_fill_manual(values = c("#F17236", "#6E687E"))+
  theme_classic()
height.x.functional.group.plot

#ggsave("./Plots/Revision.2/height.functional.group.pdf", height = 3, width = 3)
#ggsave("./Plots/Revision.2/height.functional.group.legend.pdf", height = 3, width = 3)

##### Plotting #####
#### ANNUAL Backtransform ####

# load lifespan model
annual.traits.NW.model = readRDS("./Results/lifespan.cats.imputed.traits.no_woody.rds")

attr(imputed.NW.3$mean.MAP, "scaled:scale")
attr(imputed.NW.3$mean.MAP, "scaled:center")

imputed.NW.3$leafN.bt =(imputed.NW.3$leafN.final*7.114868) + 21.16957
imputed.NW.3$height.bt = (imputed.NW.3$height.final*0.2353906) + 0.3476824
imputed.NW.3$rootN.bt = (imputed.NW.3$rootN.final*4.488773) + 9.855027
imputed.NW.3$SLA.bt = (imputed.NW.3$SLA.final*8.581554) + 19.93883
imputed.NW.3$Depth.bt = (imputed.NW.3$root.depth.final*0.5333377) + 0.5784653
imputed.NW.3$Diam.bt = (imputed.NW.3$rootDiam.final*0.1598629) + 0.3501538
imputed.NW.3$SRL.bt = (imputed.NW.3$SRL.final*74.26099) + 106.4204
imputed.NW.3$RTD.bt = (imputed.NW.3$RTD.final*0.1159369) + 0.2349392
imputed.NW.3$RMF.bt = (imputed.NW.3$RMF.final*0.1058164) + 0.4007788
imputed.NW.3$DSI.bt = (imputed.NW.3$mean.DSI*0.1950438) + -0.4604709
imputed.NW.3$MAP.bt = (imputed.NW.3$mean.MAP*487.1451) + 691.6107

imputed.NW.3.annual = imputed.NW.3 %>%
  filter(local_lifespan == "ANNUAL")

annual.effects = conditional_effects(annual.traits.NW.model)

#### Annual leafN ####
leafN.effects = annual.effects$`leafN.final:local_lifespan`
leafN.effects = leafN.effects %>%
  filter(local_lifespan == "ANNUAL")

leafN.effects$leafN.bt= (leafN.effects$leafN.final*7.114868) + 21.16957

x.value = c(-2,-1,0,1,2,3)
(x.value*7.114868) + 21.16957

annual.leafN.plot = ggplot() +
  geom_point(data = imputed.NW.3.annual, aes(x = leafN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Leaf N (mg/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  xlim(2.68951,44.19001)+
  theme_classic()
annual.leafN.plot

ggsave("./Plots/Revision.2/annual.traits.NW.leafN.pdf", height = 3, width = 3)

#### Annual height ####

height.effects = annual.effects$`height.final:local_lifespan`
height.effects = height.effects %>%
  filter(local_lifespan == "ANNUAL")
height.effects$height.bt= (height.effects$height.final*0.2353906) + 0.3476824

x.value = c(-1,0,1,2,3,4)
(x.value*0.2353906) + 0.3476824

annual.height.plot = ggplot() +
  geom_point(data = imputed.NW.3.annual, aes(x = height.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = height.effects, aes(x = height.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Height (m)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.009000064,1.322333169)+
  theme_classic()
annual.height.plot

ggsave("./Plots/Revision.2/annual.traits.NW.height.pdf", height = 3, width = 3)


#### Annual RootN ####
rootN.effects = annual.effects$`rootN.final:local_lifespan`
rootN.effects = rootN.effects %>%
  filter(local_lifespan == "ANNUAL")
rootN.effects$rootN.bt= (rootN.effects$rootN.final*4.488773) + 9.855027

x.value = c(-2,-1,0,1,2,3)
(x.value*4.488773) + 9.855027

annual.rootN.plot = ggplot() +
  geom_point(data = imputed.NW.3.annual, aes(x = rootN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Root N (mg/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.3800000,25.6283268)+
  theme_classic()
annual.rootN.plot

ggsave("./Plots/Revision.2/annual.traits.NW.rootN.pdf", height = 3, width = 3)

#### Annual SLA ####
SLA.effects = annual.effects$`SLA.final:local_lifespan`
SLA.effects = SLA.effects %>%
  filter(local_lifespan == "ANNUAL")
SLA.effects$SLA.bt= (SLA.effects$SLA.final*8.581554) + 19.93883

x.value = c(-2,-1,0,1,2,3)
(x.value*8.581554) + 19.93883

annual.SLA.plot = ggplot() +
  geom_point(data = imputed.NW.3.annual, aes(x = SLA.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = expression("Specific Leaf Area (m"^2*"/kg)"), y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.617596,46.320001)+
  theme_classic()
annual.SLA.plot

ggsave("./Plots/Revision.2/annual.traits.NW.SLA.pdf", height = 3, width = 3)

#### Annual depth ####
root.depth.effects = annual.effects$`root.depth.final:local_lifespan`
root.depth.effects = root.depth.effects %>%
  filter(local_lifespan == "ANNUAL")
root.depth.effects$Depth.bt= (root.depth.effects$root.depth.final*0.5333377) + 0.5784653

x.value = c(0,1,2,3,4,5)
(x.value*0.5333377) + 0.5784653

annual.root.depth.plot = ggplot() +
  geom_point(data = imputed.NW.3.annual, aes(x = Depth.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = Depth.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = Depth.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Rooting Depth (m)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.0479999,3.709335)+
  theme_classic()
annual.root.depth.plot

ggsave("./Plots/Revision.2/annual.traits.NW.root.depth.pdf", height = 3, width = 3)

#### Annual Diam ####

rootDiam.effects = annual.effects$`rootDiam.final:local_lifespan`
rootDiam.effects = rootDiam.effects %>%
  filter(local_lifespan == "ANNUAL")
rootDiam.effects$Diam.bt= (rootDiam.effects$rootDiam.final*0.1598629) + 0.3501538

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1598629) + 0.3501538

annual.rootDiam.plot = ggplot() +
  geom_point(data = imputed.NW.3.annual, aes(x = Diam.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = Diam.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = Diam.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Root Diameter (mm)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.07840002,1.12601235)+
  theme_classic()
annual.rootDiam.plot

ggsave("./Plots/Revision.2/annual.traits.NW.rootDiam.pdf", height = 3, width = 3)

#### Annual SRL ####

SRL.effects = annual.effects$`SRL.final:local_lifespan`
SRL.effects = SRL.effects %>%
  filter(local_lifespan == "ANNUAL")
SRL.effects$SRL.bt= (SRL.effects$SRL.final*74.26099) + 106.4204

x.value = c(-1,0,1,2,3,4)
(x.value*74.26099) + 106.4204

annual.SRL.plot = ggplot() +
  geom_point(data = imputed.NW.3.annual, aes(x = SRL.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Specific Root Length (m/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.497203,417.749988)+
  theme_classic()
annual.SRL.plot

ggsave("./Plots/Revision.2/annual.traits.NW.SRL.pdf", height = 3, width = 3)

#### Annual RTD ####

RTD.effects = annual.effects$`RTD.final:local_lifespan`
RTD.effects = RTD.effects %>%
  filter(local_lifespan == "ANNUAL")
RTD.effects$RTD.bt= (RTD.effects$RTD.final*0.1159369) + 0.2349392

x.value = c(-1,0,1,2,3)
(x.value*0.1159369) + 0.2349392

annual.RTD.plot = ggplot() +
  geom_point(data = imputed.NW.3.annual, aes(x = RTD.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = expression("Root Tissue Density (g/cm"^3*")"), y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.02325658,0.61999993)+
  theme_classic()
annual.RTD.plot

ggsave("./Plots/Revision.2/annual.traits.NW.RTD.pdf", height = 3, width = 3)

#### Annual RMF ####

RMF.effects = annual.effects$`RMF.final:local_lifespan`
RMF.effects = RMF.effects %>%
  filter(local_lifespan == "ANNUAL")
RMF.effects$RMF.bt= (RMF.effects$RMF.final*0.1058164) + 0.4007788

x.value = c(-3,-2,-1,0,1,2,3)
(x.value*0.1058164) + 0.4007788

annual.RMF.plot = ggplot() +
  geom_point(data = imputed.NW.3.annual, aes(x = RMF.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RMF.effects, aes(x = RMF.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RMF.effects, aes(x = RMF.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Root Mass Fraction (g/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.1211193,0.7328659)+
  theme_classic()
annual.RMF.plot

ggsave("./Plots/Revision.2/annual.traits.NW.RMF.pdf", height = 3, width = 3)

#### Annual DSI ####

DSI.effects = annual.effects$`mean.DSI:local_lifespan`
DSI.effects = DSI.effects %>%
  filter(local_lifespan == "ANNUAL")
DSI.effects$DSI.bt= (DSI.effects$mean.DSI*0.1950438) + -0.4604709

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1950438) + -0.4604709

annual.DSI.plot = ggplot() +
  geom_point(data = imputed.NW.3.annual, aes(x = DSI.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = DSI.effects, aes(x = DSI.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = DSI.effects, aes(x = DSI.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Drought Severity Index", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(-0.8722616,0.1242846)+
  theme_classic()
annual.DSI.plot

ggsave("./Plots/Revision.2/annual.traits.NW.DSI.pdf", height = 3, width = 3)

#### Annual MAP ####

MAP.effects = annual.effects$`mean.MAP:local_lifespan`
MAP.effects = MAP.effects %>%
  filter(local_lifespan == "ANNUAL")
MAP.effects$MAP.bt= (MAP.effects$mean.MAP*487.1451) + 691.6107

x.value = c(-1,0,1,2,3)
(x.value*487.1451) + 691.6107

annual.MAP.plot = ggplot() +
  geom_point(data = imputed.NW.3.annual, aes(x = MAP.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = MAP.effects, aes(x = MAP.bt, y = estimate__), color = "#979461", size = 1.5) +  
  geom_ribbon(data = MAP.effects, aes(x = MAP.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Mean Annual Precipitation (mm)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(132.7,2366.1)+
  theme_classic()
annual.MAP.plot

ggsave("./Plots/Revision.2/annual.traits.NW.MAP.pdf", height = 3, width = 3)


#### PERENNIAL Backtransform ####

perennial.traits.NW.model = readRDS("./Results/lifespan.cats.imputed.traits.no_woody.rds")

imputed.NW.3$leafN.bt =(imputed.NW.3$leafN.final*7.114868) + 21.16957
imputed.NW.3$height.bt = (imputed.NW.3$height.final*0.2353906) + 0.3476824
imputed.NW.3$rootN.bt = (imputed.NW.3$rootN.final*4.488773) + 9.855027
imputed.NW.3$SLA.bt = (imputed.NW.3$SLA.final*8.581554) + 19.93883
imputed.NW.3$Depth.bt = (imputed.NW.3$root.depth.final*0.5333377) + 0.5784653
imputed.NW.3$Diam.bt = (imputed.NW.3$rootDiam.final*0.1598629) + 0.3501538
imputed.NW.3$SRL.bt = (imputed.NW.3$SRL.final*74.26099) + 106.4204
imputed.NW.3$RTD.bt = (imputed.NW.3$RTD.final*0.1159369) + 0.2349392
imputed.NW.3$RMF.bt = (imputed.NW.3$RMF.final*0.1058164) + 0.4007788
imputed.NW.3$DSI.bt = (imputed.NW.3$mean.DSI*0.1950438) + -0.4604709
imputed.NW.3$MAP.bt = (imputed.NW.3$mean.MAP*487.1451) + 691.6107

imputed.NW.3.perennial = imputed.NW.3 %>%
  filter(local_lifespan == "PERENNIAL")

perennial.effects = conditional_effects(perennial.traits.NW.model)

#### Perennial LeafN ####

leafN.effects = perennial.effects$`leafN.final:local_lifespan`
leafN.effects = leafN.effects %>%
  filter(local_lifespan == "PERENNIAL")
leafN.effects$leafN.bt= (leafN.effects$leafN.final*7.114868) + 21.16957

x.value = c(-2,-1,0,1,2,3)
(x.value*7.114868) + 21.16957

perennial.leafN.plot = ggplot() +
  geom_point(data = imputed.NW.3.perennial, aes(x = leafN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Leaf N (mg/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.68951,44.19001)+
  theme_classic()
perennial.leafN.plot

ggsave("./Plots/Revision.2/perennial.traits.NW.leafN.pdf", height = 3, width = 3)

#### Perennial Height ####

height.effects = perennial.effects$`height.final:local_lifespan`
height.effects = height.effects %>%
  filter(local_lifespan == "PERENNIAL")
height.effects$height.bt= (height.effects$height.final*0.2353906) + 0.3476824

x.value = c(-1,0,1,2,3,4)
(x.value*0.2353906) + 0.3476824

perennial.height.plot = ggplot() +
  geom_point(data = imputed.NW.3.perennial, aes(x = height.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = height.effects, aes(x = height.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Height (m)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.009000064,1.322333169)+
  theme_classic()
perennial.height.plot

ggsave("./Plots/Revision.2/perennial.traits.NW.height.pdf", height = 3, width = 3)

#### Perennial RootN ####

rootN.effects = perennial.effects$`rootN.final:local_lifespan`
rootN.effects = rootN.effects %>%
  filter(local_lifespan == "PERENNIAL")
rootN.effects$rootN.bt= (rootN.effects$rootN.final*4.488773) + 9.855027

x.value = c(-2,-1,0,1,2,3)
(x.value*4.488773) + 9.855027

perennial.rootN.plot = ggplot() +
  geom_point(data = imputed.NW.3.perennial, aes(x = rootN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Root N (mg/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.3800000,25.6283268)+
  theme_classic()
perennial.rootN.plot

ggsave("./Plots/Revision.2/perennial.traits.NW.rootN.pdf", height = 3, width = 3)

#### Perennial SLA ####

SLA.effects = perennial.effects$`SLA.final:local_lifespan`
SLA.effects = SLA.effects %>%
  filter(local_lifespan == "PERENNIAL")
SLA.effects$SLA.bt= (SLA.effects$SLA.final*8.581554) + 19.93883

x.value = c(-2,-1,0,1,2,3)
(x.value*8.581554) + 19.93883

perennial.SLA.plot = ggplot() +
  geom_point(data = imputed.NW.3.perennial, aes(x = SLA.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = expression("Specific Leaf Area (m"^2*"/kg)"), y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.617596,46.320001)+
  theme_classic()
perennial.SLA.plot

ggsave("./Plots/Revision.2/perennial.traits.NW.SLA.pdf", height = 3, width = 3)

#### Perennial Depth ####

root.depth.effects = perennial.effects$`root.depth.final:local_lifespan`
root.depth.effects = root.depth.effects %>%
  filter(local_lifespan == "PERENNIAL")
root.depth.effects$Depth.bt= (root.depth.effects$root.depth.final*0.5333377) + 0.5784653

x.value = c(-1,0,1,2,3,4,5)
(x.value*0.5333377) + 0.5784653

perennial.root.depth.plot = ggplot() +
  geom_point(data = imputed.NW.3.perennial, aes(x = Depth.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = Depth.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = Depth.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Rooting Depth (m)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.0479999,3.709335)+
  theme_classic()
perennial.root.depth.plot

ggsave("./Plots/Revision.2/perennial.traits.NW.root.depth.pdf", height = 3, width = 3)

#### Perennial Diam ####

rootDiam.effects = perennial.effects$`rootDiam.final:local_lifespan`
rootDiam.effects = rootDiam.effects %>%
  filter(local_lifespan == "PERENNIAL")
rootDiam.effects$Diam.bt= (rootDiam.effects$rootDiam.final*0.1598629) + 0.3501538

x.value = c(-1,0,1,2,3,4,5)
(x.value*0.1598629) + 0.3501538

perennial.rootDiam.plot = ggplot() +
  geom_point(data = imputed.NW.3.perennial, aes(x = Diam.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = Diam.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = Diam.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Root Diameter (mm)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.07840002,1.12601235)+
  theme_classic()
perennial.rootDiam.plot

ggsave("./Plots/Revision.2/perennial.traits.NW.rootDiam.pdf", height = 3, width = 3)

#### Perennial SRL ####

SRL.effects = perennial.effects$`SRL.final:local_lifespan`
SRL.effects = SRL.effects %>%
  filter(local_lifespan == "PERENNIAL")
SRL.effects$SRL.bt= (SRL.effects$SRL.final*74.26099) + 106.4204

x.value = c(-1,0,1,2,3,4)
(x.value*74.26099) + 106.4204

perennial.SRL.plot = ggplot() +
  geom_point(data = imputed.NW.3.perennial, aes(x = SRL.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Specific Root Length (m/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.497203,417.749988)+
  theme_classic()
perennial.SRL.plot

ggsave("./Plots/Revision.2/perennial.traits.NW.SRL.pdf", height = 3, width = 3)

#### Perennial RTD ####

RTD.effects = perennial.effects$`RTD.final:local_lifespan`
RTD.effects = RTD.effects %>%
  filter(local_lifespan == "PERENNIAL")
RTD.effects$RTD.bt= (RTD.effects$RTD.final*0.1159369) + 0.2349392

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1159369) + 0.2349392

perennial.RTD.plot = ggplot() +
  geom_point(data = imputed.NW.3.perennial, aes(x = RTD.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = expression("Root Tissue Density (g/cm"^3*")"), y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.02325658,0.61999993)+
  theme_classic()
perennial.RTD.plot

ggsave("./Plots/Revision.2/perennial.traits.NW.RTD.pdf", height = 3, width = 3)

#### Perennial RMF ####

RMF.effects = perennial.effects$`RMF.final:local_lifespan`
RMF.effects = RMF.effects %>%
  filter(local_lifespan == "PERENNIAL")
RMF.effects$RMF.bt= (RMF.effects$RMF.final*0.1058164) + 0.4007788

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1058164) + 0.4007788

perennial.RMF.plot = ggplot() +
  geom_point(data = imputed.NW.3.perennial, aes(x = RMF.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RMF.effects, aes(x = RMF.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RMF.effects, aes(x = RMF.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Root Mass Fraction (g/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.1211193,0.7328659)+
  theme_classic()
perennial.RMF.plot

ggsave("./Plots/Revision.2/perennial.traits.NW.RMF.pdf", height = 3, width = 3)

#### Perennial DSI ####

DSI.effects = perennial.effects$`mean.DSI:local_lifespan`
DSI.effects = DSI.effects %>%
  filter(local_lifespan == "PERENNIAL")
DSI.effects$DSI.bt= (DSI.effects$mean.DSI*0.1950438) + -0.4604709

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1950438) + -0.4604709

perennial.DSI.plot = ggplot() +
  geom_point(data = imputed.NW.3.perennial, aes(x = DSI.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = DSI.effects, aes(x = DSI.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = DSI.effects, aes(x = DSI.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Drought Severity Index", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(-0.8722616,0.1242846)+
  theme_classic()
perennial.DSI.plot

ggsave("./Plots/Revision.2/perennial.traits.NW.DSI.pdf", height = 3, width = 3)

#### Perennial MAP ####

MAP.effects = perennial.effects$`mean.MAP:local_lifespan`
MAP.effects = MAP.effects %>%
  filter(local_lifespan == "PERENNIAL")
MAP.effects$MAP.bt= (MAP.effects$mean.MAP*487.1451) + 691.6107

x.value = c(-1,0,1,2,3)
(x.value*487.1451) + 691.6107

perennial.MAP.plot = ggplot() +
  geom_point(data = imputed.NW.3.perennial, aes(x = MAP.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = MAP.effects, aes(x = MAP.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = MAP.effects, aes(x = MAP.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Mean Annual Precipitation (mm)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(132.7,2366.1)+
  theme_classic()
perennial.MAP.plot

ggsave("./Plots/Revision.2/perennial.traits.NW.MAP.pdf", height = 3, width = 3)

#### FORB Backtransform ####

imputed.NW.functional.group = imputed.NW %>%
  select(cover.change,functional_group,leafN.final,height.final,rootN.final,SLA.final,root.depth.final,
         rootDiam.final,SRL.final,RTD.final,RMF.final,mean.DSI,mean.MAP,site_code,Taxon) %>%
  drop_na()

forb.traits.NW.model = readRDS("./Results/functional.group.cats.imputed.traits.no_woody.rds")

imputed.NW.functional.group$leafN.bt =(imputed.NW.functional.group$leafN.final*7.114868) + 21.16957
imputed.NW.functional.group$height.bt = (imputed.NW.functional.group$height.final*0.2353906) + 0.3476824
imputed.NW.functional.group$rootN.bt = (imputed.NW.functional.group$rootN.final*4.488773) + 9.855027
imputed.NW.functional.group$SLA.bt = (imputed.NW.functional.group$SLA.final*8.581554) + 19.93883
imputed.NW.functional.group$Depth.bt = (imputed.NW.functional.group$root.depth.final*0.5333377) + 0.5784653
imputed.NW.functional.group$Diam.bt = (imputed.NW.functional.group$rootDiam.final*0.1598629) + 0.3501538
imputed.NW.functional.group$SRL.bt = (imputed.NW.functional.group$SRL.final*74.26099) + 106.4204
imputed.NW.functional.group$RTD.bt = (imputed.NW.functional.group$RTD.final*0.1159369) + 0.2349392
imputed.NW.functional.group$RMF.bt = (imputed.NW.functional.group$RMF.final*0.1058164) + 0.4007788
imputed.NW.functional.group$DSI.bt = (imputed.NW.functional.group$mean.DSI*0.1950438) + -0.4604709
imputed.NW.functional.group$MAP.bt = (imputed.NW.functional.group$mean.MAP*487.1451) + 691.6107

imputed.NW.forb = imputed.NW.functional.group %>%
  filter(functional_group == "FORB")

forb.effects = conditional_effects(forb.traits.NW.model)

#### Forb leafN ####

leafN.effects = forb.effects$`leafN.final:functional_group`
leafN.effects = leafN.effects %>%
  filter(functional_group == "FORB")
leafN.effects$leafN.bt= (leafN.effects$leafN.final*7.114868) + 21.16957

x.value = c(-3,-2,-1,0,1,2,3)
(x.value*7.114868) + 21.16957

forb.leafN.plot = ggplot() +
  geom_point(data = imputed.NW.forb, aes(x = leafN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.bt, y = estimate__), color = "#F17236", size = 1.5) +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Leaf N (mg/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.68951,44.19001)+
  theme_classic()
forb.leafN.plot

ggsave("./Plots/Revision.2/forb.traits.NW.leafN.pdf", height = 3, width = 3)

#### Forb height ####

height.effects = forb.effects$`height.final:functional_group`
height.effects = height.effects %>%
  filter(functional_group == "FORB")
height.effects$height.bt= (height.effects$height.final*0.2353906) + 0.3476824

x.value = c(-1,0,1,2,3,4)
(x.value*0.2353906) + 0.3476824

forb.height.plot = ggplot() +
  geom_point(data = imputed.NW.forb, aes(x = height.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = height.effects, aes(x = height.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Height (m)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.009000064,1.322333169)+
  theme_classic()
forb.height.plot

ggsave("./Plots/Revision.2/forb.traits.NW.height.pdf", height = 3, width = 3)

#### Forb rootN ####

rootN.effects = forb.effects$`rootN.final:functional_group`
rootN.effects = rootN.effects %>%
  filter(functional_group == "FORB")
rootN.effects$rootN.bt= (rootN.effects$rootN.final*4.488773) + 9.855027

x.value = c(-2,-1,0,1,2,3)
(x.value*4.488773) + 9.855027

forb.rootN.plot = ggplot() +
  geom_point(data = imputed.NW.forb, aes(x = rootN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Root N (mg/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.3800000,25.6283268) +
  theme_classic()
forb.rootN.plot

ggsave("./Plots/Revision.2/forb.traits.NW.rootN.pdf", height = 3, width = 3)

#### Forb SLA ####

SLA.effects = forb.effects$`SLA.final:functional_group`
SLA.effects = SLA.effects %>%
  filter(functional_group == "FORB")
SLA.effects$SLA.bt= (SLA.effects$SLA.final*8.581554) + 19.93883

x.value = c(-2,-1,0,1,2,3)
(x.value*8.581554) + 19.93883

forb.SLA.plot = ggplot() +
  geom_point(data = imputed.NW.forb, aes(x = SLA.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = expression("Specific Leaf Area (m"^2*"/kg)"), y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.617596,46.320001)+
  theme_classic()
forb.SLA.plot

ggsave("./Plots/Revision.2/forb.traits.NW.SLA.pdf", height = 3, width = 3)

#### Forb Depth ####

root.depth.effects = forb.effects$`root.depth.final:functional_group`
root.depth.effects = root.depth.effects %>%
  filter(functional_group == "FORB")
root.depth.effects$Depth.bt= (root.depth.effects$root.depth.final*0.5333377) + 0.5784653

x.value = c(-1,0,1,2,3,4,5)
(x.value*0.5333377) + 0.5784653

forb.root.depth.plot = ggplot() +
  geom_point(data = imputed.NW.forb, aes(x = Depth.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = Depth.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = Depth.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Rooting Depth (m)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.0479999,3.709335)+
  theme_classic()
forb.root.depth.plot

ggsave("./Plots/Revision.2/forb.traits.NW.root.depth.pdf", height = 3, width = 3)

#### Forb Diam ####

rootDiam.effects = forb.effects$`rootDiam.final:functional_group`
rootDiam.effects = rootDiam.effects %>%
  filter(functional_group == "FORB")
rootDiam.effects$Diam.bt= (rootDiam.effects$rootDiam.final*0.1598629) + 0.3501538

x.value = c(-1,0,1,2,3,4,5)
(x.value*0.1598629) + 0.3501538

forb.rootDiam.plot = ggplot() +
  geom_point(data = imputed.NW.forb, aes(x = Diam.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = Diam.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = Diam.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Root Diameter (mm)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.07840002,1.12601235)+
  theme_classic()
forb.rootDiam.plot

ggsave("./Plots/Revision.2/forb.traits.NW.rootDiam.pdf", height = 3, width = 3)

#### Forb SRL ####

SRL.effects = forb.effects$SRL.final
SRL.effects = SRL.effects %>%
  filter(functional_group == "FORB")
SRL.effects$SRL.bt= (SRL.effects$SRL.final*74.26099) + 106.4204

x.value = c(-1,0,1,2,3,4)
(x.value*74.26099) + 106.4204

forb.SRL.plot = ggplot() +
  geom_point(data = imputed.NW.forb, aes(x = SRL.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Specific Root Length (m/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.497203,417.749988)+
  theme_classic()
forb.SRL.plot

ggsave("./Plots/Revision.2/forb.traits.NW.SRL.pdf", height = 3, width = 3)

#### Forb RTD ####

RTD.effects = forb.effects$RTD.final
RTD.effects = RTD.effects %>%
  filter(functional_group == "FORB")
RTD.effects$RTD.bt= (RTD.effects$RTD.final*0.1159369) + 0.2349392

x.value = c(-1,0,1,2,3)
(x.value*0.1159369) + 0.2349392

forb.RTD.plot = ggplot() +
  geom_point(data = imputed.NW.forb, aes(x = RTD.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = expression("Root Tissue Density (g/cm"^3*")"), y = "Cover Change (%)") + 
  ylim(-20.66667,21.24500)+
  #xlim(0.02325658,0.61999993)+
  theme_classic()
forb.RTD.plot

ggsave("./Plots/Revision.2/forb.traits.NW.RTD.pdf", height = 3, width = 3)

#### Forb RMF ####

RMF.effects = forb.effects$`RMF.final:functional_group`
RTD.effects = RTD.effects %>%
  filter(functional_group == "FORB")
RMF.effects$RMF.bt= (RMF.effects$RMF.final*0.1159929) + 0.3917262

attr(forb.imputed.NW.2$RMF.final, "scaled:scale")
attr(forb.imputed.NW.2$RMF.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1176732) + 0.3992223

forb.RMF.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = RMF.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RMF.effects, aes(x = RMF.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RMF.effects, aes(x = RMF.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Root Mass Fraction (g/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.1211193,0.7328659)+
  theme_classic()
forb.RMF.plot

ggsave("./Plots/forb.traits.NW.RMF.pdf", height = 3, width = 3)

#### Forb DSI ####

DSI.effects = forb.effects$mean.DSI
DSI.effects$DSI.bt= (DSI.effects$mean.DSI*0.1957044) + -0.4647529

attr(forb.imputed.NW.2$mean.DSI, "scaled:scale")
attr(forb.imputed.NW.2$mean.DSI, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1949549) + -0.4644806

forb.DSI.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = DSI.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = DSI.effects, aes(x = DSI.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = DSI.effects, aes(x = DSI.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Drought Severity Index", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(-0.8722616,0.1242846)+
  theme_classic()
forb.DSI.plot

ggsave("./Plots/forb.traits.NW.DSI.pdf", height = 3, width = 3)

#### Forb MAP ####

MAP.effects = forb.effects$mean.MAP
MAP.effects$MAP.bt= (MAP.effects$mean.MAP*433.6646) + 664.4671

attr(forb.imputed.NW.2$mean.MAP, "scaled:scale")
attr(forb.imputed.NW.2$mean.MAP, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*448.2603) + 657.7821

forb.MAP.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = MAP.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = MAP.effects, aes(x = MAP.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = MAP.effects, aes(x = MAP.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Mean Annual Precipitation (mm)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(132.7,2366.1)+
  theme_classic()
forb.MAP.plot

ggsave("./Plots/forb.traits.NW.MAP.pdf", height = 3, width = 3)


#### Model with lifespan groups: Relative BACI ####

priors <- c(prior(normal(0, 10), class = b))

relative.imputed.traits.NW.lifespan.model = brm(relative.cover.change ~ local_lifespan + 
                                                  leafN.final*local_lifespan + height.final*local_lifespan + 
                                                  rootN.final*local_lifespan + SLA.final*local_lifespan +
                                                  root.depth.final*local_lifespan + rootDiam.final*local_lifespan +
                                                  SRL.final*local_lifespan + RTD.final*local_lifespan + 
                                                  RMF.final*local_lifespan + mean.DSI*local_lifespan + 
                                                  mean.MAP*local_lifespan + (1|site_code) + (1|Taxon), 
                                                family = gaussian(),
                                                prior = priors,
                                                data = imputed.NW.3)

# checking model assumptions
pp_check(relative.imputed.traits.NW.lifespan.model)
pp_check(relative.imputed.traits.NW.lifespan.model, type = "error_scatter_avg")

summary(relative.imputed.traits.NW.lifespan.model)
bayes_R2(relative.imputed.traits.NW.lifespan.model)
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


