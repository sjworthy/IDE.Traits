## Script to evaluate linear mixed effects models of cover change and traits
# will run the models for:
# complete cases with and without woody species
# with and without RMF, with and without woody species
# imputed traits with and without woody species

# load libraries 

library(tidyverse)
library(lme4)
library(lmerTest)
library(MuMIn)
library(corrplot)
library(effects)
library(sjPlot)

# steps
# test for trait correlations
# fit model with only single traits
# include belowground trait interactions from the two different axes: RTD and SRL

#### final.data.CC ####
# model of trait data with complete cases for all species

final.data.CC = read.csv("./Formatted.Data/Revisions/final.data.CC.csv",row.names = 1)

# look at correlations between traits
cor.traits.cc = cor(final.data.CC[,c(4:12)],use = "pairwise") 
corrplot(cor.traits.cc, method="number",tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")
# highest correlation is -0.38 and 0.35

final.data.CC.2 = scale(final.data.CC[,c(4:12)])
final.data.CC.3 = cbind(final.data.CC[,c(1:3)],final.data.CC.2)

# 236 data points estimating 11 variables

CC.model = lmer(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                  SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + (1|site_code) + (1|Taxon), data = final.data.CC.3)
summary(CC.model)
# RTD is positive, significant
r.squaredGLMM(CC.model) # 0.06, 0.21

plot(CC.model)
output.plot = allEffects(CC.model)
plot(output.plot)
hist(resid(CC.model))

RTD.plot = plot_model(CC.model, type = "pred", terms = "RTD.groot.cahill.merge")
RTD.plot + geom_point(data = final.data.CC.3, aes(x = RTD.groot.cahill.merge, y = cover.change))

# Add interaction between RTD and SRL

CC.model.interact = lmer(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                  SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + SRL.groot.cahill.merge*RTD.groot.cahill.merge +
                    (1|site_code) + (1|Taxon), data = final.data.CC.3)
summary(CC.model.interact)
# RTD is significant
r.squaredGLMM(CC.model.interact) # 0.07, 0.21

plot(CC.model.interact)
output.plot = allEffects(CC.model.interact)
plot(output.plot)
hist(resid(CC.model.interact))

#### final.data.CC.NW ####
# model of trait data with complete cases only non-woody species

final.data.CC.NW = read.csv("./Formatted.Data/Revisions/final.data.CC.NW.csv",row.names = 1)

# look at correlations between traits
cor.traits.cc.nw = cor(final.data.CC.NW[,c(4:12)],use = "pairwise") 
corrplot(cor.traits.cc.nw, method="number",tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")
# highest correlation is -0.39 and 0.35

final.data.CC.NW.2 = scale(final.data.CC.NW[,c(4:12)])
final.data.CC.NW.3 = cbind(final.data.CC.NW[,c(1:3)],final.data.CC.NW.2)


# 224 data points estimating 11 variables

CC.NW.model = lmer(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                  SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + (1|Taxon) + (1|site_code), data = final.data.CC.NW.3)
summary(CC.NW.model)

plot(CC.NW.model)
output.plot = allEffects(CC.NW.model)
plot(output.plot)
hist(resid(CC.NW.model))

# Add interaction between RTD and SRL

CC.NW.model.interact = lmer(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                           SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + SRL.groot.cahill.merge*RTD.groot.cahill.merge +
                           (1|site_code) + (1|Taxon), data = final.data.CC.NW.3)
summary(CC.NW.model.interact)
# RTD is significant

#### final.data.CC.weiglet ####

# model of trait data with complete cases for all species, excluding RMF

final.data.CC.weigelt = read.csv("./Formatted.Data/Revisions/final.data.CC.weigelt.csv",row.names = 1)

# look at correlations between traits
final.data.cc.weigelt = cor(final.data.CC.weigelt[,c(4:11)],use = "pairwise") 
corrplot(final.data.cc.weigelt, method="number",tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")
# highest correlation is -0.36 and 0.35

final.data.cc.weigelt.2 = scale(final.data.CC.weigelt[,c(4:11)])
final.data.cc.weigelt.3 = cbind(final.data.CC.weigelt[,c(1:3)],final.data.cc.weigelt.2)


# 288 data points estimating 10 variables

CC.weigelt.model = lmer(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                  SRL.groot.cahill.merge + RTD.groot.cahill.merge + (1|site_code) + (1|Taxon), data = final.data.cc.weigelt.3)
summary(CC.weigelt.model)
# RTD is positive, significant

# Add interaction between RTD and SRL

CC.weigelt.model.interact = lmer(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                          SRL.groot.cahill.merge + RTD.groot.cahill.merge + SRL.groot.cahill.merge*RTD.groot.cahill.merge +
                            (1|site_code) + (1|Taxon), data = final.data.cc.weigelt.3)
summary(CC.weigelt.model.interact)
# RTD is positive, significant

#### final.data.CC.NW.weiglet ####
# model of trait data with complete cases only non-woody species

final.data.CC.NW.weigelt = read.csv("./Formatted.Data/Revisions/final.data.CC.NW.weigelt.csv",row.names = 1)

# look at correlations between traits
cor.traits.cc.nw.weigelt = cor(final.data.CC.NW.weigelt[,c(4:11)],use = "pairwise") 
corrplot(cor.traits.cc.nw.weigelt, method="number",tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")
# highest correlation is -0.35 and 0.34

final.data.CC.NW.weigelt.2 = scale(final.data.CC.NW.weigelt[,c(4:11)])
final.data.CC.NW.weigelt.3 = cbind(final.data.CC.NW.weigelt[,c(1:3)],final.data.CC.NW.weigelt.2)

# 249 data points estimating 11 variables

CC.NW.weigelt.model = lmer(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                     SRL.groot.cahill.merge + RTD.groot.cahill.merge+ (1|Taxon) + (1|site_code), data = final.data.CC.NW.weigelt.3)
summary(CC.NW.weigelt.model)
# nothing significant

# Add interaction between RTD and SRL

CC.NW.weigelt.model.interact = lmer(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                             SRL.groot.cahill.merge + RTD.groot.cahill.merge + SRL.groot.cahill.merge*RTD.groot.cahill.merge + 
                               (1|Taxon) + (1|site_code), data = final.data.CC.NW.weigelt.3)
summary(CC.NW.weigelt.model.interact)
# RTD is significant

#### imputed.traits.final ####
# model of trait data with complete cases for all species

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)

# look at correlations between traits
cor.traits.impute = cor(imputed.traits[,c(26:34)],use = "pairwise") 
corrplot(cor.traits.impute, method="number",tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")
# highest correlation is -0.31 and 0.35

imputed.traits.2 = scale(imputed.traits[,c(26:34)])
imputed.traits.3 = cbind(imputed.traits[,c(1:3)],imputed.traits.2)

# 907 data points estimating 11 variables

imputed.traits.model = lmer(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                  SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), data = imputed.traits.3)
summary(imputed.traits.model)
# nothing significant

plot(imputed.traits.model)
output.plot = allEffects(imputed.traits.model)
plot(output.plot)
hist(resid(imputed.traits.model))

# Add interaction between RTD and SRL

imputed.traits.model.interact = lmer(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                              SRL.final + RTD.final + RMF.final + SRL.final*RTD.final + (1|site_code) + (1|Taxon), data = imputed.traits.3)
summary(imputed.traits.model.interact)
# nothing is significant

#### imputed.traits.NW.final ####
# model of trait data with complete cases for all species

imputed.traits.NW = read.csv("./Formatted.Data/Revisions/imputed.traits.NW.final.csv", row.names = 1)

# look at correlations between traits
cor.traits.impute = cor(imputed.traits.NW[,c(26:34)],use = "pairwise") 
corrplot(cor.traits.impute, method="number",tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")
# highest correlation is -0.27 and 0.34

imputed.traits.NW.2 = scale(imputed.traits.NW[,c(26:34)])
imputed.traits.NW.3 = cbind(imputed.traits.NW[,c(1:3)],imputed.traits.NW.2)

# 794 data points estimating 11 variables

imputed.traits.model.nw = lmer(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                              SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), data = imputed.traits.NW.3)
summary(imputed.traits.model.nw)
# nothing significant
r.squaredGLMM(imputed.traits.model.nw) # 0.02, 0.13

# Add interaction between RTD and SRL

imputed.traits.model.nw.interact = lmer(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                       SRL.final + RTD.final + RMF.final + SRL.final*RTD.final + (1|site_code) + (1|Taxon), data = imputed.traits.NW.3)
summary(imputed.traits.model.nw.interact)
# nothing is significant

#### imputed.traits.final lifespan ####

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)

table(imputed.traits$local_lifespan)

annual.imputed = imputed.traits %>%
  filter(local_lifespan == "ANNUAL")
perennial.imputed = imputed.traits %>%
  filter(local_lifespan == "PERENNIAL")

annual.imputed.2 = scale(annual.imputed[,c(26:34)])
annual.imputed.3 = cbind(annual.imputed[,c(1:3)],annual.imputed.2)

annual.impute.model = lmer(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                              SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), data = annual.imputed.3)
summary(annual.impute.model)

perennial.imputed.2 = scale(perennial.imputed[,c(26:34)])
perennial.imputed.3 = cbind(perennial.imputed[,c(1:3)],perennial.imputed.2)

perennial.impute.model = lmer(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                             SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), data = perennial.imputed.3)
summary(perennial.impute.model)
# RTD is significant


#### imputed.traits.final functional group ####

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)

table(imputed.traits$functional_group)

forb.imputed = imputed.traits %>%
  filter(functional_group == "FORB")
grass.imputed = imputed.traits %>%
  filter(functional_group == "GRASS")

forb.imputed.2 = scale(forb.imputed[,c(26:34)])
forb.imputed.3 = cbind(forb.imputed[,c(1:3)],forb.imputed.2)

forb.impute.model = lmer(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                             SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), data = forb.imputed.3)
summary(forb.impute.model)

grass.imputed.2 = scale(grass.imputed[,c(26:34)])
grass.imputed.3 = cbind(grass.imputed[,c(1:3)],grass.imputed.2)

grass.impute.model = lmer(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), data = grass.imputed.3)
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

annual.forb.imputed.2 = scale(annual.forb.imputed[,c(26:34)])
annual.forb.imputed.3 = cbind(annual.forb.imputed[,c(1:3)],annual.forb.imputed.2)

annual.forb.impute.model = lmer(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                           SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), data = annual.forb.imputed.3)
summary(annual.forb.impute.model)

perennial.forb.imputed.2 = scale(perennial.forb.imputed[,c(26:34)])
perennial.forb.imputed.3 = cbind(perennial.forb.imputed[,c(1:3)],perennial.forb.imputed.2)

perennial.forb.impute.model = lmer(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                     SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), data = perennial.forb.imputed.3)
summary(perennial.forb.impute.model)

perennial.grass.imputed.2 = scale(perennial.grass.imputed[,c(26:34)])
perennial.grass.imputed.3 = cbind(perennial.grass.imputed[,c(1:3)],perennial.grass.imputed.2)

perennial.grass.impute.model = lmer(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                      SRL.final + RTD.final + RMF.final + (1|site_code) + (1|Taxon), data = perennial.grass.imputed.3)
summary(perennial.grass.impute.model)


