# Script to prep trait data for analyses

# load libraries
library(tidyverse)
library(corrplot)

# read in full data frame 
trait.data = read.csv("./Formatted.Data/Revisions/BACI.trait.meta.csv", row.names = 1) 
# 1096 individuals, 749 species

# subset data for analyses
trait.data.2 = trait.data %>%
  select(site_code,Taxon,cover.change,leafN.mg.g,height.m,rootN.mg.g,SLA_m2.kg,root.depth_m,rootDiam.mm,
         SRL.groot.cahill.merge,RTD.groot.cahill.merge,RMF.g.g,local_lifespan,local_lifeform,functional_group)

# subset out bryophytes (9), mosses (6), cactus (1), none had SLA and generally few trait data available

trait.data.3 = trait.data.2 %>%
  filter(!local_lifeform %in% c("BRYOPHYTE","MOSS","CACTUS")) 
# 1080 individuals, 739 species

# write.csv(trait.data.3, file = "./Formatted.Data/Revisions/final.data.for.impute.csv")

# filter out woody species and trees based on functional group and local_lifeform

trait.data.NW = trait.data.3 %>%
  filter(!functional_group == "WOODY") %>%
  filter(!local_lifeform %in% c("SHRUB","TREE"))
# 956 individuals, 648 species
# lost 124 individuals (11%) and 91 species (12%)

# write.csv(trait.data.NW, file = "./Formatted.Data/Revisions/final.data.NW.for.impute.csv")

#### Get complete cases of traits ####

traits.cc.nw = trait.data.NW[complete.cases(trait.data.NW), ]
# 224 individuals, 96 species

# write.csv(traits.cc.nw, file = "./Formatted.Data/Revisions/final.data.CC.NW.csv")

traits.cc = trait.data.3[complete.cases(trait.data.3), ]
# 236 individuals, 103 species

# write.csv(traits.cc.nw, file = "./Formatted.Data/Revisions/final.data.CC.csv")

# look at correlations between traits

cor.traits.cc.nw = cor(traits.cc.nw[,c(4:12)],use = "pairwise") 
corrplot(cor.traits.cc.nw, method="number",tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")
# highest correlation is -0.39 and 0.35

cor.traits.cc = cor(traits.cc[,c(4:12)],use = "pairwise") 
corrplot(cor.traits.cc, method="number",tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")
# highest correlation is -0.38 and 0.35

#### complete cases without RMF since not in Weigelt 2021 ####

traits.data.NW.2 = trait.data.NW[,c(1:11,13:15)]

traits.cc.nw.weigelt = traits.data.NW.2[complete.cases(traits.data.NW.2), ]
# 259 individuals, 114 species

# write.csv(traits.cc.nw, file = "./Formatted.Data/Revisions/final.data.CC.NW.weigelt.csv")

trait.data.3.2 = trait.data.3[,c(1:11,13:15)]

traits.cc.weigelt = trait.data.3.2[complete.cases(trait.data.3.2), ]
# 288 individuals, 129 species

# write.csv(traits.cc.nw, file = "./Formatted.Data/Revisions/final.data.CC.weigelt.csv")


#### data set split by lifespan ####

# only enough data for perennials to fit models, not annuals

table(traits.cc.nw$local_lifespan)
# 184 perennial, 31 annual
table(traits.cc$local_lifespan)
# 196 perennial, 31 annual
table(traits.cc.nw.weigelt$local_lifespan)
# 218 perennial, 32 annual
table(traits.cc.weigelt$local_lifespan)
# 247 perennial, 32 annual

#### data set split by functional group

table(traits.cc.nw$functional_group)
# 97 forb, 3 graminoid, 94 grass, 30 legume
table(traits.cc$functional_group)
# 97 forb, 3 graminoid, 94 grass, 30 legume, 12 woody
table(traits.cc.nw.weigelt$functional_group)
# 114 forb, 5 graminoid, 106 grass, 34 legume
table(traits.cc.weigelt$functional_group)
# 116 forb, 5 graminoid, 106 grass, 34 legume, 27 woody

#### checking for trait outliers from trait data without woody species ####
# not going to remove outliers since need complete cases of traits for individuals to be included in analyses using lmer

hist(trait.data.NW$leafN.mg.g)
boxplot(trait.data.NW$leafN.mg.g)
mean = mean(trait.data.NW$leafN.mg.g, na.rm = TRUE)
std = sd(trait.data.NW$leafN.mg.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.NW$leafN.mg.g[which(trait.data.NW$leafN.mg.g <Tmin | trait.data.NW$leafN.mg.g > Tmax)])
# removed leafN 54.71597
# percent removed 
table(is.na(trait.data.NW$leafN.mg.g))
958-368
2/590*100 #0.34%

hist(trait.data.NW$height.m)
boxplot(trait.data.NW$height.m)
mean = mean(trait.data.NW$height.m, na.rm = TRUE)
std = sd(trait.data.NW$height.m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.NW$height.m[which(trait.data.NW$height.m <Tmin | trait.data.NW$height.m > Tmax)])
# removed height 1.371600 - 3.80000
# percent removed 
table(is.na(trait.data.NW$height.m))
958-216
14/742*100 #1.89%

hist(trait.data.NW$rootN.mg.g)
boxplot(trait.data.NW$rootN.mg.g)
mean = mean(trait.data.NW$rootN.mg.g, na.rm = TRUE)
std = sd(trait.data.NW$rootN.mg.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.NW$rootN.mg.g[which(trait.data.NW$rootN.mg.g <Tmin | trait.data.NW$rootN.mg.g > Tmax)])
# remove rootN 31.17241-39.26274
# percent removed 
table(is.na(trait.data.NW$rootN.mg.g))
958-613
6/345*100 #1.74%

hist(trait.data.NW$SLA_m2.kg)
boxplot(trait.data.NW$SLA_m2.kg)
mean = mean(trait.data.NW$SLA_m2.kg, na.rm = TRUE)
std = sd(trait.data.NW$SLA_m2.kg, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.NW$SLA_m2.kg[which(trait.data.NW$SLA_m2.kg <Tmin | trait.data.NW$SLA_m2.kg > Tmax)])
# remove SLA 53.50 - 68.90
# percent removed 
table(is.na(trait.data.NW$SLA_m2.kg))
958-244
8/714*100 #1.12%

hist(trait.data.NW$root.depth_m)
boxplot(trait.data.NW$root.depth_m)
mean = mean(trait.data.NW$root.depth_m, na.rm = TRUE)
std = sd(trait.data.NW$root.depth_m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.NW$root.depth_m[which(trait.data.NW$root.depth_m <Tmin | trait.data.NW$root.depth_m > Tmax)])
# remove depth 2.426850 - 3.709333
# percent removed 
table(is.na(trait.data.NW$root.depth_m))
958-478
15/480*100 #3.13%

hist(trait.data.NW$RTD.groot.cahill.merge)
boxplot(trait.data.NW$RTD.groot.cahill.merge)
mean = mean(trait.data.NW$RTD.groot.cahill.merge, na.rm = TRUE)
std = sd(trait.data.NW$RTD.groot.cahill.merge, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.NW$RTD.groot.cahill.merge[which(trait.data.NW$RTD.groot.cahill.merge <Tmin | trait.data.NW$RTD.groot.cahill.merge > Tmax)])
# remove RTD 0.7000297
# percent removed 
table(is.na(trait.data.NW$RTD.groot.cahill.merge))
958-581
1/377*100 #0.27%

hist(trait.data.NW$SRL.groot.cahill.merge)
boxplot(trait.data.NW$SRL.groot.cahill.merge)
mean = mean(trait.data.NW$SRL.groot.cahill.merge, na.rm = TRUE)
std = sd(trait.data.NW$SRL.groot.cahill.merge, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.NW$SRL.groot.cahill.merge[which(trait.data.NW$SRL.groot.cahill.merge <Tmin | trait.data.NW$SRL.groot.cahill.merge > Tmax)])
# remove SRL 527.2000 - 760.8500
# percent removed 
table(is.na(trait.data.NW$SRL.groot.cahill.merge))
958-575
10/383*100 #2.61%

hist(trait.data.NW$rootDiam.mm)
boxplot(trait.data.NW$rootDiam.mm)
mean = mean(trait.data.NW$rootDiam.mm, na.rm = TRUE)
std = sd(trait.data.NW$rootDiam.mm, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.NW$rootDiam.mm[which(trait.data.NW$rootDiam.mm <Tmin | trait.data.NW$rootDiam.mm > Tmax)])
# remove diam 1.640000 - 4.830000
# percent removed 
table(is.na(trait.data.NW$rootDiam.mm))
958-556
4/402*100 #0.995%


#### read in new trait data without outliers and cover data ####

trait.data.new = read.csv("./Formatted.Data/BACI.traits.no.woody.SLA.no.outlier.csv", row.names = 1)
cover.data = read.csv("./Formatted.Data/BACI.data.final.csv", row.names = 1)

# merge trait data and cover data

all.data = merge(cover.data, trait.data.new, by="Taxon") 
# 420 species from 642 data points from 65 sites

#write.csv(all.data, file = "./Formatted.Data/BACI.all.data.csv")


#### data set split by lifespan ####
# need to read out data and fix the lifespan to particular sites since 
# some species have different lifespan at different sites

# merge all.data and metadata
all.data.ls = read.csv("./Formatted.Data/BACI.all.data.metadata.csv", row.names = 1)
table(all.data.ls$local_lifespan)
# 497 perennial, 126 annual

annual.data = subset(all.data.ls, all.data.ls$local_lifespan == "ANNUAL")
# 126 data points from 81 species

perennial.data = subset(all.data.ls, all.data.ls$local_lifespan == "PERENNIAL") 
# 497 data points of 324 species

table(all.data.ls$functional_group)
# 327 forb, 28 graminoid, 228 grass, 59 legume

grass = subset(all.data.ls, all.data.ls$functional_group == "GRASS") 
# 228 data points of 143 species

table(grass$local_lifespan)
# 49 annuals, 174 perennials

forb = subset(all.data.ls, all.data.ls$functional_group == "FORB")
# 327 data points of 221 species

table(forb$local_lifespan)
# 67 annuals, 249 perennials

# grass.annuals
grass.annual = subset(grass, grass$local_lifespan == "ANNUAL")
# 49 data points of 25 species

# grass.perennial
grass.perennial = subset(grass, grass$local_lifespan == "PERENNIAL")
# 174 data points of 114 species

# forb.annual
forb.annual = subset(forb, forb$local_lifespan == "ANNUAL")
# 67 data points of 50 species

# forb.perennial
forb.perennial = subset(forb, forb$local_lifespan == "PERENNIAL")
# 249 data points of 163 species

#### write out the files ####

#write.csv(all.data, file = "./Formatted.Data/BACI.all.data.csv")
#write.csv(annual.data, file = "./Formatted.Data/BACI.annual.data.csv")
#write.csv(forb, file = "./Formatted.Data/BACI.forb.csv")
#write.csv(grass, file = "./Formatted.Data/BACI.grass.csv")
#write.csv(perennial.data, file = "./Formatted.Data/BACI.perennial.data.csv")
#write.csv(forb.annual, file = "./Formatted.Data/BACI.forb.annual.csv")
#write.csv(grass.annual, file = "./Formatted.Data/BACI.grass.annual.csv")
#write.csv(forb.perennial, file = "./Formatted.Data/BACI.forb.perennial.csv")
#write.csv(grass.perennial, file = "./Formatted.Data/BACI.grass.perennial.csv")


