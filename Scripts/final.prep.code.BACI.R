# Script to prep trait data for analyses

# load libraries
library(tidyverse)
library(corrplot)

# read in full data frame 
trait.data = read.csv("./Raw.Data/trait.species.BACI.final.csv")

# subset traits to only those of interest
trait.data.2 = trait.data[,c(1:9,14:16)] # 1239 species 

# subset out woody species and trees from local_lifeform

trait.data.3 = subset(trait.data.2, !trait.data.2$functional_group == "WOODY") # 1117 species
trait.data.3 = subset(trait.data.3, !trait.data.3$local_lifeform == "TREE") # 1115 species

# BACI analysis:
# 244 new species, 29 removed because TREE or WOODY, 1 FERN removed, 4 CACTUS removed
# these were removed prior to adding species and traits to the end of trait.data
# 210 total new species added (if traits available)

# subset for species with SLA
trait.data.4 = subset(trait.data.3, trait.data.3$SLA_m2.kg > 0 ) # 692 taxon

#write.csv(trait.data.4, file = "./Formatted.Data/BACI.traits.no.woody.SLA.csv")

#### checking for outliers ####
# all outliers removed manually from traits.no.woody.SLA.csv and 
# made new file BACI.traits.no.woody.SLA.no.outlier.csv

hist(trait.data.4$leafN.mg.g)
boxplot(trait.data.4$leafN.mg.g)
mean = mean(trait.data.4$leafN.mg.g, na.rm = TRUE)
std = sd(trait.data.4$leafN.mg.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$leafN.mg.g[which(trait.data.4$leafN.mg.g <Tmin | trait.data.4$leafN.mg.g > Tmax)])
# removed leafN 54.82000 58.30000 67.58333
# percent removed 
table(is.na(trait.data.4$leafN.mg.g))
692-207
3/485*100 #0.62%

hist(trait.data.4$height.m)
boxplot(trait.data.4$height.m)
mean = mean(trait.data.4$height.m, na.rm = TRUE)
std = sd(trait.data.4$height.m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$height.m[which(trait.data.4$height.m <Tmin | trait.data.4$height.m > Tmax)])
# removed height 2.600709 - 10
# percent removed 
table(is.na(trait.data.4$height.m))
692-62
7/630*100 #1.11%

hist(trait.data.4$rootN.mg.g)
boxplot(trait.data.4$rootN.mg.g)
mean = mean(trait.data.4$rootN.mg.g, na.rm = TRUE)
std = sd(trait.data.4$rootN.mg.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$rootN.mg.g[which(trait.data.4$rootN.mg.g <Tmin | trait.data.4$rootN.mg.g > Tmax)])
# remove rootN 34.39000 34.84397 39.26274
# percent removed 
table(is.na(trait.data.4$rootN.mg.g))
692-490
3/202*100 #1.49%

hist(trait.data.4$SLA_m2.kg)
boxplot(trait.data.4$SLA_m2.kg)
mean = mean(trait.data.4$SLA_m2.kg, na.rm = TRUE)
std = sd(trait.data.4$SLA_m2.kg, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$SLA_m2.kg[which(trait.data.4$SLA_m2.kg <Tmin | trait.data.4$SLA_m2.kg > Tmax)])
# remove SLA 52.45601 - 98.20000
# percent removed 
table(is.na(trait.data.4$SLA_m2.kg))
13/692*100 #1.88%

hist(trait.data.4$root.depth_m)
boxplot(trait.data.4$root.depth_m)
mean = mean(trait.data.4$root.depth_m, na.rm = TRUE)
std = sd(trait.data.4$root.depth_m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$root.depth_m[which(trait.data.4$root.depth_m <Tmin | trait.data.4$root.depth_m > Tmax)])
# remove depth 2.500000 - 3.828571
# percent removed 
table(is.na(trait.data.4$root.depth_m))
692-353
9/339*100 #2.65%

hist(trait.data.4$RTD.g.cm3)
boxplot(trait.data.4$RTD.g.cm3)
mean = mean(trait.data.4$RTD.g.cm3, na.rm = TRUE)
std = sd(trait.data.4$RTD.g.cm3, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$RTD.g.cm3[which(trait.data.4$RTD.g.cm3 <Tmin | trait.data.4$RTD.g.cm3 > Tmax)])
# remove RTD 0.7839500 - 0.9167417
# percent removed 
table(is.na(trait.data.4$RTD.g.cm3))
692-458
3/234*100 #1.28%

hist(trait.data.4$SRL_m.g)
boxplot(trait.data.4$SRL_m.g)
mean = mean(trait.data.4$SRL_m.g, na.rm = TRUE)
std = sd(trait.data.4$SRL_m.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$SRL_m.g[which(trait.data.4$SRL_m.g <Tmin | trait.data.4$SRL_m.g > Tmax)])
# remove SRL 652.6000 - 929.9185
# percent removed 
table(is.na(trait.data.4$SRL_m.g))
692-380
7/312*100 #2.24%

hist(trait.data.4$rootDiam.mm)
boxplot(trait.data.4$rootDiam.mm)
mean = mean(trait.data.4$rootDiam.mm, na.rm = TRUE)
std = sd(trait.data.4$rootDiam.mm, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$rootDiam.mm[which(trait.data.4$rootDiam.mm <Tmin | trait.data.4$rootDiam.mm > Tmax)])
# remove diam 1.242171 - 2.015308
# percent removed 
table(is.na(trait.data.4$rootDiam.mm))
692-387
9/305*100 #2.95%

#### read in new trait data without outliers and cover data ####

trait.data.new = read.csv("./Formatted.Data/BACI.traits.no.woody.SLA.no.outlier.csv", row.names = 1)
cover.data = read.csv("./Formatted.Data/BACI.data.final.csv", row.names = 1)

# merge trait data and cover data

all.data = merge(cover.data, trait.data.new, by="Taxon") 
# 420 species from 642 data points from 65 sites

#write.csv(all.data, file = "./Formatted.Data/BACI.all.data.csv")

# deteremine number of woody species removed

woody.data = inner_join(cover.data, trait.data.2)
table(woody.data$functional_group) # 126 out of 1003
woody.data.2 = subset(woody.data, !woody.data$functional_group == "WOODY") #877
table(woody.data.2$local_lifeform) # 2 trees removed

# 127 out of 1003 populations, 13% of populations removed for being woody

woody.data.2 = subset(woody.data, woody.data$SLA_m2.kg > 0 ) # 745 taxon
table(woody.data.2$functional_group) # 102 out of 1003
woody.data.3 = subset(woody.data.2, !woody.data.2$functional_group == "WOODY") #643
table(woody.data.3$local_lifeform) # 1 trees removed

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


