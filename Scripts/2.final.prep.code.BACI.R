# Script to prep trait data for analyses
# The first of this script is run to generate dataframes used to impute traits.
# After imputation, this code is used to identify outlier trait values.

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

# write.csv(traits.cc, file = "./Formatted.Data/Revisions/final.data.CC.csv")

# look at correlations between traits

cor.traits.cc.nw = cor(traits.cc.nw[,c(4:12)],use = "pairwise") 
corrplot(cor.traits.cc.nw, method="number",tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")
# highest correlation is -0.39 and 0.35

cor.traits.cc = cor(traits.cc[,c(4:12)],use = "pairwise") 
corrplot(cor.traits.cc, method="number",tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")
# highest correlation is -0.38 and 0.35

#### Checking and removing outliers from complete cases data with woody species ####

final.data.CC = read.csv("./Formatted.Data/Revisions/final.data.CC.csv",row.names = 1)

hist(final.data.CC$leafN.mg.g)
boxplot(final.data.CC$leafN.mg.g)
mean = mean(final.data.CC$leafN.mg.g, na.rm = TRUE)
std = sd(final.data.CC$leafN.mg.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$leafN.mg.g[which(final.data.CC$leafN.mg.g <Tmin | final.data.CC$leafN.mg.g > Tmax)])
# removed leafN 54.71597
# percent removed 
2/236*100 #0.85%

hist(final.data.CC$height.m)
boxplot(final.data.CC$height.m)
mean = mean(final.data.CC$height.m, na.rm = TRUE)
std = sd(final.data.CC$height.m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$height.m[which(final.data.CC$height.m <Tmin | final.data.CC$height.m > Tmax)])
# removed height 12.84847 27.10186 32.57366
# percent removed 
3/236*100 #1.27%

hist(final.data.CC$rootN.mg.g)
boxplot(final.data.CC$rootN.mg.g)
mean = mean(final.data.CC$rootN.mg.g, na.rm = TRUE)
std = sd(final.data.CC$rootN.mg.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$rootN.mg.g[which(final.data.CC$rootN.mg.g <Tmin | final.data.CC$rootN.mg.g > Tmax)])
# remove rootN 31.17241 31.34076 33.34667 33.34667 33.34667
# percent removed 
5/236*100 #2.11%

hist(final.data.CC$SLA_m2.kg)
boxplot(final.data.CC$SLA_m2.kg)
mean = mean(final.data.CC$SLA_m2.kg, na.rm = TRUE)
std = sd(final.data.CC$SLA_m2.kg, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$SLA_m2.kg[which(final.data.CC$SLA_m2.kg <Tmin | final.data.CC$SLA_m2.kg > Tmax)])
# remove SLA: none
# percent removed 
0/236*100 #0%

hist(final.data.CC$root.depth_m)
boxplot(final.data.CC$root.depth_m)
mean = mean(final.data.CC$root.depth_m, na.rm = TRUE)
std = sd(final.data.CC$root.depth_m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$root.depth_m[which(final.data.CC$root.depth_m <Tmin | final.data.CC$root.depth_m > Tmax)])
# remove depth 2.426850 2.426850 2.598333 2.598333 2.682500 2.682500 2.915000
# percent removed 
7/236*100 #2.97%

hist(final.data.CC$RTD.groot.cahill.merge)
boxplot(final.data.CC$RTD.groot.cahill.merge)
mean = mean(final.data.CC$RTD.groot.cahill.merge, na.rm = TRUE)
std = sd(final.data.CC$RTD.groot.cahill.merge, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$RTD.groot.cahill.merge[which(final.data.CC$RTD.groot.cahill.merge <Tmin | final.data.CC$RTD.groot.cahill.merge > Tmax)])
# remove RTD 1.19455 1.19455 1.19455 1.19455
# percent removed 
4/236*100 #1.69%

hist(final.data.CC$SRL.groot.cahill.merge)
boxplot(final.data.CC$SRL.groot.cahill.merge)
mean = mean(final.data.CC$SRL.groot.cahill.merge, na.rm = TRUE)
std = sd(final.data.CC$SRL.groot.cahill.merge, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$SRL.groot.cahill.merge[which(final.data.CC$SRL.groot.cahill.merge <Tmin | final.data.CC$SRL.groot.cahill.merge > Tmax)])
# remove SRL 471.2364 471.2364 471.2364 601.5858 627.5050
# percent removed 
5/236*100 #2.12%

hist(final.data.CC$rootDiam.mm)
boxplot(final.data.CC$rootDiam.mm)
mean = mean(final.data.CC$rootDiam.mm, na.rm = TRUE)
std = sd(final.data.CC$rootDiam.mm, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$rootDiam.mm[which(final.data.CC$rootDiam.mm <Tmin | final.data.CC$rootDiam.mm > Tmax)])
# remove diam 1.3965
# percent removed 
1/236*100 #0.42%

hist(final.data.CC$RMF.g.g)
boxplot(final.data.CC$RMF.g.g)
mean = mean(final.data.CC$RMF.g.g, na.rm = TRUE)
std = sd(final.data.CC$RMF.g.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$RMF.g.g[which(final.data.CC$RMF.g.g <Tmin | final.data.CC$RMF.g.g > Tmax)])
# remove RMF: none
# percent removed 
0/236*100 #0%

# 25 total removed because some removed with multiple traits
final.data.CC.2 = final.data.CC %>%
  filter(round(leafN.mg.g, 5) != 54.71597) %>% # 2
  filter(round(height.m, 5) < 12.84847) %>% # 3
  filter(round(rootN.mg.g, 5) < 31.17241) %>% # 5
  filter(round(root.depth_m, 5) < 2.426850) %>% # 7
  filter(round(RTD.groot.cahill.merge, 5) < 1.19455) %>% # 4
  filter(round(SRL.groot.cahill.merge, 5) < 471.2364) %>% # 5
  filter(round(rootDiam.mm, 5) != 1.3965) # 1

# write.csv(final.data.CC.2, file = "./Formatted.Data/Revisions/Final.Data/final.data.CC.outliersRM.csv")

# split by lifespan

table(final.data.CC.2$local_lifespan)
# annuals (22), biennials (7), Both (1), perennial (180), three (1)

final.data.CC.perennial = final.data.CC.2 %>%
  filter(local_lifespan == "PERENNIAL")

# write.csv(final.data.CC.perennial, file = "./Formatted.Data/Revisions/Final.Data/final.data.CC.perennials.outliersRM.csv")

# split by functional group

table(final.data.CC.2$functional_group)
# forb (93), graminoid (3), grass (88), legume (22), woody (5)

final.data.forb = final.data.CC.2 %>%
  filter(functional_group == "FORB")

final.data.grass = final.data.CC.2 %>%
  filter(functional_group == "GRASS")

# write.csv(final.data.forb, file = "./Formatted.Data/Revisions/Final.Data/final.data.CC.forbs.outliersRM.csv")
# write.csv(final.data.grass, file = "./Formatted.Data/Revisions/Final.Data/final.data.CC.grass.outliersRM.csv")

# split by local lifespan x functional group
table(final.data.CC.2$local_lifespan,final.data.CC.2$functional_group)
# annual forb (10), perennial forb (78), perennial grass (77)

#### Checking and removing outliers from complete cases data without woody species ####

final.data.CC.NW = read.csv("./Formatted.Data/Revisions/final.data.CC.NW.csv", row.names = 1)

hist(final.data.CC.NW$leafN.mg.g)
boxplot(final.data.CC.NW$leafN.mg.g)
mean = mean(final.data.CC.NW$leafN.mg.g, na.rm = TRUE)
std = sd(final.data.CC.NW$leafN.mg.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC.NW$leafN.mg.g[which(final.data.CC.NW$leafN.mg.g <Tmin | final.data.CC.NW$leafN.mg.g > Tmax)])
# removed leafN 54.71597
# percent removed 
2/224*100 #0.89%

hist(final.data.CC.NW$height.m)
boxplot(final.data.CC.NW$height.m)
mean = mean(final.data.CC.NW$height.m, na.rm = TRUE)
std = sd(final.data.CC.NW$height.m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC.NW$height.m[which(final.data.CC.NW$height.m <Tmin | final.data.CC.NW$height.m > Tmax)])
# removed height 1.800048 2.847733
# percent removed 
2/224*100 #0.89%

hist(final.data.CC.NW$rootN.mg.g)
boxplot(final.data.CC.NW$rootN.mg.g)
mean = mean(final.data.CC.NW$rootN.mg.g, na.rm = TRUE)
std = sd(final.data.CC.NW$rootN.mg.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC.NW$rootN.mg.g[which(final.data.CC.NW$rootN.mg.g <Tmin | final.data.CC.NW$rootN.mg.g > Tmax)])
# remove rootN 31.17241 31.34076 33.34667 33.34667 33.34667
# percent removed 
5/224*100 #2.23%

hist(final.data.CC.NW$SLA_m2.kg)
boxplot(final.data.CC.NW$SLA_m2.kg)
mean = mean(final.data.CC.NW$SLA_m2.kg, na.rm = TRUE)
std = sd(final.data.CC.NW$SLA_m2.kg, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC.NW$SLA_m2.kg[which(final.data.CC.NW$SLA_m2.kg <Tmin | final.data.CC.NW$SLA_m2.kg > Tmax)])
# remove SLA: none
# percent removed 
0/224*100 #0%

hist(final.data.CC.NW$root.depth_m)
boxplot(final.data.CC.NW$root.depth_m)
mean = mean(final.data.CC.NW$root.depth_m, na.rm = TRUE)
std = sd(final.data.CC.NW$root.depth_m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC.NW$root.depth_m[which(final.data.CC.NW$root.depth_m <Tmin | final.data.CC.NW$root.depth_m > Tmax)])
# remove depth 2.426850 2.426850 2.598333 2.598333 2.682500 2.682500 2.915000
# percent removed 
7/224*100 #3.13%

hist(final.data.CC.NW$RTD.groot.cahill.merge)
boxplot(final.data.CC.NW$RTD.groot.cahill.merge)
mean = mean(final.data.CC.NW$RTD.groot.cahill.merge, na.rm = TRUE)
std = sd(final.data.CC.NW$RTD.groot.cahill.merge, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC.NW$RTD.groot.cahill.merge[which(final.data.CC.NW$RTD.groot.cahill.merge <Tmin | final.data.CC.NW$RTD.groot.cahill.merge > Tmax)])
# remove RTD: none
# percent removed 
0/224*100 #0.0%

hist(final.data.CC.NW$SRL.groot.cahill.merge)
boxplot(final.data.CC.NW$SRL.groot.cahill.merge)
mean = mean(final.data.CC.NW$SRL.groot.cahill.merge, na.rm = TRUE)
std = sd(final.data.CC.NW$SRL.groot.cahill.merge, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC.NW$SRL.groot.cahill.merge[which(final.data.CC.NW$SRL.groot.cahill.merge <Tmin | final.data.CC.NW$SRL.groot.cahill.merge > Tmax)])
# remove SRL 471.2364 471.2364 471.2364 601.5858 627.5050
# percent removed 
5/224*100 #2.23%

hist(final.data.CC.NW$rootDiam.mm)
boxplot(final.data.CC.NW$rootDiam.mm)
mean = mean(final.data.CC.NW$rootDiam.mm, na.rm = TRUE)
std = sd(final.data.CC.NW$rootDiam.mm, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC.NW$rootDiam.mm[which(final.data.CC.NW$rootDiam.mm <Tmin | final.data.CC.NW$rootDiam.mm > Tmax)])
# remove diam: none
# percent removed 
0/224*100 #0%

hist(final.data.CC.NW$RMF.g.g)
boxplot(final.data.CC.NW$RMF.g.g)
mean = mean(final.data.CC.NW$RMF.g.g, na.rm = TRUE)
std = sd(final.data.CC.NW$RMF.g.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC.NW$RMF.g.g[which(final.data.CC.NW$RMF.g.g <Tmin | final.data.CC.NW$RMF.g.g > Tmax)])
# remove RMF: none
# percent removed 
0/402*100 #0%

# 20 total removed because some removed with multiple traits
final.data.CC.NW.2 = final.data.CC.NW %>%
  filter(round(leafN.mg.g, 5) != 54.71597) %>% # 2
  filter(round(height.m, 5) < 1.800048) %>% # 3
  filter(round(rootN.mg.g, 5) < 31.17241) %>% # 5
  filter(round(root.depth_m, 5) < 2.426850) %>% # 7
  filter(round(SRL.groot.cahill.merge, 5) < 471.2364) # 5

#write.csv(final.data.CC.NW.2, file = "./Formatted.Data/Revisions/Final.Data/final.data.CC.NW.outliersRM.csv")

# split by lifespan

table(final.data.CC.NW.2$local_lifespan)
# annuals (21), biennials (7), Both (1), perennial (174), three (1)

final.data.CC.NW.perennial = final.data.CC.NW.2 %>%
  filter(local_lifespan == "PERENNIAL")

# write.csv(final.data.CC.NW.perennial, file = "./Formatted.Data/Revisions/Final.Data/final.data.CC.NW.perennials.outliersRM.csv")

# split by functional group

table(final.data.CC.NW.2$functional_group)
# forb (92), graminoid (3), grass (87), legume (22)

final.data.NW.forb = final.data.CC.NW.2 %>%
  filter(functional_group == "FORB")

final.data.NW.grass = final.data.CC.NW.2 %>%
  filter(functional_group == "GRASS")

#write.csv(final.data.NW.forb, file = "./Formatted.Data/Revisions/Final.Data/final.data.CC.NW.forbs.outliersRM.csv")
#write.csv(final.data.NW.grass, file = "./Formatted.Data/Revisions/Final.Data/final.data.CC.NW.grass.outliersRM.csv")

# split by local lifespan x functional group
table(final.data.CC.NW.2$local_lifespan,final.data.CC.NW.2$functional_group)
# annual forb (9), perennial forb (78), perennial grass (76)

#### Checking and removing outliers from imputed data with woody species ####

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.final.csv",row.names = 1)

hist(imputed.traits$leafN.final)
boxplot(imputed.traits$leafN.final)
mean = mean(imputed.traits$leafN.final, na.rm = TRUE)
std = sd(imputed.traits$leafN.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$leafN.final[which(imputed.traits$leafN.final <Tmin | imputed.traits$leafN.final > Tmax)])
# removed leafN 46.00560 54.71597 54.71597 81.93849
# percent removed 
4/907*100 #0.44%

hist(imputed.traits$height.final)
boxplot(imputed.traits$height.final)
mean = mean(imputed.traits$height.final, na.rm = TRUE)
std = sd(imputed.traits$height.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$height.final[which(imputed.traits$height.final <Tmin | imputed.traits$height.final > Tmax)])
# removed height 7.787850  9.457333 10.000000 10.000000 12.848467 19.050000 20.170667 20.501000 23.394667 27.101857 32.573656
# percent removed 
11/907*100 #1.21%

hist(imputed.traits$rootN.final)
boxplot(imputed.traits$rootN.final)
mean = mean(imputed.traits$rootN.final, na.rm = TRUE)
std = sd(imputed.traits$rootN.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$rootN.final[which(imputed.traits$rootN.final <Tmin | imputed.traits$rootN.final > Tmax)])
# remove rootN 25.55616 - 39.26274
# percent removed 
18/907*100 #1.98%

hist(imputed.traits$SLA.final)
boxplot(imputed.traits$SLA.final)
mean = mean(imputed.traits$SLA.final, na.rm = TRUE)
std = sd(imputed.traits$SLA.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$SLA.final[which(imputed.traits$SLA.final <Tmin | imputed.traits$SLA.final > Tmax)])
# remove SLA: 49.05439 53.50000 53.50000 57.70000 61.28000 63.13000 68.90000 68.90000
# percent removed 
8/907*100 #0.88%

hist(imputed.traits$root.depth.final)
boxplot(imputed.traits$root.depth.final)
mean = mean(imputed.traits$root.depth.final, na.rm = TRUE)
std = sd(imputed.traits$root.depth.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$root.depth.final[which(imputed.traits$root.depth.final <Tmin | imputed.traits$root.depth.final > Tmax)])
# remove depth 2.915000 - 6.451600
# percent removed 
15/907*100 #1.65%

hist(imputed.traits$RTD.final)
boxplot(imputed.traits$RTD.final)
mean = mean(imputed.traits$RTD.final, na.rm = TRUE)
std = sd(imputed.traits$RTD.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$RTD.final[which(imputed.traits$RTD.final <Tmin | imputed.traits$RTD.final > Tmax)])
# remove RTD 0.7460162 - 1.4757011
# percent removed 
11/907*100 #1.21%

hist(imputed.traits$SRL.final)
boxplot(imputed.traits$SRL.final)
mean = mean(imputed.traits$SRL.final, na.rm = TRUE)
std = sd(imputed.traits$SRL.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$SRL.final[which(imputed.traits$SRL.final <Tmin | imputed.traits$SRL.final > Tmax)])
# remove SRL 437.3523-808.0298
# percent removed 
21/907*100 #2.31%

hist(imputed.traits$rootDiam.final)
boxplot(imputed.traits$rootDiam.final)
mean = mean(imputed.traits$rootDiam.final, na.rm = TRUE)
std = sd(imputed.traits$rootDiam.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$rootDiam.final[which(imputed.traits$rootDiam.final <Tmin | imputed.traits$rootDiam.final > Tmax)])
# remove diam 1.200000 1.396500 1.640000 1.850000 1.963333 4.830000
# percent removed 
6/907*100 #0.66%

hist(imputed.traits$RMF.final)
boxplot(imputed.traits$RMF.final)
mean = mean(imputed.traits$RMF.final, na.rm = TRUE)
std = sd(imputed.traits$RMF.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits$RMF.final[which(imputed.traits$RMF.final <Tmin | imputed.traits$RMF.final > Tmax)])
# remove RMF: 0.7670084 0.7670084 0.7670084 0.7670084 0.7768355 0.7995028 0.7995028 0.8274671 0.8274671 0.8630137
# percent removed 
10/907*100 #1.10%

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

write.csv(imputed.traits.2, file = "./Formatted.Data/Revisions/Final.Data/imputed.traits.outliersRM.csv")

# split by lifespan

table(imputed.traits.2$local_lifespan)
# annuals (206), biennials (11), Both (9), perennial (584), three (2)

imputed.traits.perennial = imputed.traits.2 %>%
  filter(local_lifespan == "PERENNIAL")
imputed.traits.annual = imputed.traits.2 %>%
  filter(local_lifespan == "ANNUAL")

#write.csv(imputed.traits.perennial, file = "./Formatted.Data/Revisions/Final.Data/imputed.traits.perennials.outliersRM.csv")
#write.csv(imputed.traits.annual, file = "./Formatted.Data/Revisions/Final.Data/imputed.traits.annuals.outliersRM.csv")

# split by functional group

table(imputed.traits.2$functional_group)
# forb (398), graminoid (39), grass (250), legume (49), woody (76)

imputed.traits.forb = imputed.traits.2 %>%
  filter(functional_group == "FORB")

imputed.traits.grass = imputed.traits.2 %>%
  filter(functional_group == "GRASS")

#write.csv(imputed.traits.forb, file = "./Formatted.Data/Revisions/Final.Data/imputed.traits.forbs.outliersRM.csv")
#write.csv(imputed.traits.grass, file = "./Formatted.Data/Revisions/Final.Data/imputed.traits.grass.outliersRM.csv")

# split by lifespan and functional group

table(imputed.traits.2$local_lifespan,imputed.traits.2$functional_group)
# annual forb (140), perennial forb (242), perennial grass (196)

imputed.traits.annual.forb = imputed.traits.2 %>%
  filter(local_lifespan == "ANNUAL") %>%
  filter(functional_group == "FORB")

imputed.traits.perennial.forb = imputed.traits.2 %>%
  filter(local_lifespan == "PERENNIAL") %>%
  filter(functional_group == "FORB")

imputed.traits.perennial.grass = imputed.traits.2 %>%
  filter(local_lifespan == "PERENNIAL") %>%
  filter(functional_group == "GRASS")

write.csv(imputed.traits.annual.forb, file = "./Formatted.Data/Revisions/Final.Data/imputed.traits.annual.forbs.outliersRM.csv")
write.csv(imputed.traits.perennial.forb, file = "./Formatted.Data/Revisions/Final.Data/imputed.traits.perennial.forb.outliersRM.csv")
write.csv(imputed.traits.perennial.grass, file = "./Formatted.Data/Revisions/Final.Data/imputed.traits.perennial.grass.outliersRM.csv")


#### Checking and removing outliers from imputed data without woody species ####

imputed.traits.NW= read.csv("./Formatted.Data/Revisions/imputed.traits.NW.final.csv", row.names = 1)

hist(imputed.traits.NW$leafN.final)
boxplot(imputed.traits.NW$leafN.final)
mean = mean(imputed.traits.NW$leafN.final, na.rm = TRUE)
std = sd(imputed.traits.NW$leafN.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits.NW$leafN.final[which(imputed.traits.NW$leafN.final <Tmin | imputed.traits.NW$leafN.final > Tmax)])
# removed leafN 45.18532 45.32000 46.00560 54.71597 54.71597
# percent removed 
5/794*100 #0.63%

hist(imputed.traits.NW$height.final)
boxplot(imputed.traits.NW$height.final)
mean = mean(imputed.traits.NW$height.final, na.rm = TRUE)
std = sd(imputed.traits.NW$height.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits.NW$height.final[which(imputed.traits.NW$height.final <Tmin | imputed.traits.NW$height.final > Tmax)])
# removed height 1.371600 1.383400 1.392886 1.439388 1.462500 1.493520 1.500000 1.750000 1.800000 1.800048 1.828800 1.996670 2.847733 3.800000
# percent removed 
14/794*100 #1.76%

hist(imputed.traits.NW$rootN.final)
boxplot(imputed.traits.NW$rootN.final)
mean = mean(imputed.traits.NW$rootN.final, na.rm = TRUE)
std = sd(imputed.traits.NW$rootN.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits.NW$rootN.final[which(imputed.traits.NW$rootN.final <Tmin | imputed.traits.NW$rootN.final > Tmax)])
# remove rootN 27.18327 - 39.39083
# percent removed 
16/794*100 #2.0%

hist(imputed.traits.NW$SLA.final)
boxplot(imputed.traits.NW$SLA.final)
mean = mean(imputed.traits.NW$SLA.final, na.rm = TRUE)
std = sd(imputed.traits.NW$SLA.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits.NW$SLA.final[which(imputed.traits.NW$SLA.final <Tmin | imputed.traits.NW$SLA.final > Tmax)])
# remove SLA: 49.05439 53.50000 53.50000 57.70000 61.28000 63.13000 68.90000 68.90000
# percent removed 
8/794*100 #1.01%

hist(imputed.traits.NW$root.depth.final)
boxplot(imputed.traits.NW$root.depth.final)
mean = mean(imputed.traits.NW$root.depth.final, na.rm = TRUE)
std = sd(imputed.traits.NW$root.depth.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits.NW$root.depth.final[which(imputed.traits.NW$root.depth.final <Tmin | imputed.traits.NW$root.depth.final > Tmax)])
# remove depth 6.451600 - 6.451600
# percent removed 
20/794*100 #2.52%

hist(imputed.traits.NW$RTD.final)
boxplot(imputed.traits.NW$RTD.final)
mean = mean(imputed.traits.NW$RTD.final, na.rm = TRUE)
std = sd(imputed.traits.NW$RTD.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits.NW$RTD.final[which(imputed.traits.NW$RTD.final <Tmin | imputed.traits.NW$RTD.final > Tmax)])
# remove RTD 0.6525000 - 0.9708695
# percent removed 
8/794*100 #1.01%

hist(imputed.traits.NW$SRL.final)
boxplot(imputed.traits.NW$SRL.final)
mean = mean(imputed.traits.NW$SRL.final, na.rm = TRUE)
std = sd(imputed.traits.NW$SRL.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits.NW$SRL.final[which(imputed.traits.NW$SRL.final <Tmin | imputed.traits.NW$SRL.final > Tmax)])
# remove SRL 437.3523-760.8500
# percent removed 
17/794*100 #2.14%

hist(imputed.traits.NW$rootDiam.final)
boxplot(imputed.traits.NW$rootDiam.final)
mean = mean(imputed.traits.NW$rootDiam.final, na.rm = TRUE)
std = sd(imputed.traits.NW$rootDiam.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits.NW$rootDiam.final[which(imputed.traits.NW$rootDiam.final <Tmin | imputed.traits.NW$rootDiam.final > Tmax)])
# remove diam 1.126012 1.396500 1.640000 1.850000 1.963333 4.830000
# percent removed 
5/794*100 #0.63%

hist(imputed.traits.NW$RMF.final)
boxplot(imputed.traits.NW$RMF.final)
mean = mean(imputed.traits.NW$RMF.final, na.rm = TRUE)
std = sd(imputed.traits.NW$RMF.final, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(imputed.traits.NW$RMF.final[which(imputed.traits.NW$RMF.final <Tmin | imputed.traits.NW$RMF.final > Tmax)])
# remove RMF: 0.04850791- 0.8630137
# percent removed 
11/794*100 #1.39%

# 75 total removed because some removed with multiple traits
imputed.traits.NW.2 = imputed.traits.NW %>%
  filter(round(leafN.final, 5) < 45.18532) %>% # 4
  filter(round(height.final, 5) < 1.371600) %>% # 11
  filter(round(rootN.final, 5) < 27.18327) %>% # 18
  filter(round(SLA.final, 5) < 49.05439) %>% # 8
  filter(round(root.depth.final, 5) < 6.451600) %>% # 15
  filter(round(RTD.final, 5) < 0.6525000) %>% # 11
  filter(round(SRL.final, 5) < 437.3523) %>% # 21
  filter(round(rootDiam.final, 5) < 1.126012) %>% # 6
  filter(round(RMF.final, 5) > 0.04850791 & round(RMF.final, 5) < 0.76700839) # 10

# write.csv(imputed.traits.NW.2, file = "./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.outliersRM.csv")

# split by lifespan

table(imputed.traits.NW.2$local_lifespan)
# annuals (203), biennials (10), Both (9), perennial (495), three (2)

imputed.traits.NW.perennial = imputed.traits.NW.2 %>%
  filter(local_lifespan == "PERENNIAL")
imputed.traits.NW.annual = imputed.traits.NW.2 %>%
  filter(local_lifespan == "ANNUAL")

write.csv(imputed.traits.NW.perennial, file = "./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennials.outliersRM.csv")
write.csv(imputed.traits.NW.annual, file = "./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.annuals.outliersRM.csv")

# split by functional group

table(imputed.traits.NW.2$functional_group)
# forb (396), graminoid (38), grass (238), legume (47)

imputed.traits.NW.forb = imputed.traits.NW.2 %>%
  filter(functional_group == "FORB")

imputed.traits.NW.grass = imputed.traits.NW.2 %>%
  filter(functional_group == "GRASS")

#write.csv(imputed.traits.NW.forb, file = "./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.forbs.outliersRM.csv")
#write.csv(imputed.traits.NW.grass, file = "./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.grass.outliersRM.csv")

# split by lifespan and functional group

table(imputed.traits.NW.2$local_lifespan,imputed.traits.NW.2$functional_group)
# annual forb (139), perennial forb (242), perennial grass (184)

imputed.traits.NW.annual.forb = imputed.traits.NW.2 %>%
  filter(local_lifespan == "ANNUAL") %>%
  filter(functional_group == "FORB")

imputed.traits.NW.perennial.forb = imputed.traits.NW.2 %>%
  filter(local_lifespan == "PERENNIAL") %>%
  filter(functional_group == "FORB")

imputed.traits.NW.perennial.grass = imputed.traits.NW.2 %>%
  filter(local_lifespan == "PERENNIAL") %>%
  filter(functional_group == "GRASS")

write.csv(imputed.traits.NW.annual.forb, file = "./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.annual.forbs.outliersRM.csv")
write.csv(imputed.traits.NW.perennial.forb, file = "./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennial.forb.outliersRM.csv")
write.csv(imputed.traits.NW.perennial.grass, file = "./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennial.grass.outliersRM.csv")



