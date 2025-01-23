# Imputing missing trait data using Bayesian Hierarchical Probabilistic Matrix Factorization
# includes all non-woody and woody species, full dataset

# https://github.com/fisw10/BHPMF/blob/master/BHPMF_vignette.pdf
# https://github.com/klapierre/CoRRE_Traits/blob/master/08a_bhpmf.R
# https://github.com/ilamatos/venation_tradeoffs/blob/main/scripts/1_prepare_functional_data.R

library(devtools)
#install_github("fisw10/BHPMF")
library(BHPMF)
library(stringr)
library(tidyverse)
library(abind)


#### Make the hierarchy file ####

traits = read.csv("./Formatted.Data/Revisions/final.data.for.impute.csv", row.names = 1)

# get rid of species that are var/ssp

traits.2 = traits %>%
  filter(!Taxon %in% c("HELIANTHEMUM NUMMULARIUM var. GRANDIFLORUM","NAVARRETIA INTERTEXTA SSP. PROPINQUA"))

# see which species have all data missing

all.na = as.data.frame(apply(traits.2[,c(4:12)], 1, function(x) all(is.na(x))))
table(all.na$`apply(traits.2[, c(4:12)], 1, function(x) all(is.na(x)))`)
# 170 removed

traits.3 = traits.2[!apply(traits.2[c(4:12)], 1, function(x) all(is.na(x))), ]
# 907 individuals, 601 species

traits.3$plant_id = 1:907 # add plant_id

hierarchy = traits.3[,c(16,2)]
colnames(hierarchy)[2] = "species"
hierarchy$species = str_to_sentence(hierarchy$species) # make sentence case
hierarchy$species.2 = hierarchy$species

hierarchy.2 = hierarchy %>%
  separate(species.2, c("genus", "species.3")) # splitting into genus and latin binomial

hierarchy.3 = hierarchy.2[,c(1:3)]            

# get family 

data=read.csv("./Raw.Data/IDE_cover_2023-01-02.csv") %>%
  select(Family, Taxon)

colnames(data)[2] = "species"
data$species = str_to_sentence(data$species)
data.2 = distinct(data)

hierarchy.4 = left_join(hierarchy.3,data.2)

hierarchy.4[126,4] = "Compositae"
colnames(hierarchy.4)[4] = "family"

hierarchy.info = hierarchy.4

#### Make the trait file ####

# need to get rid of everything but traits
# must be in same order as hierarchy.info

trait.data.4 = traits.3[,c(4:12)]

# how much trait data is missing (out of 907)
sum(is.na(trait.data.4$leafN.mg.g))
232/907*100 # 25.58 %
sum(is.na(trait.data.4$height.m))
57/907*100 # 6.28 %
sum(is.na(trait.data.4$rootN.mg.g))
519/907*100 # 57.22 %
sum(is.na(trait.data.4$SLA_m2.kg))
97/907*100 # 10.69 %
sum(is.na(trait.data.4$root.depth_m))
361/907*100 # 39.80 %
sum(is.na(trait.data.4$rootDiam.mm))
457/907*100 # 50.39 %
sum(is.na(trait.data.4$SRL.groot.cahill.merge))
473/907*100 # 52.14 %
sum(is.na(trait.data.4$RTD.groot.cahill.merge))
492/907*100 # 54.24 %
sum(is.na(trait.data.4$RMF.g.g))
512/907*100 # 56.45 %

trait.info = as.matrix(trait.data.4)

#check if both datasets are equal

nrow(hierarchy.info) == nrow(trait.info)

#### z % log transform ####
back_trans_pars <- list()
rm_col <- c()
for(i in 1:ncol(trait.info)){
  x <- trait.info[,i] # goes through the columns
  min_x <- min(x,na.rm = T) # takes the min of each column
  if(min_x < 0.00000000001){
    x <- x - min_x + 1 # make this optional if min x is neg
  }
  logx <- log10(x)
  mlogx <- mean(logx, na.rm = T)
  slogx <- sd(logx, na.rm = T)
  x <- (logx - mlogx)/slogx # Z transformation
  back_trans_pars[[i]] <- list(min_x = min_x,
                               mlogx = mlogx,
                               slogx = slogx)
  trait.info[,i] <- x
}

write.csv(back_trans_pars, file = "./Formatted.Data/Revisions/back_trans_pars.csv")

#### Performing Gap Filling ####
# repeat imputation 50 times

for(i in 1:50){
  capture.output(GapFilling(trait.info, hierarchy.info,
             mean.gap.filled.output.path = paste0("./Formatted.Data/Revisions/impute/mean_gap_filled_",i,".csv"),
             std.gap.filled.output.path = paste0("./Formatted.Data/Revisions/impute/std_gap_filled_",i,".csv"),
             verbose=T, rmse.plot.test.data = TRUE), file=paste0("./Formatted.Data/Revisions/impute/RMSE", i,".txt"))
}

# mean_gap_filled.csv = contains the imputed mean values
# std_gap_filled.csv = contains the standard deviation for the imputed values
# RMSE = contains the average root mean square error for each imputation
  
# load imputed traits and std

mean.trait<-list()
std.trait<-list()
for(i in 1:50) {
  print(i)
  trt <- read.table(paste0("./Formatted.Data/Revisions/impute/mean_gap_filled_",i,".csv"), row.names=NULL, header=T)
  std <- read.table(paste0("./Formatted.Data/Revisions/impute/std_gap_filled_",i,".csv"), row.names=NULL, header=T)
  mean.trait[[i]] <- trt
  std.trait[[i]] <- std
}

# get mean across all means

mean.trait.2 <- abind(mean.trait, along=3)
mean.trait.3 <- apply(mean.trait.2, c(1,2), mean, na.rm=T)

# read in data for back transforming

back <- read.csv("./Formatted.Data/Revisions/back_trans_pars.csv", row.names = 1)

# don't replace original values:
trait.info.noreplacement <- as.data.frame(mean.trait.3)

o <- 1 #to select the appropriate columns:
for(i in 1:ncol(trait.info.noreplacement)){
  
  #recover values:
  min_x <- back[1,o]
  mlogx <- back[1,o+1]
  slogx <- back[1,o+2]
  
  #back transform:
  x <- trait.info.noreplacement[,i] # goes through the columns
  logx <- (x*slogx) + mlogx
  b <- 10^logx
  
  #for negative values
  if(min_x < 0.00000000001){
    b <- b + min_x - 1 # make this optional if min x is neg
  }
  
  trait.info.noreplacement[,i] <- b
  o <- o+3
}

# merge imputed traits with original data
colnames(trait.info.noreplacement) = c("impute.leafN.mg.g","impute.height.m","impute.rootN.mg.g","impute.SLA_m2.kg","impute.root.depth_m",
                                       "impute.rootDiam.mm","impute.SRL.groot.cahill.merge","impute.RTD.groot.cahill.merge","impute.RMF.g.g")

all.data = cbind(traits.3, trait.info.noreplacement)

# write.csv(all.data, "./Formatted.Data/Revisions/imputed.traits.csv", row.names=F)


#### get mean across all std ####

std.trait.2 <- abind(std.trait, along=3)
std.trait.3 <- apply(std.trait.2, c(1,2), mean, na.rm=T)

#don't replace original values:
trait.std.noreplacement <- as.data.frame(std.trait.3)

o <- 1 #to select the appropriate columns:
for(i in 1:ncol(trait.std.noreplacement)){
  
  #recover values:
  min_x <- back[1,o]
  mlogx <- back[1,o+1]
  slogx <- back[1,o+2]
  
  #back transform:
  x <- trait.std.noreplacement[,i] # goes through the columns
  logx <- (x*slogx) + mlogx
  b <- 10^logx
  
  #for negative values
  if(min_x < 0.00000000001){
    b <- b + min_x - 1 # make this optional if min x is neg
  }
  
  trait.std.noreplacement[,i] <- b
  o <- o+3
}

# write.csv(trait.std.noreplacement, "./Formatted.Data/Revisions/imputed.traits.sd.csv", row.names=F)

#### Looking at ranges of real vs. imputed data ####

range(all.data$leafN.mg.g, na.rm = TRUE)
# 2.689516 54.715972
range(all.data$impute.leafN.mg.g)
# 4.15724 81.93849; above range
range(all.data$height.m, na.rm = TRUE)
# 0.00900 32.57366
range(all.data$impute.height.m)
#  0.009586053 25.631534692
range(all.data$rootN.mg.g, na.rm = TRUE)
# 0.38000 39.26274
range(all.data$impute.rootN.mg.g)
# 0.5347225 32.9601585
range(all.data$SLA_m2.kg, na.rm = TRUE)
# 2.094693 68.900000
range(all.data$impute.SLA_m2.kg)
# 2.305738 58.579231
range(all.data$root.depth_m, na.rm = TRUE)
# 0.0070 6.4516
range(all.data$impute.root.depth_m)
# 0.007867208 6.248630389
range(all.data$rootDiam.mm, na.rm = TRUE)
# 0.0784 4.8300
range(all.data$impute.rootDiam.mm)
# 0.0829587 3.5804349
range(all.data$SRL.groot.cahill.merge, na.rm = TRUE)
# 2.49723 760.85000
range(all.data$impute.SRL.groot.cahill.merge)
# 5.174528 808.029806; above range
range(all.data$RTD.groot.cahill.merge, na.rm = TRUE)
# 0.008390413 1.194550386
range(all.data$impute.RTD.groot.cahill.merge)
# 0.01664353 1.47570114; above range
range(all.data$RMF.g.g, na.rm = TRUE)
# 0.04850791 0.86301370
range(all.data$impute.RMF.g.g)
# 0.07712438 0.81268793

### Plot original vs. imputed trait mean values ####
# Validate imputation

leafN = ggplot(all.data, aes(x = leafN.mg.g, y = impute.leafN.mg.g))+
  geom_point()+
  geom_smooth(method = "lm", color = "black")+
  stat_regline_equation(label.y = 55, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 53, aes(label = ..rr.label..), size = 3)+
  theme_classic()
leafN # R2 = 0.97

height = ggplot(all.data, aes(x = height.m, y = impute.height.m))+
  geom_point()+
  geom_smooth(method = "lm", color = "black")+
  stat_regline_equation(label.y = 30, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 28, aes(label = ..rr.label..), size = 3)+
  theme_classic()
height # R2 = 0.98

SLA = ggplot(all.data, aes(x = SLA_m2.kg, y = impute.SLA_m2.kg))+
  geom_point()+
  geom_smooth(method = "lm", color = "black")+
  stat_regline_equation(label.y = 40, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 38, aes(label = ..rr.label..), size = 3)+
  theme_classic()
SLA # R2 = 0.93

root.depth = ggplot(all.data, aes(x = root.depth_m, y = impute.root.depth_m))+
  geom_point()+
  geom_smooth(method = "lm", color = "black")+
  stat_regline_equation(label.y = 7, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 6.7, aes(label = ..rr.label..), size = 3)+
  theme_classic()
root.depth # R2 = 0.96

rootDiam = ggplot(all.data, aes(x = rootDiam.mm, y = impute.rootDiam.mm))+
  geom_point()+
  geom_smooth(method = "lm", color = "black")+
  stat_regline_equation(label.y = 5, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 4.7, aes(label = ..rr.label..), size = 3)+
  theme_classic()
rootDiam # R2 = 0.96

SRL = ggplot(all.data, aes(x = SRL.groot.cahill.merge, y = impute.SRL.groot.cahill.merge))+
  geom_point()+
  geom_smooth(method = "lm", color = "black")+
  stat_regline_equation(label.y = 760, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 740, aes(label = ..rr.label..), size = 3)+
  theme_classic()
SRL # R2 = 0.99

RTD = ggplot(all.data, aes(x = RTD.groot.cahill.merge, y = impute.RTD.groot.cahill.merge))+
  geom_point()+
  geom_smooth(method = "lm", color = "black")+
  stat_regline_equation(label.y = 1.1, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 1, aes(label = ..rr.label..), size = 3)+
  theme_classic()
RTD # R2 = 0.99

RMF = ggplot(all.data, aes(x = RMF.g.g, y = impute.RMF.g.g))+
  geom_point()+
  geom_smooth(method = "lm", color = "black")+
  stat_regline_equation(label.y = 0.85, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 0.8, aes(label = ..rr.label..), size = 3)+
  theme_classic()
RMF # R2 = 0.99

### Replace NAs in the original trait data with the imputed values

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.csv")

imput_final<-imputed.traits%>%
  mutate(leafN.final = coalesce(leafN.mg.g,impute.leafN.mg.g),
         height.final = coalesce(height.m,impute.height.m),
         rootN.final = coalesce(rootN.mg.g,impute.rootN.mg.g),
         SLA.final = coalesce(SLA_m2.kg,impute.SLA_m2.kg),
         root.depth.final = coalesce(root.depth_m,impute.root.depth_m),
         rootDiam.final = coalesce(rootDiam.mm,impute.rootDiam.mm),
         SRL.final = coalesce(SRL.groot.cahill.merge,impute.SRL.groot.cahill.merge),
         RTD.final = coalesce(RTD.groot.cahill.merge,impute.RTD.groot.cahill.merge),
         RMF.final = coalesce(RMF.g.g,impute.RMF.g.g))

sum(is.na(imput_final[,c(26:34)])) # no NAs in the final imputed trait dataset

write.csv(imput_final, file = "./Formatted.Data/Revisions/imputed.traits.final.csv")


