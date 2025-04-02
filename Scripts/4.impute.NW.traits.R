# Imputing missing trait data using Bayesian Hierarchical Probabilistic Matrix Factorization
# Includes species after woody/tree removal.

# https://github.com/fisw10/BHPMF/blob/master/BHPMF_vignette.pdf
# https://github.com/klapierre/CoRRE_Traits/blob/master/08a_bhpmf.R
# https://github.com/ilamatos/venation_tradeoffs/blob/main/scripts/1_prepare_functional_data.R

library(devtools)
#install_github("fisw10/BHPMF")
library(BHPMF)
library(stringr)
library(tidyverse)
library(abind)
library(ggpubr)
library(cowplot)


#### Make the hierarchy file ####

traits = read.csv("./Formatted.Data/Revisions/final.data.NW.for.impute.csv", row.names = 1)

# get rid of species that are var/ssp

traits.2 = traits %>%
  filter(!Taxon %in% c("NAVARRETIA INTERTEXTA SSP. PROPINQUA"))

# see which species have all data missing

all.na = as.data.frame(apply(traits.2[,c(4:12)], 1, function(x) all(is.na(x))))
table(all.na$`apply(traits.2[, c(4:12)], 1, function(x) all(is.na(x)))`)
# 161 removed

traits.3 = traits.2[!apply(traits.2[c(4:12)], 1, function(x) all(is.na(x))), ]
# 794 individuals, 519 species

traits.3$plant_id = 1:794 # add plant_id

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
colnames(hierarchy.4)[4] = "family"

hierarchy.info = hierarchy.4

#### Make the trait file ####

# need to get rid of everything but traits
# must be in same order as hierarchy.info

trait.data.4 = traits.3[,c(4:12)]

# how much trait data is missing (out of 794)
sum(is.na(trait.data.4$leafN.mg.g))
204/794*100 # 25.69 %
sum(is.na(trait.data.4$height.m))
53/794*100 # 6.68 %
sum(is.na(trait.data.4$rootN.mg.g))
449/794*100 # 56.55 %
sum(is.na(trait.data.4$SLA_m2.kg))
80/794*100 # 10.08 %
sum(is.na(trait.data.4$root.depth_m))
314/794*100 # 39.55 %
sum(is.na(trait.data.4$rootDiam.mm))
393/794*100 # 49.49 %
sum(is.na(trait.data.4$SRL.groot.cahill.merge))
413/794*100 # 52.02 %
sum(is.na(trait.data.4$RTD.groot.cahill.merge))
417/794*100 # 52.52 %
sum(is.na(trait.data.4$RMF.g.g))
440/794*100 # 55.42 %

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

write.csv(back_trans_pars, file = "./Formatted.Data/Revisions/back_trans_pars_NW.csv")

#### Performing Gap Filling ####
# repeat imputation 50 times

for(i in 1:50){
  capture.output(GapFilling(trait.info, hierarchy.info,
                            mean.gap.filled.output.path = paste0("./Formatted.Data/Revisions/impute/NW/mean_gap_filled_",i,".csv"),
                            std.gap.filled.output.path = paste0("./Formatted.Data/Revisions/impute/NW/std_gap_filled_",i,".csv"),
                            verbose=T, rmse.plot.test.data = TRUE), file=paste0("./Formatted.Data/Revisions/impute/NW/RMSE", i,".txt"))
}

# mean_gap_filled.csv = contains the imputed mean values
# std_gap_filled.csv = contains the standard deviation for the imputed values
# RMSE = contains the average root mean square error for each imputation

# load imputed traits and std

mean.trait<-list()
std.trait<-list()
for(i in 1:50) {
  print(i)
  trt <- read.table(paste0("./Formatted.Data/Revisions/impute/NW/mean_gap_filled_",i,".csv"), row.names=NULL, header=T)
  std <- read.table(paste0("./Formatted.Data/Revisions/impute/NW/std_gap_filled_",i,".csv"), row.names=NULL, header=T)
  mean.trait[[i]] <- trt
  std.trait[[i]] <- std
}

# get mean across all means

mean.trait.2 <- abind(mean.trait, along=3)
mean.trait.3 <- apply(mean.trait.2, c(1,2), mean, na.rm=T)

# read in data for back transforming

back <- read.csv("./Formatted.Data/Revisions/back_trans_pars_NW.csv", row.names = 1)

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

# write.csv(all.data, "./Formatted.Data/Revisions/imputed.traits.NW.csv", row.names=F)


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

# write.csv(trait.std.noreplacement, "./Formatted.Data/Revisions/imputed.traits.sd.NW.csv", row.names=F)

#### Looking at ranges of real vs. imputed data ####

range(all.data$leafN.mg.g, na.rm = TRUE)
# 2.689516 54.715972
range(all.data$impute.leafN.mg.g)
# 4.678917 53.870502
range(all.data$height.m, na.rm = TRUE)
# 0.009 3.800
range(all.data$impute.height.m)
#  0.009432988 2.403367027
range(all.data$rootN.mg.g, na.rm = TRUE)
# 0.38000 39.26274
range(all.data$impute.rootN.mg.g)
# 0.4166559 39.3908306, above range
range(all.data$SLA_m2.kg, na.rm = TRUE)
# 2.094693 68.900000
range(all.data$impute.SLA_m2.kg)
# 3.070733 64.607813
range(all.data$root.depth_m, na.rm = TRUE)
# 0.048000 3.709333
range(all.data$impute.root.depth_m)
# 0.08080788 2.89370345
range(all.data$rootDiam.mm, na.rm = TRUE)
# 0.0784 4.8300
range(all.data$impute.rootDiam.mm)
# 0.08439273 2.38465479
range(all.data$SRL.groot.cahill.merge, na.rm = TRUE)
# 2.49723 760.85000
range(all.data$impute.SRL.groot.cahill.merge)
# 3.052451 675.086085
range(all.data$RTD.groot.cahill.merge, na.rm = TRUE)
# 0.008390413 0.700029683
range(all.data$impute.RTD.groot.cahill.merge)
# 0.01542964 0.97086948; above range
range(all.data$RMF.g.g, na.rm = TRUE)
# 0.04850791 0.86301370
range(all.data$impute.RMF.g.g)
# 0.08000226 0.82353118

### Plot original vs. imputed trait mean values ####
# Validate imputation

all.data = read.csv("./Formatted.Data/Revisions/imputed.traits.NW.csv")

leafN = ggplot(all.data, aes(x = leafN.mg.g, y = impute.leafN.mg.g))+
  geom_point(color = "gray50")+
  geom_smooth(method = "lm", color = "black")+
  #stat_regline_equation(label.y = 55, aes(label = after_stat(eq.label)), size = 5) +
  stat_regline_equation(label.y = 53, aes(label = after_stat(rr.label)), size = 5)+
  theme_classic(base_size = 15)+
  labs(y = "Imputed Leaf N", x = "Observed Leaf N")
leafN # R2 = 0.97

height = ggplot(all.data, aes(x = height.m, y = impute.height.m))+
  geom_point(color = "gray50")+
  geom_smooth(method = "lm", color = "black")+
  #stat_regline_equation(label.y = 3.2, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 3.3, aes(label = ..rr.label..), size = 3)+
  theme_classic(base_size = 15)+
  labs(y = "Imputed Height", x = "Observed Height")
height # R2 = 0.93

rootN = ggplot(all.data, aes(x = rootN.mg.g, y = impute.rootN.mg.g))+
  geom_point(color = "gray50")+
  geom_smooth(method = "lm", color = "black")+
  #stat_regline_equation(label.y = 55, aes(label = after_stat(eq.label)), size = 5) +
  stat_regline_equation(label.y = 38, aes(label = after_stat(rr.label)), size = 5)+
  theme_classic(base_size = 15)+
  labs(y = "Imputed Root N", x = "Observed Root N")
rootN # R2 = 0.99

SLA = ggplot(all.data, aes(x = SLA_m2.kg, y = impute.SLA_m2.kg))+
  geom_point(color = "gray50")+
  geom_smooth(method = "lm", color = "black")+
  #stat_regline_equation(label.y = 40, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 65, aes(label = ..rr.label..), size = 3)+
  theme_classic(base_size = 15)+
  labs(y = "Imputed SLA", x = "Observed SLA")
SLA # R2 = 0.94

root.depth = ggplot(all.data, aes(x = root.depth_m, y = impute.root.depth_m))+
  geom_point(color = "gray50")+
  geom_smooth(method = "lm", color = "black")+
  #stat_regline_equation(label.y = 7, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 3.3, aes(label = ..rr.label..), size = 3)+
  theme_classic(base_size = 15)+
  labs(y = "Imputed Rooting Depth", x = "Observed Rooting Depth")
root.depth # R2 = 0.94

rootDiam = ggplot(all.data, aes(x = rootDiam.mm, y = impute.rootDiam.mm))+
  geom_point(color = "gray50")+
  geom_smooth(method = "lm", color = "black")+
  #stat_regline_equation(label.y = 5, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 3.3, aes(label = ..rr.label..), size = 3)+
  theme_classic(base_size = 15)+
  labs(y = "Imputed Root Diameter", x = "Observed Root Diameter")
rootDiam # R2 = 0.89

SRL = ggplot(all.data, aes(x = SRL.groot.cahill.merge, y = impute.SRL.groot.cahill.merge))+
  geom_point(color = "gray50")+
  geom_smooth(method = "lm", color = "black")+
  #stat_regline_equation(label.y = 760, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 675, aes(label = ..rr.label..), size = 3)+
  theme_classic(base_size = 15)+
  labs(y = "Imputed SRL", x = "Observed SRL")
SRL # R2 = 0.98

RTD = ggplot(all.data, aes(x = RTD.groot.cahill.merge, y = impute.RTD.groot.cahill.merge))+
  geom_point(color = "gray50")+
  geom_smooth(method = "lm", color = "black")+
  #stat_regline_equation(label.y = 1.1, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 0.85, aes(label = ..rr.label..), size = 3)+
  theme_classic(base_size = 15)+
  labs(y = "Imputed RTD", x = "Observed RTD")
RTD # R2 = 0.99

RMF = ggplot(all.data, aes(x = RMF.g.g, y = impute.RMF.g.g))+
  geom_point(color = "gray50")+
  geom_smooth(method = "lm", color = "black")+
  #stat_regline_equation(label.y = 0.85, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 0.83, aes(label = ..rr.label..), size = 3)+
  theme_classic(base_size = 15)+
  labs(y = "Imputed RMF", x = "Observed RMF")
RMF # R2 = 0.99

impute.eval.plot = plot_grid(leafN,height,rootN,SLA,RMF,root.depth,RTD,SRL,rootDiam, 
                             labels = c("A","B","C","D","E","F","G","H","I"))
impute.eval.plot

ggsave("./Plots/impute.traits.eval.plot.png", height = 11, width = 8.5)
### Replace NAs in the original trait data with the imputed values

imputed.traits = read.csv("./Formatted.Data/Revisions/imputed.traits.NW.csv")

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

write.csv(imput_final, file = "./Formatted.Data/Revisions/imputed.traits.NW.final.csv")


