# Imputing missing trait data using Bayesian Hierarchical Probabilistic Matrix Factorization
# Same method used to build CoRRE database
# https://github.com/fisw10/BHPMF/blob/master/BHPMF_vignette.pdf
# https://github.com/klapierre/CoRRE_Traits/blob/master/08a_bhpmf.R

library(devtools)
install_github("fisw10/BHPMF")
library(BHPMF)
library(stringr)
library(tidyverse)

# Set temporary directory

tmp.dir = dirname("./Formatted.Data/Revisions/tmp/")

#### Make the hierarchy file ####

traits = read.csv("./Formatted.Data/Revisions/final.data.for.impute.csv", row.names = 1)

# get rid of species that are var/ssp

traits.2 = traits %>%
  filter(!Taxon %in% c("HELIANTHEMUM NUMMULARIUM var. GRANDIFLORUM","NAVARRETIA INTERTEXTA SSP. PROPINQUA"))

# see which species have all data missing

all.na = as.data.frame(apply(traits.2[,c(4:12)], 1, function(x) all(is.na(x))))
table(all.na$`apply(traits.2[, c(4:12)], 1, function(x) all(is.na(x)))`)
# 153 removed

traits.3 = traits.2[!apply(traits.2[c(4:12)], 1, function(x) all(is.na(x))), ]
# 924 individuals, 614 species

traits.3$plant_id = 1:924 # add plant_id

hierarchy = traits.3[,c(17,2)]
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

hierarchy.4[128,4] = "Compositae"
colnames(hierarchy.4)[4] = "family"

hierarchy.info = hierarchy.4

#### Make the trait file ####

# need to get rid of everything but traits
# must be in same order as hierarchy.info

trait.data.4 = traits.3[,c(4:12)]

# how much trait data is missing (out of 924)
sum(is.na(trait.data.4$leafN.mg.g))
275/924*100 # 29.8 %
sum(is.na(trait.data.4$height.m))
43/924*100 # 4.65 %
sum(is.na(trait.data.4$rootN.mg.g))
536/924*100 # 58.0 %
sum(is.na(trait.data.4$SLA_m2.kg))
66/924*100 # 7.14 %
sum(is.na(trait.data.4$root.depth_m))
378/924*100 # 40.9 %
sum(is.na(trait.data.4$rootDiam.mm))
474/924*100 # 51.3 %
sum(is.na(trait.data.4$SRL.groot.cahill.merge))
489/924*100 # 52.9 %
sum(is.na(trait.data.4$RTD.groot.cahill.merge))
510/924*100 # 55.2 %
sum(is.na(trait.data.4$RMF.g.g))
529/924*100 # 57.3 %

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

GapFilling(trait.info, hierarchy.info,
           tuning = TRUE,
           mean.gap.filled.output.path = paste0(tmp.dir,"/mean_gap_filled_",i,".txt"),
           std.gap.filled.output.path = paste0(tmp.dir,"/std_gap_filled_",i,".txt"),
           tmp.dir = tmp.dir, verbose=F, rmse.plot.test.data = TRUE)

#### calculate cross validation RMSE ####

out1 = CalculateCvRmse(trait.info,hierarchy.info)
avg.rmse = out1$avg.rmse
std.rmse = out1$std.rmse


