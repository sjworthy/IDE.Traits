# Code to summarize and organize traits from TRY

# Traits
# 6: root rooting depth (m)
# 14: leaf nitrogen content (mg/g)
# 80: root nitrogen content (mg/g)
# 83: root diameter (mm)
# 614: fine root length per fine root dry mass (specific fine root length) (cm/g)
# 3106: Plant height vegetative (m)
# 3115: Leaf area per leaf dry mass (SLA): petiole excluded mm2/mg
# 3116: Leaf area per leaf dry mass (SLA): petiole included mm2/mg
# 3117: Leaf area per leaf dry mass (SLA): undefined if petiole in-or-excluded mm2/mg

library(stringr)
library(dplyr)
library(tidyverse)

try.data = read.csv("./Raw.Data/TRYTRAIT.RAW.csv")

# remove rows that don't include trait data - most include metadata we don't need
try.data.2 <- subset(try.data, try.data$TraitName != "")

# Select the traits we are actually working with

try.data.3 = try.data.2 %>%
  filter(TraitID %in% c(6,14,80,83,614,3106,3115,3116,3117))

# remove values that correspond to maximum, minimum, high, low
try.data.4 = try.data.3 %>%
  filter(ValueKindName %in% c("Best estimate", "Mean","Median","Single","Site specific mean","Species mean"))

#### Summarize Traits ####

summary_data <- try.data.4 %>% 
  group_by(DatasetID, AccSpeciesName, TraitName) %>% 
  summarise(mean_value=mean(StdValue, na.rm = T),
            .groups = 'drop')

# summarize means across study
summary <- summary_data %>% group_by(AccSpeciesName, TraitName) %>% 
  summarise(mean_value=mean(mean_value, na.rm = T),
            .groups = 'drop')

summary_wide <- pivot_wider(data = summary, names_from = TraitName, values_from = mean_value)

colnames(summary_wide) = c("Species","Height.m","Root.depth.m","SLA.no.petiole.mm2.mg","SLA.unknown.petiole.mm2.mg",
                           "SLA.petiole.mm2.mg","leafN.mg.g","rootN.mg.g","SRL.fine.cm.g","root.diam.mm")

summary_wide = summary_wide %>%
  group_by(Species) %>%
  mutate(SLA.mg.g = rowMeans(cbind(SLA.no.petiole.mm2.mg,SLA.unknown.petiole.mm2.mg,SLA.petiole.mm2.mg), na.rm = TRUE))

# write.csv(summary_wide, file="./Formatted.Data/Revisions/TRY.traits.summary.csv")


