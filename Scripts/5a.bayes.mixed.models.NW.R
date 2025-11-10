## Script to evaluate linear mixed effects models of cover change and traits
# cover outliers removed

# load libraries 
library(tidyverse)
library(brms)
library(emmeans)
library(performance)
library(cowplot)

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
imputed.NW.annual = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.annuals.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.NW.perennial = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennials.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.NW.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.forbs.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.NW.grass = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.grass.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.NW.annual.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.annual.forbs.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.NW.perennial.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennial.forb.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.NW.perennial.grass = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennial.grass.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)

imputed.NW.forb.legume = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.forbs.legume.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.NW.graminoid = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.graminoid.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.NW.annual.forb.legume = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.annual.forb.legume.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.NW.perennial.forb.legume = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennial.forb.legume.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)
imputed.NW.perennial.graminoid = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennial.graminoid.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)


#### imputed traits model NW ####

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

imputed.traits.NW.model = readRDS("./Results/all.imputed.traits.no_woody.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.traits.NW.model, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.traits.NW.model, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals: All Species",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.traits.NW.model) # p = 0.753
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values: All Species") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.traits.NW.model)

imputed.NW.traits.height.leafN = brm(cover.change ~ height.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                     family = gaussian(),
                                     prior = priors,
                                     data = imputed.NW)

summary(imputed.NW.traits.height.leafN)
# height x leafN

#saveRDS(imputed.NW.traits.height.leafN, file = "./Results/all.imputed.traits.NW.height.leafN.rds")
bayes_R2(imputed.NW.traits.height.leafN)
# R2 0.04979303 0.02318114 0.01594568 0.1060524

conditional_effects(imputed.NW.traits.height.leafN)
# positive cover change with taller height and higher leafN

imputed.NW.traits.height.leafN = readRDS("./Results/all.imputed.traits.NW.height.leafN.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.traits.height.leafN, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.traits.height.leafN, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.traits.height.leafN) # p = 0.676
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.traits.height.leafN)

imputed.NW.traits.depth.leafN = brm(cover.change ~ root.depth.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                    family = gaussian(),
                                    prior = priors,
                                    data = imputed.NW)

summary(imputed.NW.traits.depth.leafN)
# leafN

saveRDS(imputed.NW.traits.depth.leafN, file = "./Results/all.imputed.NW.traits.depth.leafN.rds")
bayes_R2(imputed.NW.traits.depth.leafN)
# 0.04902346 0.02481235 0.01394253 0.1094119

imputed.NW.traits.depth.leafN = readRDS("./Results/all.imputed.NW.traits.depth.leafN.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.traits.depth.leafN, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.traits.depth.leafN, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.traits.depth.leafN) # p = 0.713
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.traits.depth.leafN)

imputed.NW.traits.RTD.SRL = brm(cover.change ~ RTD.final*SRL.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                family = gaussian(),
                                prior = priors,
                                data = imputed.NW)

summary(imputed.NW.traits.RTD.SRL)
# nothing significant

saveRDS(imputed.NW.traits.RTD.SRL, file = "./Results/all.imputed.NW.traits.RTD.SRL.rds")
bayes_R2(imputed.NW.traits.RTD.SRL)
# R2 0.04620861 0.02366104 0.01308091 0.1012018

imputed.NW.traits.RTD.SRL = readRDS("./Results/all.imputed.NW.traits.RTD.SRL.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.traits.RTD.SRL, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.traits.RTD.SRL, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.traits.RTD.SRL) # p = 0.782
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.traits.RTD.SRL)

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

imputed.NW.traits.leafN.RMF = readRDS("./Results/all.imputed.NW.traits.leafN.RMF.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.traits.leafN.RMF, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.traits.leafN.RMF, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.traits.leafN.RMF) # p = 0.682
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.traits.leafN.RMF)

#### All-species Environment ####

priors <- c(prior(normal(0, 10), class = b))

imputed.NW.traits.height.MAP.DSI = brm(cover.change ~ height.final*mean.MAP + height.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                   family = gaussian(),
                                   prior = priors,
                                   data = imputed.NW)

summary(imputed.NW.traits.height.MAP.DSI)

# saveRDS(imputed.NW.traits.height.MAP.DSI, file = "./Results/all.imputed.NW.traits.height.MAP.DSI.rds")

imputed.NW.traits.height.MAP.DSI = readRDS("./Results/all.imputed.NW.traits.height.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.traits.height.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.traits.height.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.traits.height.MAP.DSI) # p = 0.791
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.traits.height.MAP.DSI)

imputed.NW.traits.leafN.MAP.DSI = brm(cover.change ~ leafN.final*mean.MAP + leafN.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                  family = gaussian(),
                                  prior = priors,
                                  data = imputed.NW)

summary(imputed.NW.traits.leafN.MAP.DSI)
# leafN 
# leafN x DSI

#saveRDS(imputed.NW.traits.leafN.MAP.DSI, file = "./Results/all.imputed.NW.traits.leafN.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.leafN.MAP.DSI)
# R2 0.06092939 0.0258172 0.02256498 0.1241577

conditional_effects(imputed.NW.traits.leafN.MAP.DSI)
# higher leafN x less drought (higher DSI) or lower leafN x more drought (lower DSI)

imputed.NW.traits.leafN.MAP.DSI = readRDS("./Results/all.imputed.NW.traits.leafN.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.traits.leafN.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.traits.leafN.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.traits.leafN.MAP.DSI) # p = 0.712
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.traits.leafN.MAP.DSI)

imputed.NW.traits.root.depth.MAP.DSI = brm(cover.change ~ root.depth.final*mean.MAP + root.depth.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                       family = gaussian(),
                                       prior = priors,
                                       data = imputed.NW)

summary(imputed.NW.traits.root.depth.MAP.DSI)
#saveRDS(imputed.NW.traits.root.depth.MAP.DSI, file = "./Results/all.imputed.NW.traits.depth.MAP.DSI.rds")

imputed.NW.traits.root.depth.MAP.DSI = readRDS("./Results/all.imputed.NW.traits.depth.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.traits.root.depth.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.traits.root.depth.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.traits.root.depth.MAP.DSI) # p = 0.815
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.traits.root.depth.MAP.DSI)

imputed.NW.traits.RTD.MAP.DSI = brm(cover.change ~ RTD.final*mean.MAP + RTD.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                family = gaussian(),
                                prior = priors,
                                data = imputed.NW)

summary(imputed.NW.traits.RTD.MAP.DSI)
# saveRDS(imputed.NW.traits.RTD.MAP.DSI, file = "./Results/all.imputed.NW.traits.RTD.MAP.DSI.rds")

imputed.NW.traits.RTD.MAP.DSI = readRDS("./Results/all.imputed.NW.traits.RTD.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.traits.RTD.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.traits.RTD.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.traits.RTD.MAP.DSI) # p = 0.779
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.traits.RTD.MAP.DSI)

imputed.NW.traits.SRL.MAP.DSI = brm(cover.change ~ SRL.final*mean.MAP + SRL.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                family = gaussian(),
                                prior = priors,
                                data = imputed.NW)

summary(imputed.NW.traits.SRL.MAP.DSI)
#saveRDS(imputed.NW.traits.SRL.MAP.DSI, file = "./Results/all.imputed.NW.traits.SRL.MAP.DSI.rds")

imputed.NW.traits.SRL.MAP.DSI = readRDS("./Results/all.imputed.NW.traits.SRL.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.traits.SRL.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.traits.SRL.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.traits.SRL.MAP.DSI) # p = 0.830
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.traits.SRL.MAP.DSI)

imputed.NW.traits.RMF.MAP.DSI = brm(cover.change ~ RMF.final*mean.MAP + RMF.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                family = gaussian(),
                                prior = priors,
                                data = imputed.NW)

summary(imputed.NW.traits.RMF.MAP.DSI)
#saveRDS(imputed.NW.traits.RMF.MAP.DSI, file = "./Results/all.imputed.NW.traits.RMF.MAP.DSI.rds")

imputed.NW.traits.RMF.MAP.DSI = readRDS("./Results/all.imputed.NW.traits.RMF.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.traits.RMF.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.traits.RMF.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.traits.RMF.MAP.DSI) # p = 0.811
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.traits.RMF.MAP.DSI)

#### imputed ANNUAL traits lifespan model NW ####

priors <- c(prior(normal(0, 10), class = b))

annual.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                               SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                             family = gaussian(),
                             prior = priors,
                             data = imputed.NW.annual)

summary(annual.traits.NW.model)
# MAP
#saveRDS(annual.traits.NW.model, file = "./Results/annual.imputed.traits.no_woody.rds")
bayes_R2(annual.traits.NW.model)
# R2 0.1542627 0.04448325 0.07546093 0.2476329

annual.traits.NW.model = readRDS("./Results/annual.imputed.traits.no_woody.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(annual.traits.NW.model, summary = TRUE))
fitted_vals <- as.data.frame(fitted(annual.traits.NW.model, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(annual.traits.NW.model) # p = 0.831
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(annual.traits.NW.model)

imputed.NW.annual.traits.height.leafN = brm(cover.change ~ height.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                            family = gaussian(),
                                            prior = priors,
                                            data = imputed.NW.annual)

summary(imputed.NW.annual.traits.height.leafN)
# MAP
#saveRDS(imputed.NW.annual.traits.height.leafN, file = "./Results/annual.imputed.traits.NW.height.leafN.rds")

imputed.NW.annual.traits.height.leafN = readRDS("./Results/annual.imputed.traits.NW.height.leafN.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.annual.traits.height.leafN, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.annual.traits.height.leafN, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.annual.traits.height.leafN) # p = 0.737
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.annual.traits.height.leafN)

imputed.NW.annual.traits.depth.leafN = brm(cover.change ~ root.depth.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                           family = gaussian(),
                                           prior = priors,
                                           data = imputed.NW.annual)

summary(imputed.NW.annual.traits.depth.leafN)
# MAP is positive, significant
#saveRDS(imputed.NW.annual.traits.depth.leafN, file = "./Results/annual.imputed.traits.NW.depth.leafN.rds")

imputed.NW.annual.traits.depth.leafN = readRDS("./Results/annual.imputed.traits.NW.depth.leafN.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.annual.traits.depth.leafN, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.annual.traits.depth.leafN, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.annual.traits.depth.leafN) # p = 0.845
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.annual.traits.depth.leafN)

imputed.NW.annual.traits.RTD.SRL = brm(cover.change ~ RTD.final*SRL.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                       family = gaussian(),
                                       prior = priors,
                                       data = imputed.NW.annual)

summary(imputed.NW.annual.traits.RTD.SRL)
# nothing significant
# saveRDS(imputed.NW.annual.traits.RTD.SRL, file = "./Results/annual.imputed.traits.NW.RTD.SRL.rds")

imputed.NW.annual.traits.RTD.SRL = readRDS("./Results/annual.imputed.traits.NW.RTD.SRL.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.annual.traits.RTD.SRL, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.annual.traits.RTD.SRL, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.annual.traits.RTD.SRL) # p = 0.889
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.annual.traits.RTD.SRL)

imputed.NW.annual.traits.leafN.RMF = brm(cover.change ~ leafN.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                         family = gaussian(),
                                         prior = priors,
                                         data = imputed.NW.annual)

summary(imputed.NW.annual.traits.leafN.RMF)
# MAP significant positive
# leafN:RMF significant positive

#saveRDS(imputed.NW.annual.traits.leafN.RMF, file = "./Results/annual.imputed.NW.traits.leafN.RMF.rds")
bayes_R2(imputed.NW.annual.traits.leafN.RMF)
# R2 0.1310284 0.04710103 0.0522819 0.2401864

conditional_effects(imputed.NW.annual.traits.leafN.RMF)
# higher cover change with high leafN and high RMF or low leafN and low RMF

imputed.NW.annual.traits.leafN.RMF = readRDS("./Results/annual.imputed.NW.traits.leafN.RMF.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.annual.traits.leafN.RMF, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.annual.traits.leafN.RMF, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.annual.traits.leafN.RMF) # p = 0.765
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.annual.traits.leafN.RMF)

#### Annuals Environment ####

priors <- c(prior(normal(0, 10), class = b))

imputed.NW.annual.traits.height.MAP.DSI = brm(cover.change ~ height.final*mean.MAP + height.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                              family = gaussian(),
                                              prior = priors,
                                              data = imputed.NW.annual)

summary(imputed.NW.annual.traits.height.MAP.DSI)
# height x MAP

#saveRDS(imputed.NW.annual.traits.height.MAP.DSI, file = "./Results/imputed.NW.annual.traits.height.MAP.DSI.rds")
bayes_R2(imputed.NW.annual.traits.height.MAP.DSI)
# R2 0.1208656 0.04510458 0.04694877 0.2220583

imputed.NW.annual.traits.height.MAP.DSI = readRDS("./Results/imputed.NW.annual.traits.height.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.annual.traits.height.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.annual.traits.height.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.annual.traits.height.MAP.DSI) # p = 0.877
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.annual.traits.height.MAP.DSI)

conditional_effects(imputed.NW.annual.traits.height.MAP.DSI)
# higher height x higher MAP

imputed.NW.annual.traits.leafN.MAP.DSI = brm(cover.change ~ leafN.final*mean.MAP + leafN.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                             family = gaussian(),
                                             prior = priors,
                                             data = imputed.NW.annual)

summary(imputed.NW.annual.traits.leafN.MAP.DSI)
# leafN x DSI

#saveRDS(imputed.NW.annual.traits.leafN.MAP.DSI, file = "./Results/imputed.NW.annual.traits.leafN.MAP.DSI.rds")
bayes_R2(imputed.NW.annual.traits.leafN.MAP.DSI)
# R2 0.1705448 0.05178855 0.07863608 0.2805161

conditional_effects(imputed.NW.annual.traits.leafN.MAP.DSI)
# higher leafN x less drought (higher DSI) or lower leafN x more drought (lower DSI)

imputed.NW.annual.traits.leafN.MAP.DSI = readRDS("./Results/imputed.NW.annual.traits.leafN.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.annual.traits.leafN.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.annual.traits.leafN.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.annual.traits.leafN.MAP.DSI) # p = 0.759
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.annual.traits.leafN.MAP.DSI)

imputed.NW.annual.traits.root.depth.MAP.DSI = brm(cover.change ~ root.depth.final*mean.MAP + root.depth.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                  family = gaussian(),
                                                  prior = priors,
                                                  data = imputed.NW.annual)

summary(imputed.NW.annual.traits.root.depth.MAP.DSI)
# MAP
#saveRDS(imputed.NW.annual.traits.root.depth.MAP.DSI, file = "./Results/imputed.NW.annual.traits.depth.MAP.DSI.rds")

imputed.NW.annual.traits.root.depth.MAP.DSI = readRDS("./Results/imputed.NW.annual.traits.depth.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.annual.traits.root.depth.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.annual.traits.root.depth.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.annual.traits.root.depth.MAP.DSI) # p = 0.939
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.annual.traits.root.depth.MAP.DSI)

imputed.NW.annual.traits.RTD.MAP.DSI = brm(cover.change ~ RTD.final*mean.MAP + RTD.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                           family = gaussian(),
                                           prior = priors,
                                           data = imputed.NW.annual)

summary(imputed.NW.annual.traits.RTD.MAP.DSI)
#saveRDS(imputed.NW.annual.traits.RTD.MAP.DSI, file = "./Results/imputed.NW.annual.traits.RTD.MAP.DSI.rds")

imputed.NW.annual.traits.RTD.MAP.DSI = readRDS("./Results/imputed.NW.annual.traits.RTD.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.annual.traits.RTD.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.annual.traits.RTD.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.annual.traits.RTD.MAP.DSI) # p = 0.824
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.annual.traits.RTD.MAP.DSI)

imputed.NW.annual.traits.SRL.MAP.DSI = brm(cover.change ~ SRL.final*mean.MAP + SRL.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                           family = gaussian(),
                                           prior = priors,
                                           data = imputed.NW.annual)

summary(imputed.NW.annual.traits.SRL.MAP.DSI)
# MAP
#saveRDS(imputed.NW.annual.traits.SRL.MAP.DSI, file = "./Results/imputed.NW.annual.traits.SRL.MAP.DSI.rds")

imputed.NW.annual.traits.SRL.MAP.DSI = readRDS("./Results/imputed.NW.annual.traits.SRL.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.annual.traits.SRL.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.annual.traits.SRL.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.annual.traits.SRL.MAP.DSI) # p = 0.915
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.annual.traits.SRL.MAP.DSI)

imputed.NW.annual.traits.RMF.MAP.DSI = brm(cover.change ~ RMF.final*mean.MAP + RMF.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                           family = gaussian(),
                                           prior = priors,
                                           data = imputed.NW.annual)

summary(imputed.NW.annual.traits.RMF.MAP.DSI)
# MAP
#saveRDS(imputed.NW.annual.traits.RMF.MAP.DSI, file = "./Results/imputed.NW.annual.traits.RMF.MAP.DSI.rds")

imputed.NW.annual.traits.RMF.MAP.DSI = readRDS("./Results/imputed.NW.annual.traits.RMF.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.annual.traits.RMF.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.annual.traits.RMF.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.annual.traits.RMF.MAP.DSI) # p = 0.876
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.annual.traits.RMF.MAP.DSI)

#### imputed PERENNIAL traits lifespan model NW ####

priors <- c(prior(normal(0, 10), class = b))

perennial.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                  SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                family = gaussian(),
                                prior = priors,
                                data = imputed.NW.perennial)

summary(perennial.traits.NW.model)
# nothing significant

#saveRDS(perennial.traits.NW.model, file = "./Results/perennial.imputed.traits.no_woody.rds")
bayes_R2(perennial.traits.NW.model)
# R2 0.081794 0.03372474 0.03279104 0.1602133

perennial.traits.NW.model = readRDS("./Results/perennial.imputed.traits.no_woody.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(perennial.traits.NW.model, summary = TRUE))
fitted_vals <- as.data.frame(fitted(perennial.traits.NW.model, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(perennial.traits.NW.model) # p = 0.790
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(perennial.traits.NW.model)

imputed.NW.perennial.traits.height.leafN = brm(cover.change ~ height.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                               family = gaussian(),
                                               prior = priors,
                                               data = imputed.NW.perennial)

summary(imputed.NW.perennial.traits.height.leafN)
# height x leafN significant positive

#saveRDS(imputed.NW.perennial.traits.height.leafN, file = "./Results/perennial.imputed.traits.NW.height.leafN.rds")
bayes_R2(imputed.NW.perennial.traits.height.leafN)
#R2 0.06295824 0.03318674 0.01677983 0.1430209

conditional_effects(imputed.NW.perennial.traits.height.leafN)
# higher cover for taller plants with higher leafN

imputed.NW.perennial.traits.height.leafN = readRDS("./Results/perennial.imputed.traits.NW.height.leafN.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.traits.height.leafN, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.traits.height.leafN, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.traits.height.leafN) # p = 0.726
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.traits.height.leafN)

imputed.NW.perennial.traits.depth.leafN = brm(cover.change ~ root.depth.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                              family = gaussian(),
                                              prior = priors,
                                              data = imputed.NW.perennial)

summary(imputed.NW.perennial.traits.depth.leafN)
# nothing significant
#saveRDS(imputed.NW.perennial.traits.depth.leafN, file = "./Results/perennial.imputed.traits.NW.depth.leafN.rds")

imputed.NW.perennial.traits.depth.leafN = readRDS("./Results/perennial.imputed.traits.NW.depth.leafN.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.traits.depth.leafN, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.traits.depth.leafN, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.traits.depth.leafN) # p = 0.711
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.traits.depth.leafN)

imputed.NW.perennial.traits.RTD.SRL = brm(cover.change ~ RTD.final*SRL.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                          family = gaussian(),
                                          prior = priors,
                                          data = imputed.NW.perennial)

summary(imputed.NW.perennial.traits.RTD.SRL)
# nothing significant
#saveRDS(imputed.NW.perennial.traits.RTD.SRL, file = "./Results/perennial.imputed.traits.NW.RTD.SRL.rds")

imputed.NW.perennial.traits.RTD.SRL = readRDS("./Results/perennial.imputed.traits.NW.RTD.SRL.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.traits.RTD.SRL, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.traits.RTD.SRL, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.traits.RTD.SRL) # p = 0.735
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.traits.RTD.SRL)

imputed.NW.perennial.traits.leafN.RMF = brm(cover.change ~ leafN.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                            family = gaussian(),
                                            prior = priors,
                                            data = imputed.NW.perennial)

summary(imputed.NW.perennial.traits.leafN.RMF)
# nothing significant
#saveRDS(imputed.NW.perennial.traits.leafN.RMF, file = "./Results/perennial.imputed.traits.NW.leafN.RMF.rds")

imputed.NW.perennial.traits.leafN.RMF = readRDS("./Results/perennial.imputed.traits.NW.leafN.RMF.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.traits.leafN.RMF, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.traits.leafN.RMF, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.traits.leafN.RMF) # p = 0.716
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.traits.leafN.RMF)

#### Perennial Environment ####

priors <- c(prior(normal(0, 10), class = b))

imputed.NW.perennial.traits.height.MAP.DSI = brm(cover.change ~ height.final*mean.MAP + height.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                 family = gaussian(),
                                                 prior = priors,
                                                 data = imputed.NW.perennial)

summary(imputed.NW.perennial.traits.height.MAP.DSI)
#saveRDS(imputed.NW.perennial.traits.height.MAP.DSI, file = "./Results/imputed.NW.perennial.traits.height.MAP.DSI.rds")

imputed.NW.perennial.traits.height.MAP.DSI = readRDS("./Results/imputed.NW.perennial.traits.height.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.traits.height.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.traits.height.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.traits.height.MAP.DSI) # p = 0.809
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.traits.height.MAP.DSI)

imputed.NW.perennial.traits.leafN.MAP.DSI = brm(cover.change ~ leafN.final*mean.MAP + leafN.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                family = gaussian(),
                                                prior = priors,
                                                data = imputed.NW.perennial)

summary(imputed.NW.perennial.traits.leafN.MAP.DSI)
#saveRDS(imputed.NW.perennial.traits.leafN.MAP.DSI, file = "./Results/imputed.NW.perennial.traits.leafN.MAP.DSI.rds")

imputed.NW.perennial.traits.leafN.MAP.DSI = readRDS("./Results/imputed.NW.perennial.traits.leafN.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.traits.leafN.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.traits.leafN.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.traits.leafN.MAP.DSI) # p = 0.719
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.traits.leafN.MAP.DSI)

imputed.NW.perennial.traits.root.depth.MAP.DSI = brm(cover.change ~ root.depth.final*mean.MAP + root.depth.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                     family = gaussian(),
                                                     prior = priors,
                                                     data = imputed.NW.perennial)

summary(imputed.NW.perennial.traits.root.depth.MAP.DSI)
#saveRDS(imputed.NW.perennial.traits.root.depth.MAP.DSI, file = "./Results/imputed.NW.perennial.traits.depth.MAP.DSI.rds")

imputed.NW.perennial.traits.root.depth.MAP.DSI = readRDS("./Results/imputed.NW.perennial.traits.depth.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.traits.root.depth.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.traits.root.depth.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.traits.root.depth.MAP.DSI) # p = 0.775
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.traits.root.depth.MAP.DSI)

imputed.NW.perennial.traits.RTD.MAP.DSI = brm(cover.change ~ RTD.final*mean.MAP + RTD.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                              family = gaussian(),
                                              prior = priors,
                                              data = imputed.NW.perennial)

summary(imputed.NW.perennial.traits.RTD.MAP.DSI)
#saveRDS(imputed.NW.perennial.traits.RTD.MAP.DSI, file = "./Results/imputed.NW.perennial.traits.RTD.MAP.DSI.rds")

imputed.NW.perennial.traits.RTD.MAP.DSI = readRDS("./Results/imputed.NW.perennial.traits.RTD.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.traits.RTD.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.traits.RTD.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.traits.RTD.MAP.DSI) # p = 0.778
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.traits.RTD.MAP.DSI)

imputed.NW.perennial.traits.SRL.MAP.DSI = brm(cover.change ~ SRL.final*mean.MAP + SRL.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                              family = gaussian(),
                                              prior = priors,
                                              data = imputed.NW.perennial)

summary(imputed.NW.perennial.traits.SRL.MAP.DSI)
#saveRDS(imputed.NW.perennial.traits.SRL.MAP.DSI, file = "./Results/imputed.NW.perennial.traits.SRL.MAP.DSI.rds")

imputed.NW.perennial.traits.SRL.MAP.DSI = readRDS("./Results/imputed.NW.perennial.traits.SRL.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.traits.SRL.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.traits.SRL.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.traits.SRL.MAP.DSI) # p = 0.857
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.traits.SRL.MAP.DSI)

imputed.NW.perennial.traits.RMF.MAP.DSI = brm(cover.change ~ RMF.final*mean.MAP + RMF.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                              family = gaussian(),
                                              prior = priors,
                                              data = imputed.NW.perennial)

summary(imputed.NW.perennial.traits.RMF.MAP.DSI)
#saveRDS(imputed.NW.perennial.traits.RMF.MAP.DSI, file = "./Results/imputed.NW.perennial.traits.RMF.MAP.DSI.rds")

imputed.NW.perennial.traits.RMF.MAP.DSI = readRDS("./Results/imputed.NW.perennial.traits.RMF.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.traits.RMF.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.traits.RMF.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.traits.RMF.MAP.DSI) # p = 0.759
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.traits.RMF.MAP.DSI)

#### impute FORB and LEGUME traits functional group NW ####

priors <- c(prior(normal(0, 10), class = b))

forb.legume.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                             SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                           family = gaussian(),
                           prior = priors,
                           data = imputed.NW.forb.legume)

summary(forb.legume.traits.NW.model)
# leafN significant
# height significant
# RTD significant

#saveRDS(forb.legume.traits.NW.model, file = "./Results/forb.legume.imputed.traits.no_woody.rds")
bayes_R2(forb.legume.traits.NW.model)
# R2 0.1206621 0.04300289 0.05401455 0.2241898

forb.legume.traits.NW.model = readRDS("./Results/forb.legume.imputed.traits.no_woody.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(forb.legume.traits.NW.model, summary = TRUE))
fitted_vals <- as.data.frame(fitted(forb.legume.traits.NW.model, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(forb.legume.traits.NW.model) # p = 0.508
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(forb.legume.traits.NW.model)

imputed.NW.forb.legume.traits.height.leafN = brm(cover.change ~ height.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                          family = gaussian(),
                                          prior = priors,
                                          data = imputed.NW.forb.legume)

summary(imputed.NW.forb.legume.traits.height.leafN)
# leafN
#saveRDS(imputed.NW.forb.legume.traits.height.leafN, file = "./Results/forb.legume.imputed.traits.NW.height.leafN.rds")

imputed.NW.forb.legume.traits.height.leafN = readRDS("./Results/forb.legume.imputed.traits.NW.height.leafN.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.forb.legume.traits.height.leafN, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.forb.legume.traits.height.leafN, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.forb.legume.traits.height.leafN) # p = 0.452
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.forb.legume.traits.height.leafN)

imputed.NW.forb.legume.traits.depth.leafN = brm(cover.change ~ root.depth.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                         family = gaussian(),
                                         prior = priors,
                                         data = imputed.NW.forb.legume)

summary(imputed.NW.forb.legume.traits.depth.leafN)
# leafN
#saveRDS(imputed.NW.forb.legume.traits.depth.leafN, file = "./Results/forb.legume.imputed.traits.NW.depth.leafN.rds")

imputed.NW.forb.legume.traits.depth.leafN = readRDS("./Results/forb.legume.imputed.traits.NW.depth.leafN.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.forb.legume.traits.depth.leafN, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.forb.legume.traits.depth.leafN, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.forb.legume.traits.depth.leafN) # p = 0.397
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.forb.legume.traits.depth.leafN)

imputed.NW.forb.legume.traits.RTD.SRL = brm(cover.change ~ RTD.final*SRL.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                     family = gaussian(),
                                     prior = priors,
                                     data = imputed.NW.forb.legume)

summary(imputed.NW.forb.legume.traits.RTD.SRL)
# RTD x SRL significant, negative

saveRDS(imputed.NW.forb.legume.traits.RTD.SRL, file = "./Results/forb.legume.imputed.traits.NW.RTD.SRL.rds")
bayes_R2(imputed.NW.forb.legume.traits.RTD.SRL)
# R2 0.08202759 0.03926922 0.02649765 0.1789433

conditional_effects(imputed.NW.forb.legume.traits.RTD.SRL)
# higher cover with higher RTD and lower SRL or lower RTD with higher SRL

imputed.NW.forb.legume.traits.RTD.SRL = readRDS("./Results/forb.legume.imputed.traits.NW.RTD.SRL.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.forb.legume.traits.RTD.SRL, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.forb.legume.traits.RTD.SRL, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.forb.legume.traits.RTD.SRL) # p = 0.532
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.forb.legume.traits.RTD.SRL)

imputed.NW.forb.legume.traits.leafN.RMF = brm(cover.change ~ leafN.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                       family = gaussian(),
                                       prior = priors,
                                       data = imputed.NW.forb.legume)

summary(imputed.NW.forb.legume.traits.leafN.RMF)
# leafN significant positive
#saveRDS(imputed.NW.forb.legume.traits.leafN.RMF, file = "./Results/forb.legume.imputed.traits.NW.leafN.RMF.rds")

imputed.NW.forb.legume.traits.leafN.RMF = readRDS("./Results/forb.legume.imputed.traits.NW.leafN.RMF.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.forb.legume.traits.leafN.RMF, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.forb.legume.traits.leafN.RMF, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.forb.legume.traits.leafN.RMF) # p = 0.391
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.forb.legume.traits.leafN.RMF)

#### Forb and Legume Environment ####

priors <- c(prior(normal(0, 10), class = b))

imputed.NW.forb.legume.traits.height.MAP.DSI = brm(cover.change ~ height.final*mean.MAP + height.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                   family = gaussian(),
                                                   prior = priors,
                                                   data = imputed.NW.forb.legume)

summary(imputed.NW.forb.legume.traits.height.MAP.DSI)
# height
#saveRDS(imputed.NW.forb.legume.traits.height.MAP.DSI, file = "./Results/imputed.NW.forb.legume.traits.height.MAP.DSI.rds")

imputed.NW.forb.legume.traits.height.MAP.DSI = readRDS("./Results/imputed.NW.forb.legume.traits.height.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.forb.legume.traits.height.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.forb.legume.traits.height.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.forb.legume.traits.height.MAP.DSI) # p = 0.559
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.forb.legume.traits.height.MAP.DSI)

imputed.NW.forb.legume.traits.leafN.MAP.DSI = brm(cover.change ~ leafN.final*mean.MAP + leafN.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                  family = gaussian(),
                                                  prior = priors,
                                                  data = imputed.NW.forb.legume)

summary(imputed.NW.forb.legume.traits.leafN.MAP.DSI)
# leafN 
# leafN x DSI

saveRDS(imputed.NW.forb.legume.traits.leafN.MAP.DSI, file = "./Results/imputed.NW.forb.legume.traits.leafN.MAP.DSI.rds")
bayes_R2(imputed.NW.forb.legume.traits.leafN.MAP.DSI)
# R2 0.1104155 0.04644593 0.04187595 0.2197352

conditional_effects(imputed.NW.forb.legume.traits.leafN.MAP.DSI)
# higher leafN x less drought (higher DSI) or lower leafN x more drought (lower DSI)

imputed.NW.forb.legume.traits.leafN.MAP.DSI = readRDS("./Results/imputed.NW.forb.legume.traits.leafN.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.forb.legume.traits.leafN.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.forb.legume.traits.leafN.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.forb.legume.traits.leafN.MAP.DSI) # p = 0.470
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.forb.legume.traits.leafN.MAP.DSI)

imputed.NW.forb.legume.traits.root.depth.MAP.DSI = brm(cover.change ~ root.depth.final*mean.MAP + root.depth.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                       family = gaussian(),
                                                       prior = priors,
                                                       data = imputed.NW.forb.legume)

summary(imputed.NW.forb.legume.traits.root.depth.MAP.DSI)
#saveRDS(imputed.NW.forb.legume.traits.root.depth.MAP.DSI, file = "./Results/imputed.NW.forb.legume.traits.depth.MAP.DSI.rds")

imputed.NW.forb.legume.traits.root.depth.MAP.DSI = readRDS("./Results/imputed.NW.forb.legume.traits.depth.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.forb.legume.traits.root.depth.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.forb.legume.traits.root.depth.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.forb.legume.traits.root.depth.MAP.DSI) # p = 0.520
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.forb.legume.traits.root.depth.MAP.DSI)

imputed.NW.forb.legume.traits.RTD.MAP.DSI = brm(cover.change ~ RTD.final*mean.MAP + RTD.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                family = gaussian(),
                                                prior = priors,
                                                data = imputed.NW.forb.legume)

summary(imputed.NW.forb.legume.traits.RTD.MAP.DSI)
#saveRDS(imputed.NW.forb.legume.traits.RTD.MAP.DSI, file = "./Results/imputed.NW.forb.legume.traits.RTD.MAP.DSI.rds")

imputed.NW.forb.legume.traits.RTD.MAP.DSI = readRDS("./Results/imputed.NW.forb.legume.traits.RTD.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.forb.legume.traits.RTD.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.forb.legume.traits.RTD.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.forb.legume.traits.RTD.MAP.DSI) # p = 0.538
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.forb.legume.traits.RTD.MAP.DSI)

imputed.NW.forb.legume.traits.SRL.MAP.DSI = brm(cover.change ~ SRL.final*mean.MAP + SRL.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                family = gaussian(),
                                                prior = priors,
                                                data = imputed.NW.forb.legume)

summary(imputed.NW.forb.legume.traits.SRL.MAP.DSI)

saveRDS(imputed.NW.forb.legume.traits.SRL.MAP.DSI, file = "./Results/imputed.NW.forb.legume.traits.SRL.MAP.DSI.rds")

imputed.NW.forb.legume.traits.SRL.MAP.DSI = readRDS("./Results/imputed.NW.forb.legume.traits.SRL.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.forb.legume.traits.SRL.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.forb.legume.traits.SRL.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.forb.legume.traits.SRL.MAP.DSI) # p = 0.515
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.forb.legume.traits.SRL.MAP.DSI)

imputed.NW.forb.legume.traits.RMF.MAP.DSI = brm(cover.change ~ RMF.final*mean.MAP + RMF.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                family = gaussian(),
                                                prior = priors,
                                                data = imputed.NW.forb.legume)

summary(imputed.NW.forb.legume.traits.RMF.MAP.DSI)
#saveRDS(imputed.NW.forb.legume.traits.RMF.MAP.DSI, file = "./Results/imputed.NW.forb.legume.traits.RMF.MAP.DSI.rds")

imputed.NW.forb.legume.traits.RMF.MAP.DSI = readRDS("./Results/imputed.NW.forb.legume.traits.RMF.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.forb.legume.traits.RMF.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.forb.legume.traits.RMF.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.forb.legume.traits.RMF.MAP.DSI) # p = 0.493
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.forb.legume.traits.RMF.MAP.DSI)

#### impute GRAMINOID traits functional group NW ####

priors <- c(prior(normal(0, 10), class = b))

graminoid.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                  SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                family = gaussian(),
                                prior = priors,
                                data = imputed.NW.graminoid)

summary(graminoid.traits.NW.model)
# nothing significant

#saveRDS(graminoid.traits.NW.model, file = "./Results/graminoid.imputed.traits.no_woody.rds")
bayes_R2(graminoid.traits.NW.model)
# R2 0.1064767 0.03930761 0.04390368 0.1977966

graminoid.traits.NW.model = readRDS("./Results/graminoid.imputed.traits.no_woody.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(graminoid.traits.NW.model, summary = TRUE))
fitted_vals <- as.data.frame(fitted(graminoid.traits.NW.model, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(graminoid.traits.NW.model) # p = 0.919
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(graminoid.traits.NW.model)

imputed.NW.graminoid.traits.height.leafN = brm(cover.change ~ height.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                               family = gaussian(),
                                               prior = priors,
                                               data = imputed.NW.graminoid)

summary(imputed.NW.graminoid.traits.height.leafN)
# nothing significant
saveRDS(imputed.NW.graminoid.traits.height.leafN, file = "./Results/graminoid.imputed.traits.NW.height.leafN.rds")

imputed.NW.graminoid.traits.height.leafN = readRDS("./Results/graminoid.imputed.traits.NW.height.leafN.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.graminoid.traits.height.leafN, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.graminoid.traits.height.leafN, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.graminoid.traits.height.leafN) # p = 0.890
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.graminoid.traits.height.leafN)

imputed.NW.graminoid.traits.depth.leafN = brm(cover.change ~ root.depth.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                              family = gaussian(),
                                              prior = priors,
                                              data = imputed.NW.graminoid)

summary(imputed.NW.graminoid.traits.depth.leafN)
# MAP is positive, significant
saveRDS(imputed.NW.graminoid.traits.depth.leafN, file = "./Results/graminoid.imputed.traits.NW.depth.leafN.rds")

imputed.NW.graminoid.traits.depth.leafN = readRDS("./Results/graminoid.imputed.traits.NW.depth.leafN.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.graminoid.traits.depth.leafN, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.graminoid.traits.depth.leafN, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.graminoid.traits.depth.leafN) # p = 0.917
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.graminoid.traits.depth.leafN)

imputed.NW.graminoid.traits.RTD.SRL = brm(cover.change ~ RTD.final*SRL.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                          family = gaussian(),
                                          prior = priors,
                                          data = imputed.NW.graminoid)

summary(imputed.NW.graminoid.traits.RTD.SRL)
# nothing significant
#saveRDS(imputed.NW.graminoid.traits.RTD.SRL, file = "./Results/graminoid.imputed.traits.NW.RTD.SRL.rds")

imputed.NW.graminoid.traits.RTD.SRL = readRDS("./Results/graminoid.imputed.traits.NW.RTD.SRL.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.graminoid.traits.RTD.SRL, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.graminoid.traits.RTD.SRL, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.graminoid.traits.RTD.SRL) # p = 0.922
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.graminoid.traits.RTD.SRL)

imputed.NW.graminoid.traits.leafN.RMF = brm(cover.change ~ leafN.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                            family = gaussian(),
                                            prior = priors,
                                            data = imputed.NW.graminoid)

summary(imputed.NW.graminoid.traits.leafN.RMF)
# nothing significant
saveRDS(imputed.NW.graminoid.traits.leafN.RMF, file = "./Results/graminoid.imputed.traits.NW.leafN.RMF.rds")

imputed.NW.graminoid.traits.leafN.RMF = readRDS("./Results/graminoid.imputed.traits.NW.leafN.RMF.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.graminoid.traits.leafN.RMF, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.graminoid.traits.leafN.RMF, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.graminoid.traits.leafN.RMF) # p = 0.891
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.graminoid.traits.leafN.RMF)

#### Graminoid Environment ####

priors <- c(prior(normal(0, 10), class = b))

imputed.NW.graminoid.traits.height.MAP.DSI = brm(cover.change ~ height.final*mean.MAP + height.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                 family = gaussian(),
                                                 prior = priors,
                                                 data = imputed.NW.graminoid)

summary(imputed.NW.graminoid.traits.height.MAP.DSI)
#saveRDS(imputed.NW.graminoid.traits.height.MAP.DSI, file = "./Results/imputed.NW.graminoid.traits.height.MAP.DSI.rds")

imputed.NW.graminoid.traits.height.MAP.DSI = readRDS("./Results/imputed.NW.graminoid.traits.height.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.graminoid.traits.height.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.graminoid.traits.height.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.graminoid.traits.height.MAP.DSI) # p = 0.901
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.graminoid.traits.height.MAP.DSI)

imputed.NW.graminoid.traits.leafN.MAP.DSI = brm(cover.change ~ leafN.final*mean.MAP + leafN.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                family = gaussian(),
                                                prior = priors,
                                                data = imputed.NW.graminoid)

summary(imputed.NW.graminoid.traits.leafN.MAP.DSI)
saveRDS(imputed.NW.graminoid.traits.leafN.MAP.DSI, file = "./Results/imputed.NW.graminoid.traits.leafN.MAP.DSI.rds")

imputed.NW.graminoid.traits.leafN.MAP.DSI = readRDS("./Results/imputed.NW.graminoid.traits.leafN.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.graminoid.traits.leafN.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.graminoid.traits.leafN.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.graminoid.traits.leafN.MAP.DSI) # p = 0.948
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.graminoid.traits.leafN.MAP.DSI)

imputed.NW.graminoid.traits.root.depth.MAP.DSI = brm(cover.change ~ root.depth.final*mean.MAP + root.depth.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                     family = gaussian(),
                                                     prior = priors,
                                                     data = imputed.NW.graminoid)

summary(imputed.NW.graminoid.traits.root.depth.MAP.DSI)
#saveRDS(imputed.NW.graminoid.traits.root.depth.MAP.DSI, file = "./Results/imputed.NW.graminoid.traits.depth.MAP.DSI.rds")

imputed.NW.graminoid.traits.root.depth.MAP.DSI = readRDS("./Results/imputed.NW.graminoid.traits.depth.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.graminoid.traits.root.depth.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.graminoid.traits.root.depth.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.graminoid.traits.root.depth.MAP.DSI) # p = 0.901
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.graminoid.traits.root.depth.MAP.DSI)

imputed.NW.graminoid.traits.RTD.MAP.DSI = brm(cover.change ~ RTD.final*mean.MAP + RTD.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                              family = gaussian(),
                                              prior = priors,
                                              data = imputed.NW.graminoid)

summary(imputed.NW.graminoid.traits.RTD.MAP.DSI)
#saveRDS(imputed.NW.graminoid.traits.RTD.MAP.DSI, file = "./Results/imputed.NW.graminoid.traits.RTD.MAP.DSI.rds")

imputed.NW.graminoid.traits.RTD.MAP.DSI = readRDS("./Results/imputed.NW.graminoid.traits.RTD.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.graminoid.traits.RTD.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.graminoid.traits.RTD.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.graminoid.traits.RTD.MAP.DSI) # p = 0.924
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.graminoid.traits.RTD.MAP.DSI)

imputed.NW.graminoid.traits.SRL.MAP.DSI = brm(cover.change ~ SRL.final*mean.MAP + SRL.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                              family = gaussian(),
                                              prior = priors,
                                              data = imputed.NW.graminoid)

summary(imputed.NW.graminoid.traits.SRL.MAP.DSI)
#saveRDS(imputed.NW.graminoid.traits.SRL.MAP.DSI, file = "./Results/imputed.NW.graminoid.traits.SRL.MAP.DSI.rds")

imputed.NW.graminoid.traits.SRL.MAP.DSI = readRDS("./Results/imputed.NW.graminoid.traits.SRL.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.graminoid.traits.SRL.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.graminoid.traits.SRL.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.graminoid.traits.SRL.MAP.DSI) # p = 0.929
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.graminoid.traits.SRL.MAP.DSI)

imputed.NW.graminoid.traits.RMF.MAP.DSI = brm(cover.change ~ RMF.final*mean.MAP + RMF.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                              family = gaussian(),
                                              prior = priors,
                                              data = imputed.NW.graminoid)

summary(imputed.NW.graminoid.traits.RMF.MAP.DSI)
#saveRDS(imputed.NW.graminoid.traits.RMF.MAP.DSI, file = "./Results/imputed.NW.graminoid.traits.RMF.MAP.DSI.rds")

imputed.NW.graminoid.traits.RMF.MAP.DSI = readRDS("./Results/imputed.NW.graminoid.traits.RMF.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.graminoid.traits.RMF.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.graminoid.traits.RMF.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.graminoid.traits.RMF.MAP.DSI) # p = 0.921
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.graminoid.traits.RMF.MAP.DSI)

#### imputed traits ANNUAL FORB LEGUME NW ####

priors <- c(prior(normal(0, 10), class = b))

annual.forb.legume.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                           SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                         family = gaussian(),
                                         prior = priors,
                                         data = imputed.NW.annual.forb.legume)

summary(annual.forb.legume.traits.NW.model)
# nothing significant

#saveRDS(annual.forb.legume.traits.NW.model, file = "./Results/annual.forb.legume.imputed.traits.no_woody.rds")
bayes_R2(annual.forb.legume.traits.NW.model)
# R2 0.2273124 0.06464005 0.1131978 0.3652143

annual.forb.legume.traits.NW.model = readRDS("./Results/annual.forb.legume.imputed.traits.no_woody.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(annual.forb.legume.traits.NW.model, summary = TRUE))
fitted_vals <- as.data.frame(fitted(annual.forb.legume.traits.NW.model, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(annual.forb.legume.traits.NW.model) # p = 0.758
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(annual.forb.legume.traits.NW.model)

imputed.NW.annual.forb.legume.traits.height.leafN = brm(cover.change ~ height.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                               family = gaussian(),
                                                               prior = priors,
                                                               data = imputed.NW.annual.forb.legume)

summary(imputed.NW.annual.forb.legume.traits.height.leafN)
# nothing significant
#saveRDS(imputed.NW.annual.forb.legume.traits.height.leafN, file = "./Results/annual.forb.legume.imputed.traits.NW.height.leafN.rds")

imputed.NW.annual.forb.legume.traits.height.leafN = readRDS("./Results/annual.forb.legume.imputed.traits.NW.height.leafN.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.annual.forb.legume.traits.height.leafN, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.annual.forb.legume.traits.height.leafN, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.annual.forb.legume.traits.height.leafN) # p = 0.652
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.annual.forb.legume.traits.height.leafN)

imputed.NW.annual.forb.legume.traits.depth.leafN = brm(cover.change ~ root.depth.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                              family = gaussian(),
                                                              prior = priors,
                                                              data = imputed.NW.annual.forb.legume)

summary(imputed.NW.annual.forb.legume.traits.depth.leafN)
# MAP is positive, significant
#saveRDS(imputed.NW.annual.forb.legume.traits.depth.leafN, file = "./Results/annual.forb.legume.imputed.traits.NW.depth.leafN.rds")

imputed.NW.annual.forb.legume.traits.depth.leafN = readRDS("./Results/annual.forb.legume.imputed.traits.NW.depth.leafN.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.annual.forb.legume.traits.depth.leafN, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.annual.forb.legume.traits.depth.leafN, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.annual.forb.legume.traits.depth.leafN) # p = 0.605
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.annual.forb.legume.traits.depth.leafN)

imputed.NW.annual.forb.legume.traits.RTD.SRL = brm(cover.change ~ RTD.final*SRL.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                          family = gaussian(),
                                                          prior = priors,
                                                          data = imputed.NW.annual.forb.legume)

summary(imputed.NW.annual.forb.legume.traits.RTD.SRL)
# RTD x SRL
#saveRDS(imputed.NW.annual.forb.legume.traits.RTD.SRL, file = "./Results/annual.forb.legume.imputed.traits.NW.RTD.SRL.rds")
bayes_R2(imputed.NW.annual.forb.legume.traits.RTD.SRL)
# R2 0.1818185 0.06325504 0.0745352 0.3182902

conditional_effects(imputed.NW.annual.forb.legume.traits.RTD.SRL)
# higher cover with higher RTD and lower SRL or lower RTD with higher SRL

imputed.NW.annual.forb.legume.traits.RTD.SRL = readRDS("./Results/annual.forb.legume.imputed.traits.NW.RTD.SRL.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.annual.forb.legume.traits.RTD.SRL, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.annual.forb.legume.traits.RTD.SRL, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.annual.forb.legume.traits.RTD.SRL) # p = 0.769
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.annual.forb.legume.traits.RTD.SRL)

imputed.NW.annual.forb.legume.traits.leafN.RMF = brm(cover.change ~ leafN.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                            family = gaussian(),
                                                            prior = priors,
                                                            data = imputed.NW.annual.forb.legume)

summary(imputed.NW.annual.forb.legume.traits.leafN.RMF)
# MAP significant positive
#saveRDS(imputed.NW.annual.forb.legume.traits.leafN.RMF, file = "./Results/annual.forb.legume.imputed.traits.NW.leafN.RMF.rds")

imputed.NW.annual.forb.legume.traits.leafN.RMF = readRDS("./Results/annual.forb.legume.imputed.traits.NW.leafN.RMF.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.annual.forb.legume.traits.leafN.RMF, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.annual.forb.legume.traits.leafN.RMF, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.annual.forb.legume.traits.leafN.RMF) # p = 0.0.580
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.annual.forb.legume.traits.leafN.RMF)

#### Annual forb and legume Environment ####

priors <- c(prior(normal(0, 10), class = b))

imputed.NW.annual.forb.legume.traits.height.MAP.DSI = brm(cover.change ~ height.final*mean.MAP + height.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                          family = gaussian(),
                                                          prior = priors,
                                                          data = imputed.NW.annual.forb.legume)

summary(imputed.NW.annual.forb.legume.traits.height.MAP.DSI)
#saveRDS(imputed.NW.annual.forb.legume.traits.height.MAP.DSI, file = "./Results/imputed.NW.annual.forb.legume.traits.height.MAP.DSI.rds")

imputed.NW.annual.forb.legume.traits.height.MAP.DSI = readRDS("./Results/imputed.NW.annual.forb.legume.traits.height.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.annual.forb.legume.traits.height.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.annual.forb.legume.traits.height.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.annual.forb.legume.traits.height.MAP.DSI) # p = 0.725
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.annual.forb.legume.traits.height.MAP.DSI)

imputed.NW.annual.forb.legume.traits.leafN.MAP.DSI = brm(cover.change ~ leafN.final*mean.MAP + leafN.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                         family = gaussian(),
                                                         prior = priors,
                                                         data = imputed.NW.annual.forb.legume)

summary(imputed.NW.annual.forb.legume.traits.leafN.MAP.DSI)
# leafN x DSI
#saveRDS(imputed.NW.annual.forb.legume.traits.leafN.MAP.DSI, file = "./Results/imputed.NW.annual.forb.legume.traits.leafN.MAP.DSI.rds")
bayes_R2(imputed.NW.annual.forb.legume.traits.leafN.MAP.DSI)
# R2 0.2487581 0.07176966 0.1201865 0.400603

conditional_effects(imputed.NW.annual.forb.legume.traits.leafN.MAP.DSI)
# higher leafN x less drought (higher DSI) or lower leafN x more drought (lower DSI)

imputed.NW.annual.forb.legume.traits.leafN.MAP.DSI = readRDS("./Results/imputed.NW.annual.forb.legume.traits.leafN.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.annual.forb.legume.traits.leafN.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.annual.forb.legume.traits.leafN.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.annual.forb.legume.traits.leafN.MAP.DSI) # p = 0.567
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.annual.forb.legume.traits.leafN.MAP.DSI)

imputed.NW.annual.forb.legume.traits.root.depth.MAP.DSI = brm(cover.change ~ root.depth.final*mean.MAP + root.depth.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                              family = gaussian(),
                                                              prior = priors,
                                                              data = imputed.NW.annual.forb.legume)

summary(imputed.NW.annual.forb.legume.traits.root.depth.MAP.DSI)
# MAP
#saveRDS(imputed.NW.annual.forb.legume.traits.root.depth.MAP.DSI, file = "./Results/imputed.NW.annual.forb.legume.traits.depth.MAP.DSI.rds")

imputed.NW.annual.forb.legume.traits.root.depth.MAP.DSI = readRDS("./Results/imputed.NW.annual.forb.legume.traits.depth.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.annual.forb.legume.traits.root.depth.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.annual.forb.legume.traits.root.depth.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.annual.forb.legume.traits.root.depth.MAP.DSI) # p = 0.667
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.annual.forb.legume.traits.root.depth.MAP.DSI)

imputed.NW.annual.forb.legume.traits.RTD.MAP.DSI = brm(cover.change ~ RTD.final*mean.MAP + RTD.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                       family = gaussian(),
                                                       prior = priors,
                                                       data = imputed.NW.annual.forb.legume)

summary(imputed.NW.annual.forb.legume.traits.RTD.MAP.DSI)
# RTD x MAP
#saveRDS(imputed.NW.annual.forb.legume.traits.RTD.MAP.DSI, file = "./Results/imputed.NW.annual.forb.legume.traits.RTD.MAP.DSI.rds")
bayes_R2(imputed.NW.annual.forb.legume.traits.RTD.MAP.DSI)
# R2 0.189431 0.0672549 0.07478275 0.3386802

conditional_effects(imputed.NW.annual.forb.legume.traits.RTD.MAP.DSI)
# higher leafN x less drought (higher DSI) or lower leafN x more drought (lower DSI)

imputed.NW.annual.forb.legume.traits.RTD.MAP.DSI = readRDS("./Results/imputed.NW.annual.forb.legume.traits.RTD.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.annual.forb.legume.traits.RTD.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.annual.forb.legume.traits.RTD.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.annual.forb.legume.traits.RTD.MAP.DSI) # p = 0.665
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.annual.forb.legume.traits.RTD.MAP.DSI)

imputed.NW.annual.forb.legume.traits.SRL.MAP.DSI = brm(cover.change ~ SRL.final*mean.MAP + SRL.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                       family = gaussian(),
                                                       prior = priors,
                                                       data = imputed.NW.annual.forb.legume)

summary(imputed.NW.annual.forb.legume.traits.SRL.MAP.DSI)
#saveRDS(imputed.NW.annual.forb.legume.traits.SRL.MAP.DSI, file = "./Results/imputed.NW.annual.forb.legume.traits.SRL.MAP.DSI.rds")

imputed.NW.annual.forb.legume.traits.SRL.MAP.DSI = readRDS("./Results/imputed.NW.annual.forb.legume.traits.SRL.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.annual.forb.legume.traits.SRL.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.annual.forb.legume.traits.SRL.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.annual.forb.legume.traits.SRL.MAP.DSI) # p = 0.678
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.annual.forb.legume.traits.SRL.MAP.DSI)

imputed.NW.annual.forb.legume.traits.RMF.MAP.DSI = brm(cover.change ~ RMF.final*mean.MAP + RMF.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                       family = gaussian(),
                                                       prior = priors,
                                                       data = imputed.NW.annual.forb.legume)

summary(imputed.NW.annual.forb.legume.traits.RMF.MAP.DSI)
# MAP
#saveRDS(imputed.NW.annual.forb.legume.traits.RMF.MAP.DSI, file = "./Results/imputed.NW.annual.forb.legume.traits.RMF.MAP.DSI.rds")

imputed.NW.annual.forb.legume.traits.RMF.MAP.DSI = readRDS("./Results/imputed.NW.annual.forb.legume.traits.RMF.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.annual.forb.legume.traits.RMF.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.annual.forb.legume.traits.RMF.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.annual.forb.legume.traits.RMF.MAP.DSI) # p = 0.631
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.annual.forb.legume.traits.RMF.MAP.DSI)

#### imputed traits PERENNIAL FORB LEGUME NW ####

priors <- c(prior(normal(0, 10), class = b))

perennial.forb.legume.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                              SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                            family = gaussian(),
                                            prior = priors,
                                            data = imputed.NW.perennial.forb.legume)

summary(perennial.forb.legume.traits.NW.model)
# RTD significant
#saveRDS(perennial.forb.legume.traits.NW.model, file = "./Results/perennial.forb.legume.imputed.traits.no_woody.rds")
bayes_R2(perennial.forb.legume.traits.NW.model)
#R2 0.1371829 0.05049297 0.05734497 0.2561353

perennial.forb.legume.traits.NW.model = readRDS("./Results/perennial.forb.legume.imputed.traits.no_woody.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(perennial.forb.legume.traits.NW.model, summary = TRUE))
fitted_vals <- as.data.frame(fitted(perennial.forb.legume.traits.NW.model, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(perennial.forb.legume.traits.NW.model) # p = 0.697
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(perennial.forb.legume.traits.NW.model)

imputed.NW.perennial.forb.legume.traits.height.leafN = brm(cover.change ~ height.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                                  family = gaussian(),
                                                                  prior = priors,
                                                                  data = imputed.NW.perennial.forb.legume)

summary(imputed.NW.perennial.forb.legume.traits.height.leafN)
# nothing significant
#saveRDS(imputed.NW.perennial.forb.legume.traits.height.leafN, file = "./Results/perennial.forb.legume.imputed.traits.NW.height.leafN.rds")

imputed.NW.perennial.forb.legume.traits.height.leafN = readRDS("./Results/perennial.forb.legume.imputed.traits.NW.height.leafN.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.forb.legume.traits.height.leafN, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.forb.legume.traits.height.leafN, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.forb.legume.traits.height.leafN) # p = 0.599
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.forb.legume.traits.height.leafN)

imputed.NW.perennial.forb.legume.traits.depth.leafN = brm(cover.change ~ root.depth.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                                 family = gaussian(),
                                                                 prior = priors,
                                                                 data = imputed.NW.perennial.forb.legume)

summary(imputed.NW.perennial.forb.legume.traits.depth.leafN)
#saveRDS(imputed.NW.perennial.forb.legume.traits.depth.leafN, file = "./Results/perennial.forb.legume.imputed.traits.NW.depth.leafN.rds")

imputed.NW.perennial.forb.legume.traits.depth.leafN = readRDS("./Results/perennial.forb.legume.imputed.traits.NW.depth.leafN.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.forb.legume.traits.depth.leafN, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.forb.legume.traits.depth.leafN, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.forb.legume.traits.depth.leafN) # p = 0.591
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.forb.legume.traits.depth.leafN)

imputed.NW.perennial.forb.legume.traits.RTD.SRL = brm(cover.change ~ RTD.final*SRL.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                             family = gaussian(),
                                                             prior = priors,
                                                             data = imputed.NW.perennial.forb.legume)

summary(imputed.NW.perennial.forb.legume.traits.RTD.SRL)
# nothing significant
#saveRDS(imputed.NW.perennial.forb.legume.traits.RTD.SRL, file = "./Results/perennial.forb.legume.imputed.traits.NW.RTD.SRL.rds")

imputed.NW.perennial.forb.legume.traits.RTD.SRL = readRDS("./Results/perennial.forb.legume.imputed.traits.NW.RTD.SRL.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.forb.legume.traits.RTD.SRL, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.forb.legume.traits.RTD.SRL, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.forb.legume.traits.RTD.SRL) # p = 0.547
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.forb.legume.traits.RTD.SRL)

imputed.NW.perennial.forb.legume.traits.leafN.RMF = brm(cover.change ~ leafN.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                               family = gaussian(),
                                                               prior = priors,
                                                               data = imputed.NW.perennial.forb.legume)

summary(imputed.NW.perennial.forb.legume.traits.leafN.RMF)
# nothing significant
#saveRDS(imputed.NW.perennial.forb.legume.traits.leafN.RMF, file = "./Results/perennial.forb.legume.imputed.traits.NW.leafN.RMF.rds")

imputed.NW.perennial.forb.legume.traits.leafN.RMF = readRDS("./Results/perennial.forb.legume.imputed.traits.NW.leafN.RMF.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.forb.legume.traits.leafN.RMF, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.forb.legume.traits.leafN.RMF, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.forb.legume.traits.leafN.RMF) # p = 0.576
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.forb.legume.traits.leafN.RMF)

#### Perennial forb and legume Environment ####

priors <- c(prior(normal(0, 10), class = b))

imputed.NW.perennial.forb.legume.traits.height.MAP.DSI = brm(cover.change ~ height.final*mean.MAP + height.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                             family = gaussian(),
                                                             prior = priors,
                                                             data = imputed.NW.perennial.forb.legume)

summary(imputed.NW.perennial.forb.legume.traits.height.MAP.DSI)
#saveRDS(imputed.NW.perennial.forb.legume.traits.height.MAP.DSI, file = "./Results/imputed.NW.perennial.forb.legume.traits.height.MAP.DSI.rds")

imputed.NW.perennial.forb.legume.traits.height.MAP.DSI = readRDS("./Results/imputed.NW.perennial.forb.legume.traits.height.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.forb.legume.traits.height.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.forb.legume.traits.height.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.forb.legume.traits.height.MAP.DSI) # p = 0.668
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.forb.legume.traits.height.MAP.DSI)

imputed.NW.perennial.forb.legume.traits.leafN.MAP.DSI = brm(cover.change ~ leafN.final*mean.MAP + leafN.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                            family = gaussian(),
                                                            prior = priors,
                                                            data = imputed.NW.perennial.forb.legume)

summary(imputed.NW.perennial.forb.legume.traits.leafN.MAP.DSI)
#saveRDS(imputed.NW.perennial.forb.legume.traits.leafN.MAP.DSI, file = "./Results/imputed.NW.perennial.forb.legume.traits.leafN.MAP.DSI.rds")

imputed.NW.perennial.forb.legume.traits.leafN.MAP.DSI = readRDS("./Results/imputed.NW.perennial.forb.legume.traits.leafN.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.forb.legume.traits.leafN.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.forb.legume.traits.leafN.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.forb.legume.traits.leafN.MAP.DSI) # p = 0.584
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.forb.legume.traits.leafN.MAP.DSI)

imputed.NW.perennial.forb.legume.traits.root.depth.MAP.DSI = brm(cover.change ~ root.depth.final*mean.MAP + root.depth.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                                 family = gaussian(),
                                                                 prior = priors,
                                                                 data = imputed.NW.perennial.forb.legume)

summary(imputed.NW.perennial.forb.legume.traits.root.depth.MAP.DSI)
#saveRDS(imputed.NW.perennial.forb.legume.traits.root.depth.MAP.DSI, file = "./Results/imputed.NW.perennial.forb.legume.traits.depth.MAP.DSI.rds")

imputed.NW.perennial.forb.legume.traits.root.depth.MAP.DSI = readRDS("./Results/imputed.NW.perennial.forb.legume.traits.depth.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.forb.legume.traits.root.depth.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.forb.legume.traits.root.depth.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.forb.legume.traits.root.depth.MAP.DSI) # p = 0.683
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.forb.legume.traits.root.depth.MAP.DSI)

imputed.NW.perennial.forb.legume.traits.RTD.MAP.DSI = brm(cover.change ~ RTD.final*mean.MAP + RTD.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                          family = gaussian(),
                                                          prior = priors,
                                                          data = imputed.NW.perennial.forb.legume)

summary(imputed.NW.perennial.forb.legume.traits.RTD.MAP.DSI)
#saveRDS(imputed.NW.perennial.forb.legume.traits.RTD.MAP.DSI, file = "./Results/imputed.NW.perennial.forb.legume.traits.RTD.MAP.DSI.rds")

imputed.NW.perennial.forb.legume.traits.RTD.MAP.DSI = readRDS("./Results/imputed.NW.perennial.forb.legume.traits.RTD.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.forb.legume.traits.RTD.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.forb.legume.traits.RTD.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.forb.legume.traits.RTD.MAP.DSI) # p = 0.686
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.forb.legume.traits.RTD.MAP.DSI)

imputed.NW.perennial.forb.legume.traits.SRL.MAP.DSI = brm(cover.change ~ SRL.final*mean.MAP + SRL.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                          family = gaussian(),
                                                          prior = priors,
                                                          data = imputed.NW.perennial.forb.legume)

summary(imputed.NW.perennial.forb.legume.traits.SRL.MAP.DSI)
#saveRDS(imputed.NW.perennial.forb.legume.traits.SRL.MAP.DSI, file = "./Results/imputed.NW.perennial.forb.legume.traits.SRL.MAP.DSI.rds")

imputed.NW.perennial.forb.legume.traits.SRL.MAP.DSI = readRDS("./Results/imputed.NW.perennial.forb.legume.traits.SRL.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.forb.legume.traits.SRL.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.forb.legume.traits.SRL.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.forb.legume.traits.SRL.MAP.DSI) # p = 0.675
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.forb.legume.traits.SRL.MAP.DSI)

imputed.NW.perennial.forb.legume.traits.RMF.MAP.DSI = brm(cover.change ~ RMF.final*mean.MAP + RMF.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                          family = gaussian(),
                                                          prior = priors,
                                                          data = imputed.NW.perennial.forb.legume)

summary(imputed.NW.perennial.forb.legume.traits.RMF.MAP.DSI)
#saveRDS(imputed.NW.perennial.forb.legume.traits.RMF.MAP.DSI, file = "./Results/imputed.NW.perennial.forb.legume.traits.RMF.MAP.DSI.rds")

imputed.NW.perennial.forb.legume.traits.RMF.MAP.DSI = readRDS("./Results/imputed.NW.perennial.forb.legume.traits.RMF.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.forb.legume.traits.RMF.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.forb.legume.traits.RMF.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.forb.legume.traits.RMF.MAP.DSI) # p = 0.622
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.forb.legume.traits.RMF.MAP.DSI)

#### imputed traits PERENNIAL GRAMINOID NW #####

priors <- c(prior(normal(0, 10), class = b))

perennial.graminoid.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                            SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                          family = gaussian(),
                                          prior = priors,
                                          data = imputed.NW.perennial.graminoid)

summary(perennial.graminoid.traits.NW.model)
# nothing significant

saveRDS(perennial.graminoid.traits.NW.model, file = "./Results/perennial.graminoid.imputed.traits.no_woody.rds")
bayes_R2(perennial.graminoid.traits.NW.model)
#R2 0.1316797 0.04890526 0.0545399 0.2443937

perennial.graminoid.traits.NW.model = readRDS("./Results/perennial.graminoid.imputed.traits.no_woody.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(perennial.graminoid.traits.NW.model, summary = TRUE))
fitted_vals <- as.data.frame(fitted(perennial.graminoid.traits.NW.model, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(perennial.graminoid.traits.NW.model) # p = 0.936
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(perennial.graminoid.traits.NW.model)

imputed.NW.perennial.graminoid.traits.height.leafN = brm(cover.change ~ height.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                         family = gaussian(),
                                                         prior = priors,
                                                         data = imputed.NW.perennial.graminoid)

summary(imputed.NW.perennial.graminoid.traits.height.leafN)
# nothing significant
#saveRDS(imputed.NW.perennial.graminoid.traits.height.leafN, file = "./Results/perennial.graminoid.imputed.traits.NW.height.leafN.rds")

imputed.NW.perennial.graminoid.traits.height.leafN = readRDS("./Results/perennial.graminoid.imputed.traits.NW.height.leafN.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.graminoid.traits.height.leafN, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.graminoid.traits.height.leafN, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.graminoid.traits.height.leafN) # p = 0.885
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.graminoid.traits.height.leafN)

imputed.NW.perennial.graminoid.traits.depth.leafN = brm(cover.change ~ root.depth.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                        family = gaussian(),
                                                        prior = priors,
                                                        data = imputed.NW.perennial.graminoid)

summary(imputed.NW.perennial.graminoid.traits.depth.leafN)
# nothing significant
saveRDS(imputed.NW.perennial.graminoid.traits.depth.leafN, file = "./Results/perennial.graminoid.imputed.traits.NW.depth.leafN.rds")

imputed.NW.perennial.graminoid.traits.depth.leafN = readRDS("./Results/perennial.graminoid.imputed.traits.NW.depth.leafN.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.graminoid.traits.depth.leafN, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.graminoid.traits.depth.leafN, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.graminoid.traits.depth.leafN) # p = 0.873
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.graminoid.traits.depth.leafN)

imputed.NW.perennial.graminoid.traits.RTD.SRL = brm(cover.change ~ RTD.final*SRL.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                    family = gaussian(),
                                                    prior = priors,
                                                    data = imputed.NW.perennial.graminoid)

summary(imputed.NW.perennial.graminoid.traits.RTD.SRL)
# nothing significant
#saveRDS(imputed.NW.perennial.graminoid.traits.RTD.SRL, file = "./Results/perennial.graminoid.imputed.traits.NW.RTD.SRL.rds")

imputed.NW.perennial.graminoid.traits.RTD.SRL = readRDS("./Results/perennial.graminoid.imputed.traits.NW.RTD.SRL.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.graminoid.traits.RTD.SRL, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.graminoid.traits.RTD.SRL, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.graminoid.traits.RTD.SRL) # p = 0.920
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.graminoid.traits.RTD.SRL)

imputed.NW.perennial.graminoid.traits.leafN.RMF = brm(cover.change ~ leafN.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                      family = gaussian(),
                                                      prior = priors,
                                                      data = imputed.NW.perennial.graminoid)

summary(imputed.NW.perennial.graminoid.traits.leafN.RMF)
saveRDS(imputed.NW.perennial.graminoid.traits.leafN.RMF, file = "./Results/perennial.graminoid.imputed.traits.NW.leafN.RMF.rds")

imputed.NW.perennial.graminoid.traits.leafN.RMF = readRDS("./Results/perennial.graminoid.imputed.traits.NW.leafN.RMF.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.graminoid.traits.leafN.RMF, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.graminoid.traits.leafN.RMF, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.graminoid.traits.leafN.RMF) # p = 0.916
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.graminoid.traits.leafN.RMF)

#### Perennial graminoid Environment ####

priors <- c(prior(normal(0, 10), class = b))

imputed.NW.perennial.graminoid.traits.height.MAP.DSI = brm(cover.change ~ height.final*mean.MAP + height.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                           family = gaussian(),
                                                           prior = priors,
                                                           data = imputed.NW.perennial.graminoid)

summary(imputed.NW.perennial.graminoid.traits.height.MAP.DSI)
#saveRDS(imputed.NW.perennial.graminoid.traits.height.MAP.DSI, file = "./Results/imputed.NW.perennial.graminoid.traits.height.MAP.DSI.rds")

imputed.NW.perennial.graminoid.traits.height.MAP.DSI = readRDS("./Results/imputed.NW.perennial.graminoid.traits.height.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.graminoid.traits.height.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.graminoid.traits.height.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.graminoid.traits.height.MAP.DSI) # p = 0.916
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.graminoid.traits.height.MAP.DSI)

imputed.NW.perennial.graminoid.traits.leafN.MAP.DSI = brm(cover.change ~ leafN.final*mean.MAP + leafN.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                          family = gaussian(),
                                                          prior = priors,
                                                          data = imputed.NW.perennial.graminoid)

summary(imputed.NW.perennial.graminoid.traits.leafN.MAP.DSI)
#saveRDS(imputed.NW.perennial.graminoid.traits.leafN.MAP.DSI, file = "./Results/imputed.NW.perennial.graminoid.traits.leafN.MAP.DSI.rds")

imputed.NW.perennial.graminoid.traits.leafN.MAP.DSI = readRDS("./Results/imputed.NW.perennial.graminoid.traits.leafN.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.graminoid.traits.leafN.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.graminoid.traits.leafN.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.graminoid.traits.leafN.MAP.DSI) # p = 0.877
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.graminoid.traits.leafN.MAP.DSI)

imputed.NW.perennial.graminoid.traits.root.depth.MAP.DSI = brm(cover.change ~ root.depth.final*mean.MAP + root.depth.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                               family = gaussian(),
                                                               prior = priors,
                                                               data = imputed.NW.perennial.graminoid)

summary(imputed.NW.perennial.graminoid.traits.root.depth.MAP.DSI)
#saveRDS(imputed.NW.perennial.graminoid.traits.root.depth.MAP.DSI, file = "./Results/imputed.NW.perennial.graminoid.traits.depth.MAP.DSI.rds")

imputed.NW.perennial.graminoid.traits.root.depth.MAP.DSI = readRDS("./Results/imputed.NW.perennial.graminoid.traits.depth.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.graminoid.traits.root.depth.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.graminoid.traits.root.depth.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.graminoid.traits.root.depth.MAP.DSI) # p = 0.888
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.graminoid.traits.root.depth.MAP.DSI)

imputed.NW.perennial.graminoid.traits.RTD.MAP.DSI = brm(cover.change ~ RTD.final*mean.MAP + RTD.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                        family = gaussian(),
                                                        prior = priors,
                                                        data = imputed.NW.perennial.graminoid)

summary(imputed.NW.perennial.graminoid.traits.RTD.MAP.DSI)
#saveRDS(imputed.NW.perennial.graminoid.traits.RTD.MAP.DSI, file = "./Results/imputed.NW.perennial.graminoid.traits.RTD.MAP.DSI.rds")

imputed.NW.perennial.graminoid.traits.RTD.MAP.DSI = readRDS("./Results/imputed.NW.perennial.graminoid.traits.RTD.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.graminoid.traits.RTD.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.graminoid.traits.RTD.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.graminoid.traits.RTD.MAP.DSI) # p = 0.921
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.graminoid.traits.RTD.MAP.DSI)

imputed.NW.perennial.graminoid.traits.SRL.MAP.DSI = brm(cover.change ~ SRL.final*mean.MAP + SRL.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                        family = gaussian(),
                                                        prior = priors,
                                                        data = imputed.NW.perennial.graminoid)

summary(imputed.NW.perennial.graminoid.traits.SRL.MAP.DSI)
#saveRDS(imputed.NW.perennial.graminoid.traits.SRL.MAP.DSI, file = "./Results/imputed.NW.perennial.graminoid.traits.SRL.MAP.DSI.rds")

imputed.NW.perennial.graminoid.traits.SRL.MAP.DSI = readRDS("./Results/imputed.NW.perennial.graminoid.traits.SRL.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.graminoid.traits.SRL.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.graminoid.traits.SRL.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.graminoid.traits.SRL.MAP.DSI) # p = 0.919
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.graminoid.traits.SRL.MAP.DSI)

imputed.NW.perennial.graminoid.traits.RMF.MAP.DSI = brm(cover.change ~ RMF.final*mean.MAP + RMF.final*mean.DSI + (1|site_code) + (1|Taxon), 
                                                        family = gaussian(),
                                                        prior = priors,
                                                        data = imputed.NW.perennial.graminoid)

summary(imputed.NW.perennial.graminoid.traits.RMF.MAP.DSI)
#saveRDS(imputed.NW.perennial.graminoid.traits.RMF.MAP.DSI, file = "./Results/imputed.NW.perennial.graminoid.traits.RMF.MAP.DSI.rds")

imputed.NW.perennial.graminoid.traits.RMF.MAP.DSI = readRDS("./Results/imputed.NW.perennial.graminoid.traits.RMF.MAP.DSI.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.NW.perennial.graminoid.traits.RMF.MAP.DSI, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.NW.perennial.graminoid.traits.RMF.MAP.DSI, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.NW.perennial.graminoid.traits.RMF.MAP.DSI) # p = 0.916
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.NW.perennial.graminoid.traits.RMF.MAP.DSI)

#### ALL Backtransform ####

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)

# imputed traits with woody species
imputed.NW = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.outliersRM.csv", row.names = 1) %>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)

imputed.traits.NW.model = readRDS("./Results/all.imputed.traits.no_woody.rds")

imputed.NW$leafN.bt =(imputed.NW$leafN.final*7.114868) + 21.16957
imputed.NW$height.bt = (imputed.NW$height.final*0.2353906) + 0.3476824
imputed.NW$rootN.bt = (imputed.NW$rootN.final*4.488773) + 9.855027
imputed.NW$SLA.bt = (imputed.NW$SLA.final*8.581554) + 19.93883
imputed.NW$Depth.bt = (imputed.NW$root.depth.final*0.5333377) + 0.5784653
imputed.NW$Diam.bt = (imputed.NW$rootDiam.final*0.1598629) + 0.3501538
imputed.NW$SRL.bt = (imputed.NW$SRL.final*74.26099) + 106.4204
imputed.NW$RTD.bt = (imputed.NW$RTD.final*0.1159369) + 0.2349392
imputed.NW$RMF.bt = (imputed.NW$RMF.final*0.1058164) + 0.4007788
imputed.NW$DSI.bt = (imputed.NW$mean.DSI*0.1950438) + -0.4604709
imputed.NW$MAP.bt = (imputed.NW$mean.MAP*487.1451) + 691.6107

imputed.NW.2 = imputed.NW %>%
  dplyr::select(cover.change,leafN.bt:MAP.bt) %>%
  drop_na()

range(imputed.NW.2$cover.change) # -20.66667  21.24500
range(imputed.NW.2$leafN.bt) # 2.68952 44.19000
range(imputed.NW.2$height.bt) # 0.009000065 1.322333168
range(imputed.NW.2$rootN.bt) # 0.3800001 25.6283267
range(imputed.NW.2$SLA.bt) # 2.617597 46.320000
range(imputed.NW.2$Depth.bt) # 0.048000 3.709334
range(imputed.NW.2$Diam.bt) # 0.07840003 1.12601234
range(imputed.NW.2$SRL.bt) # 2.497204 417.749987
range(imputed.NW.2$RTD.bt) # 0.02325659 0.61999992
range(imputed.NW.2$RMF.bt) # 0.1211194 0.7328658
range(imputed.NW.2$DSI.bt) # -0.8722616  0.1242845
range(imputed.NW.2$MAP.bt) # 132.8 2366.0

all.effects = conditional_effects(imputed.traits.NW.model)

#### ALL leafN ####

leafN.effects = all.effects$leafN.final
leafN.effects$leafN.bt = (leafN.effects$leafN.final*7.114868) + 21.16957

x.value = c(-2,-1,0,1,2,3)
(x.value*7.114868) + 21.16957

all.leafN.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = leafN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.bt, y = estimate__), color = "#769370", size = 1.5) +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Leaf N (mg/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  xlim(2.68951,44.19001)+
  theme_classic()
all.leafN.plot

ggsave("./Plots/all.traits.NW.leafN.pdf", height = 3, width = 3)

#scale_x_continuous(limits = c(-2.6,3.3),
                   #breaks = seq(-2, 3, by = 1),
                   #labels = c(6.9, 14.1, 21.2, 28.3, 35.4, 42.5))+

#### All height ####

height.effects = all.effects$height.final
height.effects$height.bt = (height.effects$height.final*0.2353906) + 0.3476824

x.value = c(-1,0,1,2,3,4)
(x.value*0.2353906) + 0.3476824

all.height.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = height.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.bt, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = height.effects, aes(x = height.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Height (m)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  xlim(0.009000064,1.322333169)+
  theme_classic()
all.height.plot

ggsave("./Plots/all.traits.NW.height.pdf", height = 3, width = 3)

#### All rootN ####

rootN.effects = all.effects$rootN.final
rootN.effects$rootN.bt = (rootN.effects$rootN.final*4.488773) + 9.855027

x.value = c(-2,-1,0,1,2,3)
(x.value*4.488773) + 9.855027

all.rootN.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = rootN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.bt, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Root N (mg/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  xlim(0.3800000,25.6283268)+
  theme_classic()
all.rootN.plot

ggsave("./Plots/all.traits.NW.rootN.pdf", height = 3, width = 3)

#### All SLA ####

SLA.effects = all.effects$SLA.final
SLA.effects$SLA.bt = (SLA.effects$SLA.final*8.581554) + 19.93883

x.value = c(-2,-1,0,1,2,3)
(x.value*8.581554) + 19.93883

all.SLA.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = SLA.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.bt, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = expression("Specific Leaf Area (m"^2*"/kg)"), y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  xlim(2.617596, 46.320001)+
  theme_classic()
all.SLA.plot

ggsave("./Plots/all.traits.NW.SLA.pdf", height = 3, width = 3)

#### All Depth ####

root.depth.effects = all.effects$root.depth.final
root.depth.effects$Depth.bt = (root.depth.effects$root.depth.final*0.5333377) + 0.5784653

x.value = c(0,1,2,3,4,5)
(x.value*0.5333377) + 0.5784653

all.root.depth.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = Depth.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = Depth.bt, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = Depth.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Rooting Depth (m)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  xlim(0.0479999,3.709335)+
  theme_classic()
all.root.depth.plot

ggsave("./Plots/all.traits.NW.root.depth.pdf", height = 3, width = 3)

#### All Diameter ####

rootDiam.effects = all.effects$rootDiam.final
rootDiam.effects$Diam.bt = (rootDiam.effects$rootDiam.final*0.1598629) + 0.3501538

x.value = c(-1,0,1,2,3,4)
(x.value*0.1598629) + 0.3501538

all.rootDiam.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = Diam.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = Diam.bt, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = Diam.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Root Diameter (mm)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  xlim(0.07840002,1.12601235)+
  theme_classic()
all.rootDiam.plot

ggsave("./Plots/all.traits.NW.rootDiam.pdf", height = 3, width = 3)

#### All SRL ####

SRL.effects = all.effects$SRL.final
SRL.effects$SRL.bt = (SRL.effects$SRL.final*74.26099) + 106.4204

x.value = c(-1,0,1,2,3,4)
(x.value*74.26099) + 106.4204

all.SRL.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = SRL.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.bt, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Specific Root Length (m/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  xlim(2.497203,417.749988)+
  theme_classic()
all.SRL.plot

ggsave("./Plots/all.traits.NW.SRL.pdf", height = 3, width = 3)

#### All RTD ####

RTD.effects = all.effects$RTD.final
RTD.effects$RTD.bt = (RTD.effects$RTD.final*0.1159369) + 0.2349392

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1159369) + 0.2349392

all.RTD.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = RTD.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.bt, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = expression("Root Tissue Density (g/cm"^3*")"), y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  xlim(0.02325658, 0.61999993)+
  theme_classic()
all.RTD.plot

ggsave("./Plots/all.traits.NW.RTD.pdf", height = 3, width = 3)

#### All RMF ####

RMF.effects = all.effects$RMF.final
RMF.effects$RMF.bt = (RMF.effects$RMF.final*0.1058164) + 0.4007788

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1058164) + 0.4007788

all.RMF.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = RMF.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RMF.effects, aes(x = RMF.bt, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RMF.effects, aes(x = RMF.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Root Mass Fraction (g/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  xlim(0.1211193,0.7328659)+
  theme_classic()
all.RMF.plot

ggsave("./Plots/all.traits.NW.RMF.pdf", height = 3, width = 3)

#### All DSI ####

DSI.effects = all.effects$mean.DSI
DSI.effects$DSI.bt = (DSI.effects$mean.DSI*0.1950438) + -0.4604709

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1950438) + -0.4604709

all.DSI.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = DSI.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = DSI.effects, aes(x = DSI.bt, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = DSI.effects, aes(x = DSI.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Drought Severity Index", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  xlim(-0.8722616,0.1242846)+
   theme_classic()
all.DSI.plot

ggsave("./Plots/all.traits.NW.DSI.pdf", height = 3, width = 3)

#### All MAP ####

MAP.effects = all.effects$mean.MAP
MAP.effects$MAP.bt= (MAP.effects$mean.MAP*487.1451) + 691.6107

x.value = c(-1,0,1,2,3)
(x.value*487.1451) + 691.6107

all.MAP.plot = ggplot() +
  geom_point(data = imputed.NW.2, aes(x = MAP.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = MAP.effects, aes(x = MAP.bt, y = estimate__), color = "#769370", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = MAP.effects, aes(x = MAP.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#769370") +  
  labs(x = "Mean Annual Precipitation (mm)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  xlim(132.7,2366.1)+
  theme_classic()
all.MAP.plot

ggsave("./Plots/all.traits.NW.MAP.pdf", height = 3, width = 3)

#### ALL interactions ####

imputed.NW.traits.height.leafN = readRDS("./Results/all.imputed.traits.NW.height.leafN.rds")

height.leafN.effect = conditional_effects(imputed.NW.traits.height.leafN, effects = "height.final:leafN.final")$`height.final:leafN.final`
height.leafN.effect$height.bt = (height.leafN.effect$height.final*0.2353906) + 0.3476824

# only plot high and low values

height.leafN.effect.2 = height.leafN.effect %>%
  filter(effect2__ %in% c(-0.98,1.03))

# height
x.value = c(-1,0,1,2,3,4)
(x.value*0.2353906) + 0.3476824

# leafN
x.value = c(-0.98,1.03)
(x.value*7.114868) + 21.16957


all.height.x.leafN.plot = ggplot(data = height.leafN.effect.2, aes(x = height.bt, y = estimate__, group = as.factor(effect2__))) +
  geom_line(aes(color = as.factor(effect2__)), size = 1.5,show.legend = FALSE) +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = as.factor(effect2__)), alpha = 0.5,show.legend = FALSE) +  
  geom_line(aes(y = upper__, color = factor(effect2__)), size = 1.2) +
  geom_line(aes(y = lower__, color = factor(effect2__)), size = 1.2) + 
  labs(x = "Height (m)", y = "Cover Change (%)", color = "Leaf N (mg/g)") +
  scale_colour_manual(values = c("black", "#769370"), labels = c("28.5","14.2"))+
  scale_fill_manual(values = c("black", "#769370"))+
  theme_classic()
all.height.x.leafN.plot

ggsave("./Plots/all.traits.NW.height.leafN.pdf", height = 3, width = 3)
ggsave("./Plots/all.traits.NW.height.leafN.legend.pdf", height = 3, width = 3)

imputed.NW.traits.leafN.RMF = readRDS("./Results/all.imputed.NW.traits.leafN.RMF.rds")

leafN.RMF.effect = conditional_effects(imputed.NW.traits.leafN.RMF, effects = "leafN.final:RMF.final")$`leafN.final:RMF.final`
leafN.RMF.effect$leafN.bt = (leafN.RMF.effect$leafN.final*7.114868) + 21.16957

# only plot high and low values

leafN.RMF.effect.2 = leafN.RMF.effect %>%
  filter(effect2__ %in% c(-1.02,0.99))
# leafN
x.value = c(-2,-1,0,1,2,3)
(x.value*7.114868) + 21.16957

# RMF
x.value = c(-1.02,0.99)
(x.value*0.1058164) + 0.4007788

all.leafN.x.RMF.plot = ggplot(data = leafN.RMF.effect.2, aes(x = leafN.bt, y = estimate__, group = as.factor(effect2__))) +
  geom_line(aes(color = as.factor(effect2__)), size = 1.5,show.legend = FALSE) +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = as.factor(effect2__)), alpha = 0.5,show.legend = FALSE) +  
  geom_line(aes(y = upper__, color = factor(effect2__)), size = 1.2) +
  geom_line(aes(y = lower__, color = factor(effect2__)), size = 1.2) + 
  labs(x = "Leaf N (mg/g)", y = "Cover Change (%)", color = "RMF (g/g)") +
  scale_colour_manual(values = c("black", "#769370"), labels = c("0.5","0.3"))+
  scale_fill_manual(values = c("black", "#769370"))+
  xlim(2.68951,44.19001)+
  theme_classic()
all.leafN.x.RMF.plot

ggsave("./Plots/all.traits.NW.leafN.RMF.pdf", height = 3, width = 3)
ggsave("./Plots/all.traits.NW.leafN.RMF.legend.pdf", height = 3, width = 3)

#### ALL environment interactions ####

all.imputed.NW.traits.leafN.DSI = readRDS("./Results/all.imputed.NW.traits.leafN.MAP.DSI.rds")

leafN.DSI.effect = conditional_effects(all.imputed.NW.traits.leafN.DSI, effects = "leafN.final:mean.DSI")$`leafN.final:mean.DSI`
leafN.DSI.effect$leafN.bt= (leafN.DSI.effect$leafN.final*7.114868) + 21.16957

# only plot high and low values

leafN.DSI.effect.2 = leafN.DSI.effect %>%
  dplyr::filter(effect2__ %in% c(-1,1))

# leafN
x.value = c(-3,0,1,2,3)
(x.value*7.150147) + 22.84208

# DSI
x.value = c(-1,1)
(x.value*0.1950438) + -0.4604709

all.leafN.x.DSI.plot = ggplot(data = leafN.DSI.effect.2, aes(x = leafN.bt, y = estimate__, group = as.factor(effect2__))) +
  geom_line(aes(color = as.factor(effect2__)), size = 1.5,show.legend = FALSE) +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = as.factor(effect2__)), alpha = 0.5,show.legend = FALSE) +  
  geom_line(aes(y = upper__, color = factor(effect2__)), size = 1.2) +
  geom_line(aes(y = lower__, color = factor(effect2__)), size = 1.2) + 
  labs(x = "Leaf N (mg/g)", y = "Cover Change (%)", color = "DSI") +
  scale_colour_manual(values = c("black", "#769370"), labels = c("-0.27","-0.66"))+
  scale_fill_manual(values = c("black", "#769370"))+
  xlim(2.68951,44.19001)+
  theme_classic()
all.leafN.x.DSI.plot

ggsave("./Plots/all.traits.NW.leafN.DSI.pdf", height = 3, width = 3)
ggsave("./Plots/all.traits.NW.leafN.DSI.legend.pdf", height = 3, width = 3)


#### ANNUAL Backtransform ####

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)

# imputed traits with woody species
imputed.NW.annual = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.annuals.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)

annual.traits.NW.model = readRDS("./Results/annual.imputed.traits.no_woody.rds")

imputed.NW.annual$leafN.bt =(imputed.NW.annual$leafN.final*7.384253) + 22.54886
imputed.NW.annual$height.bt = (imputed.NW.annual$height.final*0.2308751) + 0.3156657
imputed.NW.annual$rootN.bt = (imputed.NW.annual$rootN.final*4.33202) + 9.997142
imputed.NW.annual$SLA.bt = (imputed.NW.annual$SLA.final*9.312314) + 22.81754
imputed.NW.annual$Depth.bt = (imputed.NW.annual$root.depth.final*0.4346265) + 0.5402869
imputed.NW.annual$Diam.bt = (imputed.NW.annual$rootDiam.final*0.1367222) + 0.3407223
imputed.NW.annual$SRL.bt = (imputed.NW.annual$SRL.final*73.49115) + 125.697
imputed.NW.annual$RTD.bt = (imputed.NW.annual$RTD.final*0.1111574) + 0.2026466
imputed.NW.annual$RMF.bt = (imputed.NW.annual$RMF.final*0.1037188) + 0.3934184
imputed.NW.annual$DSI.bt = (imputed.NW.annual$mean.DSI*0.2300417) + -0.4690337
imputed.NW.annual$MAP.bt = (imputed.NW.annual$mean.MAP*293.1676) + 477.7734

annual.imputed.NW.2 = imputed.NW.annual %>%
  dplyr::select(cover.change,leafN.bt:MAP.bt) %>%
  drop_na()

annual.effects = conditional_effects(annual.traits.NW.model)

#### Annual leafN ####
leafN.effects = annual.effects$leafN.final
leafN.effects$leafN.bt= (leafN.effects$leafN.final*7.384253) + 22.54886

attr(annual.imputed.NW.2$leafN.final, "scaled:scale")
attr(annual.imputed.NW.2$leafN.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*7.384253) + 22.54886

annual.leafN.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = leafN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Leaf N (mg/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  xlim(2.68951,44.19001)+
  theme_classic()
annual.leafN.plot

ggsave("./Plots/annual.traits.NW.leafN.pdf", height = 3, width = 3)

#### Annual height ####

height.effects = annual.effects$height.final
height.effects$height.bt= (height.effects$height.final*0.2308751) + 0.3156657

attr(annual.imputed.NW.2$height.final, "scaled:scale")
attr(annual.imputed.NW.2$height.final, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*0.2308751) + 0.3156657

annual.height.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = height.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = height.effects, aes(x = height.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Height (m)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.009000064,1.322333169)+
  theme_classic()
annual.height.plot

ggsave("./Plots/annual.traits.NW.height.pdf", height = 3, width = 3)


#### Annual RootN ####
rootN.effects = annual.effects$rootN.final
rootN.effects$rootN.bt= (rootN.effects$rootN.final*4.33202) + 9.997142

attr(annual.imputed.NW.2$rootN.final, "scaled:scale")
attr(annual.imputed.NW.2$rootN.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*4.33202) + 9.997142

annual.rootN.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = rootN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Root N (mg/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.3800000,25.6283268)+
  theme_classic()
annual.rootN.plot

ggsave("./Plots/annual.traits.NW.rootN.pdf", height = 3, width = 3)

#### Annual SLA ####

SLA.effects = annual.effects$SLA.final
SLA.effects$SLA.bt= (SLA.effects$SLA.final*9.312314) + 22.81754

attr(annual.imputed.NW.2$SLA.final, "scaled:scale")
attr(annual.imputed.NW.2$SLA.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*9.312314) + 22.81754

annual.SLA.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = SLA.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = expression("Specific Leaf Area (m"^2*"/kg)"), y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.617596,46.320001)+
  theme_classic()
annual.SLA.plot

ggsave("./Plots/annual.traits.NW.SLA.pdf", height = 3, width = 3)

#### Annual depth ####

root.depth.effects = annual.effects$root.depth.final
root.depth.effects$Depth.bt= (root.depth.effects$root.depth.final*0.4346265) + 0.5402869

attr(annual.imputed.NW.2$root.depth.final, "scaled:scale")
attr(annual.imputed.NW.2$root.depth.final, "scaled:center")

x.value = c(0,1,2,3,4,5)
(x.value*0.4346265) + 0.5402869

annual.root.depth.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = Depth.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = Depth.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = Depth.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Rooting Depth (m)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.0479999,3.709335)+
  theme_classic()
annual.root.depth.plot

ggsave("./Plots/annual.traits.NW.root.depth.pdf", height = 3, width = 3)

#### Annual Diam ####

rootDiam.effects = annual.effects$rootDiam.final
rootDiam.effects$Diam.bt= (rootDiam.effects$rootDiam.final*0.1367222) + 0.3407223

attr(annual.imputed.NW.2$rootDiam.final, "scaled:scale")
attr(annual.imputed.NW.2$rootDiam.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1367222) + 0.3407223

annual.rootDiam.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = Diam.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = Diam.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = Diam.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Root Diameter (mm)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.07840002,1.12601235)+
  theme_classic()
annual.rootDiam.plot

ggsave("./Plots/annual.traits.NW.rootDiam.pdf", height = 3, width = 3)

#### Annual SRL ####

SRL.effects = annual.effects$SRL.final
SRL.effects$SRL.bt= (SRL.effects$SRL.final*73.49115) + 125.697

attr(annual.imputed.NW.2$SRL.final, "scaled:scale")
attr(annual.imputed.NW.2$SRL.final, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*73.49115) + 125.697

annual.SRL.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = SRL.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Specific Root Length (m/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.497203,417.749988)+
  theme_classic()
annual.SRL.plot

ggsave("./Plots/annual.traits.NW.SRL.pdf", height = 3, width = 3)

#### Annual RTD ####

RTD.effects = annual.effects$RTD.final
RTD.effects$RTD.bt= (RTD.effects$RTD.final*0.1111574) + 0.2026466

attr(annual.imputed.NW.2$RTD.final, "scaled:scale")
attr(annual.imputed.NW.2$RTD.final, "scaled:center")

x.value = c(-1,0,1,2,3)
(x.value*0.1111574) + 0.2026466

annual.RTD.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = RTD.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = expression("Root Tissue Density (g/cm"^3*")"), y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.02325658,0.61999993)+
  theme_classic()
annual.RTD.plot

ggsave("./Plots/annual.traits.NW.RTD.pdf", height = 3, width = 3)

#### Annual RMF ####

RMF.effects = annual.effects$RMF.final
RMF.effects$RMF.bt= (RMF.effects$RMF.final*0.1037188) + 0.3934184

attr(annual.imputed.NW.2$RMF.final, "scaled:scale")
attr(annual.imputed.NW.2$RMF.final, "scaled:center")

x.value = c(-3,-2,-1,0,1,2,3)
(x.value*0.1037188) + 0.3934184

annual.RMF.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = RMF.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RMF.effects, aes(x = RMF.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RMF.effects, aes(x = RMF.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Root Mass Fraction (g/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.1211193,0.7328659)+
  theme_classic()
annual.RMF.plot

ggsave("./Plots/annual.traits.NW.RMF.pdf", height = 3, width = 3)

#### Annual DSI ####

DSI.effects = annual.effects$mean.DSI
DSI.effects$DSI.bt= (DSI.effects$mean.DSI*0.2300417) + -0.4690337

attr(annual.imputed.NW.2$mean.DSI, "scaled:scale")
attr(annual.imputed.NW.2$mean.DSI, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.2300417) + -0.4690337

annual.DSI.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = DSI.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = DSI.effects, aes(x = DSI.bt, y = estimate__), color = "#979461", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = DSI.effects, aes(x = DSI.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Drought Severity Index", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(-0.8722616,0.1242846)+
  theme_classic()
annual.DSI.plot

ggsave("./Plots/annual.traits.NW.DSI.pdf", height = 3, width = 3)

#### Annual MAP ####

MAP.effects = annual.effects$mean.MAP
MAP.effects$MAP.bt= (MAP.effects$mean.MAP*293.1676) + 477.7734

attr(annual.imputed.NW.2$mean.MAP, "scaled:scale")
attr(annual.imputed.NW.2$mean.MAP, "scaled:center")

x.value = c(-1,0,1,2,3)
(x.value*293.1676) + 477.7734

annual.MAP.plot = ggplot() +
  geom_point(data = annual.imputed.NW.2, aes(x = MAP.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = MAP.effects, aes(x = MAP.bt, y = estimate__), color = "#979461", size = 1.5) +  
  geom_ribbon(data = MAP.effects, aes(x = MAP.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#979461") +  
  labs(x = "Mean Annual Precipitation (mm)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(132.7,2366.1)+
  theme_classic()
annual.MAP.plot

ggsave("./Plots/annual.traits.NW.MAP.pdf", height = 3, width = 3)

#### ANNUAL interactions ####

annual.imputed.NW.traits.leafN.RMF = readRDS("./Results/annual.imputed.NW.traits.leafN.RMF.rds")

leafN.RMF.effect = conditional_effects(annual.imputed.NW.traits.leafN.RMF, effects = "leafN.final:RMF.final")$`leafN.final:RMF.final`
leafN.RMF.effect$leafN.bit= (leafN.RMF.effect$leafN.final*7.384253) + 22.54886

# only plot high and low values

leafN.RMF.effect.2 = leafN.RMF.effect %>%
  filter(effect2__ %in% c(-1.05,0.97))

# leafN
x.value = c(-2,-1,0,1,2,3)
(x.value*7.384253) + 22.54886

# RMF
x.value = c(-1.05,0.97)
(x.value*0.1037188) + 0.3934184

annual.leafN.x.RMF.plot = ggplot(data = leafN.RMF.effect.2, aes(x = leafN.bit, y = estimate__, group = as.factor(effect2__))) +
  geom_line(aes(color = as.factor(effect2__)), size = 1.5,show.legend = FALSE) +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = as.factor(effect2__)), alpha = 0.5,show.legend = FALSE) +  
  geom_line(aes(y = upper__, color = factor(effect2__)), size = 1.2,show.legend = FALSE) +
  geom_line(aes(y = lower__, color = factor(effect2__)), size = 1.2,show.legend = FALSE) + 
  labs(x = "Leaf N (mg/g)", y = "Cover Change (%)", color = "RMF (g/g)") +
  scale_colour_manual(values = c("black", "#979461"), labels = c("0.49","0.28"))+
  scale_fill_manual(values = c("black", "#979461"))+
  #xlim(2.68951,44.19001)+
  theme_classic()
annual.leafN.x.RMF.plot

ggsave("./Plots/annual.traits.NW.leafN.RMF.pdf", height = 3, width = 3)
ggsave("./Plots/annual.traits.NW.leafN.RMF.legend.pdf", height = 3, width = 3)

#### ANNUAL environment interactions ####

annual.imputed.NW.traits.leafN.DSI = readRDS("./Results/imputed.NW.annual.traits.leafN.MAP.DSI.rds")

leafN.DSI.effect = conditional_effects(annual.imputed.NW.traits.leafN.DSI, effects = "leafN.final:mean.DSI")$`leafN.final:mean.DSI`
leafN.DSI.effect$leafN.bt= (leafN.DSI.effect$leafN.final*7.384253) + 22.54886

# only plot high and low values

leafN.DSI.effect.2 = leafN.DSI.effect %>%
  dplyr::filter(effect2__ %in% c(-1,1))

# leafN
x.value = c(-3,0,1,2,3)
(x.value*7.150147) + 22.84208

# DSI
x.value = c(-1,1)
(x.value*0.2300417) + -0.4690337

annual.leafN.x.DSI.plot = ggplot(data = leafN.DSI.effect.2, aes(x = leafN.bt, y = estimate__, group = as.factor(effect2__))) +
  geom_line(aes(color = as.factor(effect2__)), size = 1.5,show.legend = FALSE) +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = as.factor(effect2__)), alpha = 0.5,show.legend = FALSE) +  
  geom_line(aes(y = upper__, color = factor(effect2__)), size = 1.2,show.legend = FALSE) +
  geom_line(aes(y = lower__, color = factor(effect2__)), size = 1.2,show.legend = FALSE) + 
  labs(x = "Leaf N (mg/g)", y = "Cover Change (%)", color = "DSI") +
  scale_colour_manual(values = c("black", "#979461"), labels = c("-0.24","-0.70"))+
  scale_fill_manual(values = c("black", "#979461"))+
  #xlim(2.68951,44.19001)+
  theme_classic()
annual.leafN.x.DSI.plot

ggsave("./Plots/annual.traits.NW.leafN.DSI.pdf", height = 3, width = 3)
ggsave("./Plots/annual.traits.NW.leafN.DSI.legend.pdf", height = 3, width = 3)

annual.imputed.NW.traits.height.MAP = readRDS("./Results/imputed.NW.annual.traits.height.MAP.DSI.rds")

height.MAP.effect = conditional_effects(annual.imputed.NW.traits.height.MAP, effects = "height.final:mean.MAP")$`height.final:mean.MAP`
height.MAP.effect$height.bt= (height.MAP.effect$height.final*0.2308751) + 0.3156657

# only plot high and low values

height.MAP.effect.2 = height.MAP.effect %>%
  dplyr::filter(effect2__ %in% c(-1,1))

# height
x.value = c(-3,0,1,2,3)
(x.value*7.150147) + 22.84208

# MAP
x.value = c(-1,1)
(x.value*293.1676) + 477.7734

annual.height.x.MAP.plot = ggplot(data = height.MAP.effect.2, aes(x = height.bt, y = estimate__, group = as.factor(effect2__))) +
  geom_line(aes(color = as.factor(effect2__)), size = 1.5,show.legend = FALSE) +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = as.factor(effect2__)), alpha = 0.5,show.legend = FALSE) +  
  geom_line(aes(y = upper__, color = factor(effect2__)), size = 1.2,show.legend = FALSE) +
  geom_line(aes(y = lower__, color = factor(effect2__)), size = 1.2,show.legend = FALSE) + 
  labs(x = "Height (m)", y = "Cover Change (%)", color = "MAP (mm)") +
  scale_colour_manual(values = c("black","#979461"), labels = c("770","185"))+
  scale_fill_manual(values = c("black", "#979461"))+
  #xlim(0.0090000,1.26)+
  theme_classic()
annual.height.x.MAP.plot

ggsave("./Plots/annual.traits.NW.height.MAP.pdf", height = 3, width = 3)
ggsave("./Plots/annual.traits.NW.height.MAP.legend.pdf", height = 3, width = 3)


#### PERENNIAL Backtransform ####

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)
imputed.NW.perennial = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennials.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)

perennial.traits.NW.model = readRDS("./Results/perennial.imputed.traits.no_woody.rds")

imputed.NW.perennial$leafN.bt =(imputed.NW.perennial$leafN.final*6.812995) + 20.47022
imputed.NW.perennial$height.bt = (imputed.NW.perennial$height.final*0.2365653) + 0.3627378
imputed.NW.perennial$rootN.bt = (imputed.NW.perennial$rootN.final*4.532848) + 9.765547
imputed.NW.perennial$SLA.bt = (imputed.NW.perennial$SLA.final*7.966249) + 18.75888
imputed.NW.perennial$Depth.bt = (imputed.NW.perennial$root.depth.final*0.5706892) + 0.5903092
imputed.NW.perennial$Diam.bt = (imputed.NW.perennial$rootDiam.final*0.166311) + 0.3535401
imputed.NW.perennial$SRL.bt = (imputed.NW.perennial$SRL.final*73.08878) + 98.7131
imputed.NW.perennial$RTD.bt = (imputed.NW.perennial$RTD.final*0.1151264) + 0.249861
imputed.NW.perennial$RMF.bt = (imputed.NW.perennial$RMF.final*0.1053099) + 0.4042604
imputed.NW.perennial$DSI.bt = (imputed.NW.perennial$mean.DSI*0.1779261) + -0.4547629
imputed.NW.perennial$MAP.bt = (imputed.NW.perennial$mean.MAP*522.5359) + 772.0729

perennial.imputed.NW.2 = imputed.NW.perennial %>%
  dplyr::select(cover.change,leafN.bt:MAP.bt) %>%
  drop_na()

perennial.effects = conditional_effects(perennial.traits.NW.model)

#### Perennial LeafN ####

leafN.effects = perennial.effects$leafN.final
leafN.effects$leafN.bt= (leafN.effects$leafN.final*6.812995) + 20.47022

attr(perennial.imputed.NW.2$leafN.final, "scaled:scale")
attr(perennial.imputed.NW.2$leafN.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*6.812995) + 20.47022

perennial.leafN.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = leafN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Leaf N (mg/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.68951,44.19001)+
  theme_classic()
perennial.leafN.plot

ggsave("./Plots/perennial.traits.NW.leafN.pdf", height = 3, width = 3)

#### Perennial Height ####

height.effects = perennial.effects$height.final
height.effects$height.bt= (height.effects$height.final*0.2365653) + 0.3627378

attr(perennial.imputed.NW.2$height.final, "scaled:scale")
attr(perennial.imputed.NW.2$height.final, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*0.2365653) + 0.3627378

perennial.height.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = height.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = height.effects, aes(x = height.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Height (m)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.009000064,1.322333169)+
  theme_classic()
perennial.height.plot

ggsave("./Plots/perennial.traits.NW.height.pdf", height = 3, width = 3)

#### Perennial RootN ####

rootN.effects = perennial.effects$rootN.final
rootN.effects$rootN.bt= (rootN.effects$rootN.final*4.532848) + 9.765547

attr(perennial.imputed.NW.2$rootN.final, "scaled:scale")
attr(perennial.imputed.NW.2$rootN.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*4.532848) + 9.765547

perennial.rootN.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = rootN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Root N (mg/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.3800000,25.6283268)+
  theme_classic()
perennial.rootN.plot

ggsave("./Plots/perennial.traits.NW.rootN.pdf", height = 3, width = 3)

#### Perennial SLA ####

SLA.effects = perennial.effects$SLA.final
SLA.effects$SLA.bt= (SLA.effects$SLA.final*7.966249) + 18.75888

attr(perennial.imputed.NW.2$SLA.final, "scaled:scale")
attr(perennial.imputed.NW.2$SLA.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*7.966249) + 18.75888

perennial.SLA.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = SLA.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = expression("Specific Leaf Area (m"^2*"/kg)"), y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.617596,46.320001)+
  theme_classic()
perennial.SLA.plot

ggsave("./Plots/perennial.traits.NW.SLA.pdf", height = 3, width = 3)

#### Perennial Depth ####

root.depth.effects = perennial.effects$root.depth.final
root.depth.effects$Depth.bt= (root.depth.effects$root.depth.final*0.5706892) + 0.5903092

attr(perennial.imputed.NW.2$root.depth.final, "scaled:scale")
attr(perennial.imputed.NW.2$root.depth.final, "scaled:center")

x.value = c(-1,0,1,2,3,4,5)
(x.value*0.5706892) + 0.5903092

perennial.root.depth.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = Depth.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = Depth.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = Depth.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Rooting Depth (m)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.0479999,3.709335)+
  theme_classic()
perennial.root.depth.plot

ggsave("./Plots/perennial.traits.NW.root.depth.pdf", height = 3, width = 3)

#### Perennial Diam ####

rootDiam.effects = perennial.effects$rootDiam.final
rootDiam.effects$Diam.bt= (rootDiam.effects$rootDiam.final*0.166311) + 0.3535401

attr(perennial.imputed.NW.2$rootDiam.final, "scaled:scale")
attr(perennial.imputed.NW.2$rootDiam.final, "scaled:center")

x.value = c(-1,0,1,2,3,4,5)
(x.value*0.166311) + 0.3535401

perennial.rootDiam.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = Diam.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = Diam.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = Diam.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Root Diameter (mm)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.07840002,1.12601235)+
  theme_classic()
perennial.rootDiam.plot

ggsave("./Plots/perennial.traits.NW.rootDiam.pdf", height = 3, width = 3)

#### Perennial SRL ####

SRL.effects = perennial.effects$SRL.final
SRL.effects$SRL.bt= (SRL.effects$SRL.final*73.08878) + 98.7131

attr(perennial.imputed.NW.2$SRL.final, "scaled:scale")
attr(perennial.imputed.NW.2$SRL.final, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*73.08878) + 98.7131

perennial.SRL.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = SRL.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Specific Root Length (m/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.497203,417.749988)+
  theme_classic()
perennial.SRL.plot

ggsave("./Plots/perennial.traits.NW.SRL.pdf", height = 3, width = 3)

#### Perennial RTD ####

RTD.effects = perennial.effects$RTD.final
RTD.effects$RTD.bt= (RTD.effects$RTD.final*0.1151264) + 0.249861

attr(perennial.imputed.NW.2$RTD.final, "scaled:scale")
attr(perennial.imputed.NW.2$RTD.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1151264) + 0.249861

perennial.RTD.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = RTD.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = expression("Root Tissue Density (g/cm"^3*")"), y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.02325658,0.61999993)+
  theme_classic()
perennial.RTD.plot

ggsave("./Plots/perennial.traits.NW.RTD.pdf", height = 3, width = 3)

#### Perennial RMF ####

RMF.effects = perennial.effects$RMF.final
RMF.effects$RMF.bt= (RMF.effects$RMF.final*0.1053099) + 0.4042604

attr(perennial.imputed.NW.2$RMF.final, "scaled:scale")
attr(perennial.imputed.NW.2$RMF.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1053099) + 0.4042604

perennial.RMF.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = RMF.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RMF.effects, aes(x = RMF.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RMF.effects, aes(x = RMF.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Root Mass Fraction (g/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.1211193,0.7328659)+
  theme_classic()
perennial.RMF.plot

ggsave("./Plots/perennial.traits.NW.RMF.pdf", height = 3, width = 3)

#### Perennial DSI ####

DSI.effects = perennial.effects$mean.DSI
DSI.effects$DSI.bt= (DSI.effects$mean.DSI*0.1779261) + -0.4547629

attr(perennial.imputed.NW.2$mean.DSI, "scaled:scale")
attr(perennial.imputed.NW.2$mean.DSI, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*0.1779261) + -0.4547629

perennial.DSI.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = DSI.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = DSI.effects, aes(x = DSI.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = DSI.effects, aes(x = DSI.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Drought Severity Index", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(-0.8722616,0.1242846)+
  theme_classic()
perennial.DSI.plot

ggsave("./Plots/perennial.traits.NW.DSI.pdf", height = 3, width = 3)

#### Perennial MAP ####

MAP.effects = perennial.effects$mean.MAP
MAP.effects$MAP.bt= (MAP.effects$mean.MAP*522.5359) + 772.0729

attr(perennial.imputed.NW.2$mean.MAP, "scaled:scale")
attr(perennial.imputed.NW.2$mean.MAP, "scaled:center")

x.value = c(-1,0,1,2,3)
(x.value*522.5359) + 772.0729

perennial.MAP.plot = ggplot() +
  geom_point(data = perennial.imputed.NW.2, aes(x = MAP.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = MAP.effects, aes(x = MAP.bt, y = estimate__), color = "#F1C646", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = MAP.effects, aes(x = MAP.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F1C646") +  
  labs(x = "Mean Annual Precipitation (mm)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(132.7,2366.1)+
  theme_classic()
perennial.MAP.plot

ggsave("./Plots/perennial.traits.NW.MAP.pdf", height = 3, width = 3)

#### PERENNIAL interactions ####

perennial.imputed.NW.traits.height.leafN = readRDS("./Results/perennial.imputed.traits.NW.height.leafN.rds")

height.leafN.effect = conditional_effects(perennial.imputed.NW.traits.height.leafN, effects = "height.final:leafN.final")$`height.final:leafN.final`
height.leafN.effect$height.bt= (height.leafN.effect$height.final*0.2365653) + 0.3627378

# only plot high and low values

height.leafN.effect.2 = height.leafN.effect %>%
  filter(effect2__ %in% c(-0.97,1.02))

# height
x.value = c(-1,0,1,2,3,4)
(x.value*0.2365653) + 0.3627378

# leafN
x.value = c(-0.97,1.02)
(x.value*6.812995) + 20.47022

perennial.height.x.leafN.plot = ggplot(data = height.leafN.effect.2, aes(x = height.bt, y = estimate__, group = as.factor(effect2__))) +
  geom_line(aes(color = as.factor(effect2__)), size = 1.5,show.legend = FALSE) +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = as.factor(effect2__)), alpha = 0.5,show.legend = FALSE) +  
  geom_line(aes(y = upper__, color = factor(effect2__)), size = 1.2) +
  geom_line(aes(y = lower__, color = factor(effect2__)), size = 1.2) + 
  labs(x = "Height (m)", y = "Cover Change (%)", color = "Leaf N (mg/g)") +
  scale_colour_manual(values = c("black", "#F1C646"), labels = c("27.41","13.86"))+
  scale_fill_manual(values = c("black", "#F1C646"))+
  #xlim(0.009000064,1.322333169)+
  theme_classic()
perennial.height.x.leafN.plot

ggsave("./Plots/perennial.traits.NW.height.leafN.pdf", height = 3, width = 3)
ggsave("./Plots/perennial.traits.NW.height.leafN.legend.pdf", height = 3, width = 3)

#### FORB Backtransform ####

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)
imputed.NW.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.forbs.legume.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)

forb.traits.NW.model = readRDS("./Results/forb.legume.imputed.traits.no_woody.rds")

imputed.NW.forb$leafN.bt =(imputed.NW.forb$leafN.final*7.150147) + 22.84208
imputed.NW.forb$height.bt = (imputed.NW.forb$height.final*0.2405401) + 0.3256687
imputed.NW.forb$rootN.bt = (imputed.NW.forb$rootN.final*4.715357) + 11.3382
imputed.NW.forb$SLA.bt = (imputed.NW.forb$SLA.final*8.061883) + 20.18927
imputed.NW.forb$Depth.bt = (imputed.NW.forb$root.depth.final*0.600499) + 0.6551094
imputed.NW.forb$Diam.bt = (imputed.NW.forb$rootDiam.final*0.152751) + 0.3553141
imputed.NW.forb$SRL.bt = (imputed.NW.forb$SRL.final*70.41695) + 103.2147
imputed.NW.forb$RTD.bt = (imputed.NW.forb$RTD.final*0.1178487) + 0.2169321
imputed.NW.forb$RMF.bt = (imputed.NW.forb$RMF.final*0.1159929) + 0.3917262
imputed.NW.forb$DSI.bt = (imputed.NW.forb$mean.DSI*0.1957044) + -0.4647529
imputed.NW.forb$MAP.bt = (imputed.NW.forb$mean.MAP*433.6646) + 664.4671

forb.imputed.NW.2 = imputed.NW.forb %>%
  dplyr::select(cover.change,leafN.bt:MAP.bt) %>%
  drop_na()

forb.effects = conditional_effects(forb.traits.NW.model)

#### Forb leafN ####

leafN.effects = forb.effects$leafN.final
leafN.effects$leafN.bt= (leafN.effects$leafN.final*7.150147) + 22.84208

attr(forb.imputed.NW.2$leafN.final, "scaled:scale")
attr(forb.imputed.NW.2$leafN.final, "scaled:center")

x.value = c(-3,-2,-1,0,1,2,3)
(x.value*6.403163) + 21.66854

forb.leafN.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = leafN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.bt, y = estimate__), color = "#F17236", size = 1.5) +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Leaf N (mg/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.68951,44.19001)+
  theme_classic()
forb.leafN.plot

ggsave("./Plots/forb.traits.NW.leafN.pdf", height = 3, width = 3)

#### Forb height ####

height.effects = forb.effects$height.final
height.effects$height.bt= (height.effects$height.final*0.2405401) + 0.3256687

attr(forb.imputed.NW.2$height.final, "scaled:scale")
attr(forb.imputed.NW.2$height.final, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*0.2396416) + 0.3214792

forb.height.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = height.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.bt, y = estimate__), color = "#F17236", size = 1.5) +  
  geom_ribbon(data = height.effects, aes(x = height.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Height (m)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.009000064,1.322333169)+
  theme_classic()
forb.height.plot

ggsave("./Plots/forb.traits.NW.height.pdf", height = 3, width = 3)

#### Forb rootN ####

rootN.effects = forb.effects$rootN.final
rootN.effects$rootN.bt= (rootN.effects$rootN.final*4.715357) + 11.3382

attr(forb.imputed.NW.2$rootN.final, "scaled:scale")
attr(forb.imputed.NW.2$rootN.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*4.173927) + 10.5821

forb.rootN.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = rootN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Root N (mg/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.3800000,25.6283268) +
  theme_classic()
forb.rootN.plot

ggsave("./Plots/forb.traits.NW.rootN.pdf", height = 3, width = 3)

#### Forb SLA ####

SLA.effects = forb.effects$SLA.final
SLA.effects$SLA.bt= (SLA.effects$SLA.final*8.061883) + 20.18927

attr(forb.imputed.NW.2$SLA.final, "scaled:scale")
attr(forb.imputed.NW.2$SLA.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*8.362926) + 20.14641

forb.SLA.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = SLA.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = expression("Specific Leaf Area (m"^2*"/kg)"), y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.617596,46.320001)+
  theme_classic()
forb.SLA.plot

ggsave("./Plots/forb.traits.NW.SLA.pdf", height = 3, width = 3)

#### Forb Depth ####

root.depth.effects = forb.effects$root.depth.final
root.depth.effects$Depth.bt= (root.depth.effects$root.depth.final*0.600499) + 0.6551094

attr(forb.imputed.NW.2$root.depth.final, "scaled:scale")
attr(forb.imputed.NW.2$root.depth.final, "scaled:center")

x.value = c(-1,0,1,2,3,4,5)
(x.value*0.6095345) + 0.6503508

forb.root.depth.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = Depth.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = Depth.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = Depth.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Rooting Depth (m)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.0479999,3.709335)+
  theme_classic()
forb.root.depth.plot

ggsave("./Plots/forb.traits.NW.root.depth.pdf", height = 3, width = 3)

#### Forb Diam ####

rootDiam.effects = forb.effects$rootDiam.final
rootDiam.effects$Diam.bt= (rootDiam.effects$rootDiam.final*0.152751) + 0.3553141

attr(forb.imputed.NW.2$rootDiam.final, "scaled:scale")
attr(forb.imputed.NW.2$rootDiam.final, "scaled:center")

x.value = c(-1,0,1,2,3,4,5)
(x.value*0.1590925) + 0.3579696

forb.rootDiam.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = Diam.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = Diam.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = Diam.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Root Diameter (mm)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.07840002,1.12601235)+
  theme_classic()
forb.rootDiam.plot

ggsave("./Plots/forb.traits.NW.rootDiam.pdf", height = 3, width = 3)

#### Forb SRL ####

SRL.effects = forb.effects$SRL.final
SRL.effects$SRL.bt= (SRL.effects$SRL.final*70.41695) + 103.2147

attr(forb.imputed.NW.2$SRL.final, "scaled:scale")
attr(forb.imputed.NW.2$SRL.final, "scaled:center")

x.value = c(-1,0,1,2,3,4)
(x.value*71.25811) + 105.5323

forb.SRL.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = SRL.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.bt, y = estimate__), color = "#F17236", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = "Specific Root Length (m/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.497203,417.749988)+
  theme_classic()
forb.SRL.plot

ggsave("./Plots/forb.traits.NW.SRL.pdf", height = 3, width = 3)

#### Forb RTD ####

RTD.effects = forb.effects$RTD.final
RTD.effects$RTD.bt= (RTD.effects$RTD.final*0.1178487) + 0.2169321

attr(forb.imputed.NW.2$RTD.final, "scaled:scale")
attr(forb.imputed.NW.2$RTD.final, "scaled:center")

x.value = c(-1,0,1,2,3)
(x.value*0.1184977) + 0.227384

forb.RTD.plot = ggplot() +
  geom_point(data = forb.imputed.NW.2, aes(x = RTD.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.bt, y = estimate__), color = "#F17236", size = 1.5) +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#F17236") +  
  labs(x = expression("Root Tissue Density (g/cm"^3*")"), y = "Cover Change (%)") + 
  ylim(-20.66667,21.24500)+
  #xlim(0.02325658,0.61999993)+
  theme_classic()
forb.RTD.plot

ggsave("./Plots/forb.traits.NW.RTD.pdf", height = 3, width = 3)

#### Forb RMF ####

RMF.effects = forb.effects$RMF.final
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

#### FORB interactions ####

forb.imputed.NW.traits.RTD.SRL = readRDS("./Results/forb.legume.imputed.traits.NW.RTD.SRL.rds")

RTD.SRL.effect = conditional_effects(forb.imputed.NW.traits.RTD.SRL, effects = "RTD.final:SRL.final")$`RTD.final:SRL.final`
RTD.SRL.effect$RTD.bt= (RTD.SRL.effect$RTD.final*0.1178487) + 0.2169321

# only plot high and low values

RTD.SRL.effect.2 = RTD.SRL.effect %>%
  dplyr::filter(effect2__ %in% c(-1,1.04))

# RTD
x.value = c(-1,0,1,2,3)
(x.value*0.1184977) + 0.227384

# SRL
x.value = c(-1,1.04)
(x.value*70.41695) + 103.2147

forb.RTD.x.SRL.plot = ggplot(data = RTD.SRL.effect.2, aes(x = RTD.bt, y = estimate__, group = as.factor(effect2__))) +
  geom_line(aes(color = as.factor(effect2__)), size = 1.5,show.legend = FALSE) +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = as.factor(effect2__)), alpha = 0.5,show.legend = FALSE) +  
  geom_line(aes(y = upper__, color = factor(effect2__)), size = 1.2) +
  geom_line(aes(y = lower__, color = factor(effect2__)), size = 1.2) + 
  labs(x = expression("Root Tissue Density (g/cm"^3*")"), y = "Cover Change (%)", color = "SRL (m/g)") +
  scale_colour_manual(values = c("black", "#F17236"), labels = c("176.45","32.79"))+
  scale_fill_manual(values = c("black", "#F17236"))+
  #xlim(0,0.59)+
  theme_classic()
forb.RTD.x.SRL.plot

ggsave("./Plots/forb.traits.NW.RTD.SRL.pdf", height = 3, width = 3)
ggsave("./Plots/forb.traits.NW.RTD.SRL.legend.pdf", height = 3, width = 3)

#### FORB environment interactions ####

forb.imputed.NW.traits.leafN.DSI = readRDS("./Results/imputed.NW.forb.legume.traits.leafN.MAP.DSI.rds")

leafN.DSI.effect = conditional_effects(forb.imputed.NW.traits.leafN.DSI, effects = "leafN.final:mean.DSI")$`leafN.final:mean.DSI`
leafN.DSI.effect$leafN.bt= (leafN.DSI.effect$leafN.final*7.150147) + 22.84208

# only plot high and low values

leafN.DSI.effect.2 = leafN.DSI.effect %>%
  dplyr::filter(effect2__ %in% c(-1,1))

# leafN
x.value = c(-3,0,1,2,3)
(x.value*7.150147) + 22.84208

# DSI
x.value = c(-1,1)
(x.value*0.1957044) + -0.4647529

forb.leafN.x.DSI.plot = ggplot(data = leafN.DSI.effect.2, aes(x = leafN.bt, y = estimate__, group = as.factor(effect2__))) +
  geom_line(aes(color = as.factor(effect2__)), size = 1.5,show.legend = FALSE) +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = as.factor(effect2__)), alpha = 0.5,show.legend = FALSE) +  
  geom_line(aes(y = upper__, color = factor(effect2__)), size = 1.2) +
  geom_line(aes(y = lower__, color = factor(effect2__)), size = 1.2) + 
  labs(x = "Leaf N (mg/g)", y = "Cover Change (%)", color = "DSI") +
  scale_colour_manual(values = c("black", "#F17236"), labels = c("-0.27","-0.66"))+
  scale_fill_manual(values = c("black", "#F17236"))+
  #xlim(2.68951,44.19001)+
  theme_classic()
forb.leafN.x.DSI.plot

ggsave("./Plots/forb.traits.NW.leafN.DSI.pdf", height = 3, width = 3)
ggsave("./Plots/forb.traits.NW.leafN.DSI.legend.pdf", height = 3, width = 3)


#### Graminoid Backtransform ####
enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)

imputed.NW.graminoid = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.graminoid.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)

graminoid.traits.NW.model = readRDS("./Results/graminoid.imputed.traits.no_woody.rds")

imputed.NW.graminoid$leafN.bt =(imputed.NW.graminoid$leafN.final*6.11409) + 18.36715
imputed.NW.graminoid$height.bt = (imputed.NW.graminoid$height.final*0.2220977) + 0.3845679
imputed.NW.graminoid$rootN.bt = (imputed.NW.graminoid$rootN.final*2.612803) + 7.369866
imputed.NW.graminoid$SLA.bt = (imputed.NW.graminoid$SLA.final*9.389227) + 19.51921
imputed.NW.graminoid$Depth.bt = (imputed.NW.graminoid$root.depth.final*0.3623682) + 0.4500424
imputed.NW.graminoid$Diam.bt = (imputed.NW.graminoid$rootDiam.final*0.1710708) + 0.3415072
imputed.NW.graminoid$SRL.bt = (imputed.NW.graminoid$SRL.final*80.14288) + 111.7918
imputed.NW.graminoid$RTD.bt = (imputed.NW.graminoid$RTD.final*0.1062242) + 0.2651114
imputed.NW.graminoid$RMF.bt = (imputed.NW.graminoid$RMF.final*0.08416183) + 0.4159471
imputed.NW.graminoid$DSI.bt = (imputed.NW.graminoid$mean.DSI*0.1941463) + -0.4534765
imputed.NW.graminoid$MAP.bt = (imputed.NW.graminoid$mean.MAP*561.8211) + 735.9487

graminoid.imputed.NW.2 = imputed.NW.graminoid %>%
  dplyr::select(cover.change,leafN.bt:MAP.bt) %>%
  drop_na()

graminoid.effects = conditional_effects(graminoid.traits.NW.model)

#### graminoid leafN ####

leafN.effects = graminoid.effects$leafN.final
leafN.effects$leafN.bt= (leafN.effects$leafN.final*6.11409) + 18.36715

attr(graminoid.imputed.NW.2$leafN.final, "scaled:scale")
attr(graminoid.imputed.NW.2$leafN.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3)
(x.value*6.281442) + 18.41922

graminoid.leafN.plot = ggplot() +
  geom_point(data = graminoid.imputed.NW.2, aes(x = leafN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.bt, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Leaf N (mg/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.68951,44.19001)+
  theme_classic()
graminoid.leafN.plot

ggsave("./Plots/graminoid.traits.NW.leafN.pdf", height = 3, width = 3)

#### graminoid height ####

height.effects = graminoid.effects$height.final
height.effects$height.bt= (height.effects$height.final*0.2220977) + 0.3845679

attr(graminoid.imputed.NW.2$height.final, "scaled:scale")
attr(graminoid.imputed.NW.2$height.final, "scaled:center")

x.value = c(-1,0,1,2,3)
(x.value*0.2266876) + 0.3988779

graminoid.height.plot = ggplot() +
  geom_point(data = graminoid.imputed.NW.2, aes(x = height.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.bt, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = height.effects, aes(x = height.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Height (m)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.009000064,1.322333169)+
  theme_classic() +
  theme(plot.margin = margin(t = 5, r = 15, b = 5, l = 5))
graminoid.height.plot

ggsave("./Plots/graminoid.traits.NW.height.pdf", height = 3, width = 3)

#### graminoid rootN ####

rootN.effects = graminoid.effects$rootN.final
rootN.effects$rootN.bt= (rootN.effects$rootN.final*2.612803) + 7.369866

attr(graminoid.imputed.NW.2$rootN.final, "scaled:scale")
attr(graminoid.imputed.NW.2$rootN.final, "scaled:center")

x.value = c(-2,-1,0,1,2,3,4)
(x.value*2.712129) + 7.59045

graminoid.rootN.plot = ggplot() +
  geom_point(data = graminoid.imputed.NW.2, aes(x = rootN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.bt, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Root N (mg/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.3800000,25.6283268)+
  theme_classic()
graminoid.rootN.plot

ggsave("./Plots/graminoid.traits.NW.rootN.pdf", height = 3, width = 3)

#### graminoid SLA ####

SLA.effects = graminoid.effects$SLA.final
SLA.effects$SLA.bt= (SLA.effects$SLA.final*9.389227) + 19.51921

attr(graminoid.imputed.NW.2$SLA.final, "scaled:scale")
attr(graminoid.imputed.NW.2$SLA.final, "scaled:center")

x.value = c(-1,0,1,2,3)
(x.value*9.570923) + 20.31763

graminoid.SLA.plot = ggplot() +
  geom_point(data = graminoid.imputed.NW.2, aes(x = SLA.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.bt, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = expression("Specific Leaf Area (m"^2*"/kg)"), y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.617596,46.320001)+
  theme_classic()
graminoid.SLA.plot

ggsave("./Plots/graminoid.traits.NW.SLA.pdf", height = 3, width = 3)

#### graminoid Depth ####

root.depth.effects = graminoid.effects$root.depth.final
root.depth.effects$Depth.bt= (root.depth.effects$root.depth.final*0.3623682) + 0.4500424

graminoid.root.depth.plot = ggplot() +
  geom_point(data = graminoid.imputed.NW.2, aes(x = Depth.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = Depth.bt, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = Depth.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Rooting Depth (m)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.0479999,3.709335)+
  theme_classic()
graminoid.root.depth.plot

ggsave("./Plots/graminoid.traits.NW.root.depth.pdf", height = 3, width = 3)

#### graminoid Diam ####

rootDiam.effects = graminoid.effects$rootDiam.final
rootDiam.effects$Diam.bt= (rootDiam.effects$rootDiam.final*0.1710708) + 0.3415072

graminoid.rootDiam.plot = ggplot() +
  geom_point(data = graminoid.imputed.NW.2, aes(x = Diam.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = Diam.bt, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = Diam.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Root Diameter (mm)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.07840002,1.12601235)+
  theme_classic()
graminoid.rootDiam.plot

ggsave("./Plots/graminoid.traits.NW.rootDiam.pdf", height = 3, width = 3)

#### graminoid SRL ####

SRL.effects = graminoid.effects$SRL.final
SRL.effects$SRL.bt= (SRL.effects$SRL.final*80.14288) + 111.7918

graminoid.SRL.plot = ggplot() +
  geom_point(data = graminoid.imputed.NW.2, aes(x = SRL.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.bt, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Specific Root Length (m/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.497203,417.749988)+
  theme_classic()
graminoid.SRL.plot

ggsave("./Plots/graminoid.traits.NW.SRL.pdf", height = 3, width = 3)

#### graminoid RTD ####

RTD.effects = graminoid.effects$RTD.final
RTD.effects$RTD.bt= (RTD.effects$RTD.final*0.1062242) + 0.2651114

graminoid.RTD.plot = ggplot() +
  geom_point(data = graminoid.imputed.NW.2, aes(x = RTD.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.bt, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = expression("Root Tissue Density (g/cm"^3*")"), y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.02325658,0.61999993)+
  theme_classic()
graminoid.RTD.plot

ggsave("./Plots/graminoid.traits.NW.RTD.pdf", height = 3, width = 3)

#### graminoid RMF ####

RMF.effects = graminoid.effects$RMF.final
RMF.effects$RMF.bt= (RMF.effects$RMF.final*0.08416183) + 0.4159471

graminoid.RMF.plot = ggplot() +
  geom_point(data = graminoid.imputed.NW.2, aes(x = RMF.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RMF.effects, aes(x = RMF.bt, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RMF.effects, aes(x = RMF.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Root Mass Fraction (g/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.1211193,0.7328659)+
  theme_classic()
graminoid.RMF.plot

ggsave("./Plots/graminoid.traits.NW.RMF.pdf", height = 3, width = 3)

#### graminoid DSI ####

DSI.effects = graminoid.effects$mean.DSI
DSI.effects$DSI.bt= (DSI.effects$mean.DSI*0.1941463) + -0.4534765

graminoid.DSI.plot = ggplot() +
  geom_point(data = graminoid.imputed.NW.2, aes(x = DSI.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = DSI.effects, aes(x = DSI.bt, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = DSI.effects, aes(x = DSI.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Drought Severity Index", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(-0.8722616,0.1242846)+
  theme_classic()
graminoid.DSI.plot

ggsave("./Plots/graminoid.traits.NW.DSI.pdf", height = 3, width = 3)

#### graminoid MAP ####

MAP.effects = graminoid.effects$mean.MAP
MAP.effects$MAP.bt= (MAP.effects$mean.MAP*561.8211) + 735.9487

graminoid.MAP.plot = ggplot() +
  geom_point(data = graminoid.imputed.NW.2, aes(x = MAP.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = MAP.effects, aes(x = MAP.bt, y = estimate__), color = "#6E687E", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = MAP.effects, aes(x = MAP.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6E687E") +  
  labs(x = "Mean Annual Precipitation (mm)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(132.7,2366.1)+
  theme_classic()
graminoid.MAP.plot

ggsave("./Plots/graminoid.traits.NW.MAP.pdf", height = 3, width = 3)

#### PERENNIAL graminoid Backtransform ####

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)
imputed.NW.perennial.graminoid = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennial.graminoid.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)

perennial.graminoid.traits.NW.model = readRDS("./Results/perennial.graminoid.imputed.traits.no_woody.rds")

imputed.NW.perennial.graminoid$leafN.bt =(imputed.NW.perennial.graminoid$leafN.final*5.368206) + 17.71867
imputed.NW.perennial.graminoid$height.bt = (imputed.NW.perennial.graminoid$height.final*0.2272111) + 0.38957
imputed.NW.perennial.graminoid$rootN.bt = (imputed.NW.perennial.graminoid$rootN.final*2.803468) + 7.483601
imputed.NW.perennial.graminoid$SLA.bt = (imputed.NW.perennial.graminoid$SLA.final*8.951075) + 18.54179
imputed.NW.perennial.graminoid$Depth.bt = (imputed.NW.perennial.graminoid$root.depth.final*0.3617313) + 0.4456165
imputed.NW.perennial.graminoid$Diam.bt = (imputed.NW.perennial.graminoid$rootDiam.final*0.1696134) + 0.3351268
imputed.NW.perennial.graminoid$SRL.bt = (imputed.NW.perennial.graminoid$SRL.final*81.21342) + 108.6274
imputed.NW.perennial.graminoid$RTD.bt = (imputed.NW.perennial.graminoid$RTD.final*0.1008478) + 0.2766672
imputed.NW.perennial.graminoid$RMF.bt = (imputed.NW.perennial.graminoid$RMF.final*0.08710869) + 0.4175176
imputed.NW.perennial.graminoid$DSI.bt = (imputed.NW.perennial.graminoid$mean.DSI*0.1788675) + -0.4520457
imputed.NW.perennial.graminoid$MAP.bt = (imputed.NW.perennial.graminoid$mean.MAP*591.1718) + 800.2379

perennial.graminoid.imputed.NW.2 = imputed.NW.perennial.graminoid %>%
  dplyr::select(cover.change,leafN.bt:MAP.bt) %>%
  drop_na()

perennial.graminoid.effects = conditional_effects(perennial.graminoid.traits.NW.model)

#### Perennial graminoid leafN ####

leafN.effects = perennial.graminoid.effects$leafN.final
leafN.effects$leafN.bt= (leafN.effects$leafN.final*5.368206) + 17.71867

perennial.graminoid.leafN.plot = ggplot() +
  geom_point(data = perennial.graminoid.imputed.NW.2, aes(x = leafN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.bt, y = estimate__), color = "#6089B5", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6089B5") +  
  labs(x = "Leaf N (mg/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.68951,44.19001)+
  theme_classic()
perennial.graminoid.leafN.plot

ggsave("./Plots/perennial.graminoid.traits.NW.leafN.pdf", height = 3, width = 3)

#### Perennial graminoid height ####

height.effects = perennial.graminoid.effects$height.final
height.effects$height.bt= (height.effects$height.final*0.2272111) + 0.38957

perennial.graminoid.height.plot = ggplot() +
  geom_point(data = perennial.graminoid.imputed.NW.2, aes(x = height.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.bt, y = estimate__), color = "#6089B5", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = height.effects, aes(x = height.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6089B5") +  
  labs(x = "Height (m)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.009000064,1.322333169)+
  theme_classic() +
  theme(plot.margin = margin(t = 5, r = 15, b = 5, l = 5))
perennial.graminoid.height.plot

ggsave("./Plots/perennial.graminoid.traits.NW.height.pdf", height = 3, width = 3)

#### Perennial graminoid rootN ####

rootN.effects = perennial.graminoid.effects$rootN.final
rootN.effects$rootN.bt= (rootN.effects$rootN.final*2.803468) + 7.483601

perennial.graminoid.rootN.plot = ggplot() +
  geom_point(data = perennial.graminoid.imputed.NW.2, aes(x = rootN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.bt, y = estimate__), color = "#6089B5", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6089B5") +  
  labs(x = "Root N (mg/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.3800000,25.6283268)+
  theme_classic()
perennial.graminoid.rootN.plot

ggsave("./Plots/perennial.graminoid.traits.NW.rootN.pdf", height = 3, width = 3)

#### Perennial graminoid SLA ####

SLA.effects = perennial.graminoid.effects$SLA.final
SLA.effects$SLA.bt= (SLA.effects$SLA.final*8.951075) + 18.54179

perennial.graminoid.SLA.plot = ggplot() +
  geom_point(data = perennial.graminoid.imputed.NW.2, aes(x = SLA.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.bt, y = estimate__), color = "#6089B5", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6089B5") +  
  labs(x = expression("Specific Leaf Area (m"^2*"/kg)"), y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.617596,46.320001)+
  theme_classic()
perennial.graminoid.SLA.plot

ggsave("./Plots/perennial.graminoid.traits.NW.SLA.pdf", height = 3, width = 3)

#### Perennial graminoid depth ####

root.depth.effects = perennial.graminoid.effects$root.depth.final
root.depth.effects$Depth.bt= (root.depth.effects$root.depth.final*0.3617313) + 0.4456165

perennial.graminoid.root.depth.plot = ggplot() +
  geom_point(data = perennial.graminoid.imputed.NW.2, aes(x = Depth.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = Depth.bt, y = estimate__), color = "#6089B5", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = Depth.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6089B5") +  
  labs(x = "Rooting Depth (m)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.0479999,3.709335)+
  theme_classic()
perennial.graminoid.root.depth.plot

ggsave("./Plots/perennial.graminoid.traits.NW.root.depth.pdf", height = 3, width = 3)

#### Perennial graminoid Diam ####

rootDiam.effects = perennial.graminoid.effects$rootDiam.final
rootDiam.effects$Diam.bt= (rootDiam.effects$rootDiam.final*0.1696134) + 0.3351268

perennial.graminoid.rootDiam.plot = ggplot() +
  geom_point(data = perennial.graminoid.imputed.NW.2, aes(x = Diam.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = Diam.bt, y = estimate__), color = "#6089B5", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = Diam.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6089B5") +  
  labs(x = "Root Diameter (mm)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.07840002,1.12601235)+
  theme_classic()
perennial.graminoid.rootDiam.plot

ggsave("./Plots/perennial.graminoid.traits.NW.rootDiam.pdf", height = 3, width = 3)

#### Perennial graminoid SRL ####

SRL.effects = perennial.graminoid.effects$SRL.final
SRL.effects$SRL.bt= (SRL.effects$SRL.final*81.21342) + 108.6274

perennial.graminoid.SRL.plot = ggplot() +
  geom_point(data = perennial.graminoid.imputed.NW.2, aes(x = SRL.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.bt, y = estimate__), color = "#6089B5", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6089B5") +  
  labs(x = "Specific Root Length (m/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.497203,417.749988)+
  theme_classic()
perennial.graminoid.SRL.plot

ggsave("./Plots/perennial.graminoid.traits.NW.SRL.pdf", height = 3, width = 3)

#### Perennial graminoid RTD ####

RTD.effects = perennial.graminoid.effects$RTD.final
RTD.effects$RTD.bt= (RTD.effects$RTD.final*0.1008478) + 0.2766672

perennial.graminoid.RTD.plot = ggplot() +
  geom_point(data = perennial.graminoid.imputed.NW.2, aes(x = RTD.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.bt, y = estimate__), color = "#6089B5", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6089B5") +  
  labs(x = expression("Root Tissue Density (g/cm"^3*")"), y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.02325658,0.61999993)+
  theme_classic()
perennial.graminoid.RTD.plot

ggsave("./Plots/perennial.graminoid.traits.NW.RTD.pdf", height = 3, width = 3)

#### Perennial graminoid RMF ####

RMF.effects = perennial.graminoid.effects$RMF.final
RMF.effects$RMF.bt= (RMF.effects$RMF.final*0.08710869) + 0.4175176

perennial.graminoid.RMF.plot = ggplot() +
  geom_point(data = perennial.graminoid.imputed.NW.2, aes(x = RMF.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RMF.effects, aes(x = RMF.bt, y = estimate__), color = "#6089B5", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RMF.effects, aes(x = RMF.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6089B5") +  
  labs(x = "Root Mass Fraction (g/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.1211193,0.7328659)+
  theme_classic()
perennial.graminoid.RMF.plot

ggsave("./Plots/perennial.graminoid.traits.NW.RMF.pdf", height = 3, width = 3)

#### Perennial graminoid DSI ####

DSI.effects = perennial.graminoid.effects$mean.DSI
DSI.effects$DSI.bt= (DSI.effects$mean.DSI*0.1788675) + -0.4520457

perennial.graminoid.DSI.plot = ggplot() +
  geom_point(data = perennial.graminoid.imputed.NW.2, aes(x = DSI.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = DSI.effects, aes(x = DSI.bt, y = estimate__), color = "#6089B5", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = DSI.effects, aes(x = DSI.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6089B5") +  
  labs(x = "Drought Severity Index", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(-0.8722616,0.1242846)+
  theme_classic()
perennial.graminoid.DSI.plot

ggsave("./Plots/perennial.graminoid.traits.NW.DSI.pdf", height = 3, width = 3)

#### Perennial graminoid MAP ####

MAP.effects = perennial.graminoid.effects$mean.MAP
MAP.effects$MAP.bt= (MAP.effects$mean.MAP*591.1718) + 800.2379

perennial.graminoid.MAP.plot = ggplot() +
  geom_point(data = perennial.graminoid.imputed.NW.2, aes(x = MAP.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = MAP.effects, aes(x = MAP.bt, y = estimate__), color = "#6089B5", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = MAP.effects, aes(x = MAP.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6089B5") +  
  labs(x = "Mean Annual Precipitation (mm)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(132.7,2366.1)+
  theme_classic()
perennial.graminoid.MAP.plot

ggsave("./Plots/perennial.graminoid.traits.NW.MAP.pdf", height = 3, width = 3)

#### PERENNIAL forb Backtransform ####

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)
imputed.NW.perennial.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.perennial.forb.legume.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)

perennial.forb.traits.NW.model = readRDS("./Results/perennial.forb.imputed.traits.no_woody.rds")

imputed.NW.perennial.forb$leafN.bt =(imputed.NW.perennial.forb$leafN.final*7.053831) + 22.60696
imputed.NW.perennial.forb$height.bt = (imputed.NW.perennial.forb$height.final*0.241945) + 0.3419011
imputed.NW.perennial.forb$rootN.bt = (imputed.NW.perennial.forb$rootN.final*4.822226) + 11.53761
imputed.NW.perennial.forb$SLA.bt = (imputed.NW.perennial.forb$SLA.final*7.12074) + 18.92746
imputed.NW.perennial.forb$Depth.bt = (imputed.NW.perennial.forb$root.depth.final*0.6701013) + 0.7026712
imputed.NW.perennial.forb$Diam.bt = (imputed.NW.perennial.forb$rootDiam.final*0.1625762) + 0.367839
imputed.NW.perennial.forb$SRL.bt = (imputed.NW.perennial.forb$SRL.final*65.21507) + 91.01407
imputed.NW.perennial.forb$RTD.bt = (imputed.NW.perennial.forb$RTD.final*0.1212438) + 0.2290445
imputed.NW.perennial.forb$RMF.bt = (imputed.NW.perennial.forb$RMF.final*0.1166379) + 0.3939654
imputed.NW.perennial.forb$DSI.bt = (imputed.NW.perennial.forb$mean.DSI*0.1774911) + -0.4569303
imputed.NW.perennial.forb$MAP.bt = (imputed.NW.perennial.forb$mean.MAP*460.4944) + 749.6066

perennial.forb.imputed.NW.2 = imputed.NW.perennial.forb %>%
  dplyr::select(cover.change,leafN.bt:MAP.bt) %>%
  drop_na()

perennial.forb.effects = conditional_effects(perennial.forb.traits.NW.model)

#### Perennial forb leafN ####

leafN.effects = perennial.forb.effects$leafN.final
leafN.effects$leafN.bt= (leafN.effects$leafN.final*7.053831) + 22.60696

perennial.forb.leafN.plot = ggplot() +
  geom_point(data = perennial.forb.imputed.NW.2, aes(x = leafN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.bt, y = estimate__), color = "black", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "black") +  
  labs(x = "Leaf N (mg/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.68951,44.19001)+
  theme_classic()
perennial.forb.leafN.plot

ggsave("./Plots/perennial.forb.traits.NW.leafN.pdf", height = 3, width = 3)

#### Perennial forb height ####

height.effects = perennial.forb.effects$height.final
height.effects$height.bt= (height.effects$height.final*0.241945) + 0.3419011

perennial.forb.height.plot = ggplot() +
  geom_point(data = perennial.forb.imputed.NW.2, aes(x = height.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.bt, y = estimate__), color = "black", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = height.effects, aes(x = height.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "black") +  
  labs(x = "Height (m)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.009000064,1.322333169)+
  theme_classic()
perennial.forb.height.plot

ggsave("./Plots/perennial.forb.traits.NW.height.pdf", height = 3, width = 3)

#### Perennial forb rootN ####

rootN.effects = perennial.forb.effects$rootN.final
rootN.effects$rootN.bt= (rootN.effects$rootN.final*4.822226) + 11.53761

perennial.forb.rootN.plot = ggplot() +
  geom_point(data = perennial.forb.imputed.NW.2, aes(x = rootN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.bt, y = estimate__), color = "black", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "black") +  
  labs(x = "Root N (mg/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.3800000,25.6283268)+
  theme_classic()
perennial.forb.rootN.plot

ggsave("./Plots/perennial.forb.traits.NW.rootN.pdf", height = 3, width = 3)

#### Perennial forb SLA ####

SLA.effects = perennial.forb.effects$SLA.final
SLA.effects$SLA.bt= (SLA.effects$SLA.final*7.12074) + 18.92746

perennial.forb.SLA.plot = ggplot() +
  geom_point(data = perennial.forb.imputed.NW.2, aes(x = SLA.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.bt, y = estimate__), color = "black", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "black") +  
  labs(x = expression("Specific Leaf Area (m"^2*"/kg)"), y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.617596,46.320001)+
  theme_classic()
perennial.forb.SLA.plot

ggsave("./Plots/perennial.forb.traits.NW.SLA.pdf", height = 3, width = 3)

#### Perennial forb depth ####

root.depth.effects = perennial.forb.effects$root.depth.final
root.depth.effects$Depth.bt= (root.depth.effects$root.depth.final*0.6701013) + 0.7026712

perennial.forb.root.depth.plot = ggplot() +
  geom_point(data = perennial.forb.imputed.NW.2, aes(x = Depth.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = Depth.bt, y = estimate__), color = "black", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = Depth.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "black") +  
  labs(x = "Rooting Depth (m)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.0479999,3.709335)+
  theme_classic()
perennial.forb.root.depth.plot

ggsave("./Plots/perennial.forb.traits.NW.root.depth.pdf", height = 3, width = 3)

#### Perennial forb Diam ####

rootDiam.effects = perennial.forb.effects$rootDiam.final
rootDiam.effects$Diam.bt= (rootDiam.effects$rootDiam.final*0.1625762) + 0.367839

perennial.forb.rootDiam.plot = ggplot() +
  geom_point(data = perennial.forb.imputed.NW.2, aes(x = Diam.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = Diam.bt, y = estimate__), color = "black", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = Diam.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "black") +  
  labs(x = "Root Diameter (mm)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.07840002,1.12601235)+
  theme_classic()
perennial.forb.rootDiam.plot

ggsave("./Plots/perennial.forb.traits.NW.rootDiam.pdf", height = 3, width = 3)

#### Perennial forb SRL ####

SRL.effects = perennial.forb.effects$SRL.final
SRL.effects$SRL.bt= (SRL.effects$SRL.final*65.21507) + 91.01407

perennial.forb.SRL.plot = ggplot() +
  geom_point(data = perennial.forb.imputed.NW.2, aes(x = SRL.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.bt, y = estimate__), color = "black", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "black") +  
  labs(x = "Specific Root Length (m/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.497203,417.749988)+
  theme_classic()
perennial.forb.SRL.plot

ggsave("./Plots/perennial.forb.traits.NW.SRL.pdf", height = 3, width = 3)

#### Perennial forb RTD ####

RTD.effects = perennial.forb.effects$RTD.final
RTD.effects$RTD.bt= (RTD.effects$RTD.final*0.1212438) + 0.2290445

perennial.forb.RTD.plot = ggplot() +
  geom_point(data = perennial.forb.imputed.NW.2, aes(x = RTD.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.bt, y = estimate__), color = "black", size = 1.5) +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "black") +  
  labs(x = expression("Root Tissue Density (g/cm"^3*")"), y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.02325658,0.61999993)+
  theme_classic()
perennial.forb.RTD.plot

ggsave("./Plots/perennial.forb.traits.NW.RTD.pdf", height = 3, width = 3)

#### Perennial forb RMF ####

RMF.effects = perennial.forb.effects$RMF.final
RMF.effects$RMF.bt= (RMF.effects$RMF.final*0.1166379) + 0.3939654

perennial.forb.RMF.plot = ggplot() +
  geom_point(data = perennial.forb.imputed.NW.2, aes(x = RMF.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RMF.effects, aes(x = RMF.bt, y = estimate__), color = "black", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RMF.effects, aes(x = RMF.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "black") +  
  labs(x = "Root Mass Fraction (g/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.1211193,0.7328659)+
  theme_classic()
perennial.forb.RMF.plot

ggsave("./Plots/perennial.forb.traits.NW.RMF.pdf", height = 3, width = 3)

#### Perennial forb DSI ####

DSI.effects = perennial.forb.effects$mean.DSI
DSI.effects$DSI.bt= (DSI.effects$mean.DSI*0.1774911) + -0.4569303

perennial.forb.DSI.plot = ggplot() +
  geom_point(data = perennial.forb.imputed.NW.2, aes(x = DSI.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = DSI.effects, aes(x = DSI.bt, y = estimate__), color = "black", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = DSI.effects, aes(x = DSI.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "black") +  
  labs(x = "Drought Severity Index", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(-0.8722616,0.1242846)+
  theme_classic()
perennial.forb.DSI.plot

ggsave("./Plots/perennial.forb.traits.NW.DSI.pdf", height = 3, width = 3)

#### Perennial forb MAP ####

MAP.effects = perennial.forb.effects$mean.MAP
MAP.effects$MAP.bt= (MAP.effects$mean.MAP*460.4944) + 749.6066

perennial.forb.MAP.plot = ggplot() +
  geom_point(data = perennial.forb.imputed.NW.2, aes(x = MAP.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = MAP.effects, aes(x = MAP.bt, y = estimate__), color = "black", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = MAP.effects, aes(x = MAP.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "black") +  
  labs(x = "Mean Annual Precipitation (mm)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(132.7,2366.1)+
  theme_classic()
perennial.forb.MAP.plot

ggsave("./Plots/perennial.forb.traits.NW.MAP.pdf", height = 3, width = 3)

#### ANNUAL FORB Backtransform ####

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)
imputed.NW.annual.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.NW.annual.forb.legume.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) %>%
  mutate_at(vars(26:36),scale)

annual.forb.traits.NW.model = readRDS("./Results/annual.forb.imputed.traits.no_woody.rds")

imputed.NW.annual.forb$leafN.bt =(imputed.NW.annual.forb$leafN.final*7.176246) + 22.88502
imputed.NW.annual.forb$height.bt = (imputed.NW.annual.forb$height.final*0.2374835) + 0.2992857
imputed.NW.annual.forb$rootN.bt = (imputed.NW.annual.forb$rootN.final*4.499487) + 10.98607
imputed.NW.annual.forb$SLA.bt = (imputed.NW.annual.forb$SLA.final*9.131592) + 22.54775
imputed.NW.annual.forb$Depth.bt = (imputed.NW.annual.forb$root.depth.final*0.4510935) + 0.5646051
imputed.NW.annual.forb$Diam.bt = (imputed.NW.annual.forb$rootDiam.final*0.1283094) + 0.3362157
imputed.NW.annual.forb$SRL.bt = (imputed.NW.annual.forb$SRL.final*73.50296) + 124.7039
imputed.NW.annual.forb$RTD.bt = (imputed.NW.annual.forb$RTD.final*0.1100227) + 0.1997673
imputed.NW.annual.forb$RMF.bt = (imputed.NW.annual.forb$RMF.final*0.1124093) + 0.3868801
imputed.NW.annual.forb$DSI.bt = (imputed.NW.annual.forb$mean.DSI*0.2225145) + -0.4758213
imputed.NW.annual.forb$MAP.bt = (imputed.NW.annual.forb$mean.MAP*303.7785) + 484.4501

annual.forb.imputed.NW.2 = imputed.NW.annual.forb %>%
  dplyr::select(cover.change,leafN.bt:MAP.bt) %>%
  drop_na()

annual.forb.effects = conditional_effects(annual.forb.traits.NW.model)

#### annual forb leafN ####

leafN.effects = annual.forb.effects$leafN.final
leafN.effects$leafN.bt= (leafN.effects$leafN.final*7.176246) + 22.88502

annual.forb.leafN.plot = ggplot() +
  geom_point(data = annual.forb.imputed.NW.2, aes(x = leafN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = leafN.effects, aes(x = leafN.bt, y = estimate__), color = "#B50200", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = leafN.effects, aes(x = leafN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#B50200") +  
  labs(x = "Leaf N (mg/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.68951,44.19001)+
  theme_classic()
annual.forb.leafN.plot

ggsave("./Plots/annual.forb.traits.NW.leafN.pdf", height = 3, width = 3)

#### annual forb height ####

height.effects = annual.forb.effects$height.final
height.effects$height.bt= (height.effects$height.final*0.2374835) + 0.2992857

annual.forb.height.plot = ggplot() +
  geom_point(data = annual.forb.imputed.NW.2, aes(x = height.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = height.effects, aes(x = height.bt, y = estimate__), color = "#B50200", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = height.effects, aes(x = height.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#B50200") +  
  labs(x = "Height (m)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.009000064,1.322333169)+
  theme_classic()
annual.forb.height.plot

ggsave("./Plots/annual.forb.traits.NW.height.pdf", height = 3, width = 3)

#### annual forb rootN ####

rootN.effects = annual.forb.effects$rootN.final
rootN.effects$rootN.bt= (rootN.effects$rootN.final*4.499487) + 10.98607

annual.forb.rootN.plot = ggplot() +
  geom_point(data = annual.forb.imputed.NW.2, aes(x = rootN.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootN.effects, aes(x = rootN.bt, y = estimate__), color = "#B50200", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootN.effects, aes(x = rootN.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#B50200") +  
  labs(x = "Root N (mg/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.3800000,25.6283268)+
  theme_classic()
annual.forb.rootN.plot

ggsave("./Plots/annual.forb.traits.NW.rootN.pdf", height = 3, width = 3)

#### annual forb SLA ####

SLA.effects = annual.forb.effects$SLA.final
SLA.effects$SLA.bt= (SLA.effects$SLA.final*9.131592) + 22.54775

annual.forb.SLA.plot = ggplot() +
  geom_point(data = annual.forb.imputed.NW.2, aes(x = SLA.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SLA.effects, aes(x = SLA.bt, y = estimate__), color = "#B50200", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SLA.effects, aes(x = SLA.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#B50200") +  
  labs(x = expression("Specific Leaf Area (m"^2*"/kg)"), y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.617596,46.320001)+
  theme_classic()
annual.forb.SLA.plot

ggsave("./Plots/annual.forb.traits.NW.SLA.pdf", height = 3, width = 3)

#### annual forb depth ####

root.depth.effects = annual.forb.effects$root.depth.final
root.depth.effects$Depth.bt= (root.depth.effects$root.depth.final*0.4510935) + 0.5646051

annual.forb.root.depth.plot = ggplot() +
  geom_point(data = annual.forb.imputed.NW.2, aes(x = Depth.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = root.depth.effects, aes(x = Depth.bt, y = estimate__), color = "#B50200", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = root.depth.effects, aes(x = Depth.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#B50200") +  
  labs(x = "Rooting Depth (m)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.0479999,3.709335)+
  theme_classic()
annual.forb.root.depth.plot

ggsave("./Plots/annual.forb.traits.NW.root.depth.pdf", height = 3, width = 3)

#### annual forb Diam ####

rootDiam.effects = annual.forb.effects$rootDiam.final
rootDiam.effects$Diam.bt= (rootDiam.effects$rootDiam.final*0.1283094) + 0.3362157

annual.forb.rootDiam.plot = ggplot() +
  geom_point(data = annual.forb.imputed.NW.2, aes(x = Diam.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = rootDiam.effects, aes(x = Diam.bt, y = estimate__), color = "#B50200", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = rootDiam.effects, aes(x = Diam.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#B50200") +  
  labs(x = "Root Diameter (mm)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.07840002,1.12601235)+
  theme_classic()
annual.forb.rootDiam.plot

ggsave("./Plots/annual.forb.traits.NW.rootDiam.pdf", height = 3, width = 3)

#### annual forb SRL ####

SRL.effects = annual.forb.effects$SRL.final
SRL.effects$SRL.bt= (SRL.effects$SRL.final*73.50296) + 124.7039

annual.forb.SRL.plot = ggplot() +
  geom_point(data = annual.forb.imputed.NW.2, aes(x = SRL.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = SRL.effects, aes(x = SRL.bt, y = estimate__), color = "#B50200", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = SRL.effects, aes(x = SRL.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#B50200") +  
  labs(x = "Specific Root Length (m/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(2.497203,417.749988)+
  theme_classic()
annual.forb.SRL.plot

ggsave("./Plots/annual.forb.traits.NW.SRL.pdf", height = 3, width = 3)

#### annual forb RTD ####

RTD.effects = annual.forb.effects$RTD.final
RTD.effects$RTD.bt= (RTD.effects$RTD.final*0.1100227) + 0.1997673

annual.forb.RTD.plot = ggplot() +
  geom_point(data = annual.forb.imputed.NW.2, aes(x = RTD.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RTD.effects, aes(x = RTD.bt, y = estimate__), color = "#B50200", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RTD.effects, aes(x = RTD.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#B50200") +  
  labs(x = expression("Root Tissue Density (g/cm"^3*")"), y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.02325658,0.61999993)+
  theme_classic()
annual.forb.RTD.plot

ggsave("./Plots/annual.forb.traits.NW.RTD.pdf", height = 3, width = 3)

#### annual forb RMF ####

RMF.effects = annual.forb.effects$RMF.final
RMF.effects$RMF.bt= (RMF.effects$RMF.final*0.1124093) + 0.3868801

annual.forb.RMF.plot = ggplot() +
  geom_point(data = annual.forb.imputed.NW.2, aes(x = RMF.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = RMF.effects, aes(x = RMF.bt, y = estimate__), color = "#B50200", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = RMF.effects, aes(x = RMF.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#B50200") +  
  labs(x = "Root Mass Fraction (g/g)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(0.1211193,0.7328659)+
  theme_classic()
annual.forb.RMF.plot

ggsave("./Plots/annual.forb.traits.NW.RMF.pdf", height = 3, width = 3)

#### annual forb DSI ####

DSI.effects = annual.forb.effects$mean.DSI
DSI.effects$DSI.bt= (DSI.effects$mean.DSI*0.2225145) + -0.4758213

annual.forb.DSI.plot = ggplot() +
  geom_point(data = annual.forb.imputed.NW.2, aes(x = DSI.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = DSI.effects, aes(x = DSI.bt, y = estimate__), color = "#B50200", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = DSI.effects, aes(x = DSI.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#B50200") +  
  labs(x = "Drought Severity Index", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(-0.8722616,0.1242846)+
  theme_classic()
annual.forb.DSI.plot

ggsave("./Plots/annual.forb.traits.NW.DSI.pdf", height = 3, width = 3)

#### annual forb MAP ####

MAP.effects = annual.forb.effects$mean.MAP
MAP.effects$MAP.bt= (MAP.effects$mean.MAP*303.7785) + 484.4501

annual.forb.MAP.plot = ggplot() +
  geom_point(data = annual.forb.imputed.NW.2, aes(x = MAP.bt, y = cover.change), color = "black", alpha = 0.5) +
  geom_line(data = MAP.effects, aes(x = MAP.bt, y = estimate__), color = "#B50200", size = 1.5, linetype = "dashed") +  
  geom_ribbon(data = MAP.effects, aes(x = MAP.bt, ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#B50200") +  
  labs(x = "Mean Annual Precipitation (mm)", y = "Cover Change (%)") +
  ylim(-20.66667,21.24500)+
  #xlim(132.7,2366.1)+
  theme_classic()
annual.forb.MAP.plot

ggsave("./Plots/annual.forb.traits.NW.MAP.pdf", height = 3, width = 3)

#### Dot Plots ####
imputed.traits.NW.model = readRDS("./Results/all.imputed.traits.no_woody.rds")
sum = summary(imputed.traits.NW.model)
coef = sum$fixed
row.names(coef) = c("Intercept","Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter",
                    "SRL", "RTD","RMF","DSI","Precipitation")
coef.2 = coef %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#769370"))
  
All.species.plot = ggplot(data = coef.2,
       aes(y = factor(row.names(coef.2), levels = rev(row.names(coef.2))),
           x = Estimate,
           xmin = `l-95% CI`,
           xmax = `u-95% CI`,
           fill = resp.fill)) +
  geom_pointrange(size = 1, shape = 21, color = "#769370") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "All-Species")
All.species.plot

annual.traits.NW.model = readRDS("./Results/annual.imputed.traits.no_woody.rds")

sum = summary(annual.traits.NW.model)
coef = sum$fixed
row.names(coef) = c("Intercept","Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter",
                    "SRL", "RTD","RMF","DSI","Precipitation")
coef.annual = coef %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#979461"))

annual.plot = ggplot(data = coef.annual,
                          aes(y = factor(row.names(coef.annual), levels = rev(row.names(coef.annual))),
                              x = Estimate,
                              xmin = `l-95% CI`,
                              xmax = `u-95% CI`,
                              fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#979461") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Annuals")
annual.plot

perennial.traits.NW.model = readRDS("./Results/perennial.imputed.traits.no_woody.rds")

sum = summary(perennial.traits.NW.model)
coef = sum$fixed
row.names(coef) = c("Intercept","Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter",
                    "SRL", "RTD","RMF","DSI","Precipitation")
coef.perennial = coef %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#F1C646"))

perennial.plot = ggplot(data = coef.perennial,
                     aes(y = factor(row.names(coef.perennial), levels = rev(row.names(coef.perennial))),
                         x = Estimate,
                         xmin = `l-95% CI`,
                         xmax = `u-95% CI`,
                         fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#F1C646") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Perennials")
perennial.plot

forb.traits.NW.model = readRDS("./Results/forb.legume.imputed.traits.no_woody.rds")

sum = summary(forb.traits.NW.model)
coef = sum$fixed
row.names(coef) = c("Intercept","Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter",
                    "SRL", "RTD","RMF","DSI","Precipitation")
coef.forb = coef %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#F17236"))

forb.plot = ggplot(data = coef.forb,
                        aes(y = factor(row.names(coef.forb), levels = rev(row.names(coef.forb))),
                            x = Estimate,
                            xmin = `l-95% CI`,
                            xmax = `u-95% CI`,
                            fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#F17236") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Forbs")
forb.plot

graminoid.traits.NW.model = readRDS("./Results/graminoid.imputed.traits.no_woody.rds")

sum = summary(graminoid.traits.NW.model)
coef = sum$fixed
row.names(coef) = c("Intercept","Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter",
                    "SRL", "RTD","RMF","DSI","Precipitation")
coef.grass = coef %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#6E687E"))

graminoid.plot = ggplot(data = coef.grass,
                   aes(y = factor(row.names(coef.grass), levels = rev(row.names(coef.grass))),
                       x = Estimate,
                       xmin = `l-95% CI`,
                       xmax = `u-95% CI`,
                       fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#6E687E") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Graminoids")
graminoid.plot

annual.forb.traits.NW.model = readRDS("./Results/annual.forb.legume.imputed.traits.no_woody.rds")

sum = summary(annual.forb.traits.NW.model)
coef = sum$fixed
row.names(coef) = c("Intercept","Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter",
                    "SRL", "RTD","RMF","DSI","Precipitation")
coef.annual.forb = coef %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#B50200"))

annual.forb.plot = ggplot(data = coef.annual.forb,
                    aes(y = factor(row.names(coef.annual.forb), levels = rev(row.names(coef.annual.forb))),
                        x = Estimate,
                        xmin = `l-95% CI`,
                        xmax = `u-95% CI`,
                        fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#B50200") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Annual Forbs")
annual.forb.plot

perennial.graminoid.traits.NW.model = readRDS("./Results/perennial.graminoid.imputed.traits.no_woody.rds")

sum = summary(perennial.graminoid.traits.NW.model)
coef = sum$fixed
row.names(coef) = c("Intercept","Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter",
                    "SRL", "RTD","RMF","DSI","Precipitation")
coef.perennial.grass = coef %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#6089B5"))

perennial.graminoid.plot = ggplot(data = coef.perennial.grass,
                    aes(y = factor(row.names(coef.perennial.grass), levels = rev(row.names(coef.perennial.grass))),
                        x = Estimate,
                        xmin = `l-95% CI`,
                        xmax = `u-95% CI`,
                        fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#6089B5") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Perennial Graminoids")
perennial.graminoid.plot

perennial.forb.traits.NW.model = readRDS("./Results/perennial.forb.legume.imputed.traits.no_woody.rds")

sum = summary(perennial.forb.traits.NW.model)
coef = sum$fixed
row.names(coef) = c("Intercept","Leaf N","Height","Root N","SLA","Rooting Depth","Root Diameter",
                    "SRL", "RTD","RMF","DSI","Precipitation")
coef.perennial.forb = coef %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "black"))

perennial.forb.plot = ggplot(data = coef.perennial.forb,
                              aes(y = factor(row.names(coef.perennial.forb), levels = rev(row.names(coef.perennial.forb))),
                                  x = Estimate,
                                  xmin = `l-95% CI`,
                                  xmax = `u-95% CI`,
                                  fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "black") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Perennial Forbs")
perennial.forb.plot

png(filename = "./Plots/coef.plot.png",  width = 11, 
    height = 11, units = "in", res = 600)

plot_grid(All.species.plot,annual.plot,perennial.plot,graminoid.plot,forb.plot,
                     annual.forb.plot,perennial.forb.plot,perennial.graminoid.plot, labels = NULL)

dev.off()

#### Trait Interaction Dot Plots ####
imputed.NW.traits.height.leafN = readRDS("./Results/all.imputed.traits.NW.height.leafN.rds")
bayes_R2(imputed.NW.traits.height.leafN)
# R2 0.04620861 0.02366104 0.01308091 0.1012018
imputed.NW.traits.depth.leafN = readRDS("./Results/all.imputed.NW.traits.depth.leafN.rds")
bayes_R2(imputed.NW.traits.depth.leafN)
#R2 0.04902346 0.02481235 0.01394253 0.1094119
imputed.NW.traits.RMF.leafN = readRDS("./Results/all.imputed.NW.traits.leafN.RMF.rds")
bayes_R2(imputed.NW.traits.RMF.leafN)
#R2 0.05056307 0.02389004 0.01587269 0.1095646
imputed.NW.traits.RTD.SRL = readRDS("./Results/all.imputed.NW.traits.RTD.SRL.rds")
bayes_R2(imputed.NW.traits.RTD.SRL)
# R2 0.04620861 0.02366104 0.01308091 0.1012018

sum = summary(imputed.NW.traits.height.leafN)$fixed
sum.2 = summary(imputed.NW.traits.depth.leafN)$fixed
sum.3 = summary(imputed.NW.traits.RMF.leafN)$fixed
sum.4 = summary(imputed.NW.traits.RTD.SRL)$fixed

coefs = as.data.frame(sum[6,c(1,3,4)])
coefs[2,1:3] = sum.2[6,c(1,3,4)]
coefs[3,1:3] = sum.3[6,c(1,3,4)]
coefs[4,1:3] = sum.4[6,c(1,3,4)]

row.names(coefs) = c("Leaf N*Height", "Leaf N*Rooting Depth", "Leaf N*RMF", "RTD*SRL")

coefs.2 = coefs %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#769370"))

All.species.plot = ggplot(data = coefs.2,
                          aes(y = factor(row.names(coefs.2), levels = rev(row.names(coefs.2))),
                              x = Estimate,
                              xmin = `l-95% CI`,
                              xmax = `u-95% CI`,
                              fill = resp.fill)) +
  geom_pointrange(size = 1, shape = 21, color = "#769370") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "All-Species")
All.species.plot

## Annual

imputed.NW.traits.height.leafN = readRDS("./Results/annual.imputed.traits.NW.height.leafN.rds")
bayes_R2(imputed.NW.traits.height.leafN)
# R2 0.1082466 0.04328998 0.03755913 0.2048849
imputed.NW.traits.depth.leafN = readRDS("./Results/annual.imputed.traits.NW.depth.leafN.rds")
bayes_R2(imputed.NW.traits.depth.leafN)
#R2 0.1193793 0.04317675 0.04622809 0.2148366
imputed.NW.traits.RMF.leafN = readRDS("./Results/annual.imputed.NW.traits.leafN.RMF.rds")
bayes_R2(imputed.NW.traits.RMF.leafN)
#R2 0.1310284 0.04710103 0.0522819 0.2401864
imputed.NW.traits.RTD.SRL = readRDS("./Results/annual.imputed.traits.NW.RTD.SRL.rds")
bayes_R2(imputed.NW.traits.RTD.SRL)
# R2 0.1037567 0.0421373 0.03488088 0.1980414

sum = summary(imputed.NW.traits.height.leafN)$fixed
sum.2 = summary(imputed.NW.traits.depth.leafN)$fixed
sum.3 = summary(imputed.NW.traits.RMF.leafN)$fixed
sum.4 = summary(imputed.NW.traits.RTD.SRL)$fixed

coefs = as.data.frame(sum[6,c(1,3,4)])
coefs[2,1:3] = sum.2[6,c(1,3,4)]
coefs[3,1:3] = sum.3[6,c(1,3,4)]
coefs[4,1:3] = sum.4[6,c(1,3,4)]

row.names(coefs) = c("Leaf N*Height", "Leaf N*Rooting Depth", "Leaf N*RMF", "RTD*SRL")

coef.annual = coefs %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#979461"))

annual.plot = ggplot(data = coef.annual,
                     aes(y = factor(row.names(coef.annual), levels = rev(row.names(coef.annual))),
                         x = Estimate,
                         xmin = `l-95% CI`,
                         xmax = `u-95% CI`,
                         fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#979461") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Annuals")
annual.plot

## perennial

imputed.NW.traits.height.leafN = readRDS("./Results/perennial.imputed.traits.NW.height.leafN.rds")
bayes_R2(imputed.NW.traits.height.leafN)
# R2 0.06295824 0.03318674 0.01677983 0.1430209
imputed.NW.traits.depth.leafN = readRDS("./Results/perennial.imputed.traits.NW.depth.leafN.rds")
bayes_R2(imputed.NW.traits.depth.leafN)
#R2 0.05999619 0.03515212 0.01279108 0.1446982
imputed.NW.traits.RMF.leafN = readRDS("./Results/perennial.imputed.traits.NW.leafN.RMF.rds")
bayes_R2(imputed.NW.traits.RMF.leafN)
#R2 0.05835542 0.03305554 0.01291057 0.1367004
imputed.NW.traits.RTD.SRL = readRDS("./Results/perennial.imputed.traits.NW.RTD.SRL.rds")
bayes_R2(imputed.NW.traits.RTD.SRL)
# R2 0.06112001 0.03264814 0.01455853 0.1403494

sum = summary(imputed.NW.traits.height.leafN)$fixed
sum.2 = summary(imputed.NW.traits.depth.leafN)$fixed
sum.3 = summary(imputed.NW.traits.RMF.leafN)$fixed
sum.4 = summary(imputed.NW.traits.RTD.SRL)$fixed

coefs = as.data.frame(sum[6,c(1,3,4)])
coefs[2,1:3] = sum.2[6,c(1,3,4)]
coefs[3,1:3] = sum.3[6,c(1,3,4)]
coefs[4,1:3] = sum.4[6,c(1,3,4)]

row.names(coefs) = c("Leaf N*Height", "Leaf N*Rooting Depth", "Leaf N*RMF", "RTD*SRL")

coef.perennial = coefs %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#F1C646"))

perennial.plot = ggplot(data = coef.perennial,
                        aes(y = factor(row.names(coef.perennial), levels = rev(row.names(coef.perennial))),
                            x = Estimate,
                            xmin = `l-95% CI`,
                            xmax = `u-95% CI`,
                            fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#F1C646") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Perennials")
perennial.plot

## forb

imputed.NW.traits.height.leafN = readRDS("./Results/forb.legume.imputed.traits.NW.height.leafN.rds")
bayes_R2(imputed.NW.traits.height.leafN)
# R2 0.09349774 0.04404204 0.02967335 0.1955086
imputed.NW.traits.depth.leafN = readRDS("./Results/forb.legume.imputed.traits.NW.depth.leafN.rds")
bayes_R2(imputed.NW.traits.depth.leafN)
#R2 0.09146164 0.04923645 0.02604119 0.2098301
imputed.NW.traits.RMF.leafN = readRDS("./Results/forb.legume.imputed.traits.NW.leafN.RMF.rds")
bayes_R2(imputed.NW.traits.RMF.leafN)
#R2 0.08990255 0.04675457 0.02423598 0.2022656
imputed.NW.traits.RTD.SRL = readRDS("./Results/forb.legume.imputed.traits.NW.RTD.SRL.rds")
bayes_R2(imputed.NW.traits.RTD.SRL)
# R2 0.08990255 0.04675457 0.02423598 0.2022656

sum = summary(imputed.NW.traits.height.leafN)$fixed
sum.2 = summary(imputed.NW.traits.depth.leafN)$fixed
sum.3 = summary(imputed.NW.traits.RMF.leafN)$fixed
sum.4 = summary(imputed.NW.traits.RTD.SRL)$fixed

coefs = as.data.frame(sum[6,c(1,3,4)])
coefs[2,1:3] = sum.2[6,c(1,3,4)]
coefs[3,1:3] = sum.3[6,c(1,3,4)]
coefs[4,1:3] = sum.4[6,c(1,3,4)]

row.names(coefs) = c("Leaf N*Height", "Leaf N*Rooting Depth", "Leaf N*RMF", "RTD*SRL")

coef.forb = coefs %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#F17236"))

forb.plot = ggplot(data = coef.forb,
                   aes(y = factor(row.names(coef.forb), levels = rev(row.names(coef.forb))),
                       x = Estimate,
                       xmin = `l-95% CI`,
                       xmax = `u-95% CI`,
                       fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#F17236") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Forbs")
forb.plot

## graminoid

imputed.NW.traits.height.leafN = readRDS("./Results/graminoid.imputed.traits.NW.height.leafN.rds")
bayes_R2(imputed.NW.traits.height.leafN)
# R2 0.07519414 0.03773959 0.02012363 0.1623326
imputed.NW.traits.depth.leafN = readRDS("./Results/graminoid.imputed.traits.NW.depth.leafN.rds")
bayes_R2(imputed.NW.traits.depth.leafN)
#R2 0.07524507 0.03920688 0.01940286 0.1673162
imputed.NW.traits.RMF.leafN = readRDS("./Results/graminoid.imputed.traits.NW.leafN.RMF.rds")
bayes_R2(imputed.NW.traits.RMF.leafN)
#R2 0.08081603 0.03934892 0.02171268 0.1741477
imputed.NW.traits.RTD.SRL = readRDS("./Results/graminoid.imputed.traits.NW.RTD.SRL.rds")
bayes_R2(imputed.NW.traits.RTD.SRL)
# R2 0.07539285 0.03883973 0.01926746 0.1642017

sum = summary(imputed.NW.traits.height.leafN)$fixed
sum.2 = summary(imputed.NW.traits.depth.leafN)$fixed
sum.3 = summary(imputed.NW.traits.RMF.leafN)$fixed
sum.4 = summary(imputed.NW.traits.RTD.SRL)$fixed

coefs = as.data.frame(sum[6,c(1,3,4)])
coefs[2,1:3] = sum.2[6,c(1,3,4)]
coefs[3,1:3] = sum.3[6,c(1,3,4)]
coefs[4,1:3] = sum.4[6,c(1,3,4)]

row.names(coefs) = c("Leaf N*Height", "Leaf N*Rooting Depth", "Leaf N*RMF", "RTD*SRL")

coef.graminoid = coefs %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#6E687E"))

graminoid.plot = ggplot(data = coef.graminoid,
                    aes(y = factor(row.names(coef.graminoid), levels = rev(row.names(coef.graminoid))),
                        x = Estimate,
                        xmin = `l-95% CI`,
                        xmax = `u-95% CI`,
                        fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#6E687E") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Graminoids")
graminoid.plot

## annual.forb.legume

imputed.NW.traits.height.leafN = readRDS("./Results/annual.forb.legume.imputed.traits.NW.height.leafN.rds")
bayes_R2(imputed.NW.traits.height.leafN)
# R2 0.1879269 0.06705478 0.07328748 0.3352544
imputed.NW.traits.depth.leafN = readRDS("./Results/annual.forb.legume.imputed.traits.NW.depth.leafN.rds")
bayes_R2(imputed.NW.traits.depth.leafN)
#R2 0.19025 0.07039582 0.07139912 0.3456307
imputed.NW.traits.RMF.leafN = readRDS("./Results/annual.forb.legume.imputed.traits.NW.leafN.RMF.rds")
bayes_R2(imputed.NW.traits.RMF.leafN)
#R2 0.2055234 0.07505162 0.07897277 0.367604
imputed.NW.traits.RTD.SRL = readRDS("./Results/annual.forb.legume.imputed.traits.NW.RTD.SRL.rds")
bayes_R2(imputed.NW.traits.RTD.SRL)
# R2 0.1818185 0.06325504 0.0745352 0.3182902

sum = summary(imputed.NW.traits.height.leafN)$fixed
sum.2 = summary(imputed.NW.traits.depth.leafN)$fixed
sum.3 = summary(imputed.NW.traits.RMF.leafN)$fixed
sum.4 = summary(imputed.NW.traits.RTD.SRL)$fixed

coefs = as.data.frame(sum[6,c(1,3,4)])
coefs[2,1:3] = sum.2[6,c(1,3,4)]
coefs[3,1:3] = sum.3[6,c(1,3,4)]
coefs[4,1:3] = sum.4[6,c(1,3,4)]

row.names(coefs) = c("Leaf N*Height", "Leaf N*Rooting Depth", "Leaf N*RMF", "RTD*SRL")

coef.annual.forb = coefs %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#B50200"))

annual.forb.plot = ggplot(data = coef.annual.forb,
                          aes(y = factor(row.names(coef.annual.forb), levels = rev(row.names(coef.annual.forb))),
                              x = Estimate,
                              xmin = `l-95% CI`,
                              xmax = `u-95% CI`,
                              fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#B50200") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Annual Forbs")
annual.forb.plot

## perennial.graminoid

imputed.NW.traits.height.leafN = readRDS("./Results/perennial.graminoid.imputed.traits.NW.height.leafN.rds")
bayes_R2(imputed.NW.traits.height.leafN)
# R2 0.09954704 0.05277377 0.02366747 0.2251619
imputed.NW.traits.depth.leafN = readRDS("./Results/perennial.graminoid.imputed.traits.NW.depth.leafN.rds")
bayes_R2(imputed.NW.traits.depth.leafN)
#R2 0.09243592 0.04775324 0.02259472 0.2015672
imputed.NW.traits.RMF.leafN = readRDS("./Results/perennial.graminoid.imputed.traits.NW.leafN.RMF.rds")
bayes_R2(imputed.NW.traits.RMF.leafN)
#R2 0.1030003 0.05091094 0.0268059 0.2202317
imputed.NW.traits.RTD.SRL = readRDS("./Results/perennial.graminoid.imputed.traits.NW.RTD.SRL.rds")
bayes_R2(imputed.NW.traits.RTD.SRL)
# R2 0.09920834 0.04852001 0.02694776 0.2127798

sum = summary(imputed.NW.traits.height.leafN)$fixed
sum.2 = summary(imputed.NW.traits.depth.leafN)$fixed
sum.3 = summary(imputed.NW.traits.RMF.leafN)$fixed
sum.4 = summary(imputed.NW.traits.RTD.SRL)$fixed

coefs = as.data.frame(sum[6,c(1,3,4)])
coefs[2,1:3] = sum.2[6,c(1,3,4)]
coefs[3,1:3] = sum.3[6,c(1,3,4)]
coefs[4,1:3] = sum.4[6,c(1,3,4)]

row.names(coefs) = c("Leaf N*Height", "Leaf N*Rooting Depth", "Leaf N*RMF", "RTD*SRL")

coef.perennial.graminoid = coefs %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#6089B5"))

perennial.graminoid.plot = ggplot(data = coef.perennial.graminoid,
                              aes(y = factor(row.names(coef.perennial.graminoid), levels = rev(row.names(coef.perennial.graminoid))),
                                  x = Estimate,
                                  xmin = `l-95% CI`,
                                  xmax = `u-95% CI`,
                                  fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#6089B5") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Perennial Graminoids")
perennial.graminoid.plot

## perennial.forb.legume

imputed.NW.traits.height.leafN = readRDS("./Results/perennial.forb.legume.imputed.traits.NW.height.leafN.rds")
bayes_R2(imputed.NW.traits.height.leafN)
# R2 0.08953688 0.04958045 0.02233369 0.2128223
imputed.NW.traits.depth.leafN = readRDS("./Results/perennial.forb.legume.imputed.traits.NW.depth.leafN.rds")
bayes_R2(imputed.NW.traits.depth.leafN)
#R2 0.08383131 0.04823494 0.01869921 0.198834
imputed.NW.traits.RMF.leafN = readRDS("./Results/perennial.forb.legume.imputed.traits.NW.leafN.RMF.rds")
bayes_R2(imputed.NW.traits.RMF.leafN)
#R2 0.08396547 0.05044191 0.01737544 0.2058153
imputed.NW.traits.RTD.SRL = readRDS("./Results/perennial.forb.legume.imputed.traits.NW.RTD.SRL.rds")
bayes_R2(imputed.NW.traits.RTD.SRL)
# R2 0.09020453 0.04950654 0.0222804 0.211576

sum = summary(imputed.NW.traits.height.leafN)$fixed
sum.2 = summary(imputed.NW.traits.depth.leafN)$fixed
sum.3 = summary(imputed.NW.traits.RMF.leafN)$fixed
sum.4 = summary(imputed.NW.traits.RTD.SRL)$fixed

coefs = as.data.frame(sum[6,c(1,3,4)])
coefs[2,1:3] = sum.2[6,c(1,3,4)]
coefs[3,1:3] = sum.3[6,c(1,3,4)]
coefs[4,1:3] = sum.4[6,c(1,3,4)]

row.names(coefs) = c("Leaf N*Height", "Leaf N*Rooting Depth", "Leaf N*RMF", "RTD*SRL")

coef.perennial.forb = coefs %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "black"))

perennial.forb.plot = ggplot(data = coef.perennial.forb,
                             aes(y = factor(row.names(coef.perennial.forb), levels = rev(row.names(coef.perennial.forb))),
                                 x = Estimate,
                                 xmin = `l-95% CI`,
                                 xmax = `u-95% CI`,
                                 fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "black") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Perennial Forbs")
perennial.forb.plot

png(filename = "./Plots/trait.interaction.coef.plot.png",  width = 13, 
    height = 11, units = "in", res = 600)

plot_grid(All.species.plot,annual.plot,perennial.plot,graminoid.plot,forb.plot,
          annual.forb.plot,perennial.forb.plot,perennial.graminoid.plot)

dev.off()

#### Environment Interaction Dot Plots ####
imputed.NW.traits.height.MAP.DSI= readRDS("./Results/all.imputed.NW.traits.height.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.height.MAP.DSI)
# R2 0.03778474 0.02205149 0.008083502 0.09329935
imputed.NW.traits.leafN.MAP.DSI = readRDS("./Results/all.imputed.NW.traits.leafN.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.leafN.MAP.DSI)
#R2 0.06092939 0.0258172 0.02256498 0.1241577
imputed.NW.traits.depth.MAP.DSI = readRDS("./Results/all.imputed.NW.traits.depth.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.depth.MAP.DSI)
#R2 0.04165983 0.02175759 0.01154696 0.09635656
imputed.NW.traits.RTD.MAP.DSI = readRDS("./Results/all.imputed.NW.traits.RTD.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.RTD.MAP.DSI)
# R2 0.04314836 0.02442081 0.0106199 0.103889
imputed.NW.traits.SRL.MAP.DSI = readRDS("./Results/all.imputed.NW.traits.SRL.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.SRL.MAP.DSI)
# R2 0.04075575 0.02239753 0.009678036 0.09585312
imputed.NW.traits.RMF.MAP.DSI = readRDS("./Results/all.imputed.NW.traits.RMF.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.RMF.MAP.DSI)
# R2 0.03828717 0.02266265 0.008483114 0.09699009

sum = summary(imputed.NW.traits.height.MAP.DSI)$fixed
sum.2 = summary(imputed.NW.traits.leafN.MAP.DSI)$fixed
sum.3 = summary(imputed.NW.traits.depth.MAP.DSI)$fixed
sum.4 = summary(imputed.NW.traits.RTD.MAP.DSI)$fixed
sum.5 = summary(imputed.NW.traits.SRL.MAP.DSI)$fixed
sum.6 = summary(imputed.NW.traits.RMF.MAP.DSI)$fixed

coefs = as.data.frame(sum[c(5,6),c(1,3,4)])
coefs[c(3,4),1:3] = sum.2[c(5,6),c(1,3,4)]
coefs[c(5,6),1:3] = sum.3[c(5,6),c(1,3,4)]
coefs[c(7,8),1:3] = sum.4[c(5,6),c(1,3,4)]
coefs[c(9,10),1:3] = sum.5[c(5,6),c(1,3,4)]
coefs[c(11,12),1:3] = sum.6[c(5,6),c(1,3,4)]

row.names(coefs) = c("Height*MAP","Height*DSI","Leaf N*MAP","Leaf N*DSI","Rooting Depth*MAP","Rooting Depth*DSI",
                     "RTD*MAP","RTD*DSI","SRL*MAP","SRL*DSI","RMF*MAP","RMF*DSI")

coefs.2 = coefs %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#769370"))

All.species.plot = ggplot(data = coefs.2,
                          aes(y = factor(row.names(coefs.2), levels = rev(row.names(coefs.2))),
                              x = Estimate,
                              xmin = `l-95% CI`,
                              xmax = `u-95% CI`,
                              fill = resp.fill)) +
  geom_pointrange(size = 1, shape = 21, color = "#769370") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "All-Species")
All.species.plot

## Annuals

imputed.NW.traits.height.MAP.DSI= readRDS("./Results/imputed.NW.annual.traits.height.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.height.MAP.DSI)
# R2 0.1208656 0.04510458 0.04694877 0.2220583
imputed.NW.traits.leafN.MAP.DSI = readRDS("./Results/imputed.NW.annual.traits.leafN.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.leafN.MAP.DSI)
#R2 0.1705448 0.05178855 0.07863608 0.2805161
imputed.NW.traits.depth.MAP.DSI = readRDS("./Results/imputed.NW.annual.traits.depth.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.depth.MAP.DSI)
#R2 0.1090259 0.04274418 0.04011061 0.2032013
imputed.NW.traits.RTD.MAP.DSI = readRDS("./Results/imputed.NW.annual.traits.RTD.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.RTD.MAP.DSI)
# R2 0.111892 0.04404864 0.03903641 0.2092014
imputed.NW.traits.SRL.MAP.DSI = readRDS("./Results/imputed.NW.annual.traits.SRL.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.SRL.MAP.DSI)
# R2 0.1068335 0.04325949 0.03579992 0.2052009
imputed.NW.traits.RMF.MAP.DSI = readRDS("./Results/imputed.NW.annual.traits.RMF.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.RMF.MAP.DSI)
# R2 0.1103115 0.04439953 0.03779541 0.2071463

sum = summary(imputed.NW.traits.height.MAP.DSI)$fixed
sum.2 = summary(imputed.NW.traits.leafN.MAP.DSI)$fixed
sum.3 = summary(imputed.NW.traits.depth.MAP.DSI)$fixed
sum.4 = summary(imputed.NW.traits.RTD.MAP.DSI)$fixed
sum.5 = summary(imputed.NW.traits.SRL.MAP.DSI)$fixed
sum.6 = summary(imputed.NW.traits.RMF.MAP.DSI)$fixed

coefs = as.data.frame(sum[c(5,6),c(1,3,4)])
coefs[c(3,4),1:3] = sum.2[c(5,6),c(1,3,4)]
coefs[c(5,6),1:3] = sum.3[c(5,6),c(1,3,4)]
coefs[c(7,8),1:3] = sum.4[c(5,6),c(1,3,4)]
coefs[c(9,10),1:3] = sum.5[c(5,6),c(1,3,4)]
coefs[c(11,12),1:3] = sum.6[c(5,6),c(1,3,4)]

row.names(coefs) = c("Height*MAP","Height*DSI","Leaf N*MAP","Leaf N*DSI","Rooting Depth*MAP","Rooting Depth*DSI",
                     "RTD*MAP","RTD*DSI","SRL*MAP","SRL*DSI","RMF*MAP","RMF*DSI")

coef.annual = coefs %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#979461"))

annual.plot = ggplot(data = coef.annual,
                     aes(y = factor(row.names(coef.annual), levels = rev(row.names(coef.annual))),
                         x = Estimate,
                         xmin = `l-95% CI`,
                         xmax = `u-95% CI`,
                         fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#979461") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Annuals")
annual.plot

## perennials

imputed.NW.traits.height.MAP.DSI= readRDS("./Results/imputed.NW.perennial.traits.height.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.height.MAP.DSI)
# R2 0.05332797 0.03160344 0.01146005 0.13111
imputed.NW.traits.leafN.MAP.DSI = readRDS("./Results/imputed.NW.perennial.traits.leafN.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.leafN.MAP.DSI)
#R2 0.05879655 0.03407371 0.0125315 0.141046
imputed.NW.traits.depth.MAP.DSI = readRDS("./Results/imputed.NW.perennial.traits.depth.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.depth.MAP.DSI)
#R2 0.05988846 0.03295912 0.01457783 0.1382209
imputed.NW.traits.RTD.MAP.DSI = readRDS("./Results/imputed.NW.perennial.traits.RTD.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.RTD.MAP.DSI)
# R2 0.0614589 0.03363447 0.01443582 0.1431149
imputed.NW.traits.SRL.MAP.DSI = readRDS("./Results/imputed.NW.perennial.traits.SRL.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.SRL.MAP.DSI)
# R2 0.05748926 0.03270784 0.0125087 0.1337565
imputed.NW.traits.RMF.MAP.DSI = readRDS("./Results/imputed.NW.perennial.traits.RMF.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.RMF.MAP.DSI)
# R2 0.05922796 0.03369339 0.01331522 0.1395166

sum = summary(imputed.NW.traits.height.MAP.DSI)$fixed
sum.2 = summary(imputed.NW.traits.leafN.MAP.DSI)$fixed
sum.3 = summary(imputed.NW.traits.depth.MAP.DSI)$fixed
sum.4 = summary(imputed.NW.traits.RTD.MAP.DSI)$fixed
sum.5 = summary(imputed.NW.traits.SRL.MAP.DSI)$fixed
sum.6 = summary(imputed.NW.traits.RMF.MAP.DSI)$fixed

coefs = as.data.frame(sum[c(5,6),c(1,3,4)])
coefs[c(3,4),1:3] = sum.2[c(5,6),c(1,3,4)]
coefs[c(5,6),1:3] = sum.3[c(5,6),c(1,3,4)]
coefs[c(7,8),1:3] = sum.4[c(5,6),c(1,3,4)]
coefs[c(9,10),1:3] = sum.5[c(5,6),c(1,3,4)]
coefs[c(11,12),1:3] = sum.6[c(5,6),c(1,3,4)]

row.names(coefs) = c("Height*MAP","Height*DSI","Leaf N*MAP","Leaf N*DSI","Rooting Depth*MAP","Rooting Depth*DSI",
                     "RTD*MAP","RTD*DSI","SRL*MAP","SRL*DSI","RMF*MAP","RMF*DSI")

coef.perennial = coefs %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#F1C646"))

perennial.plot = ggplot(data = coef.perennial,
                        aes(y = factor(row.names(coef.perennial), levels = rev(row.names(coef.perennial))),
                            x = Estimate,
                            xmin = `l-95% CI`,
                            xmax = `u-95% CI`,
                            fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#F1C646") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Perennials")
perennial.plot

## forb.legumes

imputed.NW.traits.height.MAP.DSI= readRDS("./Results/imputed.NW.forb.legume.traits.height.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.height.MAP.DSI)
# R2 0.07520765 0.0402359 0.01865661 0.172068
imputed.NW.traits.leafN.MAP.DSI = readRDS("./Results/imputed.NW.forb.legume.traits.leafN.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.leafN.MAP.DSI)
#R2 0.1104155 0.04644593 0.04187595 0.2197352
imputed.NW.traits.depth.MAP.DSI = readRDS("./Results/imputed.NW.forb.legume.traits.depth.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.depth.MAP.DSI)
#R2  0.07414734 0.04162549 0.01818971 0.1759422
imputed.NW.traits.RTD.MAP.DSI = readRDS("./Results/imputed.NW.forb.legume.traits.RTD.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.RTD.MAP.DSI)
# R2 0.07357648 0.04117196 0.01686191 0.1750251
imputed.NW.traits.SRL.MAP.DSI = readRDS("./Results/imputed.NW.forb.legume.traits.SRL.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.SRL.MAP.DSI)
# R2 0.06303201 0.03884063 0.01470552 0.1649587
imputed.NW.traits.RMF.MAP.DSI = readRDS("./Results/imputed.NW.forb.legume.traits.RMF.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.RMF.MAP.DSI)
# R2 0.06886695 0.04259294 0.0134793 0.1755777

sum = summary(imputed.NW.traits.height.MAP.DSI)$fixed
sum.2 = summary(imputed.NW.traits.leafN.MAP.DSI)$fixed
sum.3 = summary(imputed.NW.traits.depth.MAP.DSI)$fixed
sum.4 = summary(imputed.NW.traits.RTD.MAP.DSI)$fixed
sum.5 = summary(imputed.NW.traits.SRL.MAP.DSI)$fixed
sum.6 = summary(imputed.NW.traits.RMF.MAP.DSI)$fixed

coefs = as.data.frame(sum[c(5,6),c(1,3,4)])
coefs[c(3,4),1:3] = sum.2[c(5,6),c(1,3,4)]
coefs[c(5,6),1:3] = sum.3[c(5,6),c(1,3,4)]
coefs[c(7,8),1:3] = sum.4[c(5,6),c(1,3,4)]
coefs[c(9,10),1:3] = sum.5[c(5,6),c(1,3,4)]
coefs[c(11,12),1:3] = sum.6[c(5,6),c(1,3,4)]

row.names(coefs) = c("Height*MAP","Height*DSI","Leaf N*MAP","Leaf N*DSI","Rooting Depth*MAP","Rooting Depth*DSI",
                     "RTD*MAP","RTD*DSI","SRL*MAP","SRL*DSI","RMF*MAP","RMF*DSI")

coef.forb = coefs %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#F17236"))

forb.plot = ggplot(data = coef.forb,
                   aes(y = factor(row.names(coef.forb), levels = rev(row.names(coef.forb))),
                       x = Estimate,
                       xmin = `l-95% CI`,
                       xmax = `u-95% CI`,
                       fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#F17236") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Forbs")
forb.plot

## graminoids

imputed.NW.traits.height.MAP.DSI= readRDS("./Results/imputed.NW.graminoid.traits.height.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.height.MAP.DSI)
# R2 0.07257722 0.03807303 0.01785619 0.1639368
imputed.NW.traits.leafN.MAP.DSI = readRDS("./Results/imputed.NW.graminoid.traits.leafN.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.leafN.MAP.DSI)
#R2 0.08168753 0.04139105 0.02202018 0.1794329
imputed.NW.traits.depth.MAP.DSI = readRDS("./Results/imputed.NW.graminoid.traits.depth.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.depth.MAP.DSI)
#R2  0.06667337 0.03725247 0.01497841 0.153952
imputed.NW.traits.RTD.MAP.DSI = readRDS("./Results/imputed.NW.graminoid.traits.RTD.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.RTD.MAP.DSI)
# R2 0.07071695 0.0389249 0.01555127 0.163423
imputed.NW.traits.SRL.MAP.DSI = readRDS("./Results/imputed.NW.graminoid.traits.SRL.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.SRL.MAP.DSI)
# R2 0.08167805 0.03947151 0.02165621 0.171412
imputed.NW.traits.RMF.MAP.DSI = readRDS("./Results/imputed.NW.graminoid.traits.RMF.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.RMF.MAP.DSI)
# R2 0.07350117 0.03865609 0.01877043 0.1640502

sum = summary(imputed.NW.traits.height.MAP.DSI)$fixed
sum.2 = summary(imputed.NW.traits.leafN.MAP.DSI)$fixed
sum.3 = summary(imputed.NW.traits.depth.MAP.DSI)$fixed
sum.4 = summary(imputed.NW.traits.RTD.MAP.DSI)$fixed
sum.5 = summary(imputed.NW.traits.SRL.MAP.DSI)$fixed
sum.6 = summary(imputed.NW.traits.RMF.MAP.DSI)$fixed

coefs = as.data.frame(sum[c(5,6),c(1,3,4)])
coefs[c(3,4),1:3] = sum.2[c(5,6),c(1,3,4)]
coefs[c(5,6),1:3] = sum.3[c(5,6),c(1,3,4)]
coefs[c(7,8),1:3] = sum.4[c(5,6),c(1,3,4)]
coefs[c(9,10),1:3] = sum.5[c(5,6),c(1,3,4)]
coefs[c(11,12),1:3] = sum.6[c(5,6),c(1,3,4)]

row.names(coefs) = c("Height*MAP","Height*DSI","Leaf N*MAP","Leaf N*DSI","Rooting Depth*MAP","Rooting Depth*DSI",
                     "RTD*MAP","RTD*DSI","SRL*MAP","SRL*DSI","RMF*MAP","RMF*DSI")

coef.graminoid = coefs %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#6E687E"))

graminoid.plot = ggplot(data = coef.graminoid,
                        aes(y = factor(row.names(coef.graminoid), levels = rev(row.names(coef.graminoid))),
                            x = Estimate,
                            xmin = `l-95% CI`,
                            xmax = `u-95% CI`,
                            fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#6E687E") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Graminoids")
graminoid.plot

## annual.forb.legumes

imputed.NW.traits.height.MAP.DSI= readRDS("./Results/imputed.NW.annual.forb.legume.traits.height.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.height.MAP.DSI)
# R2  0.193955 0.06927718 0.07696936 0.3484241
imputed.NW.traits.leafN.MAP.DSI = readRDS("./Results/imputed.NW.annual.forb.legume.traits.leafN.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.leafN.MAP.DSI)
#R2 0.2487581 0.07176966 0.1201865 0.400603
imputed.NW.traits.depth.MAP.DSI = readRDS("./Results/imputed.NW.annual.forb.legume.traits.depth.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.depth.MAP.DSI)
#R2  0.1750077 0.06859499 0.06367586 0.3316238
imputed.NW.traits.RTD.MAP.DSI = readRDS("./Results/imputed.NW.annual.forb.legume.traits.RTD.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.RTD.MAP.DSI)
# R2 0.189431 0.0672549 0.07478275 0.3386802
imputed.NW.traits.SRL.MAP.DSI = readRDS("./Results/imputed.NW.annual.forb.legume.traits.SRL.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.SRL.MAP.DSI)
# R2 0.1601988 0.06443742 0.05627058 0.3042982
imputed.NW.traits.RMF.MAP.DSI = readRDS("./Results/imputed.NW.annual.forb.legume.traits.RMF.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.RMF.MAP.DSI)
# R2  0.1818622 0.07089922 0.06765313 0.3469543

sum = summary(imputed.NW.traits.height.MAP.DSI)$fixed
sum.2 = summary(imputed.NW.traits.leafN.MAP.DSI)$fixed
sum.3 = summary(imputed.NW.traits.depth.MAP.DSI)$fixed
sum.4 = summary(imputed.NW.traits.RTD.MAP.DSI)$fixed
sum.5 = summary(imputed.NW.traits.SRL.MAP.DSI)$fixed
sum.6 = summary(imputed.NW.traits.RMF.MAP.DSI)$fixed

coefs = as.data.frame(sum[c(5,6),c(1,3,4)])
coefs[c(3,4),1:3] = sum.2[c(5,6),c(1,3,4)]
coefs[c(5,6),1:3] = sum.3[c(5,6),c(1,3,4)]
coefs[c(7,8),1:3] = sum.4[c(5,6),c(1,3,4)]
coefs[c(9,10),1:3] = sum.5[c(5,6),c(1,3,4)]
coefs[c(11,12),1:3] = sum.6[c(5,6),c(1,3,4)]

row.names(coefs) = c("Height*MAP","Height*DSI","Leaf N*MAP","Leaf N*DSI","Rooting Depth*MAP","Rooting Depth*DSI",
                     "RTD*MAP","RTD*DSI","SRL*MAP","SRL*DSI","RMF*MAP","RMF*DSI")

coef.annual.forb = coefs %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#B50200"))

annual.forb.plot = ggplot(data = coef.annual.forb,
                          aes(y = factor(row.names(coef.annual.forb), levels = rev(row.names(coef.annual.forb))),
                              x = Estimate,
                              xmin = `l-95% CI`,
                              xmax = `u-95% CI`,
                              fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#B50200") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Annual Forbs")
annual.forb.plot

## perennial.graminoids

imputed.NW.traits.height.MAP.DSI= readRDS("./Results/imputed.NW.perennial.graminoid.traits.height.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.height.MAP.DSI)
# R2  0.09514321 0.05076139 0.02167901 0.2094033
imputed.NW.traits.leafN.MAP.DSI = readRDS("./Results/imputed.NW.perennial.graminoid.traits.leafN.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.leafN.MAP.DSI)
#R2 0.09001695 0.05015101 0.02085408 0.2087146
imputed.NW.traits.depth.MAP.DSI = readRDS("./Results/imputed.NW.perennial.graminoid.traits.depth.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.depth.MAP.DSI)
#R2  0.09170852 0.04970511 0.02077364 0.2091623
imputed.NW.traits.RTD.MAP.DSI = readRDS("./Results/imputed.NW.perennial.graminoid.traits.RTD.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.RTD.MAP.DSI)
# R2 0.09152688 0.0483271 0.02238516 0.2055615
imputed.NW.traits.SRL.MAP.DSI = readRDS("./Results/imputed.NW.perennial.graminoid.traits.SRL.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.SRL.MAP.DSI)
# R2 0.1011934 0.04889843 0.02797906 0.2149136
imputed.NW.traits.RMF.MAP.DSI = readRDS("./Results/imputed.NW.perennial.graminoid.traits.RMF.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.RMF.MAP.DSI)
# R2  0.1039792 0.05003459 0.02738343 0.221052

sum = summary(imputed.NW.traits.height.MAP.DSI)$fixed
sum.2 = summary(imputed.NW.traits.leafN.MAP.DSI)$fixed
sum.3 = summary(imputed.NW.traits.depth.MAP.DSI)$fixed
sum.4 = summary(imputed.NW.traits.RTD.MAP.DSI)$fixed
sum.5 = summary(imputed.NW.traits.SRL.MAP.DSI)$fixed
sum.6 = summary(imputed.NW.traits.RMF.MAP.DSI)$fixed

coefs = as.data.frame(sum[c(5,6),c(1,3,4)])
coefs[c(3,4),1:3] = sum.2[c(5,6),c(1,3,4)]
coefs[c(5,6),1:3] = sum.3[c(5,6),c(1,3,4)]
coefs[c(7,8),1:3] = sum.4[c(5,6),c(1,3,4)]
coefs[c(9,10),1:3] = sum.5[c(5,6),c(1,3,4)]
coefs[c(11,12),1:3] = sum.6[c(5,6),c(1,3,4)]

row.names(coefs) = c("Height*MAP","Height*DSI","Leaf N*MAP","Leaf N*DSI","Rooting Depth*MAP","Rooting Depth*DSI",
                     "RTD*MAP","RTD*DSI","SRL*MAP","SRL*DSI","RMF*MAP","RMF*DSI")

coef.perennial.graminoid = coefs %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#6089B5"))

perennial.graminoid.plot = ggplot(data = coef.perennial.graminoid,
                                  aes(y = factor(row.names(coef.perennial.graminoid), levels = rev(row.names(coef.perennial.graminoid))),
                                      x = Estimate,
                                      xmin = `l-95% CI`,
                                      xmax = `u-95% CI`,
                                      fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "#6089B5") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Perennial Graminoids")
perennial.graminoid.plot

## perennial.forb.legumes

imputed.NW.traits.height.MAP.DSI= readRDS("./Results/imputed.NW.perennial.forb.legume.traits.height.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.height.MAP.DSI)
# R2  0.09176607 0.05293895 0.02110338 0.223103
imputed.NW.traits.leafN.MAP.DSI = readRDS("./Results/imputed.NW.perennial.forb.legume.traits.leafN.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.leafN.MAP.DSI)
#R2 0.08968503 0.05285282 0.01864453 0.2213102
imputed.NW.traits.depth.MAP.DSI = readRDS("./Results/imputed.NW.perennial.forb.legume.traits.depth.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.depth.MAP.DSI)
#R2  0.09036111 0.04854217 0.02345889 0.2111724
imputed.NW.traits.RTD.MAP.DSI = readRDS("./Results/imputed.NW.perennial.forb.legume.traits.RTD.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.RTD.MAP.DSI)
# R2 0.09019763 0.05142541 0.02031147 0.2162557
imputed.NW.traits.SRL.MAP.DSI = readRDS("./Results/imputed.NW.perennial.forb.legume.traits.SRL.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.SRL.MAP.DSI)
# R2 0.07678711 0.04832567 0.0155648 0.1989975
imputed.NW.traits.RMF.MAP.DSI = readRDS("./Results/imputed.NW.perennial.forb.legume.traits.RMF.MAP.DSI.rds")
bayes_R2(imputed.NW.traits.RMF.MAP.DSI)
# R2  0.08799672 0.05499119 0.01528511 0.2228772

sum = summary(imputed.NW.traits.height.MAP.DSI)$fixed
sum.2 = summary(imputed.NW.traits.leafN.MAP.DSI)$fixed
sum.3 = summary(imputed.NW.traits.depth.MAP.DSI)$fixed
sum.4 = summary(imputed.NW.traits.RTD.MAP.DSI)$fixed
sum.5 = summary(imputed.NW.traits.SRL.MAP.DSI)$fixed
sum.6 = summary(imputed.NW.traits.RMF.MAP.DSI)$fixed

coefs = as.data.frame(sum[c(5,6),c(1,3,4)])
coefs[c(3,4),1:3] = sum.2[c(5,6),c(1,3,4)]
coefs[c(5,6),1:3] = sum.3[c(5,6),c(1,3,4)]
coefs[c(7,8),1:3] = sum.4[c(5,6),c(1,3,4)]
coefs[c(9,10),1:3] = sum.5[c(5,6),c(1,3,4)]
coefs[c(11,12),1:3] = sum.6[c(5,6),c(1,3,4)]

row.names(coefs) = c("Height*MAP","Height*DSI","Leaf N*MAP","Leaf N*DSI","Rooting Depth*MAP","Rooting Depth*DSI",
                     "RTD*MAP","RTD*DSI","SRL*MAP","SRL*DSI","RMF*MAP","RMF*DSI")

coef.perennial.forb = coefs %>%
  mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "black"))

perennial.forb.plot = ggplot(data = coef.perennial.forb,
                             aes(y = factor(row.names(coef.perennial.forb), levels = rev(row.names(coef.perennial.forb))),
                                 x = Estimate,
                                 xmin = `l-95% CI`,
                                 xmax = `u-95% CI`,
                                 fill = resp.fill)) +
  geom_pointrange(size = 1,shape = 21,color = "black") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Perennial Forbs")
perennial.forb.plot

png(filename = "./Plots/enviro.interaction.coef.plot.png",  width = 13, 
    height = 11, units = "in", res = 600)

plot_grid(All.species.plot,annual.plot,perennial.plot,graminoid.plot,forb.plot,
          annual.forb.plot,perennial.forb.plot,perennial.graminoid.plot)

dev.off()

#### Plotting random effects ####

imputed.traits.NW.model = readRDS("./Results/all.imputed.traits.no_woody.rds")

rf = ranef(imputed.traits.NW.model)
species.rf = as.data.frame(rf$Taxon)
row.names(species.rf) = stringr::str_to_sentence(row.names(species.rf))
site.rf = as.data.frame(rf$site_code)


species.rf.2 = species.rf %>%
  mutate(coef.zero = Q2.5.Intercept/Q97.5.Intercept  < 0,
         resp.fill = if_else(coef.zero == TRUE, true = "white", false = "#769370"))

species.rf.plot = ggplot(data = species.rf.2,
                          aes(y = row.names(species.rf.2),
                              x = Estimate.Intercept,
                              xmin = Q2.5.Intercept,
                              xmax = Q97.5.Intercept,
                              fill = resp.fill)) +
  geom_pointrange(size = 1, shape = 21, color = "#769370") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic()+
  labs(y = "",x = "Mean with 95% CI", title = "All-Species")
species.rf.plot

#### Load imputed data for interaction models ####

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
#all.cover = read.csv("./Formatted.Data/Revisions/BACI.data.final.csv") %>%
  #select(site_code,Taxon,relative.cover.change)

# remove whitespace to help merging
#imputed.NW$site_code <- str_trim(imputed.NW$site_code)
#all.cover$site_code <- str_trim(all.cover$site_code)

# join together
#imputed.NW = left_join(imputed.NW,all.cover)

# remove all rows that are not annual or perennial
imputed.NW.2 = imputed.NW %>%
  filter(local_lifespan %in% c("ANNUAL","PERENNIAL"))
# lose 21 data points
# make them factors
imputed.NW.2$local_lifespan = as.factor(imputed.NW.2$local_lifespan)
imputed.NW.2$functional_group = as.factor(imputed.NW.2$functional_group)

# subset for columns we need 
imputed.NW.3 = imputed.NW.2 %>%
  select(cover.change,local_lifespan,functional_group,leafN.final,height.final,rootN.final,SLA.final,root.depth.final,
         rootDiam.final,SRL.final,RTD.final,RMF.final,mean.DSI,mean.MAP,site_code,Taxon) %>%
  drop_na()

###### Model with lifespan groups ###### 

priors <- c(prior(normal(0, 10), class = b))

imputed.traits.NW.lifespan.model = brm(cover.change ~ local_lifespan + mean.MAP + 
                                         mean.MAP*local_lifespan + (1|site_code) + (1|Taxon), 
                                       family = gaussian(),
                                       prior = priors,
                                       data = imputed.NW.3)

summary(imputed.traits.NW.lifespan.model)
bayes_R2(imputed.traits.NW.lifespan.model)
# R2 0.03875425 0.02160649 0.009283615 0.09443166
#saveRDS(imputed.traits.NW.lifespan.model, file = "./Results/lifespan.cats.imputed.traits.no_woody.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.traits.NW.lifespan.model, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.traits.NW.lifespan.model, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals: All Species",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.traits.NW.lifespan.model) # p = 0.749
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.traits.NW.lifespan.model)

imputed.traits.NW.lifespan.model = readRDS("./Results/lifespan.cats.imputed.traits.no_woody.rds")

# test for differences between annuals and perennials for MAP
emtrends(imputed.traits.NW.lifespan.model, pairwise ~ local_lifespan, var = "mean.MAP")
# significant difference b/t annual and perennial. Annual higher.
# annual significant, positive

#### Plotting Lifespan Groups ####

MAP.lifespan.group.effect = conditional_effects(imputed.traits.NW.lifespan.model, 
                                                     effects = "mean.MAP:local_lifespan")$`mean.MAP:local_lifespan`

attr(imputed.NW.3$mean.MAP, "scaled:scale")
attr(imputed.NW.3$mean.MAP, "scaled:center")

MAP.lifespan.group.effect$MAP.bt = (MAP.lifespan.group.effect$mean.MAP*487.1451) + 691.6107

MAP.x.lifespan.plot = ggplot(data = MAP.lifespan.group.effect, aes(x = MAP.bt, y = estimate__, group = as.factor(effect2__))) +
  geom_line(aes(color = as.factor(effect2__)), size = 1.5,show.legend = FALSE) +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = as.factor(effect2__)), alpha = 0.5,show.legend = FALSE) +  
  geom_line(aes(y = upper__, color = factor(effect2__)), size = 1.2) +
  geom_line(aes(y = lower__, color = factor(effect2__)), size = 1.2) + 
  labs(x = "Mean Annual Precipitation (mm)", y = "Cover Change (%)", color = "Lifespan") +
  scale_colour_manual(values = c("#979461", "#F1C646"), labels = c("Annuals","Perennials"))+
  scale_fill_manual(values = c("#979461", "#F1C646"))+
  theme_classic()
MAP.x.lifespan.plot

#ggsave("./Plots/Revisions.3/MAP.lifespan.pdf", height = 3, width = 3)
#ggsave("./Plots/Revisions.3/MAP.lifespan.legend.pdf", height = 3, width = 3)


###### Model with functional groups###### 

imputed.NW.functional.group = imputed.NW %>%
  select(cover.change,functional_group,leafN.final,height.final,rootN.final,SLA.final,root.depth.final,
         rootDiam.final,SRL.final,RTD.final,RMF.final,mean.DSI,mean.MAP,site_code,Taxon) %>%
  drop_na()

priors <- c(prior(normal(0, 10), class = b))

imputed.traits.NW.functional.group.model = brm(cover.change ~ functional_group + leafN.final + 
                                                 height.final + RTD.final + 
                                                 leafN.final*functional_group + height.final*functional_group + 
                                                 RTD.final*functional_group + (1|site_code) + (1|Taxon), 
                                               family = gaussian(),
                                               prior = priors,
                                               data = imputed.NW.functional.group)

summary(imputed.traits.NW.functional.group.model)
bayes_R2(imputed.traits.NW.functional.group.model)
#R2 0.05917052 0.02456162 0.02174739 0.1162923
#saveRDS(imputed.traits.NW.functional.group.model, file = "./Results/functional.group.cats.imputed.traits.no_woody.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.traits.NW.functional.group.model, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.traits.NW.functional.group.model, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.traits.NW.functional.group.model) # p = 0.628
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.traits.NW.functional.group.model)

imputed.traits.NW.functional.group.model = readRDS("./Results/functional.group.cats.imputed.traits.no_woody.rds")

# test for differences between functional groups differ
emtrends(imputed.traits.NW.functional.group.model, pairwise ~ functional_group, var = "leafN.final")
# contrasts not significant
# forb significant, positive

emtrends(imputed.traits.NW.functional.group.model, pairwise ~ functional_group, var = "height.final")
# forb and graminoid significantly different
# forb significant, positive

emtrends(imputed.traits.NW.functional.group.model, pairwise ~ functional_group, var = "RTD.final")
# not significant

#### Plotting Functional Groups ####

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

#ggsave("./Plots/Revisions.3/height.functional.group.pdf", height = 3, width = 3)
#ggsave("./Plots/Revisions.3/height.functional.group.legend.pdf", height = 3, width = 3)

###### Model with lifespan x functional group ###### 

# remove annual graminoids
imputed.NW.4 = imputed.NW.3 %>%
  filter(!(local_lifespan %in% c("ANNUAL") & functional_group %in% c("GRAMINOID")))
# group names together
imputed.NW.4$all.cats = paste(imputed.NW.4$local_lifespan, imputed.NW.4$functional_group, sep = "_")

priors <- c(prior(normal(0, 10), class = b))

imputed.traits.NW.all.cats.model = brm(cover.change ~ all.cats + RTD.final + 
                                         RTD.final*all.cats + (1|site_code) + (1|Taxon), 
                                       family = gaussian(),
                                       prior = priors,
                                       data = imputed.NW.4)

summary(imputed.traits.NW.all.cats.model)
bayes_R2(imputed.traits.NW.all.cats.model)
# 0.05418933 0.03061809 0.01208369 0.1269852
# saveRDS(imputed.traits.NW.all.cats.model, file = "./Results/all.cats.imputed.traits.no_woody.rds")

# checking model assumptions
# plots of residuals, normality of errors
model_residuals = as.data.frame(residuals(imputed.traits.NW.all.cats.model, summary = TRUE))
fitted_vals <- as.data.frame(fitted(imputed.traits.NW.all.cats.model, summary = TRUE))
ggplot(model_residuals, aes(x = Estimate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Density of Model Residuals",
       x = "Residuals",
       y = "Density") +
  theme_minimal()
# variance homogeneity
check_heteroscedasticity(imputed.traits.NW.all.cats.model) # p = 0.591
ggplot(data.frame(fitted = fitted_vals$Estimate, resid = model_residuals$Estimate),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted values") +
  theme_minimal()

# variance inflation
check_collinearity(imputed.traits.NW.all.cats.model)

imputed.traits.NW.all.cats.model = readRDS("./Results/all.cats.imputed.traits.no_woody.rds")

# test for differences between combinations of lifespan and functional group for each trait
# also see if any particular group is significant for each trait

emtrends(imputed.traits.NW.all.cats.model, pairwise ~ all.cats, var = "RTD.final")
# perennial forb alone significant, positive
# contrasts not significant

#### Group Interaction Dot Plots ####
imputed.traits.NW.lifespan.model = readRDS("./Results/lifespan.cats.imputed.traits.no_woody.rds")
imputed.traits.NW.functional.group.model = readRDS("./Results/functional.group.cats.imputed.traits.no_woody.rds")
imputed.traits.NW.all.cats.model = readRDS("./Results/all.cats.imputed.traits.no_woody.rds")

sum = summary(imputed.traits.NW.lifespan.model)$fixed
row.names(sum) = c("Intercept","Lifespan\n(Perennial)", "MAP","Lifespan*MAP\n(Perennial)")
sum.2 = sum %>%
  dplyr::mutate(sum.zero = `l-95% CI`/`u-95% CI` < 0,
         resp.fill = if_else(sum.zero == TRUE, true = "white", false = "black"))

lifespan.plot = ggplot(data = sum.2,
                          aes(y = factor(row.names(sum.2), levels = rev(row.names(sum.2))),
                              x = Estimate,
                              xmin = `l-95% CI`,
                              xmax = `u-95% CI`,
                              fill = resp.fill)) +
  geom_pointrange(size = 1, shape = 21, color = "black") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Lifespans")
lifespan.plot

coef = summary(imputed.traits.NW.functional.group.model)$fixed
row.names(coef) = c("Intercept","Functional Group\n(Graminoid)", "Leaf N","Height","RTD",
                    "Functional Group*Leaf N\n(Graminoid)","Functional Group*Height\n(Graminoid)",
                    "Functional Group*RTD\n(Graminoid)")
coef.2 = coef %>%
  dplyr::mutate(coef.zero = `l-95% CI`/`u-95% CI` < 0,
                resp.fill = if_else(coef.zero == TRUE, true = "white", false = "black"))

functional.plot = ggplot(data = coef.2,
                       aes(y = factor(row.names(coef.2), levels = rev(row.names(coef.2))),
                           x = Estimate,
                           xmin = `l-95% CI`,
                           xmax = `u-95% CI`,
                           fill = resp.fill)) +
  geom_pointrange(size = 1, shape = 21, color = "black") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Functional Groups")
functional.plot

cat = summary(imputed.traits.NW.all.cats.model)$fixed
row.names(cat) = c("Intercept","All Groups\n(Perennial Forb)", "All Groups\n(Perennial Graminoid)","RTD",
                    "All Groups*RTD\n(Perennial Forb)","All Groups*RTD\n(Perennial Graminoid)")
cat.2 = cat %>%
  dplyr::mutate(cat.zero = `l-95% CI`/`u-95% CI` < 0,
                resp.fill = if_else(cat.zero == TRUE, true = "white", false = "black"))

cat.plot = ggplot(data = cat.2,
                         aes(y = factor(row.names(cat.2), levels = rev(row.names(cat.2))),
                             x = Estimate,
                             xmin = `l-95% CI`,
                             xmax = `u-95% CI`,
                             fill = resp.fill)) +
  geom_pointrange(size = 1, shape = 21, color = "black") +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 15)+
  labs(y = "",x = "Mean with 95% CI", title = "Lifespans * Functional Groups")
cat.plot

png(filename = "./Plots/group.interaction.coef.plot.png",  width = 13, 
    height = 11, units = "in", res = 600)

plot_grid(lifespan.plot,functional.plot,cat.plot)

dev.off()

#### UNUSED Below Here ####
###### impute FORB traits functional group NW###### 

priors <- c(prior(normal(0, 10), class = b))

forb.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                             SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                           family = gaussian(),
                           prior = priors,
                           data = imputed.NW.forb)

summary(forb.traits.NW.model)
# leafN significant

saveRDS(forb.traits.NW.model, file = "./Results/forb.imputed.traits.no_woody.rds")
bayes_R2(forb.traits.NW.model)
# R2 0.1148149 0.04304049 0.04828436 0.2129581

forb.traits.NW.model = readRDS("./Results/forb.imputed.traits.no_woody.rds")

imputed.NW.forb.traits.height.leafN = brm(cover.change ~ height.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                          family = gaussian(),
                                          prior = priors,
                                          data = imputed.NW.forb)

summary(imputed.NW.forb.traits.height.leafN)
# leafN

imputed.NW.forb.traits.depth.SLA= brm(cover.change ~ root.depth.final*SLA.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                      family = gaussian(),
                                      prior = priors,
                                      data = imputed.NW.forb)

summary(imputed.NW.forb.traits.depth.SLA)
# nothing significant

imputed.NW.forb.traits.depth.leafN = brm(cover.change ~ root.depth.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                         family = gaussian(),
                                         prior = priors,
                                         data = imputed.NW.forb)

summary(imputed.NW.forb.traits.depth.leafN)
# leafN

imputed.NW.forb.traits.RTD.SRL = brm(cover.change ~ RTD.final*SRL.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                     family = gaussian(),
                                     prior = priors,
                                     data = imputed.NW.forb)

summary(imputed.NW.forb.traits.RTD.SRL)
# RTD x SRL significant, negative

saveRDS(imputed.NW.forb.traits.RTD.SRL, file = "./Results/forb.imputed.traits.NW.RTD.SRL.rds")
bayes_R2(imputed.NW.forb.traits.RTD.SRL)
# R2 0.0804646 0.04240653 0.02177166 0.1872374

conditional_effects(imputed.NW.forb.traits.RTD.SRL)
# higher cover with higher RTD and lower SRL or lower RTD with higher SRL

imputed.NW.forb.traits.SLA.RMF = brm(cover.change ~ SLA.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                     family = gaussian(),
                                     prior = priors,
                                     data = imputed.NW.forb)

summary(imputed.NW.forb.traits.SLA.RMF)
# nothing significant

imputed.NW.forb.traits.leafN.RMF = brm(cover.change ~ leafN.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                       family = gaussian(),
                                       prior = priors,
                                       data = imputed.NW.forb)

summary(imputed.NW.forb.traits.leafN.RMF)
# leafN significant positive

###### impute GRASS traits functional group NW ###### 
grass.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                              SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                            family = gaussian(),
                            prior = priors,
                            data = imputed.NW.grass)

summary(grass.traits.NW.model)
# nothing significant

saveRDS(grass.traits.NW.model, file = "./Results/grass.imputed.traits.no_woody.rds")
bayes_R2(grass.traits.NW.model)
# R2 0.1131223 0.04113157 0.04815195 0.2051698

imputed.NW.grass.traits.height.leafN = brm(cover.change ~ height.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                           family = gaussian(),
                                           prior = priors,
                                           data = imputed.NW.grass)

summary(imputed.NW.grass.traits.height.leafN)
# nothing significant

imputed.NW.grass.traits.depth.SLA= brm(cover.change ~ root.depth.final*SLA.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                       family = gaussian(),
                                       prior = priors,
                                       data = imputed.NW.grass)

summary(imputed.NW.grass.traits.depth.SLA)
# nothing significant

imputed.NW.grass.traits.depth.leafN = brm(cover.change ~ root.depth.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                          family = gaussian(),
                                          prior = priors,
                                          data = imputed.NW.grass)

summary(imputed.NW.grass.traits.depth.leafN)
# MAP is positive, significant

imputed.NW.grass.traits.RTD.SRL = brm(cover.change ~ RTD.final*SRL.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                      family = gaussian(),
                                      prior = priors,
                                      data = imputed.NW.grass)

summary(imputed.NW.grass.traits.RTD.SRL)
# nothing significant

imputed.NW.grass.traits.SLA.RMF = brm(cover.change ~ SLA.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                      family = gaussian(),
                                      prior = priors,
                                      data = imputed.NW.grass)

summary(imputed.NW.grass.traits.SLA.RMF)
# nothing significant

imputed.NW.grass.traits.leafN.RMF = brm(cover.change ~ leafN.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                        family = gaussian(),
                                        prior = priors,
                                        data = imputed.NW.grass)

summary(imputed.NW.grass.traits.leafN.RMF)
# nothing significant

###### imputed traits ANNUAL FORB NW###### 

priors <- c(prior(normal(0, 10), class = b))

annual.forb.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                    SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                  family = gaussian(),
                                  prior = priors,
                                  data = imputed.NW.annual.forb)

summary(annual.forb.traits.NW.model)
# MAP significant

saveRDS(annual.forb.traits.NW.model, file = "./Results/annual.forb.imputed.traits.no_woody.rds")
bayes_R2(annual.forb.traits.NW.model)
# R2 0.2252887 0.0659826 0.1116487 0.3650369

imputed.NW.annual.forb.traits.height.leafN = brm(cover.change ~ height.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                 family = gaussian(),
                                                 prior = priors,
                                                 data = imputed.NW.annual.forb)

summary(imputed.NW.annual.forb.traits.height.leafN)
# nothing significant

imputed.NW.annual.forb.traits.depth.SLA= brm(cover.change ~ root.depth.final*SLA.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                             family = gaussian(),
                                             prior = priors,
                                             data = imputed.NW.annual.forb)

summary(imputed.NW.annual.forb.traits.depth.SLA)
# MAP

imputed.NW.annual.forb.traits.depth.leafN = brm(cover.change ~ root.depth.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                family = gaussian(),
                                                prior = priors,
                                                data = imputed.NW.annual.forb)

summary(imputed.NW.annual.forb.traits.depth.leafN)
# MAP is positive, significant

imputed.NW.annual.forb.traits.RTD.SRL = brm(cover.change ~ RTD.final*SRL.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                            family = gaussian(),
                                            prior = priors,
                                            data = imputed.NW.annual.forb)

summary(imputed.NW.annual.forb.traits.RTD.SRL)
# nothing significant

imputed.NW.annual.forb.traits.SLA.RMF = brm(cover.change ~ SLA.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                            family = gaussian(),
                                            prior = priors,
                                            data = imputed.NW.annual.forb)

summary(imputed.NW.annual.forb.traits.SLA.RMF)
# MAP

imputed.NW.annual.forb.traits.leafN.RMF = brm(cover.change ~ leafN.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                              family = gaussian(),
                                              prior = priors,
                                              data = imputed.NW.annual.forb)

summary(imputed.NW.annual.forb.traits.leafN.RMF)
# MAP significant positive

###### imputed traits PERENNIAL FORB NW ###### 

perennial.forb.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                       SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                     family = gaussian(),
                                     prior = priors,
                                     data = imputed.NW.perennial.forb)

summary(perennial.forb.traits.NW.model)
# nothing significant

saveRDS(perennial.forb.traits.NW.model, file = "./Results/perennial.forb.imputed.traits.no_woody.rds")
bayes_R2(perennial.forb.traits.NW.model)
#R2 0.1475479 0.05794628 0.06079698 0.2874462

imputed.NW.perennial.forb.traits.height.leafN = brm(cover.change ~ height.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                    family = gaussian(),
                                                    prior = priors,
                                                    data = imputed.NW.perennial.forb)

summary(imputed.NW.perennial.forb.traits.height.leafN)
# leafN

imputed.NW.perennial.forb.traits.depth.SLA= brm(cover.change ~ root.depth.final*SLA.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                family = gaussian(),
                                                prior = priors,
                                                data = imputed.NW.perennial.forb)

summary(imputed.NW.perennial.forb.traits.depth.SLA)
# nothing significant

imputed.NW.perennial.forb.traits.depth.leafN = brm(cover.change ~ root.depth.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                   family = gaussian(),
                                                   prior = priors,
                                                   data = imputed.NW.perennial.forb)

summary(imputed.NW.perennial.forb.traits.depth.leafN)
# leafN

imputed.NW.perennial.forb.traits.RTD.SRL = brm(cover.change ~ RTD.final*SRL.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                               family = gaussian(),
                                               prior = priors,
                                               data = imputed.NW.perennial.forb)

summary(imputed.NW.perennial.forb.traits.RTD.SRL)
# nothing significant

imputed.NW.perennial.forb.traits.SLA.RMF = brm(cover.change ~ SLA.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                               family = gaussian(),
                                               prior = priors,
                                               data = imputed.NW.perennial.forb)

summary(imputed.NW.perennial.forb.traits.SLA.RMF)
# nothing significant

imputed.NW.perennial.forb.traits.leafN.RMF = brm(cover.change ~ leafN.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                 family = gaussian(),
                                                 prior = priors,
                                                 data = imputed.NW.perennial.forb)

summary(imputed.NW.perennial.forb.traits.leafN.RMF)
# leafN

###### imputed traits PERENNIAL GRASS NW ###### 

perennial.grass.traits.NW.model = brm(cover.change ~ leafN.final + height.final + rootN.final + SLA.final + root.depth.final + rootDiam.final +
                                        SRL.final + RTD.final + RMF.final + mean.DSI + mean.MAP +(1|site_code) + (1|Taxon), 
                                      family = gaussian(),
                                      prior = priors,
                                      data = imputed.NW.perennial.grass)

summary(perennial.grass.traits.NW.model)
# nothing significant

saveRDS(perennial.grass.traits.NW.model, file = "./Results/perennial.grass.imputed.traits.no_woody.rds")
bayes_R2(perennial.grass.traits.NW.model)
#R2 0.1499307 0.06219888 0.05914333 0.2792502

imputed.NW.perennial.grass.traits.height.leafN = brm(cover.change ~ height.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                     family = gaussian(),
                                                     prior = priors,
                                                     data = imputed.NW.perennial.grass)

summary(imputed.NW.perennial.grass.traits.height.leafN)
# nothing significant

imputed.NW.perennial.grass.traits.depth.SLA= brm(cover.change ~ root.depth.final*SLA.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                 family = gaussian(),
                                                 prior = priors,
                                                 data = imputed.NW.perennial.grass)

summary(imputed.NW.perennial.grass.traits.depth.SLA)
# nothing significant

imputed.NW.perennial.grass.traits.depth.leafN = brm(cover.change ~ root.depth.final*leafN.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                    family = gaussian(),
                                                    prior = priors,
                                                    data = imputed.NW.perennial.grass)

summary(imputed.NW.perennial.grass.traits.depth.leafN)
# nothing significant

imputed.NW.perennial.grass.traits.RTD.SRL = brm(cover.change ~ RTD.final*SRL.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                family = gaussian(),
                                                prior = priors,
                                                data = imputed.NW.perennial.grass)

summary(imputed.NW.perennial.grass.traits.RTD.SRL)
# nothing significant

imputed.NW.perennial.grass.traits.SLA.RMF = brm(cover.change ~ SLA.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                family = gaussian(),
                                                prior = priors,
                                                data = imputed.NW.perennial.grass)

summary(imputed.NW.perennial.grass.traits.SLA.RMF)
# nothing significant

imputed.NW.perennial.grass.traits.leafN.RMF = brm(cover.change ~ leafN.final*RMF.final + mean.DSI + mean.MAP + (1|site_code) + (1|Taxon), 
                                                  family = gaussian(),
                                                  prior = priors,
                                                  data = imputed.NW.perennial.grass)

summary(imputed.NW.perennial.grass.traits.leafN.RMF)

