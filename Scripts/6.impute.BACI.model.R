# Script to perform analyses where % cover calculated using BACI method
# drought after-drought before) - (control after - control before)

# Using the imputed trait data

library(tidyverse)
library(dismo)
library(gbm)
library(ggBRT)
library(multcomp)
library(cowplot)
library(corrplot)
#library(treezy)

# edited functions from ggBRT needed for plotting
# ggPD_boot_edit.R
# ggPD_boot_test_2.R, use because issues with xlim from ggPD_boot_edit.R
# ggInfluence_edit.R
# ggInteract_2D_edit.R

#### read in all the data frames needs for the analyses ####

# combine data with enviro and remove cover change outlier values

enviro = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)

# imputed traits with woody species
imputed = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.outliersRM.csv", row.names = 1) %>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.annual = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.annuals.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.perennial = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.perennials.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.forbs.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.grass = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.grass.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.annual.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.annual.forbs.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000)
imputed.perennial.forb = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.perennial.forb.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) 
imputed.perennial.grass = read.csv("./Formatted.Data/Revisions/Final.Data/imputed.traits.perennial.grass.outliersRM.csv", row.names = 1)%>%
  left_join(., enviro, by = "site_code") %>%
  filter(round(cover.change, 5) > -24.50000 & round(cover.change, 5) < 23.52000) 

#### Remove NAs ####

# remove NAs from climate data so dataset is the same as the one used in lmer models

all.data = imputed %>%
  dplyr::select(cover.change,leafN.final:mean.MAP) %>%
  drop_na()
annual.data = imputed.annual %>%
  dplyr::select(cover.change,leafN.final:mean.MAP) %>%
  drop_na()
perennial.data = imputed.perennial %>%
  dplyr::select(cover.change,leafN.final:mean.MAP) %>%
  drop_na()
grass = imputed.grass %>%
  dplyr::select(cover.change,leafN.final:mean.MAP) %>%
  drop_na()
forb = imputed.forb %>%
  dplyr::select(cover.change,leafN.final:mean.MAP) %>%
  drop_na()
grass.perennial = imputed.perennial.grass %>%
  dplyr::select(cover.change,leafN.final:mean.MAP) %>%
  drop_na()
forb.annual = imputed.annual.forb %>%
  dplyr::select(cover.change,leafN.final:mean.MAP) %>%
  drop_na()
forb.perennial = imputed.perennial.forb %>%
  dplyr::select(cover.change,leafN.final:mean.MAP) %>%
  drop_na()

#### determining best parameter combination to generate 1000 trees ####
# bag fraction of 0.50 and 0.75, step.size of 25 and 50, tc = 10, 

all.data.map.dsi=gbm.step(data=all.data, gbm.x = c(2:12), gbm.y=1,
                          family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                          bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(all.data.map.dsi)

annual.data.map.dsi=gbm.step(data=annual.data, gbm.x = c(2:12), gbm.y=1,
                             family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                             bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(annual.data.map.dsi)

perennial.data.map.dsi=gbm.step(data=perennial.data, gbm.x = c(2:12), gbm.y=1,
                                family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                                bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(perennial.data.map.dsi)

grass.data.map.dsi=gbm.step(data=grass, gbm.x = c(2:12), gbm.y=1,
                            family = "gaussian", tree.complexity = 10, learning.rate = 0.0000001,
                            bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(grass.data.map.dsi)
# won't add trees

forb.data.map.dsi=gbm.step(data=forb, gbm.x = c(2:12), gbm.y=1,
                           family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                           bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(forb.data.map.dsi)

grass.perennial.data.map.dsi=gbm.step(data=grass.perennial, gbm.x = c(2:12), gbm.y=1,
                                   family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                                   bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(grass.perennial.data.map.dsi)

forb.annual.data.map.dsi=gbm.step(data=forb.annual, gbm.x = c(2:12), gbm.y=1,
                                  family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                                  bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(forb.annual.data.map.dsi)

forb.perennial.data.map.dsi=gbm.step(data=forb.perennial, gbm.x = c(2:12), gbm.y=1,
                                     family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                                     bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(forb.perennial.data.map.dsi)

#### determining best tree complexity to use ####

# code edited to run for each of the data sets
# when loops stopped prematurely, code was run manually to generate 10 reps

# all variables
R2Obs.all.variables <- list()
importancePred.all.variables <- list()
nreps <- 10 #number of simulations
for (tcomp in 10:10) {
  R2Obs.all.variables[[tcomp]] <- numeric(nreps)
  importancePred.all.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:11),
                                                             ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    
    BRT.all.variables <- gbm.step(data=forb.perennial,
                                  gbm.x = c(2:12),
                                  gbm.y = 1,
                                  family = "gaussian",
                                  tree.complexity = tcomp,
                                  learning.rate = 0.0001,
                                  bag.fraction = 0.75,
                                  n.trees = 50,
                                  step.size = 50,
                                  plot.main=F, plot.folds=F)
    
    
    #R2 adj:
    R2Obs.all.variables[[tcomp]][i] <- 1 - (BRT.all.variables$self.statistics$mean.resid /
                                              BRT.all.variables$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.all.variables[[tcomp]]) <- sort(rownames(summary(BRT.all.variables)))
    }
    importancePred.all.variables[[tcomp]][, i] <-
      summary(BRT.all.variables)[rownames(importancePred.all.variables[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity

means <- sapply(R2Obs.all.variables, mean)
sds <- sapply(R2Obs.all.variables, sd)
plot(1:length(R2Obs.all.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.all.variables)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

tcFactor <- as.factor(rep(1:10, each=nreps))
R2Vector <- unlist(R2Obs.all.variables)
model <- lm(R2Vector~tcFactor)
TukeyModel<-glht(model, linfct = mcp(tcFactor="Tukey"))
TukeyLetters <- cld(TukeyModel)$mcletters$Letters

# selected the lowest tc value that did not show significant differences 
# compared to the largest tc (10) value
# all.data tc = 9
# annual tc = 3
# perennial tc = 9
# grass tc = NA
# forb tc = 8
# grass.annual = NA, sample size too low
# grass.perennial = 2
# forb.annual = 2
# forb.perennial = 6

#### all data ####

all.data.map.dsi = gbm.step(data=all.data, gbm.x = c(2:12), gbm.y=1,
                          family = "gaussian", tree.complexity = 9, learning.rate = 0.0001,
                          bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(all.data.map.dsi)
# 2550 trees Per.Expl = 4.93%

#saveRDS(all.data.map.dsi, file = "./Formatted.Data/Revisions/all.species.model.rds")

ggInfluence(all.data.map.dsi)
# MAP, leafN, rootN, SLA are significant

all.data.map.dsi$contributions$var = c("Precipitation","LeafN","RootN","SLA",
                                       "RTD","Height","DSI","Root Diameter",
                                       "SRL","Rooting Depth","RMF")

all.data.dsi.influence.plot = ggInfluence_test(all.data.map.dsi, main = expression("All Species (n = 750), R"^2*" = 6.62%"), 
            col.bar = c("#769370","#769370","#769370","#769370",
                        "gray70","gray70","gray70","gray70","gray70","gray70"), 
            col.signif = "#B50200")

# saved as pdf from device (4 x 5)

# get data to plot partial dependency plots

all.data.map.dsi$contributions$var = c("mean.MAP","leafN.final","rootN.final","SLA.final",
                                       "RTD.final","height.final","mean.DSI","rootDiam.final",
                                       "SRL.final","root.depth.final","RMF.final")
                                        

all.data.dsi.prerun<- plot.gbm.4list(all.data.map.dsi)

all.data.dsi.boot <- gbm.bootstrap.functions(all.data.map.dsi, list.predictors=all.data.dsi.prerun, n.reps=1000)

all.data.map.dsi$contributions$var = c("Precipitation","LeafN","RootN","SLA",
                                       "RTD","Height","DSI","Root Diameter",
                                       "SRL","Rooting Depth","RMF")

all.data.map.dsi$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","Root Diameter","SRL","RTD","RMF",
                                              "DSI","Precipitation")
colnames(all.data.map.dsi$gbm.call$dataframe)[2:12] = c("LeafN","Height","RootN","SLA","Rooting Depth","Root Diameter","SRL","RTD","RMF",
                                                        "DSI","Precipitation")

# significant for SLA, height, leafN, DSI, Precip, rootDaim

ggPD_boot_test_2(all.data.map.dsi,predictor = "Precipitation",list.4.preds=all.data.dsi.prerun, col.line="#769370",
                 booted.preds=all.data.dsi.boot$function.preds, cex.line=2, col.ci="#769370",
                 alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
                 y.label = "Percent Cover Change")

ggPD_boot_test_2(all.data.map.dsi,predictor = "LeafN",list.4.preds=all.data.dsi.prerun, col.line="#769370",
               booted.preds=all.data.dsi.boot$function.preds, cex.line=2, col.ci="#769370",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")

ggPD_boot_test_2(all.data.map.dsi,predictor = "RootN",list.4.preds=all.data.dsi.prerun, col.line="#769370",
                 booted.preds=all.data.dsi.boot$function.preds, cex.line=2, col.ci="#769370",
                 alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
                 y.label = "Percent Cover Change")

ggPD_boot_test_2(all.data.map.dsi,predictor = "SLA",list.4.preds=all.data.dsi.prerun, col.line="#769370",
                 booted.preds=all.data.dsi.boot$function.preds, cex.line=2, col.ci="#769370",
                 alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
                 y.label = "Percent Cover Change")


# output individual plots as 3x3

ggPD(all.data.map.dsi, col.line="#769370", common.scale = FALSE, y.label = "Percent Cover Change")
# output all plots as one 8x8

# investigation of interactions
all.data.map.dsi$contributions$var = c("mean.MAP","leafN.final","rootN.final","SLA.final",
                                       "RTD.final","height.final","mean.DSI","rootDiam.final",
                                       "SRL.final","root.depth.final","RMF.final")

all.data.map.dsi$gbm.call$predictor.names = c("leafN.final","height.final","rootN.final","SLA.final","root.depth.final","rootDiam.final",
                                              "SRL.final","RTD.final","RMF.final","mean.DSI","mean.MAP")
colnames(all.data.map.dsi$gbm.call$dataframe)[2:12] = c("leafN.final","height.final","rootN.final","SLA.final","root.depth.final","rootDiam.final",
                                                        "SRL.final","RTD.final","RMF.final","mean.DSI","mean.MAP")

gbm.interactions(all.data.map.dsi)$interactions
ggInteract_list(all.data.map.dsi, index = T)
# 11 MAP x 1 leafN 6.54
# 11 MAP x 4 SLA 3.99
# 10 DSI x 1 leafN 1.13
# 11 MAP x 3 rootN 0.84
# 3 rootN x 2 height 0.57
# 10 DSI x 3 rootN 0.56


# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model
Interact_boot_dsi<-ggInteract_boot(c('mean.MAP','leafN.final'),c('mean.MAP','SLA.final'),
                                   c('mean.DSI','leafN.final'),c('mean.MAP','rootN.final'),
                                   c('rootN.final','height.final'),c('mean.DSI','rootN.final'),
                                        nboots = 500, data=all.data, predictors =  c(2:12), 
                                        response="cover.change",
                                        family = "gaussian", tc = 9, lr = 0.0001, bf= 0.75, global.env=F)

# Significance histogram p-value< 0.05
ggInteract_boot_hist(data = Interact_boot_dsi, column = 2,obs = 6.54) # significant
# 11 MAP x 1 leafN 6.54
ggInteract_boot_hist(data = Interact_boot_dsi, column = 3,obs = 3.99) # significant
# 11 MAP x 4 SLA 3.99
ggInteract_boot_hist(data = Interact_boot_dsi, column = 4,obs = 1.13) 
ggInteract_boot_hist(data = Interact_boot_dsi, column = 5,obs = 0.84) 
ggInteract_boot_hist(data = Interact_boot_dsi, column = 6,obs = 0.57)
ggInteract_boot_hist(data = Interact_boot_dsi, column = 6,obs = 0.56)

ggInteract_2D_test(gbm.object = all.data.map.dsi, x="precip",y="height.m",col.gradient = c("white","#769370"),
                   show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = F,
                   col.contour = "#254376",show.axis = T,legend = T, x.label = "Precipitation", y.label = "Height",
                   z.range = c(-0.81, 1.38), z.label = "% Cover Change", smooth = "average")

ggInteract_3D(gbm.object = all.data.map.dsi, x="precip",y="height.m")

ggInteract_2D_test(gbm.object = all.data.map.dsi, x="mean.DSI",y="height.m",col.gradient = c("white","#769370"),
                   show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = F,
                   col.contour = "#254376",show.axis = T,legend = T, x.label = "DSI", y.label = "Height",
                   z.range = c(-0.48, 1.57), z.label = "% Cover Change", smooth = "average")

ggInteract_3D(gbm.object = all.data.map.dsi, x="mean.DSI",y="height.m")

ggInteract_2D_test(gbm.object = all.data.map.dsi, x="rootDiam.mm",y="height.m",col.gradient = c("white","#769370"),
                   show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = F,
                   col.contour = "#254376",show.axis = T,legend = T, x.label = "Root Diameter", y.label = "Height",
                   z.range = c(-0.75, 1.35), z.label = "% Cover Change", smooth = "average")

ggInteract_3D(gbm.object = all.data.map.dsi,x="rootDiam.mm",y="height.m")

ggInteract_2D_test(gbm.object = all.data.map.dsi, x="leafN.mg.g",y="height.m",col.gradient = c("white","#769370"),
                   show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = F,
                   col.contour = "#254376",show.axis = T,legend = T, x.label = "LeafN", y.label = "Height",
                   z.range = c(-0.54, 1.64), z.label = "% Cover Change", smooth = "average")

ggInteract_3D(gbm.object = all.data.map.dsi, x="leafN.mg.g",y="height.m")

ggInteract_2D_test(gbm.object = all.data.map.dsi, x="precip",y="leafN.mg.g",col.gradient = c("white","#769370"),
                   show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = F,
                   col.contour = "#254376",show.axis = T,legend = T, x.label = "Precipitation", y.label = "LeafN",
                   z.range = c(-1.14, 0.11), z.label = "% Cover Change", smooth = "average")

ggInteract_3D(gbm.object = all.data.map.dsi, x="precip",y="leafN.mg.g", z.range = c(-1.14, 0.11))


#### annual data ####
annual.no.site.map.dsi = gbm.step(data=annual.data, gbm.x = c(2:12), gbm.y=1,
                                family = "gaussian", tree.complexity = 3, learning.rate = 0.0001,
                                bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(annual.no.site.map.dsi)
# 1000 trees Per.Expl = 1.91%

#saveRDS(annual.no.site.map.dsi, file = "./Formatted.Data/Revisions/annual.species.model.rds")

ggInfluence(annual.no.site.map.dsi)

annual.no.site.map.dsi$contributions$var = c("LeafN","RTD","Height",
                                             "Root Diameter","RootN","Precipitation",
                                             "DSI","SLA","SRL",
                                              "Rooting Depth")

annual.influence.plot = ggInfluence_test(annual.no.site.map.dsi, main = expression("Annuals (n = 120), R"^2*" = 27.79%"), 
                                    col.bar = c("#979461","#979461","#979461","gray70",
                                                "gray70","gray70","gray70","gray70","gray70","gray70"),
                                    col.signif = "#B50200")
# saved as pdf from device (4 x 5)

# get data to plot partial dependency plots
annual.no.site.map.dsi$contributions$var = c("leafN.mg.g","RTD.g.cm3","height.m",
                                             "rootDiam.mm","rootN.mg.g","precip",
                                             "mean.DSI","SLA_m2.kg","SRL_m.g",
                                             "root.depth_m")

annual.dsi.prerun<- plot.gbm.4list(annual.no.site.map.dsi)

annual.dsi.boot <- gbm.bootstrap.functions(annual.no.site.map.dsi, list.predictors=annual.dsi.prerun, n.reps=1000)

annual.no.site.map.dsi$contributions$var = c("LeafN","RTD","Height",
                                             "Root Diameter","RootN","Precipitation",
                                             "DSI","SLA","SRL",
                                             "Rooting Depth")

annual.no.site.map.dsi$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation","DSI")
colnames(annual.no.site.map.dsi$gbm.call$dataframe)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(annual.no.site.map.dsi$gbm.call$dataframe)[22:23] = c("Precipitation","DSI")

# significant for DSI: LeafN, RTD, height

ggPD_boot_test_2(annual.no.site.map.dsi,predictor = "LeafN",list.4.preds=annual.dsi.prerun, col.line="#979461",
          booted.preds=annual.dsi.boot$function.preds, cex.line=2, col.ci="#979461",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = TRUE,
          y.label = "Percent Cover Change")

ggPD_boot_test_2(annual.no.site.map.dsi,predictor = "RTD",list.4.preds=annual.dsi.prerun, col.line="#979461",
               booted.preds=annual.dsi.boot$function.preds, cex.line=2, col.ci="#979461",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = TRUE,
               y.label = "Percent Cover Change")

ggPD_boot_test_2(annual.no.site.map.dsi,predictor = "Height",list.4.preds=annual.dsi.prerun, col.line="#979461",
                 booted.preds=annual.dsi.boot$function.preds, cex.line=2, col.ci="#979461",
                 alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = TRUE,
                 y.label = "Percent Cover Change")


ggPD(annual.no.site.map.dsi, col.line="#979461", common.scale = FALSE, y.label = "Percent Cover Change")

# output all plots as one 8x8

gbm.interactions(annual.no.site.map.dsi)$interactions
ggInteract_list(annual.no.site.map.dsi, index = T)
# 10 DSI x 1 leafN 0.19
# 11 MAP x 1 leafN 0.05
# 11 MAP x 2 height 0.03
# 11 MAP x 6 rootDiam 0.02
# 3 rootN x 1 leafN 0.02
# 11 MAP x 8 RTD 0.01


# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model
Interact_boot_dsi<-ggInteract_boot(c('mean.DSI','leafN.final'),c('mean.MAP','leafN.final'),
                                   c('mean.MAP','height.final'),c('mean.MAP','rootDiam.final'),
                                   c('rootN.final','leafN.final'),c('mean.MAP','RTD.final'),
                                   nboots = 500, data=annual.data, predictors =  c(2:12), 
                                   response="cover.change",
                                   family = "gaussian", tc = 3, lr = 0.0001, bf= 0.50, global.env=F)

# Significance histogram p-value< 0.05
ggInteract_boot_hist(data = Interact_boot_dsi, column = 2,obs = 0.19)
ggInteract_boot_hist(data = Interact_boot_dsi, column = 3,obs = 0.05)
ggInteract_boot_hist(data = Interact_boot_dsi, column = 4,obs = 0.03) 
ggInteract_boot_hist(data = Interact_boot_dsi, column = 5,obs = 0.02) 
ggInteract_boot_hist(data = Interact_boot_dsi, column = 6,obs = 0.02)
ggInteract_boot_hist(data = Interact_boot_dsi, column = 6,obs = 0.01)



#### Perennial without woody ####
perennial.data.map.dsi = gbm.step(data=perennial.data, gbm.x = c(2:12), gbm.y=1,
                                family = "gaussian", tree.complexity = 9, learning.rate = 0.0001,
                                bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)


ggPerformance(perennial.data.map.dsi)
# 1000 trees Per.Expl = 2.36%

#saveRDS(perennial.data.map.dsi, file = "./Formatted.Data/Revisions/perennial.species.model.rds")

ggInfluence(perennial.data.map.dsi)

perennial.data.map.dsi$contributions$var = c("SLA","LeafN","DSI","Height",
                                         "SRL","Precipitation",
                                         "Root Diameter","RTD","Rooting Depth","RootN")

perennial.influence.plot = ggInfluence_test(perennial.data.map.dsi, main = expression("Perennials (n = 490), R"^2*" = 6.77%"), 
                                       col.bar = c("#F1C646","#F1C646","#F1C646","#F1C646",
                                                   "gray70","gray70","gray70","gray70","gray70","gray70"
                                       ), col.signif = "#B50200")

# get data to plot partial dependency plots

perennial.data.map.dsi$contributions$var = c("SLA_m2.kg","leafN.mg.g","mean.DSI","height.m",
                                             "SRL_m.g","precip",
                                             "rootDiam.mm","RTD.g.cm3","root.depth_m","rootN.mg.g")

perennial.dsi.prerun<- plot.gbm.4list(perennial.data.map.dsi)

perennial.dsi.boot <- gbm.bootstrap.functions(perennial.data.map.dsi, list.predictors=perennial.dsi.prerun, n.reps=1000)

perennial.data.map.dsi$contributions$var = c("SLA","LeafN","DSI","Height",
                                             "SRL","Precipitation",
                                             "Root Diameter","RTD","Rooting Depth","RootN")

perennial.data.map.dsi$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation","DSI")
colnames(perennial.data.map.dsi$gbm.call$dataframe)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(perennial.data.map.dsi$gbm.call$dataframe)[22:23] = c("Precipitation","DSI")

# significant for both: SLA, LeafN, DSI, Height

ggPD_boot_test_2(perennial.data.map.dsi,predictor = "SLA",list.4.preds=perennial.dsi.prerun, col.line="#F1C646",
               booted.preds=perennial.dsi.boot$function.preds, cex.line=2, col.ci="#F1C646",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = TRUE,
               y.label = "Percent Cover Change")

ggPD_boot_test_2(perennial.data.map.dsi,predictor = "LeafN",list.4.preds=perennial.dsi.prerun, col.line="#F1C646",
                 booted.preds=perennial.dsi.boot$function.preds, cex.line=2, col.ci="#F1C646",
                 alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = TRUE,
                 y.label = "Percent Cover Change")

ggPD_boot_test_2(perennial.data.map.dsi,predictor = "DSI",list.4.preds=perennial.dsi.prerun, col.line="#F1C646",
                 booted.preds=perennial.dsi.boot$function.preds, cex.line=2, col.ci="#F1C646",
                 alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = TRUE,
                 y.label = "Percent Cover Change")

ggPD_boot_test_2(perennial.data.map.dsi,predictor = "Height",list.4.preds=perennial.dsi.prerun, col.line="#F1C646",
                 booted.preds=perennial.dsi.boot$function.preds, cex.line=2, col.ci="#F1C646",
                 alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = TRUE,
                 y.label = "Percent Cover Change")

ggPD(perennial.data.map.dsi, col.line="#F1C646", common.scale = FALSE, y.label = "Percent Cover Change")

# output all plots as one 8x8

# investigation of interactions

perennial.data.map.dsi$contributions$var = c("SLA_m2.kg","leafN.mg.g","mean.DSI","height.m",
                                             "SRL_m.g","precip",
                                             "rootDiam.mm","RTD.g.cm3","root.depth_m","rootN.mg.g")

perennial.data.map.dsi$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm", "precip","mean.DSI")
colnames(perennial.data.map.dsi$gbm.call$dataframe)[10:17] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(perennial.data.map.dsi$gbm.call$dataframe)[22:23] = c("precip","mean.DSI")

gbm.interactions(perennial.data.map.dsi)$rank.list
ggInteract_list(perennial.data.map.dsi)
# 8 RTD x 7 SRL 0.25
# 11 MAP x 8 RTD 0.24
# 11 MAP x 1 leafN 0.21
# 11 MAP x 4 SLA 0.12
# 9 RMF x 1 leafN 0.12
# 11 MAP x 3 rootN 0.11

# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model

Interact_boot_perennial_dsi<-ggInteract_boot(c('RTD.final','SRL.final'),c('mean.MAP','RTD.final'),
                                             c('mean.MAP','leafN.final'),c('mean.MAP','SLA.final'),
                                         c('RMF.final','leafN.final'),c('mean.MAP','rootN.final'),
                                         nboots = 500,data=perennial.data, predictors =  c(2:12), 
                                         response="cover.change",
                                         family = "gaussian", tc = 9, lr = 0.0001, bf= 0.75, global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_perennial_dsi, column = 2,obs = 0.25) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial_dsi, column = 3,obs = 0.24) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial_dsi, column = 4,obs = 0.21) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial_dsi, column = 5,obs = 0.12) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial_dsi, column = 6,obs = 0.12) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial_dsi, column = 7,obs = 0.11) # non-significant

#### Grass ####
grass.map.dsi = gbm.step(data=grass, gbm.x = c(2:12), gbm.y=1,
                       family = "gaussian", tree.complexity = 1, learning.rate = 0.001,
                       bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 25)

ggPerformance(grass.map.dsi)
# 2400 trees Per.Expl = 12.23%

saveRDS(grass.map.dsi, file = "./Formatted.Data/grass.model.rds")

ggInfluence(grass.map.dsi)

grass.map.dsi$contributions$var = c("SLA","Height","DSI",
                                    "Precipitation", 
                                    "Root Diameter","RootN",
                                    "Rooting Depth","LeafN",
                                    "RTD","SRL")
                                
grass.influence.plot = ggInfluence_test(grass.map.dsi, main = expression("Grasses (n = 218), R"^2*" = 12.23%"), 
                                   col.bar = c("#6E687E","#6E687E","#6E687E","#6E687E",
                                               "gray70","gray70","gray70","gray70","gray70","gray70"
                                   ), col.signif = "#B50200")

# get data to plot partial dependency plots

grass.map.dsi$contributions$var = c("SLA_m2.kg","height.m","mean.DSI",
                                    "precip",
                                    "rootDiam.mm","rootN.mg.g",
                                    "root.depth_m","leafN.mg.g",
                                    "RTD.g.cm3","SRL_m.g")

grass.dsi.prerun<- plot.gbm.4list(grass.map.dsi)

grass.dsi.boot <- gbm.bootstrap.functions(grass.map.dsi, list.predictors=grass.dsi.prerun, n.reps=1000)

grass.map.dsi$contributions$var = c("SLA","Height","DSI",
                                    "Precipitation", 
                                    "Root Diameter","RootN",
                                    "Rooting Depth","LeafN",
                                    "RTD","SRL")

grass.map.dsi$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation","DSI")
colnames(grass.map.dsi$gbm.call$dataframe)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(grass.map.dsi$gbm.call$dataframe)[22:23] = c("Precipitation","DSI")

# significant for both: SLA, Height, DSI, Precip

ggPD_boot_test_2(grass.map.dsi,predictor = "SLA",list.4.preds=grass.dsi.prerun, col.line="#6E687E",
               booted.preds=grass.dsi.boot$function.preds, cex.line=2, col.ci="#6E687E",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")

ggPD_boot_test_2(grass.map.dsi,predictor = "Height",list.4.preds=grass.dsi.prerun, col.line="#6E687E",
               booted.preds=grass.dsi.boot$function.preds, cex.line=2, col.ci="#6E687E",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")

ggPD_boot_test_2(grass.map.dsi,predictor = "DSI",list.4.preds=grass.dsi.prerun, col.line="#6E687E",
               booted.preds=grass.dsi.boot$function.preds, cex.line=2, col.ci="#6E687E",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")

ggPD_boot_test_2(grass.map.dsi,predictor = "Precipitation",list.4.preds=grass.dsi.prerun, col.line="#6E687E",
                 booted.preds=grass.dsi.boot$function.preds, cex.line=2, col.ci="#6E687E",
                 alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
                 y.label = "Percent Cover Change")

ggPD(grass.map.dsi, col.line="#6E687E", common.scale = FALSE, y.label = "Percent Cover Change")

# output all plots as one 8x8

#### Forb ####
forb.map.dsi = gbm.step(data=forb, gbm.x = c(2:12), gbm.y=1,
                      family = "gaussian", tree.complexity = 8, learning.rate = 0.0005,
                      bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(forb.map.dsi)
# 1000 trees Per.Expl = 13.12%

#saveRDS(forb.map.dsi, file = "./Formatted.Data/Revisions/forb.species.model.rds")


ggInfluence(forb.map.dsi)

forb.map.dsi$contributions$var = c("LeafN","Root Diameter","Height",
                                   "SRL","DSI",
                                   "Precipitation","Rooting Depth","SLA",
                                   "RTD", "RootN")

forb.influence.plot = ggInfluence_test(forb.map.dsi, main = expression("Forbs (n = 325), R"^2*" = 3.57%"),
                                  col.bar = c("#F17236","#F17236","gray70","gray70",
                                              "gray70","gray70","gray70","gray70","gray70","gray70"
                                  ), col.signif = "#B50200")

# get data to plot partial dependency plots

forb.map.dsi$contributions$var = c("leafN.mg.g","rootDiam.mm","height.m",
                                   "SRL_m.g","mean.DSI",
                                   "precip","root.depth_m","SLA_m2.kg",
                                   "RTD.g.cm3", "rootN.mg.g")

forb.dsi.prerun<- plot.gbm.4list(forb.map.dsi)

forb.dsi.boot <- gbm.bootstrap.functions(forb.map.dsi, list.predictors=forb.dsi.prerun, n.reps=1000)

forb.map.dsi$contributions$var = c("LeafN","Root Diameter","Height",
                                   "SRL","DSI",
                                   "Precipitation","Rooting Depth","SLA",
                                   "RTD", "RootN")

forb.map.dsi$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation","DSI")
colnames(forb.map.dsi$gbm.call$dataframe)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(forb.map.dsi$gbm.call$dataframe)[22:23] = c("Precipitation","DSI")

# significant: LeafN, Diameter

ggPD_boot_test_2(forb.map.dsi,predictor = "LeafN",list.4.preds=forb.dsi.prerun, col.line="#F17236",
               booted.preds=forb.dsi.boot$function.preds, cex.line=2, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")

ggPD_boot_test_2(forb.map.dsi,predictor = "Root Diameter",list.4.preds=forb.dsi.prerun, col.line="#F17236",
               booted.preds=forb.dsi.boot$function.preds, cex.line=2, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")


ggPD(forb.map.dsi, col.line="#F17236", common.scale = FALSE, y.label = "Percent Cover Change")

# output all plots as one 8x8

gbm.interactions(forb.map.dsi)$interactions
ggInteract_list(forb.map.dsi, index = T)
# 8 RTD x 7 SRL 0.19
# 2 height x 1 leafN 0.05
# 8 RTD x 1 leafN 0.03
# 10 DSI x 1 leafN 0.02
# 11 MAP x 4 SLA 0.02
# 11 MAP x 1 leafN 0.01


# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model
Interact_boot_dsi<-ggInteract_boot(c('RTD.final','SRL.final'),c('height.final','leafN.final'),
                                   c('RTD.final','leafN.final'),c('mean.DSI','leafN.final'),
                                   c('mean.MAP','SLA.final'),c('mean.MAP','leafN.final'),
                                   nboots = 500, data=forb, predictors =  c(2:12), 
                                   response="cover.change",
                                   family = "gaussian", tc = 8, lr = 0.0005, bf= 0.75, global.env=F)

# Significance histogram p-value< 0.05
ggInteract_boot_hist(data = Interact_boot_dsi, column = 2,obs = 0.19)
ggInteract_boot_hist(data = Interact_boot_dsi, column = 3,obs = 0.05)
ggInteract_boot_hist(data = Interact_boot_dsi, column = 4,obs = 0.03) 
ggInteract_boot_hist(data = Interact_boot_dsi, column = 5,obs = 0.02) 
ggInteract_boot_hist(data = Interact_boot_dsi, column = 6,obs = 0.02)
ggInteract_boot_hist(data = Interact_boot_dsi, column = 6,obs = 0.01)


#### Grass Perennials ####
grass.perennial.map.dsi = gbm.step(data=grass.perennial, gbm.x = c(2:12), gbm.y=1,
                                 family = "gaussian", tree.complexity = 2, learning.rate = 0.0001,
                                 bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(grass.perennial.map.dsi)
# 1000 trees Per.Expl = 1.38%

#saveRDS(grass.perennial.map.dsi, file = "./Formatted.Data/Revisions/grass.perennial.model.rds")

ggInfluence(grass.perennial.map.dsi)

grass.perennial.map.dsi$contributions$var = c("SLA", "Height", 
                                              "Root Diameter","DSI",
                                              "Precipitation","Rooting Depth","LeafN", 
                                              "SRL","RootN","RTD")

grass.perennial.influence.plot = ggInfluence_test(grass.perennial.map.dsi, main = expression("Perennnial Grasses (n = 169), R"^2*" = 19.51%"), 
                                        col.bar = c("#6089B5","#6089B5","#6089B5","#6089B5",
                                                    "gray70","gray70","gray70","gray70","gray70","gray70"
                                        ), col.signif = "#B50200")

# saved as pdf from device (4 x 6)

# get data to plot partial dependency plots

grass.perennial.map.dsi$contributions$var = c("SLA_m2.kg","height.m","rootDiam.mm",
                                              "mean.DSI", "precip","root.depth_m","leafN.mg.g", 
                                    "SRL_m.g","rootN.mg.g","RTD.g.cm3")
                                    

grass.perennial.dsi.prerun<- plot.gbm.4list(grass.perennial.map.dsi)

grass.perennial.dsi.boot <- gbm.bootstrap.functions(grass.perennial.map.dsi, list.predictors=grass.perennial.dsi.prerun, n.reps=1000)

grass.perennial.map.dsi$contributions$var = c("SLA", "Height", 
                                              "Root Diameter","DSI",
                                              "Precipitation","Rooting Depth","LeafN", 
                                              "SRL","RootN","RTD")

grass.perennial.map.dsi$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation","DSI")
colnames(grass.perennial.map.dsi$gbm.call$dataframe)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(grass.perennial.map.dsi$gbm.call$dataframe)[22:23] = c("Precipitation","DSI")

# significant for both: SLA, Height, Diam, DSI

ggPD_boot_test_2(grass.perennial.map.dsi,predictor = "SLA",list.4.preds=grass.perennial.dsi.prerun, col.line="#6089B5",
                 booted.preds=grass.perennial.dsi.boot$function.preds, cex.line=2, col.ci="#6089B5",
                 alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
                 y.label = "Percent Cover Change")

ggPD_boot_test_2(grass.perennial.map.dsi,predictor = "Height",list.4.preds=grass.perennial.dsi.prerun, col.line="#6089B5",
                 booted.preds=grass.perennial.dsi.boot$function.preds, cex.line=2, col.ci="#6089B5",
                 alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
                 y.label = "Percent Cover Change")

ggPD_boot_test_2(grass.perennial.map.dsi,predictor = "Root Diameter",list.4.preds=grass.perennial.dsi.prerun, col.line="#6089B5",
                 booted.preds=grass.perennial.dsi.boot$function.preds, cex.line=2, col.ci="#6089B5",
                 alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
                 y.label = "Percent Cover Change")

ggPD_boot_test_2(grass.perennial.map.dsi,predictor = "DSI",list.4.preds=grass.perennial.dsi.prerun, col.line="#6089B5",
                 booted.preds=grass.perennial.dsi.boot$function.preds, cex.line=2, col.ci="#6089B5",
                 alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
                 y.label = "Percent Cover Change")

# output individual plots as 3x3

ggPD(grass.perennial.map.dsi, col.line="#6089B5", common.scale = FALSE, y.label = "Percent Cover Change")

# output all plots as one 8x8

gbm.interactions(grass.perennial.map.dsi)$interactions
ggInteract_list(grass.perennial.map.dsi, index = T)
# 11 MAP x 7 SRL 0.06



# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model
Interact_boot_dsi<-ggInteract_boot(c('mean.MAP','SRL.final'),
                                   nboots = 500, data=grass.perennial, predictors =  c(2:12), 
                                   response="cover.change",
                                   family = "gaussian", tc = 2, lr = 0.0001, bf= 0.50, global.env=F)

# Significance histogram p-value< 0.05
ggInteract_boot_hist(data = Interact_boot_dsi, column = 2,obs = 0.06)



#### Forb Perennial ####
forb.perennial.map.dsi = gbm.step(data=forb.perennial, gbm.x = c(2:12), gbm.y=1,
                                  family = "gaussian", tree.complexity = 6, learning.rate = 0.0001,
                                  bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(forb.perennial.map.dsi)
# 1000 trees Per.Expl = 3.01%

#saveRDS(forb.perennial.map.dsi, file = "./Formatted.Data/Revisions/forb.perennial.model.rds")

ggInfluence(forb.perennial.map.dsi)

forb.perennial.map.dsi$contributions$var = c("LeafN","Root Diameter","SRL",
                                             "Rooting Depth","Precipitation",
                                             "DSI","Height","RTD","RootN","SLA")

forb.perennial.influence.plot = ggInfluence_test(forb.perennial.map.dsi, main = expression("Perennial Forbs (n = 248), R"^2*" = 2.99%"),
                                                 col.bar = c("black","black","black","gray70",
                                                             "gray70","gray70","gray70","gray70","gray70","gray70"
                                                 ), col.signif = "#B50200")

# get data to plot partial dependency plots

forb.perennial.map.dsi$contributions$var = c("leafN.mg.g","rootDiam.mm","SRL_m.g",
                                             "root.depth_m","precip",
                                             "mean.DSI",
                                             "height.m","RTD.g.cm3", 
                                             "rootN.mg.g","SLA_m2.kg")

forb.perennial.dsi.prerun<- plot.gbm.4list(forb.perennial.map.dsi)

forb.perennial.dsi.boot <- gbm.bootstrap.functions(forb.perennial.map.dsi, list.predictors=forb.perennial.dsi.prerun, n.reps=1000)

forb.perennial.map.dsi$contributions$var = c("LeafN","Root Diameter","SRL",
                                             "Rooting Depth","Precipitation",
                                             "DSI","Height","RTD","RootN","SLA")

forb.perennial.map.dsi$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation","DSI")
colnames(forb.perennial.map.dsi$gbm.call$dataframe)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(forb.perennial.map.dsi$gbm.call$dataframe)[22:23] = c("Precipitation","DSI")

# significant for DSI: leafN, Root diameter, SRL

ggPD_boot_test_2(forb.perennial.map.dsi,predictor = "LeafN",list.4.preds=forb.perennial.dsi.prerun, col.line="black",
                 booted.preds=forb.perennial.dsi.boot$function.preds, cex.line=2, col.ci="black",
                 alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
                 y.label = "Percent Cover Change")

ggPD_boot_test_2(forb.perennial.map.dsi,predictor = "Root Diameter",list.4.preds=forb.perennial.dsi.prerun, col.line="black",
                 booted.preds=forb.perennial.dsi.boot$function.preds, cex.line=2, col.ci="black",
                 alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
                 y.label = "Percent Cover Change")

ggPD_boot_test_2(forb.perennial.map.dsi,predictor = "SRL",list.4.preds=forb.perennial.dsi.prerun, col.line="black",
                 booted.preds=forb.perennial.dsi.boot$function.preds, cex.line=2, col.ci="black",
                 alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
                 y.label = "Percent Cover Change")

ggPD(forb.perennial.map.dsi, col.line="black", common.scale = FALSE, y.label = "Percent Cover Change")

# output all plots as one 8x8

gbm.interactions(forb.perennial.map.dsi)$interactions
ggInteract_list(forb.perennial.map.dsi, index = T)
# 4 SLA x 2 height 0.11
# 10 DSI x 7 SRL 0.10
# 10 DSI x 1 leafN 0.07
# 5 depth x 4 SLA 0.07
# 11 MAP x 2 height 0.06
# 11 MAP x 8 RTD 0.05


# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model
Interact_boot_dsi<-ggInteract_boot(c('SLA.final','height.final'),c('mean.DSI','SRL.final'),
                                   c('mean.DSI','leafN.final'),c('root.depth.final','SLA.final'),
                                   c('mean.MAP','height.final'),c('mean.MAP','RTD.final'),
                                   nboots = 500, data=forb.perennial, predictors =  c(2:12), 
                                   response="cover.change",
                                   family = "gaussian", tc = 6, lr = 0.0001, bf= 0.75, global.env=F)

# Significance histogram p-value< 0.05
ggInteract_boot_hist(data = Interact_boot_dsi, column = 2,obs = 0.11)
ggInteract_boot_hist(data = Interact_boot_dsi, column = 3,obs = 0.10)
ggInteract_boot_hist(data = Interact_boot_dsi, column = 4,obs = 0.07)
ggInteract_boot_hist(data = Interact_boot_dsi, column = 5,obs = 0.07)
ggInteract_boot_hist(data = Interact_boot_dsi, column = 6,obs = 0.06)
ggInteract_boot_hist(data = Interact_boot_dsi, column = 7,obs = 0.05)


#### Table S1 ####
mean(all.data$leafN.mg.g, na.rm = TRUE)
mean(all.data$height.m, na.rm = TRUE)
mean(all.data$rootN.mg.g, na.rm = TRUE)
mean(all.data$SLA_m2.kg, na.rm = TRUE)
mean(all.data$root.depth_m, na.rm = TRUE)
mean(all.data$RTD.g.cm3, na.rm = TRUE)
mean(all.data$SRL_m.g, na.rm = TRUE)
mean(all.data$rootDiam.mm, na.rm = TRUE)
mean(all.data$precip, na.rm = TRUE)
mean(all.data$mean.DSI, na.rm = TRUE)
mean(all.data$cover.change)

sd(all.data$leafN.mg.g, na.rm = TRUE)
sd(all.data$height.m, na.rm = TRUE)
sd(all.data$rootN.mg.g, na.rm = TRUE)
sd(all.data$SLA_m2.kg, na.rm = TRUE)
sd(all.data$root.depth_m, na.rm = TRUE)
sd(all.data$RTD.g.cm3, na.rm = TRUE)
sd(all.data$SRL_m.g, na.rm = TRUE)
sd(all.data$rootDiam.mm, na.rm = TRUE)
sd(all.data$precip, na.rm = TRUE)
sd(all.data$mean.DSI, na.rm = TRUE)
sd(all.data$cover.change)

range(all.data$leafN.mg.g, na.rm = TRUE)
range(all.data$height.m, na.rm = TRUE)
range(all.data$rootN.mg.g, na.rm = TRUE)
range(all.data$SLA_m2.kg, na.rm = TRUE)
range(all.data$root.depth_m, na.rm = TRUE)
range(all.data$RTD.g.cm3, na.rm = TRUE)
range(all.data$SRL_m.g, na.rm = TRUE)
range(all.data$rootDiam.mm, na.rm = TRUE)
range(all.data$precip, na.rm = TRUE)
range(all.data$mean.DSI, na.rm = TRUE)
range(all.data$cover.change)

sum(is.na(all.data$leafN.mg.g))/629*100
sum(is.na(all.data$height.m))/629*100
sum(is.na(all.data$rootN.mg.g))/629*100
sum(is.na(all.data$SLA_m2.kg))/629*100
sum(is.na(all.data$root.depth_m))/629*100
sum(is.na(all.data$RTD.g.cm3))/629*100
sum(is.na(all.data$SRL_m.g))/629*100
sum(is.na(all.data$rootDiam.mm))/629*100
sum(is.na(all.data$precip))/629*100
sum(is.na(all.data$mean.DSI))/629*100
sum(is.na(all.data$cover.change))/629*100

sd(all.data$leafN.mg.g, na.rm = TRUE)/mean(all.data$leafN.mg.g, na.rm = TRUE)
sd(all.data$height.m, na.rm = TRUE)/mean(all.data$height.m, na.rm = TRUE)
sd(all.data$rootN.mg.g, na.rm = TRUE)/mean(all.data$rootN.mg.g, na.rm = TRUE)
sd(all.data$SLA_m2.kg, na.rm = TRUE)/mean(all.data$SLA_m2.kg, na.rm = TRUE)
sd(all.data$root.depth_m, na.rm = TRUE)/mean(all.data$root.depth_m, na.rm = TRUE)
sd(all.data$RTD.g.cm3, na.rm = TRUE)/mean(all.data$RTD.g.cm3, na.rm = TRUE)
sd(all.data$SRL_m.g, na.rm = TRUE)/mean(all.data$SRL_m.g, na.rm = TRUE)
sd(all.data$rootDiam.mm, na.rm = TRUE)/mean(all.data$rootDiam.mm, na.rm = TRUE)
sd(all.data$precip, na.rm = TRUE)/mean(all.data$precip, na.rm = TRUE)
sd(all.data$mean.DSI, na.rm = TRUE)/mean(all.data$mean.DSI, na.rm = TRUE)
sd(all.data$cover.change)/mean(all.data$cover.change)

mean(annual.data$leafN.mg.g, na.rm = TRUE)
mean(annual.data$height.m, na.rm = TRUE)
mean(annual.data$rootN.mg.g, na.rm = TRUE)
mean(annual.data$SLA_m2.kg, na.rm = TRUE)
mean(annual.data$root.depth_m, na.rm = TRUE)
mean(annual.data$RTD.g.cm3, na.rm = TRUE)
mean(annual.data$SRL_m.g, na.rm = TRUE)
mean(annual.data$rootDiam.mm, na.rm = TRUE)
mean(annual.data$precip, na.rm = TRUE)
mean(annual.data$mean.DSI, na.rm = TRUE)
mean(annual.data$cover.change)

sd(annual.data$leafN.mg.g, na.rm = TRUE)
sd(annual.data$height.m, na.rm = TRUE)
sd(annual.data$rootN.mg.g, na.rm = TRUE)
sd(annual.data$SLA_m2.kg, na.rm = TRUE)
sd(annual.data$root.depth_m, na.rm = TRUE)
sd(annual.data$RTD.g.cm3, na.rm = TRUE)
sd(annual.data$SRL_m.g, na.rm = TRUE)
sd(annual.data$rootDiam.mm, na.rm = TRUE)
sd(annual.data$precip, na.rm = TRUE)
sd(annual.data$mean.DSI, na.rm = TRUE)
sd(annual.data$cover.change)

range(annual.data$leafN.mg.g, na.rm = TRUE)
range(annual.data$height.m, na.rm = TRUE)
range(annual.data$rootN.mg.g, na.rm = TRUE)
range(annual.data$SLA_m2.kg, na.rm = TRUE)
range(annual.data$root.depth_m, na.rm = TRUE)
range(annual.data$RTD.g.cm3, na.rm = TRUE)
range(annual.data$SRL_m.g, na.rm = TRUE)
range(annual.data$rootDiam.mm, na.rm = TRUE)
range(annual.data$precip, na.rm = TRUE)
range(annual.data$mean.DSI, na.rm = TRUE)
range(annual.data$cover.change)

mean(perennial.data$leafN.mg.g, na.rm = TRUE)
mean(perennial.data$height.m, na.rm = TRUE)
mean(perennial.data$rootN.mg.g, na.rm = TRUE)
mean(perennial.data$SLA_m2.kg, na.rm = TRUE)
mean(perennial.data$root.depth_m, na.rm = TRUE)
mean(perennial.data$RTD.g.cm3, na.rm = TRUE)
mean(perennial.data$SRL_m.g, na.rm = TRUE)
mean(perennial.data$rootDiam.mm, na.rm = TRUE)
mean(perennial.data$precip, na.rm = TRUE)
mean(perennial.data$mean.DSI, na.rm = TRUE)
mean(perennial.data$cover.change)

sd(perennial.data$leafN.mg.g, na.rm = TRUE)
sd(perennial.data$height.m, na.rm = TRUE)
sd(perennial.data$rootN.mg.g, na.rm = TRUE)
sd(perennial.data$SLA_m2.kg, na.rm = TRUE)
sd(perennial.data$root.depth_m, na.rm = TRUE)
sd(perennial.data$RTD.g.cm3, na.rm = TRUE)
sd(perennial.data$SRL_m.g, na.rm = TRUE)
sd(perennial.data$rootDiam.mm, na.rm = TRUE)
sd(perennial.data$precip, na.rm = TRUE)
sd(perennial.data$mean.DSI, na.rm = TRUE)
sd(perennial.data$cover.change)

range(perennial.data$leafN.mg.g, na.rm = TRUE)
range(perennial.data$height.m, na.rm = TRUE)
range(perennial.data$rootN.mg.g, na.rm = TRUE)
range(perennial.data$SLA_m2.kg, na.rm = TRUE)
range(perennial.data$root.depth_m, na.rm = TRUE)
range(perennial.data$RTD.g.cm3, na.rm = TRUE)
range(perennial.data$SRL_m.g, na.rm = TRUE)
range(perennial.data$rootDiam.mm, na.rm = TRUE)
range(perennial.data$precip, na.rm = TRUE)
range(perennial.data$mean.DSI, na.rm = TRUE)
range(perennial.data$cover.change)

mean(grass$leafN.mg.g, na.rm = TRUE)
mean(grass$height.m, na.rm = TRUE)
mean(grass$rootN.mg.g, na.rm = TRUE)
mean(grass$SLA_m2.kg, na.rm = TRUE)
mean(grass$root.depth_m, na.rm = TRUE)
mean(grass$RTD.g.cm3, na.rm = TRUE)
mean(grass$SRL_m.g, na.rm = TRUE)
mean(grass$rootDiam.mm, na.rm = TRUE)
mean(grass$precip, na.rm = TRUE)
mean(grass$mean.DSI, na.rm = TRUE)
mean(grass$cover.change)

sd(grass$leafN.mg.g, na.rm = TRUE)
sd(grass$height.m, na.rm = TRUE)
sd(grass$rootN.mg.g, na.rm = TRUE)
sd(grass$SLA_m2.kg, na.rm = TRUE)
sd(grass$root.depth_m, na.rm = TRUE)
sd(grass$RTD.g.cm3, na.rm = TRUE)
sd(grass$SRL_m.g, na.rm = TRUE)
sd(grass$rootDiam.mm, na.rm = TRUE)
sd(grass$precip, na.rm = TRUE)
sd(grass$mean.DSI, na.rm = TRUE)
sd(grass$cover.change)

range(grass$leafN.mg.g, na.rm = TRUE)
range(grass$height.m, na.rm = TRUE)
range(grass$rootN.mg.g, na.rm = TRUE)
range(grass$SLA_m2.kg, na.rm = TRUE)
range(grass$root.depth_m, na.rm = TRUE)
range(grass$RTD.g.cm3, na.rm = TRUE)
range(grass$SRL_m.g, na.rm = TRUE)
range(grass$rootDiam.mm, na.rm = TRUE)
range(grass$precip, na.rm = TRUE)
range(grass$mean.DSI, na.rm = TRUE)
range(grass$cover.change)

mean(forb$leafN.mg.g, na.rm = TRUE)
mean(forb$height.m, na.rm = TRUE)
mean(forb$rootN.mg.g, na.rm = TRUE)
mean(forb$SLA_m2.kg, na.rm = TRUE)
mean(forb$root.depth_m, na.rm = TRUE)
mean(forb$RTD.g.cm3, na.rm = TRUE)
mean(forb$SRL_m.g, na.rm = TRUE)
mean(forb$rootDiam.mm, na.rm = TRUE)
mean(forb$precip, na.rm = TRUE)
mean(forb$mean.DSI, na.rm = TRUE)
mean(forb$cover.change)

sd(forb$leafN.mg.g, na.rm = TRUE)
sd(forb$height.m, na.rm = TRUE)
sd(forb$rootN.mg.g, na.rm = TRUE)
sd(forb$SLA_m2.kg, na.rm = TRUE)
sd(forb$root.depth_m, na.rm = TRUE)
sd(forb$RTD.g.cm3, na.rm = TRUE)
sd(forb$SRL_m.g, na.rm = TRUE)
sd(forb$rootDiam.mm, na.rm = TRUE)
sd(forb$precip, na.rm = TRUE)
sd(forb$mean.DSI, na.rm = TRUE)
sd(forb$cover.change)

range(forb$leafN.mg.g, na.rm = TRUE)
range(forb$height.m, na.rm = TRUE)
range(forb$rootN.mg.g, na.rm = TRUE)
range(forb$SLA_m2.kg, na.rm = TRUE)
range(forb$root.depth_m, na.rm = TRUE)
range(forb$RTD.g.cm3, na.rm = TRUE)
range(forb$SRL_m.g, na.rm = TRUE)
range(forb$rootDiam.mm, na.rm = TRUE)
range(forb$precip, na.rm = TRUE)
range(forb$mean.DSI, na.rm = TRUE)
range(forb$cover.change)

mean(grass.annual$leafN.mg.g, na.rm = TRUE)
mean(grass.annual$height.m, na.rm = TRUE)
mean(grass.annual$rootN.mg.g, na.rm = TRUE)
mean(grass.annual$SLA_m2.kg, na.rm = TRUE)
mean(grass.annual$root.depth_m, na.rm = TRUE)
mean(grass.annual$RTD.g.cm3, na.rm = TRUE)
mean(grass.annual$SRL_m.g, na.rm = TRUE)
mean(grass.annual$rootDiam.mm, na.rm = TRUE)
mean(grass.annual$precip, na.rm = TRUE)
mean(grass.annual$mean.DSI, na.rm = TRUE)
mean(grass.annual$cover.change)

sd(grass.annual$leafN.mg.g, na.rm = TRUE)
sd(grass.annual$height.m, na.rm = TRUE)
sd(grass.annual$rootN.mg.g, na.rm = TRUE)
sd(grass.annual$SLA_m2.kg, na.rm = TRUE)
sd(grass.annual$root.depth_m, na.rm = TRUE)
sd(grass.annual$RTD.g.cm3, na.rm = TRUE)
sd(grass.annual$SRL_m.g, na.rm = TRUE)
sd(grass.annual$rootDiam.mm, na.rm = TRUE)
sd(grass.annual$precip, na.rm = TRUE)
sd(grass.annual$mean.DSI, na.rm = TRUE)
sd(grass.annual$cover.change)

range(grass.annual$leafN.mg.g, na.rm = TRUE)
range(grass.annual$height.m, na.rm = TRUE)
range(grass.annual$rootN.mg.g, na.rm = TRUE)
range(grass.annual$SLA_m2.kg, na.rm = TRUE)
range(grass.annual$root.depth_m, na.rm = TRUE)
range(grass.annual$RTD.g.cm3, na.rm = TRUE)
range(grass.annual$SRL_m.g, na.rm = TRUE)
range(grass.annual$rootDiam.mm, na.rm = TRUE)
range(grass.annual$precip, na.rm = TRUE)
range(grass.annual$mean.DSI, na.rm = TRUE)
range(grass.annual$cover.change)

mean(grass.perennial$leafN.mg.g, na.rm = TRUE)
mean(grass.perennial$height.m, na.rm = TRUE)
mean(grass.perennial$rootN.mg.g, na.rm = TRUE)
mean(grass.perennial$SLA_m2.kg, na.rm = TRUE)
mean(grass.perennial$root.depth_m, na.rm = TRUE)
mean(grass.perennial$RTD.g.cm3, na.rm = TRUE)
mean(grass.perennial$SRL_m.g, na.rm = TRUE)
mean(grass.perennial$rootDiam.mm, na.rm = TRUE)
mean(grass.perennial$precip, na.rm = TRUE)
mean(grass.perennial$mean.DSI, na.rm = TRUE)
mean(grass.perennial$cover.change)

sd(grass.perennial$leafN.mg.g, na.rm = TRUE)
sd(grass.perennial$height.m, na.rm = TRUE)
sd(grass.perennial$rootN.mg.g, na.rm = TRUE)
sd(grass.perennial$SLA_m2.kg, na.rm = TRUE)
sd(grass.perennial$root.depth_m, na.rm = TRUE)
sd(grass.perennial$RTD.g.cm3, na.rm = TRUE)
sd(grass.perennial$SRL_m.g, na.rm = TRUE)
sd(grass.perennial$rootDiam.mm, na.rm = TRUE)
sd(grass.perennial$precip, na.rm = TRUE)
sd(grass.perennial$mean.DSI, na.rm = TRUE)
sd(grass.perennial$cover.change)

range(grass.perennial$leafN.mg.g, na.rm = TRUE)
range(grass.perennial$height.m, na.rm = TRUE)
range(grass.perennial$rootN.mg.g, na.rm = TRUE)
range(grass.perennial$SLA_m2.kg, na.rm = TRUE)
range(grass.perennial$root.depth_m, na.rm = TRUE)
range(grass.perennial$RTD.g.cm3, na.rm = TRUE)
range(grass.perennial$SRL_m.g, na.rm = TRUE)
range(grass.perennial$rootDiam.mm, na.rm = TRUE)
range(grass.perennial$precip, na.rm = TRUE)
range(grass.perennial$mean.DSI, na.rm = TRUE)
range(grass.perennial$cover.change)

mean(forb.annual$leafN.mg.g, na.rm = TRUE)
mean(forb.annual$height.m, na.rm = TRUE)
mean(forb.annual$rootN.mg.g, na.rm = TRUE)
mean(forb.annual$SLA_m2.kg, na.rm = TRUE)
mean(forb.annual$root.depth_m, na.rm = TRUE)
mean(forb.annual$RTD.g.cm3, na.rm = TRUE)
mean(forb.annual$SRL_m.g, na.rm = TRUE)
mean(forb.annual$rootDiam.mm, na.rm = TRUE)
mean(forb.annual$precip, na.rm = TRUE)
mean(forb.annual$mean.DSI, na.rm = TRUE)
mean(forb.annual$cover.change)

sd(forb.annual$leafN.mg.g, na.rm = TRUE)
sd(forb.annual$height.m, na.rm = TRUE)
sd(forb.annual$rootN.mg.g, na.rm = TRUE)
sd(forb.annual$SLA_m2.kg, na.rm = TRUE)
sd(forb.annual$root.depth_m, na.rm = TRUE)
sd(forb.annual$RTD.g.cm3, na.rm = TRUE)
sd(forb.annual$SRL_m.g, na.rm = TRUE)
sd(forb.annual$rootDiam.mm, na.rm = TRUE)
sd(forb.annual$precip, na.rm = TRUE)
sd(forb.annual$mean.DSI, na.rm = TRUE)
sd(forb.annual$cover.change)

range(forb.annual$leafN.mg.g, na.rm = TRUE)
range(forb.annual$height.m, na.rm = TRUE)
range(forb.annual$rootN.mg.g, na.rm = TRUE)
range(forb.annual$SLA_m2.kg, na.rm = TRUE)
range(forb.annual$root.depth_m, na.rm = TRUE)
range(forb.annual$RTD.g.cm3, na.rm = TRUE)
range(forb.annual$SRL_m.g, na.rm = TRUE)
range(forb.annual$rootDiam.mm, na.rm = TRUE)
range(forb.annual$precip, na.rm = TRUE)
range(forb.annual$mean.DSI, na.rm = TRUE)
range(forb.annual$cover.change)

mean(forb.perennial$leafN.mg.g, na.rm = TRUE)
mean(forb.perennial$height.m, na.rm = TRUE)
mean(forb.perennial$rootN.mg.g, na.rm = TRUE)
mean(forb.perennial$SLA_m2.kg, na.rm = TRUE)
mean(forb.perennial$root.depth_m, na.rm = TRUE)
mean(forb.perennial$RTD.g.cm3, na.rm = TRUE)
mean(forb.perennial$SRL_m.g, na.rm = TRUE)
mean(forb.perennial$rootDiam.mm, na.rm = TRUE)
mean(forb.perennial$precip, na.rm = TRUE)
mean(forb.perennial$mean.DSI, na.rm = TRUE)
mean(forb.perennial$cover.change)

sd(forb.perennial$leafN.mg.g, na.rm = TRUE)
sd(forb.perennial$height.m, na.rm = TRUE)
sd(forb.perennial$rootN.mg.g, na.rm = TRUE)
sd(forb.perennial$SLA_m2.kg, na.rm = TRUE)
sd(forb.perennial$root.depth_m, na.rm = TRUE)
sd(forb.perennial$RTD.g.cm3, na.rm = TRUE)
sd(forb.perennial$SRL_m.g, na.rm = TRUE)
sd(forb.perennial$rootDiam.mm, na.rm = TRUE)
sd(forb.perennial$precip, na.rm = TRUE)
sd(forb.perennial$mean.DSI, na.rm = TRUE)
sd(forb.perennial$cover.change)

range(forb.perennial$leafN.mg.g, na.rm = TRUE)
range(forb.perennial$height.m, na.rm = TRUE)
range(forb.perennial$rootN.mg.g, na.rm = TRUE)
range(forb.perennial$SLA_m2.kg, na.rm = TRUE)
range(forb.perennial$root.depth_m, na.rm = TRUE)
range(forb.perennial$RTD.g.cm3, na.rm = TRUE)
range(forb.perennial$SRL_m.g, na.rm = TRUE)
range(forb.perennial$rootDiam.mm, na.rm = TRUE)
range(forb.perennial$precip, na.rm = TRUE)
range(forb.perennial$mean.DSI, na.rm = TRUE)
range(forb.perennial$cover.change)

table(annual.data$functional_group)

table(perennial.data$functional_group)
# 230 forbs, 180 grasses
table(annual.data$functional_group)
# 76 forbs, 39 grasses


