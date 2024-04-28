# Script to peform analyses where % cover calculated using BACI method
# drought after-drought before) - (control after - control before)

# Use the same data files just removing outlier cover change tested for in All Data

library(tidyverse)
library(dismo)
library(gbm)
library(ggBRT)
library(multcomp)
library(cowplot)
library(corrplot)
library(treezy)

# edited functions from ggBRT needed for plotting
# ggPD_boot_edit.R
# ggInfluence_edit.R
# ggInteract_2D_edit.R

#### read in all the data frames needs for the analyses ####

# all data
all.data = read.csv("./Formatted.Data/BACI.all.data.csv", row.names = 1) # 1234 data points
# annual data
annual.data = read.csv("./Formatted.Data/BACI.annual.data.csv", row.names = 1) # 292 data points
# all perennial data
perennial.data = read.csv("./Formatted.Data/BACI.perennial.data.csv", row.names = 1) # 890 data points
# grass
grass = read.csv("./Formatted.Data/BACI.grass.csv", row.names = 1) # 392 data points
# forbs
forb = read.csv("./Formatted.Data/BACI.forb.csv", row.names = 1) # 675 data points
# grass perennial

#### Testing for outliers in all data ####

hist(all.data$cover.change)
boxplot(all.data$cover.change)
mean = mean(all.data$cover.change, na.rm = TRUE)
std = sd(all.data$cover.change, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(all.data$cover.change[which(all.data$cover.change <Tmin | all.data$cover.change > Tmax)])
# removed cover.change -79.47619 to -20.58333
# removed cover.change 20.00000 to 43.00000

all.data.otl.rm = subset(all.data, all.data$cover.change >-20.58333 & all.data$cover.change < 20.00000)
# lose 30 data points, lose 13 species
annual.data.otl.rm = subset(annual.data, annual.data$cover.change >-20.58333 & annual.data$cover.change < 20.00000)
# lose 8 data points, lose 1 species
perennial.data.otl.rm = subset(perennial.data, perennial.data$cover.change >-20.58333 & perennial.data$cover.change < 20.00000)
# lose 21 data points, lose 11 species
grass.data.otl.rm = subset(grass, grass$cover.change >-20.58333 & grass$cover.change < 20.00000)
# lose 19 data points, lose 7 species
forb.data.otl.rm = subset(forb, forb$cover.change >-20.58333 & forb$cover.change < 20.00000)
# lose 9 data points, lose 4 species

#### change site code to numeric, continuous vector ####
all.data.otl.rm$site.id = as.numeric(as.factor(all.data.otl.rm$site_code))
annual.data.otl.rm$site.id = as.numeric(as.factor(annual.data.otl.rm$site_code))
perennial.data.otl.rm$site.id = as.numeric(as.factor(perennial.data.otl.rm$site_code))
grass.data.otl.rm$site.id = as.numeric(as.factor(grass.data.otl.rm$site_code))
forb.data.otl.rm$site.id = as.numeric(as.factor(forb.data.otl.rm$site_code))

#### merge with environmental data ####
# mean annual precipitation data (MAP)

Site.info = read.csv("./Raw.Data/Site_Elev-Disturb.csv")
site.info.map = Site.info[,c(2,13)]

# drought severity index

dsi = read.csv("./Raw.Data/site.drt.dev.index.csv", row.names = 1)

# merge MAP and DSI

enviro.data = left_join(site.info.map,dsi, by = "site_code")

# override original data names so can use same code
all.data = left_join(all.data.otl.rm, enviro.data, by="site_code") # 83 unique sites
annual.data = left_join(annual.data.otl.rm, enviro.data, by="site_code")
perennial.data = left_join(perennial.data.otl.rm, enviro.data, by="site_code")
grass = left_join(grass.data.otl.rm, enviro.data, by="site_code")
forb = left_join(forb.data.otl.rm, enviro.data, by="site_code")

#### testing for correlation among traits ####

# saved as pdf 7 x 7

# all.data
all.data.2 = all.data
colnames(all.data.2)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(all.data.2)[22] = "Precipitation"
colnames(all.data.2)[23] = "DSI"

cor.mat.all = cor(all.data.2[,c(10:17,22,23)],use = "pairwise") 
corrplot(cor.mat.all, method="number",tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")

annual.data.2 = annual.data
colnames(annual.data.2)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(annual.data.2)[22] = "Precipitation"
colnames(annual.data.2)[23] = "DSI"

cor.mat.annual = cor(annual.data.2[,c(10:17,22,23)],use = "pairwise") 
corrplot(cor.mat.annual, method="number", tl.col = "black",bg = "gray70", is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")

perennial.data.2 = perennial.data
colnames(perennial.data.2)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(perennial.data.2)[22] = "Precipitation"
colnames(perennial.data.2)[23] = "DSI"

cor.mat.perennial = cor(perennial.data.2[,c(10:17,22,23)],use = "pairwise") 
corrplot(cor.mat.perennial, method="number", tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")

grass.2 = grass
colnames(grass.2)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(grass.2)[22] = "Precipitation"
colnames(grass.2)[23] = "DSI"

cor.mat.grass = cor(grass.2[,c(10:17,22,23)],use = "pairwise") 
corrplot(cor.mat.grass, method="number", tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")

forb.2 = forb
colnames(forb.2)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(forb.2)[22] = "Precipitation"
colnames(forb.2)[23] = "DSI"

cor.mat.forb = cor(forb.2[,c(10:17,22,23)],use = "pairwise") 
corrplot(cor.mat.forb, method="number", tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")

#### determining best parameter combination to generate 1000 trees ####
# bag fraction of 0.50 and 0.75, step.size of 25 and 50, tc = 10, 

all.data.map.dsi=gbm.step(data=all.data, gbm.x = c(10:17,22,23), gbm.y=9,
                          family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                          bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(all.data.map.dsi)

annual.data.map.dsi=gbm.step(data=annual.data, gbm.x = c(10:17,22,23), gbm.y=9,
                             family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                             bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(annual.data.map.dsi)

perennial.data.map.dsi=gbm.step(data=perennial.data, gbm.x = c(10:17,22,23), gbm.y=9,
                                family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                                bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(perennial.data.map.dsi)

grass.data.map.dsi=gbm.step(data=grass, gbm.x = c(10:17,22,23), gbm.y=9,
                            family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                            bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(grass.data.map.dsi)

forb.data.map.dsi=gbm.step(data=forb, gbm.x = c(10:17,22,23), gbm.y=9,
                           family = "gaussian", tree.complexity = 10, learning.rate = 0.00001,
                           bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(forb.data.map.dsi)


#### determining best tree complexity to use ####

# code edited to run for each of the data sets
# when loops stopped prematurely, code was run manually to generate 10 reps

# all variables
R2Obs.all.variables <- list()
importancePred.all.variables <- list()
nreps <- 10 #number of simulations
for (tcomp in 1:10) {
  R2Obs.all.variables[[tcomp]] <- numeric(nreps)
  importancePred.all.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:10),
                                                             ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    
    BRT.all.variables <- gbm.step(data=perennial.data,
                                  gbm.x = c(10:17,22,23),
                                  gbm.y = 9,
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
# all.data tc = 5
# annual tc = 1
# perennial tree tc = 10
# grass tc = 6
# forb tc = 7

#### all data ####

all.data.map.dsi=gbm.step(data=all.data, gbm.x = c(10:17,22,23), gbm.y=9,
                          family = "gaussian", tree.complexity = 5, learning.rate = 0.0005,
                          bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(all.data.map.dsi)
# 1800 trees Per.Expl = 6.93%

ggInfluence(all.data.map.dsi)

all.data.map.dsi$contributions$var = c("LeafN","Height","SLA","RTD",
                                       "Precipitation","RootN","Rooting Depth",
                                       "Root Diameter","DSI","SRL")
all.data.dsi.influence.plot = ggInfluence_test(all.data.map.dsi, main = expression("All Species (n = 1204), R"^2*" = 6.93%"), 
            col.bar = c("#769370","#769370","#769370","gray70",
                        "gray70","gray70","gray70","gray70","gray70","gray70"), 
            col.signif = "#B50200")

# saved as pdf from device (4 x 5)

# get data to plot partial dependency plots

all.data.map.dsi$contributions$var = c("leafN.mg.g","height.m","SLA_m2.kg","RTD.g.cm3","precip",
                                       "rootN.mg.g","root.depth_m",
                                       "rootDiam.mm", "mean.drt.sev.index","SRL_m.g")

all.data.dsi.prerun<- plot.gbm.4list(all.data.map.dsi)

all.data.dsi.boot <- gbm.bootstrap.functions(all.data.map.dsi, list.predictors=all.data.dsi.prerun, n.reps=1000)

all.data.map.dsi$contributions$var = c("LeafN","Height","SLA","RTD","Precipitation",
                                       "RootN","Rooting Depth",
                                       "Root Diameter","DSI","SRL")

all.data.map.dsi$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation","DSI")
colnames(all.data.map.dsi$gbm.call$dataframe)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(all.data.map.dsi$gbm.call$dataframe)[22:23] = c("Precipitation","DSI")

# significant for both: leafN, height, SLA

ggPD_boot_test(all.data.map.dsi,predictor = "LeafN",list.4.preds=all.data.dsi.prerun, col.line="#769370",
          booted.preds=all.data.dsi.boot$function.preds, cex.line=2, col.ci="#769370",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
          y.label = "Percent Cover Change")
ggPD_boot_test(all.data.map.dsi,predictor = "Height",list.4.preds=all.data.dsi.prerun, col.line="#769370",
               booted.preds=all.data.dsi.boot$function.preds, cex.line=2, col.ci="#769370",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")
ggPD_boot_test(all.data.map.dsi,predictor = "SLA",list.4.preds=all.data.dsi.prerun, col.line="#769370",
               booted.preds=all.data.dsi.boot$function.preds, cex.line=2, col.ci="#769370",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = FALSE,
               y.label = "Percent Cover Change")

# output individual plots as 3x3

ggPD(all.data.map.dsi, col.line="#769370", common.scale = FALSE, y.label = "Percent Cover Change")
# output all plots as one 8x8

# investigation of interactions
all.data.map.dsi$contributions$var = c("leafN.mg.g","height.m","SLA_m2.kg",
                                       "RTD.g.cm3","precip","rootN.mg.g",
                                       "root.depth_m","rootDiam.mm","mean.drt.sev.index","SRL_m.g")

all.data.map.dsi$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm","precip","mean.drt.sev.index")
colnames(all.data.map.dsi$gbm.call$dataframe)[10:17] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(all.data.map.dsi$gbm.call$dataframe)[22:23] = c("precip","mean.drt.sev.index")

gbm.interactions(all.data.map.dsi)$interactions
ggInteract_list(all.data.map.dsi, index = T)
# 6 RTD x 2 height 5.62
# 6 RTD x 1 leafN 3.17
# 9 precip x 2 height 2.19
# 9 precip x 1 leafN 1.66
# 10 DSI x 2 height 1.53

# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model
Interact_boot_dsi<-ggInteract_boot(c('RTD.g.cm3','height.m'),c('RTD.g.cm3','leafN.mg.g'),
                                        c('precip','height.m'),c('precip','leafN.mg.g'),
                                   c('mean.drt.sev.index','height.m'),
                                        nboots = 500, data=all.data, predictors =  c(10:17,22,23), 
                                        response="cover.change",
                                        family = "gaussian", tc = 5, lr = 0.0005, bf= 0.50, global.env=F)

# Significance histogram p-value< 0.05
ggInteract_boot_hist(data = Interact_boot_dsi, column = 2,obs = 5.62) # significant
# lower RTD x higher height
ggInteract_boot_hist(data = Interact_boot_dsi, column = 3,obs = 3.17) # significant
# lower RTD x lower leafN 
ggInteract_boot_hist(data = Interact_boot_dsi, column = 4,obs = 2.19) # significant
# higher precip x higher height
ggInteract_boot_hist(data = Interact_boot_dsi, column = 5,obs = 1.66) # significant
# higher DSI x higher height
ggInteract_boot_hist(data = Interact_boot_dsi, column = 6,obs = 1.53) # significant
# higher precip x higher leafN

ggInteract_2D_test(gbm.object = all.data.map.dsi, x="RTD.g.cm3",y="height.m",col.gradient = c("white","#769370"),
                   show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = F,
                   col.contour = "#254376",show.axis = T,legend = T, x.label = "RTD", y.label = "Height",
                   z.range = c(-0.66, 0.84), z.label = "% Cover Change", smooth = "average")
ggInteract_2D_test(gbm.object = all.data.map.dsi, x="RTD.g.cm3",y="leafN.mg.g",col.gradient = c("white","#769370"),
                   show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = F,
                   col.contour = "#254376",show.axis = T,legend = T, x.label = "RTD", y.label = "LeafN",
                   z.range = c(-0.99, 0.13), z.label = "% Cover Change", smooth = "average")
ggInteract_2D_test(gbm.object = all.data.map.dsi, x="precip",y="height.m",col.gradient = c("white","#769370"),
                   show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = F,
                   col.contour = "#254376",show.axis = T,legend = T, x.label = "Precipitation", y.label = "Height",
                   z.range = c(-0.46, 0.80), z.label = "% Cover Change", smooth = "average")
ggInteract_2D_test(gbm.object = all.data.map.dsi, x="precip",y="leafN.mg.g",col.gradient = c("white","#769370"),
                   show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = F,
                   col.contour = "#254376",show.axis = T,legend = T, x.label = "Precipitation", y.label = "LeafN",
                   z.range = c(-0.84, 0.11), z.label = "% Cover Change", smooth = "average")
ggInteract_2D_test(gbm.object = all.data.map.dsi, x="mean.drt.sev.index",y="height.m",col.gradient = c("white","#769370"),
                   show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = F,
                   col.contour = "#254376",show.axis = T,legend = T, x.label = "DSI", y.label = "Height",
                   z.range = c(-0.25, 1.05), z.label = "% Cover Change", smooth = "average")

#### annual data ####
annual.no.site.map.dsi=gbm.step(data=annual.data, gbm.x = c(10:17,22,23), gbm.y=9,
                                family = "gaussian", tree.complexity = 1, learning.rate = 0.0005,
                                bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(annual.no.site.map.dsi)
# 1350 trees Per.Expl = 4.08%

ggInfluence(annual.no.site.map.dsi)

annual.no.site.map.dsi$contributions$var = c( "RootN","RTD","DSI","Precipitation",
                                              "LeafN","SRL","Rooting Depth",
                                              "Height", "SLA","Root Diameter")
annual.influence.plot = ggInfluence_test(annual.no.site.map.dsi, main = expression("Annuals (n = 284), R"^2*" = 4.08%"), 
                                    col.bar = c("#979461","#979461","#979461","gray70",
                                                "gray70","gray70","gray70","gray70","gray70","gray70"),
                                    col.signif = "#B50200")
# saved as pdf from device (4 x 5)

# get data to plot partial dependency plots
annual.no.site.map.dsi$contributions$var = c("rootN.mg.g","RTD.g.cm3","mean.drt.sev.index",
                                             "precip","leafN.mg.g","SRL_m.g","root.depth_m",
                                             "height.m","SLA_m2.kg","rootDiam.mm")

annual.dsi.prerun<- plot.gbm.4list(annual.no.site.map.dsi)

annual.dsi.boot <- gbm.bootstrap.functions(annual.no.site.map.dsi, list.predictors=annual.dsi.prerun, n.reps=1000)

annual.no.site.map.dsi$contributions$var = c( "RootN","RTD","DSI","Precipitation",
                                              "LeafN","SRL","Rooting Depth",
                                              "Height", "SLA","Root Diameter")

annual.no.site.map.dsi$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation","DSI")
colnames(annual.no.site.map.dsi$gbm.call$dataframe)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(annual.no.site.map.dsi$gbm.call$dataframe)[22:23] = c("Precipitation","DSI")

# significant for DSI: RootN, RTD, DSI

ggPD_boot_test(annual.no.site.map.dsi,predictor = "RootN",list.4.preds=annual.dsi.prerun, col.line="#979461",
          booted.preds=annual.dsi.boot$function.preds, cex.line=2, col.ci="#979461",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = FALSE,
          y.label = "Percent Cover Change")
ggPD_boot_test(annual.no.site.map.dsi,predictor = "RTD",list.4.preds=annual.dsi.prerun, col.line="#979461",
               booted.preds=annual.dsi.boot$function.preds, cex.line=2, col.ci="#979461",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = FALSE,
               y.label = "Percent Cover Change")
# had to change x-axis in ggPD_boot_edit.R for this graph
ggPD_boot_test(annual.no.site.map.dsi,predictor = "DSI",list.4.preds=annual.dsi.prerun, col.line="#979461",
               booted.preds=annual.dsi.boot$function.preds, cex.line=2, col.ci="#979461",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = TRUE,
               y.label = "Percent Cover Change")

ggPD(annual.no.site.map.dsi, col.line="#979461", common.scale = FALSE, y.label = "Percent Cover Change")

# output all plots as one 6x6

#### Perennial without woody ####
perennial.data.map.dsi=gbm.step(data=perennial.data, gbm.x = c(10:17,22,23), gbm.y=9,
                                family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                                bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)


ggPerformance(perennial.data.map.dsi)
# 1000 trees Per.Expl = 1.82%

gInfluence(perennial.data.map.dsi)

perennial.data.map.dsi$contributions$var = c("LeafN","Height","SLA",
                                         "SRL","RootN","Precipitation",
                                         "Rooting Depth","Root Diameter","RTD","DSI")
perennial.influence.plot = ggInfluence_test(perennial.data.map.dsi, main = expression("Perennials (n = 869), R"^2*" = 1.82%"), 
                                       col.bar = c("#F1C646","#F1C646","#F1C646","gray70",
                                                   "gray70","gray70","gray70","gray70","gray70","gray70"
                                       ), col.signif = "#B50200")

# get data to plot partial dependency plots

perennial.data.map.dsi$contributions$var = c("leafN.mg.g","height.m","SLA_m2.kg",
                                             "SRL_m.g","rootN.mg.g","precip",
                                             "root.depth_m","rootDiam.mm","RTD.g.cm3","mean.drt.sev.index")

perennial.dsi.prerun<- plot.gbm.4list(perennial.data.map.dsi)

perennial.dsi.boot <- gbm.bootstrap.functions(perennial.data.map.dsi, list.predictors=perennial.dsi.prerun, n.reps=1000)

perennial.data.map.dsi$contributions$var = c("LeafN","Height","SLA",
                                             "SRL","RootN","Precipitation",
                                             "Rooting Depth","Root Diameter","RTD","DSI")

erennial.data.map.dsi$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation","DSI")
colnames(perennial.data.map.dsi$gbm.call$dataframe)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(perennial.data.map.dsi$gbm.call$dataframe)[22:23] = c("Precipitation","DSI")

# significant for both: leafN, SLA, height

ggPD_boot_test(perennial.data.map.dsi,predictor = "LeafN",list.4.preds=perennial.dsi.prerun, col.line="#F1C646",
          booted.preds=perennial.dsi.boot$function.preds, cex.line=2, col.ci="#F1C646",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = FALSE,
          y.label = "Percent Cover Change")
ggPD_boot_test(perennial.data.map.dsi,predictor = "SLA",list.4.preds=perennial.dsi.prerun, col.line="#F1C646",
               booted.preds=perennial.dsi.boot$function.preds, cex.line=2, col.ci="#F1C646",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = TRUE,
               y.label = "Percent Cover Change")
ggPD_boot_test(perennial.data.map.dsi,predictor = "Height",list.4.preds=perennial.dsi.prerun, col.line="#F1C646",
               booted.preds=perennial.dsi.boot$function.preds, cex.line=2, col.ci="#F1C646",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = TRUE,
               y.label = "Percent Cover Change")

ggPD(perennial.data.map.dsi, col.line="#F1C646", common.scale = FALSE, y.label = "Percent Cover Change")

# output all plots as one 6x6

# investigation of interactions

perennial.data.map.dsi$contributions$var = c("leafN.mg.g","height.m","SLA_m2.kg",
                                             "SRL_m.g","rootN.mg.g","precip",
                                             "root.depth_m","rootDiam.mm","RTD.g.cm3","mean.drt.sev.index")

perennial.data.map.dsi$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm", "precip","mean.drt.sev.index")
colnames(perennial.data.map.dsi$gbm.call$dataframe)[10:17] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(perennial.data.map.dsi$gbm.call$dataframe)[22:23] = c("precip","mean.drt.sev.index")

gbm.interactions(perennial.data.map.dsi)$rank.list
ggInteract_list(perennial.data.map.dsi)
# 7 SRL x 2 height 0.24
# 7 SRL x 3 rootN 0.09
# 6 RTD x 1 leafN 0.08
# 9 precip x 1 leafN 0.06
# 6 SRL x 1 leafN 0.04

# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model

Interact_boot_perennial_dsi<-ggInteract_boot(c('SRL_m.g','height.m'),c('SRL_m.g','rootN.mg.g'),
                                             c('RTD.g.cm3','leafN.mg.g'),c('precip','leafN.mg.g'),
                                         c('SRL_m.g','leafN.mg.g'),
                                         nboots = 500,data=perennial.data, predictors =  c(10:17,22,23), 
                                         response="cover.change",
                                         family = "gaussian", tc = 10, lr = 0.0001, bf= 0.75, global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_perennial_dsi, column = 2,obs = 0.24) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial_dsi, column = 3,obs = 0.09) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial_dsi, column = 4,obs = 0.08) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial_dsi, column = 5,obs = 0.06) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial_dsi, column = 6,obs = 0.04) # non-significant

#### Grass ####
grass.map.dsi=gbm.step(data=grass, gbm.x = c(10:17,22,23), gbm.y=9,
                       family = "gaussian", tree.complexity = 6, learning.rate = 0.0001,
                       bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(grass.map.dsi)
# 1850 trees Per.Expl = 3.51%

ggInfluence(grass.map)

grass.map.dsi$contributions$var = c("LeafN", "Precipitation","Height",
                                    "SLA","DSI","Root Diameter","RootN",
                                    "RTD","SRL","Rooting Depth")
                                
grass.influence.plot = ggInfluence_test(grass.map.dsi, main = expression("Grasses (n = 373), R"^2*" = 3.51%"), 
                                   col.bar = c("#6E687E","#6E687E","#6E687E","#6E687E",
                                               "gray70","gray70","gray70","gray70","gray70","gray70"
                                   ), col.signif = "#B50200")

# get data to plot partial dependency plots

grass.map.dsi$contributions$var = c("leafN.mg.g", "precip","height.m",
                                    "SLA_m2.kg","mean.drt.sev.index","rootDiam.mm","rootN.mg.g",
                                    "RTD.g.cm3","SRL_m.g","root.depth_m")

grass.dsi.prerun<- plot.gbm.4list(grass.map.dsi)

grass.dsi.boot <- gbm.bootstrap.functions(grass.map.dsi, list.predictors=grass.dsi.prerun, n.reps=1000)

grass.map.dsi$contributions$var = c("LeafN", "Precipitation","Height",
                                    "SLA","DSI","Root Diameter","RootN",
                                    "RTD","SRL","Rooting Depth")

grass.map.dsi$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation","DSI")
colnames(grass.map.dsi$gbm.call$dataframe)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(grass.map.dsi$gbm.call$dataframe)[22:23] = c("Precipitation","DSI")

# significant for both: leafN, precip, height, SLA

ggPD_boot_test(grass.map.dsi,predictor = "LeafN",list.4.preds=grass.dsi.prerun, col.line="#6E687E",
               booted.preds=grass.dsi.boot$function.preds, cex.line=2, col.ci="#6E687E",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")
ggPD_boot_test(grass.map.dsi,predictor = "Precipitation",list.4.preds=grass.dsi.prerun, col.line="#6E687E",
               booted.preds=grass.dsi.boot$function.preds, cex.line=2, col.ci="#6E687E",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")
ggPD_boot_test(grass.map.dsi,predictor = "Height",list.4.preds=grass.dsi.prerun, col.line="#6E687E",
               booted.preds=grass.dsi.boot$function.preds, cex.line=2, col.ci="#6E687E",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")
ggPD_boot_test(grass.map.dsi,predictor = "SLA",list.4.preds=grass.dsi.prerun, col.line="#6E687E",
               booted.preds=grass.dsi.boot$function.preds, cex.line=2, col.ci="#6E687E",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")

ggPD(grass.map.dsi, col.line="#6E687E", common.scale = FALSE, y.label = "Percent Cover Change")

# output all plots as one 6x6

# investigation of interactions

grass.map.dsi$contributions$var = c("leafN.mg.g", "precip","height.m",
                                    "SLA_m2.kg","mean.drt.sev.index","rootDiam.mm","rootN.mg.g",
                                    "RTD.g.cm3","SRL_m.g","root.depth_m")

grass.map.dsi$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm","precip","mean.drt.sev.index")
colnames(grass.map.dsi$gbm.call$dataframe)[10:17] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(grass.map.dsi$gbm.call$dataframe)[22:23] = c("precip","mean.drt.sev.index")

gbm.interactions(grass.map.dsi)$interactions
ggInteract_list(grass.map.dsi, index = T)
# 9 precip x 2 height 0.29
# 10 DSI x 2 height 0.23
# 9 precip x 8 rootDiam 0.12
# 8 rootDiam x 1 leafN 0.11
# 8 rootDiam x 3 rootN 0.08

# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model
Interact_boot_dsi.grass<-ggInteract_boot(c('precip','height.m'),c('mean.drt.sev.index','height.m'),
                                          c('precip','rootDiam.mm  '),c('rootDiam.mm ','leafN.mg.g'),
                                         c('rootDiam.mm ','rootN.mg.g'),
                                          nboots = 500, data=grass, predictors =  c(10:17,22,23), 
                                          response="cover.change",
                                          family = "gaussian", tc = 6, lr = 0.0001, bf= 0.50, global.env=F)
# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_dsi.grass, column = 2,obs = 0.29) # non-significant
ggInteract_boot_hist(data = Interact_boot_dsi.grass, column = 3,obs = 0.23) # non-significant
ggInteract_boot_hist(data = Interact_boot_dsi.grass, column = 4,obs = 0.12) # non-significant
ggInteract_boot_hist(data = Interact_boot_dsi.grass, column = 5,obs = 0.11) # non-significant
ggInteract_boot_hist(data = Interact_boot_dsi.grass, column = 6,obs = 0.08) # non-significant

#### Forb ####
forb.map.dsi=gbm.step(data=forb, gbm.x = c(10:17,22,23), gbm.y=9,
                      family = "gaussian", tree.complexity = 7, learning.rate = 0.0001,
                      bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(forb.map.dsi)
# 1000 trees Per.Expl = 1.48%

forb.map.dsi$contributions$var = c("LeafN","SLA","Rooting Depth",
                                   "DSI","Height","SRL",
                                   "Precipitation","Root Diameter",
                                   "RTD", "RootN")

forb.influence.plot = ggInfluence_test(forb.map.dsi, main = expression("Forbs (n = 666), R"^2*" = 1.48%"),
                                  col.bar = c("#F17236","#F17236","#F17236","#F17236",
                                              "gray70","gray70","gray70","gray70","gray70","gray70"
                                  ), col.signif = "#B50200")

# get data to plot partial dependency plots

forb.map.dsi$contributions$var = c("leafN.mg.g","SLA_m2.kg","root.depth_m",
                                   "mean.drt.sev.index","height.m","SRL_m.g",
                                   "precip","rootDiam.mm",
                                   "RTD.g.cm3", "rootN.mg.g")

forb.dsi.prerun<- plot.gbm.4list(forb.map.dsi)

forb.dsi.boot <- gbm.bootstrap.functions(forb.map.dsi, list.predictors=forb.dsi.prerun, n.reps=1000)

forb.map.dsi$contributions$var = c("LeafN","SLA","Rooting Depth",
                                   "DSI","Height","SRL",
                                   "Precipitation","Root Diameter",
                                   "RTD", "RootN")

forb.map.dsi$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation","DSI")
colnames(forb.map.dsi$gbm.call$dataframe)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(forb.map.dsi$gbm.call$dataframe)[22:23] = c("Precipitation","DSI")

# significant for DSI: leafN, SLA, depth, DSI

ggPD_boot_test(forb.map.dsi,predictor = "LeafN",list.4.preds=forb.dsi.prerun, col.line="#F17236",
               booted.preds=forb.dsi.boot$function.preds, cex.line=2, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")
ggPD_boot_test(forb.map.dsi,predictor = "SLA",list.4.preds=forb.dsi.prerun, col.line="#F17236",
               booted.preds=forb.dsi.boot$function.preds, cex.line=2, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")
ggPD_boot_test(forb.map.dsi,predictor = "Rooting Depth",list.4.preds=forb.dsi.prerun, col.line="#F17236",
               booted.preds=forb.dsi.boot$function.preds, cex.line=2, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")
ggPD_boot_test(forb.map.dsi,predictor = "DSI",list.4.preds=forb.dsi.prerun, col.line="#F17236",
               booted.preds=forb.dsi.boot$function.preds, cex.line=2, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")

ggPD(forb.map.dsi, col.line="#F17236", common.scale = FALSE, y.label = "Percent Cover Change")

# output all plots as one 6x6

# investigation of interactions

forb.map.dsi$contributions$var = c("leafN.mg.g","SLA_m2.kg","root.depth_m",
                                   "mean.drt.sev.index","height.m","SRL_m.g",
                                   "precip","rootDiam.mm",
                                   "RTD.g.cm3", "rootN.mg.g")

forb.map.dsi$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm","precip","mean.drt.sev.index")
colnames(forb.map.dsi$gbm.call$dataframe)[10:17] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(forb.map.dsi$gbm.call$dataframe)[22:23] = c("precip","mean.drt.sev.index")

gbm.interactions(forb.map.dsi)$interactions
ggInteract_list(forb.map.dsi, index = T)
# 4 SLA x 1 leafN 0.06
# 5 depth x 4 SLA 0.04
# 10 DSI x 5 depth 0.01
# 10 DSI x 4 SLA 0.01
# 10 DSI x 2 height 0.01

# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model

Interact_boot_forb_dsi<-ggInteract_boot(c('SLA_m2.kg','leafN.mg.g'),c('root.depth_m ','SLA_m2.kg'),
                                    c('mean.drt.sev.index','root.depth_m'),c('mean.drt.sev.index','SLA_m2.kg'),
                                    c('mean.drt.sev.index','height.m'),
                                    nboots = 500, data=forb, predictors =  c(10:17,22,23), 
                                    response="cover.change",
                                    family = "gaussian", tc = 7, lr = 0.0001, bf= 0.50, global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_forb_dsi, column = 2,obs = 0.06) # non-significant
ggInteract_boot_hist(data = Interact_boot_forb_dsi, column = 3,obs = 0.04) # non-significant
ggInteract_boot_hist(data = Interact_boot_forb_dsi, column = 4,obs = 0.01) # non-significant
ggInteract_boot_hist(data = Interact_boot_forb_dsi, column = 5,obs = 0.01) # non-significant
ggInteract_boot_hist(data = Interact_boot_forb_dsi, column = 6,obs = 0.01) # non-significant

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
mean(all.data$mean.drt.sev.index, na.rm = TRUE)
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
sd(all.data$mean.drt.sev.index, na.rm = TRUE)
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
range(all.data$mean.drt.sev.index, na.rm = TRUE)
range(all.data$cover.change)

sum(is.na(all.data$leafN.mg.g))/1204*100
sum(is.na(all.data$height.m))/1204*100
sum(is.na(all.data$rootN.mg.g))/1204*100
sum(is.na(all.data$SLA_m2.kg))/1204*100
sum(is.na(all.data$root.depth_m))/1204*100
sum(is.na(all.data$RTD.g.cm3))/1204*100
sum(is.na(all.data$SRL_m.g))/1204*100
sum(is.na(all.data$rootDiam.mm))/1204*100
sum(is.na(all.data$precip))/1204*100
sum(is.na(all.data$mean.drt.sev.index))/1204*100
sum(is.na(all.data$cover.change))/1204*100

mean(annual.data$leafN.mg.g, na.rm = TRUE)
mean(annual.data$height.m, na.rm = TRUE)
mean(annual.data$rootN.mg.g, na.rm = TRUE)
mean(annual.data$SLA_m2.kg, na.rm = TRUE)
mean(annual.data$root.depth_m, na.rm = TRUE)
mean(annual.data$RTD.g.cm3, na.rm = TRUE)
mean(annual.data$SRL_m.g, na.rm = TRUE)
mean(annual.data$rootDiam.mm, na.rm = TRUE)
mean(annual.data$precip, na.rm = TRUE)
mean(annual.data$mean.drt.sev.index, na.rm = TRUE)
mean(annual.data$cover.change, na.rm = TRUE)

sd(annual.data$leafN.mg.g, na.rm = TRUE)
sd(annual.data$height.m, na.rm = TRUE)
sd(annual.data$rootN.mg.g, na.rm = TRUE)
sd(annual.data$SLA_m2.kg, na.rm = TRUE)
sd(annual.data$root.depth_m, na.rm = TRUE)
sd(annual.data$RTD.g.cm3, na.rm = TRUE)
sd(annual.data$SRL_m.g, na.rm = TRUE)
sd(annual.data$rootDiam.mm, na.rm = TRUE)
sd(annual.data$precip, na.rm = TRUE)
sd(annual.data$mean.drt.sev.index, na.rm = TRUE)
sd(annual.data$cover.change, na.rm = TRUE)

range(annual.data$leafN.mg.g, na.rm = TRUE)
range(annual.data$height.m, na.rm = TRUE)
range(annual.data$rootN.mg.g, na.rm = TRUE)
range(annual.data$SLA_m2.kg, na.rm = TRUE)
range(annual.data$root.depth_m, na.rm = TRUE)
range(annual.data$RTD.g.cm3, na.rm = TRUE)
range(annual.data$SRL_m.g, na.rm = TRUE)
range(annual.data$rootDiam.mm, na.rm = TRUE)
range(annual.data$precip, na.rm = TRUE)
range(annual.data$mean.drt.sev.index, na.rm = TRUE)
range(annual.data$cover.change, na.rm = TRUE)

mean(perennial.data$leafN.mg.g, na.rm = TRUE)
mean(perennial.data$height.m, na.rm = TRUE)
mean(perennial.data$rootN.mg.g, na.rm = TRUE)
mean(perennial.data$SLA_m2.kg, na.rm = TRUE)
mean(perennial.data$root.depth_m, na.rm = TRUE)
mean(perennial.data$RTD.g.cm3, na.rm = TRUE)
mean(perennial.data$SRL_m.g, na.rm = TRUE)
mean(perennial.data$rootDiam.mm, na.rm = TRUE)
mean(perennial.data$precip, na.rm = TRUE)
mean(perennial.data$mean.drt.sev.index, na.rm = TRUE)
mean(perennial.data$cover.change, na.rm = TRUE)

sd(perennial.data$leafN.mg.g, na.rm = TRUE)
sd(perennial.data$height.m, na.rm = TRUE)
sd(perennial.data$rootN.mg.g, na.rm = TRUE)
sd(perennial.data$SLA_m2.kg, na.rm = TRUE)
sd(perennial.data$root.depth_m, na.rm = TRUE)
sd(perennial.data$RTD.g.cm3, na.rm = TRUE)
sd(perennial.data$SRL_m.g, na.rm = TRUE)
sd(perennial.data$rootDiam.mm, na.rm = TRUE)
sd(perennial.data$precip, na.rm = TRUE)
sd(perennial.data$mean.drt.sev.index, na.rm = TRUE)
sd(perennial.data$cover.change, na.rm = TRUE)

range(perennial.data$leafN.mg.g, na.rm = TRUE)
range(perennial.data$height.m, na.rm = TRUE)
range(perennial.data$rootN.mg.g, na.rm = TRUE)
range(perennial.data$SLA_m2.kg, na.rm = TRUE)
range(perennial.data$root.depth_m, na.rm = TRUE)
range(perennial.data$RTD.g.cm3, na.rm = TRUE)
range(perennial.data$SRL_m.g, na.rm = TRUE)
range(perennial.data$rootDiam.mm, na.rm = TRUE)
range(perennial.data$precip, na.rm = TRUE)
range(perennial.data$mean.drt.sev.index, na.rm = TRUE)
range(perennial.data$cover.change, na.rm = TRUE)

mean(grass$leafN.mg.g, na.rm = TRUE)
mean(grass$height.m, na.rm = TRUE)
mean(grass$rootN.mg.g, na.rm = TRUE)
mean(grass$SLA_m2.kg, na.rm = TRUE)
mean(grass$root.depth_m, na.rm = TRUE)
mean(grass$RTD.g.cm3, na.rm = TRUE)
mean(grass$SRL_m.g, na.rm = TRUE)
mean(grass$rootDiam.mm, na.rm = TRUE)
mean(grass$precip, na.rm = TRUE)
mean(grass$mean.drt.sev.index, na.rm = TRUE)
mean(grass$cover.change, na.rm = TRUE)

sd(grass$leafN.mg.g, na.rm = TRUE)
sd(grass$height.m, na.rm = TRUE)
sd(grass$rootN.mg.g, na.rm = TRUE)
sd(grass$SLA_m2.kg, na.rm = TRUE)
sd(grass$root.depth_m, na.rm = TRUE)
sd(grass$RTD.g.cm3, na.rm = TRUE)
sd(grass$SRL_m.g, na.rm = TRUE)
sd(grass$rootDiam.mm, na.rm = TRUE)
sd(grass$precip, na.rm = TRUE)
sd(grass$mean.drt.sev.index, na.rm = TRUE)
sd(grass$cover.change, na.rm = TRUE)

range(grass$leafN.mg.g, na.rm = TRUE)
range(grass$height.m, na.rm = TRUE)
range(grass$rootN.mg.g, na.rm = TRUE)
range(grass$SLA_m2.kg, na.rm = TRUE)
range(grass$root.depth_m, na.rm = TRUE)
range(grass$RTD.g.cm3, na.rm = TRUE)
range(grass$SRL_m.g, na.rm = TRUE)
range(grass$rootDiam.mm, na.rm = TRUE)
range(grass$precip, na.rm = TRUE)
range(grass$mean.drt.sev.index, na.rm = TRUE)
range(grass$cover.change, na.rm = TRUE)

mean(forb$leafN.mg.g, na.rm = TRUE)
mean(forb$height.m, na.rm = TRUE)
mean(forb$rootN.mg.g, na.rm = TRUE)
mean(forb$SLA_m2.kg, na.rm = TRUE)
mean(forb$root.depth_m, na.rm = TRUE)
mean(forb$RTD.g.cm3, na.rm = TRUE)
mean(forb$SRL_m.g, na.rm = TRUE)
mean(forb$rootDiam.mm, na.rm = TRUE)
mean(forb$precip, na.rm = TRUE)
mean(forb$mean.drt.sev.index, na.rm = TRUE)
mean(forb$cover.change, na.rm = TRUE)

sd(forb$leafN.mg.g, na.rm = TRUE)
sd(forb$height.m, na.rm = TRUE)
sd(forb$rootN.mg.g, na.rm = TRUE)
sd(forb$SLA_m2.kg, na.rm = TRUE)
sd(forb$root.depth_m, na.rm = TRUE)
sd(forb$RTD.g.cm3, na.rm = TRUE)
sd(forb$SRL_m.g, na.rm = TRUE)
sd(forb$rootDiam.mm, na.rm = TRUE)
sd(forb$precip, na.rm = TRUE)
sd(forb$mean.drt.sev.index, na.rm = TRUE)
sd(forb$cover.change, na.rm = TRUE)

range(forb$leafN.mg.g, na.rm = TRUE)
range(forb$height.m, na.rm = TRUE)
range(forb$rootN.mg.g, na.rm = TRUE)
range(forb$SLA_m2.kg, na.rm = TRUE)
range(forb$root.depth_m, na.rm = TRUE)
range(forb$RTD.g.cm3, na.rm = TRUE)
range(forb$SRL_m.g, na.rm = TRUE)
range(forb$rootDiam.mm, na.rm = TRUE)
range(forb$precip, na.rm = TRUE)
range(forb$mean.drt.sev.index, na.rm = TRUE)
range(forb$cover.change, na.rm = TRUE)

table(perennial.data$functional_group)
# 230 forbs, 180 grasses
table(annual.data$functional_group)
# 76 forbs, 39 grasses
