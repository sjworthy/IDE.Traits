# removing outlier traits

#### final.data.CC ####
# model of trait data with complete cases for all species

final.data.CC = read.csv("./Formatted.Data/Revisions/final.data.CC.csv",row.names = 1)

#### checking for trait outliers from trait data with woody species ####
# not going to remove outliers since need complete cases of traits for individuals to be included in analyses using lmer

hist(final.data.CC$leafN.mg.g)
boxplot(final.data.CC$leafN.mg.g)
mean = mean(final.data.CC$leafN.mg.g, na.rm = TRUE)
std = sd(final.data.CC$leafN.mg.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$leafN.mg.g[which(final.data.CC$leafN.mg.g <Tmin | final.data.CC$leafN.mg.g > Tmax)])
# removed leafN 54.71597
# percent removed 
table(is.na(final.data.CC$leafN.mg.g))
958-368
2/590*100 #0.34%

hist(final.data.CC$height.m)
boxplot(final.data.CC$height.m)
mean = mean(final.data.CC$height.m, na.rm = TRUE)
std = sd(final.data.CC$height.m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$height.m[which(final.data.CC$height.m <Tmin | final.data.CC$height.m > Tmax)])
# removed height 12.84847 27.10186 32.57366
# percent removed 
table(is.na(final.data.CC$height.m))
958-216
14/742*100 #1.89%

hist(final.data.CC$rootN.mg.g)
boxplot(final.data.CC$rootN.mg.g)
mean = mean(final.data.CC$rootN.mg.g, na.rm = TRUE)
std = sd(final.data.CC$rootN.mg.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$rootN.mg.g[which(final.data.CC$rootN.mg.g <Tmin | final.data.CC$rootN.mg.g > Tmax)])
# remove rootN 31.17241 31.34076 33.34667 33.34667 33.34667
# percent removed 
table(is.na(final.data.CC$rootN.mg.g))
958-613
6/345*100 #1.74%

hist(final.data.CC$SLA_m2.kg)
boxplot(final.data.CC$SLA_m2.kg)
mean = mean(final.data.CC$SLA_m2.kg, na.rm = TRUE)
std = sd(final.data.CC$SLA_m2.kg, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$SLA_m2.kg[which(final.data.CC$SLA_m2.kg <Tmin | final.data.CC$SLA_m2.kg > Tmax)])
# remove SLA: none
# percent removed 
table(is.na(final.data.CC$SLA_m2.kg))
958-244
8/714*100 #1.12%

hist(final.data.CC$root.depth_m)
boxplot(final.data.CC$root.depth_m)
mean = mean(final.data.CC$root.depth_m, na.rm = TRUE)
std = sd(final.data.CC$root.depth_m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$root.depth_m[which(final.data.CC$root.depth_m <Tmin | final.data.CC$root.depth_m > Tmax)])
# remove depth 2.426850 2.426850 2.598333 2.598333 2.682500 2.682500 2.915000
# percent removed 
table(is.na(final.data.CC$root.depth_m))
958-478
15/480*100 #3.13%

hist(final.data.CC$RTD.groot.cahill.merge)
boxplot(final.data.CC$RTD.groot.cahill.merge)
mean = mean(final.data.CC$RTD.groot.cahill.merge, na.rm = TRUE)
std = sd(final.data.CC$RTD.groot.cahill.merge, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$RTD.groot.cahill.merge[which(final.data.CC$RTD.groot.cahill.merge <Tmin | final.data.CC$RTD.groot.cahill.merge > Tmax)])
# remove RTD 1.19455 1.19455 1.19455 1.19455
# percent removed 
table(is.na(final.data.CC$RTD.groot.cahill.merge))
958-581
1/377*100 #0.27%

hist(final.data.CC$SRL.groot.cahill.merge)
boxplot(final.data.CC$SRL.groot.cahill.merge)
mean = mean(final.data.CC$SRL.groot.cahill.merge, na.rm = TRUE)
std = sd(final.data.CC$SRL.groot.cahill.merge, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$SRL.groot.cahill.merge[which(final.data.CC$SRL.groot.cahill.merge <Tmin | final.data.CC$SRL.groot.cahill.merge > Tmax)])
# remove SRL 471.2364 471.2364 471.2364 601.5858 627.5050
# percent removed 
table(is.na(final.data.CC$SRL.groot.cahill.merge))
958-575
10/383*100 #2.61%

hist(final.data.CC$rootDiam.mm)
boxplot(final.data.CC$rootDiam.mm)
mean = mean(final.data.CC$rootDiam.mm, na.rm = TRUE)
std = sd(final.data.CC$rootDiam.mm, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$rootDiam.mm[which(final.data.CC$rootDiam.mm <Tmin | final.data.CC$rootDiam.mm > Tmax)])
# remove diam 1.3965
# percent removed 
table(is.na(final.data.CC$rootDiam.mm))
958-556
4/402*100 #0.995%

hist(final.data.CC$RMF.g.g)
boxplot(final.data.CC$RMF.g.g)
mean = mean(final.data.CC$RMF.g.g, na.rm = TRUE)
std = sd(final.data.CC$RMF.g.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$RMF.g.g[which(final.data.CC$RMF.g.g <Tmin | final.data.CC$RMF.g.g > Tmax)])
# remove RMF: none
# percent removed 
table(is.na(final.data.CC$RMF.g.g))
958-556
4/402*100 #0.995%


# 25 total removed because some removed with multiple traits
final.data.CC.2 = final.data.CC %>%
  filter(round(leafN.mg.g, 5) != 54.71597) %>% # 2
  filter(round(height.m, 5) < 12.84847) %>% # 3
  filter(round(rootN.mg.g, 5) < 31.17241) %>% # 5
  filter(round(root.depth_m, 5) < 2.426850) %>% # 7
  filter(round(RTD.groot.cahill.merge, 5) < 1.19455) %>% # 4
  filter(round(SRL.groot.cahill.merge, 5) < 471.2364) %>% # 5
  filter(round(rootDiam.mm, 5) != 1.3965) # 1

# look at correlations between traits
cor.traits.cc = cor(final.data.CC.2[,c(4:12)],use = "pairwise") 
corrplot(cor.traits.cc, method="number",tl.col = "black", bg = "gray70",is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")
# highest correlation is -0.38 and 0.35

final.data.CC.3 = scale(final.data.CC.2[,c(4:12)])
final.data.CC.4 = cbind(final.data.CC.2[,c(1:3)],final.data.CC.3)

# 211 data points estimating 11 variables

CC.model = lmer(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                  SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + (1|site_code) + (1|Taxon), data = final.data.CC.4)
summary(CC.model)
# RTD is positive, significant
r.squaredGLMM(CC.model) # 0.06, 0.21

plot(CC.model)
output.plot = allEffects(CC.model)
plot(output.plot)
hist(resid(CC.model))

RTD.plot = plot_model(CC.model, type = "pred", terms = "RTD.groot.cahill.merge")
RTD.plot + geom_point(data = final.data.CC.3, aes(x = RTD.groot.cahill.merge, y = cover.change))

# Add interaction between RTD and SRL

CC.model.interact = lmer(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                           SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + SRL.groot.cahill.merge*RTD.groot.cahill.merge +
                           (1|site_code) + (1|Taxon), data = final.data.CC.4)
summary(CC.model.interact)
# RTD is significant
r.squaredGLMM(CC.model.interact) # 0.07, 0.21

plot(CC.model.interact)
output.plot = allEffects(CC.model.interact)
plot(output.plot)
hist(resid(CC.model.interact))

#### Bayes model ###

get_prior(cover.change ~ 1 + leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm + SRL.groot.cahill.merge + 
            RTD.groot.cahill.merge + RMF.g.g + (1|site_code) + (1|Taxon),
          family = gaussian(),
          data = final.data.CC.4)

priors <- c(prior(normal(0, 10), class = b))


CC.model = brm(cover.change ~ leafN.mg.g + height.m + rootN.mg.g + SLA_m2.kg + root.depth_m + rootDiam.mm +
                 SRL.groot.cahill.merge + RTD.groot.cahill.merge + RMF.g.g + (1|site_code) + (1|Taxon), 
               family = gaussian(),
               prior = priors,
               data = final.data.CC.4)

summary(CC.model)
plot(CC.model, nvariables=6, ask=FALSE)
pp_check(CC.model, ndraws = 500)

conditional_effects(CC.model, re_formula = NULL)

#### checking for trait outliers from trait data without woody species ####
# not going to remove outliers since need complete cases of traits for individuals to be included in analyses using lmer

final.data.CC = read.csv("./Formatted.Data/Revisions/final.data.CC.NW.csv")

hist(final.data.CC$leafN.mg.g)
boxplot(final.data.CC$leafN.mg.g)
mean = mean(final.data.CC$leafN.mg.g, na.rm = TRUE)
std = sd(final.data.CC$leafN.mg.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$leafN.mg.g[which(final.data.CC$leafN.mg.g <Tmin | final.data.CC$leafN.mg.g > Tmax)])
# removed leafN 54.71597
# percent removed 
table(is.na(final.data.CC$leafN.mg.g))
958-368
2/590*100 #0.34%

hist(final.data.CC$height.m)
boxplot(final.data.CC$height.m)
mean = mean(final.data.CC$height.m, na.rm = TRUE)
std = sd(final.data.CC$height.m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$height.m[which(final.data.CC$height.m <Tmin | final.data.CC$height.m > Tmax)])
# removed height 1.800048 2.847733
# percent removed 
table(is.na(final.data.CC$height.m))
958-216
14/742*100 #1.89%

hist(final.data.CC$rootN.mg.g)
boxplot(final.data.CC$rootN.mg.g)
mean = mean(final.data.CC$rootN.mg.g, na.rm = TRUE)
std = sd(final.data.CC$rootN.mg.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$rootN.mg.g[which(final.data.CC$rootN.mg.g <Tmin | final.data.CC$rootN.mg.g > Tmax)])
# remove rootN 31.17241 31.34076 33.34667 33.34667 33.34667
# percent removed 
table(is.na(final.data.CC$rootN.mg.g))
958-613
6/345*100 #1.74%

hist(final.data.CC$SLA_m2.kg)
boxplot(final.data.CC$SLA_m2.kg)
mean = mean(final.data.CC$SLA_m2.kg, na.rm = TRUE)
std = sd(final.data.CC$SLA_m2.kg, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$SLA_m2.kg[which(final.data.CC$SLA_m2.kg <Tmin | final.data.CC$SLA_m2.kg > Tmax)])
# remove SLA: none
# percent removed 
table(is.na(final.data.CC$SLA_m2.kg))
958-244
8/714*100 #1.12%

hist(final.data.CC$root.depth_m)
boxplot(final.data.CC$root.depth_m)
mean = mean(final.data.CC$root.depth_m, na.rm = TRUE)
std = sd(final.data.CC$root.depth_m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$root.depth_m[which(final.data.CC$root.depth_m <Tmin | final.data.CC$root.depth_m > Tmax)])
# remove depth 2.426850 2.426850 2.598333 2.598333 2.682500 2.682500 2.915000
# percent removed 
table(is.na(final.data.CC$root.depth_m))
958-478
15/480*100 #3.13%

hist(final.data.CC$RTD.groot.cahill.merge)
boxplot(final.data.CC$RTD.groot.cahill.merge)
mean = mean(final.data.CC$RTD.groot.cahill.merge, na.rm = TRUE)
std = sd(final.data.CC$RTD.groot.cahill.merge, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$RTD.groot.cahill.merge[which(final.data.CC$RTD.groot.cahill.merge <Tmin | final.data.CC$RTD.groot.cahill.merge > Tmax)])
# remove RTD 1.19455 1.19455 1.19455 1.19455
# percent removed 
table(is.na(final.data.CC$RTD.groot.cahill.merge))
958-581
1/377*100 #0.27%

hist(final.data.CC$SRL.groot.cahill.merge)
boxplot(final.data.CC$SRL.groot.cahill.merge)
mean = mean(final.data.CC$SRL.groot.cahill.merge, na.rm = TRUE)
std = sd(final.data.CC$SRL.groot.cahill.merge, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$SRL.groot.cahill.merge[which(final.data.CC$SRL.groot.cahill.merge <Tmin | final.data.CC$SRL.groot.cahill.merge > Tmax)])
# remove SRL 471.2364 471.2364 471.2364 601.5858 627.5050
# percent removed 
table(is.na(final.data.CC$SRL.groot.cahill.merge))
958-575
10/383*100 #2.61%

hist(final.data.CC$rootDiam.mm)
boxplot(final.data.CC$rootDiam.mm)
mean = mean(final.data.CC$rootDiam.mm, na.rm = TRUE)
std = sd(final.data.CC$rootDiam.mm, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$rootDiam.mm[which(final.data.CC$rootDiam.mm <Tmin | final.data.CC$rootDiam.mm > Tmax)])
# remove diam 1.3965
# percent removed 
table(is.na(final.data.CC$rootDiam.mm))
958-556
4/402*100 #0.995%

hist(final.data.CC$RMF.g.g)
boxplot(final.data.CC$RMF.g.g)
mean = mean(final.data.CC$RMF.g.g, na.rm = TRUE)
std = sd(final.data.CC$RMF.g.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(final.data.CC$RMF.g.g[which(final.data.CC$RMF.g.g <Tmin | final.data.CC$RMF.g.g > Tmax)])
# remove RMF: none
# percent removed 
table(is.na(final.data.CC$RMF.g.g))
958-556
4/402*100 #0.995%


# 25 total removed because some removed with multiple traits
final.data.CC.2 = final.data.CC %>%
  filter(round(leafN.mg.g, 5) != 54.71597) %>% # 2
  filter(round(height.m, 5) < 1.800048) %>% # 3
  filter(round(rootN.mg.g, 5) < 31.17241) %>% # 5
  filter(round(root.depth_m, 5) < 2.426850) %>% # 7
  filter(round(SRL.groot.cahill.merge, 5) < 471.2364) # 5



