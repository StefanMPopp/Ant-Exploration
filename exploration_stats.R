#### Ant exploration: (G)LMMs for exploration, inter & inter day comparisons ####
# For "Ant colonies explore novel environments with increased activity and slower, curvier walks."
# Stefan Popp & Anna Dornhaus
# March 2023

list.of.packages <- c("lme4", "lmerTest","readr","emmeans")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(lme4) # GLMMs
library(lmerTest) # adds p values to lme4 output
library(emmeans) # Tukey's test on LMMs
library(dplyr) # For Exploration analysis (antmeters walked by col & day)
library(readr) # Reading in data file easily

fileDir = "/home/stefan/Documents/MATLAB/AntExploration"
datWalls <- read_csv(paste(fileDir,"/dataWithWalls.txt", sep=''))
datWalls = datWalls[c("day","col","s","nestDisp","v","st","t","idUni")]
dat <- read_csv(paste(fileDir,"/dataClean.txt",sep=''))
dat = dat[c("day","col","s","nestDisp","v","st","t","idUni")]

#### Exploration (ant meters walked) #####
mmWalked <- datWalls %>%
  group_by(col, day) %>%
  summarize(s = sum(s,na.rm=TRUE))
activ = lmer(s ~ day + (1|col), data=mmWalked) # needs sum of walked dist by col&day
summary(activ)
activFactor = lmer(s ~ factor(day) + (1|col), data=mmWalked) # needs sum of walked dist by col&day
summary(emmeans(activFactor, specs = pairwise ~ day, adjust = "tukey"))

#### Dispersion (nestDisp) ####
expl = lmer(log(nestDisp) ~ day + (1|col), data=datWalls)
summary(expl)
explFactor = lmer(log(nestDisp) ~ factor(day) + (1|col), data=datWalls)
summary(emmeans(explFactor, specs = pairwise ~ day, adjust = "tukey")) #




#### Speed & Straightness BETWEEN days & nest distances ####
### Between days
## Both near & far from nest
# Speed
vBetween = lmer(v ~ day + (1|col), data=dat)
summary(vBetween)
vBetweenFactor = lmer(v ~ factor(day) + (1|col), data=dat) # For Tukey's, as.factor needed
summary(emmeans(vBetweenFactor, specs = pairwise ~ day, adjust = "tukey")) # gives NAs

# Straightness
data = dat[seq(1,nrow(dat), 50),]
stBetween = glmer(st ~ day + (1|col), data=dat, family=binomial())
summary(stBetween)
stBetweenFactor = glmer(st ~ factor(day) + (1|col), data=data, family=binomial()) # For Tukey's, as.factor needed
summary(emmeans(stBetweenFactor, specs = pairwise ~ day, adjust = "tukey"))


## Only near the nest
# Speed
vBetweenNear = lmer(v ~ day + (1|col), data=dat[dat$nestDisp < 1000,])
summary(vBetweenNear)
vBetweenNearFactor = lmer(v ~ factor(day) + (1|col), data=dat[dat$nestDisp < 1000,]) # For Tukey's, as.factor needed
summary(emmeans(vBetweenNearFactor, specs = pairwise ~ day, adjust = "tukey")) # gives NAs

# Straightness
stBetween = lmer(st ~ day + (1|col), data=dat[dat$nestDisp < 1000,])
summary(stBetween)
stBetweenNearFactor = lmer(v ~ factor(day) + (1|col), data=dat[dat$nestDisp < 1000,]) # For Tukey's, as.factor needed
summary(emmeans(stBetweenNearFactor, specs = pairwise ~ day, adjust = "tukey")) # gives NAs


## Only far from the nest
# Speed
vBetween = lmer(v ~ day + (1|col), data=dat[dat$nestDisp >= 1000,])
summary(vBetween)
vBetweenFarFactor = lmer(v ~ factor(day) + (1|col), data=dat[dat$nestDisp > 1000,]) # For Tukey's, as.factor needed
summary(emmeans(vBetweenFarFactor, specs = pairwise ~ day, adjust = "tukey")) # gives NAs

# Straightness
stBetween = lmer(st ~ day + (1|col), data=dat[dat$nestDisp >= 1000,])
summary(stBetween)
stBetweenFarFactor = lmer(v ~ factor(day) + (1|col), data=dat[dat$nestDisp > 1000,]) # For Tukey's, as.factor needed
summary(emmeans(stBetweenFarFactor, specs = pairwise ~ day, adjust = "tukey")) # gives NAs


### Between NEAR & FAR within days
dat$nearFar = 'near'
dat$nearFar[dat$nestDisp >= 1000] = 'far'
## Day 1
# Speed
vNearFar = lmer(v ~ nearFar + (1|col), data=dat[dat$day==1,])
summary(vNearFar)

# Straightness
stNearFar1 = lmer(st ~ nearFar + (1|col), data=dat[dat$day==1,])
summary(stNearFar1)


## Day 2
# Speed
vNearFar2 = lmer(v ~ nearFar + (1|col), data=dat[dat$day==2,])
summary(vNearFar2)

# Straightness
stNearFar3 = lmer(st ~ nearFar + (1|col), data=dat[dat$day==2,])
summary(stNearFar3)


## Day 3
# Speed
vNearFar = lmer(v ~ nearFar + (1|col), data=dat[dat$day==3,])
summary(vNearFar)

# Straightness
stNearFar = lmer(st ~ nearFar + (1|col), data=dat[dat$day==3,])
summary(stNearFar)




#### Speed & Straightness WITHIN days ####
# Speed
vWithin1 = lmer(v ~ t + (1|col), data=dat[dat$day==1,])
summary(vWithin1)

vWithin2 = lmer(v ~ t + (1|col), data=dat[dat$day==2,])
summary(vWithin2)

vWithin3 = lmer(v ~ t + (1|col), data=dat[dat$day==3,])
summary(vWithin3)

# Straightness
stWithin1 = glmer(st ~ t + (1|col), data=dat[dat$day==1,], family=binomial())
summary(stWithin1)

stWithin2 = glmer(st ~ t + (1|col), data=dat[dat$day==2,], family=binomial())
summary(stWithin2)

stWithin3 = glmer(st ~ t + (1|col), data=dat[dat$day==3,], family=binomial())
summary(stWithin3)

