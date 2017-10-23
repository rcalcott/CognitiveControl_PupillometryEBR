outputPath <- '~/Dropbox/Studies/Dissertation/output/Study1b/'


excludeRTmax <- 2000
excludeRTmin <- 100
excludeRTmaxSDs <- 3

library(itrackR)
library(edfR)
library(dplyr)
library(dtplyr)
library(ggplot2)
library(lme4)
library(haven)
library(lmerTest)
library(effects)
library(psych)

load(paste(outputPath, 'cleanDataBehAnalyses.Rda', sep='')) 


#### RT Analyses - Trial-by-Trial RTs ####
cleanData <- filter(cleanData_behAnalyses, Error==0)

# Create cue variable
cleanData$Cue <- NA
cleanData$Cue[cleanData$Trial==1 & cleanData$Set!=1 & cleanData$Blocktype!=2] <- 1
cleanData$Cue[cleanData$Trial>=2 & cleanData$Blocktype!=2] <- 2


findOutliers <- group_by(cleanData, ID) %>%
  summarise(meanRT=mean(RT, na.rm=TRUE),
            sdRT=sd(RT, na.rm=TRUE))

findOutliers$exclude <- findOutliers$meanRT + excludeRTmaxSDs*findOutliers$sdRT
findOutliers$excludeMin <- findOutliers$meanRT - excludeRTmaxSDs*findOutliers$sdRT
cleanData <- merge(cleanData, findOutliers, by='ID', sort=FALSE)
cleanData$outlierMax <- cleanData$exclude-cleanData$RT
cleanData$outlierMin <- cleanData$excludeMin-cleanData$RT

cleanData$switchSet <- factor(cleanData$switchSet,
                              levels=c(0,1),
                              labels=c('Non-Switch Set','Switch Set'))

cleanData$Cue <- factor(cleanData$Cue,
                        levels=c(2,1),
                        labels=c('Uncued','Cued'))

cleanData$Congruent <- factor(cleanData$Congruent,
                              levels=c(1,2),
                              labels=c('Congruent','Incongruent'))

cleanData$Blocktype <- factor(cleanData$Blocktype,
                              levels=c(3,4),
                              labels=c('Pure Updating','Perseveration-Inhibition'))
cleanData$ID <- factor(cleanData$ID)


cleanData_outliersRemoved <- filter(cleanData, 
                                    RT>=excludeRTmin & RT<=excludeRTmax & outlierMax>0 & outlierMin<0)

model0 <- lmer(RT ~ 1 +
                 (1|ID),
               data=cleanData_outliersRemoved,
               na.action=na.exclude)

# calculate ICC
varcom <- as.data.frame(VarCorr(model0))
L2var <- varcom$vcov[1]
L1var <- varcom$vcov[2]
icc <- L2var/(L2var+L1var)


model1 <- lmer(RT ~ Cue*switchSet*Congruent*Blocktype + 
                 (1|ID),
               data=cleanData_outliersRemoved,
               na.action=na.exclude)

eff1 <- effect('Cue:switchSet:Blocktype', model1)
x1 <- as.data.frame(eff1)
x1$Cue <- relevel(x1$Cue,'Uncued')

limits = aes(ymax = fit + (se), ymin=fit - (se))
dodge = position_dodge(width=0.9)

s=ggplot(x1, aes(x = switchSet, y = fit, fill = Cue))+
  facet_grid(Blocktype~.) +
  geom_bar(stat='identity', position=dodge)+
  geom_errorbar(limits, position=dodge, width=0.25)+
  scale_fill_manual(values=c('orangered3','dodgerblue4'),
                    name="Cue",
                    breaks=c("Uncued", "Cued"),
                    labels=c("Uncued", "Cued")) +
  ylab('RT (ms)')+
  xlab('Set Type')

s + coord_cartesian(ylim=c(650,950))


#### Mean RT Analyses - By Condition ####

summaryStats <- group_by(cleanData_outliersRemoved, ID, Cue, switchSet, Congruent, Blocktype) %>%
  summarise(meanRT = mean(RT, na.rm=TRUE))

model0.3 <- lmer(meanRT ~ 1 +
                   (1|ID),
                 data=summaryStats,
                 na.action=na.exclude)

varcom <- as.data.frame(VarCorr(model0.3))
L2var <- varcom$vcov[1]
L1var <- varcom$vcov[2]
icc <- L2var/(L2var+L1var)


model3 <- lmer(meanRT ~ Cue*switchSet*Congruent*Blocktype + 
                 (1|ID),
               data=summaryStats,
               na.action=na.exclude)

#### Median and RTCV Analyses - By Condition ####
summaryStats <- group_by(cleanData, ID, Cue, switchSet, Congruent, Blocktype) %>%
  summarise(medianRT = median(RT, na.rm=TRUE),
            meanRT = mean(RT, na.rm = TRUE),
            sdRT = sd(RT, na.rm=TRUE))

# Compute RT CV
summaryStats$RTCV <-  summaryStats$sdRT/summaryStats$meanRT


model0.4 <- lmer(medianRT ~ 1 +
                   (1|ID),
                 data=summaryStats,
                 na.action=na.exclude)

model0.5 <- lmer(RTCV ~ 1 +
                   (1|ID),
                 data=summaryStats,
                 na.action=na.exclude)

varcom <- as.data.frame(VarCorr(model0.5))
L2var <- varcom$vcov[1]
L1var <- varcom$vcov[2]
icc <- L2var/(L2var+L1var)



model4 <- lmer(medianRT ~ Cue*switchSet*Congruent*Blocktype + 
                 (1|ID),
               data=summaryStats,
               na.action=na.exclude)

model5 <- lmer(RTCV ~ Cue*switchSet*Congruent*Blocktype + 
                 (1|ID),
               data=summaryStats,
               na.action=na.exclude)


eff1 <- effect('Cue:switchSet:Blocktype', model5)
x1 <- as.data.frame(eff1)
x1$Cue <- relevel(x1$Cue,'Uncued')

limits = aes(ymax = fit + (se), ymin=fit - (se))
dodge = position_dodge(width=0.9)

s=ggplot(x1, aes(x = switchSet, y = fit, fill = Cue))+
  facet_grid(Blocktype~.) +
  geom_bar(stat='identity', position=dodge)+
  geom_errorbar(limits, position=dodge, width=0.25)+
  scale_fill_manual(values=c('orangered3','dodgerblue4'),
                    name="Cue",
                    breaks=c("Uncued", "Cued"),
                    labels=c("Uncued", "Cued")) +
  ylab('RT (ms)')+
  xlab('Set Type')

s + coord_cartesian(ylim=c(650,950))
s + coord_cartesian(ylim=c(0.2,0.3))


eff1 <- effect('Cue:switchSet:Congruent:Blocktype', model5)
x1 <- as.data.frame(eff1)
x1$Cue <- relevel(x1$Cue,'Uncued')

limits = aes(ymax = fit + (se), ymin=fit - (se))
dodge = position_dodge(width=0.9)

s=ggplot(x1, aes(x = switchSet, y = fit, fill = Cue))+
  facet_grid(Blocktype~Congruent) +
  geom_bar(stat='identity', position=dodge)+
  geom_errorbar(limits, position=dodge, width=0.25)+
  scale_fill_manual(values=c('orangered3','dodgerblue4'),
                    name="Cue",
                    breaks=c("Uncued", "Cued"),
                    labels=c("Uncued", "Cued")) +
  ylab('RT (ms)')+
  xlab('Set Type')

s + coord_cartesian(ylim=c(0.2,0.3))


#### Error Rate By-Condition Analyses ####
# Looks at error rates, mean RT, median RT, and RT variability

cleanDataErrors <- cleanData_behAnalyses

# Create cue variable
cleanDataErrors$Cue <- NA
cleanDataErrors$Cue[cleanDataErrors$Trial==1 & cleanDataErrors$Set!=1 & cleanDataErrors$Blocktype!=2] <- 1
cleanDataErrors$Cue[cleanDataErrors$Trial>=2 & cleanDataErrors$Blocktype!=2] <- 2


# create factors
cleanDataErrors$switchSet <- factor(cleanDataErrors$switchSet,
                              levels=c(0,1),
                              labels=c('Non-Switch Set','Switch Set'))

cleanDataErrors$Cue <- factor(cleanDataErrors$Cue,
                        levels=c(2,1),
                        labels=c('uncued','cued'))

cleanDataErrors$Congruent <- factor(cleanDataErrors$Congruent,
                              levels=c(1,2),
                              labels=c('Congruent','Incongruent'))

cleanDataErrors$Blocktype <- factor(cleanDataErrors$Blocktype,
                              levels=c(3,4),
                              labels=c('Pure Updating','Perseveration-Inhibition'))
cleanDataErrors$ID <- factor(cleanDataErrors$ID)


summaryStats <- group_by(cleanDataErrors, ID, Cue, switchSet, Congruent, Blocktype) %>%
  summarise(eRate = mean(Error, na.rm=TRUE))


# Null models
model0.1 <- lmer(eRate ~ 1 +
                   (1|ID),
                 data=summaryStats,
                 na.action=na.exclude)

varcom <- as.data.frame(VarCorr(model0.1))
L2var <- varcom$vcov[1]
L1var <- varcom$vcov[2]
icc <- L2var/(L2var+L1var)


model2 <- lmer(eRate ~ Cue*switchSet*Congruent*Blocktype + 
                 (1|ID),
               data=summaryStats,
               na.action=na.exclude)





