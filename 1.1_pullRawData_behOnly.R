# This script will:
# - pull the raw behavioral data 
# - initial cleaning of the data (remove practice blocks, remove subjects because of too many errors and "other")


# Adjust paths based on data being analyzed
behPath <- '~/Dropbox/Studies/Dissertation/Study1/Study1b/data/beh/txt/'
surveyPath <- '~/Dropbox/Studies/Dissertation/Study1/Study1b/data/surveyMeasures/'
outputPath <- '~/Dropbox/Studies/Dissertation/output/Study1b/'
codePath <- '~/Dropbox/Studies/Dissertation/code/reAnalysis2017/'


#### Adjust parameters ####
numPracBlocks <- 3 # practice blocks are removed from analyses

# select relevant experiment
# excludeOther <- c(279, 291, 307, 316) # study 1a
excludeOther <- 617 # study 1b
# excludeOther <- c(708,714) # study 1c

# excludeErrorsPercent <- .15 #threshold for excluding subjects due to poor accuracy (only use if want to set threshold other than 3SDs)

###############

library(itrackR)
library(edfR)
library(dtplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(plyr)
library(dplyr)
source(paste(codePath, 'load_files.R', sep=''))

### Loading in samples ###
beh <- load_files(path=behPath,pattern='*.txt')

save(beh, file=paste(outputPath, 'behByTrial.Rda', sep=''))
load(paste(outputPath, 'behByTrial.Rda', sep='')) 

#load and add demographic data
demog <- read.csv(paste(surveyPath, 'demographicsCaffeineSleep.csv', sep=''), header=T)
allData <- join(beh, demog, by='ID')

detach(package:plyr)
Errors <- dplyr::group_by(allData, ID) %>%
  summarise(errors=mean(Error, na.rm=TRUE))

excludeErrors <- Errors$ID[Errors$errors>excludeErrorsPercent]

meanErrorCalc <- filter(Errors, !(ID %in% c(excludeErrors)) & !(ID %in% c(excludeOther)))
meanErrors <- mean(meanErrorCalc$errors)
sdErrors <- sd(meanErrorCalc$errors)

maxErrors <- meanErrors + 3*sdErrors
excludeMaxErrors <- Errors$ID[Errors$errors>maxErrors]



# includes everything except ineligible participants and practice trials
cleanData_behAnalyses <- filter(allData,
                                  Block>numPracBlocks &
                                  #!(ID %in% c(excludeErrors)) &
                                  !(ID %in% c(excludeOther)))
# cleanData_behAnalyses_eRate <- filter(allData,
#                                         Block>numPracBlocks &
#                                         #!(ID %in% c(excludeErrors)) &
#                                         !(ID %in% c(excludeOther)) &
#                                         !(ID %in% c(excludeMaxErrors)))




save(cleanData_behAnalyses, file=paste(outputPath, 'cleanDataBehAnalyses.Rda', sep='')) 
save(cleanData_behAnalyses_eRate, file=paste(outputPath, 'cleanDataBehAnalyses_eRate.Rda', sep='')) 




## Examining Data
beh_noErrors <- filter(cleanData_behAnalyses, Error==0)
hist(beh_noErrors$RT)
