# This script will:
# - pull the raw behavioral data and eyetracking data
# - initial cleaning of the data (remove practice blocks, remove subjects because of too many errors and "other")
## Note: this script should be used only for those experiments with Non-Switch cues (i.e. not for Study 1a)

#### Adjust parameters ####
numPracBlocks <- 3 # practice blocks are removed from analyses

# offset stim changes depending on experiment
offsetTime <- 1000 # Study1b

# select relevant experiment
# excludeOther <- c(279, 291, 307, 316) # study 1a
excludeOther <- 617 # study 1b
# excludeOther <- c(708,714) # study 1c

# Exclude pupil outliers or other issues with pupils
#excludePupilOutliers <- c(250, 287) #study 1a
excludePupilOutliers <- c(607, 609, 613, 614) #study 1b
#excludePupilOutliers <- 434 #study 2


# Adjust paths based on data being analyzed
behPath <- '~/Dropbox/Studies/Dissertation/Study1/Study1b/data/beh/txt/'
edfPath <- '~/Dropbox/Studies/Dissertation/Study1/Study1b/data/edf/'
surveyPath <- '~/Dropbox/Studies/Dissertation/Study1/Study1b/data/surveyMeasures/'
outputPath <- '~/Dropbox/Studies/Dissertation/output/Study1b/'
codePath <- '~/Dropbox/Studies/Dissertation/code/reAnalysis2017/'

library('itrackR')
library('edfR')
library('dplyr')
library('dtplyr')
library('ggplot2')
source(paste(codePath, 'load_files.R', sep=''))


### Loading in samples ###
setwd(edfPath)
edf_files <- list.files(path=edfPath, pattern='*.edf')

z <- itrackr(edfs=edf_files)

z <- find_messages(z, c('STIMONSET','RESPONSE', 'BASELINE_START', 'BASELINE_END'), c('STIMONSET','RESPONSE', 'BASELINE_START', 'BASELINE_END'),timestamp = T)
z <- set_index(z,c('Block','Set','Trial'),c('^BLOCK [0-9]*','^SET [0-9]*','TRIAL [0-9]*'),numeric.only = T)
beh <- load_files(path=behPath,pattern='*.txt')
z <- add_behdata(z,beh)

innercoords <- radialCoords(x=512,y=384, numpoints=12, radius=200)
angles <- roiFlower(12)

z <- makeROIs(z, innercoords, shape='ellipse', xradius=45, yradius=90,
              angles=angles[c(1,2,3,4,5,6,7,8,9,10,11,12)])
centercoords <- matrix(c(512,384),nrow=1)
z <- makeROIs(z,centercoords,shapes='circle',radius=65, names=13, append=T)

z <- drift_correct(z,vars='Block',threshold=15)

z <- calcHits(z)
z <- mapROIs(z,names=c('target','distractor','neutral'),indicators=c('Targetpos','Distpos','Neutralpos'))

fixes <- eyemerge(z,'fixations', all.rois=TRUE)
saccs <- eyemerge(z,'saccades', all.rois=TRUE)


# remove fixations that occurred prior to stimulus onset and after stimulus offset, and those that occur after participants have responded
fixes$relOnsetTime <- fixes$sttime - fixes$STIMONSET
fixes$relOffsetTime <- fixes$sttime-(fixes$STIMONSET + 1000) #In Study 1b, stimulus offset occurs 1000ms after onset.
fixes$dur <- fixes$entime-fixes$sttime
fixesClean <- filter(fixes, Block>numPracBlocks, Blocktype>=3, relOnsetTime>0, relOffsetTime<0, sttime<RESPONSE)


# Calculate saccade velocity
saccs$xdistance <- abs(saccs$gstx-saccs$genx)
saccs$ydistance <- abs(saccs$gsty-saccs$geny)
saccs$distance <- sqrt((saccs$xdistance)^2 + (saccs$ydistance)^2)
saccs$velocity <- saccs$distance/(saccs$entime-saccs$sttime)

# remove saccades that occur prior to target onset and after target offset
saccs$relOnsetTime <- saccs$sttime-saccs$STIMONSET
saccs$relOffsetTime <- saccs$sttime-(saccs$STIMONSET + 1000) 
saccsClean <- filter(saccs, Block>numPracBlocks, Blocktype>=3, relOnsetTime>0, relOffsetTime<0, sttime<RESPONSE)

# pull indices of first saccade in trial
getfirstSaccade <- group_by(saccsClean, ID, eyetrial)
getfirstSaccade <- subset(getfirstSaccade, !duplicated(getfirstSaccade[,1:4]))
firstSaccade <- select(getfirstSaccade, ID, Block, Set, Trial, eyetrial, target_end_hit, distractor_end_hit, neutral_end_hit, roi_start_13, roi_end_13,
                       sttime, entime, distance, velocity, relOnsetTime) %>%
  rename(firstSacc_tar=target_end_hit, firstSacc_dist=distractor_end_hit, firstSacc_neut=neutral_end_hit, firstSacc_startFix=roi_start_13,
           firstSacc_fix=roi_end_13, firstSacc_sttime=sttime, firstSacc_entime=entime, firstSacc_distance=distance, firstSacc_velocity=velocity, 
           firstSacc_relOnsetTime=relOnsetTime)

# pull indices of first saccade to target in trial
getfirstSaccadeToTarget <- group_by(saccsClean, ID, eyetrial) %>%
  filter(target_end_hit==1)
getfirstSaccadeToTarget <- subset(getfirstSaccadeToTarget, !duplicated(getfirstSaccadeToTarget[,1:4]))
firstSaccadeToTarget <- select(getfirstSaccadeToTarget, ID, Block, Set, Trial, eyetrial, target_end_hit,
                       sttime, entime, distance, velocity, relOnsetTime) %>%
  rename(firstTarSacc_tar=target_end_hit, firstTarSacc_sttime=sttime, firstTarSacc_entime=entime, firstTarSacc_distance=distance, firstTarSacc_velocity=velocity, 
         firstTarSacc_relOnsetTime=relOnsetTime)

# pull indices of first saccade to distractor in each trial
getfirstSaccadeToDistractor <- group_by(saccsClean, ID, eyetrial) %>%
  filter(distractor_end_hit==1)
getfirstSaccadeToDistractor <- subset(getfirstSaccadeToDistractor, !duplicated(getfirstSaccadeToDistractor[,1:4]))
firstSaccadeToDistractor <- select(getfirstSaccadeToDistractor, ID, Block, Set, Trial, eyetrial, distractor_end_hit,
                               sttime, entime, distance, velocity, relOnsetTime) %>%
  rename(firstDistSacc_dist=distractor_end_hit, firstDistSacc_sttime=sttime, firstDistSacc_entime=entime, firstDistSacc_distance=distance, firstDistSacc_velocity=velocity, 
         firstDistSacc_relOnsetTime=relOnsetTime)

# pull indices of first saccade to neutral in trial
getfirstSaccadeToNeutral <- group_by(saccsClean, ID, eyetrial) %>%
  filter(neutral_end_hit==1)
getfirstSaccadeToNeutral <- subset(getfirstSaccadeToNeutral, !duplicated(getfirstSaccadeToNeutral[,1:4]))
firstSaccadeToNeutral <- select(getfirstSaccadeToNeutral, ID, Block, Set, Trial, eyetrial, neutral_end_hit,
                                   sttime, entime, distance, velocity, relOnsetTime) %>%
  rename(firstNeutSacc_neut=neutral_end_hit, firstNeutSacc_sttime=sttime, firstNeutSacc_entime=entime, firstNeutSacc_distance=distance, firstNeutSacc_velocity=velocity, 
         firstNeutSacc_relOnsetTime=relOnsetTime)

# merge all saccade data together
allSaccs <- merge(firstSaccade, firstSaccadeToTarget, by=c('ID','Block','Set','Trial','eyetrial'), all=TRUE)
allSaccs <- merge(allSaccs, firstSaccadeToDistractor, by=c('ID','Block','Set','Trial','eyetrial'), all=TRUE)
allSaccs <- merge(allSaccs, firstSaccadeToNeutral, by=c('ID','Block','Set','Trial','eyetrial'), all=TRUE)


# Create similar dataframe for fixations by trial (first fixation, first to target, distractor, and neutral)
getfirstFix <- group_by(fixesClean, ID, eyetrial)
getfirstFix <- subset(getfirstFix, !duplicated(getfirstFix[,1:4]))
firstFix <- select(getfirstFix, ID, Block, Set, Trial, eyetrial, target_hit, distractor_hit, neutral_hit, roi_13,
                       sttime, entime, relOnsetTime, dur) %>%
  rename(firstFix_target=target_hit, firstFix_distractor=distractor_hit, firstFix_neutral=neutral_hit, firstFix_fix=roi_13, firstFix_sttime=sttime,
         firstFix_entime=entime, firstFix_relOnsetTime=relOnsetTime, firstFix_dur=dur)

getfirstTarFix <- group_by(fixesClean, ID, eyetrial) %>%
  filter(target_hit==1)
getfirstTarFix <- subset(getfirstTarFix, !duplicated(getfirstTarFix[,1:4]))
firstTarFix <- select(getfirstTarFix, ID, Block, Set, Trial, eyetrial, target_hit,
                   sttime, entime, relOnsetTime, dur) %>%
  rename(firstTarFix_tar=target_hit, firstTarFix_sttime=sttime,
         firstTarFix_entime=entime, firstTarFix_relOnsetTime=relOnsetTime, firstTarFix_dur=dur)

getfirstDistFix <- group_by(fixesClean, ID, eyetrial) %>%
  filter(distractor_hit==1)
getfirstDistFix <- subset(getfirstDistFix, !duplicated(getfirstDistFix[,1:4]))
firstDistFix <- select(getfirstDistFix, ID, Block, Set, Trial, eyetrial, distractor_hit,
                   sttime, entime, relOnsetTime, dur) %>%
  rename(firstDistFix_dist=distractor_hit, firstDistFix_sttime=sttime,
         firstDistFix_entime=entime, firstDistFix_relOnsetTime=relOnsetTime, firstDistFix_dur=dur)

getfirstNeutFix <- group_by(fixesClean, ID, eyetrial) %>%
  filter(neutral_hit==1)
getfirstNeutFix <- subset(getfirstNeutFix, !duplicated(getfirstNeutFix[,1:4]))
firstNeutFix <- select(getfirstNeutFix, ID, Block, Set, Trial, eyetrial, neutral_hit, 
                   sttime, entime, relOnsetTime, dur) %>%
  rename(firstNeutFix_neut=neutral_hit, firstNeutFix_sttime=sttime,
         firstNeutFix_entime=entime, firstNeutFix_relOnsetTime=relOnsetTime, firstNeutFix_dur=dur)

allFixes <- merge(firstFix, firstTarFix, by=c('ID','Block','Set','Trial','eyetrial'), all=TRUE)
allFixes <- merge(allFixes, firstDistFix, by=c('ID','Block','Set','Trial','eyetrial'), all=TRUE)
allFixes <- merge(allFixes, firstNeutFix, by=c('ID','Block','Set','Trial','eyetrial'), all=TRUE)


# find subjects with too many errors:
Errors <- dplyr::group_by(beh, ID) %>%
  filter(Block>numPracBlocks) %>%
  summarise(errors=mean(Error, na.rm=TRUE))

meanErrorCalc <- filter(Errors, !(ID %in% c(excludeOther)))
meanErrors <- mean(meanErrorCalc$errors)
sdErrors <- sd(meanErrorCalc$errors)

maxErrors <- meanErrors + 3*sdErrors
excludeMaxErrors <- Errors$ID[Errors$errors>maxErrors]


# Merge files and exclude subjects with pupil issues, other issues, and too many errors
allEyes <- merge(allSaccs, allFixes, by=c('ID','Block','Set','Trial','eyetrial'), all=TRUE) %>%
  filter(!(ID %in% c(excludePupilOutliers)) &
             !(ID %in% c(excludeOther)) &
           !(ID %in% c(excludeMaxErrors)))

allEyes$firstTarSacc_tar[is.na(allEyes$firstTarSacc_tar)] <- 0
allEyes$firstDistSacc_dist[is.na(allEyes$firstDistSacc_dist)] <- 0
allEyes$firstNeutSacc_neut[is.na(allEyes$firstNeutSacc_neut)] <- 0
allEyes$firstTarFix_tar[is.na(allEyes$firstTarFix_tar)] <- 0
allEyes$firstDistFix_dist[is.na(allEyes$firstDistFix_dist)] <- 0
allEyes$firstNeutFix_neut[is.na(allEyes$firstNeutFix_neut)] <- 0



save(allEyes, file=paste(outputPath, 'eyetrackingData_byTrial.Rda', sep=''))
