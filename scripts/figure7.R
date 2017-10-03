# This script will be used to produce figure 7 for Rosen et al. T1QA
# It will explore univariate relationships between data quality estimates and 
# local FS CT estimates

## Load the data
source('/home/adrose/T1QA/scripts/galton/loadGo1Data.R')
detachAllPackages()
set.seed(16)
tbvData <- read.csv('/home/adrose/dataPrepForHiLoPaper/data/preRaw/t1/n1601_antsCtVol.csv')
## Now load all of the freesurfer values
fsVol <- read.csv('/home/adrose/qapQA/data/n1601_freesurferVol_20161220.csv')
fsCt <- read.csv('/home/adrose/qapQA/data/n1601_freesurferCt_20161220.csv')
fsVals <- merge(fsVol, fsCt, by=c('bblid', 'scanid'))
## Source all functions
source('/home/adrose/RosenT1QA/scripts/figures78Functions.R')

## Load library(s) we will need
install_load('caret', 'lme4', 'bda', 'ggplot2', 'R.matlab')

## Now lets prep the data
## Now create the training data set and create the outcomes for all of the training data sets
raw.lme.data <- merge(isolatedVars, manualQAData2, by='bblid')
raw.lme.data <- merge(raw.lme.data, tbvData, by='bblid')
raw.lme.data$averageRating.x <- as.numeric(as.character(raw.lme.data$averageRating.x))
raw.lme.data$averageRating.x[raw.lme.data$averageRating.x>1] <- 1
load('/home/adrose/RosenT1QA/data/folds/foldsToUse.RData')
index <- unlist(folds[1])
trainingData <- raw.lme.data[index,]
validationData <- raw.lme.data[-index,]

## Now create our outcomes
# Now make sure our data is organized for the functions 
trainingData$oneVsTwoOutcome <- trainingData$mean_euler
validationData$oneVsTwoOutcome <- validationData$mean_euler

## Now merge our scaled data values with the original data values
all.train.data <- merge(mergedQAP, trainingData, by='bblid')
all.train.data <- merge(all.train.data, fsVals, by='bblid')
all.valid.data <- merge(mergedQAP, validationData, by='bblid')
all.valid.data <- merge(all.valid.data, fsVals, by='bblid')

# Remove 0 scans 
all.train.data <- all.train.data[which(all.train.data$averageRating!=0),]
all.valid.data <- all.valid.data[which(all.valid.data$averageRating!=0),]

# Now create our z scores
tmp <- all.train.data[,-seq(2862, 2997)[1:38]]
fsCTVals <- pvalLoop('_thickness', tmp, correct=TRUE)
fsCTVals <- fsCTVals[-grep('ean', fsCTVals[,1]),]
rm(tmp)

## Now create our color values to export to ITK snap
ctColors <- returnPosNegAndNeuColorScale(fsCTVals[,2], colorScaleNeg=c('blue', 'light blue'),colorScalePos=c('yellow', 'red'))[-1,]
ctColors[,8] <- fsCTVals[,1]
ctColors <- cbind(ctColors, fsCTVals[,2])

# Now I need to save these color scales and the other thing
writeMat('fsctColorScale.mat', vals=ctColors)

# Now do the validation data down here
static <- all.train.data
all.train.data <- all.valid.data

tmp <- all.train.data[,-seq(2862, 2997)[1:38]]
fsCTVals <- pvalLoop('_thickness', tmp, correct=TRUE)
fsCTVals <- fsCTVals[-grep('ean', fsCTVals[,1]),]
rm(tmp)

## Now create our color values to export to ITK snap
ctColors <- returnPosNegAndNeuColorScale(fsCTVals[,2], colorScaleNeg=c('blue', 'light blue'),colorScalePos=c('yellow', 'red'))[-1,]
ctColors[,8] <- fsCTVals[,1]
ctColors <- cbind(ctColors, fsCTVals[,2])

# Now I need to save these color scales and the other thing
writeMat('fsctvalidColorScale.mat', vals=ctColors)

# Now do MGI
source('/home/adrose/T1QA/scripts/galton/loadMgiData.R')
all.train.data <- mergedQAP
all.train.data <- all.train.data[which(all.train.data$averageRating!=0),]
all.train.data$oneVsTwoOutcome <- all.train.data$mean_euler
all.train.data$ageAtGo1Scan <- all.train.data$age
all.train.data$sex <- all.train.data$Gender

fsCTVals <- pvalLoop('_thickness', all.train.data, correct=TRUE)
fsCTVals <- fsCTVals[-grep('ean', fsCTVals[,1]),]
## Now create our color values to export to ITK snap
ctColors <- returnPosNegAndNeuColorScale(fsCTVals[,2], colorScaleNeg=c('blue', 'light blue'),colorScalePos=c('yellow', 'red'))[-1,]
ctColors[,8] <- fsCTVals[,1]
ctColors <- cbind(ctColors, fsCTVals[,2])
# Now I need to save these color scales and the other thing
writeMat('fsctColorScaleMGI.mat', vals=ctColors)
