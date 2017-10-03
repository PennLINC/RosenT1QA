# This script will be used to produce figure 8 for Rosen et al. T1QA
# It will explore cortical thickness and age relationships mediated by 
# data quality estimates

## Load the data
source('/home/adrose/T1QA/scripts/galton/loadGo1Data.R')
detachAllPackages()
set.seed(16)
tbvData <- read.csv('/home/adrose/dataPrepForHiLoPaper/data/preRaw/t1/n1601_antsCtVol.csv')
fsVol <- read.csv('/home/adrose/qapQA/data/n1601_freesurferVol_20161220.csv')
fsCt <- read.csv('/home/adrose/qapQA/data/n1601_freesurferCt_20161220.csv')
fsVals <- merge(fsVol, fsCt, by=c('bblid', 'scanid'))

## Source all functions
source('/home/adrose/RosenT1QA/scripts/figures78Functions.R')

## Load library(s) we will need
install_load('caret', 'lme4', 'bda', 'ggplot2', 'R.matlab')

## Now lets prep the data
raw.lme.data <- merge(isolatedVars, manualQAData2, by='bblid')
raw.lme.data <- merge(raw.lme.data, tbvData, by='bblid')
raw.lme.data$averageRating.x <- as.numeric(as.character(raw.lme.data$averageRating.x))
raw.lme.data$averageRating.x[raw.lme.data$averageRating.x>1] <- 1
load('/home/adrose/RosenT1QA/data/folds/foldsToUse.RData')
index <- unlist(folds[1])
trainingData <- raw.lme.data[index,]
validationData <- raw.lme.data[-index,]

## Now merge our scaled data values with the original data values
all.train.data <- merge(mergedQAP, trainingData, by='bblid')
all.train.data <- merge(all.train.data, fsVals, by='bblid')
all.valid.data <- merge(mergedQAP, validationData, by='bblid')
all.valid.data <- merge(all.valid.data, fsVals, by='bblid')

## Now remove the 0 values because they do weird things
all.train.data <- all.train.data[which(all.train.data$averageRating!=0),]
all.valid.data <- all.valid.data[which(all.valid.data$averageRating!=0),]
vals <- grep('_thickness', names(all.train.data))
vals <- vals[-1]
zScoreCT <- NULL
binVals <- NULL
for(i in vals){
    foo <- mediation.test(mv=all.train.data$mean_euler, iv=all.train.data$ageAtGo1Scan, dv= all.train.data[,i])$Sobel[1]
    bar <- mediation.test(mv=all.train.data$mean_euler, iv=all.train.data$ageAtGo1Scan, dv= all.train.data[,i])$Sobel[2]
    toAppend <- c(names(all.train.data)[i], foo)
    toAppend2 <- c(names(all.train.data)[i], bar)
    zScoreCT <- rbind(zScoreCT, toAppend)
    binVals <- rbind(binVals, toAppend2)
}
binValsApply <- rep(0, length(vals))
binValsApply[which(p.adjust(binVals[,2], method='fdr')<.05)] <- 1
zScoreCT[,2] <- as.numeric(zScoreCT[,2]) * binValsApply
zScoreCT <- zScoreCT[order(as.numeric(zScoreCT[,2])),]
zScoreCT <- zScoreCT[-c(1, grep('ean', zScoreCT[,1])),]

## Now create our color values to export to ITK snap
ctColors <- returnPosNegAndNeuColorScale(zScoreCT[,2], colorScaleNeg=c('blue', 'light blue'),colorScalePos=c('yellow', 'pink'))[-1,]
ctColors[,8] <- zScoreCT[,1]
ctColors <- cbind(ctColors, c(zScoreCT[,2]))
ctColors[ctColors=="NaN"] <- 190
# Now I need to save these color scales and the other thing
writeMat('fsctColorScale.mat', vals=ctColors)

# Now do the validation data down here
static <- all.train.data
all.train.data <- all.valid.data

zScoreCT <- NULL
binVals <- NULL
for(i in vals){
    foo <- mediation.test(mv=all.train.data$mean_euler, iv=all.train.data$ageAtGo1Scan, dv= all.train.data[,i])$Sobel[1]
    bar <- mediation.test(mv=all.train.data$mean_euler, iv=all.train.data$ageAtGo1Scan, dv= all.train.data[,i])$Sobel[2]
    toAppend <- c(names(all.train.data)[i], foo)
    toAppend2 <- c(names(all.train.data)[i], bar)
    zScoreCT <- rbind(zScoreCT, toAppend)
    binVals <- rbind(binVals, toAppend2)
}
binValsApply <- rep(0, length(vals))
binValsApply[which(p.adjust(binVals[,2], method='fdr')<.05)] <- 1
zScoreCT[,2] <- as.numeric(zScoreCT[,2]) * binValsApply
zScoreCT <- zScoreCT[order(as.numeric(zScoreCT[,2])),]
zScoreCT <- zScoreCT[-c(1, grep('ean', zScoreCT[,1])),]

## Now create our color values to export to ITK snap
ctColors <- returnPosNegAndNeuColorScale(zScoreCT[,2], colorScalePos=c('blue', 'light blue'),colorScaleNeg=c('yellow', 'pink'))[-1,]
ctColors[,8] <- zScoreCT[,1]
ctColors <- cbind(ctColors, zScoreCT[,2])
ctColors[ctColors=="NaN"] <- 190
# Now I need to save these color scales and the other thing
writeMat('fsctColorScaleTest.mat', vals=ctColors)
