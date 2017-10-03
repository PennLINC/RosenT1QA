# This script will be used to produce figure 6 for Rosen et al. T1QA
# This script will explore the drawbacks of out of sequence motion estimates 
# when comapred to data-driven quality estimates
# This will include information lost from out out sequence motion estimates, differences amongst estimates, and AUC

## Load the data
source('/home/adrose/T1QA/scripts/galton/loadGo1Data.R')
detachAllPackages()
set.seed(16)
motionData <- read.csv('/home/adrose/qapQA/data/qaData/allMergedQAData.csv')

## Now load the library(s)
install_load('caret', 'ggplot2', 'scales', 'pROC')

# Now load the models
raw.lme.data <- merge(isolatedVars, manualQAData2, by='bblid')
raw.lme.data$averageRating.x <- as.numeric(as.character(raw.lme.data$averageRating.x))
raw.lme.data$averageRating.x[raw.lme.data$averageRating.x>1] <- 1
raw.lme.data <- merge(raw.lme.data, mergedQAP, by=intersect(names(raw.lme.data),names(mergedQAP)))

# Now create our full DF
raw.lme.data <- merge(raw.lme.data, motionData, by=intersect(names(raw.lme.data), names(motionData)))
all.train.data <- raw.lme.data
cols <- grep('Meanrelrms', names(all.train.data))
colsCorrect <- cols[c(3, 2, 1, 4)]

val1 <- summarySE(data=all.train.data[,colsCorrect],  measurevar=names(all.train.data[,colsCorrect])[1], na.rm=T)
colnames(val1)[3] <- 'mean'
val2 <- summarySE(data=all.train.data[,colsCorrect],  measurevar=names(all.train.data[,colsCorrect])[2], na.rm=T)
colnames(val2)[3] <- 'mean'
val3 <- summarySE(data=all.train.data[,colsCorrect],  measurevar=names(all.train.data[,colsCorrect])[3], na.rm=T)
colnames(val3)[3] <- 'mean'
val4 <- summarySE(data=all.train.data[,colsCorrect],  measurevar=names(all.train.data[,colsCorrect])[4], na.rm=T)
colnames(val4)[3] <- 'mean'

motionValues <- rbind(val1, val2, val3, val4)
motionValues$.id <- c('2:46', '14:51', '20:19', '43:01')
motionValues$seq <- c('PCASL', 'tfMRI 1', 'tfMRI 2', 'rsfMRI')
motionValues$.id <- factor(motionValues$.id, levels=c('2:46', '14:51', '20:19', '43:01'))
motionValues$seq <- factor(motionValues$seq, levels=c("T1", "PCASL", "tfMRI 1", "tfMRI 2", "rsfMRI"))

# Create the data to plot a survivor
motionCols <- c(2945, 2953, 2980, 2967, 2943)
survivorPercents <- matrix(NA, nrow=nrow(raw.lme.data), ncol=length(motionCols))
z <- 1
for(i in motionCols){
  survivorPercents[which(raw.lme.data[,i]!=1),z] <- raw.lme.data$averageRating[which(raw.lme.data[,i]!=1)]
  z <- z+1
}

# Now produce percents for each
percentVals <- NULL
for(i in 1:length(motionCols)){
    centVals <- table(survivorPercents[,i]) / table(survivorPercents[,5])
    tmpOut <- c(names(raw.lme.data)[motionCols][i] , centVals)
    percentVals <- rbind(percentVals, tmpOut)
}
rownames(percentVals) <- NULL
centValsToPlot <- as.data.frame(percentVals)
centValsToPlot[,2:6] <- apply(centValsToPlot[,2:6], 2, function(x) as.numeric(as.character(x)))
centValsToPlot$V1 <- c("tfMRI 1", "tfMRI 2", "PCASL", "rsfMRI", "T1")
centValsToPlot$V1 <- factor(centValsToPlot$V1, levels=c("T1", "PCASL", "tfMRI 1", "tfMRI 2", "rsfMRI"))
plotData <- melt(centValsToPlot, id.vars='V1')


# Now make the plot
survivorPlot <- ggplot(plotData, aes(x=V1, y=value*100, group=variable)) +
  geom_line(aes(size=1.5)) +
  labs(title='', x='Sequence', y="% Remaining") +
  theme_bw()+
  theme(legend.position="none",
    axis.text.x = element_text(angle=90, size=16, face='bold'),
    axis.text.y = element_text(angle=0, size=16, face='bold'),
    axis.ticks.x=element_blank(),
    axis.text.y = element_text(size=16, face="bold"),
    axis.title=element_text(size=30,face="bold")) + 
  annotate("text", x=5.2, y=45, label=0, parse=T, hjust=c(1), size=8) + 
  annotate("text", x=5.2, y=63, label=1, parse=T, hjust=c(1), size=8) + 
  annotate("text", x=5.5, y=79, label=1.3, parse=T, hjust=c(1), size=8) +
  annotate("text", x=5.5, y=85, label=1.6, parse=T, hjust=c(1), size=8) +
  annotate("text", x=5.2, y=92, label=2, parse=T, hjust=c(1), size=8)

# Now produce the AUC vals 
motionCols <- c(307, 2954, 2981, 2968, 2943)
raw.lme.data <- raw.lme.data[complete.cases(raw.lme.data[,motionCols]),]
aucVals <- NULL
for(i in motionCols){
  roc.tmp <- roc(raw.lme.data$averageRating.x ~ raw.lme.data[,i])
  auc.tmp <- auc(roc.tmp)
  outTmp <- c(names(raw.lme.data)[i], auc.tmp, '0 vs !0', 'Training')
  aucVals <- rbind(aucVals, outTmp)
}

# Produce a mean motion vs time bar plot
meanMotionPlot <- ggplot(motionValues, aes(x=seq, y=mean, fill=seq)) +
  geom_bar(stat='identity', position=position_dodge(), size=.1) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
    width = .2, position=position_dodge(.9)) +
  theme_bw() +
  labs(title='', y='Mean Relative RMS', x='Sequence') +
  coord_cartesian(ylim=c(.06,.17)) +
  theme(text=element_text(size=20), axis.text.x = element_text(angle = 0)) + 
  theme(legend.position="none",
    axis.text.x = element_text(angle=90, size=16, face='bold'),
    axis.ticks.x=element_blank(),
    axis.text.y = element_text(size=16, face="bold"),
    axis.title=element_text(size=30,face="bold"),
    strip.text.y = element_text(size = 16, angle = 270, face="bold"),
    strip.text.x = element_text(size = 16, angle = 0, face="bold")) + 
    scale_fill_manual(values=c("#888888", "#AEAEAE", "#CCCCCC", "#E6E6E6")) + 
  geom_text(aes(label=.id), position=position_dodge(width=0.9), vjust=4, size=10)

# Now produce the bar plot!
row.names(aucVals) <- NULL
prettyNames <- c('tfMRI 1', 'tfMRI 2', 'PCASL', 'rsfMRI', 'T1/Euler')
aucVals[,1] <- prettyNames
plotData <- as.data.frame(aucVals)
plotData$V2 <- as.numeric(as.character(plotData$V2))
plotData$V1 <- factor(plotData$V1, levels=c('T1/Euler', 'PCASL', 'tfMRI 1', 'tfMRI 2', 'rsfMRI'))
aucPlot <- ggplot(plotData, aes(x=V1, y=V2, fill=V1)) + 
                 geom_bar(stat='identity', position=position_dodge(), size=.1) + 
                 labs(title='', x='Sequence', y="AUC") +
                 theme_bw() + 
                 theme(legend.position="none",
                 axis.text.x = element_text(angle=90, size=16, face='bold'),
        	 axis.ticks.x=element_blank(),
                 axis.text.y = element_text(size=16, face="bold"),
                 axis.title=element_text(size=30,face="bold"),
                 strip.text.y = element_text(size = 16, angle = 270, face="bold"),
                 strip.text.x = element_text(size = 16, angle = 0, face="bold")) + 
                 coord_cartesian(ylim=c(.75,1)) + 
                 scale_fill_manual(values=c("#4D4D4D", "#888888", "#AEAEAE", "#CCCCCC", "#E6E6E6"))


# Now plot our data
png('fig6.png', height=6, width=18, units='in', res=300)
multiplot(meanMotionPlot, survivorPlot, aucPlot, cols=3)
dev.off()
