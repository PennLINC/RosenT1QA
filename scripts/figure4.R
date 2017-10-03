# This script will be used to produce figure 4 for Rosen et al. T1QA
# This script is going to be used to plot the differences amongst the training and testing sets 
# mean values from the quantitative variables and their relation to the qualititative manual ratings.
# It is going to plot the 9 values from the outcome of the 1 vs 2 octavariate model
# These include:
#	qi1 wm.skewness cnr bg.kurtosis efc bg.skewness fber snr euler_number

## Load data
# Start with Go1
source('/home/adrose/T1QA/scripts/galton/loadMgiData.R')
source('/home/adrose/T1QA/scripts/galton/loadGo1Data.R')
set.seed(16)

# Modify summarySE to make it easier to work with
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)
    
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

# load library(s)
install_load('caret', 'ggplot2', 'grid', 'gridExtra')

# Now lets create our train and validation sets
raw.lme.data <- merge(isolatedVars, manualQAData2, by='bblid')
raw.lme.data$averageRating.x <- as.numeric(as.character(raw.lme.data$averageRating.x))
raw.lme.data$averageRating.x[raw.lme.data$averageRating.x>1] <- 1
raw.lme.data[,2:32] <- scale(raw.lme.data[,2:32], center=T, scale=T)
# Load our folds
load('/home/adrose/RosenT1QA/data/folds/foldsToUse.RData')
index <- unlist(folds[1])
trainingData <- raw.lme.data[index,]
validationData <- raw.lme.data[-index,] 

## Now merge our scaled data values with the original data values 
all.train.data <- merge(mergedQAP, trainingData, by='bblid')
all.valid.data <- merge(mergedQAP, validationData, by='bblid')
names(all.train.data) <- gsub(pattern='.x', x=names(all.train.data), replacement='')
names(all.valid.data) <- gsub(pattern='.x', x=names(all.valid.data), replacement='')
trainingData$mean_euler <- scale(trainingData$mean_euler)
validationData$mean_euler <- scale(validationData$mean_euler)

# Now ensure complete cases
all.train.data <- all.train.data[complete.cases(all.train.data$mean_euler),]
all.valid.data <- all.valid.data[complete.cases(all.valid.data$mean_euler),]
all.mgi.data <- all.mgi.data[complete.cases(all.mgi.data$mean_euler),]

# Now lets declare our variables of interest 
varsOfInterest <- c('bg.kurtosis', 'bg.skewness', 'cnr', 'efc', 'fber', 'qi1', 'snr', 'wm.skewness', 'mean_euler')
prettyNames <- c('BG Kurtosis', 'BG Skewness', 'CNR', 'EFC', 'FBER', 'QI1', 'SNR', 'WM Skewness', 'Mean Euler')

# Now lets create our training values
trainingValues <- NULL 
i <- 1
for(qapVal in varsOfInterest){
  valsToAppend <- summarySE(data=trainingData, measurevar=qapVal, groupvars='averageRating.y', na.rm=T)
  qapValue <- rep(prettyNames[i], nrow(valsToAppend))
  Dataset <- rep('Training', nrow(valsToAppend))
  valsToAppend <- cbind(valsToAppend, qapValue, Dataset)
  trainingValues <- rbind(trainingValues, valsToAppend)
  i <- i + 1
}
i <- 1
for(qapVal in varsOfInterest){
  valsToAppend <- summarySE(data=validationData, measurevar=qapVal, groupvars='averageRating.y', na.rm=T)
  qapValue <- rep(prettyNames[i], nrow(valsToAppend))
  Dataset <- rep('Testing: Internal', nrow(valsToAppend))
  valsToAppend <- cbind(valsToAppend, qapValue, Dataset)
  trainingValues <- rbind(trainingValues, valsToAppend)
  i <- i + 1
}
all.mgi.data$averageRating <- all.mgi.data$rawAverageRating
all.mgi.data$averageRating[all.mgi.data$averageRating < 1] <- 0.00
i <- 1
for(qapVal in varsOfInterest){
    all.mgi.data[qapVal] <- scale(unlist(all.mgi.data[qapVal]), center=T, scale=T)
    valsToAppend <- summarySE(data=all.mgi.data, measurevar=qapVal, groupvars='averageRating')
    qapValue <- rep(prettyNames[i], nrow(valsToAppend))
    Dataset <- rep('Testing: External', nrow(valsToAppend))
    valsToAppend <- cbind(valsToAppend, qapValue, Dataset)
    colnames(valsToAppend)[1] <- 'averageRating.y'
    trainingValues <- rbind(trainingValues, valsToAppend)
    i <- i + 1
}
trainingValues$averageRating.y[which(trainingValues$averageRating.y==1.670)] <- 1.667
trainingValues$averageRating.y[which(trainingValues$averageRating.y==1.330)] <- 1.333
# Now produce all of our partial corellations 
corVals <- NULL
i <- 1
all.train.data$ageSquared <- all.train.data$ageAtGo1Scan^2
residAverageRating <- lm(averageRating.y ~ ageAtGo1Scan + s + ageSquared, data=all.train.data)$residuals
for(qapVal in varsOfInterest){
  form1 <- as.formula(paste(qapVal,' ~ ageAtGo1Scan + s + ageSquared', sep=''))
  residQuant <- lm(form1, data=all.train.data)$residuals
  corVal <- cor(residAverageRating, residQuant, method='spearman')
  corPVal <- cor.test(residAverageRating, residQuant, method='spearman')$p.value
  qapValue <- prettyNames[i]
  Dataset <- 'Training'
  valsToAppend <- cbind(corVal, qapValue, Dataset, corPVal)
  corVals <- rbind(corVals, valsToAppend)
  i <- i + 1
}
i <- 1
all.valid.data$ageSquared <- all.valid.data$ageAtGo1Scan^2
residAverageRating <- lm(averageRating.y ~ ageAtGo1Scan + s + ageSquared, data=all.valid.data)$residuals
for(qapVal in varsOfInterest){
  form1 <- as.formula(paste(qapVal,' ~ ageAtGo1Scan + s + ageSquared', sep=''))
  residQuant <- lm(form1, data=all.valid.data)$residuals
  corVal <- cor(residAverageRating, residQuant, method='spearman')
  corPVal <- cor.test(residAverageRating, residQuant, method='spearman')$p.value
  qapValue <- prettyNames[i]
  Dataset <- 'Testing: Internal'
  valsToAppend <- cbind(corVal, qapValue, Dataset, corPVal)
  corVals <- rbind(corVals, valsToAppend)
  i <- i + 1
}
i <- 1
all.mgi.data$ageSquared <- all.mgi.data$age^2
residAverageRating <- lm(averageRating ~ age + Gender + ageSquared, data=all.mgi.data)$residuals
for(qapVal in varsOfInterest){
    form1 <- as.formula(paste(qapVal,' ~ age + Gender + ageSquared', sep=''))
    residQuant <- lm(form1, data=all.mgi.data)$residuals
    corVal <- cor(residAverageRating, residQuant, method='spearman')
    corPVal <- cor.test(residAverageRating, residQuant, method='spearman')$p.value
    qapValue <- prettyNames[i]
    Dataset <- 'Testing: External'
    valsToAppend <- cbind(corVal, qapValue, Dataset, corPVal)
    corVals <- rbind(corVals, valsToAppend)
    i <- i + 1
}



# Now make our plot
corVals <- corVals[order(as.numeric(corVals[,1])),]
corVals <- as.data.frame(corVals)
corVals$qapValue <- factor(corVals$qapValue, levels=unique(as.character(corVals$qapValue)))
corVals$Dataset <- factor(corVals$Datase, levels=c('Training', 'Testing: Internal', 'Testing: External'))
corPlot <- ggplot(corVals, 
                 aes(x=factor(qapValue), y=as.numeric(as.character(corVal)), fill=factor(Dataset))) + 
                 geom_bar(stat='identity', position=position_dodge(), size=.1) + 
                 labs(title='', x='Quantitative Metric', y=expression("Correlation with Manual Quality Rating (partial"~rho*")")) +
                 theme_bw() + 
                 facet_grid(Dataset ~ qapValue, space = "free", scales='free_x') + 
                 theme(legend.position="none",
                 axis.text.x = element_blank(),
        	 axis.ticks.x=element_blank(),
                 axis.text.y = element_text(size=16, face="bold"),
                 axis.title=element_text(size=30,face="bold"),
                 strip.text.y = element_text(size = 16, angle = 270, face="bold"),
                 strip.text.x = element_text(size = 16, angle = 90, face="bold")) +
		 scale_fill_manual(values=c("black", "black", "black"))

trainingValues$qapValue <- factor(trainingValues$qapValue, levels=levels(corVals$qapValue))
allPlot <- ggplot(trainingValues, 
                 aes(x=factor(averageRating.y), y=as.numeric(as.character(mean)), fill=factor(averageRating.y))) + 
                 geom_bar(stat='identity', position=position_dodge(), size=.1) + 
                 labs(title='', x='Manual Quality Rating (mean value)', y='Quality Metric (mean z-score)') +
                 geom_errorbar(aes(ymin=as.numeric(as.character(mean))-se, ymax=as.numeric(as.character(mean))+se), 
                       width = .1, position=position_dodge(.9)) + 
                 theme_bw() + 
                 facet_grid(Dataset ~ qapValue) + 
                 theme(legend.position="none",
                 axis.text.x = element_blank(),
                 axis.text.y = element_text(size=16, face="bold"),
        	 axis.ticks.x=element_blank(),
                 axis.title=element_text(size=30,face="bold"),
                 strip.text.y = element_text(size = 16, angle = 270, face="bold"),
                 strip.text.x = element_text(size = 16, angle = 90, face="bold"),
                 panel.margin = unit(2, "lines")) + scale_fill_grey()

png('fig4.png', height=12, width=20, units='in', res=300)
grid.arrange(allPlot, corPlot, ncol = 2, layout_matrix = cbind(c(1,1,1,1,1,1),c(1,1,1,1,1,1),c(2,2,2,2,2,2)))
dev.off()
