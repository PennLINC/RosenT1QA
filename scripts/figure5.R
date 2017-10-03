# This script will be used to produce figure 5 for Rosen et al. T1QA
# It will explore the predictive ability from each of the quantitative metrics 
# ability to bin 0 vs !0 images. Tjis will be done individually in each 
# dataset

## Load the data
source('/home/adrose/T1QA/scripts/galton/loadMgiData.R')
source('/home/adrose/T1QA/scripts/galton/loadGo1Data.R')
detachAllPackages()
set.seed(16)

# Declare functions 
rocdata <- function(grp, pred){
    # Produces x and y co-ordinates for ROC curve plot
    # Arguments: grp - labels classifying subject status
    #            pred - values of each observation
    # Output: List with 2 components:
    #         roc = data.frame with x and y co-ordinates of plot
    #         stats = data.frame containing: area under ROC curve, p value, upper and lower 95% confidence interval
    
    grp <- as.factor(grp)
    if (length(pred) != length(grp)) {
        stop("The number of classifiers must match the number of data points")
    }
    
    if (length(levels(grp)) != 2) {
        stop("There must only be 2 values for the classifier")
    }
    
    cut <- unique(pred)
    tp <- sapply(cut, function(x) length(which(pred > x & grp == levels(grp)[2])))
    fn <- sapply(cut, function(x) length(which(pred < x & grp == levels(grp)[2])))
    fp <- sapply(cut, function(x) length(which(pred > x & grp == levels(grp)[1])))
    tn <- sapply(cut, function(x) length(which(pred < x & grp == levels(grp)[1])))
    tpr <- tp / (tp + fn)
    fpr <- fp / (fp + tn)
    roc = data.frame(x = fpr, y = tpr)
    roc <- roc[order(roc$x, roc$y),]
    
    i <- 2:nrow(roc)
    auc <- (roc$x[i] - roc$x[i - 1]) %*% (roc$y[i] + roc$y[i - 1])/2
    
    pos <- pred[grp == levels(grp)[2]]
    neg <- pred[grp == levels(grp)[1]]
    q1 <- auc/(2-auc)
    q2 <- (2*auc^2)/(1+auc)
    se.auc <- sqrt(((auc * (1 - auc)) + ((length(pos) -1)*(q1 - auc^2)) + ((length(neg) -1)*(q2 - auc^2)))/(length(pos)*length(neg)))
    ci.upper <- auc + (se.auc * 0.96)
    ci.lower <- auc - (se.auc * 0.96)
    
    se.auc.null <- sqrt((1 + length(pos) + length(neg))/(12*length(pos)*length(neg)))
    z <- (auc - 0.5)/se.auc.null
    p <- 2*pnorm(-abs(z))
    
    stats <- data.frame (auc = auc,
    p.value = p,
    ci.upper = ci.upper,
    ci.lower = ci.lower
    )
    
    return (list(roc = roc, stats = stats))
}

rocplot.single <- function(grp, pred, title = "ROC Plot", p.value = FALSE){
    require(ggplot2)
    plotdata <- rocdata(grp, pred)

    p <- ggplot(plotdata$roc, aes(x = x, y = y)) +
    geom_line(aes(colour = ""), size=3) +
    geom_abline (intercept = 0, slope = 1) +
    theme_bw() +
    scale_x_continuous("1-Specificity") +
    scale_y_continuous("Sensitivity") +
    scale_colour_manual(values = "#000000") +
    ggtitle(title) +
    theme_bw() +
    theme(legend.position=c(1,0)) +
    theme(legend.justification=c(1,0)) +
    theme(legend.title=element_blank(),
    text = element_text(size=30)) + theme(legend.position="none")
    
    return(p)
}

## Load Library(s)
install_load('pROC', 'ggplot2', 'caret', 'lme4', 'grid', 'gridExtra')

# Now split data into train and traning 
raw.lme.data <- merge(isolatedVars, manualQAData2, by='bblid')
raw.lme.data$averageRating.x <- as.numeric(as.character(raw.lme.data$averageRating.x))
raw.lme.data$averageRating.x[raw.lme.data$averageRating.x>1] <- 1
# Load our folds
load('/home/adrose/RosenT1QA/data/folds/foldsToUse.RData')
index <- unlist(folds[1])
trainingData <- raw.lme.data[index,]
trainingData <- trainingData[complete.cases(trainingData),]
validationData <- raw.lme.data[-index,]
validationData <- validationData[complete.cases(validationData),]

# Now get the measure vars
measureVars <- names(raw.lme.data)[1:33]
# Now get the id.vars
idVars <- names(raw.lme.data)[35:37]

raw.lme.data <- melt(trainingData, id.vars=measureVars, measure.vars=idVars)
raw.lme.data$value[raw.lme.data$value > 1] <- 1
raw.lme.data.test <- melt(validationData, id.vars=measureVars, measure.vars=idVars)
raw.lme.data.test$value[raw.lme.data.test$value > 1] <- 1

# Now run through each variable of interest and build an ROC curve for it
# first limit it to QAP variables of interest
qapValNamesUse <- qapValNames[-c(1:3, 5:7, 10:12, 14:19, 22:26,33:34)]
qapValNamesUse <- c('bg.kurtosis', 'bg.skewness', 'cnr', 'efc', 'fber', 'qi1', 'snr', 'wm.skewness', 'mean_euler')
qapValNamesUse <- c('qi1', 'efc', 'wm.skewness', 'snr', 'cnr', 'fber', 'bg.skewness', 'bg.kurtosis', 'mean_euler')
aucVals <- NULL
# Now create a series of loops which will train in each of the data sets and validate in the other remaining
# The output will be the output AUC value
# The first for loop will train in the test and test and validate in the remaining

# First prepare our temporary data sets to predict in
testTmpSet <- validationData
validTmpSet <- all.mgi.data[c(qapValNames, 'averageRating')]
validTmpSet$averageRating[validTmpSet$averageRating>=1] <- 1
for(qapVal in qapValNamesUse){
  # Now get our auc vals
  form <- as.formula(paste("averageRating.x ~ ", paste(qapVal)))
  trainAUC <- auc(roc(form, data=trainingData))
  # Now get our test value auc
  testAUC <- auc(roc(form, data=validationData))
  form <- as.formula(paste("averageRating ~ ", paste(qapVal)))  
  validAUC <- auc(roc(form, data=all.mgi.data))
  # Now prepare the output
  outputData <- rbind(cbind(trainAUC, 'Training', qapVal), cbind(testAUC, 'Testing: Internal', qapVal), cbind(validAUC, 'Testing: External', qapVal))
  aucVals <- rbind(aucVals, outputData)  
}
aucValsAll <- cbind(aucVals, rep('Training', 27))

# Now create our data frame to plot
aucVals <- as.data.frame(aucValsAll)
aucVals$trainAUC <- as.numeric(as.character(aucVals$trainAUC))
aucVals$V2 <- factor(aucVals$V2, levels=c('Training', 'Testing: Internal', 'Testing: External'))
aucVals$V4 <- factor(aucVals$V4, levels=c('Training', 'Testing: Internal', 'Testing: External'))
aucVals$BG <- 0
aucVals$BG[which(aucVals$V2==aucVals$V4)] <- 1
aucVals$prettyQap <- rep(rep(c('QI1', 'EFC', 'WM Skewness', 'SNR', 'CNR', 'FBER', 'BG Skewness', 'BG Kurtosis', 'Mean Euler'), each=3), 1)
aucVals$prettyQap <- factor(aucVals$prettyQap, levels=c('QI1', 'EFC', 'WM Skewness', 'SNR', 'CNR', 'FBER', 'BG Skewness', 'BG Kurtosis', 'Mean Euler'))
aucValPlot <- ggplot(aucVals, aes(x=prettyQap, y=trainAUC)) +
  geom_bar(stat="identity", width=0.4, position=position_dodge(width=0.5)) +
  coord_cartesian(ylim=c(.5,1)) +
  facet_grid(. ~ V2, scales="free", space="free_x") +
  theme(axis.text.x = element_text(angle=90,hjust=1, size=20),
    axis.title.x = element_text(size=30),
    axis.title.y = element_text(size=30),
    text = element_text(size=30),
    panel.margin = unit(1, "lines")) +
  ggtitle("") +
  xlab("") +
  ylab("AUC")

# Now create the ROC curves
install_load('pROC', 'ggplot2', 'caret', 'lme4', 'grid', 'gridExtra')
raw.lme.data <- merge(isolatedVars, manualQAData2, by='bblid')
raw.lme.data$averageRating.x <- as.numeric(as.character(raw.lme.data$averageRating.x))
raw.lme.data$averageRating.x[raw.lme.data$averageRating.x>1] <- 1
load('/home/adrose/RosenT1QA/data/folds/foldsToUse.RData')
index <- unlist(folds[1])

# Now produce our data sets
raw.lme.data.test <- raw.lme.data[-index,]
raw.lme.data <- raw.lme.data[index,]

raw.lme.data <- raw.lme.data[complete.cases(raw.lme.data$zeroVsNotZero),]
roc.tmp <- roc(averageRating.x ~ mean_euler, data=trainingData)
trainText2 <- paste("AUC =  ", round(auc(roc.tmp), digits=2), sep='')
trainText <- c(trainText2)
trainZeroPlot <- rocplot.single(pred=trainingData$mean_euler, grp=trainingData$averageRating.x, title="")
trainZeroPlot <- trainZeroPlot + annotate("text", x=c(Inf), y=c(-Inf), label=trainText, vjust=c(-3.4), hjust="inward", size=8) + theme(axis.text.x=element_text(color='black'), axis.title.x=element_text(color='black'))
threshold <- coords(roc.tmp, 'best')[1]

# Now produce the testing ROC plot using the training data set
roc.tmp <- roc(averageRating.x ~ mean_euler, data=validationData)
trainText2 <- paste("AUC =  ", round(auc(roc.tmp), digits=2), sep='')
trainText <- c(trainText2)
testZeroPlot <- rocplot.single(pred=validationData$mean_euler, grp=validationData$averageRating.x, title="")
testZeroPlot <- testZeroPlot  + annotate("text", x=c(Inf), y=c(-Inf), label=trainText, vjust=c(-3.4), hjust="inward", size=8) + 
  theme(axis.text.x=element_text(color='black'), axis.title.x=element_text(color='black'), axis.title.y=element_text(color='white'), axis.text.y=element_text(color='white'), axis.ticks.y=element_blank())

# Now do the external testing data set 
all.mgi.data <- all.mgi.data[complete.cases(all.mgi.data$mean_euler),]
all.mgi.data$averageRating.x <- 1
all.mgi.data$averageRating.x[which(all.mgi.data$averageRating == 0)] <- 0
roc.tmp <- roc(averageRating.x ~ mean_euler, data=all.mgi.data)
trainText2 <- paste("AUC =  ", round(auc(roc.tmp), digits=2), sep='')
trainText <- c(trainText2)
validZeroPlot <- rocplot.single(pred=all.mgi.data$mean_euler, grp=all.mgi.data$averageRating.x, title="")
validZeroPlot <- validZeroPlot + annotate("text", x=c(Inf), y=c(-Inf), label=trainText, vjust=c(-3.4), hjust="inward", size=8) + theme(axis.text.x=element_text(color='black'), axis.title.x=element_text(color='black'), axis.title.y=element_text(color='white'), axis.text.y=element_text(color='white'), axis.ticks.y=element_blank())

png('fig5.png', width=20, height=16, units='in', res=300)
grid.arrange(aucValPlot, trainZeroPlot, testZeroPlot, validZeroPlot, ncol = 3, layout_matrix = rbind(c(1, 1, 1), c(2, 3, 4)))
dev.off()
