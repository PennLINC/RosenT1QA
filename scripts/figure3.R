# This script will be used to produce figure 3 for Rosen et al. T1QA
# The figure will give details of sex and age differences amongst the mean manual rating
# The end product will be a 2x3 image, cols will be data set, rows will be demographic trends

## Load data
# Start with Go1
source('/home/adrose/T1QA/scripts/galton/loadMgiData.R')
source('/home/adrose/T1QA/scripts/galton/loadGo1Data.R')
set.seed(16)

## Now declare any functions 
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
          mean = median   (xx[[col]], na.rm=na.rm),
          sd   = IQR     (xx[[col]], na.rm=na.rm)/2
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

# load library(s)
install_load('caret', 'ggplot2', 'lme4', 'car', 'visreg', 'scales', 'MASS')

# Now lets create our train and validation sets
raw.lme.data <- merge(isolatedVars, manualQAData2, by='bblid')
raw.lme.data$averageRating.x <- as.numeric(as.character(raw.lme.data$averageRating.x))
raw.lme.data$averageRating.x[raw.lme.data$averageRating.x>1] <- 1
raw.lme.data <- raw.lme.data[complete.cases(raw.lme.data$mean_euler),]
all.mgi.data <- all.mgi.data[complete.cases(all.mgi.data$mean_euler),]
# Load our folds
load('/home/adrose/RosenT1QA/data/folds/foldsToUse.RData')
index <- unlist(folds[1])
trainingData <- raw.lme.data[index,]
validationData <- raw.lme.data[-index,] 
manualQAData$age <- (manualQAData$ageAtGo1Scan / 12)

# Now prep our individual data sets
all.train.data <- merge(trainingData, manualQAData, by='bblid')
all.valid.data <- merge(validationData, manualQAData, by='bblid')

# Now age reg the ratings
all.train.data$averageRatingAR <- scale(residuals(lm(averageRating ~ ageAtGo1Scan, data = all.train.data)))
all.valid.data$averageRatingAR <- scale(residuals(lm(averageRating ~ ageAtGo1Scan, data = all.valid.data)))
all.mgi.data$averageRatingAR <- scale(residuals(lm(averageRating ~ age, data = all.mgi.data)))

# Now sex reg the ratings
all.train.data$averageRatingSR <- residuals(lm(averageRating ~ sex, data = all.train.data))
all.valid.data$averageRatingSR <- residuals(lm(averageRating ~ sex, data = all.valid.data))
all.mgi.data$averageRatingSR <- residuals(lm(averageRating ~ Gender, data = all.mgi.data))

# Now age reg 
all.train.data$ageSR <- residuals(lm(age ~ sex, data=all.train.data))
all.valid.data$ageSR <- residuals(lm(age ~ sex, data=all.valid.data))
all.mgi.data$ageSR <- residuals(lm(age ~ Gender, data=all.mgi.data))

# Now prepare our values
bg2.vals.train <- summarySE(data=all.train.data, measurevar='averageRatingAR', groupvars='sex')
bg2.vals.train$Dataset <- rep('Training', nrow(bg2.vals.train))
bg2.vals.valid <- summarySE(data=all.valid.data, measurevar='averageRatingAR', groupvars='sex')
bg2.vals.valid$Dataset <- rep('Testing', nrow(bg2.vals.valid))
bg2.vals.mgi <- summarySE(data=all.mgi.data, measurevar='averageRatingAR', groupvars='Gender')
colnames(bg2.vals.mgi)[1] <- 'sex'
bg2.vals.mgi$Dataset <- 'Validation'
bg2.vals <- rbind(bg2.vals.train, bg2.vals.valid, bg2.vals.mgi)
bg2.vals$Dataset <- factor(bg2.vals$Dataset)
bg2.vals$sex <- c('Male', 'Female', 'Male', 'Female', 'Male', 'Female')

# Now lets plot our values
# Grab a p value from a mann-whitney
pValue <- wilcox.test(all.train.data$averageRatingAR ~ all.train.data$sex)
# Find the min value to plot 
minVal <- round(min(bg2.vals$averageRatingAR), digits=2)-.7
maxVal <- round(max(bg2.vals$averageRatingAR), digits=2)+1
bg1 <- ggplot(bg2.vals[which(bg2.vals$Dataset=='Training'),], aes(x=factor(sex), y=as.numeric(as.character(averageRatingAR)), group=Dataset)) + 
                geom_bar(stat='identity', position=position_dodge(), width=.5) + 
                labs(title='Training', x='Sex', y='Manual Quality Rating (z-score)') +
                theme_bw() + 
                coord_cartesian(ylim=c(minVal,maxVal)) + 
                       geom_errorbar(aes(ymin=as.numeric(as.character(averageRatingAR))-se, 
                                         ymax=as.numeric(as.character(averageRatingAR))+se), 
                width = .1, position=position_dodge(.9)) +
                #facet_grid(Dataset ~ .) +
                theme(legend.position="none",
                axis.text=element_text(size=20),
                axis.title=element_text(size=30),
                strip.text.y = element_text(size = 16, angle = 270, face="bold"),
                title=element_text(size=30)) + 
		geom_path(aes(x=factor(sex), y=c(maxVal-.2,maxVal-.2))) +
		geom_path(aes(x=factor(sex)[1], y=c(maxVal-.4, maxVal-.2))) +
		geom_path(aes(x=factor(sex)[2], y=c(maxVal-.4, maxVal-.2))) +
		geom_text(aes(x=factor(sex)[1], y=.19), label='',angle=90, size=10) +
		scale_y_continuous(limits=c(minVal, maxVal), 
                           breaks=round(seq(minVal, maxVal, .5), digits=1), oob=rescale_none) + 
		annotate("text", x=c(Inf), y=c(Inf), label="p > 0.1", hjust=c(2.1), vjust=c(1.1), size=8, parse=T)
		
pValue <- wilcox.test(all.valid.data$averageRatingAR ~ all.valid.data$sex)
bg2 <- ggplot(bg2.vals[which(bg2.vals$Dataset=='Testing'),], aes(x=factor(sex), y=as.numeric(as.character(averageRatingAR)), group=Dataset)) +
                geom_bar(stat='identity', position=position_dodge(), width=.5) + 
                labs(title='Testing: Internal', x='Sex', y='Manual Quality Rating (mean value)') +
                theme_bw() + 
                       geom_errorbar(aes(ymin=as.numeric(as.character(averageRatingAR))-se, 
                                         ymax=as.numeric(as.character(averageRatingAR))+se), 
                width = .1, position=position_dodge(.9)) +
                #facet_grid(Dataset ~ .) +
                theme(legend.position="none",
                axis.text.y=element_text(size=20, color='white'),
                axis.title.y=element_text(size=30, color='white'),
                axis.text=element_text(size=20),
                axis.title=element_text(size=30),
		axis.ticks.y=element_blank(),
                strip.text.y = element_text(size = 16, angle = 270, face="bold"),
                title=element_text(size=30)) + 
		geom_path(aes(x=factor(sex), y=c(maxVal-.2,maxVal-.2))) +
		geom_path(aes(x=factor(sex)[1], y=c(maxVal-.4, maxVal-.2))) +
		geom_path(aes(x=factor(sex)[2], y=c(maxVal-.4, maxVal-.2))) +
		geom_text(aes(x=factor(sex)[1], y=.19), label='',angle=90, size=10) +
		scale_y_continuous(limits=c(minVal, maxVal), 
                           breaks=round(seq(minVal, maxVal, .5), digits=1), oob=rescale_none) + 
		annotate("text", x=c(Inf), y=c(Inf), label="p > 0.1", hjust=c(2.1), vjust=c(1.1), size=8, parse=T)
		
pValue <- wilcox.test(all.mgi.data$averageRatingAR~ all.mgi.data$Gender)
bg3 <- ggplot(bg2.vals[which(bg2.vals$Dataset=='Validation'),], aes(x=factor(sex), y=as.numeric(as.character(averageRatingAR)), group=Dataset)) +
		geom_bar(stat='identity', position=position_dodge(), width=.5) +
		labs(title='Testing: External', x='Sex', y='Manual Quality Rating (mean value)') +
		theme_bw() +
		geom_errorbar(aes(ymin=as.numeric(as.character(averageRatingAR))-se,
		ymax=as.numeric(as.character(averageRatingAR))+se),
		width = .1, position=position_dodge(.9)) +
		#facet_grid(Dataset ~ .) +
		theme(legend.position="none",
		axis.text.y=element_text(size=20, color='white'),
		axis.title.y=element_text(size=30, color='white'),
		axis.text=element_text(size=20),
		axis.title=element_text(size=30),
		axis.ticks.y=element_blank(),
		strip.text.y = element_text(size = 16, angle = 270, face="bold"),
		title=element_text(size=30)) +
		geom_path(aes(x=factor(sex), y=c(maxVal-.2,maxVal-.2))) +
		geom_path(aes(x=factor(sex)[1], y=c(maxVal-.4, maxVal-.2))) +
		geom_path(aes(x=factor(sex)[2], y=c(maxVal-.4, maxVal-.2))) +
		geom_text(aes(x=factor(sex)[1], y=.19), label='',angle=90, size=10) +
		scale_y_continuous(limits=c(minVal, maxVal), 
                           breaks=round(seq(minVal, maxVal, .5), digits=1), oob=rescale_none) + 
		annotate("text", x=c(Inf), y=c(Inf), label="p > 0.1", hjust=c(2.1), vjust=c(1.1), size=8, parse=T)

# Now build our models to show general age trends
corVal <- cor(all.train.data$averageRatingSR, all.train.data$ageSR, method='spearman')
corSig <- cor.test(all.train.data$averageRatingSR, all.train.data$ageSR, method='spearman')$p.value
corText1 <- expression(~rho == .14)
corText2 <- paste("p < 0.001")
mod1 <- ggplot(all.train.data, aes(y=scale(averageRatingSR), x=age)) +
   geom_smooth(method=lm, color='black') +
   theme_bw() +
   coord_cartesian(xlim=c(8,22), ylim=c(-.4,.6)) +
   labs(title='', y='Manual Quality Rating (z-score)', x='Age (years)') +
   theme(
    axis.text=element_text(size=20),
    axis.title=element_text(size=30)) +
   scale_y_continuous(breaks=c(-.4,-.2,0,.2,.4,.6)) +
   annotate("text", x=c(Inf, Inf), y=c(-Inf, -Inf), label=c(as.character(corText2), as.character(corText1)), hjust=c(1, 1), vjust=c(-.5, -2.5), size=8, parse=T)

corVal <- cor(all.valid.data$averageRatingSR, all.valid.data$age, method='spearman')
corSig <- cor.test(all.valid.data$averageRatingSR, all.valid.data$age, method='spearman')$p.value
corText1 <- expression(~rho == paste(0.12))
corText2 <- paste("p < 0.01")
mod2 <- ggplot(all.valid.data, aes(y=scale(averageRatingSR), x=age)) +
   geom_smooth(method=lm, color='black') +
   theme_bw() +
   coord_cartesian(xlim=c(8,22), ylim=c(-.4,.6)) +
   labs(title='', y='Mean Manual Quality Rating', x='Age (years)') +
   theme(
    axis.text=element_text(size=20),
    axis.title=element_text(size=30),
    axis.title.y=element_text(size=20, color='white'),
    axis.text.y=element_text(size=30, color='white'),
    axis.ticks.y=element_blank()) +
    scale_y_continuous(breaks=c(-.4,-.2,0,.2,.4,.6)) +
   annotate("text", x=c(Inf, Inf), y=c(-Inf, -Inf), label=c(as.character(corText2), as.character(corText1)), hjust=c(1, 1), vjust=c(-.5, -2.5), size=8, parse=T)

corVal <- cor(all.mgi.data$averageRatingSR, all.mgi.data$age, method='spearman')
corSig <- cor.test(all.mgi.data$averageRatingSR, all.mgi.data$age, method='spearman')$p.value
corText1 <-expression(~rho == paste(-0.15))
corText2 <- paste("p < 0.05")
mod3 <- ggplot(all.mgi.data, aes(y=scale(averageRatingSR), x=age)) +
	geom_smooth(method=lm, color='black') +
	theme_bw() +
	coord_cartesian(xlim=c(20,80), ylim=c(-.4,.6)) +
	labs(title='', y='Manual Quality Rating (mean value)', x='Age (years)') +
	theme(
	axis.text=element_text(size=20),
	axis.title=element_text(size=30),
	axis.title.y=element_text(size=20, color='white'),
	axis.text.y=element_text(size=30, color='white'),
	axis.ticks.y=element_blank()) +
	scale_y_continuous(breaks=c(0,.33,.66,1,1.33,1.66,2)) +
	scale_x_continuous(breaks=c(20,40,60,80)) +
	annotate("text", x=c(Inf, Inf), y=c(-Inf, -Inf), label=c(as.character(corText2), as.character(corText1)), hjust=c(1, 1), vjust=c(-.5, -2.5), size=8, parse=T)


png('fig3.png', width=24, height=16, units='in', res=300)
multiplot(bg1, mod1, bg2, mod2, bg3, mod3, cols=3)
dev.off()

# Now look at beta weights in these demographics in lm's
# Now build a lm model in the training data 
m1 <- lm(rawAverageRating.y ~ age + sex, data = all.train.data)
sigValsTrain <- summary(m1)

# Now do the same for the internal testing data
m1 <- lm(averageRating ~ age + sex, data=all.valid.data)
sigValsTest <- Anova(m1)

# Now do the external testing dataset
m1 <- lm(rawAverageRating ~ age + Gender, data=all.mgi.data)
sigValsValid <- Anova(m1)
