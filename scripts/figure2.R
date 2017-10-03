# This script will be used to produce figure 2 for Rosen et al. T1QA
# The figure will give details of mean rating bin count, weighted kappa, and polychoric cor vals
# The end product will be a 3x3 image, rows will be data set, cols will be concordance metrics

## Load data
source("/home/adrose/T1QA/scripts/galton/loadMgiData.R")
source("/home/adrose/T1QA/scripts/galton/loadGo1Data.R")

## Now load any library(s) we need
install_load("corrplot", "caret", "ggplot2", "irr", "grid", "polycor")
set.seed(16)

## Now prepare our data
raw.lme.data <- merge(isolatedVars, manualQAData2, by = "bblid")
raw.lme.data$averageRating.x <- as.numeric(as.character(raw.lme.data$averageRating.x))
raw.lme.data$averageRating.x[raw.lme.data$averageRating.x > 1] <- 1
# Load our folds
load('/home/adrose/RosenT1QA/data/folds/foldsToUse.RData')
index <- unlist(folds[1])
trainingData <- raw.lme.data[index, ]
validationData <- raw.lme.data[-index, ]
all.train.data <- merge(trainingData, manualQAData, by = "bblid")
all.valid.data <- merge(validationData, manualQAData, by = "bblid")

# Ensure all subjects have Euler numbers
all.train.data <- all.train.data[complete.cases(all.train.data$mean_euler),]
all.valid.data <- all.valid.data[complete.cases(all.valid.data$mean_euler),]
all.mgi.data <- all.mgi.data[complete.cases(all.mgi.data$mean_euler),]


## Now lets produce our ICC values 
attach(all.train.data)
trainValue <- matrix(NA, nrow=3, ncol=3)

# Start with jason's values
trainValue[1,1] <- kappa2(cbind(ratingJB.x, ratingJB.x),weight='squared')$value
trainValue[2,1] <- kappa2(cbind(ratingJB.x, ratingKS.x),weight='squared')$value
trainValue[3,1] <- kappa2(cbind(ratingJB.x, ratingLV.x),weight='squared')$value

# Now do Kevin's column
trainValue[1,2] <- kappa2(cbind(ratingKS.x, ratingJB.x),weight='squared')$value
trainValue[2,2] <- kappa2(cbind(ratingKS.x, ratingKS.x),weight='squared')$value
trainValue[3,2] <- kappa2(cbind(ratingKS.x, ratingLV.x),weight='squared')$value

# And now Prayosha's
trainValue[1,3] <- kappa2(cbind(ratingLV.x, ratingJB.x),weight='squared')$value
trainValue[2,3] <- kappa2(cbind(ratingLV.x, ratingKS.x),weight='squared')$value
trainValue[3,3] <- kappa2(cbind(ratingLV.x, ratingLV.x),weight='squared')$value

# Now fix the column and row names
colnames(trainValue) <- c('Rater 1', 'Rater 2', 'Rater 3')
rownames(trainValue) <- c('Rater 1', 'Rater 2', 'Rater 3')

# All done with the train data
detach(all.train.data)
trainValueDone <- trainValue

## Now do the validation data
attach(all.valid.data)
trainValue <- matrix(NA, nrow=3, ncol=3)

# Start with jason's values
trainValue[1,1] <- kappa2(cbind(ratingJB.x, ratingJB.x),weight='squared')$value
trainValue[2,1] <- kappa2(cbind(ratingJB.x, ratingKS.x),weight='squared')$value
trainValue[3,1] <- kappa2(cbind(ratingJB.x, ratingLV.x),weight='squared')$value

# Now do Kevin's column
trainValue[1,2] <- kappa2(cbind(ratingKS.x, ratingJB.x),weight='squared')$value
trainValue[2,2] <- kappa2(cbind(ratingKS.x, ratingKS.x),weight='squared')$value
trainValue[3,2] <- kappa2(cbind(ratingKS.x, ratingLV.x),weight='squared')$value

# And now Prayosha's
trainValue[1,3] <- kappa2(cbind(ratingLV.x, ratingJB.x),weight='squared')$value
trainValue[2,3] <- kappa2(cbind(ratingLV.x, ratingKS.x),weight='squared')$value
trainValue[3,3] <- kappa2(cbind(ratingLV.x, ratingLV.x),weight='squared')$value

# Now fix the column and row names
colnames(trainValue) <- c('Rater 1', 'Rater 2', 'Rater 3')
rownames(trainValue) <- c('Rater 1', 'Rater 2', 'Rater 3')

# All done with the validation data set
detach(all.valid.data)
validValueDone <- trainValue

# Now do the MGI data which we will call our validation data
attach(all.mgi.data)
trainValue <- matrix(NA, nrow=3, ncol=3)
# Start with jason's values
trainValue[1,1] <- kappa2(cbind(ratingJB, ratingJB),weight='squared')$value
trainValue[2,1] <- kappa2(cbind(ratingJB, ratingKS),weight='squared')$value
trainValue[3,1] <- kappa2(cbind(ratingJB, ratingLV),weight='squared')$value

# Now do Kevin's column
trainValue[1,2] <- kappa2(cbind(ratingKS, ratingJB),weight='squared')$value
trainValue[2,2] <- kappa2(cbind(ratingKS, ratingKS),weight='squared')$value
trainValue[3,2] <- kappa2(cbind(ratingKS, ratingLV),weight='squared')$value

# And now Prayosha's
trainValue[1,3] <- kappa2(cbind(ratingLV, ratingJB),weight='squared')$value
trainValue[2,3] <- kappa2(cbind(ratingLV, ratingKS),weight='squared')$value
trainValue[3,3] <- kappa2(cbind(ratingLV, ratingLV),weight='squared')$value

# Now fix the column and row names
colnames(trainValue) <- c('Rater 1', 'Rater 2', 'Rater 3')
rownames(trainValue) <- c('Rater 1', 'Rater 2', 'Rater 3')

# All done with MGI
detach(all.mgi.data)
mgiValid <- trainValue

# Now create our cor matrices plots
trainData <- melt(trainValueDone)
trainCor <- ggplot(data = trainData, aes(x=Var1, y=Var2, fill=value)) +
	geom_tile() +
	scale_fill_gradient(low = "black", high = "white",
		limit = c(.5,1), space = "Lab") +
	geom_text(aes(Var2, Var1, label = round(value, digits=2)), color = "white", size = 16) +
	theme(
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		panel.grid.major = element_blank(),
		panel.border = element_blank(),
		panel.background = element_blank(),
		axis.ticks = element_blank(),
		legend.justification = c(1, 0),
		legend.position = c(0.6, 0.7),
		legend.direction = "horizontal",
		plot.title=element_text(size=40),
		axis.text.x=element_text(size=30, angle=90, color='white'),
		axis.text.y=element_text(size=30)) +
	guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
		title.position = "top", title.hjust = 0.5)) +
	theme(legend.position="none") + 
		ggtitle(expression(paste("Weighted-", kappa))) + 
	coord_equal()

validData <- melt(validValueDone)
validCor <- ggplot(data = validData, aes(x=Var1, y=Var2, fill=value)) +
	geom_tile() +
	scale_fill_gradient(low = "black", high = "white",
		limit = c(.5,1), space = "Lab") +
	geom_text(aes(Var2, Var1, label = round(value, digits=2)), color = "white", size = 16) +
	theme(
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		panel.grid.major = element_blank(),
		panel.border = element_blank(),
		panel.background = element_blank(),
		axis.ticks = element_blank(),
		legend.justification = c(1, 0),
		legend.position = c(0.6, 0.7),
		legend.direction = "horizontal",
		plot.title=element_text(size=40, color="white"),
		axis.text.x=element_text(size=30, angle=90, color='white'),
		axis.text.y=element_text(size=30, color='black')) +
	theme(legend.position="none") + 
	ggtitle(expression(paste("Testing: Internal", kappa))) + 
	coord_equal()

mgiData <- melt(mgiValid)
mgiCor <- ggplot(data = mgiData, aes(x=Var1, y=Var2, fill=value)) +
	geom_tile() +
	scale_fill_gradient(low = "black", high = "white",
		limit = c(.5,1), space = "Lab") +
	geom_text(aes(Var2, Var1, label = round(value, digits=2)), color = "white", size = 16) +
	theme(
		axis.title.x = element_blank(),
		panel.grid.major = element_blank(),
		panel.border = element_blank(),
		panel.background = element_blank(),
		axis.ticks = element_blank(),
		legend.justification = c(1, 0),
		legend.position = c(0.6, 0.7),
		legend.direction = "horizontal",
		plot.title=element_text(size=40, color="white"),
		axis.text.x=element_text(size=30, angle=90, color='black'),
		axis.text.y=element_text(size=30, color='black')) +
		theme(legend.position="none") +
	ggtitle(expression(paste("Testing: External weighted-", kappa))) +
	coord_equal()

# Now we need to create our bar graphs
dataQaDfTrain <- as.data.frame(table(round(all.train.data$rawAverageRating.x, digits=2)))
dataQaDfTrain <- cbind(dataQaDfTrain, c(0,0,0,1,1.33,1.67,2))
colnames(dataQaDfTrain)[3] <- 'color'
trainBG <- ggplot(dataQaDfTrain, aes(x=Var1, y=Freq, fill=factor(color))) +
	geom_bar(stat='identity') +
	labs(title='Distribution of Manual Quality Ratings', x='', y='Training') +
	geom_text(data=dataQaDfTrain,aes(x=Var1,y=Freq,label=Freq),vjust=0, size=12) +
	theme_bw() +
	theme(legend.position="none",
		axis.title.x = element_blank(),
		axis.ticks = element_blank(),
		axis.text.x=element_text(size=30, color='white'),
		axis.text.y=element_text(size=30),
		axis.title.y=element_text(size=40, angle=90),
		axis.title.x=element_text(size=30),
		plot.title=element_text(size=40)) +
		scale_y_continuous(limits=c(0,1000), breaks=round(seq(0, 1000, 200), digits=2)) + 
	scale_fill_grey()

# Now do the validation data
dataQaDfValid <- as.data.frame(table(round(all.valid.data$rawAverageRating.x, digits=2)))
dataQaDfValid <- cbind(dataQaDfValid, c(0,0,0,1,1.33,1.67,2))
colnames(dataQaDfValid)[3] <- 'color'
validBG <- ggplot(dataQaDfValid, aes(x=Var1, y=Freq, fill=factor(color))) +
	geom_bar(stat='identity') +
	labs(title='', x='', y='Testing: Internal') +
	geom_text(data=dataQaDfValid,aes(x=Var1,y=Freq,label=Freq),vjust=0, size=12) +
	theme_bw() +
	theme(legend.position="none",
		axis.text.x=element_blank(),
		axis.text.y=element_text(size=30),
		axis.title.x=element_text(size=30),
		axis.title.y=element_text(size=40, angle=90),
		plot.title=element_text(size=40),
		axis.title.x=element_blank(), 
		axis.ticks.x=element_blank()) +
	scale_y_continuous(limits=c(0,500), breaks=round(seq(0, 500, 100), digits=2)) + 
	scale_fill_grey()

# Now do MGI data
# Now do the validation data
dataQaDfMGI <- as.data.frame(table(round(all.mgi.data$rawAverageRating, digits=2)))
dataQaDfMGI <- cbind(dataQaDfMGI, c(0,0,0,1,1.33,1.67,2))
colnames(dataQaDfMGI)[3] <- 'color'
mgiBG <- ggplot(dataQaDfMGI, aes(x=Var1, y=Freq, fill=factor(color))) +
geom_bar(stat='identity') +
labs(title='', x='Mean Manual Quality Rating', y='Testing: External') +
geom_text(data=dataQaDfMGI,aes(x=Var1,y=Freq,label=Freq),vjust=0, size=12) +
theme_bw() +
theme(legend.position="none",
axis.text.x=element_text(size=30),
axis.text.y=element_text(size=30),
axis.title.x=element_text(size=30),
axis.title.y=element_text(size=40, angle=90),
plot.title=element_text(size=40)) +
scale_y_continuous(limits=c(0,250), breaks=round(seq(0, 250, 50), digits=2)) + scale_fill_grey()

# Now do the polychoric cor's down here
attach(all.train.data)
trainValue <- matrix(NA, nrow=3, ncol=3)
# Start with jason's values
trainValue[1,1] <- polychor(x=ratingJB.x, y=ratingJB.x, ML=TRUE, maxcor=1)
trainValue[2,1] <- polychor(x=ratingJB.x, y=ratingKS.x, ML=TRUE, maxcor=1)
trainValue[3,1] <- polychor(x=ratingJB.x, y=ratingLV.x, ML=TRUE, maxcor=1)

# Now do Kevin's column
trainValue[1,2] <- polychor(x=ratingKS.x, y=ratingJB.x, ML=TRUE, maxcor=1)
trainValue[2,2] <- polychor(x=ratingKS.x, y=ratingKS.x, ML=TRUE, maxcor=1)
trainValue[3,2] <- polychor(x=ratingKS.x, y=ratingLV.x, ML=TRUE, maxcor=1)

# And now prayosha's
trainValue[1,3] <- polychor(x=ratingLV.x, y=ratingJB.x, ML=TRUE, maxcor=1)
trainValue[2,3] <- polychor(x=ratingLV.x, y=ratingKS.x, ML=TRUE, maxcor=1)
trainValue[3,3] <- polychor(x=ratingLV.x, y=ratingLV.x, ML=TRUE, maxcor=1)

# Now fix the column and row names
diag(trainValue) <- 1
colnames(trainValue) <- c('Rater 1', 'Rater 2', 'Rater 3')
rownames(trainValue) <- c('Rater 1', 'Rater 2', 'Rater 3')

# All done with the train data
detach(all.train.data)
trainValueDone <- trainValue

## Now do the validation data
attach(all.valid.data)
trainValue <- matrix(NA, nrow=3, ncol=3)
# Start with jason's values
trainValue[1,1] <- polychor(x=ratingJB.x, y=ratingJB.x, ML=TRUE, maxcor=1)
trainValue[2,1] <- polychor(x=ratingJB.x, y=ratingKS.x, ML=TRUE, maxcor=1)
trainValue[3,1] <- polychor(x=ratingJB.x, y=ratingLV.x, ML=TRUE, maxcor=1)

# Now do Kevin's column
trainValue[1,2] <- polychor(x=ratingKS.x, y=ratingJB.x, ML=TRUE, maxcor=1)
trainValue[2,2] <- polychor(x=ratingKS.x, y=ratingKS.x, ML=TRUE, maxcor=1)
trainValue[3,2] <- polychor(x=ratingKS.x, y=ratingLV.x, ML=TRUE, maxcor=1)

# And now prayosha's
trainValue[1,3] <- polychor(x=ratingLV.x, y=ratingJB.x, ML=TRUE, maxcor=1)
trainValue[2,3] <- polychor(x=ratingLV.x, y=ratingKS.x, ML=TRUE, maxcor=1)
trainValue[3,3] <- polychor(x=ratingLV.x, y=ratingLV.x, ML=TRUE, maxcor=1)

# Now fix the column and row names
diag(trainValue) <- 1
colnames(trainValue) <- c('Rater 1', 'Rater 2', 'Rater 3')
rownames(trainValue) <- c('Rater 1', 'Rater 2', 'Rater 3')

# All done with the validation data set
detach(all.valid.data)
validValueDone <- trainValue

## Now do the MGI data
attach(all.mgi.data)
trainValue <- matrix(NA, nrow=3, ncol=3)
# Start with jason's values
trainValue[1,1] <- polychor(x=ratingJB, y=ratingJB, ML=TRUE, maxcor=1)
trainValue[2,1] <- polychor(x=ratingJB, y=ratingKS, ML=TRUE, maxcor=1)
trainValue[3,1] <- polychor(x=ratingJB, y=ratingLV, ML=TRUE, maxcor=1)

# Now do Kevin's column
trainValue[1,2] <- polychor(x=ratingKS, y=ratingJB, ML=TRUE, maxcor=1)
trainValue[2,2] <- polychor(x=ratingKS, y=ratingKS, ML=TRUE, maxcor=1)
trainValue[3,2] <- polychor(x=ratingKS, y=ratingLV, ML=TRUE, maxcor=1)

# And now prayosha's
trainValue[1,3] <- polychor(x=ratingLV, y=ratingJB, ML=TRUE, maxcor=1)
trainValue[2,3] <- polychor(x=ratingLV, y=ratingKS, ML=TRUE, maxcor=1)
trainValue[3,3] <- polychor(x=ratingLV, y=ratingLV, ML=TRUE, maxcor=1)

# Now fix the column and row names
diag(trainValue) <- 1
colnames(trainValue) <- c('Rater 1', 'Rater 2', 'Rater 3')
rownames(trainValue) <- c('Rater 1', 'Rater 2', 'Rater 3')

# All done with the validation data set
detach(all.mgi.data)
mgiValueDone <- trainValue

# Now create our cor matrices plots
my_y_title <- expression(paste("Polychoric ", italic("r")))
trainDataPoly <- melt(trainValueDone)
trainCorPoly <- ggplot(data = trainDataPoly, aes(x=Var1, y=Var2, fill=value)) +
	geom_tile() +
		scale_fill_gradient(low = "black", high = "white",
		limit = c(.5,1), space = "Lab") +
	geom_text(aes(Var2, Var1, label = round(value, digits=2)), color = "black", size = 16) +
	theme(
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		panel.grid.major = element_blank(),
		panel.border = element_blank(),
		panel.background = element_blank(),
		axis.ticks = element_blank(),
		legend.justification = c(1, 0),
		legend.position = c(0.6, 0.7),
		legend.direction = "horizontal",
		plot.title=element_text(size=40),
		axis.text.x=element_text(size=30, angle=90, color='white'),
		axis.text.y=element_text(size=30, color='white')) +
	guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
	title.position = "top", title.hjust = 0.5)) +
	theme(legend.position="none") + 
	labs(title=my_y_title) + 
	coord_equal()

validDataPoly <- melt(validValueDone)
validCorPoly <- ggplot(data = validDataPoly, aes(x=Var1, y=Var2, fill=value)) +
	geom_tile() +
	scale_fill_gradient(low = "black", high = "white",
		limit = c(.5,1), space = "Lab") +
	geom_text(aes(Var2, Var1, label = round(value, digits=2)), color = "black", size = 16) +
	theme(
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		panel.grid.major = element_blank(),
		panel.border = element_blank(),
		panel.background = element_blank(),
		axis.ticks = element_blank(),
		legend.justification = c(1, 0),
		legend.position = c(0.6, 0.7),
		legend.direction = "horizontal",
		plot.title=element_text(size=40, color='white'),
		axis.text.x=element_text(size=30, angle=90, color='white'),
		axis.text.y=element_text(size=30, color='white')) +
	theme(legend.position="none") + 
	labs(title=my_y_title) + 
	coord_equal()

mgiDataPoly <- melt(mgiValueDone)
mgiCorPoly <- ggplot(data = mgiDataPoly, aes(x=Var1, y=Var2, fill=value)) +
	geom_tile() +
	scale_fill_gradient(low = "black", high = "white",
		limit = c(.5,1), space = "Lab") +
	geom_text(aes(Var2, Var1, label = round(value, digits=2)), color = "black", size = 16) +
	theme(
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		panel.grid.major = element_blank(),
		panel.border = element_blank(),
		panel.background = element_blank(),
		axis.ticks = element_blank(),
		legend.justification = c(1, 0),
		legend.position = c(0.6, 0.7),
		legend.direction = "horizontal",
		plot.title=element_text(size=40, color='white'),
		axis.text.x=element_text(size=30, angle=90),
		axis.text.y=element_text(size=30, color='white')) +
	theme(legend.position="none") +
	labs(title=my_y_title) +
	coord_equal() 

# Now create our plot
png('fig2.png', height=24, width=30, units='in', res=300)
multiplot(trainBG,  validBG, mgiBG, trainCor, validCor, mgiCor, trainCorPoly, validCorPoly, mgiCorPoly, cols=3)
dev.off()

# Now run a repeated effects anova down here
raw.lme.data <- merge(isolatedVars, manualQAData2, by='bblid')
raw.lme.data$averageRating.x <- as.numeric(as.character(raw.lme.data$averageRating.x))
raw.lme.data$averageRating.x[raw.lme.data$averageRating.x>1] <- 1
load('/home/adrose/qapQA/data/foldsToUse.RData')
raw.lme.data[,2:32] <- scale(raw.lme.data[,2:32], center=T, scale=T)
index <- unlist(folds[1])
trainingData <- raw.lme.data[index,]
validationData <- raw.lme.data[-index,] 

names(all.mgi.data)[1] <- 'bblid'
raw.lme.data.mgi <- melt(all.mgi.data, id.vars=names(raw.lme.data)[1:32], measure.vars=names(raw.lme.data)[35:37])
raw.lme.data.train <- melt(trainingData, id.vars=names(raw.lme.data)[1:32], measure.vars=names(raw.lme.data)[35:37])
raw.lme.data.valid <- melt(validationData, id.vars=names(raw.lme.data)[1:32], measure.vars=names(raw.lme.data)[35:37])

# Now run the repated effects AOV's
aov.train <- aov(value ~ variable, data=raw.lme.data.train)
aov.valid <- aov(value ~ variable, data=raw.lme.data.valid)
aov.mgi <- aov(value ~ variable, data=raw.lme.data.mgi)

