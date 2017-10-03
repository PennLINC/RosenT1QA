# Now create a function which will return the p val for the prediction of the imaging values
# based on the quality index with age and sex regressed values
returnPVal <- function(imagingVal, qualityVal, regVals, df, regressAgeBOO=TRUE, regressSexBOO=TRUE, regressTBV=FALSE){
    if(regressAgeBOO == 'TRUE'){
        imagingVal <- regressAge(df, imagingVal)
        qualityVal <- regressAge(df, qualityVal)
    }
    if(regressSexBOO == 'TRUE'){
        imagingVal <- regressSex(df, imagingVal)
        qualityVal <- regressSex(df, qualityVal)
    }
    if(regressTBV == 'TRUE'){
        regVals <- paste(regVals, 'df$mprage_antsCT_vol_TBV', sep='+')
    }
    form <- as.formula(paste('imagingVal~qualityVal', paste(regVals), sep=''))
    outputPVal <- summary(lm(form))$coefficients[2,4]
    #outputPVal <- cor.test(imagingVal, qualityVal, method='kendall')$p.value
    
    return(outputPVal)
}

# Now same function for t value
returnTVal <- function(imagingVal, qualityVal, regVals, df, regressAgeBOO=TRUE, regressSexBOO=TRUE, regressTBV=FALSE){
    if(regressAgeBOO == 'TRUE'){
        imagingVal <- regressAge(df, imagingVal)
        qualityVal <- regressAge(df, qualityVal)
    }
    if(regressSexBOO == 'TRUE'){
        imagingVal <- regressSex(df, imagingVal)
        qualityVal <- regressSex(df, qualityVal)
    }
    if(regressTBV == 'TRUE'){
        regVals <- paste(regVals, 'df$mprage_antsCT_vol_TBV', sep='+')
    }
    form <- as.formula(paste('imagingVal~qualityVal', paste(regVals), sep=''))
    outputPVal <- summary(lm(form))$coefficients[2,3]
    #outputPVal <- cor.test(imagingVal, qualityVal, method='kendall')$p.value
    
    return(outputPVal)
}

# Now create some functions which will be used to create the output colored labels for itksnap
returnHeatMapITKSnapVals <- function(inputZScores, lowColor='blue', hiColor='red'){
    # Create some functions this function will call... yeesh
    range01 <- function(x){
        # Now make sure we have some standard deviation
        # If no standard deviation return 1
        if( is.na(sd(x)) == 'TRUE'){
            output <- rep(1, length(x))
            return(output)
        }
        else if (sd(x) < 0 ){
            output <- rep(1, length(x))
            return(output)
        }
        else
        (x-min(x))/diff(range(x))
    }
    cRamp <- function(x){
        cols <- colorRamp(c(lowColor, hiColor))(range01(as.numeric(x)))
    }
    # Output values
    outputValues <- matrix(0, nrow=(length(inputZScores)+1), ncol=8)
    
    # Now cretae our rgb values
    redValues <- round(cRamp(inputZScores)[,1], digits=0)
    greenValues <- round(cRamp(inputZScores)[,2], digits=0)
    blueValues <- round(cRamp(inputZScores)[,3], digits=0)
    
    # Now create our index column
    outputValues[,1] <- seq(0, length(inputZScores))
    
    # Now put the proper values in the correct place
    outputValues[2:(length(inputZScores)+1),2] <- redValues
    outputValues[2:(length(inputZScores)+1),3] <- greenValues
    outputValues[2:(length(inputZScores)+1),4] <- blueValues
    
    # Now we need to do the Transperancy column
    outputValues[,5] <- c(0, rep(1, length(inputZScores)))
    
    # Now the visibility column
    outputValues[,6] <- c(0, rep(1, length(inputZScores)))
    
    # Now the mesh visibility
    outputValues[,7] <- c(0, rep(1, length(inputZScores)))
    
    # Now the label indicies
    labelIndexNames <- c('Clear Label', paste('Label ', inputZScores, sep=''))
    labelIndexNames <- paste('"', labelIndexNames, '"', sep='')
    outputValues[,8] <- labelIndexNames
    
    # Now return our output
    return(outputValues)
}

# Now I need to create a function which will combine my output colormaps
returnPosNegAndNeuColorScale <- function(outputZScores, colorScaleNeg=c('light blue', 'blue'), colorScalePos=c('yellow', 'red'), colorScaleNeu=c('gray'), sigThreshold=.05){
    # MAKE SURE WE ARE DEALING WITH NUMERICS!!!!
    outputZScores <- as.numeric(as.character(outputZScores))
    
    # First convert our sig threshold into a z score to find our cut off value
    cutOff <- abs(qnorm(sigThreshold))
    
    # Now we need to make our seperate our data into neutral, positive, and negative values
    # We are going to order these just so it is easier to match the labesl to the output ROI
    # when working with the ouput of this function
    negativeValues <- outputZScores[which(outputZScores < 0 & abs(outputZScores) >= cutOff)]
    negativeValues <- negativeValues[order(negativeValues)]
    positiveValues <- outputZScores[which(outputZScores >= cutOff)]
    positiveValues <- positiveValues[order(positiveValues)]
    neutralValues <- outputZScores[which(abs(outputZScores) < cutOff )]
    neutralValues <- neutralValues[order(neutralValues)]
    
    # Create our blank label row first
    values <- rep(0, 7)
    blankRow <- append(values, paste('"', 'Clear Label' ,'"', sep=''))
    
    # Now we need to create our individual color scales
    #startPoint <- NULL
    output <- blankRow
    if(length(negativeValues) > 0 ){
        negativeColors <- returnHeatMapITKSnapVals(negativeValues, lowColor=colorScaleNeg[1], hiColor=colorScaleNeg[2])[2:(length(negativeValues)+1),]
        #negIndex <- max(as.numeric(as.character(negativeColors[,1])))
        #startPoint <- cbind(startPoint, negIndex)
        output <- rbind(output, negativeColors)
    }
    if(length(neutralValues) > 0){
        neutralColors <- returnHeatMapITKSnapVals(neutralValues, lowColor=colorScaleNeu[1], hiColor=colorScaleNeu[1])[2:(length(neutralValues)+1),]
        #neuIndex <- max(as.numeric(as.character(neutralColors[,1])))
        #startPoint <- cbind(startPoint, neuIndex)
        output <- rbind(output, neutralColors)
    }
    if(length(positiveValues) > 0 ){
        positiveColors <- returnHeatMapITKSnapVals(positiveValues, lowColor=colorScalePos[1], hiColor=colorScalePos[2])[2:(length(positiveValues)+1),]
        #posIndex <- max(as.numeric(as.character(positiveColors[,1])))
        #startPoint <- cbind(startPoint, posIndex)
        output <- rbind(output, positiveColors)
    }
    # Now I need to make sure that the index column doesn't have any repeats
    # This will be done by running an an index thorugh the first column
    output[,1] <- seq(0, length(outputZScores))
    
    # Now we are all set! just need to return our output
    return(output)
}


# Now create a function which will return all of the pVals for a specific grep pattern
# which will be the prefix for an imaging value
pvalLoop <- function(grepPattern, dataFrame, TBV=FALSE, correct=TRUE){
    # First lets find the values that we need to loop through
    colVals <- grep(grepPattern, names(dataFrame))
    
    # Now lets compute our p vals
    outputPVals <- apply(dataFrame[,colVals], 2, function(x) returnPVal(x ,dataFrame$oneVsTwoOutcome, '+df$ageAtGo1Scan+df$sex', all.train.data, regressAgeBOO=FALSE, regressSexBOO=FALSE))
    if(TBV=='TRUE'){
        outputPVals <- apply(dataFrame[,colVals], 2, function(x) returnPVal(x ,dataFrame$oneVsTwoOutcome, '+df$ageAtGo1Scan+df$sex', all.train.data, regressAgeBOO=FALSE, regressSexBOO=FALSE, regressTBV=TRUE))
    }
    # Now fdr correct these suckers
    if(correct==TRUE){
      outputPVals.fdr <- p.adjust(outputPVals, method='fdr')
    }
    if(correct==FALSE){
      outputPVals.fdr <- outputPVals
    }
    # Now append the T values to the output
    outputTVals <- apply(dataFrame[,colVals], 2, function(x) returnTVal(x ,dataFrame$oneVsTwoOutcome, '+df$ageAtGo1Scan+df$sex', all.train.data, regressAgeBOO=FALSE, regressSexBOO=FALSE))
    
    # I need to add some direction to my z scores
    # These directions will come from the t values
    outputZScores <- qnorm(outputPVals.fdr, lower.tail=F) * sign(outputTVals)
    
    # Now return the output
    output <- cbind(names(outputZScores), outputZScores)
    output <- output[order(as.numeric(output[,2])),]
    rownames(output) <- NULL
    return(output)
}
