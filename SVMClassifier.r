#MuSICC.Classifier.r @ Bahiyah Nor, Aug 18th 2015
#
#Usage: Rscript MuSICC.Classifier.r

##-----------------------------------------------------------
#Install "e1071" package for SVM
if(!is.element("e1071",installed.packages()[,1]))
  install.packages("e1071")
library("e1071")

##-----------------------------------------------------------
#Load the Pre-trained SVM model
SVM.model <- readRDS("SVMModel/MuSICC.SVMrmodel.rds")

#Load the new data to predict
pseAAC <- read.table("Output/PseAAC.csv", sep="\t",header=TRUE)
data <- pseAAC[,2:476]
data.predict <- predict(SVM.model,data,probability=TRUE)

#Combine the prediction values with Sequence ID
tmp <- matrix(data=NA,nrow=nrow(data),ncol=3)
fullclass <- data.frame(tmp)
rm(tmp)
colnames(fullclass) <- c("SequenceId","SVMClassification","Probability")
for(i in 1:nrow(data))
{
  fullclass$SequenceId[i] <- as.character(pseAAC[i,1])
  fullclass$SVMClasssification[i] <- as.character(data.predict[i])
  fullclass$Probability[i] <- attr(data.predict,"probabilities")[i,data.predict[i]]
}
#Sort it according to the probability values
sorted <- fullclass[order(-fullclass$Probability),]

filename <- paste("MuSICC.Results",Sys.Date(),".txt")
write.table(sorted,filename,sep="\t")
