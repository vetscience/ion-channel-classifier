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
#AJS: I have changed the read.table command so that the first column is read as the rownames: row.names=1
pseAAC <- read.table("Output/PseAAC.csv", sep="\t",header=TRUE, row.names=1)
#AJS: Here you then don't need to skip the first column
data <- pseAAC
data.predict <- predict(SVM.model,data,probability=TRUE)

#Combine the prediction values with Sequence ID
tmp <- matrix(data=NA,nrow=nrow(data),ncol=3)
fullclass <- data.frame(tmp)
rm(tmp)
colnames(fullclass) <- c("SequenceId","SVMClassification","Probability")
for(i in 1:nrow(data))
{
#AJS: Here I have changed as.character(pseAAC[i,1]) to row.names(pseAAC[i,]), we now just use the row.names command
#now that the IDs are defined as row names
  fullclass$SequenceId[i] <- row.names(pseAAC[i,])
#There was a typo here: fullclass$SVMClasssification with three "s"
  fullclass$SVMClassification[i] <- as.character(data.predict[i])
  fullclass$Probability[i] <- attr(data.predict,"probabilities")[i,data.predict[i]]
}
#Sort it according to the probability values
sorted <- fullclass[order(-fullclass$Probability),]

#AJS: Here I have changed the separator to "" so that there is no "_" between the date and the ".txt"
filename <- paste("MuSICC.Results_",Sys.Date(),".txt",sep="")
#AJS: Writing the table, we need to make sure that the row names (which will be 1,2,3 etc. if not defined otherwise) are not included
write.table(sorted,filename,sep="\t", row.names=F)
