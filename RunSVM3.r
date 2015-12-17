##Bahiyah Nor @ May 23 2015
# For identification and classification of ion channels
# Build, validate and test SVM using e1071
# Data consists of 475 predictors
# Three types: 1. frequency of amino acids [,2:21]
#              2. dipeptide composition [,22:421]
#              3. hydrophilicity, hydrophobicity, side chain mass [,422,476]
# Construct 4 models: 1 for each type of predictors, 1 using all of them
# Using 5-fold cross validation to evaluate
# Based on prediction and training accuracy
##-----------------------------------------------------------
#Install "e1071" package for SVM
if(!is.element("e1071",installed.packages()[,1]))
  install.packages("e1071")

##-----------------------------------------------------------
#Define initial data
attach(TrainingPseAAC)
genes <- factor(TrainingPseAAC$Gene)

##-----------------------------------------------------------
#Build my SVM
#Using "e1071" svm
library("e1071")
#First model - Amino acid frequency
aaFrequency <- TrainingPseAAC[,2:21]
aa <- as.matrix(aaFrequency)
aa.data <- cbind(Gene,aaFrequency)
rm(aaFrequency)
aa.tune <- tune.svm(genes~aa,kernel="radial",cross=5)
aa.best <- aa.tune$best.model
rm(aa.tune)
attach(aa.data)
aa.model <- svm(Gene~.,data=aa.data,type="C-classification",kernel="radial",gamma=aa.best$gamma,cost=aa.best$cost,cross=5, probability=TRUE)
rm(aa)
#Second model - dipeptide
diPeptide <- TrainingPseAAC[,22:421]
diPep <- as.matrix(diPeptide)
diPep.data <- cbind(Gene,diPeptide)
rm(diPeptide)
diPep.tune <- tune.svm(genes~diPep,kernel="radial",cross=5)
diPep.best <- diPep.tune$best.model
rm(diPep.tune)
attach(diPep.data)
diPep.model <- svm(Gene~.,data=diPep.data,type="C-classification",kernel="radial",gamma=diPep.best$gamma,cost=diPep.best$cost,cross=5, probability=TRUE)
rm(diPep)
#Third model - hydrophilicity, hydrophobicity, side chain mass
chemistry <- TrainingPseAAC[,422:476]
chem <- as.matrix(chemistry)
chem.data <- cbind(Gene,chemistry)
rm(chemistry)
chem.tune <- tune.svm(genes~chem,kernel="radial",cross=5)
chem.best <- chem.tune$best.model
rm(chem.tune)
attach(chem.data)
chem.model <- svm(Gene~.,data=chem.data,type="C-classification",kernel="radial",gamma=chem.best$gamma,cost=chem.best$cost,cross=5, probability=TRUE)
rm(chem)
#Fourth model - all the predictors
traindata <- TrainingPseAAC[,2:ncol(TrainingPseAAC)]
train <- as.matrix(traindata)
rm(traindata)
train.tune <- tune.svm(genes~train,kernel="radial", cross=5)
train.best <- train.tune$best.model
attach(TrainingPseAAC)
train.model <- svm(Gene~.,data=TrainingPseAAC,type="C-classification",kernel="radial",gamma=train.best$gamma,cost=train.best$cost,cross=5, probability=TRUE)
rm(train)
#Fifth Model - Chou's
aa <- TrainingPseAAC[,2:21]
chem <- TrainingPseAAC[,422:476]
chou.data <- cbind(Gene,aa,chem)
rm(aa,chem)
chou <- as.matrix(chou.data[,2:ncol(chou.data)])
chou.tune <- tune.svm(genes~chou,kernel="radial",cross=5)
chou.best <- chou.tune$best.model
attach(chou.data)
chou.model <- svm(Gene~.,data=chou.data,type="C-classification",kernel="radial",gamma=chou.best$gamma,cost=chou.best$cost,cross=5,probability=TRUE)
rm(chou)

##-----------------------------------------------------------
#Test the SVM
testdata <- BLASTp.e45.PseAAC0522[,2:476]
predicted.aa <- predict(aa.model,testdata,probability=TRUE)
predicted.chem <- predict(chem.model,testdata,probability=TRUE)
predicted.chou <- predict(chou.model,testdata,probability=TRUE)
predicted.diPep <- predict(diPep.model,testdata,probability=TRUE)
predicted.train <- predict(train.model,testdata,probability=TRUE)

##-----------------------------------------------------------
#Get the test probabilities for 'train' model
pvals <- matrix(0,nrow=768,ncol=5)
roc.data <- data.frame(pvals)
rm(pvals)
colnames(roc.data) <- c("SequenceID","Pvals","ActualClass","Predicted","True")
c <- 1

for(i in 1:768)
{
  actual.label <- testdata.labels$Label[i]
  #if(actual.label != 'NIOC')
  #{
    roc.data$SequenceID[c] <- as.character(testdata.labels$Gene.True[i])
    roc.data$Pvals[c] <- attr(predicted.train,"probabilities")[i,predicted.train[i]]
    roc.data$ActualClass[c] <- actual.label
    roc.data$Predicted[c] <- as.character(predicted.train[i])
    if(roc.data$ActualClass[c] == roc.data$Predicted[c])
      roc.data$True[c] <- 1
    c <- c+1
  #}
}
sorted.rocdata <- roc.data[order(-roc.data$Pvals),]

#Get the test probabilities for 'dipep' model
pvals <- matrix(0,nrow=389,ncol=4)
dipep.roc <- data.frame(pvals)
rm(pvals)
colnames(dipep.roc) <- c("Pvals","ActualClass","Predicted","True")
c <- 1

for(i in 1:768)
{
  actual.label <- testdata.labels$Label[i]
  if(actual.label != 'NIOC')
  {
    dipep.roc$Pvals[c] <- attr(predicted.diPep,"probabilities")[i,predicted.diPep[i]]
    dipep.roc$ActualClass[c] <- actual.label
    dipep.roc$Predicted[c] <- as.character(predicted.diPep[i])
    if(dipep.roc$ActualClass[c] == dipep.roc$Predicted[c])
      dipep.roc$True[c] <- 1
    c <- c+1
  }
}

sorted.dipeproc <- dipep.roc[order(-dipep.roc$Pvals),]

#Plot the ROC curce
if(!is.element("pROC",installed.packages()[,1]))
  install.packages("pROC")
library("pROC")
train.roc <- roc(sorted.rocdata$True ~ sorted.rocdata$Pvals, data=sorted.rocdata)
dipep.roc <- roc(sorted.dipeproc$True ~ sorted.dipeproc$Pvals, data=sorted.dipeproc)
train.auc <- plot(train.roc)$auc
dipep.auc <- plot(dipep.roc)$auc

## Get the OV prediction results
OVpvals <- matrix(0,nrow=123,ncol=3)
OV.results <-data.frame(OVpvals)
colnames(OV.results) <- c("SeqID","Classification","Pvals")
rm(OVpvals)

for(i in 1:123)
{
  OV.results$SeqID[i] <- as.character(OV.PseAAC$Gene[i])
  OV.results$Classification[i] <- as.character(OV.predicted[i])
  OV.results$Pvals[i] <- attr(OV.predicted,"probabilities")[i,OV.predicted[i]]
}

OV.sorted <- OV.results[order(-OV.results$Pvals),]
write.table(OV.results, "C:/Users/Bahiyah Nor/Documents/MSc Bioinformatics/BINF90001/OV_results.txt", sep="\t")
