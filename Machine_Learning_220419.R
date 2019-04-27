ML
#310319
#05.04.19 - set functions to iterate over files
#220419 - evaluation functions moved to new script, model and prediction saving introduced
#J Baxter

library(nnet)
library(randomForest)

library(kernlab)
library(caret)
library(MLmetrics)
library(ModelMetrics)
library(dplyr)
library(multiROC)
library(mltools)
library(data.table)
library(ggpubr)

#set current working directory
setwd('~/Desktop/fastas_for_raxml/test/')

#define functions

formatfeatures = function(x){
  feature_set = list()
  for (i in 1:length(x)){
    set = x[[i]]
    Host = make.names(set$Host)
    set$Host = NULL
    feature_set[[i]] = cbind.data.frame(Host , set)
    colnames(feature_set[[i]]) = gsub('X' , '' , colnames(feature_set[[i]]))
  }
  return(feature_set)
}

names = unlist(lapply(filenames , function(x) {
  split = unlist(strsplit(x , '[/,_]+'))
  name = split[1]
  return(name)
}))

mixclasses = function(x){
  hostdf= data.frame(x$Host)
  x$Host = NULL
  redistributedhosts = hostdf[sample(1:nrow(hostdf)), ]
  nulldata = cbind(redistributedhosts , x)
  colnames(nulldata) =c('Host' , colnames(x))
  return(nulldata)
}

setFolds = function(x){
  require(dplyr)
  index = 1:nrow(x)
  indexed = cbind.data.frame(x , index)
  splitdf = split.data.frame(indexed , indexed$Host)
  minhost = min(unlist(lapply(splitdf , nrow)))
  folds = list()
  for (i in 1:10){
    t = bind_rows(lapply(splitdf , function(i) sample_n(i, (minhost), replace=FALSE)))
    folds[[i]] = t$index
  }
  return(folds)}


setCV = function(x){
  cv.cntrl = trainControl(method = "repeatedcv", number = 10, repeats=10, index = x , classProbs = TRUE , 
                          summaryFunction = multiClassSummary , search = 'random')
}

trainModels = function(data, cv){
  rf = caret::train(Host ~ ., data = data, 
                    method = "rf", 
                    trControl = cv, 
                    verbose = TRUE, 
                    metric = "AUC" , 
                    importance = TRUE)
  
  svm =caret::train(Host ~ ., data = data, 
                    method = "svmPoly", 
                    trControl = cv, 
                    verbose = FALSE, 
                    metric = "AUC")
  
  nnet= caret::train(Host ~ ., data = data, 
                     method = "nnet", 
                     trControl = cv, 
                     verbose = FALSE,
                     MaxNWts = 25000,
                     metric = "AUC")
  
  knn= caret::train(Host ~ ., data = data, 
                    method = "knn", 
                    trControl = cv,
                    metric = "AUC")
  
  trained.models = list(rf , knn , nnet , svm)
  names(trained.models) = c('RF' , 'KNN' , 'NNet' , 'SVM')
  return(trained.models)
}

predictProbs = function(x,y){
  probs = lapply(x , function(i) predict(i , y , type = 'prob'))
  names(probs) = c('RF.probs' , 'GBM.probs' , 'NNet.probs' , 'SVM.probs')
  return(probs)
}

predictClass = function(x,y){
  class = lapply(x , function(i) predict(i , y , type = 'raw'))
  names(class) = c('RF' , 'KNN' , 'NNET' , 'SVM')
  return(class)
}


##
#import datasets
filenames = list.files(path = '.')
filenames = filenames[grep('featureset' , filenames)]
importfiles = lapply(filenames[-c(5,12,13,14,18,17)], read.csv)

virus.names = names[-c(5,12,13,14,18,17)]

names(importfiles) = virus.names 

feature_set = formatfeatures(importfiles)

#split test and training sets
ind = lapply(feature_set , function(x) createDataPartition(x[,1], times = 1 , p= 0.7 , list = FALSE))

trainData = list()
for (i in 1:length(feature_set)){
  trainData[[i]] = feature_set[[i]][ind[[i]],]
}


testData = list()
for (i in 1:length(feature_set)){
  testData[[i]] = feature_set[[i]][-ind[[i]],]
}

#define null training set
trainNull = lapply(trainData , mixclasses)

#set cv folds
Kfolds = lapply(trainData , setFolds)
Nullfolds = lapply(trainNull , setFolds)

cv.cntrl = lapply(Kfolds , setCV)
null.cntrl = lapply(Nullfolds , setCV)

#train models
library(doSNOW)
library(parallel)

start = Sys.time()
print(start)
cl = makeCluster(6)
clusterExport(cl = cl, c('trainModels' , 'trainData' , 'cv.cntrl'), envir = .GlobalEnv)

trainedModels = parallel::clusterMap(cl = cl , trainModels , trainData, cv.cntrl , SIMPLIFY = FALSE)

stopCluster(cl)
remove(cl)
duration = Sys.time() - start
print(duration)

#define test data
test_nolabels= lapply(testData , function(x) x[,-1])
test_labels = lapply(testData , function(x) x[,1])

#predict
prob.predict = mapply(predictProbs , x = trainedModels , y = test_nolabels , SIMPLIFY = FALSE)
class.predict = mapply(predictClass , x = trainedModels , y = test_nolabels , SIMPLIFY = FALSE)

#run Null Models
start = Sys.time()
print(start)
cl = makeCluster(6)
clusterExport(cl = cl, c('trainModels' , 'trainNull' , 'null.cntrl'), envir = .GlobalEnv)

nullModels = parallel::clusterMap(cl = cl , trainModels , trainNull, null.cntrl , SIMPLIFY = FALSE)

stopCluster(cl)
remove(cl)
duration = Sys.time() - start
print(duration)

#predict
nullprob.predict = mapply(predictProbs , x = nullModels , y = test_nolabels , SIMPLIFY = FALSE)
nullclass.predict = mapply(predictClass , x = nullModels , y = test_nolabels , SIMPLIFY = FALSE)


#Set filenames
models = c('RF' , 'SVM' , 'NNET' , 'KNN')
modelnames = lapply(virus.names , function(x) {
  name = list()
  for (i in 1:length(models)){
    name[[i]] = paste(x , models[i] , sep = '.')
  }
  return(unlist(name))
})

modelnames = unlist(modelnames)

#Save Models
trainedModels.f = purrr::flatten(trainedModels)
lapply(1:length(trainedModels.f), function(i) saveRDS(trainedModels.f[[i]], 
       file = paste0(modelnames[i], '_test_' , Sys.Date(), ".RDS")))

nullModels.f = purrr::flatten(nullModels)
lapply(1:length(nullModels.f), function(i) saveRDS(nullModels.f[[i]], 
       file = paste0(modelnames[i], '_null_' , Sys.Date(), ".RDS")))

#Save Predictions
prob.predict.f = purrr::flatten(prob.predict)
lapply(1:length(prob.predict.f), function(i) write.csv(prob.predict.f[[i]],
                                                     file = paste0(modelnames[i],'_testprobs_' , Sys.Date(), ".csv"), row.names = TRUE))

class.predict.f = purrr::flatten(class.predict)
lapply(1:length(class.predict.f), function(i) write.csv(class.predict.f[[i]], 
                                                   file = paste0(modelnames[i], '_testclass_' , Sys.Date(), ".csv"), row.names = TRUE))

nullprob.predict.f = purrr::flatten(nullprob.predict)
lapply(1:length(nullprob.predict.f), function(i) write.csv(nullprob.predict.f[[i]], 
                                                     file = paste0(modelnames[i], '_nullprobs_' , Sys.Date(), ".csv"), row.names = TRUE))

nullclass.predict.f = purrr::flatten(nullclass.predict)
lapply(1:length(nullclass.predict.f), function(i) write.csv(nullclass.predict.f[[i]], 
                                                      file = paste0(modelnames[i], '_nullclass_' , Sys.Date(), ".csv"), row.names = TRUE))
###END###
