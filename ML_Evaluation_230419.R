#ML Model Evaluation
#JBaxter

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
groupbyvirus = function(x, name = name){
  grouped = list(list(x[[1]],x[[2]],x[[3]],x[[4]]),
                 list(x[[5]],x[[6]],x[[7]],x[[8]]),
                 list(x[[9]],x[[10]],x[[11]],x[[12]]),
                 list(x[[13]],x[[14]],x[[15]],x[[16]]),
                 list(x[[17]],x[[18]],x[[19]],x[[20]]),
                 list(x[[21]],x[[22]],x[[23]],x[[24]]),
                 list(x[[25]],x[[26]],x[[27]],x[[28]]),
                 list(x[[29]],x[[30]],x[[31]],x[[32]]),
                 list(x[[33]],x[[34]],x[[35]],x[[36]]),
                 list(x[[37]],x[[38]],x[[39]],x[[40]]),
                 list(x[[41]],x[[42]],x[[43]],x[[44]]),
                 list(x[[45]],x[[46]],x[[47]],x[[48]]))
  for(i in 1:length(grouped)){
    names(grouped[[i]]) = name
  }
  return(grouped)
}

prep_multiroc = function(modelpredictions , trueclasses){
  classOHE = one_hot(data.table(trueclasses))
  colnames(classOHE) = c(unlist(lapply(levels(trueclasses) , function(x) paste(x ,'true',sep = '_'))))
  classpredict = data.frame(modelpredictions)
  colnames(classpredict) =  c(unlist(lapply(levels(trueclasses) , function(x) paste(x ,'pred', 'model', sep = '_'))))
  data = cbind.data.frame(classOHE , classpredict)
  return(data)
}

mc_roc = function(x){
  t = multi_roc(x)
  roc_res_df =plot_roc_data(t)
  groups = split.data.frame(roc_res_df , roc_res_df$Group)
  micro = groups$Micro
  micro = as.data.frame(cbind((1-micro[,1]) , micro[,2]))
  colnames(micro) = c( 'FPR','Sensitivity')
  return(micro)
}


eval.MCC = function(x,y){
  MCC = lapply(x, function(i) mltools::mcc(i , y))
}

eval.ACC = function(x,y){
  ACC = lapply(x, function(i) MLmetrics::Accuracy(i , y))
} 

eval.AUC = function(x,y){
  ACC = lapply(x, function(i) MLmetrics::AUC(i , y))
} 

metric.dfs = function(x , metric = metric){
  bymodel = reshape2::melt(x)
  met.df = cbind.data.frame(bymodel , 'Test Model')
  colnames(met.df) = c(paste(metric), 'Algorithm' , 'Virus' , 'Group')
  return(met.df )
}

modelmeans = function(x){
  splitbymodel = split.data.frame(x , x$Algorithm)
  modelmeans = lapply(splitbymodel , function(x) mean(x$MCC))
  return(modelmeans)
}

#Import test labels
test_Data = lapply(filenames[grep('testData' , filenames)] , function(x) read.csv(x , stringsAsFactors = FALSE)[,-1])
test_labels = lapply(test_Data , function(x) x[,1])

#Import predictions as matrices
filenames = list.files(path = '.')
filenames = filenames[grep('.csv' , filenames)]
prob.predict.init = lapply(filenames[grep('testprobs' , filenames)] , function(x) read.csv(x, stringsAsFactors = FALSE)[,-1])
class.predict.init =  lapply(filenames[grep('testclass' , filenames)] , function(x) read.csv(x, stringsAsFactors = FALSE)[,-1])
nullprob.predict.init =  lapply(filenames[grep('nullprobs' , filenames)] , function(x) read.csv(x, stringsAsFactors = FALSE)[,-1])
nullclass.predict.init =  lapply(filenames[grep('nullclass' , filenames)] , function(x) read.csv(x, stringsAsFactors = FALSE)[,-1])

names = unlist(lapply(filenames , function(x) {
  split = unlist(strsplit(x , '[.,_]+'))
  name = split[1]
  return(name)
}))

names = unique(names)
names = names[-c(5,12,13,14,18,17)]

models = c('KNN' , 'NNET' ,'RF' , 'SVM')

#Group by Model
prob.predict = groupbyvirus(prob.predict.init, name = models)
names(prob.predict) = names

class.predict = groupbyvirus(class.predict.init, name = models)
names(prob.predict) = names

nullprob.predict = groupbyvirus(nullprob.predict.init , name = models)
names(prob.predict) = names

nullclass.predict = groupbyvirus(nullclass.predict.init , name = models)
names(nullclass.predict) = names


#Matthews Correlation Coefficient
test.MCC = mapply(eval.MCC , x = class.predict , y = test_labels, SIMPLIFY = FALSE)
names(test.MCC) = names
null.MCC = mapply(eval.MCC , x = nullclass.predict , y = test_labels, SIMPLIFY = FALSE)
names(null.MCC) = names

testmcc.df = metric.dfs(test.MCC , 'MCC')
testbymodel = reshape2::melt(test.MCC)
testmcc.df = cbind.data.frame(testbymodel , 'Test Model')
colnames(testmcc.df) = c('MCC' , 'Algorithm' , 'Virus' , 'Group')
nullbymodel = reshape2::melt(null.MCC)
nullmcc.df = cbind.data.frame(nullbymodel  , 'Null Model')
colnames(nullmcc.df) = c('MCC' , 'Algorithm' , 'Virus' , 'Group')

ggpubr::ggbarplot(data = testmcc.df , x = 'Virus' , y = 'MCC' ,  add = c("mean_se") , color = 'Group' , palette = 'npg') + rotate_x_text(45) 


#Accuracy
test.ACC = mapply(eval.ACC , x = class.predict , y = test_labels, SIMPLIFY = FALSE)
names(test.ACC) = names
null.ACC = mapply(eval.ACC , x = nullclass.predict , y = test_labels, SIMPLIFY = FALSE)
names(null.ACC) = names

testACC.df = metric.dfs(test.ACC , 'Accuracy')

#ROC
multirocd = list()
rocplots = list()

for (i in 1:12){
  data_roc = lapply(prob.predict[[i]], function(x) prep_multiroc(x , as.factor(test_labels[[i]])))
  micro_roc = lapply(data_roc , mc_roc)
  names(micro_roc) = c('RF' , 'KNN' , 'NNet' , 'SVM')
  multiroc_data = bind_rows(micro_roc, .id = "column_label")
  colnames(multiroc_data) = c('Model','FPR' , 'TPR' )
  multirocd[[i]] = multiroc_data
  rocplots[[i]] = ggline(multiroc_data , 'FPR', 'TPR',plot_type = 'l' , group = 'Model' , color = 'Model' , size = 1,palette = 'npg' , numeric.x.axis = TRUE )+geom_abline(intercept = 0, slope = 1)
}

ggarrange(rocplots[[1]], rocplots[[2]] , 
          rocplots[[3]],rocplots[[4]],
          rocplots[[5]],rocplots[[6]],
          rocplots[[7]],rocplots[[8]],
          rocplots[[9]],rocplots[[10]],
          rocplots[[11]],rocplots[[12]],
          ncol = 3 , nrow = 4,labels='AUTO' , common.legend = TRUE)

split.rates = lapply(multirocd , function(x) split.data.frame(x , x$Model))

auc = lapply(split.rates , function(x) {
  auc = list()
  for (i in 1:4){
    model = x[[i]]
    auc[[i]] = cal_auc(model[,2] , model[,3])
  }
  return(auc)
})

names(auc) = names
auc = lapply(auc , function(x) {
  names(x) = c('KNN' , 'NNET' , 'RF' , 'SVM')
  return(x)})

testAUC.df = reshape2::melt(auc)
colnames(testAUC.df) = c ('AUC', 'Algorithm' , 'Virus' )

#Evaluate suitability of data for statistical test (para vs non-para)
#hist(combinedmetrics$Value)

#Fitted vs Residuals
#res.aov = aov(Value ~ Metric, data = combinedmetrics)
#plot(res.aov, 1)
#plot(res.aov, 2)

#Compare Models and Controls Panel
p1 = ggboxplot(testmcc.df , x = 'Algorithm' , y = 'MCC' , color = 'black' , fill = 'Algorithm' , palette = 'npg', label = 'Virus', repel = TRUE, label.select = list(criteria = "`y` < 0.2"))+
  stat_compare_means(method = "kruskal.test" , label.y = 2)  +
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "KNN" , inherit.aes = TRUE , label.y = 1.6)

combinedmcc.df = rbind.data.frame(testmcc.df , nullmcc.df)

my_comparisons = compare_means(MCC~Group, data = combinedmcc.df, 
                               group.by = "Algorithm")

p2 = ggpubr::ggbarplot(data = combinedmcc.df , x = 'Algorithm' , y = 'MCC' ,  add = c("mean_ci"),
                       color = 'Group'  , palette = 'npg',
                       position = position_dodge(0.9)) +  stat_compare_means(aes(group = Group),label = "p.signif")

ggarrange(p1 , p2 ,ncol = 1 , nrow = 2, labels = 'AUTO')


#Compare Metrics Panel
metcomp1 = ggpubr::ggbarplot(data = testmcc.df , x = 'Virus' , y = 'MCC' ,  add = c("mean_se") , remove = '') + theme(axis.text.x = element_blank())

metcomp2 = ggpubr::ggbarplot(data = testACC.df, x = 'Virus' , y = 'Accuracy' ,  add = c("mean_se")) + theme(axis.text.x = element_blank())

metcomp3 = ggpubr::ggbarplot(data = testAUC.df, x = 'Virus' , y = 'AUC' ,  add = c("mean_se"))+ theme(axis.text.x = element_blank())


mcc2 = cbind.data.frame(testmcc.df , 'MCC')
colnames(mcc2) = colnames(mcc2) = c('Value' , 'Algorithm' , 'Virus' , 'Group' , 'Metric')

auc2 = cbind.data.frame(testAUC.df , 'AUC')
colnames(auc2) = colnames(auc2) = c('Value' , 'Algorithm' , 'Virus' , 'Metric')

acc2 = cbind.data.frame(testACC.df , 'Accuracy')
colnames(acc2) = colnames(acc2) = c('Value' , 'Algorithm' , 'Virus' , 'Group' , 'Metric')

combinedmetrics = rbind.data.frame(mcc2[,c(1,5)] , auc2[,c(1,4)] , acc2[,c(1,5)])

means = ggpubr::ggbarplot(data = combinedmetrics, x = 'Metric' , y = 'Value' ,  add = c("mean_se"), color = 'black' , 
                          fill = 'Metric' , palette = 'npg' , legend = 'none' , xlab = FALSE)+ stat_compare_means(method = "kruskal.test" , label.y = 1.1 , label.x =0.75)  +
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "MCC" , inherit.aes = TRUE , label.y = 0.9)

ggarrange(metcomp1 , metcomp2 , metcomp3 , means ,labels = 'AUTO')





###END###
