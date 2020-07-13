library(tidyverse)
library(ADNIMERGE) # load phenotypes from ADNI data
library(caret)
library(plyr) # FOR multiple merging

## load datasets from ADNIMERGE

uaspmvbm <- (ADNIMERGE::uaspmvbm) #UA (Gene Alexander) MRI SPM voxel based morphometry (VBM) analysis
ucsdvol <- (ADNIMERGE::ucsdvol) # Anders Dale Lab (UCSD) - Derived Volumes
ucsffsx51final <- (ADNIMERGE::ucsffsx51final)
upenn_roi_mars <- (ADNIMERGE::upenn_roi_mars)
picslashs <- (ADNIMERGE::picslashs)

## LOAD PHENOTYPES

label <- read.csv(file = "labels.csv")
label <-label[!duplicated(label[,c("RID")]),]
label <- label[,c("RID" , "PTID", "DX_bl")]
label$DX_bl <- str_replace_all(label$DX_bl, "LMCI", "MCI")
label$DX_bl <- str_replace_all(label$DX_bl, "SMC", "MCI")
label$DX_bl <- str_replace_all(label$DX_bl, "EMCI", "MCI")
label$DX_bl <- str_replace_all(label$DX_bl, "AD", "Dementia")
colnames(label) <- c("RID" ,"PTID", "Pheno")

# remove duplicates
uaspmvbm <- uaspmvbm[!duplicated(uaspmvbm[,c("RID")]),]
ucsdvol <- ucsdvol[!duplicated(ucsdvol[,c("RID")]),]
ucsffsx51final <- ucsffsx51final[!duplicated(ucsffsx51final[,c("RID")]),]
upenn_roi_mars <- upenn_roi_mars[!duplicated(upenn_roi_mars[,c("RID")]),]
picslashs <- picslashs[!duplicated(picslashs[,c("RID")]),]

## match with labeled file for names

uaspmvbm <- uaspmvbm[uaspmvbm$RID %in% label$RID,]
ucsdvol <- ucsdvol[ucsdvol$RID %in% label$RID,]
ucsffsx51final <- ucsffsx51final[ucsffsx51final$RID %in% label$RID,]
upenn_roi_mars <- upenn_roi_mars[upenn_roi_mars$RID %in% label$RID,]
picslashs <- picslashs[picslashs$RID %in% label$RID,]

uaspmvbm <- merge(uaspmvbm, label, by="RID")
ucsdvol <- merge(ucsdvol, label, by="RID")
ucsffsx51final <- merge(ucsffsx51final, label, by="RID")
upenn_roi_mars <- merge(upenn_roi_mars, label, by="RID")
picslashs <- merge(picslashs, label, by="RID")


## exclude those that have high missingess

uaspmvbm <- uaspmvbm[,colMeans(is.na(uaspmvbm)) <= 0.05]
ucsdvol <- ucsdvol[,colMeans(is.na(ucsdvol)) <= 0.05]
ucsffsx51final <- ucsffsx51final[,colMeans(is.na(ucsffsx51final)) <= 0.05]
upenn_roi_mars <- upenn_roi_mars[,colMeans(is.na(upenn_roi_mars)) <= 0.05]
picslashs <- picslashs[,colMeans(is.na(picslashs)) <= 0.05]


## preprocessing data

rownames(uaspmvbm) <- uaspmvbm$PTID
uaspmvbm <- uaspmvbm[-c(1:9, 127:128)]
table(uaspmvbm$Pheno)

rownames(ucsdvol) <- ucsdvol$PTID
ucsdvol <- ucsdvol[-c(1:7, 23:24)]
table(ucsdvol$Pheno)

rownames(ucsffsx51final) <- ucsffsx51final$PTID
ucsffsx51final <- ucsffsx51final[-c(1:20, 361:362)]
table(ucsffsx51final$Pheno)

rownames(upenn_roi_mars) <- upenn_roi_mars$PTID
upenn_roi_mars <- upenn_roi_mars[-c(1:10, 270:271)]
table(upenn_roi_mars$Pheno)

rownames(picslashs) <- picslashs$PTID
picslashs <- picslashs[-c(1:5, 56:59)]
table(picslashs$Pheno)

# Normalized

all_datasets <- vector("list")
dataset <- vector("list")

all_datasets[["uaspmvbm"]] <- uaspmvbm
all_datasets[["ucsdvol"]] <- ucsdvol
all_datasets[["upenn_roi_mars"]] <- upenn_roi_mars
all_datasets[["picslashs"]] <- picslashs
all_datasets[["ucsffsx51final"]] <- ucsffsx51final


datalist <- c("uaspmvbm", "ucsdvol", "upenn_roi_mars", 
              "picslashs", "ucsffsx51final")

for (e in datalist) {
  data <- all_datasets[[e]]
  collist <- colnames(data)
  collist <- collist[collist != "Pheno"]
  
  for (k in collist) {
    data[[k]] <- as.numeric(data[[k]])
  }
  dataset[[e]] <- data
}

## Normalised

myNormalize <- caret::preProcess(dataset[["uaspmvbm"]], methods = c("center", "scale"))
uaspmvbm_norm <- predict(myNormalize,uaspmvbm)

myNormalize <- caret::preProcess(dataset[["ucsdvol"]], methods = c("center", "scale"))
ucsdvol_norm <- predict(myNormalize,ucsdvol)

myNormalize <- caret::preProcess(dataset[["upenn_roi_mars"]], methods = c("center", "scale"))
upenn_roi_mars_norm <- predict(myNormalize, upenn_roi_mars)

myNormalize <- caret::preProcess(dataset[["picslashs"]], methods = c("center", "scale"))
picslashs_norm <- predict(myNormalize,picslashs)

myNormalize <- caret::preProcess(dataset[["ucsffsx51final"]], methods = c("center", "scale"))
ucsffsx51final_norm <- predict(myNormalize,ucsffsx51final)

## combine data

merge1 <- merge(uaspmvbm_norm, ucsdvol_norm, by=0, all=TRUE)
rownames(merge1) <- merge1$Row.names
merge1 <- merge1[,-1]
merge2 <- merge(merge1, upenn_roi_mars_norm, by=0, all= TRUE)
rownames(merge2) <- merge2$Row.names
merge2 <- merge2[,-1]
merge3 <- merge(merge2, picslashs_norm, by=0, all= TRUE)
rownames(merge3) <- merge3$Row.names
merge3 <- merge3[,-1]
merge4 <- merge(merge3, ucsffsx51final_norm, by=0, all= TRUE)
rownames(merge4) <- merge4$Row.names
merge4 <- merge4[,-1]
merged_data <- merge4
library(Amelia)
#missmap(merged_data)
dim(merged_data)


## Imputation 

vector <- colnames(merged_data)

for (k in vector) {
  merged_data[[k]] <- as.numeric(merged_data[[k]])
}

myimputation <- caret::preProcess(merged_data, "knnImpute")
library(RANN)
merged_data_imputed <- predict(myimputation, newdata = merged_data, na.action = na.pass)
#missmap(merged_data_imputed)
dim(merged_data_imputed)

## prepare for machine learning

merged_data_imputed$PTID <- rownames(merged_data_imputed)

label <- label[,-1]

merged_data_imputed <- merge(merged_data_imputed, label, by="PTID")
rownames(merged_data_imputed) <- merged_data_imputed$PTID
merged_data_imputed <- merged_data_imputed[,-1]

a <- merged_data_imputed
a <- filter(a, Pheno == "Dementia" |
              Pheno == "CN" |
              Pheno == "MCI")

table(a$Pheno)



filter_MCI <- rownames(a)[a$Pheno != "MCI"] #keep MCI samples
filter_CTL <- rownames(a)[a$Pheno != "CN"] #keep CTL samples
filter_AD <- rownames(a)[a$Pheno != "Dementia"] #keep AD samples

ExpressionSet_CTL_AD <- a[filter_MCI,] # make an ExpressionSet only with CTL_AD
ExpressionSet_MCI_AD <- a[filter_CTL,] # make an ExpressionSet only with MCI_AD
ExpressionSet_CTL_MCI <- a[filter_AD,] # make an ExpressionSet only with CTL_MCI

colnames(ExpressionSet_CTL_AD) <- make.names(as.vector(colnames(ExpressionSet_CTL_AD)), unique = FALSE, allow_ = TRUE)
colnames(ExpressionSet_MCI_AD) <- make.names(as.vector(colnames(ExpressionSet_MCI_AD)), unique = FALSE, allow_ = TRUE)
colnames(ExpressionSet_CTL_MCI) <- make.names(as.vector(colnames(ExpressionSet_CTL_MCI)), unique = FALSE, allow_ = TRUE)

a_CTL_AD <- ExpressionSet_CTL_AD
a_CTL_MCI <- ExpressionSet_CTL_MCI
a_MCI_AD <- ExpressionSet_MCI_AD

datalist <- colnames(a_CTL_AD)
datalist <- datalist[datalist != "Pheno"]
for (z in datalist) {
  a_CTL_AD[[z]] <- as.numeric(a_CTL_AD[[z]])
  a_CTL_MCI[[z]] <- as.numeric(a_CTL_MCI[[z]])
  a_MCI_AD[[z]] <- as.numeric(a_MCI_AD[[z]])
}


a_CTL_AD$Pheno <- as.factor(a_CTL_AD$Pheno)
a_CTL_MCI$Pheno <- as.factor(a_CTL_MCI$Pheno)
a_MCI_AD$Pheno <- as.factor(a_MCI_AD$Pheno)

#save(a_CTL_AD, a_CTL_MCI, a_MCI_AD, file = "raw.RData")

obj_rf_CTL_AD_5folds <- vector("list")
obj_rf_CTL_MCI_5folds <- vector("list")
obj_rf_MCI_AD_5folds <- vector("list")


obj_svm_CTL_AD_5folds <- vector("list")
obj_svm_MCI_AD_5folds <- vector("list")
obj_svm_CTL_MCI_5folds <- vector("list")

final_model_rf_CTL_AD_5folds <- vector("list")
final_model_rf_MCI_AD_5folds <- vector("list")
final_model_rf_CTL_MCI_5folds <- vector("list")

final_model_svm_CTL_AD_5folds <- vector("list")
final_model_svm_MCI_AD_5folds <- vector("list")
final_model_svm_CTL_MCI_5folds <- vector("list")

myROC_rf_CTL_AD_5folds <- vector("list")
myROC_rf_CTL_MCI_5folds <- vector("list")
myROC_rf_MCI_AD_5folds <- vector("list")

myROC_svm_CTL_AD_5folds <- vector("list")
myROC_svm_CTL_MCI_5folds <- vector("list")
myROC_svm_MCI_AD_5folds <- vector("list")


set.seed(1333)

##split partitions

a_CTL_AD_5folds <- caret::createFolds(a_CTL_AD$Pheno, k=5)
a_CTL_MCI_5folds <- caret::createFolds(a_CTL_MCI$Pheno, k=5)
a_MCI_AD_5folds <- caret::createFolds(a_MCI_AD$Pheno, k=5)

for (i in seq_along(a_CTL_AD_5folds)) {
  my_training_index_CTL_AD <- a_CTL_AD_5folds[[i]]
  a_CTL_AD_train <- a_CTL_AD[-my_training_index_CTL_AD,]
  a_CTL_AD_test <- a_CTL_AD[my_training_index_CTL_AD,]
  
  my_training_index_CTL_MCI <- a_CTL_MCI_5folds[[i]]
  a_CTL_MCI_train <- a_CTL_MCI[-my_training_index_CTL_MCI,]
  a_CTL_MCI_test <- a_CTL_MCI[my_training_index_CTL_MCI,]
  
  my_training_index_MCI_AD <- a_MCI_AD_5folds[[i]]
  a_MCI_AD_train <- a_MCI_AD[-my_training_index_MCI_AD,]
  a_MCI_AD_test <- a_MCI_AD[my_training_index_MCI_AD,]
  
  print(i)
  #fit models and save results
  modFit_rf_CTL_AD <- caret::train(Pheno ~ ., data = a_CTL_AD_train,
                                   method = "rf", ## svmLinear, etc
                                   preProcess = c("center", "scale", "nzv"),
                                   trControl = trainControl(method = "cv", number = 10,
                                                            summaryFunction=twoClassSummary, 
                                                            classProbs=T,
                                                            savePredictions = T),
                                   metric = "ROC")
  
  modFit_rf_CTL_MCI <- caret::train(Pheno ~ ., data = a_CTL_MCI_train,
                                    method = "rf", ## svmLinear, etc
                                    preProcess = c("center", "scale", "nzv"),
                                    trControl = trainControl(method = "cv", number = 10,
                                                             summaryFunction=twoClassSummary, 
                                                             classProbs=T,
                                                             savePredictions = T),metric = "ROC")
  
  modFit_rf_MCI_AD <- caret::train(Pheno ~ ., data = a_MCI_AD_train,
                                   method = "rf", ## svmLinear, etc
                                   preProcess = c("center", "scale", "nzv"),
                                   trControl = trainControl(method = "cv", number = 10,
                                                            summaryFunction=twoClassSummary, 
                                                            classProbs=T,
                                                            savePredictions = T),metric = "ROC")
  
  
  
  modFit_svm_CTL_AD <- caret::train(Pheno ~ ., data = a_CTL_AD_train,
                                    method = "svmRadial", ## svmLinear, etc
                                    preProcess = c("center", "scale", "nzv"),
                                    trControl = trainControl(method = "cv", number = 10,
                                                             summaryFunction=twoClassSummary, 
                                                             classProbs=T,
                                                             savePredictions = T),metric = "ROC")
  
  modFit_svm_CTL_MCI <- caret::train(Pheno ~ ., data = a_CTL_MCI_train,
                                     method = "svmRadial", ## svmLinear, etc
                                     preProcess = c("center", "scale", "nzv"),
                                     trControl = trainControl(method = "cv", number = 10,
                                                              summaryFunction=twoClassSummary, 
                                                              classProbs=T,
                                                              savePredictions = T),metric = "ROC")
  
  modFit_svm_MCI_AD <- caret::train(Pheno ~ ., data = a_MCI_AD_train,
                                    method = "svmRadial", ## svmLinear, etc
                                    preProcess = c("center", "scale", "nzv"),
                                    trControl = trainControl(method = "cv", number = 10,
                                                             summaryFunction=twoClassSummary, 
                                                             classProbs=T,
                                                             savePredictions = T),metric = "ROC")
  
  myVarImp_rf_CTL_AD <- varImp(modFit_rf_CTL_AD, scale = FALSE)
  myVarImp_rf_CTL_MCI <- varImp(modFit_rf_CTL_MCI, scale = FALSE)
  myVarImp_rf_MCI_AD <- varImp(modFit_rf_MCI_AD, scale = FALSE)
  
  myVarImp_svm_CTL_AD <- varImp(modFit_svm_CTL_AD, scale = FALSE)
  myVarImp_svm_CTL_MCI <- varImp(modFit_svm_CTL_MCI, scale = FALSE)
  myVarImp_svm_MCI_AD <- varImp(modFit_svm_MCI_AD, scale = FALSE)
  
  myROC_rf_CTL_AD <- pROC::roc(response = modFit_rf_CTL_AD$pred$obs, predictor = modFit_rf_CTL_AD$pred$AD, auc = TRUE, ci = TRUE)
  myROC_rf_CTL_MCI <- pROC::roc(response = modFit_rf_CTL_MCI$pred$obs, predictor = modFit_rf_CTL_MCI$pred$MCI, auc = TRUE, ci = TRUE)
  myROC_rf_MCI_AD <- pROC::roc(response = modFit_rf_MCI_AD$pred$obs, predictor = modFit_rf_MCI_AD$pred$AD, auc = TRUE, ci = TRUE)
  
  myROC_svm_CTL_AD <- pROC::roc(response = modFit_svm_CTL_AD$pred$obs, predictor = modFit_svm_CTL_AD$pred$AD, auc = TRUE, ci = TRUE)
  myROC_svm_CTL_MCI <- pROC::roc(response = modFit_svm_CTL_MCI$pred$obs, predictor = modFit_svm_CTL_MCI$pred$MCI, auc = TRUE, ci = TRUE)
  myROC_svm_MCI_AD <- pROC::roc(response = modFit_svm_MCI_AD$pred$obs, predictor = modFit_svm_MCI_AD$pred$AD, auc = TRUE, ci = TRUE)
  
  
  final_model_rf_CTL_AD_5folds[[i]] <- modFit_rf_CTL_AD
  final_model_rf_MCI_AD_5folds[[i]] <- modFit_rf_MCI_AD
  final_model_rf_CTL_MCI_5folds[[i]] <- modFit_rf_CTL_MCI
  
  final_model_svm_CTL_AD_5folds[[i]] <- modFit_svm_CTL_AD
  final_model_svm_MCI_AD_5folds[[i]] <- modFit_svm_MCI_AD
  final_model_svm_CTL_MCI_5folds[[i]] <- modFit_svm_CTL_MCI
  
  obj_rf_CTL_AD_5folds[[i]] <- myVarImp_rf_CTL_AD
  obj_rf_MCI_AD_5folds[[i]] <- myVarImp_rf_MCI_AD
  obj_rf_CTL_MCI_5folds[[i]] <- myVarImp_rf_CTL_MCI
  
  obj_svm_CTL_AD_5folds[[i]] <- myVarImp_svm_CTL_AD
  obj_svm_MCI_AD_5folds[[i]] <- myVarImp_svm_MCI_AD
  obj_svm_CTL_MCI_5folds[[i]] <- myVarImp_svm_CTL_MCI
  
  myROC_rf_CTL_AD_5folds[[i]] <- myROC_rf_CTL_AD
  myROC_rf_CTL_MCI_5folds[[i]] <- myROC_rf_CTL_MCI
  myROC_rf_MCI_AD_5folds[[i]] <- myROC_rf_MCI_AD
  
  myROC_svm_CTL_AD_5folds[[i]] <- myROC_svm_CTL_AD
  myROC_svm_CTL_MCI_5folds[[i]] <- myROC_svm_CTL_MCI
  myROC_svm_MCI_AD_5folds[[i]] <- myROC_svm_MCI_AD
  
  print("END")
}

## roc plot

myPlottedROC_rf_1_CTL_AD <- pROC::plot.roc(myROC_rf_CTL_AD_5folds[[1]], main = "rf_CTL_AD_1", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_rf_2_CTL_AD <- pROC::plot.roc(myROC_rf_CTL_AD_5folds[[2]], main = "rf_CTL_AD_2", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_rf_3_CTL_AD <- pROC::plot.roc(myROC_rf_CTL_AD_5folds[[3]], main = "rf_CTL_AD_3", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_rf_4_CTL_AD <- pROC::plot.roc(myROC_rf_CTL_AD_5folds[[4]], main = "rf_CTL_AD_4", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_rf_5_CTL_AD <- pROC::plot.roc(myROC_rf_CTL_AD_5folds[[5]], main = "rf_CTL_AD_5", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)

myPlottedROC_rf_1_MCI_AD <- pROC::plot.roc(myROC_rf_MCI_AD_5folds[[1]], main = "rf_MCI_AD_1", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_rf_2_MCI_AD <- pROC::plot.roc(myROC_rf_MCI_AD_5folds[[2]], main = "rf_MCI_AD_2", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_rf_3_MCI_AD <- pROC::plot.roc(myROC_rf_MCI_AD_5folds[[3]], main = "rf_MCI_AD_3", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_rf_4_MCI_AD <- pROC::plot.roc(myROC_rf_MCI_AD_5folds[[4]], main = "rf_MCI_AD_4", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_rf_5_MCI_AD <- pROC::plot.roc(myROC_rf_MCI_AD_5folds[[5]], main = "rf_MCI_AD_5", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)

myPlottedROC_rf_1_CTL_MCI <- pROC::plot.roc(myROC_rf_CTL_MCI_5folds[[1]], main = "rf_CTL_MCI_1", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_rf_2_CTL_MCI <- pROC::plot.roc(myROC_rf_CTL_MCI_5folds[[2]], main = "rf_CTL_MCI_2", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_rf_3_CTL_MCI <- pROC::plot.roc(myROC_rf_CTL_MCI_5folds[[3]], main = "rf_CTL_MCI_3", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_rf_4_CTL_MCI <- pROC::plot.roc(myROC_rf_CTL_MCI_5folds[[4]], main = "rf_CTL_MCI_4", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_rf_5_CTL_MCI <- pROC::plot.roc(myROC_rf_CTL_MCI_5folds[[5]], main = "rf_CTL_MCI_5", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)


##for SMV

myPlottedROC_svm_1_CTL_AD <- pROC::plot.roc(myROC_svm_CTL_AD_5folds[[1]], main = "svm_CTL_AD_1", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_svm_2_CTL_AD <- pROC::plot.roc(myROC_svm_CTL_AD_5folds[[2]], main = "svm_CTL_AD_2", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_svm_3_CTL_AD <- pROC::plot.roc(myROC_svm_CTL_AD_5folds[[3]], main = "svm_CTL_AD_3", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_svm_4_CTL_AD <- pROC::plot.roc(myROC_svm_CTL_AD_5folds[[4]], main = "svm_CTL_AD_4", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_svm_5_CTL_AD <- pROC::plot.roc(myROC_svm_CTL_AD_5folds[[5]], main = "svm_CTL_AD_5", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)

myPlottedROC_svm_1_MCI_AD <- pROC::plot.roc(myROC_svm_MCI_AD_5folds[[1]], main = "svm_MCI_AD_1", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_svm_2_MCI_AD <- pROC::plot.roc(myROC_svm_MCI_AD_5folds[[2]], main = "svm_MCI_AD_2", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_svm_3_MCI_AD <- pROC::plot.roc(myROC_svm_MCI_AD_5folds[[3]], main = "svm_MCI_AD_3", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_svm_4_MCI_AD <- pROC::plot.roc(myROC_svm_MCI_AD_5folds[[4]], main = "svm_MCI_AD_4", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_svm_5_MCI_AD <- pROC::plot.roc(myROC_svm_MCI_AD_5folds[[5]], main = "svm_MCI_AD_5", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)

myPlottedROC_svm_1_CTL_MCI <- pROC::plot.roc(myROC_svm_CTL_MCI_5folds[[1]], main = "svm_CTL_MCI_1", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_svm_2_CTL_MCI <- pROC::plot.roc(myROC_svm_CTL_MCI_5folds[[2]], main = "svm_CTL_MCI_2", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_svm_3_CTL_MCI <- pROC::plot.roc(myROC_svm_CTL_MCI_5folds[[3]], main = "svm_CTL_MCI_3", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_svm_4_CTL_MCI <- pROC::plot.roc(myROC_svm_CTL_MCI_5folds[[4]], main = "svm_CTL_MCI_4", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)
myPlottedROC_svm_5_CTL_MCI <- pROC::plot.roc(myROC_svm_CTL_MCI_5folds[[5]], main = "svm_CTL_MCI_5", print.auc=TRUE, print.thres=TRUE, ci=TRUE,asp=FALSE)


total_rf_important <- vector("list")
total_svm_important <- vector("list")

# RANDOM FOREST
rf_CTL_AD_1 <- data.frame(obj_rf_CTL_AD_5folds[[1]]$importance)
rf_CTL_AD_2 <- data.frame(obj_rf_CTL_AD_5folds[[2]]$importance)
rf_CTL_AD_3 <- data.frame(obj_rf_CTL_AD_5folds[[3]]$importance)
rf_CTL_AD_4 <- data.frame(obj_rf_CTL_AD_5folds[[4]]$importance)
rf_CTL_AD_5 <- data.frame(obj_rf_CTL_AD_5folds[[5]]$importance)

rf_MCI_AD_1 <- data.frame(obj_rf_MCI_AD_5folds[[1]]$importance)
rf_MCI_AD_2 <- data.frame(obj_rf_MCI_AD_5folds[[2]]$importance)
rf_MCI_AD_3 <- data.frame(obj_rf_MCI_AD_5folds[[3]]$importance)
rf_MCI_AD_4 <- data.frame(obj_rf_MCI_AD_5folds[[4]]$importance)
rf_MCI_AD_5 <- data.frame(obj_rf_MCI_AD_5folds[[5]]$importance)

rf_CTL_MCI_1 <- data.frame(obj_rf_CTL_MCI_5folds[[1]]$importance)
rf_CTL_MCI_2 <- data.frame(obj_rf_CTL_MCI_5folds[[2]]$importance)
rf_CTL_MCI_3 <- data.frame(obj_rf_CTL_MCI_5folds[[3]]$importance)
rf_CTL_MCI_4 <- data.frame(obj_rf_CTL_MCI_5folds[[4]]$importance)
rf_CTL_MCI_5 <- data.frame(obj_rf_CTL_MCI_5folds[[5]]$importance)

# SVM

svm_CTL_AD_1 <- data.frame(obj_svm_CTL_AD_5folds[[1]]$importance)
svm_CTL_AD_2 <- data.frame(obj_svm_CTL_AD_5folds[[2]]$importance)
svm_CTL_AD_3 <- data.frame(obj_svm_CTL_AD_5folds[[3]]$importance)
svm_CTL_AD_4 <- data.frame(obj_svm_CTL_AD_5folds[[4]]$importance)
svm_CTL_AD_5 <- data.frame(obj_svm_CTL_AD_5folds[[5]]$importance)

svm_MCI_AD_1 <- data.frame(obj_svm_MCI_AD_5folds[[1]]$importance)
svm_MCI_AD_2 <- data.frame(obj_svm_MCI_AD_5folds[[2]]$importance)
svm_MCI_AD_3 <- data.frame(obj_svm_MCI_AD_5folds[[3]]$importance)
svm_MCI_AD_4 <- data.frame(obj_svm_MCI_AD_5folds[[4]]$importance)
svm_MCI_AD_5 <- data.frame(obj_svm_MCI_AD_5folds[[5]]$importance)

svm_CTL_MCI_1 <- data.frame(obj_svm_CTL_MCI_5folds[[1]]$importance)
svm_CTL_MCI_2 <- data.frame(obj_svm_CTL_MCI_5folds[[2]]$importance)
svm_CTL_MCI_3 <- data.frame(obj_svm_CTL_MCI_5folds[[3]]$importance)
svm_CTL_MCI_4 <- data.frame(obj_svm_CTL_MCI_5folds[[4]]$importance)
svm_CTL_MCI_5 <- data.frame(obj_svm_CTL_MCI_5folds[[5]]$importance)

# CTL - AD

rf_CTL_AD_1$GENES <- rownames(rf_CTL_AD_1)
rf_CTL_AD_2$GENES <- rownames(rf_CTL_AD_2)
rf_CTL_AD_3$GENES <- rownames(rf_CTL_AD_3)
rf_CTL_AD_4$GENES <- rownames(rf_CTL_AD_4)
rf_CTL_AD_5$GENES <- rownames(rf_CTL_AD_5)

colnames(rf_CTL_AD_1) <- c("Values_1", "GENES")
colnames(rf_CTL_AD_2) <- c("Values_2", "GENES")
colnames(rf_CTL_AD_3) <- c("Values_3", "GENES")
colnames(rf_CTL_AD_4) <- c("Values_4", "GENES")
colnames(rf_CTL_AD_5) <- c("Values_5", "GENES")

total_rf_CTL_AD <- join_all(list(rf_CTL_AD_1,rf_CTL_AD_2,rf_CTL_AD_3,rf_CTL_AD_4,rf_CTL_AD_5), by = 'GENES')
row.names(total_rf_CTL_AD) <- total_rf_CTL_AD$GENES
total_rf_CTL_AD <- total_rf_CTL_AD[-2]
total_rf_CTL_AD$Total <- rowSums(total_rf_CTL_AD)
total_rf_CTL_AD <- total_rf_CTL_AD[order(-total_rf_CTL_AD$Total),]
total_rf_important[["CTL_AD"]] <- total_rf_CTL_AD

# CTL - MCI

rf_CTL_MCI_1$GENES <- rownames(rf_CTL_MCI_1)
rf_CTL_MCI_2$GENES <- rownames(rf_CTL_MCI_2)
rf_CTL_MCI_3$GENES <- rownames(rf_CTL_MCI_3)
rf_CTL_MCI_4$GENES <- rownames(rf_CTL_MCI_4)
rf_CTL_MCI_5$GENES <- rownames(rf_CTL_MCI_5)

colnames(rf_CTL_MCI_1) <- c("Values_1", "GENES")
colnames(rf_CTL_MCI_2) <- c("Values_2", "GENES")
colnames(rf_CTL_MCI_3) <- c("Values_3", "GENES")
colnames(rf_CTL_MCI_4) <- c("Values_4", "GENES")
colnames(rf_CTL_MCI_5) <- c("Values_5", "GENES")

total_rf_CTL_MCI <- join_all(list(rf_CTL_MCI_1,rf_CTL_MCI_2,rf_CTL_MCI_3,rf_CTL_MCI_4,rf_CTL_MCI_5), by = 'GENES')
row.names(total_rf_CTL_MCI) <- total_rf_CTL_MCI$GENES
total_rf_CTL_MCI <- total_rf_CTL_MCI[-2]
total_rf_CTL_MCI$Total <- rowSums(total_rf_CTL_MCI)
total_rf_CTL_MCI <- total_rf_CTL_MCI[order(-total_rf_CTL_MCI$Total),]
total_rf_important[["CTL_MCI"]] <- total_rf_CTL_MCI

# MCI - AD

rf_MCI_AD_1$GENES <- rownames(rf_MCI_AD_1)
rf_MCI_AD_2$GENES <- rownames(rf_MCI_AD_2)
rf_MCI_AD_3$GENES <- rownames(rf_MCI_AD_3)
rf_MCI_AD_4$GENES <- rownames(rf_MCI_AD_4)
rf_MCI_AD_5$GENES <- rownames(rf_MCI_AD_5)

colnames(rf_MCI_AD_1) <- c("Values_1", "GENES")
colnames(rf_MCI_AD_2) <- c("Values_2", "GENES")
colnames(rf_MCI_AD_3) <- c("Values_3", "GENES")
colnames(rf_MCI_AD_4) <- c("Values_4", "GENES")
colnames(rf_MCI_AD_5) <- c("Values_5", "GENES")

total_rf_MCI_AD <- join_all(list(rf_MCI_AD_1,rf_MCI_AD_2,rf_MCI_AD_3,rf_MCI_AD_4,rf_MCI_AD_5), by = 'GENES')
row.names(total_rf_MCI_AD) <- total_rf_MCI_AD$GENES
total_rf_MCI_AD <- total_rf_MCI_AD[-2]
total_rf_MCI_AD$Total <- rowSums(total_rf_MCI_AD)
total_rf_MCI_AD <- total_rf_MCI_AD[order(-total_rf_MCI_AD$Total),]
total_rf_important[["MCI_AD"]] <- total_rf_MCI_AD


## SUPPORT VECTOR MACHINE

svm_CTL_AD_1 <- svm_CTL_AD_1[1]
svm_CTL_AD_2 <- svm_CTL_AD_2[1]
svm_CTL_AD_3 <- svm_CTL_AD_3[1]
svm_CTL_AD_4 <- svm_CTL_AD_4[1]
svm_CTL_AD_5 <- svm_CTL_AD_5[1]

svm_CTL_AD_1$GENES <- rownames(svm_CTL_AD_1)
svm_CTL_AD_2$GENES <- rownames(svm_CTL_AD_2)
svm_CTL_AD_3$GENES <- rownames(svm_CTL_AD_3)
svm_CTL_AD_4$GENES <- rownames(svm_CTL_AD_4)
svm_CTL_AD_5$GENES <- rownames(svm_CTL_AD_5)

colnames(svm_CTL_AD_1) <- c("Values_1", "GENES")
colnames(svm_CTL_AD_2) <- c("Values_2", "GENES")
colnames(svm_CTL_AD_3) <- c("Values_3", "GENES")
colnames(svm_CTL_AD_4) <- c("Values_4", "GENES")
colnames(svm_CTL_AD_5) <- c("Values_5", "GENES")

total_svm_CTL_AD <- join_all(list(svm_CTL_AD_1,svm_CTL_AD_2,svm_CTL_AD_3,svm_CTL_AD_4,svm_CTL_AD_5), by = 'GENES')
row.names(total_svm_CTL_AD) <- total_svm_CTL_AD$GENES
total_svm_CTL_AD <- total_svm_CTL_AD[-2]
total_svm_CTL_AD$Total <- rowSums(total_svm_CTL_AD)
total_svm_CTL_AD <- total_svm_CTL_AD[order(-total_svm_CTL_AD$Total),]
total_svm_important[["CTL_AD"]] <- total_svm_CTL_AD

# CTL - MCI

svm_CTL_MCI_1 <- svm_CTL_MCI_1[1]
svm_CTL_MCI_2 <- svm_CTL_MCI_2[1]
svm_CTL_MCI_3 <- svm_CTL_MCI_3[1]
svm_CTL_MCI_4 <- svm_CTL_MCI_4[1]
svm_CTL_MCI_5 <- svm_CTL_MCI_5[1]

svm_CTL_MCI_1$GENES <- rownames(svm_CTL_MCI_1)
svm_CTL_MCI_2$GENES <- rownames(svm_CTL_MCI_2)
svm_CTL_MCI_3$GENES <- rownames(svm_CTL_MCI_3)
svm_CTL_MCI_4$GENES <- rownames(svm_CTL_MCI_4)
svm_CTL_MCI_5$GENES <- rownames(svm_CTL_MCI_5)

colnames(svm_CTL_MCI_1) <- c("Values_1", "GENES")
colnames(svm_CTL_MCI_2) <- c("Values_2", "GENES")
colnames(svm_CTL_MCI_3) <- c("Values_3", "GENES")
colnames(svm_CTL_MCI_4) <- c("Values_4", "GENES")
colnames(svm_CTL_MCI_5) <- c("Values_5", "GENES")

total_svm_CTL_MCI <- join_all(list(svm_CTL_MCI_1,svm_CTL_MCI_2,svm_CTL_MCI_3,svm_CTL_MCI_4,svm_CTL_MCI_5), by = 'GENES')
row.names(total_svm_CTL_MCI) <- total_svm_CTL_MCI$GENES
total_svm_CTL_MCI <- total_svm_CTL_MCI[-2]
total_svm_CTL_MCI$Total <- rowSums(total_svm_CTL_MCI)
total_svm_CTL_MCI <- total_svm_CTL_MCI[order(-total_svm_CTL_MCI$Total),]
total_svm_important[["CTL_MCI"]] <- total_svm_CTL_MCI

# MCI - AD

svm_MCI_AD_1 <- svm_MCI_AD_1[1]
svm_MCI_AD_2 <- svm_MCI_AD_2[1]
svm_MCI_AD_3 <- svm_MCI_AD_3[1]
svm_MCI_AD_4 <- svm_MCI_AD_4[1]
svm_MCI_AD_5 <- svm_MCI_AD_5[1]

svm_MCI_AD_1$GENES <- rownames(svm_MCI_AD_1)
svm_MCI_AD_2$GENES <- rownames(svm_MCI_AD_2)
svm_MCI_AD_3$GENES <- rownames(svm_MCI_AD_3)
svm_MCI_AD_4$GENES <- rownames(svm_MCI_AD_4)
svm_MCI_AD_5$GENES <- rownames(svm_MCI_AD_5)

colnames(svm_MCI_AD_1) <- c("Values_1", "GENES")
colnames(svm_MCI_AD_2) <- c("Values_2", "GENES")
colnames(svm_MCI_AD_3) <- c("Values_3", "GENES")
colnames(svm_MCI_AD_4) <- c("Values_4", "GENES")
colnames(svm_MCI_AD_5) <- c("Values_5", "GENES")

total_svm_MCI_AD <- join_all(list(svm_MCI_AD_1,svm_MCI_AD_2,svm_MCI_AD_3,svm_MCI_AD_4,svm_MCI_AD_5), by = 'GENES')
row.names(total_svm_MCI_AD) <- total_svm_MCI_AD$GENES
total_svm_MCI_AD <- total_svm_MCI_AD[-2]
total_svm_MCI_AD$Total <- rowSums(total_svm_MCI_AD)
total_svm_MCI_AD <- total_svm_MCI_AD[order(-total_svm_MCI_AD$Total),]
total_svm_important[["MCI_AD"]] <- total_svm_MCI_AD

barplot(total_rf_important[["MCI_AD"]]$Total, main = "RF_MCI_AD")
barplot(total_rf_important[["CTL_MCI"]]$Total, main = "RF_CTL_MCI")
barplot(total_rf_important[["CTL_AD"]]$Total, main = "RF_CTL_AD")

barplot(total_svm_important[["MCI_AD"]]$Total, main = "SVM_MCI_AD")
barplot(total_svm_important[["CTL_MCI"]]$Total, main = "SVM_CTL_MCI")
barplot(total_svm_important[["CTL_AD"]]$Total, main = "SVM_CTL_AD")

## random forest

CTL_AD_RF <- total_rf_important[["CTL_AD"]]
CTL_MCI_RF <- total_rf_important[["CTL_MCI"]]
MCI_AD_RF <- total_rf_important[["MCI_AD"]]

CTL_AD_RF$GENES<- row.names(CTL_AD_RF)
CTL_MCI_RF$GENES<- row.names(CTL_MCI_RF)
MCI_AD_RF$GENES<- row.names(MCI_AD_RF)

CTL_AD_RF_filt <- filter(CTL_AD_RF, Total >= 4) # 60
CTL_MCI_RF_filt <- filter(CTL_MCI_RF, Total >= 6) # 64
MCI_AD_RF_filt <- filter(MCI_AD_RF, Total >= 5) # 67

dim(CTL_AD_RF_filt)
dim(CTL_MCI_RF_filt)
dim(MCI_AD_RF_filt)

## support vector maachine

CTL_AD_SVM <- total_svm_important[["CTL_AD"]]
CTL_MCI_SVM <- total_svm_important[["CTL_MCI"]]
MCI_AD_SVM <- total_svm_important[["MCI_AD"]]

CTL_AD_SVM$GENES<- row.names(CTL_AD_SVM)
CTL_MCI_SVM$GENES<- row.names(CTL_MCI_SVM)
MCI_AD_SVM$GENES<- row.names(MCI_AD_SVM)

CTL_AD_SVM_filt <- filter(CTL_AD_SVM, Total >= 3.7) # 22
CTL_MCI_SVM_filt <- filter(CTL_MCI_SVM, Total >= 2.8) # 23
MCI_AD_SVM_filt <- filter(MCI_AD_SVM, Total >= 3.5) # 26

dim(CTL_AD_SVM_filt)
dim(CTL_MCI_SVM_filt)
dim(MCI_AD_SVM_filt)

## FIND ALL FEATURES FROM RF & SVM

test1 <- CTL_AD_SVM_filt[!(CTL_AD_SVM_filt$GENES %in% CTL_AD_RF_filt$GENES),]
test2 <- CTL_MCI_SVM_filt[!(CTL_MCI_SVM_filt$GENES %in% CTL_MCI_RF_filt$GENES),]
test3 <- MCI_AD_SVM_filt[!(MCI_AD_SVM_filt$GENES %in% MCI_AD_RF_filt$GENES),]

## 9 FEATURES MORE IN CTL_MCI_SVM ADD THEM IN EXTRACTION

vector_CTL_MCI_test2 <- test2$GENES

## EXTRACT SIGNIFICANT VARIANTS

vector_CTL_AD <- CTL_AD_RF_filt$GENES
vector_CTL_MCI <- CTL_MCI_RF_filt$GENES
vector_CTL_MCI <- c(vector_CTL_MCI, vector_CTL_MCI_test2 )
vector_MCI_AD <- MCI_AD_RF_filt$GENES

CTL_AD_final <- a_CTL_AD[,vector_CTL_AD]
CTL_MCI_final <- a_CTL_MCI[,vector_CTL_MCI]
MCI_AD_final <- a_MCI_AD[,vector_MCI_AD]

dim(CTL_AD_final)
dim(CTL_MCI_final)
dim(MCI_AD_final)

save(CTL_AD_final, CTL_MCI_final, MCI_AD_final, file = "MRI_results.RData")