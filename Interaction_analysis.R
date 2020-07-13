## load libraries
library(tidyverse)
library(caret)
library(glmnet)

## convert GWAS DATA into appropirate file format

vcf <- read.vcfR("CTL_AD_total.vcf", verbose = FALSE )
vcf_snps <- as.data.frame(extract.gt(vcf, element = 'GT', as.numeric = FALSE))

vcf_snps[vcf_snps == "0/0"] <- 0
vcf_snps[is.na(vcf_snps)] <- 0
vcf_snps[vcf_snps == "0/1"] <- 1
vcf_snps[vcf_snps == "1/1"] <- 2

colnames(vcf_snps) <- test

genotype <- as.data.frame(t(vcf_snps))
write.table(genotype, file ="GWAS_data.txt", row.names=F, quote = F)

## loading and find common samples

#MRI data
load("MRI_results.RData")
MRI_CTL_AD <- CTL_AD_final
#GWAS data
GWAS_CTL_AD <- read.csv(file = "GWAS_data.csv", header =  TRUE)
rownames(GWAS_CTL_AD) <- GWAS_CTL_AD[,1]
GWAS_CTL_AD <- GWAS_CTL_AD[,-1]
#transcriptomic data
load("ADNI.RData")
Transcriptomic_CTL_AD <- as.data.frame(Biobase::exprs(ExpressionSet_adni))
Transcriptomic_CTL_AD <- as.data.frame(t(Transcriptomic_CTL_AD))
load("total_genes.RData")
vector_CTL_AD <- CTL_AD_FINAL$GENES
Transcriptomic_CTL_AD <- Transcriptomic_CTL_AD[,vector_CTL_AD]

save(Transcriptomic_CTL_AD, GWAS_CTL_AD, MRI_CTL_AD, file = "Laura_MRI_GWAS_EXPR.RData")

## find common samples

load("Laura_MRI_GWAS_EXPR.RData")

myNormalize <- caret::preProcess(Transcriptomic_CTL_AD, methods = c("center", "scale"))
Transcriptomic_CTL_AD_norm <- predict(myNormalize,Transcriptomic_CTL_AD)

myNormalize <- caret::preProcess(GWAS_CTL_AD, methods = c("center", "scale"))
GWAS_CTL_AD_norm <- predict(myNormalize,GWAS_CTL_AD)

myNormalize <- caret::preProcess(MRI_CTL_AD, methods = c("center", "scale"))
MRI_CTL_AD_norm <- predict(myNormalize,MRI_CTL_AD)

## ploting normalised data

boxplot(Transcriptomic_CTL_AD_norm)
boxplot(GWAS_CTL_AD_norm)
boxplot(MRI_CTL_AD_norm)


## preprocessing for merge

rownames(Transcriptomic_CTL_AD_norm) <- toupper(rownames(Transcriptomic_CTL_AD_norm))
rownames(GWAS_CTL_AD_norm) <- toupper(rownames(GWAS_CTL_AD_norm))
rownames(MRI_CTL_AD_norm) <- toupper(rownames(MRI_CTL_AD_norm))

round1 <- GWAS_CTL_AD_norm[rownames(GWAS_CTL_AD_norm) %in% rownames(MRI_CTL_AD_norm),]
round2 <- MRI_CTL_AD_norm[rownames(MRI_CTL_AD_norm) %in% rownames(GWAS_CTL_AD_norm),]

all(rownames(round1)==rownames(round2)) # should be true

round_total <- Transcriptomic_CTL_AD_norm[rownames(Transcriptomic_CTL_AD_norm) %in% rownames(round1),]

common_samples <- rownames(round_total)

GWAS_CTL_AD <- GWAS_CTL_AD_norm[rownames(GWAS_CTL_AD_norm) %in% common_samples,]
Transcriptomic_CTL_AD <- Transcriptomic_CTL_AD_norm[rownames(Transcriptomic_CTL_AD_norm) %in% common_samples,]
MRI_CTL_AD <- MRI_CTL_AD_norm[rownames(MRI_CTL_AD_norm) %in% common_samples,]

merge1 <- merge(MRI_CTL_AD, Transcriptomic_CTL_AD, by = 0)
rownames(merge1) <- merge1$Row.names
merge1 <- merge1[,-1]
merged_data <- merge(merge1, GWAS_CTL_AD, by = 0)
rownames(merged_data) <- merged_data$Row.names
merged_data <- merged_data[,-1]

label <- read.csv(file = "labels.csv")
label <-label[!duplicated(label[,c("RID")]),]
label <- label[,c("RID" , "PTID", "DX_bl")]
label$DX_bl <- str_replace_all(label$DX_bl, "LMCI", "MCI")
label$DX_bl <- str_replace_all(label$DX_bl, "SMC", "MCI")
label$DX_bl <- str_replace_all(label$DX_bl, "EMCI", "MCI")
label$DX_bl <- str_replace_all(label$DX_bl, "AD", "Dementia")
colnames(label) <- c("RID" ,"PTID", "Pheno")

merged_data$PTID <- rownames(merged_data)
merged_data <- merge(merged_data, label, by = "PTID")
merged_data <- merged_data[,-1]

save(merged_data, file = "interaction_matrix.RData")


## start LASSO workflow

load("interaction_matrix.RData")

n<-100

ListError<-vector(mode="double", length=n)
BetasFin<-vector(mode="character", length=n)
LambdaFin<-vector(mode="double", length=n)
BNumFin<-vector(mode="double", length=n)
see2<-data.frame(All="All")
LauCoef1L<-data.frame(Coeff="See",stringsAsFactors=FALSE)
BetasTodoL<-data.frame(Features="Name",Coefficients=1)

set.seed(1333)

for(i in 1:n) {
  
  partisions <- caret::createDataPartition(merged_data$Pheno, p = 0.80, list = FALSE, times = 1)
  train <- merged_data[partisions, ]
  test <- merged_data[-partisions, ]
  
  xtrain <- train[ , !(names(train) %in% "Pheno")] %>% data.matrix()
  ytrain = train$Pheno
  xtest <- test[ , !(names(test) %in% "Pheno")] %>% data.matrix()
  ytest = test$Pheno
  
  library(glmnet)
  CV=cv.glmnet(xtrain,ytrain,family="binomial",type.measure="class",alpha=1,nlambda=100)
  
  attach(CV)
  
  Lambda.Best<-CV$lambda.min
  CVFin=glmnet(xtrain,ytrain,family="binomial",alpha=1,lambda=Lambda.Best)
  Coef<-coef(CVFin)
  Intercept<-Coef@x[1]
  Betas<-CVFin$beta
  Betas2<-data.frame(Features=Betas@Dimnames[[1]][Betas@i+1], Coefficients=Betas@x)
  CVPred1 = predict(CVFin, family="binomial", s=Lambda.Best, newx = xtest,type="class")
  
  
  
  #Calculate error for categorical values
  ytest2<-as.factor(ytest)
  
  Results<-table(CVPred1,ytest)
  Accuracy<-(Results[1]+Results[4])/sum(Results[1:4])
  Error<-1-Accuracy
  
  #visual display of for
  
  BetasTodoL<-rbind(BetasTodoL,Betas2)
  see<-Betas2$Features
  see1<-data.frame(All=see)
  see2<-rbind(see2,see1)
  m<-count(see2)
  
  see3<-toString(see)
  ListError[i]<-Error
  BetasFin[i]<-see3
  BNumFin[i]<-length(see)
  LambdaFin[i]<-Lambda.Best
  detach(CV)
  
  
}

All_info<-data.frame(Error=ListError, BetasNames=BetasFin, BetasNum=BNumFin, Lambda=LambdaFin) 
m<-count(see2, var = "All")
m<-m[-1,]
Final_LASSO<-m[order(-m$freq),]
Final_LASSO1<-filter(Final_LASSO,freq>50)

### RUN RANDOM FOREST

library(caret)
library(tidyverse)

load("interaction_matrix.RData")

a_CTL_AD <- merged_data

a_CTL_AD$Pheno <- as.factor(a_CTL_AD$Pheno)

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

a_CTL_AD_5folds <- caret::createFolds(a_CTL_AD$Pheno, k=5)

for (i in seq_along(a_CTL_AD_5folds)) {
  print(i)
  my_training_index_CTL_AD <- a_CTL_AD_5folds[[i]]
  a_CTL_AD_train <- a_CTL_AD[-my_training_index_CTL_AD,]
  a_CTL_AD_test <- a_CTL_AD[my_training_index_CTL_AD,]
  
  modFit_rf_CTL_AD <- caret::train(Pheno ~ .^2, data = a_CTL_AD_train,
                                   method = "rf", ## svmLinear, etc
                                   preProcess = c("center", "scale", "nzv"),
                                   trControl = trainControl(method = "cv", number = 10,
                                                            summaryFunction=twoClassSummary,
                                                            classProbs=T,
                                                            savePredictions = T),
                                   metric = "ROC")
  
  myVarImp_rf_CTL_AD <- varImp(modFit_rf_CTL_AD, scale = FALSE)
  
  myROC_rf_CTL_AD <- pROC::roc(response = modFit_rf_CTL_AD$pred$obs,
                               predictor = modFit_rf_CTL_AD$pred$Dementia, auc = TRUE, ci = TRUE)
  
  final_model_rf_CTL_AD_5folds[[i]] <- modFit_rf_CTL_AD
  
  obj_rf_CTL_AD_5folds[[i]] <- myVarImp_rf_CTL_AD
  
  myROC_rf_CTL_AD_5folds[[i]] <- myROC_rf_CTL_AD
  
}

save(final_model_rf_CTL_AD_5folds, final_model_svm_CTL_AD_5folds,
     obj_rf_CTL_AD_5folds, obj_svm_CTL_AD_5folds,
     myROC_rf_CTL_AD_5folds, myROC_svm_CTL_AD_5folds, file ="rf_interaction_output.RData")

## filtering from rf _ svm 

load("rf_interaction_output.RData")

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

## check alongside with lasso, rf, svm

gwas <- read.csv("GWAS_data.csv", header = T)
gwas <- colnames(gwas) 
expr <- read.table("expr_output.txt", header = T)
expr <- colnames(expr)
mri <- read.table("MRI_output.txt", header = T)
mri <- colnames(mri)
## CHECK HOW MANY EACH CATEGORY

test <- filter(total_rf_CTL_AD, Total > 0.05)
vector_test <- rownames(test)
length(vector_test) # 484 interactions with 0.05
fix(vector_test) ## replace : with "," by ctrl h
vector_test <-unique(vector_test)
length(vector_test)

snps_rf <- gwas[gwas %in% vector_test]
imag_rf <- mri[mri %in% vector_test]
expr_rf <- expr[expr %in% vector_test]

length(snps_rf) # 147 0.05
length(imag_rf) # 33 0.05
length(expr_rf) # 7 0.05

## END