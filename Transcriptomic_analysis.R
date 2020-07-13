## Libraries
library(tidyverse) # multitools
library(ggplot2) # for plotting
library(Biobase) # for ExpressionSet
library(oligo) # for oligoplotting
library(RColorBrewer) # for colourploting
library(tidymodels) # for non-linear analysis
library(ADNIMERGE) # load phenotypes from ADNI data
library(biomaRt) # download names for genes mapping
library(WGCNA) # for multi probs to one gene
library(impute) # for multi probs to one gene, it is needed for WGCNA
library(limma) # Run differencial expression with limma
library(a4Base) # combine ExpressionSets
library(sva) # batch correction
library(caret) # data spliting function for bootstraping
library(svMisc) # library necessary for looping (progress)
library(plyr) # FOR multiple merging
library(rminer) # Load significant variants from SVM
library(pROC) # for ROC ploting


#CREATE EXPRESSIONSET FROM ADNI

# load metadata table from installed ADNIMerge package
adni <- ADNIMERGE::adnimerge
# load transcript data
transcripts <- data.table::fread("ADNI_Gene_Expression_Profile.csv")
# get subject IDs from the transcripts
TranSubjects <- transcripts[3,]
# subset the adni metadata by subject IDs
adniTX <- adni[adni$PTID %in% TranSubjects,]
# select the subject ID and the diagnosis. There are repeated fields (data is in long format)
# so unique them, then look at the diagnosis field and table it to count them. You can use
# the PTID and DX fields with a unique() to get a list of which subjects are normal and which are
# diagnosed. CN == cognitive normal, MCI == mild cognitive impairment, and Dementia == Alzheimer's Disease
adniTX %>% dplyr::select(PTID,DX) %>% unique %>% dplyr::select(DX) %>% table

## --> PREPATE DATA
##############input ADNI_Gene_Expression_Profile.csv

ADNI_Gene_Expression_Profile <- read.csv(file="ADNI_Gene_Expression_Profile.csv", header = FALSE, sep=",")

############## ---- process phenofile
phenotype<- subset(head(ADNI_Gene_Expression_Profile, n=8), select=-c(V2,V3,V748))
rownames(phenotype) <- phenotype[,1]
phenotype <- phenotype[,-1]
phenotype <- as.data.frame(t(phenotype))
phenotype <- phenotype[order(phenotype$SubjectID),]

adni_meta <- adniTX
adni_meta<-adni_meta[!duplicated(adni_meta[,c("PTID")]),]
colnames(adni_meta)[4] <- "SubjectID"
adni_meta <- adni_meta[order(adni_meta$SubjectID),]
phenofile <- merge(phenotype, adni_meta, by="SubjectID")


rownames(phenofile) <- phenofile[,1]
phenofile <- as.data.frame(t(phenofile))


############## --- process featuresfile

featuresfile <- ADNI_Gene_Expression_Profile[-c(1:8), c("V1","V2","V3","V748")]
names(featuresfile) <- as.matrix(featuresfile[1,])
rownames(featuresfile) <- featuresfile[,1]
featuresfile <- featuresfile[-1,]
names(featuresfile)[4] <- "ANNO"
name_vector<-rownames(featuresfile)
featuresfile <- featuresfile %>% mutate_if(is.factor, as.character)
rownames(featuresfile) <- name_vector


############## ---- process expressionfile
expressionfile <- ADNI_Gene_Expression_Profile[-c(1,2,4:9), -c(2,3,748)]
rownames(expressionfile) <- expressionfile[,1]
expressionfile <- expressionfile[,-1]
expressionfile <- as.data.frame(t(expressionfile))
expressionfile <- expressionfile[order(expressionfile$SubjectID),]
rownames(expressionfile) <- expressionfile[,1]
expressionfile <- expressionfile[,-1]
vector_name <- as.vector(colnames(expressionfile))
expressionfile <- as.data.frame(t(expressionfile))

#### change the table's variants into numeric

expressionfile <- expressionfile %>% mutate_if(is.factor, as.character)
expressionfile <- expressionfile %>% mutate_if(is.character, as.numeric)
rownames(expressionfile) <- vector_name

exprs<- as.matrix(expressionfile)

##############create from phenotype
pData <- phenofile
pData<- t(pData)
pData<-as.data.frame(pData)

summary(pData)

##############check for similarities -- essential for ExpressionSet

all(rownames(pData)==colnames(exprs)) # should be TRUE

##############convert AnnotatedDataFrame 

phenoData <- new("AnnotatedDataFrame", data=pData)

##############create from features and convert AnnotatedDataFrame

fData <- featuresfile

featureData <- new("AnnotatedDataFrame", data=fData)

############## check for similarities -- essential ExpressionSet

all(rownames(featureData)==rownames(exprs)) # should be TRUE

##############  Create common ExpressionSet #####################

ExpressionSet <- ExpressionSet(assayData = exprs, phenoData = phenoData, featureData =featureData)


######## LOAD FILES FOR PROCESSING

#Load files
ADNI <- ExpressionSet
GSE63063 <- readRDS('GSE63063.Rds')
NCBI1 <- GSE63063[[1]]
NCBI2 <- GSE63063[[2]]

#### MULTIPROBES TO ONE GENE
# Prepare for WGCNA

# ADNI

ADNI_datET <- Biobase::exprs(ADNI) # counts data.frame
ADNI_rowGroup <- Biobase::fData(ADNI) # gene names (common)
ADNI_rowID <- row.names(Biobase::fData(ADNI)) #probes names (unique)
ADNI_rowGroup$Symbol <- gsub('10-Mar', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup$Symbol <- gsub('10-Sep', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup$Symbol <- gsub('11-Mar', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup$Symbol <- gsub('11-Sep', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup$Symbol <- gsub('12-Sep', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup$Symbol <- gsub('14-Sep', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup$Symbol <- gsub('1-Dec', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup$Symbol <- gsub('1-Mar', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup$Symbol <- gsub('1-Sep', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup$Symbol <- gsub('2-Mar', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup$Symbol <- gsub('2-Sep', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup$Symbol <- gsub('3-Mar', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup$Symbol <- gsub('3-Sep', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup$Symbol <- gsub('4-Mar', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup$Symbol <- gsub('4-Sep', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup$Symbol <- gsub('5-Mar', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup$Symbol <- gsub('6-Mar', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup$Symbol <- gsub('6-Sep', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup$Symbol <- gsub('7-Mar', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup$Symbol <- gsub('7-Sep', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup$Symbol <- gsub('8-Mar', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup$Symbol <- gsub('8-Sep', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup$Symbol <- gsub('9-Mar', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup$Symbol <- gsub('9-Sep', '', ADNI_rowGroup$Symbol)
ADNI_rowGroup <- ADNI_rowGroup$Symbol

# NCBI1

NCBI1_datET <- Biobase::exprs(NCBI1) # counts data.frame
NCBI1_rowGroup <- Biobase::fData(NCBI1)$`Gene symbol` # gene names (common)
NCBI1_rowID <- row.names(Biobase::fData(NCBI1)) #probes names (unique)

# NCBI2

NCBI2_datET <- Biobase::exprs(NCBI2) # counts data.frame
NCBI2_rowGroup <- Biobase::fData(NCBI2)$`Gene symbol` # gene names (common)
NCBI2_rowID <- row.names(Biobase::fData(NCBI2)) #probes names (unique)

## keep and merge probes into genes

adni_filt <- WGCNA::collapseRows(ADNI_datET, ADNI_rowGroup, ADNI_rowID, 
                                 method="maxRowVariance", connectivityBasedCollapsing=FALSE,
                                 methodFunction=NULL, connectivityPower=8,
                                 selectFewestMissing=TRUE, thresholdCombine=NA)

ncbi1_filt <- WGCNA::collapseRows(NCBI1_datET, NCBI1_rowGroup, NCBI1_rowID, 
                                  method="maxRowVariance", connectivityBasedCollapsing=FALSE,
                                  methodFunction=NULL, connectivityPower=8,
                                  selectFewestMissing=TRUE, thresholdCombine=NA)

ncbi2_filt <- WGCNA::collapseRows(NCBI2_datET, NCBI2_rowGroup, NCBI2_rowID, 
                                  method="maxRowVariance", connectivityBasedCollapsing=FALSE,
                                  methodFunction=NULL, connectivityPower=8,
                                  selectFewestMissing=TRUE, thresholdCombine=NA)


# Make new ExpressionSet

# ADNI

phenotype_adni <- Biobase::pData(ADNI)
feature_adni <- as.data.frame(adni_filt$group2row)
exprs_adni <- adni_filt$datETcollapsed
all(rownames(phenotype_adni)==colnames(exprs_adni))# TRUE
all(rownames(feature_adni)==rownames(exprs_adni)) # TRUE

phenoData_adni <- new("AnnotatedDataFrame", data=phenotype_adni)
featureData_adni <- new("AnnotatedDataFrame", data=feature_adni)

ExpressionSet_adni <- ExpressionSet(assayData = exprs_adni, phenoData = phenoData_adni, featureData =featureData_adni)

# NCBI1


phenotype_ncbi1 <- Biobase::pData(NCBI1)
feature_ncbi1 <- as.data.frame(ncbi1_filt$group2row)
exprs_ncbi1 <- ncbi1_filt$datETcollapsed
all(rownames(phenotype_ncbi1)==colnames(exprs_ncbi1))# TRUE
all(rownames(feature_ncbi1)==rownames(exprs_ncbi1)) # TRUE

phenoData_ncbi1 <- new("AnnotatedDataFrame", data=phenotype_ncbi1)
featureData_ncbi1 <- new("AnnotatedDataFrame", data=feature_ncbi1)

ExpressionSet_ncbi1 <- ExpressionSet(assayData = exprs_ncbi1, phenoData = phenoData_ncbi1, featureData =featureData_ncbi1)

# NCBI2

phenotype_ncbi2 <- Biobase::pData(NCBI2)
feature_ncbi2 <- as.data.frame(ncbi2_filt$group2row)
exprs_ncbi2 <- ncbi2_filt$datETcollapsed
all(rownames(phenotype_ncbi2)==colnames(exprs_ncbi2))# TRUE
all(rownames(feature_ncbi2)==rownames(exprs_ncbi2)) # TRUE

phenoData_ncbi2 <- new("AnnotatedDataFrame", data=phenotype_ncbi2)
featureData_ncbi2 <- new("AnnotatedDataFrame", data=feature_ncbi2)

ExpressionSet_ncbi2 <- ExpressionSet(assayData = exprs_ncbi2, phenoData = phenoData_ncbi2, featureData =featureData_ncbi2)

#### Filtering paper for ADNI due to low intensity and low RIN  1)

ADNI <- ExpressionSet_adni

#Filtering samples 1) 2)

#take the RIN_values from pData
RIN_values <- pData(ADNI)[,c("SubjectID","RIN")]
RIN_values <- transform(RIN_values, RIN = as.character(RIN))
RIN_values <- transform(RIN_values, RIN = as.numeric(RIN))
RIN_values$RIN[is.na(RIN_values$RIN)] <- 0

#Histogram
threshold <- 6.9 # recommended

hist_res <- hist(RIN_values$RIN, 100, col = "cornsilk", freq = FALSE, xlim=c(4.0,max(RIN_values$RIN)),
                 main = "Histogram of the RIN frequency",
                 border = "antiquewhite4",
                 xlab = "RIN")
abline(v = threshold, col = "coral4", lwd = 2)

RIN_filtered <- filter(RIN_values, RIN >= threshold)
ADNI_total <- ADNI[,RIN_filtered$SubjectID]

no_of_samples_1 <- table(pData(ADNI_total)$DX)
no_of_samples_1

# Filtering genes


ADNI_medians <- rowMedians(Biobase::exprs(ADNI_total))
man_threshold <- median(ADNI_medians) #3.886 

hist_res <- hist(ADNI_medians, 100, col = "cornsilk", freq = FALSE, 
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")
abline(v = man_threshold, col = "coral4", lwd = 2)


no_of_samples <- table(pData(ADNI_total)$DX)
no_of_samples

samples_cutoff <- min(no_of_samples)
idx_man_threshold <- apply(Biobase::exprs(ADNI_total), 1,
                           function(x){
                             sum(x > man_threshold) >= samples_cutoff})
idx_man_threshold
table(idx_man_threshold)

# exclude low median intensity probes
ADNI_filtered <- subset(ADNI_total, idx_man_threshold)
ExpressionSet_adni <-ADNI_filtered

## proting results

oligo::boxplot(ExpressionSet_adni, xaxt="n", main = "Normalised ADNI", xlab = "samples", ylab="exptADNI intensity")
oligo::boxplot(ExpressionSet_ncbi1, xaxt="n", main = "Normalised ADNI", xlab = "samples", ylab="exptADNI intensity")
oligo::boxplot(ExpressionSet_ncbi2, xaxt="n", main = "Normalised ADNI", xlab = "samples", ylab="exptADNI intensity")


###################### RUN LIMMA #############################
# ADNI

pheno_ADNI <- as.character(Biobase::pData(ExpressionSet_adni)$DX)
pheno_ADNI

counts <- Biobase::exprs(ExpressionSet_adni)
mm <- model.matrix(~0 + pheno_ADNI)
fit <- lmFit(counts, mm)
head(coef(fit))
contr_Dem_MCI <- makeContrasts(pheno_ADNIDementia - pheno_ADNIMCI, levels = colnames(coef(fit)))
contr_CN_MCI <- makeContrasts(pheno_ADNICN - pheno_ADNIMCI, levels = colnames(coef(fit)))
contr_CN_Dem <- makeContrasts(pheno_ADNICN - pheno_ADNIDementia, levels = colnames(coef(fit)))
tmp_Dem_MCI <- eBayes(contrasts.fit(fit, contr_Dem_MCI))
tmp_CN_MCI <- eBayes(contrasts.fit(fit, contr_CN_MCI))
tmp_CN_Dem <- eBayes(contrasts.fit(fit, contr_CN_Dem))
top.table_Dem_MCI <- limma::topTable(tmp_Dem_MCI, sort.by = "P", n = Inf)
top.table_CN_MCI <- limma::topTable(tmp_CN_MCI, sort.by = "P", n = Inf)
top.table_CN_Dem <- limma::topTable(tmp_CN_Dem, sort.by = "P", n = Inf)

hist(top.table_Dem_MCI$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "AD-MCI", xlab = "p-values")
hist(top.table_CN_MCI$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "CN-MCI", xlab = "p-values")
hist(top.table_CN_Dem$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "CN-AD", xlab = "p-values")

nrow(subset(top.table_Dem_MCI, P.Value < 0.001)) #33
nrow(subset(top.table_CN_MCI, P.Value < 0.001)) #33
nrow(subset(top.table_CN_Dem, P.Value < 0.001)) #88

## Statistically significant differential exptNCBIressed genes
length(which(top.table_Dem_MCI$adj.P.Val < 0.05)) #0
length(which(top.table_CN_MCI$adj.P.Val < 0.05)) #0
length(which(top.table_CN_Dem$adj.P.Val < 0.05)) #0

hist_res <- hist(top.table_Dem_MCI$adj.P.Val, 100, col = "cornsilk", freq = FALSE, xlim=c(0, 1),
                 main = "Histogram of AD - MCI",
                 border = "antiquewhite4",
                 xlab = "P-adj value")
abline(v = 0.05, col = "coral4", lwd = 2)

hist_res <- hist(top.table_CN_MCI$adj.P.Val, 100, col = "cornsilk", freq = FALSE, xlim=c(0, 1),
                 main = "Histogram of CN - MCI",
                 border = "antiquewhite4",
                 xlab = "P-adj value")
abline(v = 0.05, col = "coral4", lwd = 2)

hist_res <- hist(top.table_CN_Dem$adj.P.Val, 100, col = "cornsilk", freq = FALSE, xlim=c(0, 1),
                 main = "Histogram of CN - AD",
                 border = "antiquewhite4",
                 xlab = "P-adj value")
abline(v = 0.05, col = "coral4", lwd = 2)

# NCBI1

ncbi1_samples <- table(pData(ExpressionSet_ncbi1)$characteristics_ch1)
ncbi1_samples

pheno_NCBI1 <- str_replace_all(Biobase::pData(ExpressionSet_ncbi1)$characteristics_ch1,"status: ", "")
pheno_NCBI1 <- str_replace_all(pheno_NCBI1," ", "_")
pheno_NCBI1

counts1 <- Biobase::exprs(ExpressionSet_ncbi1)
mm1 <- model.matrix(~0 + pheno_NCBI1)
fit1 <- lmFit(counts1, mm1)
head(coef(fit1))

contr_AD_MCI1 <- makeContrasts(pheno_NCBI1AD - pheno_NCBI1MCI, levels = colnames(coef(fit1)))
contr_CTL_MCI1 <- makeContrasts(pheno_NCBI1CTL - pheno_NCBI1MCI, levels = colnames(coef(fit1)))
contr_CTL_AD1 <- makeContrasts(pheno_NCBI1CTL - pheno_NCBI1AD, levels = colnames(coef(fit1)))
tmp_AD_MCI1 <- eBayes(contrasts.fit(fit1, contr_AD_MCI1))
tmp_CTL_MCI1 <- eBayes(contrasts.fit(fit1, contr_CTL_MCI1))
tmp_CTL_AD1 <- eBayes(contrasts.fit(fit1, contr_CTL_AD1))
top.table_AD_MCI1 <- limma::topTable(tmp_AD_MCI1, sort.by = "P", n = Inf)
top.table_CTL_MCI1 <- limma::topTable(tmp_CTL_MCI1, sort.by = "P", n = Inf)
top.table_CTL_AD1 <- limma::topTable(tmp_CTL_AD1, sort.by = "P", n = Inf)

hist(top.table_AD_MCI1$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "AD-MCI", xlab = "p-values")
hist(top.table_CTL_MCI1$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "CN-MCI", xlab = "p-values")
hist(top.table_CTL_AD1$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "CN-AD", xlab = "p-values")

nrow(subset(top.table_AD_MCI1, P.Value < 0.001)) #49
nrow(subset(top.table_CTL_MCI1, P.Value < 0.001)) #1035
nrow(subset(top.table_CTL_AD1, P.Value < 0.001)) #1219

## Statistically significant differential exptNCBIressed genes

length(which(top.table_AD_MCI1$adj.P.Val < 0.05)) #0
length(which(top.table_CTL_MCI1$adj.P.Val < 0.05)) #2150
length(which(top.table_CTL_AD1$adj.P.Val < 0.05)) #2051

length(which(abs(top.table_AD_MCI1$logFC) > 1)) #0
length(which(abs(top.table_CTL_MCI1$logFC) > 1)) #0
length(which(abs(top.table_CTL_AD1$logFC) > 1)) #0

## Statistically significant and high differential expressed genes (PADJ < 0.05 & ABSLOG2FOLD > 1)

length(which(top.table_AD_MCI1$adj.P.Val < 0.05 & abs(top.table_AD_MCI1$logFC) > 1)) #0
length(which(top.table_CTL_MCI1$adj.P.Val < 0.05 & abs(top.table_CTL_MCI1$logFC) > 1)) #0
length(which(top.table_CTL_AD1$adj.P.Val < 0.05 & abs(top.table_CTL_AD1$logFC) > 1)) #0

hist_res <- hist(top.table_AD_MCI1$adj.P.Val, 100, col = "cornsilk", freq = FALSE, xlim=c(0, 1),
                 main = "Histogram of AD - MCI",
                 border = "antiquewhite4",
                 xlab = "P-adj value")
abline(v = 0.05, col = "coral4", lwd = 2)

hist_res <- hist(top.table_CTL_MCI1$adj.P.Val, 100, col = "cornsilk", freq = FALSE, xlim=c(0, 1),
                 main = "Histogram of CN - MCI",
                 border = "antiquewhite4",
                 xlab = "P-adj value")
abline(v = 0.05, col = "coral4", lwd = 2)

hist_res <- hist(top.table_CTL_AD1$adj.P.Val, 100, col = "cornsilk", freq = FALSE, xlim=c(0, 1),
                 main = "Histogram of CN - AD",
                 border = "antiquewhite4",
                 xlab = "P-adj value")
abline(v = 0.05, col = "coral4", lwd = 2)

# NCBI 2

ncbi2_samples <- table(pData(ExpressionSet_ncbi2)$characteristics_ch1)
ncbi2_samples

pheno_NCBI2 <- str_replace_all(Biobase::pData(ExpressionSet_ncbi2)$characteristics_ch1,"status: ", "")
pheno_NCBI2 <- str_replace_all(pheno_NCBI2," ", "_")
pheno_NCBI2

counts2 <- Biobase::exprs(ExpressionSet_ncbi2)
mm2 <- model.matrix(~0 + pheno_NCBI2)
fit2 <- lmFit(counts2, mm2)
head(coef(fit2))

contr_AD_MCI2 <- makeContrasts(pheno_NCBI2AD - pheno_NCBI2CTL, levels = colnames(coef(fit2)))
contr_CTL_MCI2 <- makeContrasts(pheno_NCBI2CTL - pheno_NCBI2MCI, levels = colnames(coef(fit2)))
contr_CTL_AD2 <- makeContrasts(pheno_NCBI2CTL - pheno_NCBI2AD, levels = colnames(coef(fit2)))
tmp_AD_MCI2 <- eBayes(contrasts.fit(fit2, contr_AD_MCI2))
tmp_CTL_MCI2 <- eBayes(contrasts.fit(fit2, contr_CTL_MCI2))
tmp_CTL_AD2 <- eBayes(contrasts.fit(fit2, contr_CTL_AD2))
top.table_AD_MCI2 <- limma::topTable(tmp_AD_MCI2, sort.by = "P", n = Inf)
top.table_CTL_MCI2 <- limma::topTable(tmp_CTL_MCI2, sort.by = "P", n = Inf)
top.table_CTL_AD2 <- limma::topTable(tmp_CTL_AD2, sort.by = "P", n = Inf)

hist(top.table_AD_MCI2$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "AD-MCI", xlab = "p-values")
hist(top.table_CTL_MCI2$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "CN-MCI", xlab = "p-values")
hist(top.table_CTL_AD2$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "CN-AD", xlab = "p-values")

nrow(subset(top.table_AD_MCI2, P.Value < 0.001)) #2099
nrow(subset(top.table_CTL_MCI2, P.Value < 0.001)) #3616
nrow(subset(top.table_CTL_AD2, P.Value < 0.001)) #2099

## Statistically significant differential expressed genes (PADJ < 0.05)
length(which(top.table_AD_MCI2$adj.P.Val < 0.05)) #3162
length(which(top.table_CTL_MCI2$adj.P.Val < 0.05)) #5537
length(which(top.table_CTL_AD2$adj.P.Val < 0.05)) #3162

## Statistically significant and high differential expressed genes (PADJ < 0.05 & ABSLOG2FOLD > 1)

length(which(top.table_AD_MCI2$adj.P.Val < 0.05 & abs(top.table_AD_MCI2$logFC) > 1)) #5
length(which(top.table_CTL_MCI2$adj.P.Val < 0.05 & abs(top.table_CTL_MCI2$logFC) > 1)) #9
length(which(top.table_CTL_AD2$adj.P.Val < 0.05 & abs(top.table_CTL_AD2$logFC) > 1)) #5

## HISTOGRAMMS

hist_res <- hist(top.table_AD_MCI2$adj.P.Val, 100, col = "cornsilk", freq = FALSE, xlim=c(0, 1),
                 main = "Histogram of AD - MCI",
                 border = "antiquewhite4",
                 xlab = "P-adj value")
abline(v = 0.05, col = "coral4", lwd = 2)

hist_res <- hist(top.table_CTL_MCI2$adj.P.Val, 100, col = "cornsilk", freq = FALSE, xlim=c(0, 1),
                 main = "Histogram of CN - MCI",
                 border = "antiquewhite4",
                 xlab = "P-adj value")
abline(v = 0.05, col = "coral4", lwd = 2)

hist_res <- hist(top.table_CTL_AD2$adj.P.Val, 100, col = "cornsilk", freq = FALSE, xlim=c(0, 1),
                 main = "Histogram of CN - AD",
                 border = "antiquewhite4",
                 xlab = "P-adj value")
abline(v = 0.05, col = "coral4", lwd = 2)

rm(list=setdiff(ls(), c("ExpressionSet_adni", "ExpressionSet_ncbi1", "ExpressionSet_ncbi2")))

##features
adni_genes <- Biobase::fData(ExpressionSet_adni)
ncbi1_genes <- Biobase::fData(ExpressionSet_ncbi1)
ncbi2_genes <- Biobase::fData(ExpressionSet_ncbi2)
adni_ncbi1_genes <- adni_genes[adni_genes$group %in% ncbi1_genes$group,]
common_genes <- as.data.frame(adni_ncbi1_genes[adni_ncbi1_genes$group %in% ncbi2_genes$group,]) #9643 common genes
common_genes <- new("AnnotatedDataFrame", data=common_genes)

##expression
adni_table <- Biobase::exprs(ExpressionSet_adni)
ncbi1_table <- Biobase::exprs(ExpressionSet_ncbi1)
ncbi2_table <- Biobase::exprs(ExpressionSet_ncbi2)


adni_total <- adni_table[rownames(adni_table) %in% common_genes$group,]
ncbi1_total <- ncbi1_table[rownames(ncbi1_table) %in% common_genes$group,]
ncbi2_total <- ncbi2_table[rownames(ncbi2_table) %in% common_genes$group,]
all(rownames(adni_total)==rownames(ncbi1_total)) # should be TRUE
all(rownames(adni_total)==rownames(ncbi2_total)) # should be TRUE
all(rownames(ncbi1_total)==rownames(ncbi2_total)) # should be TRUE
dim(adni_total) # 9643 431
dim(ncbi1_total) # 9643 388
dim(ncbi2_total) # 9643 329

#pheno

adni_pheno <- Biobase::pData(ExpressionSet_adni)
ncbi1_pheno <- Biobase::pData(ExpressionSet_ncbi1)
ncbi2_pheno <- Biobase::pData(ExpressionSet_ncbi2)

SubjectID <- rownames(adni_pheno)
Condition <- adni_pheno$DX
adni_pheno <- data.frame(SubjectID,Condition)
adni_pheno$Condition <- str_replace_all(adni_pheno$Condition, "Dementia", "AD")
adni_pheno$Condition <- str_replace_all(adni_pheno$Condition, "CN", "CTL")


SubjectID <- rownames(ncbi1_pheno)
Condition <- ncbi1_pheno$characteristics_ch1
ncbi1_pheno <- data.frame(SubjectID,Condition)
rownames(ncbi1_pheno) <- ncbi1_pheno$SubjectID
ncbi1_pheno$Condition <- str_replace_all(ncbi1_pheno$Condition, "status: ", "")
ncbi1_pheno$Condition <- str_replace_all(ncbi1_pheno$Condition, "OTHER", "CTL")
ncbi1_pheno$Condition <- str_replace_all(ncbi1_pheno$Condition, "borderline ", "")
ncbi1_pheno$Condition <- str_replace_all(ncbi1_pheno$Condition, "CTL to AD", "AD")
ncbi1_pheno$Condition <- str_replace_all(ncbi1_pheno$Condition, "MCI to CTL", "MCI")


SubjectID <- rownames(ncbi2_pheno)
Condition <- ncbi2_pheno$characteristics_ch1
ncbi2_pheno <- data.frame(SubjectID,Condition)
rownames(ncbi2_pheno) <- ncbi2_pheno$SubjectID
ncbi2_pheno$Condition <- str_replace_all(ncbi2_pheno$Condition, "status: ", "")

adni_pheno$batch <- 1
ncbi1_pheno$batch <- 2
ncbi2_pheno$batch <- 3


adni_pheno <- new("AnnotatedDataFrame", data=adni_pheno)
ncbi1_pheno <- new("AnnotatedDataFrame", data=ncbi1_pheno)
ncbi2_pheno <- new("AnnotatedDataFrame", data=ncbi2_pheno)


ExpressionSet_adni <- ExpressionSet(assayData = adni_total, phenoData = adni_pheno, featureData = common_genes)
ExpressionSet_ncbi1 <- ExpressionSet(assayData = ncbi1_total, phenoData = ncbi1_pheno, featureData = common_genes)
ExpressionSet_ncbi2 <- ExpressionSet(assayData = ncbi2_total, phenoData = ncbi2_pheno, featureData = common_genes)


###########################################################################################
#combine data & batch correction
## COMBINE
part <- combineTwoExpressionSet(ExpressionSet_ncbi1,ExpressionSet_ncbi2)
common <- combineTwoExpressionSet(part,ExpressionSet_adni)
## BATCH CORRECTION
#pre-filtering
oligo::boxplot(common, xaxt="n", main = "Non-Normalized", xlab="samples", ylab="Features Intensity")

#process

pheno <- Biobase::pData(common)
feature <- Biobase::fData(common)
edata <- Biobase::exprs(common)
batch<-pheno$batch
modcombat<-model.matrix(~1, data=pheno)
combat_mydata= ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)


phenodata <- new("AnnotatedDataFrame", data=pheno)
feature <- new("AnnotatedDataFrame", data=feature)


TotalExpressionSet <-ExpressionSet(assayData = combat_mydata, phenoData = phenodata, featureData =feature)

#Pro-filtering
oligo::boxplot(TotalExpressionSet, xaxt="n", main = "Non-Normalized", xlab="samples", ylab="Features Intensity")

table(pData(TotalExpressionSet)$Condition)

saveRDS(TotalExpressionSet, 'Total_ExpressionSet.Rds')

################### END OF FIRST PART ##################################################
## START TIDYMODELS AND CARET NON-LINEAR METHODS
## INPUT ExpressionSet that contain the combined dataset (ADNI-NCBI1-NCBI2) COMMON GENES, SAME PHENOTYPES NAME
ExpressionSet <- readRDS("Total_ExpressionSet.Rds")
## Non-linear method

##Separating into CN-AD, CN-MCI, AD-MCI

filter_MCI <- colnames(ExpressionSet)[ExpressionSet@phenoData@data$Condition != "MCI"] #keep MCI samples
filter_CTL <- colnames(ExpressionSet)[ExpressionSet@phenoData@data$Condition != "CTL"] #keep CTL samples
filter_AD <- colnames(ExpressionSet)[ExpressionSet@phenoData@data$Condition != "AD"] #keep AD samples

ExpressionSet_CTL_AD <- ExpressionSet[,filter_MCI] # make an ExpressionSet only with CTL_AD
ExpressionSet_MCI_AD <- ExpressionSet[,filter_CTL] # make an ExpressionSet only with MCI_AD
ExpressionSet_CTL_MCI <- ExpressionSet[,filter_AD] # make an ExpressionSet only with CTL_MCI

### SPECIFY THE MODEL ###

svm_model <- 
  # specify that the model is svm use either _poly or _rbf
  svm_rbf() %>%
  # select the engine/package that underlies the model
  set_engine("kernlab") %>%
  # choose either the continuous regression or binary classification mode
  set_mode("classification")

rf_model <- 
  # specify that the model is a random forest
  rand_forest() %>%
  # select the engine/package that underlies the model
  set_engine("ranger", importance = "impurity") %>%
  # choose either the continuous regression or binary classification mode
  set_mode("classification") 
#Pre-processing

##Pre-processing file -- create for each comparison (CTL-AD, CTL-MCI, MCI-AD) the exprs table
#that contains the expression data (outcomes) with a last column of Phenotype (predictions)

#take outcomes table as reversed to have genes at columns
a_CTL_AD <- as.data.frame(t(Biobase::exprs(ExpressionSet_CTL_AD)))
a_MCI_AD <- as.data.frame(t(Biobase::exprs(ExpressionSet_MCI_AD)))
a_CTL_MCI <- as.data.frame(t(Biobase::exprs(ExpressionSet_CTL_MCI)))

# take phenotypes as vector from expressionSet
X_CTL_AD <- Biobase::pData(ExpressionSet_CTL_AD)$Condition #phenotype
X_MCI_AD <- Biobase::pData(ExpressionSet_MCI_AD)$Condition #phenotype
X_CTL_MCI <- Biobase::pData(ExpressionSet_CTL_MCI)$Condition #phenotype

#Take the phenotype from ExpressionSet and put it as predictors in last columns
a_CTL_AD$Pheno <- as.factor(X_CTL_AD) #input phenotype into expression
a_MCI_AD$Pheno <- as.factor(X_MCI_AD) #input phenotype into expression
a_CTL_MCI$Pheno <- as.factor(X_CTL_MCI) #input phenotype into expression

#make.names the columns, otherwise does not run in analysis
colnames(a_CTL_AD) <- make.names(as.vector(colnames(a_CTL_AD)), unique = FALSE, allow_ = TRUE)
colnames(a_MCI_AD) <- make.names(as.vector(colnames(a_MCI_AD)), unique = FALSE, allow_ = TRUE)
colnames(a_CTL_MCI) <- make.names(as.vector(colnames(a_CTL_MCI)), unique = FALSE, allow_ = TRUE)

## Spliting data into train-test
# I used initial_split() factor to separate them into train and test parts
#For CTL_AD

set.seed(1333)

##split partitions

a_CTL_AD_5folds <- caret::createFolds(a_CTL_AD$Pheno, k=5)
a_CTL_MCI_5folds <- caret::createFolds(a_CTL_MCI$Pheno, k=5)
a_MCI_AD_5folds <- caret::createFolds(a_MCI_AD$Pheno, k=5)

#This step occurs because of the file format needed in analysis

a_CTL_AD_splited <- initial_split(a_CTL_AD)
a_CTL_MCI_splited <- initial_split(a_CTL_MCI)
a_MCI_AD_splited <- initial_split(a_MCI_AD)


#make lists to save results

test_performance_rf_CTL_AD_5folds <- vector("list")
test_performance_rf_MCI_AD_5folds <- vector("list")
test_performance_rf_CTL_MCI_5folds <- vector("list")

test_performance_svm_CTL_AD_5folds <- vector("list")
test_performance_svm_MCI_AD_5folds <- vector("list")
test_performance_svm_CTL_MCI_5folds <- vector("list")

test_predictions_rf_CTL_AD_5folds <- vector("list")
test_predictions_rf_CTL_MCI_5folds <- vector("list")
test_predictions_rf_MCI_AD_5folds <- vector("list")

test_predictions_svm_CTL_AD_5folds <- vector("list")
test_predictions_svm_MCI_AD_5folds <- vector("list") 
test_predictions_svm_CTL_MCI_5folds <- vector("list")

final_model_rf_CTL_AD_5folds <- vector("list")
final_model_rf_MCI_AD_5folds <- vector("list")
final_model_rf_CTL_MCI_5folds <- vector("list")

final_model_svm_CTL_AD_5folds <- vector("list")
final_model_svm_MCI_AD_5folds <- vector("list")
final_model_svm_CTL_MCI_5folds <- vector("list")

myROC_rf_CTL_AD_5folds <- vector("list")
myROC_rf_MCI_AD_5folds <- vector("list")
myROC_rf_CTL_MCI_5folds <- vector("list")

myROC_svm_CTL_AD_5folds <- vector("list")
myROC_svm_MCI_AD_5folds <- vector("list")
myROC_svm_CTL_MCI_5folds <- vector("list")

obj_rf_CTL_AD_5folds <- vector("list")
obj_rf_MCI_AD_5folds <- vector("list")
obj_rf_CTL_MCI_5folds <- vector("list")

obj_svm_CTL_AD_5folds <- vector("list")
obj_svm_MCI_AD_5folds <- vector("list")
obj_svm_CTL_MCI_5folds <- vector("list")

for (i in seq_along(a_CTL_AD_5folds)) {
  my_training_index_CTL_AD <- a_CTL_AD_5folds[[i]]
  a_CTL_AD_train <- a_CTL_AD[-my_training_index_CTL_AD,]
  a_CTL_AD_test <- a_CTL_AD[my_training_index_CTL_AD,]
  all_integer_CTL_AD <-sort(c(a_CTL_AD_5folds[[1]],a_CTL_AD_5folds[[2]],
                              a_CTL_AD_5folds[[3]],a_CTL_AD_5folds[[4]],
                              a_CTL_AD_5folds[[5]]))
  integer_CTL_AD <- all_integer_CTL_AD[-a_CTL_AD_5folds[[i]]]
  a_CTL_AD_splited$in_id <- integer_CTL_AD
  
  
  my_training_index_CTL_MCI <- a_CTL_MCI_5folds[[i]]
  a_CTL_MCI_train <- a_CTL_MCI[-my_training_index_CTL_MCI,]
  a_CTL_MCI_test <- a_CTL_MCI[my_training_index_CTL_MCI,]
  
  all_integer_CTL_MCI <-sort(c(a_CTL_MCI_5folds[[1]],a_CTL_MCI_5folds[[2]],
                               a_CTL_MCI_5folds[[3]],a_CTL_MCI_5folds[[4]],
                               a_CTL_MCI_5folds[[5]]))
  integer_CTL_MCI <- all_integer_CTL_MCI[-a_CTL_MCI_5folds[[i]]]
  a_CTL_MCI_splited$in_id <- integer_CTL_MCI
  
  my_training_index_MCI_AD <- a_MCI_AD_5folds[[i]]
  a_MCI_AD_train <- a_MCI_AD[-my_training_index_MCI_AD,]
  a_MCI_AD_test <- a_MCI_AD[my_training_index_MCI_AD,]
  
  all_integer_MCI_AD <-sort(c(a_MCI_AD_5folds[[1]],a_MCI_AD_5folds[[2]],
                              a_MCI_AD_5folds[[3]],a_MCI_AD_5folds[[4]],
                              a_MCI_AD_5folds[[5]]))
  integer_MCI_AD <- all_integer_MCI_AD[-a_MCI_AD_5folds[[i]]]
  a_MCI_AD_splited$in_id <- integer_MCI_AD
  
  ## define the recipe in each comparison FOR RANDOM FOREST
  
  #CTL_AD
  a_CTL_AD_rf_recipe <-
    # which consists of the formula (outcome ~ predictors)
    recipe(Pheno~ ., # . means that take all outcomes
           data = a_CTL_AD_train) %>%
    #remove variables that are same (or almost same) for all numeric (our outcomes)
    step_nzv(all_numeric())
  #MCI_AD
  a_MCI_AD_rf_recipe <-
    # which consists of the formula (outcome ~ predictors)
    recipe(Pheno~ ., # . means that take all outcomes
           data = a_MCI_AD_train) %>%
    #remove variables that are same (or almost same) for all numeric (our outcomes)
    step_nzv(all_numeric())
  #CTL_MCI
  a_CTL_MCI_rf_recipe <-
    # which consists of the formula (outcome ~ predictors)
    recipe(Pheno~ ., # . means that take all outcomes
           data = a_CTL_MCI_train) %>%
    #remove variables that are same (or almost same) for all numeric (our outcomes)
    step_nzv(all_numeric())
  
  ## define the recipe in each comparison FOR SUPPORT VECTOR MACHINE
  a_CTL_AD_svm_recipe <-
    # which consists of the formula (outcome ~ predictors)
    recipe(Pheno ~ ., # . means that take all columns except Pheno(phenotype that I input above)
           data = a_CTL_AD_train) %>%
    # normalize numeric data to be within a pre-defined range of values
    #I remember that you told me it is necessary for SVM
    step_range(all_numeric()) %>%
    #remove variables that are same (or almost same) for all numeric (our outcomes)
    step_nzv(all_numeric())
  a_MCI_AD_svm_recipe <-
    # which consists of the formula (outcome ~ predictors)
    recipe(Pheno ~ ., # . means that take all columns except Pheno(phenotype that I input above)
           data = a_MCI_AD_train) %>%
    #step_range : normalize numeric data to be within a pre-defined range of values
    step_range(all_numeric()) %>%
    #remove variables that are same (or almost same) for all numeric (our outcomes)
    step_nzv(all_numeric())
  a_CTL_MCI_svm_recipe <-
    # which consists of the formula (outcome ~ predictors)
    recipe(Pheno ~ ., # . means that take all columns except Pheno(phenotype that I input above)
           data = a_CTL_MCI_train) %>%
    # normalize numeric data to be within a pre-defined range of values 
    step_range(all_numeric()) %>%
    #remove variables that are same (or almost same) for all numeric (our outcomes)
    step_nzv(all_numeric())
  
  
  ##### set the workflow for each combarison for RANDOM FOREST
  # set the workflow
  rf_workflow_CTL_AD <- workflow() %>%
    # add the recipe
    add_recipe(a_CTL_AD_rf_recipe) %>%
    # add the model
    add_model(rf_model)
  # set the workflow
  rf_workflow_MCI_AD <- workflow() %>%
    # add the recipe
    add_recipe(a_MCI_AD_rf_recipe) %>%
    # add the model
    add_model(rf_model)
  # set the workflow
  rf_workflow_CTL_MCI <- workflow() %>%
    # add the recipe
    add_recipe(a_CTL_MCI_rf_recipe) %>%
    # add the model
    add_model(rf_model)
  
  ##### set the workflow for each combarison for SUPPORT VECTOR MACHINE
  # set the workflow
  svm_workflow_CTL_AD <- workflow() %>%
    # add the recipe
    add_recipe(a_CTL_AD_svm_recipe) %>%
    # add the model
    add_model(svm_model)
  # set the workflow
  svm_workflow_MCI_AD <- workflow() %>%
    # add the recipe
    add_recipe(a_MCI_AD_svm_recipe) %>%
    # add the model
    add_model(svm_model)
  # set the workflow
  svm_workflow_CTL_MCI <- workflow() %>%
    # add the recipe
    add_recipe(a_CTL_MCI_svm_recipe) %>%
    # add the model
    add_model(svm_model)
  #FIT RESAMPLES BOTH RANDOM FOREST AND SUPPORT VECTOR MACHINE
  #make either a bootstraping or simple cross validation file
  #for CTL_AD
  a_CTL_AD_boot <- bootstraps(a_CTL_AD_train, times= 5, apparent = FALSE)
  #for MCI_AD
  a_MCI_AD_boot <- bootstraps(a_MCI_AD_train, times= 5, apparent = FALSE)
  #for CTL_MCI
  a_CTL_MCI_boot <- bootstraps(a_CTL_MCI_train, times= 5, apparent = FALSE) 
  
  # FIT FOR RANDOM FOREST
  rf_results_CTL_AD <- rf_workflow_CTL_AD %>%
    #use either cross_Validation or bootstrap
    fit_resamples(resamples = a_CTL_AD_boot,
                  #use only roc_auc
                  metrics = metric_set(roc_auc))
  
  rf_results_MCI_AD <- rf_workflow_MCI_AD %>%
    #use either cross_Validation or bootstrap
    fit_resamples(resamples = a_MCI_AD_boot,
                  #use only roc_auc
                  metrics = metric_set(roc_auc))
  
  rf_results_CTL_MCI <- rf_workflow_CTL_MCI %>%
    #use either cross_Validation or bootstrap
    fit_resamples(resamples = a_CTL_MCI_boot,
                  #use only roc_auc
                  metrics = metric_set(roc_auc))
  
  # FIT FOR SUPPORT VECTOR MACHINE
  svm_results_CTL_AD <- svm_workflow_CTL_AD %>%
    #use either cross_Validation or bootstrap
    fit_resamples(resamples = a_CTL_AD_boot,
                  #use only roc_auc
                  metrics = metric_set(roc_auc))
  
  svm_results_MCI_AD <- svm_workflow_MCI_AD %>%
    fit_resamples(resamples = a_MCI_AD_boot,
                  #use only roc_auc
                  metrics = metric_set(roc_auc))
  
  svm_results_CTL_MCI <- svm_workflow_CTL_MCI %>%
    fit_resamples(resamples = a_CTL_MCI_boot,
                  #use only roc_auc
                  metrics = metric_set(roc_auc))
  
  # WANT TO TAKE THE BEST roc_auc BUT I HAVE A TRIBBLE 0 X 1 output
  param_final_rf_CTL_AD <- rf_results_CTL_AD %>%
    #from fit_resamples output take the best roc_auc
    select_best(metric = "roc_auc")
  param_final_rf_MCI_AD <- rf_results_MCI_AD %>%
    #from fit_resamples output take the best roc_auc
    select_best(metric = "roc_auc")
  param_final_rf_CTL_MCI <- rf_results_CTL_MCI %>%
    #from fit_resamples output take the best roc_auc
    select_best(metric = "roc_auc")
  
  
  param_final_svm_CTL_AD <- svm_results_CTL_AD %>%
    #from fit_resamples output take the best roc_auc
    select_best(metric = "roc_auc")
  param_final_svm_MCI_AD <- svm_results_MCI_AD %>%
    #from fit_resamples output take the best roc_auc
    select_best(metric = "roc_auc")
  param_final_svm_CTL_MCI <- svm_results_CTL_MCI %>%
    #from fit_resamples output take the best roc_auc
    select_best(metric = "roc_auc")
  ## 
  
  # In case we have a significant roc_auc from above we input it into our workflow
  rf_workflow_CTL_AD <- rf_workflow_CTL_AD %>%
    finalize_workflow(param_final_rf_CTL_AD)
  rf_workflow_MCI_AD <- rf_workflow_MCI_AD %>%
    finalize_workflow(param_final_rf_MCI_AD)
  rf_workflow_CTL_MCI <- rf_workflow_CTL_MCI %>%
    finalize_workflow(param_final_rf_CTL_MCI)
  
  svm_workflow_CTL_AD <- svm_workflow_CTL_AD %>%
    finalize_workflow(param_final_svm_CTL_AD)
  svm_workflow_MCI_AD <- svm_workflow_MCI_AD %>%
    finalize_workflow(param_final_svm_MCI_AD)
  svm_workflow_CTL_MCI <- svm_workflow_CTL_MCI %>%
    finalize_workflow(param_final_svm_CTL_MCI)
  
  #Evaluate the model on the test set on RANDOM FOREST
  
  rf_fit_CTL_AD <- rf_workflow_CTL_AD %>%
    # fit on the training set and evaluate on test set
    last_fit(a_CTL_AD_splited)
  rf_fit_MCI_AD <- rf_workflow_MCI_AD %>%
    # fit on the training set and evaluate on test set
    last_fit(a_MCI_AD_splited)
  rf_fit_CTL_MCI <- rf_workflow_CTL_MCI %>%
    # fit on the training set and evaluate on test set
    last_fit(a_CTL_MCI_splited)
  
  #Evaluate the model on the test set on SUPPORT VECTOR MACHINE
  
  svm_fit_CTL_AD <- svm_workflow_CTL_AD %>%
    # fit on the training set and evaluate on test set
    last_fit(a_CTL_AD_splited)
  svm_fit_MCI_AD <- svm_workflow_MCI_AD %>%
    # fit on the training set and evaluate on test set
    last_fit(a_MCI_AD_splited)
  svm_fit_CTL_MCI <- svm_workflow_CTL_MCI %>%
    # fit on the training set and evaluate on test set
    last_fit(a_CTL_MCI_splited)
  
  ## Collect metrics from above rf
  test_performance_rf_CTL_AD <- rf_fit_CTL_AD %>% collect_metrics()
  test_performance_rf_MCI_AD <- rf_fit_MCI_AD %>% collect_metrics()
  test_performance_rf_CTL_MCI <- rf_fit_CTL_MCI %>% collect_metrics()
  
  ## Collect metrics from above svm
  test_performance_svm_CTL_AD <- svm_fit_CTL_AD %>% collect_metrics()
  test_performance_svm_MCI_AD <- svm_fit_MCI_AD %>% collect_metrics()
  test_performance_svm_CTL_MCI <- svm_fit_CTL_MCI %>% collect_metrics()
  
  ## Save them into list for saving
  test_performance_rf_CTL_AD_5folds[[i]] <- test_performance_rf_CTL_AD
  test_performance_rf_MCI_AD_5folds[[i]] <- test_performance_rf_MCI_AD
  test_performance_rf_CTL_MCI_5folds[[i]] <- test_performance_rf_CTL_MCI
  
  test_performance_svm_CTL_AD_5folds[[i]] <- test_performance_svm_CTL_AD
  test_performance_svm_MCI_AD_5folds[[i]] <- test_performance_svm_MCI_AD
  test_performance_svm_CTL_MCI_5folds[[i]] <- test_performance_svm_CTL_MCI
  
  ## extract prediction metrics from above rf
  # generate predictions from the test set
  test_predictions_rf_CTL_AD <- rf_fit_CTL_AD %>% collect_predictions()
  test_predictions_rf_MCI_AD <- rf_fit_MCI_AD %>% collect_predictions()
  test_predictions_rf_CTL_MCI <- rf_fit_CTL_MCI %>% collect_predictions()
  
  ## extract prediction metrics from above svm
  # generate predictions from the test set
  test_predictions_svm_CTL_AD <- svm_fit_CTL_AD %>% collect_predictions()
  test_predictions_svm_MCI_AD <- svm_fit_MCI_AD %>% collect_predictions()
  test_predictions_svm_CTL_MCI <- svm_fit_CTL_MCI %>% collect_predictions()
  
  # generate a confusion matrix for rf
  test_predictions_rf_CTL_AD %>% 
    conf_mat(truth = Pheno, estimate = .pred_class)
  test_predictions_rf_MCI_AD %>% 
    conf_mat(truth = Pheno, estimate = .pred_class)
  test_predictions_rf_CTL_MCI %>% 
    conf_mat(truth = Pheno, estimate = .pred_class)
  
  # generate a confusion matrix for svm
  test_predictions_svm_CTL_AD %>% 
    conf_mat(truth = Pheno, estimate = .pred_class)
  test_predictions_svm_MCI_AD %>% 
    conf_mat(truth = Pheno, estimate = .pred_class)
  test_predictions_svm_CTL_MCI %>% 
    conf_mat(truth = Pheno, estimate = .pred_class)
  
  ## save them into list for saving
  test_predictions_rf_CTL_AD_5folds[[i]] <- test_predictions_rf_CTL_AD
  test_predictions_rf_CTL_MCI_5folds[[i]] <- test_predictions_rf_MCI_AD
  test_predictions_rf_MCI_AD_5folds[[i]] <- test_predictions_rf_CTL_MCI
  
  test_predictions_svm_CTL_AD_5folds[[i]] <- test_predictions_svm_CTL_AD
  test_predictions_svm_CTL_MCI_5folds[[i]] <- test_predictions_svm_MCI_AD
  test_predictions_svm_MCI_AD_5folds[[i]] <- test_predictions_svm_CTL_MCI
  
  ## Prepare them for plotting Roc
  ## FOR RANDOM FOREST
  
  pred_results_rf_CTL_AD <- test_predictions_rf_CTL_AD
  pred_results_rf_CTL_AD$Orig_Class <- ifelse(pred_results_rf_CTL_AD$Pheno == "AD",1,0) %>% as.factor()
  confusionMatrix(data = pred_results_rf_CTL_AD$.pred_class, reference = pred_results_rf_CTL_AD$Pheno, positive = "AD", mode = "everything")
  myROC_rf_CTL_AD <- pROC::roc(response = pred_results_rf_CTL_AD$Orig_Class, predictor = pred_results_rf_CTL_AD$.pred_AD, auc = TRUE, ci = TRUE)
  
  pred_results_rf_MCI_AD <- test_predictions_rf_MCI_AD
  pred_results_rf_MCI_AD$Orig_Class <- ifelse(pred_results_rf_MCI_AD$Pheno == "AD",1,0) %>% as.factor()
  confusionMatrix(data = pred_results_rf_MCI_AD$.pred_class, reference = pred_results_rf_MCI_AD$Pheno, positive = "AD", mode = "everything")
  myROC_rf_MCI_AD <- pROC::roc(response = pred_results_rf_MCI_AD$Orig_Class, predictor = pred_results_rf_MCI_AD$.pred_AD, auc = TRUE, ci = TRUE)
  
  pred_results_rf_CTL_MCI <- test_predictions_rf_CTL_MCI
  pred_results_rf_CTL_MCI$Orig_Class <- ifelse(pred_results_rf_CTL_MCI$Pheno == "MCI",1,0) %>% as.factor()
  confusionMatrix(data = pred_results_rf_CTL_MCI$.pred_class, reference = pred_results_rf_CTL_MCI$Pheno, positive = "MCI", mode = "everything")
  myROC_rf_CTL_MCI <- pROC::roc(response = pred_results_rf_CTL_MCI$Orig_Class, predictor = pred_results_rf_CTL_MCI$.pred_MCI, auc = TRUE, ci = TRUE)
  
  ## ROC PLOTTING FOR SUPPORT VECTOR MACHINE 
  
  pred_results_svm_CTL_AD <- test_predictions_svm_CTL_AD
  pred_results_svm_CTL_AD$Orig_Class <- ifelse(pred_results_svm_CTL_AD$Pheno == "AD",1,0) %>% as.factor()
  confusionMatrix(data = pred_results_svm_CTL_AD$.pred_class, reference = pred_results_svm_CTL_AD$Pheno, positive = "AD", mode = "everything")
  myROC_svm_CTL_AD <- pROC::roc(response = pred_results_svm_CTL_AD$Orig_Class, predictor = pred_results_svm_CTL_AD$.pred_AD, auc = TRUE, ci = TRUE)
  
  pred_results_svm_MCI_AD <- test_predictions_svm_MCI_AD
  pred_results_svm_MCI_AD$Orig_Class <- ifelse(pred_results_svm_MCI_AD$Pheno == "AD",1,0) %>% as.factor()
  confusionMatrix(data = pred_results_svm_MCI_AD$.pred_class, reference = pred_results_svm_MCI_AD$Pheno, positive = "AD", mode = "everything")
  myROC_svm_MCI_AD <- pROC::roc(response = pred_results_svm_MCI_AD$Orig_Class, predictor = pred_results_svm_MCI_AD$.pred_AD, auc = TRUE, ci = TRUE)
  
  pred_results_svm_CTL_MCI <- test_predictions_svm_CTL_MCI
  pred_results_svm_CTL_MCI$Orig_Class <- ifelse(pred_results_svm_CTL_MCI$Pheno == "MCI",1,0) %>% as.factor()
  confusionMatrix(data = pred_results_svm_CTL_MCI$.pred_class, reference = pred_results_svm_CTL_MCI$Pheno, positive = "MCI", mode = "everything")
  myROC_svm_CTL_MCI <- pROC::roc(response = pred_results_svm_CTL_MCI$Orig_Class, predictor = pred_results_svm_CTL_MCI$.pred_MCI, auc = TRUE, ci = TRUE)
  
  myROC_rf_CTL_AD_5folds[[i]] <- myROC_rf_CTL_AD
  myROC_rf_CTL_MCI_5folds[[i]] <- myROC_rf_CTL_MCI
  myROC_rf_MCI_AD_5folds[[i]] <- myROC_rf_MCI_AD
  
  myROC_svm_CTL_AD_5folds[[i]] <- myROC_svm_CTL_AD
  myROC_svm_CTL_MCI_5folds[[i]] <- myROC_svm_CTL_MCI
  myROC_svm_MCI_AD_5folds[[i]] <- myROC_svm_MCI_AD
  
  ####Fitting and using your final model
  #For RANDOM FOREST
  
  final_model_rf_CTL_AD <- fit(rf_workflow_CTL_AD, a_CTL_AD)
  final_model_rf_MCI_AD <- fit(rf_workflow_MCI_AD, a_MCI_AD)
  final_model_rf_CTL_MCI <- fit(rf_workflow_CTL_MCI, a_CTL_MCI)
  
  #For SUPPORT VECTOR MACHINE
  final_model_svm_CTL_AD <- fit(svm_workflow_CTL_AD, a_CTL_AD)
  final_model_svm_MCI_AD <- fit(svm_workflow_MCI_AD, a_MCI_AD)
  final_model_svm_CTL_MCI <- fit(svm_workflow_CTL_MCI, a_CTL_MCI)
  
  final_model_rf_CTL_AD_5folds[[i]] <- final_model_rf_CTL_AD
  final_model_rf_MCI_AD_5folds[[i]] <- final_model_rf_MCI_AD
  final_model_rf_CTL_MCI_5folds[[i]] <- final_model_rf_CTL_MCI
  
  final_model_svm_CTL_AD_5folds[[i]] <- final_model_svm_CTL_AD
  final_model_svm_MCI_AD_5folds[[i]] <- final_model_svm_MCI_AD
  final_model_svm_CTL_MCI_5folds[[i]] <- final_model_svm_CTL_MCI
  
  ## PREPARE FOR IMPORTANT VARIANTS
  
  obj_rf_CTL_AD <- pull_workflow_fit(final_model_rf_CTL_AD)$fit
  obj_rf_MCI_AD <- pull_workflow_fit(final_model_rf_MCI_AD)$fit
  obj_rf_CTL_MCI <- pull_workflow_fit(final_model_rf_CTL_MCI)$fit
  
  obj_svm_CTL_AD <- pull_workflow_fit(final_model_svm_CTL_AD)$fit
  obj_svm_MCI_AD <- pull_workflow_fit(final_model_svm_MCI_AD)$fit
  obj_svm_CTL_MCI <- pull_workflow_fit(final_model_svm_CTL_MCI)$fit
  
  obj_rf_CTL_AD_5folds[[i]] <- obj_rf_CTL_AD
  obj_rf_MCI_AD_5folds[[i]] <- obj_rf_MCI_AD
  obj_rf_CTL_MCI_5folds[[i]] <- obj_rf_CTL_MCI
  
  obj_svm_CTL_AD_5folds[[i]] <- obj_svm_CTL_AD
  obj_svm_MCI_AD_5folds[[i]] <- obj_svm_MCI_AD
  obj_svm_CTL_MCI_5folds[[i]] <- obj_svm_CTL_MCI
  
  print(i)
  if (i == 5) cat("Done!")
}

## PLOTING ROC_AUC



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

## RUN CARET FOR SVM AGAIN

## INPUT ExpressionSet that contain the combined dataset (ADNI-NCBI1-NCBI2) COMMON GENES, SAME PHENOTYPES NAME
ExpressionSet <- readRDS("~/transcriptomics/Total_ExpressionSet.Rds")
## Non-linear method

##Separating into CN-AD, CN-MCI, AD-MCI

filter_MCI <- colnames(ExpressionSet)[ExpressionSet@phenoData@data$Condition != "MCI"] #keep MCI samples
filter_CTL <- colnames(ExpressionSet)[ExpressionSet@phenoData@data$Condition != "CTL"] #keep CTL samples
filter_AD <- colnames(ExpressionSet)[ExpressionSet@phenoData@data$Condition != "AD"] #keep AD samples

ExpressionSet_CTL_AD <- ExpressionSet[,filter_MCI] # make an ExpressionSet only with CTL_AD
ExpressionSet_MCI_AD <- ExpressionSet[,filter_CTL] # make an ExpressionSet only with MCI_AD
ExpressionSet_CTL_MCI <- ExpressionSet[,filter_AD] # make an ExpressionSet only with CTL_MCI

rm(list=setdiff(ls(), c("ExpressionSet_CTL_AD", "ExpressionSet_MCI_AD", "ExpressionSet_CTL_MCI")))

#take outcomes table as reversed to have genes at columns
a_CTL_AD <- as.data.frame(t(Biobase::exprs(ExpressionSet_CTL_AD)))
a_MCI_AD <- as.data.frame(t(Biobase::exprs(ExpressionSet_MCI_AD)))
a_CTL_MCI <- as.data.frame(t(Biobase::exprs(ExpressionSet_CTL_MCI)))

# take phenotypes as vector from expressionSet
X_CTL_AD <- Biobase::pData(ExpressionSet_CTL_AD)$Condition #phenotype
X_MCI_AD <- Biobase::pData(ExpressionSet_MCI_AD)$Condition #phenotype
X_CTL_MCI <- Biobase::pData(ExpressionSet_CTL_MCI)$Condition #phenotype

#Take the phenotype from ExpressionSet and put it as predictors in last columns
a_CTL_AD$Pheno <- as.factor(X_CTL_AD) #input phenotype into expression
a_MCI_AD$Pheno <- as.factor(X_MCI_AD) #input phenotype into expression
a_CTL_MCI$Pheno <- as.factor(X_CTL_MCI) #input phenotype into expression

#make.names the columns, otherwise does not run in analysis
colnames(a_CTL_AD) <- make.names(as.vector(colnames(a_CTL_AD)), unique = FALSE, allow_ = TRUE)
colnames(a_MCI_AD) <- make.names(as.vector(colnames(a_MCI_AD)), unique = FALSE, allow_ = TRUE)
colnames(a_CTL_MCI) <- make.names(as.vector(colnames(a_CTL_MCI)), unique = FALSE, allow_ = TRUE)

## Spliting data into train-test
# I used initial_split() factor to separate them into train and test parts
#For CTL_AD

obj_svm_CTL_AD_5folds <- vector("list")
obj_svm_MCI_AD_5folds <- vector("list")
obj_svm_CTL_MCI_5folds <- vector("list")

final_model_svm_CTL_AD_5folds <- vector("list")
final_model_svm_MCI_AD_5folds <- vector("list")
final_model_svm_CTL_MCI_5folds <- vector("list")

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
  a_CTL_MCI_train <- a_CTL_MCI[!my_training_index_CTL_MCI,]
  a_CTL_MCI_test <- a_CTL_MCI[my_training_index_CTL_MCI,]
  
  my_training_index_MCI_AD <- a_MCI_AD_5folds[[i]]
  a_MCI_AD_train <- a_MCI_AD[!my_training_index_MCI_AD,]
  a_MCI_AD_test <- a_MCI_AD[my_training_index_MCI_AD,]
  
  #fit models and save results
  
  modFit_CTL_AD <- caret::train(Pheno ~ ., data = a_CTL_AD_train,
                                method = "svmRadial", ## svmLinear, etc
                                preProcess = c("range", "nzv"),
                                trControl = trainControl(method = "boot", number = 10))
  
  modFit_CTL_MCI <- caret::train(Pheno ~ ., data = a_CTL_MCI_train,
                                 method = "svmRadial", ## svmLinear, etc
                                 preProcess = c("range", "nzv"),
                                 trControl = trainControl(method = "boot", number = 10))
  
  modFit_MCI_AD <- caret::train(Pheno ~ ., data = a_MCI_AD_train,
                                method = "svmRadial", ## svmLinear, etc
                                preProcess = c("range", "nzv"),
                                trControl = trainControl(method = "boot", number = 10))
  
  myVarImp_CTL_AD <- varImp(modFit_CTL_AD, scale = TRUE)
  myVarImp_CTL_MCI <- varImp(modFit_CTL_MCI, scale = TRUE)
  myVarImp_MCI_AD <- varImp(modFit_MCI_AD, scale = TRUE)
  
  
  final_model_svm_CTL_AD_5folds[[i]] <- modFit_CTL_AD
  final_model_svm_MCI_AD_5folds[[i]] <- modFit_MCI_AD
  final_model_svm_CTL_MCI_5folds[[i]] <- modFit_CTL_MCI
  
  obj_svm_CTL_AD_5folds[[i]] <- myVarImp_CTL_AD
  obj_svm_MCI_AD_5folds[[i]] <- myVarImp_MCI_AD
  obj_svm_CTL_MCI_5folds[[i]] <- myVarImp_CTL_MCI
  
  print{i}
}

############## END OF PART 2################################################################
####### START PART 3 ###################################################################
## FIND SIGNIFICANT VARIANTS 
# RANDOM FOREST

load("rf_obj.RData")

total_rf_important <- vector("list")

#CTL_AD
rf_CTL_AD_1 <- data.frame(obj_rf_CTL_AD_5folds[[1]]$variable.importance)
rf_CTL_AD_2 <- data.frame(obj_rf_CTL_AD_5folds[[2]]$variable.importance)
rf_CTL_AD_3 <- data.frame(obj_rf_CTL_AD_5folds[[3]]$variable.importance)
rf_CTL_AD_4 <- data.frame(obj_rf_CTL_AD_5folds[[4]]$variable.importance)
rf_CTL_AD_5 <- data.frame(obj_rf_CTL_AD_5folds[[5]]$variable.importance)

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
total_rf_important[[1]] <- total_rf_CTL_AD

# CTL - MCI

rf_CTL_MCI_1 <- data.frame(obj_rf_CTL_MCI_5folds[[1]]$variable.importance)
rf_CTL_MCI_2 <- data.frame(obj_rf_CTL_MCI_5folds[[2]]$variable.importance)
rf_CTL_MCI_3 <- data.frame(obj_rf_CTL_MCI_5folds[[3]]$variable.importance)
rf_CTL_MCI_4 <- data.frame(obj_rf_CTL_MCI_5folds[[4]]$variable.importance)
rf_CTL_MCI_5 <- data.frame(obj_rf_CTL_MCI_5folds[[5]]$variable.importance)

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
total_rf_important[[2]] <- total_rf_CTL_MCI

# MCI - AD

rf_MCI_AD_1 <- data.frame(obj_rf_MCI_AD_5folds[[1]]$variable.importance)
rf_MCI_AD_2 <- data.frame(obj_rf_MCI_AD_5folds[[2]]$variable.importance)
rf_MCI_AD_3 <- data.frame(obj_rf_MCI_AD_5folds[[3]]$variable.importance)
rf_MCI_AD_4 <- data.frame(obj_rf_MCI_AD_5folds[[4]]$variable.importance)
rf_MCI_AD_5 <- data.frame(obj_rf_MCI_AD_5folds[[5]]$variable.importance)

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
total_rf_important[[3]] <- total_rf_MCI_AD

CTL_AD_rf <- total_rf_important[[1]]
CTL_MCI_rf <- total_rf_important[[2]]
MCI_AD_rf <- total_rf_important[[3]]

CTL_AD_rf$GENES<- row.names(CTL_AD_rf)
CTL_MCI_rf$GENES<- row.names(CTL_MCI_rf)
MCI_AD_rf$GENES<- row.names(MCI_AD_rf)

CTL_AD_rf_filt <- filter(CTL_AD_rf, Total >= 1) # 40
CTL_MCI_rf_filt <- filter(CTL_MCI_rf, Total >= 1) # 45
MCI_AD_rf_filt <- filter(MCI_AD_rf, Total >= 0.7) # 42

rm(list=setdiff(ls(), c("total_rf_important","CTL_AD_rf_filt", "CTL_MCI_rf_filt", "MCI_AD_rf_filt")))
save(total_rf_important, CTL_AD_rf_filt, CTL_MCI_rf_filt, MCI_AD_rf_filt, file = "rf_filtered.RData")
load("rf_filtered.RData")
#save("total_rf_important", file = "total_rf_significant.RData")

##TAKE THE MOST IMPORTANT VARIANTS

barplot(total_rf_important[[3]]$Total, xlim = c(1, 200), main = "Barplot MCI_AD")
barplot(total_rf_important[[2]]$Total, xlim = c(1, 200), main = "Barplot CTL_MCI")
barplot(total_rf_important[[1]]$Total, xlim = c(1, 200), main = "Barplot CTL_AD")

## ggplots

## keep first 200

y_CTL_AD <- head(total_rf_important[[1]] , n=100)
y_CTL_MCI <- head(total_rf_important[[2]] , n=200)
y_MCI_AD <- head(total_rf_important[[3]] , n=200)

y_CTL_AD$rownames <- rownames(y_CTL_AD)
y_CTL_MCI$rownames <- rownames(y_CTL_MCI)
y_MCI_AD$rownames <- rownames(y_MCI_AD)

p <- ggplot(data = y_CTL_AD, aes(x=reorder(rownames, -Total),y=Total)) +
  geom_bar(stat="identity")
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## SUPPORT VECTOR MACHINE

#load("svm_obj.RData")
load("svm_obj.RData")

total_svm_important <- vector("list")

#CTL_AD
svm_CTL_AD_1 <- data.frame(obj_svm_CTL_AD_5folds[[1]]$importance)
svm_CTL_AD_2 <- data.frame(obj_svm_CTL_AD_5folds[[2]]$importance)
svm_CTL_AD_3 <- data.frame(obj_svm_CTL_AD_5folds[[3]]$importance)
svm_CTL_AD_4 <- data.frame(obj_svm_CTL_AD_5folds[[4]]$importance)
svm_CTL_AD_5 <- data.frame(obj_svm_CTL_AD_5folds[[5]]$importance)

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
total_svm_important[[1]] <- total_svm_CTL_AD

# CTL - MCI

svm_CTL_MCI_1 <- data.frame(obj_svm_CTL_MCI_5folds[[1]]$importance)
svm_CTL_MCI_2 <- data.frame(obj_svm_CTL_MCI_5folds[[2]]$importance)
svm_CTL_MCI_3 <- data.frame(obj_svm_CTL_MCI_5folds[[3]]$importance)
svm_CTL_MCI_4 <- data.frame(obj_svm_CTL_MCI_5folds[[4]]$importance)
svm_CTL_MCI_5 <- data.frame(obj_svm_CTL_MCI_5folds[[5]]$importance)

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
total_svm_important[[2]] <- total_svm_CTL_MCI

# MCI - AD

svm_MCI_AD_1 <- data.frame(obj_svm_MCI_AD_5folds[[1]]$importance)
svm_MCI_AD_2 <- data.frame(obj_svm_MCI_AD_5folds[[2]]$importance)
svm_MCI_AD_3 <- data.frame(obj_svm_MCI_AD_5folds[[3]]$importance)
svm_MCI_AD_4 <- data.frame(obj_svm_MCI_AD_5folds[[4]]$importance)
svm_MCI_AD_5 <- data.frame(obj_svm_MCI_AD_5folds[[5]]$importance)

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
total_svm_important[[3]] <- total_svm_MCI_AD

## PLOTTING RESULTS

barplot(total_svm_important[[3]]$Total, xlim = c(1, 100), main = "Barplot MCI_AD")
barplot(total_svm_important[[2]]$Total, xlim = c(1, 600), main = "Barplot CTL_MCI")
barplot(total_svm_important[[1]]$Total, xlim = c(1, 600), main = "Barplot CTL_AD")

## FILTERING ACCORDING TO BARPLOT

CTL_AD_svm <- total_svm_important[[1]]
CTL_MCI_svm <- total_svm_important[[2]]
MCI_AD_svm <- total_svm_important[[3]]

CTL_AD_svm$GENES<- row.names(CTL_AD_svm)
CTL_MCI_svm$GENES<- row.names(CTL_MCI_svm)
MCI_AD_svm$GENES<- row.names(MCI_AD_svm)


CTL_AD_svm_filt <- filter(CTL_AD_svm, Total >= 3.2) # 186
CTL_MCI_svm_filt <- filter(CTL_MCI_svm, Total >= 3.1) # 131
MCI_AD_svm_filt <- filter(MCI_AD_svm, Total >= 2.85) # 114

#COMMON GENES FROM RF SVM

CTL_AD_TOTAL_GENES <- CTL_AD_rf_filt[CTL_AD_rf_filt$GENES %in% CTL_AD_svm_filt$GENES,]
MCI_AD_TOTAL_GENES <- MCI_AD_rf_filt[MCI_AD_rf_filt$GENES %in% MCI_AD_svm_filt$GENES,]
CTL_MCI_TOTAL_GENES <- CTL_MCI_rf_filt[CTL_MCI_rf_filt$GENES %in% CTL_MCI_svm_filt$GENES,]

CTL_AD_TOTAL_GENES2 <- CTL_AD_svm_filt[CTL_AD_svm_filt$GENES %in% CTL_AD_rf_filt$GENES,]
MCI_AD_TOTAL_GENES2 <- MCI_AD_svm_filt[MCI_AD_svm_filt$GENES %in% MCI_AD_rf_filt$GENES,]
CTL_MCI_TOTAL_GENES2 <- CTL_MCI_svm_filt[CTL_MCI_svm_filt$GENES %in% CTL_MCI_rf_filt$GENES,]

colnames(CTL_AD_TOTAL_GENES2) <- c("val1","val2","val3","val4","val5", "total2", "GENES")
colnames(CTL_MCI_TOTAL_GENES2) <- c("val1","val2","val3","val4","val5", "total2", "GENES")
colnames(MCI_AD_TOTAL_GENES2) <- c("val1","val2","val3","val4","val5", "total2", "GENES")

CTL_AD_FINAL <- merge(CTL_AD_TOTAL_GENES, CTL_AD_TOTAL_GENES2, by="GENES")
CTL_MCI_FINAL <- merge(CTL_MCI_TOTAL_GENES, CTL_MCI_TOTAL_GENES2, by="GENES")
MCI_AD_FINAL <- merge(MCI_AD_TOTAL_GENES, MCI_AD_TOTAL_GENES2, by="GENES")

CTL_AD_FINAL$FINAL <- CTL_AD_FINAL$Total + CTL_AD_FINAL$total2
CTL_MCI_FINAL$FINAL <- CTL_MCI_FINAL$Total + CTL_MCI_FINAL$total2
MCI_AD_FINAL$FINAL <- MCI_AD_FINAL$Total + MCI_AD_FINAL$total2

CTL_AD_FINAL <- CTL_AD_FINAL[order(-CTL_AD_FINAL$FINAL),]
CTL_MCI_FINAL <- CTL_MCI_FINAL[order(-CTL_MCI_FINAL$FINAL),] 
MCI_AD_FINAL <- MCI_AD_FINAL[order(-MCI_AD_FINAL$FINAL),] 

CTL_AD_GENES <- CTL_AD_FINAL[["GENES"]]
CTL_MCI_GENES <- CTL_MCI_FINAL[["GENES"]]
MCI_AD_GENES <- MCI_AD_FINAL[["GENES"]]


save(CTL_AD_FINAL, CTL_MCI_FINAL, MCI_AD_FINAL, file = "total_genes.RData")
load("total_genes.RData")

CTL_AD_FINAL <- CTL_AD_FINAL[,c("GENES","FINAL")]
CTL_MCI_FINAL <- CTL_MCI_FINAL[,c("GENES","FINAL")]
MCI_AD_FINAL <- MCI_AD_FINAL[,c("GENES","FINAL")]

## PLOTING RESULTS GENES - AUC - 

p <- ggplot(data = CTL_AD_FINAL, aes(x=reorder(GENES, FINAL),y=FINAL)) +
  geom_bar(stat="identity") +
  xlab("Genes") +
  ylab("Frequency") +
  ggtitle("Control - Alzheimer's Disease") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
p + theme(axis.text.y = element_text(hjust = 0.4, size = 6),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
) + coord_flip()

p <- ggplot(data = MCI_AD_FINAL, aes(x=reorder(GENES, FINAL),y=FINAL)) +
  geom_bar(stat="identity") +
  xlab("Genes") +
  ylab("Frequency") +
  ggtitle("Mild Cognitive Impairment - Alzheimer's Disease") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
p + theme(axis.text.y = element_text(hjust = 0.4, size = 6),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
) + coord_flip()

## GENES BARPLOT

p <- ggplot(data = CTL_MCI_FINAL, aes(x=reorder(GENES, FINAL),y=FINAL)) +
  geom_bar(stat="identity") +
  xlab("Genes") +
  ylab("Frequency") +
  ggtitle("Control - Mild Cognitive Impairment") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
p + theme(axis.text.y = element_text(hjust = 0.4, size = 6),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
) + coord_flip()

## AUC BOXPLOTS

p <-ggplot(df, aes(x=condition, y=notes, fill=machine)) + 
  geom_boxplot() +
  xlab("Conditions") +
  ylab("AUC") +
  ggtitle("5 Partitions AUC Range") +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
p + theme(axis.text.y = element_text(hjust = 0.4, size = 6),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15))

## TAKE SNPs FROM ENSEMBLE ANNOTATION

MCI_AD_rs <- read.table("MCI_AD_ANNO.txt", header = FALSE)
CTL_AD_rs <- read.table("CTL_AD_ANNO.txt", header = FALSE)
CTL_MCI_rs <- read.table("CTL_MCI_ANNO.txt", header = FALSE)

MCI_AD_rs <- MCI_AD_rs[c(1,6)]
CTL_AD_rs <- CTL_AD_rs[c(1,6)]
CTL_MCI_rs <- CTL_MCI_rs[c(1,6)]

colnames(CTL_AD_rs) <- c("rs_Name","GENES")
colnames(MCI_AD_rs) <- c("rs_Name","GENES")
colnames(CTL_MCI_rs) <- c("rs_Name","GENES")

CTL_AD_rs<-CTL_AD_rs[!duplicated(CTL_AD_rs[,c("rs_Name")]),]
MCI_AD_rs<-MCI_AD_rs[!duplicated(MCI_AD_rs[,c("rs_Name")]),]
CTL_MCI_rs<-CTL_MCI_rs[!duplicated(CTL_MCI_rs[,c("rs_Name")]),]

## CHECK IF THERE ARE MATCHES IN TRANSCRIPT AND GENOMICS
#RANDOM FOREST

CTL_AD_rf_match <- CTL_AD_rs[CTL_AD_rs$GENES %in% rownames(total_rf_CTL_AD),] # 0
CTL_MCI_rf_match <- CTL_MCI_rf_filt[CTL_MCI_rf_filt$GENES %in% CTL_MCI_rs$GENES,] # 0
MCI_AD_rf_match <- MCI_AD_rf_filt[MCI_AD_rf_filt$GENES %in% MCI_AD_rs$GENES,] # 0


#SUPPORT VECTOR MACHINE

CTL_AD_svm_match <- CTL_AD_svm_filt[CTL_AD_svm_filt$GENES %in% CTL_AD_rs$GENES,] # 0
CTL_MCI_svm_match <- CTL_MCI_svm_filt[CTL_MCI_svm_filt$GENES %in% CTL_MCI_rs$GENES,] # 0
MCI_AD_svm_match <- MCI_AD_svm_filt[MCI_AD_svm_filt$GENES %in% MCI_AD_rs$GENES,] # 0


x<-CTL_AD_rf_match[!duplicated(CTL_AD_rf_match[,c("GENES")]),]

y <- filter(x, Gene.name == vector)

## USE ENSEMBL FROM BIOMART

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

it2 <-getBM(attributes = c("external_gene_name", "phenotype_description"), 
            filters="external_gene_name",
            values = genes_names,
            mart = ensembl)

### END ###