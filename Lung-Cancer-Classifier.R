###########################################
#  - Lung Cancer DEG analysis and classifier#
#  - TCGA-LUAD project - RNAseq data      #                        #    
#  - 2023-3-10                          #
#  - Copyright: @Youssef Fahim           #
###########################################
# R version 4.2.0 (2022-04-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: windows 10
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

####load libraries ########  
library(readr)
library(DESeq2)
library(ComplexHeatmap)
library(ggplot2)
library(scatterplot3d)
library(dplyr)
library(ggfortify)
library(rgl)
library(GOplot)
library(clusterProfiler)
library(pathview)
library(org.Hs.eg.db)
library(plotly)
library(caret)
library(caretEnsemble)
library(doParallel)
library(randomForest)
library(randomForestSRC)
library(pROC)
library(e1071)
library(class)
load(".RDATA")
data.path="Dataset/TCGA-LUAD-RNA_SEQ"
pheno.path="Dataset/Sample_Sheet.tsv"

##### load the mRNA-Seq data #####
files <- list.files(path=data.path,recursive=T, pattern = "tsv")

# read the first file for the first time
file=files[1]
filepath=file.path(data.path,files[1])
colnames= unlist( strsplit ( readLines(filepath ,n=2) [2] ,"\t"))

temp <- read.table(filepath, header=F,skip = 6)
names(temp)=colnames
genetype_mapper=temp[,c(2,3)]

#create a storing object exp to save the whole TPMs of each file TPMs read in an iteration
exp=temp[temp$gene_type == "protein_coding" ,c(1,2)]

for(i in 1: length(files))
{
  ## refer to the next file (note that we start from index 2, bec we already read the first file)
  file=files[i]
  file.id=strsplit(file,"/")[[1]][1]
  filepath=file.path(data.path,files[i])
  
  # read the next file  
  temp <- read.table(filepath, header=F,skip = 6)
  temp=temp [ temp[,3] == "protein_coding" , ]
  
  ## change the colname of the tpm with sample name 
  exp=cbind(exp,temp[,4])
  colnames(exp)[dim(exp)[2]]=file.id
}
View(exp)

#########################  Assignment #########################
## load the data using the two other approaches in the following link
# a professional example to follow is here using two ways ( reading from files , TCGAbiolink package) >>>  https://www.biostars.org/p/9500223/

### General Reminder
# don't forget to master the tidyvers packages https://www.tidyverse.org/packages/


# check duplciation of of gene symbols?
x=duplicated(exp$gene_name)
sum(x)

### yes .. why ? transcripts?  solutions : aggregation
exp.data=exp[ , 3:dim(exp)[2]]
exp.data=apply(exp.data,2, as.numeric)

#### remove duplication by aggregation
exp.data.agg= aggregate(exp.data, list(exp$gene_name),FUN=max)
genes=exp.data.agg$Group.1
exp.data.agg=exp.data.agg[-1]
exp.data.agg=apply(exp.data.agg,2, as.numeric)
exp.data.agg=round(exp.data.agg)
rownames(exp.data.agg)=genes
dim(data)
file.ids=colnames(exp.data.agg)

###### load the mrna sample sheets  # sample sheets
pheno <- read_delim(pheno.path,"\t", escape_double = FALSE, trim_ws = TRUE)
table(pheno$`Sample Type`)
names(pheno)=sub (" ", "_",names(pheno)) # rename the names of the pheno to add "_" instead of " "

### Excluding recurrent tumor cases
exp.data.agg <- subset(exp.data.agg, select = -c(`347db78b-ddfc-46ed-9fe4-0ca52ac432ef`,`eb29ba2f-e237-403c-ab72-a5dfae705b05`))
pheno <- pheno %>% filter(Sample_Type!='Recurrent Tumor')

#we will rename the columns of our exp data with the sample ids columns of the pheno file
#however we need to match the file ids
file.ids.pheno=pheno$File_ID
index.files=match(file.ids,file.ids.pheno)
colnames(exp.data.agg) = pheno$Sample_ID[index.files] %>% na.omit()


#### Exploratory analysis + filteration process  ####################

#### Do differential analysis using Deseq2 package as it works on readcounts  ##########
## or use Limma package instead >> you can follow this tutorial here https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html


table(pheno$Sample_Type)

cond1="Solid Tissue Normal" 
cond2="Primary Tumor"

dds = DESeqDataSetFromMatrix( countData = exp.data.agg, colData = pheno , design = ~ Sample_Type)
dds.run = DESeq(dds)

### direct results or specifying the contrast (to make a res object based on two specific conditions/treatment)
#res=results(dds.run)
res=results(dds.run, contrast = c("Sample_Type",cond1 ,cond2))

# remove nulls
res=res[complete.cases(res), ]

res.df=as.data.frame(res)
write.table(res.df, file = "res.txt")

plotMA(res, ylim=c(-1,1)) 

res.degs=res.df[res.df$padj< 0.05 & abs(res.df$log2FoldChange)>log2(2),]
res.degs=res.degs[order(res.degs$padj), ]

# add a column for expression level classification
res.degs$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
res.degs$diffexpressed[res.degs$log2FoldChange > 1] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res.degs$diffexpressed[res.degs$log2FoldChange < -1] <- "DOWN"

# add a column for expression level classification
res.df$diffexpressed <- "Normal"
# if log2Foldchange > 1 and pvalue < 0.05, set as "UP" 
res.df$diffexpressed[res.df$log2FoldChange > 1 & res.df$padj < 0.05] <- "UP"
# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
res.df$diffexpressed[res.df$log2FoldChange < -1 & res.df$padj < 0.05] <- "DOWN"

degs.genes=rownames(res.degs)
genes = row.names(res.df)

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
res.degs$delabel <- NA
res.degs$delabel[res.degs$diffexpressed != "NO"] <- degs.genes[res.degs$diffexpressed != "NO"]

res.df$delabel <- NA
res.df$delabel[res.df$diffexpressed != "Normal"] <- genes[res.df$diffexpressed != "Normal"]

# or export them to a gene list to start the gene set enrichment analysis.
write.table(degs.genes, file="res.degs.txt", quote = F, col.names = F, row.names = F)
write.csv(res.degs, "final.csv")

#### get the normalized and logged transformed values of all exp data
#using the the variance stabilizing transformation. vsn package

ntd = normTransform(dds)
exp.norm = assay(ntd)
exp.norm = exp.norm[, order(colnames(exp.norm))]

## 1- creating a heatmap for the top 100 DEG genes
#  get the expression profiles of the top 100 degs only and create heatmap
top100_degs.exp <- degs.exp[1:100, ]
top100_degs.exp <- top100_degs.exp[, order(colnames(top100_degs.exp))]

pheno_ordered_matched <- pheno[match(colnames(top100_degs.exp), pheno$Sample_ID),]
pheno_ordered_sorted <- pheno_ordered_matched[order(pheno_ordered_matched$Sample_Type),]
top100_degs.exp_sorted <- top100_degs.exp[, pheno_ordered_sorted$Sample_ID]

sample_type <- pheno_ordered_sorted$Sample_Type  # Replace with your actual column name from pheno
colors <- c("Solid Tissue Normal" = "blue", "Primary Tumor" = "red")  # Choose your preferred colors
annotation <- HeatmapAnnotation(df = data.frame(SampleType = sample_type),
                                col = list(SampleType = colors),
                                show_annotation_name = FALSE)
scaled_data <- t(scale(t(as.matrix(top100_degs.exp_sorted))))

Heatmap(as.matrix(scaled_data),top_annotation = annotation ,cluster_columns = F, name = "expression", row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 2))

### 2- creating 2D PCA for all degs and look how it segregates the healthy and tumor samples ### 2- creating 2D PCA foTRUEr all degs and look how it segregates the healthy and tumor samples 

constant_columns <- apply(t(exp.norm), 2, var) == 0
# Filter out the constant columns
exp.norm.filtered <- t(exp.norm)[, !constant_columns]
######### get the normalized expression levels of the degs ###################
degs.exp = exp.norm.filtered[,degs.genes]

pca_result <- prcomp(exp.norm.filtered, center = TRUE, scale. = TRUE)
pheno_ordered <- pheno[order(pheno$Sample_ID), ]##order the samples in pheno to concatenate the sample type to the scores_3D

plot_ly(colors=c('red', 'green'), x = pca_result$x[,1], y = pca_result$x[,2], type = "scatter", mode = "markers", hoverinfo = "text", color=pheno_ordered$Sample_Type)
plot_ly(colors=c('red', 'green'), x = pca_result$x[,1], y = pca_result$x[,2], z = pca_result$x[,3], type = "scatter3d", mode = "markers", hoverinfo = "text", color=pheno_ordered$Sample_Type)

####################################################################################
save(dds,dds.run,exp.data.agg,genetype_mapper,  pheno, res.df,res.degs, degs.exp, exp.norm, pca_result, file="data.RDATA") 
####################################################################################

### 4- creating  Volcano Plot using the advanced methods

#Assuming 'res' is a dataframe of your DESeq2 results including log2FoldChange and padj columns

ggplot(res.df, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() + 
  theme(axis.text = element_text(size = 12)) + 
  scale_color_manual(values=c("blue", "black", "red")) + 
  geom_vline(xintercept=c(-1, 1), col="red") + 
  geom_hline(yintercept=-log10(0.002), col="red") + geom_text(aes(label = delabel), size = 3, vjust = -1) 

### Machine learning using Caret (Degs only)
dataForTraining <- degs.exp
dataForTraining <- t(scale(t(as.matrix(degs.exp))))##testing only
dataForTraining <- data.frame(dataForTraining,stringsAsFactors=FALSE)
dataForTraining <- cbind(dataForTraining,pheno_ordered$Sample_Type)
dataForTraining <- dataForTraining %>% rename(Sample_Type = `pheno_ordered$Sample_Type`)

#preparing the data set

deg.train <- dataForTraining %>% mutate(class = dataForTraining$sample_type)
deg.train$Sample_Type <- gsub(" ", "_", deg.train$Sample_Type)
constant_vars <- sapply(deg.train, function(x) length(unique(x)) == 1)
print(names(constant_vars)[constant_vars])
deg.train <- deg.train[, !constant_vars]
write.csv(deg.train, file = "dataForTraining.csv")

# Assuming deg.train is your training dataset
char_cols <- sapply(deg.train, is.character)
deg.train[char_cols] <- lapply(deg.train[char_cols], as.factor)
model <- rfsrc(Sample_Type ~ ., data = deg.train)

#training and testing data sets
set.seed(123)
is.factor(deg.train$Sample_Type) ##If not Factor convert it to Factor
deg.train$Sample_Type <- as.factor(deg.train$Sample_Type)
trainIndex <- createDataPartition(y = deg.train$Sample_Type, p = 0.7, list = FALSE)
training_data <- deg.train[trainIndex, ]
testing_data <- deg.train[-trainIndex, ]

#training models
algorithms <- c('rf', 'svmLinear', 'glm')

numCores <- detectCores()
cl <- makeCluster(numCores - 2)  # leave one core free for other tasks
registerDoParallel(cl)

train_control <- trainControl(method = "repeatedcv", repeats = 3, classProbs = TRUE, summaryFunction = twoClassSummary)

##training 
#model <- caretList(Sample_Type ~ ., data = training_data, methodList = algorithms, trControl = train_control)

model.rf <- train(Sample_Type ~ ., data = training_data, method="rf", trControl=train_control, tuneLength = 15, metric="ROC")

model.svm <- train(Sample_Type ~ ., data = training_data, method="svmLinear", trControl=train_control, tuneLength = 15, metric="ROC")

model.glm <- train(Sample_Type ~ ., data = training_data, method="glm", trControl=train_control, tuneLength = 15, metric="ROC")

##Accuracy
pred.rf <- predict(model.rf, testing_data)
accuracy_rf <- sum(pred.rf == testing_data$Sample_Type) / nrow(testing_data)
print(paste("rf Accuracy:", accuracy_rf))

pred.svm <- predict(model_svm, testing_data)
accuracy_svm <- sum(pred.svm == testing_data$Sample_Type) / nrow(testing_data)
print(paste("svm Accuracy:", accuracy_svm))

pred.glm <- predict(model.glm, testing_data)
accuracy_glm <- sum(pred.glm == testing_data$Sample_Type) / nrow(testing_data)
print(paste("glm Accuracy:", accuracy_glm))

##importance of genes after training with random forest
imp.rf <- varImp(model.rf, scale= FALSE)
imp.svm <- varImp(model.svm, scale= FALSE)

importance_rf <- imp.rf$importance
importance_svm <- imp.svm$importance

merged_importance <- merge(importance_rf, importance_svm, by = "row.names", all = TRUE)
merged_importance$SumImportance <- merged_importance[order(-merged_importance$Overall),]

##String network and scaling
string_degrees <- `string_node_degrees.(2)` 

##Precision
actual <- testing_data$Sample_Type
true_positives_rf <- sum(pred.rf == "Primary_Tumor" & actual == "Primary_Tumor")
predicted_positives_rf <- sum(pred.rf == "Primary_Tumor")
precision_rf <- true_positives_rf / predicted_positives_rf

true_positives_svm <- sum(pred.svm == "Primary_Tumor" & actual == "Primary_Tumor")
predicted_positives_svm <- sum(pred.svm == "Primary_Tumor")
precision_svm <- true_positives_svm / predicted_positives_svm

true_positives_glm <- sum(pred.glm == "Primary_Tumor" & actual == "Primary_Tumor")
predicted_positives_glm <- sum(pred.glm == "Primary_Tumor")
precision_glm <- true_positives_glm / predicted_positives_glm

##Recall 
false_negatives_rf <- sum(pred.rf == "Solid_Tissue_Normal" & actual == "Primary_Tumor")
recall_rf <- true_positives_rf / (true_positives_rf + false_negatives_rf)

false_negatives_svm <- sum(pred.svm == "Solid_Tissue_Normal" & actual == "Primary_Tumor")
recall_svm <- true_positives_svm / (true_positives_svm + false_negatives_svm)

false_negatives_glm <- sum(pred.glm == "Solid_Tissue_Normal" & actual == "Primary_Tumor")
recall_glm <- true_positives_glm / (true_positives_glm + false_negatives_glm)

##F1 Score
f1_score_rf <- 2 * ((precision_rf * recall_rf) / (precision_rf + recall_rf))

f1_score_svm <- 2 * ((precision_svm * recall_svm) / (precision_svm + recall_svm))

f1_score_glm <- 2 * ((precision_glm * recall_glm) / (precision_glm + recall_glm))

## Plotting results
# preparing data for plotting
metrics_df <- data.frame(
  Model = rep(c("Random Forest", "SVM", "Logistic Regression"), each = 4),
  Metric = rep(c("Accuracy", "Precision", "Recall", "F1 Score"), times = 3),
  Value = c(accuracy_rf, precision_rf, recall_rf, f1_score_rf,   # Values for Model 1
            accuracy_svm, precision_svm, recall_svm, f1_score_svm,   # Values for Model 2
            accuracy_glm, precision_glm, recall_glm, f1_score_glm)   # Values for Model 3
)

metrics_df$Metric <- factor(metrics_df$Metric, levels = c("Accuracy", "Precision", "Recall", "F1 Score"))


# Creating the bar plot
p <- ggplot(metrics_df, aes(x = Model, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  labs(title = "Comparison of Model Metrics",
       x = "Models",
       y = "Metric Value") +
  theme_minimal() +
  theme(text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))

# Display the plot
print(p)

##plotting barplots
metrics_long <- reshape2::melt(metrics_df, id.vars = "Model")

#Accuracy
ggplot(metrics_df, aes(x = Model, y = Accuracy, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(x = "Model", y = "Accuracy", title = "Model Accuracy Comparison") +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 12))

#precision
ggplot(metrics_df, aes(x = Model, y = Precision, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(x = "Model", y = "Precision", title = "Model Precision Comparison") +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 12))

#Recall
ggplot(metrics_df, aes(x = Model, y = Recall, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(x = "Model", y = "Recall", title = "Model Recall Comparison") +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 12))

#F1 Score
ggplot(metrics_df, aes(x = Model, y = F1_Score, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(x = "Model", y = "F1 Score", title = "Model F1 Score Comparison") +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 12))

##Plotting Confusion Matrix
testing_data$Sample_Type <- factor(testing_data$Sample_Type)
preds.rf <- factor(preds.rf)
conf_matrix_rf <- confusionMatrix(data = pred.rf, reference = testing_data$Sample_Type)
fourfoldplot(conf_matrix_rf$table, main = "ROC Curve Random Forest")

conf_matrix_svm <- confusionMatrix(data = pred.svm, reference = testing_data$Sample_Type)
fourfoldplot(conf_matrix_svm$table, main = "ROC Curve Support Vector Machines")

conf_matrix_glm <- confusionMatrix(data = pred.glm, reference = testing_data$Sample_Type)
fourfoldplot(conf_matrix_glm$table, main = "ROC Curve Logistic Regression")

##ROC Curve
testing_data$Sample_Type_numeric <- as.numeric(testing_data$Sample_Type == "Primary_Tumor")
pred.rf.numeric <- as.numeric(pred.rf == "Primary_Tumor")
pred.svm.numeric <- as.numeric(pred.svm == "Primary_Tumor")
pred.glm.numeric <- as.numeric(pred.glm == "Primary_Tumor")

roc_curve_rf <- roc(testing_data$Sample_Type_numeric, pred.rf.numeric)
plot(roc_curve_rf, main = "ROC Curve Random Forest", col = "#1c61b6")
abline(a = 0, b = 1, lty = 2, col = "red")

roc_curve_svm <- roc(testing_data$Sample_Type_numeric, pred.svm.numeric)
plot(roc_curve_svm, main = "ROC Curve Support Vector Machines", col = "#1c61b6")
abline(a = 0, b = 1, lty = 2, col = "red")

roc_curve_glm <- roc(testing_data$Sample_Type_numeric, pred.glm.numeric)
plot(roc_curve_glm, main = "ROC Curve Logistic Regression", col = "#1c61b6")
abline(a = 0, b = 1, lty = 2, col = "red")





# Plot the first ROC curve
plot(roc_curve_rf, col = "blue", main = "ROC Curves Comparison", xlim = c(1, 0), ylim = c(0, 1))

# Add the second and third ROC curves
lines(roc_curve_svm, col = "green", )
lines(roc_curve_glm, col = "red")

# Adding a legend
legend("bottomright", legend = c("Random Forest", "Support Vector Machines", "Logistic Regression"), col = c("blue", "green", "red"), lwd = 2)

