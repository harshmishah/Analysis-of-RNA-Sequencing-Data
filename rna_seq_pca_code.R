#!/usr/bin/R

library(DESeq2)
library(ggplot2)
library(openxlsx)
library(tidyverse)
library(ggfortify)
library(sva)
library(pcaExplorer)
library(edgeR)
library(rafalib)

#Read Count Data
countData <- read.csv('Exp6_Exp3_Counts_Test.csv', header = TRUE, sep = ",")

#Read Metadata
metaData <- read.xlsx('Exp6_Exp3_Counts_Test_Metadata.xlsx')

#Format Data
coldata = metaData[-1]
row.names(coldata) <- metaData$Sample_Filenames

coldata <- coldata[,c("Experiment","Time","Donor")]
coldata$Experiment <- factor(coldata$Experiment)
coldata$Time <- factor(coldata$Time)
coldata$Donor <- factor(coldata$Donor)
all(rownames(coldata) == colnames(countData))


countdata <- countData %>%
	    column_to_rownames("Gene") %>% # turn the Gene column into rownames
	        as.matrix()
head(countdata)

#Build DESeq Data Set
dds <- DESeqDataSetFromMatrix(countData = countdata,
			      colData = coldata,
			      design = ~ Experiment + Time)

#Normalize Data
ddsObj <- estimateSizeFactors(dds)
ddsObj@colData$sizeFactor
normalized_count <- counts(ddsObj, normalized=TRUE)

dds <- DESeqDataSetFromMatrix(countData = round(normalized_count),
			      colData = coldata,
			      design = ~ Experiment + Time)
 
#Transform Data
dds <- DESeq(dds)
vsd <- vst(dds,blind=TRUE)

#Plot PCA Results for top 5000 genes
jpeg('pca_time_exp_5000.jpg')
plotPCA(vsd,intgroup=c("Time","Experiment"),ntop=5000)
dev.off()

#Plot PCA Results for top 500 genes
jpeg('pca_time_exp_500.jpg')
plotPCA(vsd,intgroup=c("Time","Experiment"),ntop=500)
dev.off()

#Batch Effect Model
#Provide 2 model matrices to SVA
mm <- model.matrix(~ Donor, colData(dds))
mm0 <- model.matrix(~ 1, colData(dds))
normalized_count <- normalized_count[rowSums(normalized_count) > 0,]
fit <- svaseq(normalized_count, mod=mm, mod0=mm0, n.sv=2)

#Plot estimated surrogate variables
bigpar()
dds$Donor.int <- as.integer(dds$Donor) + 15
dds$Donor.int

jpeg('batch_effect.jpg')
plot(fit$sv[,1:2], col=dds$Experiment, pch=dds$Donor.int, cex=2,
          xlab="SV1", ylab="SV2")
legend("top", levels(dds$Experiment
		     ), pch=16,
              col=1:3, cex=.8, ncol=3, title="batch")
dev.off()

#Removing Batch Effect
mat <- assay(vsd)
mat <- limma::removeBatchEffect(mat, vsd$Experiment)
assay(vsd) <- mat

#Re-plotting PCA Results for top 500 genes
jpeg('pca_time_exp_500_post_batch.jpg')
plotPCA(vsd,intgroup=c("Time","Experiment"),ntop=500)
dev.off()

#Re-plotting PCA Results for top 5000 genes
jpeg('pca_time_exp_5000_post_batch.jpg')
plotPCA(vsd,intgroup=c("Time","Experiment"),ntop=5000)
dev.off()
