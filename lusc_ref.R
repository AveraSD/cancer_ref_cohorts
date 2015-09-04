# LUSC Reference Cohort
# original Author: Tobias Meissner

rm(list=ls())

library(plyr)
library(mygene)
library(reshape2)
library(DESeq2)
library(BiocParallel)
library(parallel)
library(EDASeq)
library(RUVSeq)
library(RColorBrewer)

source('src/func.R')

## Read in TCGA data
datTumor <- readTCGA('/home/tobias/AWS/s3/averaprojects/tcga/GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FeatureCounts.txt.gz', 'LUSC', 'tumor')
datNormal <- readTCGA('/home/tobias/AWS/s3/averaprojects/tcga/GSM1697009_06_01_15_TCGA_24.normal_Rsubread_FeatureCounts.txt.gz', 'LUSC', 'normal')

## get clinical data
clinical <- clinicalTCGA(colnames(datTumor))
clinicalN <- clinicalTCGA(colnames(datNormal))

## read in GTEX normal samples
datGTEX <- readGTEX('Lung')

# combine data
commonHGNC <- intersect(rownames(datGTEX), rownames(datTumor))
df <- cbind(datGTEX[commonHGNC,], datNormal[commonHGNC,], datTumor[commonHGNC,])

# preprocess data
colData <- data.frame(condition=factor(rep(c('NORMAL', 'MNORMAL', 'TUMOR'),
                                           c(dim(datGTEX)[2], dim(datNormal)[2], dim(datTumor)[2])))
)
dds <- DESeqDataSetFromMatrix(df, colData, formula(~ condition))
rownames(colData(dds)) <- colnames(df)

workers <- 8
register(MulticoreParam(workers))
dds <- DESeq(dds, parallel=T)

sfDESeq <- sizeFactors(dds) # save the size factors of the ref cohort
loggeomeansRef <- rowMeans(log(df))

save(loggeomeansRef, file='ref_data/loggeoameansLUSC.Rdata')

reference <- list(reference_dds=dds,
                  reference_raw_count=df,
                  reference_sf=sfDESeq,
                  group=colData$condition)

save(reference, file='ref_data/DESeq_TCGA_GTEX_LUSC.Rdata')


