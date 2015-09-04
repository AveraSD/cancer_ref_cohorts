# OV Reference Cohort
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
# note: there are no matched normal samples for OV in TCGA
datTumor <- readTCGA('/home/tobias/AWS/s3/averaprojects/tcga/GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FeatureCounts.txt.gz', 'OV', 'tumor')

## get clinical data
clinical <- clinicalTCGA(colnames(datTumor))
clinicalN <- clinicalTCGA(colnames(datNormal))

## read in GTEX normal samples
datGTEX <- readGTEX('Breast')

# combine data
hg19syms <- read.csv2('/home/tobias/AWS/database/Homo_sapiens/UCSC/hg19/Annotation/Genes/geneid.txt', 
                      header=F,
                      stringsAsFactor=F)$V1

commonHGNC <- intersect(intersect(rownames(datGTEX), rownames(datTumor)), hg19syms)
df <- cbind(datGTEX[commonHGNC,], datTumor[commonHGNC,])


# preprocess data
colData <- data.frame(condition=factor(rep(c('NORMAL', 'TUMOR'),
                                          c(dim(datGTEX)[2], dim(datTumor)[2])))
)
dds <- DESeqDataSetFromMatrix(df, colData, formula(~ condition))
rownames(colData(dds)) <- colnames(df)

workers <- 32
register(MulticoreParam(workers))
dds <- DESeq(dds, parallel=F)

sfDESeq <- sizeFactors(dds) # save the size factors of the ref cohort
loggeomeansRef <- rowMeans(log(df))
vsdRef <- varianceStabilizingTransformation(dds)

save(loggeomeansRef, file='ref_data/loggeoameansOV.Rdata')

reference <- list(reference_dds=dds,
                  reference_raw_count=df,
                  reference_sf=sfDESeq,
                  reference_vsd=vsdRef,
                  group=colData$condition)

save(reference, file='ref_data/DESeq_TCGA_GTEX_OV.Rdata')

