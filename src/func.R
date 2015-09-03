# return counts matrix based on TCGA tissue and sample type 
readTCGA <- function(countsFile, tissue='BRCA', type='tumor') {
  zz <- gzfile(countsFile,
               'rt')  
  dat <- read.csv(zz, sep='\t')
  close(zz)
  colnames(dat) <- gsub('\\.', '-', colnames(dat))
  gnames <- dat$X
  
  if(type=='tumor') {
    all <- read.csv2('/home/tobias/AWS/s3/averaprojects/tcga/GSE62944_06_01_15_TCGA_24_CancerType_Samples.txt', 
                          header=F, 
                          stringsAsFactor=F, 
                          sep='\t')
  }
  
  if(type=='normal') {
    all <- read.csv2('/home/tobias/AWS/s3/averaprojects/tcga/GSE62944_06_01_15_TCGA_24_Normal_CancerType_Samples.txt', 
                           header=F, 
                           stringsAsFactor=F, 
                           sep='\t')
  }
  id <- all$V1[all$V2==tissue]
  dat <- dat[,id]
  rownames(dat) <- gnames
  return(dat)
}

# return tcga clinical data based on vector of sample ids
clinicalTCGA <- function(id) {
  allClinical <- t(read.csv2('/home/tobias/AWS/s3/averaprojects/tcga/GSE62944_06_01_15_TCGA_24_548_Clinical_Variables_9264_Samples.txt', 
                             sep='\t'))
  colnames(allClinical) <- allClinical[1,]
  allClinical <- allClinical[-c(1:3),]
  rownames(allClinical) <- gsub('\\.', '-', rownames(allClinical))
  clinical <- allClinical[match(id, rownames(allClinical)), ]
}

readGTEX <- function(tissue) {
  dat <- read.csv2('/home/tobias/AWS/s3/averaprojects/GTEX/GTEx_Analysis_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct', 
                    sep='\t', 
                    skip=2)
  colnames(dat) <- gsub('\\.', '-', colnames(dat))
  dat0 <- aggregate(. ~ Description, sum, data=dat)
  gtexSamples <- read.csv2('/home/tobias/AWS/s3/averaprojects/GTEX/GTEx_Data_V4_Annotations_SampleAttributesDS.txt', 
                           sep='\t', 
                           stringsAsFactor=F)
  id <- gtexSamples$SAMPID[gtexSamples$SMTS==tissue]
  dat1 <- dat0[,colnames(dat0) %in% id]
  rownames(dat1) <- dat0$Description
  return(dat1)
}