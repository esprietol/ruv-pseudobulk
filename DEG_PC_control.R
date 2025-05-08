rm(list=ls())
setwd("/mnt")

BiocManager::install(c('EDASeq','iCOBRA','RUVSeq'),update = F)
install.packages(c('rlist','ruv','factoextra'))
system("R CMD INSTALL /domino/datasets/local/RUV/swapper-master/")

library(BiocManager)
library(tidyverse)
library(DESeq2)
library(EDASeq)
library(scuttle)
library(scater)
library(uwot)
library(edgeR)
library(Seurat)
library(rlist)
library(cluster)
library(factoextra)
library(UpSetR)  
library(foreach)
library(doParallel)

numCores <- detectCores()

# Create a parallel backend with the number of cores
cl <- makeCluster(8)
registerDoParallel(cl)
# sc data -----------------------------------------------------------------

source("/mnt/auxf.R")
ds2 <- readRDS("/domino/datasets/local/RUV/sclupusds2.rds")
ds2 <- ds2[,ds2$Sex=='Female']
Idents(object = ds2) <- "cg_cov"
celltypes <- c('B', 'NK', 'T4', 'T8', 'cDC', 'cM', 'ncM', 'pDC')
ds2s <- subset(x = ds2, idents = celltypes)
rm(list = 'ds2')

pdata <- droplevels(ds2s@meta.data)
pdata$Age <- as.numeric(as.character(pdata$Age))

ds2s <- SetIdent(ds2s, value = "sample_cell") # change ident to group later
ds2s <- FindVariableFeatures(ds2s) #### adds metafeatures on ds2s@assays[["originalexp"]]@meta.features
ds2s <- ScaleData(ds2s) #### adds scaledata on ds2s@assays[["originalexp"]]@data

# pseudobulk --------------------------------------------------------------

PBC <- AggregateExpression(ds2s, group_by = "sample_cell", assays = 'originalexp',slot='counts', return.seurat = F)
PBC <- PBC$originalexp # same result than manual procedure
colnames(PBC) <- unique(pdata$sample_cell)
metacovs <- colnames(pdata)[c(4:6,8,10:16, 18:19)] # variables at type cell-sample level
pbc.metaD <- unique(pdata[,metacovs])
bc.metaD <- unique(pbc.metaD[,!colnames(pbc.metaD) %in% c("cg_cov" ,"sample_cell" )])
rownames(pbc.metaD) <- pbc.metaD$sample_cell 
pbc.metaD$Age <- as.numeric(as.character(pbc.metaD$Age))



path <- ('/domino/datasets/local/RUV')

# Global filter and DGEList -----------------------------------------------------------

PBC <- PBC[rowSums(PBC >= 10) >= 5,]
gene.D <- data.frame(gene=rownames(PBC))
rownames(gene.D) <- gene.D$gene

# hyperparameters -----------------------------------------------------------


DEG <- foreach(ct = celltypes) %do% {
  
  pbc.ct <- PBC[,pbc.metaD$cg_cov==ct]
  pbc.ct <-  pbc.ct[rowSums(pbc.ct >= 10) >= 5,]
  gene.D.ct <- gene.D[ rownames(pbc.ct),]
  
  y.ct <- DGEList(counts=pbc.ct )
  pmin.ct <- find_p(pbc.ct )
  y.ct <- calcNormFactors(y.ct, method="upperquartile", p=pmin.ct)
  #nf.ct <- y.ct$samples$norm.factors
  #logy.ct <- voom(y.ct, plot=F)
  
  design <- model.matrix(~ 0 + pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct] ) 
  colnames(design) <- c('PC1','PC2', 'PC3', 'PC4')
  v <- voom(y.ct, design, plot=F) 
  vfit <- lmFit(v, design)
  
  
  contrast.matrix <- makeContrasts(
    PC2_vs_PC1 = PC2 - PC1,
    PC3_vs_PC1 = PC3 - PC1,
    PC4_vs_PC1 = PC4 - PC1,
    PC3_vs_PC2 = PC3 - PC2,
    PC4_vs_PC2 = PC4 - PC2,
    PC4_vs_PC3 = PC4 - PC3,
    levels = colnames(design)
  )
  
  contrast.matrix <- makeContrasts(
    PC2 - PC1, PC3 - PC1, PC4 - PC1, PC3 - PC2, PC4 - PC2, PC4 - PC3,
    levels = colnames(design)
  )
  
  # Fit contrasts
  vfit <- contrasts.fit(vfit, contrast.matrix)
  efit <- eBayes(vfit)
  topTable(efit, number=dim(y.ct)[1], adjust="BH", p.value = 0.01)
  
}

names(DEG)<- celltypes 

saveRDS(DEG,paste0(path,'/DiffGPCs_popE.rds'))





