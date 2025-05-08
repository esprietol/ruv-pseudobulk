rm(list=ls())
setwd("/mnt")

BiocManager::install(c('EDASeq','iCOBRA','RUVSeq'),update = F)
install.packages(c('rlist','ruv','factoextra'))
system("R CMD INSTALL /domino/datasets/local/RUV/swapper-master/")

library(tidyverse)
library(BiocManager)
library(DESeq2)
library(EDASeq)
library(scuttle)
library(scater)
library(uwot)
library(edgeR)
library(ruv)
library(Seurat)
library(swapper)
library(rlist)
library(iCOBRA)
library(cluster)
library(factoextra)

# Previously all cells with more than 100 genes, all cells with more than 400 counts, cellls with <15 percentage of mitochondrial genes


# sc data -----------------------------------------------------------------

source("/mnt/auxf.R")
ds2 <- readRDS("/domino/datasets/local/RUV/sclupusds2.rds")
ds2 <- ds2[,ds2$Sex=='Female']
Idents(object = ds2) <- "cg_cov"

celltypes <- c('B', 'NK', 'T4', 'T8', 'cDC', 'cM', 'ncM', 'pDC')
ds2s <- subset(x = ds2, idents = celltypes)
rm(list = 'ds2')
ds2s@meta.data <- droplevels(ds2s@meta.data)
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
rownames(pbc.metaD) <- pbc.metaD$sample_cell 
pbc.metaD$Age <- as.numeric(as.character(pbc.metaD$Age))


# NCG ---------------------------------------------------------------------

path <- ('/domino/datasets/local/RUV')
hkGagnon <- read.csv(paste0(path,"/Genes_Gagnon.txt"), sep="")
hkGagnon <- unlist(hkGagnon)
filter.vst<-FindVariableFeatures(PBC,selection.method = 'vst')
high_varg <- which(filter.vst[hkGagnon,]$vst.variance.standardized>1.8)


# Fake treatment ----------------------------------------------------------

subs <- unique(pbc.metaD$ind_cov)

iter <- subs %in% pbc.metaD[pbc.metaD$Processing_Cohort=='4.0','ind_cov']
set.seed(1)
treatment <- NULL
treatment[iter] <- sample(c("A","B"),sum(iter), replace=T,prob=c(0.9,0.1))
treatment[!iter] <- sample(c("A","B"),sum(!iter), replace=T,prob=c(0.5,0.5))

pbc.metaD <- pbc.metaD %>% left_join(bind_cols(ind_cov = subs, fk.tr = treatment), by='ind_cov')
pdata <- pdata %>% left_join(bind_cols(ind_cov = subs, fk.tr = treatment), by='ind_cov')

rownames(pbc.metaD) <- pbc.metaD$sample_cell
rownames(pdata) <- pdata$cell_id

# Global filter and DGEList -----------------------------------------------------------

PBC <- PBC[rowSums(PBC >= 10) >= 5,]
gene.D <- data.frame(gene=rownames(PBC))
rownames(gene.D) <- gene.D$gene
hk.ind <- rownames(PBC) %in% hkGagnon[-high_varg]

# PBPS -----------------------------------------------------------------

# deffine biological subgroup and group
n=10

groups_info <- group_by(pbc.metaD,fk.tr,cg_cov,pop_cov,Processing_Cohort)%>%
  summarise(nsample=n())%>%
  mutate(bsg=paste(fk.tr,cg_cov,pop_cov,sep='_'),
         group=paste(bsg,Processing_Cohort,sep='_'))

bsg.keep <- names(table(groups_info$bsg)[table(groups_info$bsg)>1]) # remove bsg with only 1 group

pbc.metaD <- mutate(pbc.metaD, bsg=paste(fk.tr,cg_cov,pop_cov,sep='_'),
                    group=paste(bsg,Processing_Cohort,sep='_')) # add bsg and group info

sc.metaD <- pdata %>% mutate(bsg=paste(fk.tr,cg_cov,pop_cov,sep='_'),
                             group=paste(bsg,Processing_Cohort,sep='_')) %>%
  filter(bsg %in% bsg.keep)

cells_to_pbps <- cellspbps(subsample = sc.metaD, seed = 2, n = n)

pbpscounts <- sapply(names(cells_to_pbps),function(x) pbps(cells_to_pbps,Gr=x,ds=ds2s))

pbps_info <- sapply(colnames(pbpscounts), function (x) regmatches(x, regexpr("\\.", x), invert = TRUE))
pa <- sapply(pbps_info,function (x) x[1])
pa_group <- sapply(pbps_info,function (x) x[2])

pbps_info <- data.frame(sample_cell =colnames(pbpscounts), pa = pa, group = pa_group) %>%
  left_join(groups_info, by='group') %>% select(-nsample)
miss_var <- colnames(pbc.metaD)[!colnames(pbc.metaD) %in% colnames(pbps_info) ]
otherv <- sapply(miss_var,function(x) pbps_info[,x ]<- rep(NA,dim(pbps_info)[1]))
pbps_info <- cbind(pbps_info,otherv) %>% mutate(pbps=1)
sifull <- mutate(pbc.metaD,pbps=0,pa='orig') %>% rbind(pbps_info)

sifull$ind_cov <- as.character(sifull$ind_cov)
sifull[sifull$pbps==1,'ind_cov'] <- paste(sifull$pa[sifull$pbps==1] ,sifull$fk.tr[sifull$pbps==1],sifull$pop_cov[sifull$pbps==1] ,sep="_")
sifull$ind_cov_cg_cov <- paste(sifull$ind_cov,sifull$cg_cov,sep='_')

fullcount <- cbind(PBC,pbpscounts[rownames(PBC),])
pbpsds <- list(sifull,fullcount)

write_rds(pbpsds,'/domino/datasets/local/RUV/pbps10rep_ds2.rds')