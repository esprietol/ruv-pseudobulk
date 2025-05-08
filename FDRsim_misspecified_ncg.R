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
library(ggpubr)
library(cluster)
library(factoextra)


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


# NCG ---------------------------------------------------------------------

path <- ('/domino/datasets/local/RUV')
hkGagnon <- read.csv(paste0(path,"/Genes_Gagnon.txt"), sep="")
hkGagnon <- unlist(hkGagnon)
filter.vst<-FindVariableFeatures(PBC,selection.method = 'vst')
high_varg <- which(filter.vst[hkGagnon,]$vst.variance.standardized>1.8)


# Fake treatment ----------------------------------------------------------

subs <- unique(pbc.metaD$ind_cov)
# set.seed(1)
# treatment <- sample(c("A","B"),length(subs), replace=T)

iter <- subs %in% pbc.metaD[pbc.metaD$Processing_Cohort=='4.0','ind_cov']
set.seed(1)
treatment <- NULL
treatment[iter] <- sample(c("A","B"),sum(iter), replace=T,prob=c(0.9,0.1))
treatment[!iter] <- sample(c("A","B"),sum(!iter), replace=T,prob=c(0.5,0.5))

pbc.metaD <- pbc.metaD %>% left_join(bind_cols(ind_cov = subs, fk.tr = treatment), by='ind_cov')
pdata <- pdata %>% left_join(bind_cols(ind_cov = subs, fk.tr = treatment), by='ind_cov')
bc.metaD <- bc.metaD %>% left_join(bind_cols(ind_cov = subs, fk.tr = treatment), by='ind_cov')




# Global filter and DGEList -----------------------------------------------------------

PBC <- PBC[rowSums(PBC >= 10) >= 5,]
gene.D <- data.frame(gene=rownames(PBC))
rownames(gene.D) <- gene.D$gene
hk.ind <- rownames(PBC) %in% hkGagnon[-high_varg]


# hyperparameters -----------------------------------------------------------

# ct = 'T4'
# i=1
seeds <- 1:100*1000
k=5
sw2its <- list()
swpits <- list()
pvalits <- list()
truthits <- list()

for (i in 1:length(seeds)) {
  
  
  
  # In silico DEG -----------------------------------------------------------
  
  
  samp_to_swap <- pbc.metaD$fk.tr == "A"
  
  sim <- PBC
  
  for(j in 1:length(celltypes)){
    set.seed(seeds[i]+j-1)
    sim.a <- simulateDE(PBC[,pbc.metaD$cg_cov==celltypes[j] ], which_cols = samp_to_swap[pbc.metaD$cg_cov==celltypes[j]], prop_DE = 0.1) # I also swap NCG
    gene.D[,  paste0('trueDE_',celltypes[j])] <- sim.a@elementMetadata@listData[["is_DE"]]
    sim.a <- assays(sim.a)$counts
    sim.a <- sim.a[gene.D$gene,]
    sim[,pbc.metaD$cg_cov==celltypes[j] ] <- sim.a    
  }
  
  #Misspecified NCG
  sim.vst <- FindVariableFeatures(sim,selection.method = 'vst')
  hk.ind <- rownames(sim) %in% rownames(slice_min(sim.vst,vst.variance.standardized,n=100))  
  
  
  UQDE <- list()
  UQ_sw2 <- list()
  UQ_swp<- list()
  
  t2ruv2DE <- list()
  t2ruv2_sw2 <- list()
  t2ruv2_swp<- list()
  
  t2ruv3DE <- list()
  t2ruv3_sw2 <- list()
  t2ruv3_swp<- list()
  
  t2ruv3nocDE <- list()
  t2ruv3noc_sw2 <- list()
  t2ruv3noc_swp<- list()
  
  truth <- list()
  
  
  # W Type 2 --------------------------------------------------------------
  
  for (ct in celltypes) {
    
    pbc.ct <- sim[,pbc.metaD$cg_cov==ct]
    pbc.ct <-  pbc.ct[rowSums(pbc.ct >= 10) >= 5,]
    gene.D.ct <- gene.D[ rownames(pbc.ct),]
    gene.D.ct$ind <- paste0(ct,'_',gene.D.ct$gene,'_',i)
    truth[[ct]] <- gene.D.ct
    
    hk.ind <- rownames(pbc.ct) %in% rownames(slice_min(sim.vst,vst.variance.standardized,n=100))
    y.ct <- DGEList(counts=pbc.ct )
    pmin.ct <- find_p(pbc.ct )
    y.ct <- calcNormFactors(y.ct, method="upperquartile", p=pmin.ct)
    nf.ct <- y.ct$samples$norm.factors
    logy.ct <- voom(y.ct, plot = F)$E
    
    # UQ ----------------------------------------------------------------------
    
    pca <- calculatePCA(logy.ct,ncomponents=10)
    pco <- as.integer(pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    UQ_sw2[[ct]] <- mean(sil2[,3])
    UQ_swp[[ct]] <- mean(silp[,3])
    
    design <- model.matrix(~ pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct] ) 
    colnames(design) <- c('(Intercept)','tr')
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    #etable$trueDEG <-  gene.D.ct[rownames(etable),2]
    # etable <- mutate(etable, TP = if_else(DEG+trueDEG==2,T,F), TN =if_else(DEG+trueDEG==0,T,F) )
    # 
    # sum(etable$TP)/sum(etable$trueDEG)
    UQDE[[ct]] <- etable
    
    
    # RUV2 --------------------------------------------------------------
    
    
    ruv2 <- RUV2mod(Y = t(logy.ct), X = pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct], ctl = hk.ind, k=k, Z =pbc.metaD$pop_cov[pbc.metaD$cg_cov==ct])
    
    pca <- calculatePCA(t(ruv2$newY),ncomponents=10)
    pco <- as.integer(pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t2ruv2_sw2[[ct]] <- mean(sil2[,3])
    t2ruv2_swp[[ct]] <- mean(silp[,3])
    
    design <- model.matrix(~ pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct]+ ruv2$W + pbc.metaD$pop_cov[pbc.metaD$cg_cov==ct] ) 
    colnames(design) <- c('(Intercept)','tr',paste0('W', 1:k), 'pop')
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    # etable$trueDEG <- gene.D.ct[rownames(etable),2]
    # 
    # etable <- mutate(etable, TP = if_else(DEG+trueDEG==2,T,F), TN =if_else(DEG+trueDEG==0,T,F) )
    
    # sum(etable$TP)/sum(etable$trueDEG)
    t2ruv2DE[[ct]] <- etable
    
    
    # RUVIII ------------------------------------------------------------------
    
    
    
    Mct <- replicate.matrix(pbc.metaD$ind_cov[pbc.metaD$cg_cov==ct])
    rownames(Mct) <- pbc.metaD$ind_cov[pbc.metaD$cg_cov==ct]
    
    ruv3 <- RUVIIIW(Y = t(logy.ct), M = Mct, ctl=hk.ind,  k = k,
                    return.info = T)
    
    pca <- calculatePCA(t(ruv3$newY),ncomponents=10)
    pco <- as.integer(pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t2ruv3_sw2[[ct]] <- mean(sil2[,3])
    t2ruv3_swp[[ct]] <- mean(silp[,3])
    
    design <- model.matrix(~ pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct]+ ruv3$W +pbc.metaD$pop_cov[pbc.metaD$cg_cov==ct] ) 
    colnames(design) <- c('(Intercept)','tr',paste0('W', 1:k), 'pop')
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct[rownames(etable),2]
    
    t2ruv3DE[[ct]] <- etable
    
    
    ruv3 <- RUVIIIW(Y = t(logy.ct), M = Mct, ctl=rep(T,nrow(logy.ct)),  k = k,
                    return.info = T)
    
    pca <- calculatePCA(t(ruv3$newY),ncomponents=10)
    pco <- as.integer(pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t2ruv3noc_sw2[[ct]] <- mean(sil2[,3])
    t2ruv3noc_swp[[ct]] <- mean(silp[,3])
    
    design <- model.matrix(~ pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct]+ ruv3$W + pbc.metaD$pop_cov[pbc.metaD$cg_cov==ct]) 
    colnames(design) <- c('(Intercept)','tr',paste0('W', 1:k), 'pop')
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct[rownames(etable),2]
    
    t2ruv3nocDE[[ct]] <- etable
    
    
  }
  
  
  pval2 <- NULL 
  
  
  sw2 <- NULL
  swp <- NULL
  
  
  
  pval2[['ind']] <- mapply(function(x,ct) paste0(ct,'_',row.names(x),'_',i),UQDE,names(UQDE),SIMPLIFY = F)
  pval2[['UQ']] <- sapply(UQDE, function(x) x$P.Value, simplify = F)
  pval2[['T22']] <- sapply(t2ruv2DE, function(x) x$P.Value, simplify = F)
  pval2[['T23']] <- sapply(t2ruv3DE, function(x) x$P.Value, simplify = F)
  pval2[['T23noc']] <- sapply(t2ruv3nocDE, function(x) x$P.Value, simplify = F)
  
  sw2[['UQ']] <- UQ_sw2
  sw2[['ruv2T2']] <- t2ruv2_sw2
  sw2[['ruv3T2']] <- t2ruv3_sw2
  sw2[['ruv3T2noc']] <- t2ruv3noc_sw2
  
  swp[['UQ']] <- UQ_swp
  swp[['ruv2T2']] <- t2ruv2_swp
  swp[['ruv3T2']] <- t2ruv3_swp
  swp[['ruv3T2noc']] <- t2ruv3noc_swp
  
  
  sw2its[[i]] <- list_transpose(sw2) #%>% bind_cols()
  swpits[[i]] <- list_transpose(swp)
  
  pvalits[[i]] <- list_transpose(pval2) %>% lapply(bind_cols)
  
  truthits[[i]] <- truth
  
  print(i)
  
}


pvals <- list(pvalits,truthits)
saveRDS(pvals,paste0(path,'/pvals_misspecifiednc_tr_imbalance_popE.rds'))


