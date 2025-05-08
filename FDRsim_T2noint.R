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


# pseudobulk  data -----------------------------------------------------------------

source("/mnt/auxf.R")

PBC <- readRDS("/domino/datasets/local/RUV/pblupusds2.rds")

PBC_split <- split(PBC$samples$sample_cell,PBC$samples$ind_cov_batch_cov)
PBC_S <- DGEList(counts=sapply(PBC_split, function(x) rowSums(PBC$counts[,x])),
                 samples=unique(PBC$samples[,!colnames(PBC$samples) %in% c("cg_cov" ,"sample_cell","lib.size","norm.factors" )]))  

celltypes <- c('B', 'NK', 'T4', 'T8', 'cDC', 'cM', 'ncM', 'pDC')



# NCG ---------------------------------------------------------------------

path <- ('/domino/datasets/local/RUV')
hkGagnon <- read.csv(paste0(path,"/Genes_Gagnon.txt"), sep="")
hkGagnon <- unlist(hkGagnon)
filter.vst <- FindVariableFeatures(PBC,selection.method = 'vst')
top_high_varg <- rownames(arrange(filter.vst,-vst.variance.standardized))[1:500]
high_varg <- which(filter.vst[hkGagnon,]$vst.variance.standardized>1.8)



# Fake treatment ----------------------------------------------------------

subs <- unique(PBC$samples$ind_cov)
set.seed(1)
treatment <- sample(c("A","B"),length(subs), replace=T)

PBC$samples <- PBC$samples %>% left_join(bind_cols(ind_cov = subs, fk.tr = treatment), by='ind_cov')
PBC_S$samples <- PBC_S$samples %>% left_join(bind_cols(ind_cov = subs, fk.tr = treatment), by='ind_cov')




# Global filter and DGEList -----------------------------------------------------------

PBC <- PBC[rowSums(PBC$counts >= 10) >= 5,]
gene.D <- data.frame(gene=rownames(PBC))
rownames(gene.D) <- gene.D$gene



# hyperparameters -----------------------------------------------------------

# ct = 'T4'
# i=1
seeds <- 1:100*1000
k=5
sw2its <- list()
swpits <- list()

sw2itstr <- list()
swpitstr <- list()

pvalitsruv2 <- list()
pvalitsruv3 <- list()
pvalitsruv4 <- list()
truthits <- list()
padjits <- list()
for (i in 1:length(seeds)) {
  
  
  
  # In silico DEG -----------------------------------------------------------
  hk.ind <- rownames(PBC) %in% hkGagnon[-high_varg]
  
  samp_to_swap <- PBC$samples$fk.tr == "A"
  
  sim <- PBC$counts
  
  for(j in 1:length(celltypes)){
    
    set.seed(seeds[i]+j-1)
    sim.a <- simulateDE(PBC$counts[!hk.ind,PBC$samples$cg_cov==celltypes[j] ], which_cols = samp_to_swap[PBC$samples$cg_cov==celltypes[j]], prop_DE = 0.1)
    gene.D[, paste0('trueDE_',celltypes[j])] <- FALSE
    gene.D[!hk.ind,  paste0('trueDE_',celltypes[j])] <- sim.a@elementMetadata@listData[["is_DE"]]
    sim.a <- rbind(assays(sim.a)$counts,PBC$counts[hk.ind,PBC$samples$cg_cov==celltypes[j] ])
    sim.a <- sim.a[gene.D$gene,]
    sim[,PBC$samples$cg_cov==celltypes[j] ] <- sim.a    
  }
  
  sim <- DGEList(counts=sim, samples = PBC$samples)
  sim.t3 <- split(PBC$samples$sample_cell,PBC$samples$ind_cov_batch_cov)
  sim.b <- DGEList(counts=sapply(sim.t3, function(x) rowSums(sim$counts[,x])),
                   samples=unique(PBC$samples[,!colnames(PBC$samples) %in% c("cg_cov" ,"sample_cell","lib.size","norm.factors" )]))  
  
  
  UQDE <- list()
  UQ_sw2 <- list()
  UQ_swp<- list()
  UQ_sw2tr <- list()
  UQ_swptr<- list()
  UQBDE <- list()
  UQB_sw2 <- list()
  UQB_swp<- list()
  UQB_sw2tr <- list()
  UQB_swptr<- list()
  
  t1ruv2DE <- list()
  t1ruv3DE <- list()
  t1ruv4DE <- list()
  t1ruv2_sw2 <- list()
  t1ruv2_swp<- list()
  t1ruv3_sw2 <- list()
  t1ruv3_swp<- list()
  t1ruv4_sw2 <- list()
  t1ruv4_swp<- list()
  t1ruv2_sw2tr <- list()
  t1ruv2_swptr<- list()
  t1ruv3_sw2tr <- list()
  t1ruv3_swptr<- list()
  t1ruv4_sw2tr <- list()
  t1ruv4_swptr<- list()
  
  t2ruv2DE <- list()
  t2ruv3DE <- list()
  t2ruv4DE <- list()
  t2ruv2_sw2 <- list()
  t2ruv2_swp<- list()
  t2ruv3_sw2 <- list()
  t2ruv3_swp<- list()
  t2ruv4_sw2 <- list()
  t2ruv4_swp<- list()
  t2ruv2_sw2tr <- list()
  t2ruv2_swptr<- list()
  t2ruv3_sw2tr <- list()
  t2ruv3_swptr<- list()
  t2ruv4_sw2tr <- list()
  t2ruv4_swptr<- list()
  
  t3ruv2DE <- list()
  t3ruv3DE <- list()
  t3ruv4DE <- list()
  t3ruv2_sw2 <- list()
  t3ruv2_swp<- list()
  t3ruv3_sw2 <- list()
  t3ruv3_swp<- list()
  t3ruv4_sw2 <- list()
  t3ruv4_swp<- list()
  t3ruv2_sw2tr <- list()
  t3ruv2_swptr<- list()
  t3ruv3_sw2tr <- list()
  t3ruv3_swptr<- list()
  t3ruv4_sw2tr <- list()
  t3ruv4_swptr<- list()
  
  truth <- list()
  
  
  # W Type 1 ------------------------------------------------------------------
  
  pbc.all <- sim
  pbc.all <-  pbc.all[rowSums(pbc.all$counts >= 10) >= 5,]
  gene.D.all <- gene.D[ rownames(pbc.all),]
  
  hk.ind <- rownames(pbc.all) %in% hkGagnon[-high_varg]
  #y.all <- DGEList(counts=pbc.all )
  pmin.all <- find_p(pbc.all )
  pbc.all <- calcNormFactors(pbc.all, method="upperquartile", p=pmin.all)
  #nf.all <- pbc.all$samples$norm.factors
  v.all <- voom(pbc.all, plot=T) 
  
  pbc.all$samples$ind_covT1 <- factor(paste(pbc.all$samples$ind_cov, pbc.all$samples$cg_cov,sep="."))
  
  M <- replicate.matrix(pbc.all$samples$ind_covT1)
  rownames(M) <- pbc.all$samples$ind_covT1
  
  trT1 <- factor(paste(pbc.all$samples$fk.tr, pbc.all$samples$cg_cov,sep="."))
  
  ruv2T1 <- RUV2mod(Y = t(v.all$E), X = trT1, ctl = hk.ind, k=k,Z=pbc.all$samples$pop_cov, include.intercept = F )
  ruv3T1 <- RUVIIIW(Y = t(v.all$E), M = M, ctl = hk.ind,  k = k, return.info = T)
  ruv4T1 <- RUV4mod(Y = t(v.all$E), X = trT1, ctl = hk.ind, k = k,Z=pbc.all$samples$pop_cov, include.intercept = F)
  
  
  # W Type 3 ------------------------------------------------------------------
  
  bc <- sim.b
  bc <-  bc[rowSums(bc$counts >= 10) >= 5,]
  gene.D.bc <- gene.D[ rownames(bc),]
  
  hk.ind <- rownames(bc) %in% hkGagnon[-high_varg]
  pmin.bc <- find_p(bc )
  bc <- calcNormFactors(bc, method="upperquartile", p=pmin.bc)
  v.bc <- voom(bc, plot=T)
  
  
  M.bc <- replicate.matrix(bc$samples$ind_cov)
  rownames(M.bc) <- bc$samples$ind_cov
  
  
  
  ruv2T3 <- RUV2mod(Y = t(v.bc$E), X = bc$samples$fk.tr, ctl = hk.ind, k=k,Z=bc$samples$pop_cov)
  ruv3T3 <- RUVIIIW(Y = t(v.bc$E), M = M.bc, ctl = hk.ind,  k = k, return.info = T)
  ruv4T3 <- RUV4mod(Y = t(v.bc$E), X = bc$samples$fk.tr, ctl = hk.ind, k=k,Z=bc$samples$pop_cov)
  
  
  # W Type 2 --------------------------------------------------------------
  
  for (ct in celltypes) {
    indsamp.ct <- bc$samples$ind_cov_batch_cov %in% PBC$samples$ind_cov_batch_cov[PBC$samples$cg_cov==ct]
    pbc.ct <- sim[,PBC$samples$cg_cov==ct]
    pbc.ct <-  pbc.ct[rowSums(pbc.ct$counts >= 10) >= 5,]
    gene.D.ct <- gene.D[ rownames(pbc.ct),]
    gene.D.ct$ind <- paste0(ct,'_',gene.D.ct$gene,'_',i)
    truth[[ct]] <- gene.D.ct
    
    hk.ind <- rownames(pbc.ct) %in% hkGagnon[-high_varg]
    y.ct <- pbc.ct 
    pmin.ct <- find_p(pbc.ct )
    y.ct <- calcNormFactors(y.ct, method="upperquartile", p=pmin.ct)
    v.ct <- voom(y.ct, plot=F)
    
    indg.ct <- rownames(bc) %in% rownames(y.ct)
    indg.ct2 <- rownames(y.ct) %in% rownames(bc) 
    
    # UQ ----------------------------------------------------------------------
    
    top_high_varg.ct <- FindVariableFeatures(v.ct$E)
    top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]
    
    pca <- calculatePCA(v.ct$E,ncomponents=10,subset_row=top_high_varg.ct)
    
    pco <- as.integer(PBC$samples$Processing_Cohort[PBC$samples$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    UQ_sw2[[ct]] <- mean(sil2[,3])
    UQ_swp[[ct]] <- mean(silp[,3])
    
    pco <- as.integer(PBC$samples$fk.tr[PBC$samples$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    UQ_sw2tr[[ct]] <- mean(sil2[,3])
    UQ_swptr[[ct]] <- mean(silp[,3])
    
    design <- model.matrix(~ PBC$samples$fk.tr[PBC$samples$cg_cov==ct] ) 
    colnames(design) <- c('(Intercept)','tr')
    v <- voom(y.ct, design, plot=F) # to recompute weights
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct$trueDE
    
    UQDE[[ct]] <- etable
    
    # UQB ----------------------------------------------------------------------
    
    design <- model.matrix(~ PBC$samples$fk.tr[PBC$samples$cg_cov==ct]+ PBC$samples$pop_cov[PBC$samples$cg_cov==ct]+ PBC$samples$Processing_Cohort[PBC$samples$cg_cov==ct] ) 
    colnames(design) <- c('(Intercept)','tr','pop',paste0('Proc', 2:4))
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    alpha <- vfit$coefficients[,4:6] # same coefficients as efit and returned by toptable
    newY <- t(v$E) - design[,4:6]%*%t(alpha)
    
    top_high_varg.ct <- FindVariableFeatures(t(newY))
    top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]
    pca <- calculatePCA(t(newY),ncomponents=10, subset_row=top_high_varg.ct)
    pco <- as.integer(PBC$samples$Processing_Cohort[PBC$samples$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    UQB_sw2[[ct]] <- mean(sil2[,3])
    UQB_swp[[ct]] <- mean(silp[,3])
    
    pco <- as.integer(PBC$samples$fk.tr[PBC$samples$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    UQB_sw2tr[[ct]] <- mean(sil2[,3])
    UQB_swptr[[ct]] <- mean(silp[,3])
    
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct$trueDE
    
    UQBDE[[ct]] <- etable
    
    
    # RUV2 --------------------------------------------------------------
    top_high_varg.ct <- FindVariableFeatures(t(ruv2T1$newY[PBC$samples$cg_cov==ct,]))
    top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]
    pca <- calculatePCA(t(ruv2T1$newY[PBC$samples$cg_cov==ct,]),ncomponents=10,subset_row=top_high_varg.ct)
    pco <- as.integer(PBC$samples$Processing_Cohort[PBC$samples$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t1ruv2_sw2[[ct]] <- mean(sil2[,3])
    t1ruv2_swp[[ct]] <- mean(silp[,3])
    
    pco <- as.integer(PBC$samples$fk.tr[PBC$samples$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t1ruv2_sw2tr[[ct]] <- mean(sil2[,3])
    t1ruv2_swptr[[ct]] <- mean(silp[,3])
    
    
    design <- model.matrix(~ PBC$samples$fk.tr[PBC$samples$cg_cov==ct]+PBC$samples$pop_cov[PBC$samples$cg_cov==ct]+ ruv2T1$W[PBC$samples$cg_cov==ct,] ) 
    colnames(design) <- c('(Intercept)','tr','pop',paste0('W', 1:k))
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct[rownames(etable),2]
    
    t1ruv2DE[[ct]] <- etable
    
    
    ruv2 <- RUV2mod(Y = t(v.ct$E), X = PBC$samples$fk.tr[PBC$samples$cg_cov==ct], ctl = hk.ind, k=k, Z=PBC$samples$pop_cov[PBC$samples$cg_cov==ct])
    
    top_high_varg.ct <- FindVariableFeatures(t(ruv2$newY))
    top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]
    pca <- calculatePCA(t(ruv2$newY),ncomponents=10,subset_row=top_high_varg.ct )
    pco <- as.integer(PBC$samples$Processing_Cohort[PBC$samples$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t2ruv2_sw2[[ct]] <- mean(sil2[,3])
    t2ruv2_swp[[ct]] <- mean(silp[,3])
    
    pco <- as.integer(PBC$samples$fk.tr[PBC$samples$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t2ruv2_sw2tr[[ct]] <- mean(sil2[,3])
    t2ruv2_swptr[[ct]] <- mean(silp[,3])
    
    design <- model.matrix(~ PBC$samples$fk.tr[PBC$samples$cg_cov==ct]+ ruv2$W +PBC$samples$pop_cov[PBC$samples$cg_cov==ct]) 
    colnames(design) <- c('(Intercept)','tr',paste0('W', 1:k),'pop')
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct[rownames(etable),2]
    
    t2ruv2DE[[ct]] <- etable
    
    
    newY <- t(v.ct$E[indg.ct2,]) - ruv2T3$W[indsamp.ct,] %*% ruv2T3$fullalpha[,indg.ct]
    top_high_varg.ct <- FindVariableFeatures(t(newY))
    top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]
    pca <- calculatePCA(t(newY),ncomponents=10,subset_row=top_high_varg.ct)
    pco <- as.integer(PBC$samples$Processing_Cohort[PBC$samples$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t3ruv2_sw2[[ct]] <- mean(sil2[,3])
    t3ruv2_swp[[ct]] <- mean(silp[,3])
    
    pco <- as.integer(PBC$samples$fk.tr[PBC$samples$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t3ruv2_sw2tr[[ct]] <- mean(sil2[,3])
    t3ruv2_swptr[[ct]] <- mean(silp[,3])
    
    design <- model.matrix(~ PBC$samples$fk.tr[PBC$samples$cg_cov==ct]+PBC$samples$pop_cov[PBC$samples$cg_cov==ct]+ ruv2T3$W[indsamp.ct,] ) 
    colnames(design) <- c('(Intercept)','tr','pop',paste0('W', 1:k))
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct[rownames(etable),2]
    
    t3ruv2DE[[ct]] <- etable
    
    
    
    
    # RUV3 --------------------------------------------------------------------
    top_high_varg.ct <- FindVariableFeatures(t(ruv3T1$newY[PBC$samples$cg_cov==ct,]))
    top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]
    pca <- calculatePCA(t(ruv3T1$newY[PBC$samples$cg_cov==ct,]),ncomponents=10, subset_row=top_high_varg.ct)
    pco <- as.integer(PBC$samples$Processing_Cohort[PBC$samples$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t1ruv3_sw2[[ct]] <- mean(sil2[,3])
    t1ruv3_swp[[ct]] <- mean(silp[,3])
    
    pco <- as.integer(PBC$samples$fk.tr[PBC$samples$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t1ruv3_sw2tr[[ct]] <- mean(sil2[,3])
    t1ruv3_swptr[[ct]] <- mean(silp[,3])
    
    design <- model.matrix(~ PBC$samples$fk.tr[PBC$samples$cg_cov==ct]+PBC$samples$pop_cov[PBC$samples$cg_cov==ct]+ ruv3T1$W[PBC$samples$cg_cov==ct,] ) 
    colnames(design) <- c('(Intercept)','tr','pop',paste0('W', 1:k))
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct[rownames(etable),2]
    
    t1ruv3DE[[ct]] <- etable
    
    
    Mct <- replicate.matrix(PBC$samples$ind_cov[PBC$samples$cg_cov==ct])
    rownames(Mct) <- PBC$samples$ind_cov[PBC$samples$cg_cov==ct]
    
    ruv3 <- RUVIIIW(Y = t(v.ct$E), M = Mct, ctl=hk.ind,  k = k,
                    return.info = T)
    
    top_high_varg.ct <- FindVariableFeatures(t(ruv3$newY))
    top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]
    pca <- calculatePCA(t(ruv3$newY),ncomponents=10)
    pco <- as.integer(PBC$samples$Processing_Cohort[PBC$samples$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t2ruv3_sw2[[ct]] <- mean(sil2[,3])
    t2ruv3_swp[[ct]] <- mean(silp[,3])
    
    pco <- as.integer(PBC$samples$fk.tr[PBC$samples$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t2ruv3_sw2tr[[ct]] <- mean(sil2[,3])
    t2ruv3_swptr[[ct]] <- mean(silp[,3])
    
    design <- model.matrix(~ PBC$samples$fk.tr[PBC$samples$cg_cov==ct]+PBC$samples$pop_cov[PBC$samples$cg_cov==ct]+ ruv3$W ) 
    colnames(design) <- c('(Intercept)','tr','pop',paste0('W', 1:k))
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct[rownames(etable),2]
    
    t2ruv3DE[[ct]] <- etable
    
    newY <- t(v.ct$E[indg.ct2,]) - ruv3T3$W[indsamp.ct,] %*% ruv3T3$fullalpha[,indg.ct]
    top_high_varg.ct <- FindVariableFeatures(t(newY))
    top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]
    pca <- calculatePCA(t(newY),ncomponents=10)
    pco <- as.integer(PBC$samples$Processing_Cohort[PBC$samples$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t3ruv3_sw2[[ct]] <- mean(sil2[,3])
    t3ruv3_swp[[ct]] <- mean(silp[,3])
    
    pco <- as.integer(PBC$samples$fk.tr[PBC$samples$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t3ruv3_sw2tr[[ct]] <- mean(sil2[,3])
    t3ruv3_swptr[[ct]] <- mean(silp[,3])
    
    design <- model.matrix(~ PBC$samples$fk.tr[PBC$samples$cg_cov==ct]+PBC$samples$pop_cov[PBC$samples$cg_cov==ct]+ ruv3T3$W[indsamp.ct,] ) 
    colnames(design) <- c('(Intercept)','tr','pop',paste0('W', 1:k))
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct[rownames(etable),2]
    
    t3ruv3DE[[ct]] <- etable
    
    # RUV4 --------------------------------------------------------------
    top_high_varg.ct <- FindVariableFeatures(t(ruv4T1$newY[PBC$samples$cg_cov==ct,]))
    top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]
    pca <- calculatePCA(t(ruv4T1$newY[PBC$samples$cg_cov==ct,]),ncomponents=10, subset_row=top_high_varg.ct)
    pco <- as.integer(PBC$samples$Processing_Cohort[PBC$samples$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t1ruv4_sw2[[ct]] <- mean(sil2[,3])
    t1ruv4_swp[[ct]] <- mean(silp[,3])
    
    pco <- as.integer(PBC$samples$fk.tr[PBC$samples$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t1ruv4_sw2tr[[ct]] <- mean(sil2[,3])
    t1ruv4_swptr[[ct]] <- mean(silp[,3])
    
    design <- model.matrix(~ PBC$samples$fk.tr[PBC$samples$cg_cov==ct]+PBC$samples$pop_cov[PBC$samples$cg_cov==ct]+ ruv4T1$W[PBC$samples$cg_cov==ct,] ) 
    colnames(design) <- c('(Intercept)','tr','pop',paste0('W', 1:k))
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct[rownames(etable),2]
    
    t1ruv4DE[[ct]] <- etable
    
    
    ruv4 <- RUV4mod(Y = t(v.ct$E), X = PBC$samples$fk.tr[PBC$samples$cg_cov==ct], ctl = hk.ind, k=k, Z=PBC$samples$pop_cov[PBC$samples$cg_cov==ct])
    top_high_varg.ct <- FindVariableFeatures(t(ruv4$newY))
    top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]
    pca <- calculatePCA(t(ruv4$newY),ncomponents=10,subset_row=top_high_varg.ct)
    pco <- as.integer(PBC$samples$Processing_Cohort[PBC$samples$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t2ruv4_sw2[[ct]] <- mean(sil2[,3])
    t2ruv4_swp[[ct]] <- mean(silp[,3])
    
    pco <- as.integer(PBC$samples$fk.tr[PBC$samples$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t2ruv4_sw2tr[[ct]] <- mean(sil2[,3])
    t2ruv4_swptr[[ct]] <- mean(silp[,3])
    
    design <- model.matrix(~ PBC$samples$fk.tr[PBC$samples$cg_cov==ct]+PBC$samples$pop_cov[PBC$samples$cg_cov==ct]+ ruv4$W ) 
    colnames(design) <- c('(Intercept)','tr','pop',paste0('W', 1:k))
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct[rownames(etable),2]
    
    t2ruv4DE[[ct]] <- etable
    
    newY <- t(v.ct$E[indg.ct2,]) - ruv4T3$W[indsamp.ct,] %*% ruv4T3$fullalpha[,indg.ct]
    top_high_varg.ct <- FindVariableFeatures(t(newY))
    top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]
    pca <- calculatePCA(t(newY),ncomponents=10, subset_row=top_high_varg.ct)
    pco <- as.integer(PBC$samples$Processing_Cohort[PBC$samples$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t3ruv4_sw2[[ct]] <- mean(sil2[,3])
    t3ruv4_swp[[ct]] <- mean(silp[,3])
    
    pco <- as.integer(PBC$samples$fk.tr[PBC$samples$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t3ruv4_sw2tr[[ct]] <- mean(sil2[,3])
    t3ruv4_swptr[[ct]] <- mean(silp[,3])
    
    design <- model.matrix(~ PBC$samples$fk.tr[PBC$samples$cg_cov==ct]+PBC$samples$pop_cov[PBC$samples$cg_cov==ct]+ ruv4T3$W[indsamp.ct,] ) 
    colnames(design) <- c('(Intercept)','tr','pop',paste0('W', 1:k))
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct[rownames(etable),2]
    
    t3ruv4DE[[ct]] <- etable
    
    
  }
  
  
  pval2 <- NULL 
  pval3 <- NULL
  pval4 <- NULL
  sw2 <- NULL
  swp <- NULL
  
  
  
  pval2[['ind']] <- mapply(function(x,ct) paste0(ct,'_',row.names(x),'_',i),UQDE,names(UQDE),SIMPLIFY = F)
  pval2[['UQ']] <- sapply(UQDE, function(x) x$P.Value, simplify = F)
  pval2[['UQ_Batch']] <- sapply(UQBDE, function(x) x$P.Value, simplify = F)
  pval2[['T1']] <- sapply(t1ruv2DE, function(x) x$P.Value, simplify = F)
  pval2[['T2']] <- sapply(t2ruv2DE, function(x) x$P.Value, simplify = F)
  pval2[['T3']] <- sapply(t3ruv2DE, function(x) x$P.Value, simplify = F)
  
  pval3[['ind']] <- mapply(function(x,ct) paste0(ct,'_',row.names(x),'_',i),UQDE,names(UQDE),SIMPLIFY = F)
  pval3[['UQ']] <- sapply(UQDE, function(x) x$P.Value, simplify = F)
  pval3[['UQ_Batch']] <- sapply(UQBDE, function(x) x$P.Value, simplify = F)
  pval3[['T1']] <- sapply(t1ruv3DE, function(x) x$P.Value, simplify = F)
  pval3[['T2']] <- sapply(t2ruv3DE, function(x) x$P.Value, simplify = F)
  pval3[['T3']] <- sapply(t3ruv3DE, function(x) x$P.Value, simplify = F)
  
  
  pval4[['ind']] <- mapply(function(x,ct) paste0(ct,'_',row.names(x),'_',i),UQDE,names(UQDE),SIMPLIFY = F)
  pval4[['UQ']] <- sapply(UQDE, function(x) x$P.Value, simplify = F)
  pval4[['UQ_Batch']] <- sapply(UQBDE, function(x) x$P.Value, simplify = F)
  pval4[['T1']] <- sapply(t1ruv4DE, function(x) x$P.Value, simplify = F)  
  pval4[['T2']] <- sapply(t2ruv4DE, function(x) x$P.Value, simplify = F)
  pval4[['T3']] <- sapply(t3ruv4DE, function(x) x$P.Value, simplify = F)
  
  sw2[['UQ']] <- UQ_sw2
  sw2[['UQ_Batch']] <- UQB_sw2
  sw2[['ruv2T1']] <- t1ruv2_sw2
  sw2[['ruv2T2']] <- t2ruv2_sw2
  sw2[['ruv2T3']] <- t3ruv2_sw2
  sw2[['ruv3T1']] <- t1ruv3_sw2
  sw2[['ruv3T2']] <- t2ruv3_sw2
  sw2[['ruv3T3']] <- t3ruv3_sw2
  sw2[['ruv4T1']] <- t1ruv4_sw2
  sw2[['ruv4T2']] <- t2ruv4_sw2
  sw2[['ruv4T3']] <- t3ruv4_sw2
  
  swp[['UQ']] <- UQ_swp
  swp[['UQ_Batch']] <- UQB_swp
  swp[['ruv2T1']] <- t1ruv2_swp
  swp[['ruv2T2']] <- t2ruv2_swp
  swp[['ruv2T3']] <- t3ruv2_swp
  swp[['ruv3T1']] <- t1ruv3_swp
  swp[['ruv3T2']] <- t2ruv3_swp
  swp[['ruv3T3']] <- t3ruv3_swp
  swp[['ruv4T1']] <- t1ruv4_swp
  swp[['ruv4T2']] <- t2ruv4_swp
  swp[['ruv4T3']] <- t3ruv4_swp
  
  sw2its[[i]] <- list_transpose(sw2) #%>% bind_cols()
  swpits[[i]] <- list_transpose(swp)
  
  sw2 <- NULL
  swp <- NULL
  
  sw2[['UQ']] <- UQ_sw2tr
  sw2[['UQ_Batch']] <- UQB_sw2tr
  sw2[['ruv2T1']] <- t1ruv2_sw2tr
  sw2[['ruv2T2']] <- t2ruv2_sw2tr
  sw2[['ruv2T3']] <- t3ruv2_sw2tr
  sw2[['ruv3T1']] <- t1ruv3_sw2tr
  sw2[['ruv3T2']] <- t2ruv3_sw2tr
  sw2[['ruv3T3']] <- t3ruv3_sw2tr
  sw2[['ruv4T1']] <- t1ruv4_sw2tr
  sw2[['ruv4T2']] <- t2ruv4_sw2tr
  sw2[['ruv4T3']] <- t3ruv4_sw2tr
  
  swp[['UQ']] <- UQ_swptr
  swp[['UQ_Batch']] <- UQB_swptr
  swp[['ruv2T1']] <- t1ruv2_swptr
  swp[['ruv2T2']] <- t2ruv2_swptr
  swp[['ruv2T3']] <- t3ruv2_swptr
  swp[['ruv3T1']] <- t1ruv3_swptr
  swp[['ruv3T2']] <- t2ruv3_swptr
  swp[['ruv3T3']] <- t3ruv3_swptr
  swp[['ruv4T1']] <- t1ruv4_swptr
  swp[['ruv4T2']] <- t2ruv4_swptr
  swp[['ruv4T3']] <- t3ruv4_swptr
  
  sw2itstr[[i]] <- list_transpose(sw2) #%>% bind_cols()
  swpitstr[[i]] <- list_transpose(swp)
  
  pvalitsruv2[[i]] <- list_transpose(pval2) %>% lapply(bind_cols)
  pvalitsruv3[[i]] <- list_transpose(pval3) %>% lapply(bind_cols)
  pvalitsruv4[[i]] <- list_transpose(pval4) %>% lapply(bind_cols)
  
  truthits[[i]] <- truth
  
  print(i)
  
}


pvals <- list(pvalitsruv2,pvalitsruv3,pvalitsruv4,truthits)
saveRDS(pvals,paste0(path,'/pvals_ruvdg_diffg_pop_E.rds'))

saveRDS(sw2its,paste0(path,'/asw2_ruvdg_diffg_pop_E.rds'))
saveRDS(swpits,paste0(path,'/aswp_ruvdg_diffg_pop_E.rds'))

