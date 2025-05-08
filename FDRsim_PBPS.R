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

PBC <- PBC[rowSums(PBC >= 5) >= 5,]
gene.D <- data.frame(gene=rownames(PBC))
rownames(gene.D) <- gene.D$gene
hk.ind <- rownames(PBC) %in% hkGagnon[-high_varg]

# PBPR -----------------------------------------------------------------

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

# pbps2 <- readRDS('/domino/datasets/local/RUV/pbps10rep_ds2.rds')
# 
# fullcount <- pbps2$counts
# sifull <- pbps2$samples


# graphs ------------------------------------------------------------------
path2 <- ('/domino/datasets/local/RUV/jul24')

group_by(sifull,Processing_Cohort,ind_cov,pop_cov,Age,fk.tr) %>%
  summarise(count=n()) %>% 
  ggplot(aes(y=ind_cov,x=Processing_Cohort, group=ind_cov, color=Age, shape=pop_cov)) +
  geom_point(size=2) + geom_line() + scale_color_viridis_b() +
  ggtitle('Replicates and Pseudo Replicates')

ggsave(paste0(path2,'pbpsreps.png'))

# hyperparameters -----------------------------------------------------------

# ct = 'T4'
# i=1
seeds <- 1:100*1000
k=5
pvalits <- list()
truthits <- list()
aswits <- list()
for (i in 1:length(seeds)) {
  
  # In silico DEG PS -----------------------------------------------------------
  
  hk.ind <- rownames(PBC) %in% hkGagnon[-high_varg]
  samp_to_swap <- sifull$fk.tr == "A"
  sim.ps <- fullcount
  
  for(j in 1:length(celltypes)){
    set.seed(seeds[i]+j-1)
    sim.ct <- simulateDE(fullcount[!hk.ind, sifull$cg_cov==celltypes[j]], which_cols = samp_to_swap[sifull$cg_cov==celltypes[j]], prop_DE = 0.1)
    gene.D[, paste0('trueDE_',celltypes[j])] <- FALSE
    gene.D[!hk.ind,  paste0('trueDE_',celltypes[j])] <- sim.ct@elementMetadata@listData[["is_DE"]]
    sim.ct <- rbind(assays(sim.ct)$counts,fullcount[hk.ind, sifull$cg_cov==celltypes[j]])
    sim.ct <- sim.ct[gene.D$gene,]
    sim.ps[,sifull$cg_cov==celltypes[j] ] <- sim.ct
  }
  
  
  
  # In silico DEG -----------------------------------------------------------
  
  sim <- sim.ps[,colnames(PBC)]
  
  UQBDE <- list()
  ruv3DE <- list()
  ruv3BDE <- list()
  ruv3psDE <- list()
  ruv3psallDE <- list()
  UQB_asw <- list()
  ruv3_asw <- list()
  ruv3B_asw <- list()
  ruv3ps_asw <- list()
  ruv3psall_asw <- list()
  truth <- list()
  
  # ct deglist --------------------------------------------------------------
  for (ct in celltypes) {
    
    pbc.ct <- sim[,pbc.metaD$cg_cov==ct]
    pbc.ct <-  pbc.ct[rowSums(pbc.ct >= 10) >= 5,]
    gene.D.ct <- gene.D[ rownames(pbc.ct),]
    gene.D.ct$ind <- paste0(ct,'_',gene.D.ct$gene,'_',i)
    #rownames(gene.D.ct) <-NULL
    truth[[ct]] <- gene.D.ct
    #rownames(gene.D.ct) <- gene.D.ct$gene
    
    hk.ind <- rownames(pbc.ct) %in% hkGagnon[-high_varg]
    y.ct <- DGEList(counts=pbc.ct )
    pmin.ct <- find_p(pbc.ct )
    y.ct <- calcNormFactors(y.ct, method="upperquartile", p=pmin.ct)
    #nf.ct <- y.ct$samples$norm.factors
    logy.ct <- voom(y.ct,plot= F)$E
    
    # UQB ----------------------------------------------------------------------
    
    design <- model.matrix(~ pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct]+ 
                             pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct]+ 
                             pbc.metaD$pop_cov[pbc.metaD$cg_cov==ct] ) 
    colnames(design) <- c('(Intercept)','tr',paste0('Proc', 2:4),'pop')
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    
    alphaUQB <- vfit$coefficients[,3:5] # same coefficients as efit and returned by toptable
    newY <- t(logy.ct) - design[,3:5]%*%t(alphaUQB)
    
    top_high_varg.ct <- FindVariableFeatures(t(newY))
    top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]
    pca <- calculatePCA(t(newY),ncomponents=10, subset_row=top_high_varg.ct)
    
    aswPC <- ASW(pca,pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
    aswTr<- ASW(pca,pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct])
    aswEt<- ASW(pca,pbc.metaD$pop_cov[pbc.metaD$cg_cov==ct])
    
    UQB_asw[[ct]] <- rbind(aswPC,aswTr,aswEt)
    
    UQBDE[[ct]] <- limma.matrix(y=y.ct, design=design, gene.D=gene.D.ct)
    
    
    # RUV3 --------------------------------------------------------------------
    
    Mct <- replicate.matrix(pbc.metaD$ind_cov[pbc.metaD$cg_cov==ct])
    rownames(Mct) <- pbc.metaD$ind_cov[pbc.metaD$cg_cov==ct]
    
    #ruv3 <- RUVIIIW(Y = t(logy.ct), M = Mct, ctl=rep(T,nrow(logy.ct)),  k = k,
    #               return.info = T)
    
    ruv3 <- RUVIIIW(Y = t(logy.ct), M = Mct, ctl = hk.ind,  k = k,
                    return.info = T)
    
    design <- model.matrix(~ pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct]+ ruv3$W + pbc.metaD$pop_cov[pbc.metaD$cg_cov==ct] ) 
    colnames(design) <- c('(Intercept)','tr',paste0('W', 1:k),'pop')
    
    top_high_varg.ct <- FindVariableFeatures(t(ruv3$newY))
    top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]
    pca <- calculatePCA(t(ruv3$newY),ncomponents=10, subset_row=top_high_varg.ct)
    aswPC <- ASW(pca,pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
    aswTr<- ASW(pca,pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct])
    aswEt<- ASW(pca,pbc.metaD$pop_cov[pbc.metaD$cg_cov==ct])
    
    ruv3_asw[[ct]] <- rbind(aswPC,aswTr,aswEt)
    
    ruv3DE[[ct]] <- limma.matrix(y=y.ct, design=design, gene.D=gene.D.ct)
    
    # RUV3 + Batch --------------------------------------------------------------------
    
    
    design <- model.matrix(~ pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct]+ ruv3$W + 
                             pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct]+
                             pbc.metaD$pop_cov[pbc.metaD$cg_cov==ct]) 
    colnames(design) <- c('(Intercept)','tr',paste0('W', 1:k), paste0('PCohort', 2:4),'pop')
    
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    
    alpharuv3 <- vfit$coefficients[,3:10] # same coefficients as efit and returned by toptable
    newY <- t(logy.ct) - design[,3:10]%*%t(alpharuv3)
    
    top_high_varg.ct <- FindVariableFeatures(t(newY))
    top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]
    pca <- calculatePCA(t(newY),ncomponents=10, subset_row=top_high_varg.ct)
    aswPC <- ASW(pca,pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
    aswTr<- ASW(pca,pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct])
    aswEt<- ASW(pca,pbc.metaD$pop_cov[pbc.metaD$cg_cov==ct])
    
    ruv3B_asw[[ct]] <- rbind(aswPC,aswTr,aswEt)
    
    ruv3BDE[[ct]] <- limma.matrix(y=y.ct, design=design, gene.D=gene.D.ct)
    
    # ct ps deglist --------------------------------------------------------------
    
    ps.pbc.ct <- sim.ps[,sifull$cg_cov==ct]
    ps.pbc.ct <-  ps.pbc.ct[rownames(pbc.ct),]
    ps.y.ct <- DGEList(counts=ps.pbc.ct )
    pmin.ct <- find_p(ps.pbc.ct )
    ps.y.ct <- calcNormFactors(ps.y.ct, method="upperquartile", p=pmin.ct)
    # ps.nf.ct <- ps.y.ct$samples$norm.factors
    ps.logy.ct <- voom(ps.y.ct,plot=F)$E
    
    
    # RUV3 PS ---------------------------------------------------------------
    
    ps.M.ct <- replicate.matrix(sifull$ind_cov_cg_cov[sifull$cg_cov==ct])
    rownames(ps.M.ct) <- sifull$ind_cov_cg_cov[sifull$cg_cov==ct]
    
    ruv3ps <- RUVIIIW(Y = t(ps.logy.ct), M = ps.M.ct, ctl = hk.ind,  k = k,
                      return.info = T)
    
    # ruv3ps <- RUVIIIW(Y = t(ps.logy.ct), M = ps.M.ct, ctl = rep(T,nrow(logy.ct)),  k = k,
    #                   return.info = T)
    # 
    orig.s <- sifull$sample_cell[sifull$cg_cov==ct & sifull$pbps==0]
    design <- model.matrix(~ sifull$fk.tr[sifull$sample_cell%in%orig.s]+ ruv3ps$W[orig.s,]+ sifull$pop_cov[sifull$sample_cell%in%orig.s]) 
    colnames(design) <-c('(Intercept)','tr',paste0('W', 1:k),'pop')
    
    top_high_varg.ct <- FindVariableFeatures(t(ruv3ps$newY[orig.s,]))
    top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]
    pca <- calculatePCA(t(ruv3ps$newY[orig.s,]),ncomponents=10, subset_row=top_high_varg.ct)
    aswPC <- ASW(pca,pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
    aswTr<- ASW(pca,pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct])
    aswEt<- ASW(pca,pbc.metaD$pop_cov[pbc.metaD$cg_cov==ct])
    
    ruv3ps_asw[[ct]] <- rbind(aswPC,aswTr,aswEt)
    
    ruv3psDE[[ct]] <- limma.matrix(y=ps.y.ct[,orig.s], design=design, gene.D=gene.D.ct)
    
    # RUV3 PS all---------------------------------------------------------------
    
    ps.M.ct <- replicate.matrix(sifull$ind_cov_cg_cov[sifull$cg_cov==ct])
    rownames(ps.M.ct) <- sifull$ind_cov_cg_cov[sifull$cg_cov==ct]
    
    ruv3ps <- RUVIIIW(Y = t(ps.logy.ct), M = ps.M.ct, ctl = rep(T,nrow(ps.logy.ct)),  k = k,
                      return.info = T)
    
    # ruv3ps <- RUVIIIW(Y = t(ps.logy.ct), M = ps.M.ct, ctl = rep(T,nrow(logy.ct)),  k = k,
    #                   return.info = T)
    # 
    orig.s <- sifull$sample_cell[sifull$cg_cov==ct & sifull$pbps==0]
    design <- model.matrix(~ sifull$fk.tr[sifull$sample_cell%in%orig.s]+ ruv3ps$W[orig.s,]+ sifull$pop_cov[sifull$sample_cell%in%orig.s]) 
    colnames(design) <-c('(Intercept)','tr',paste0('W', 1:k),'pop')
    
    top_high_varg.ct <- FindVariableFeatures(t(ruv3ps$newY[orig.s,]))
    top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]
    pca <- calculatePCA(t(ruv3ps$newY[orig.s,]),ncomponents=10, subset_row=top_high_varg.ct)
    aswPC <- ASW(pca,pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
    aswTr<- ASW(pca,pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct])
    aswEt<- ASW(pca,pbc.metaD$pop_cov[pbc.metaD$cg_cov==ct])
    
    ruv3psall_asw[[ct]] <- rbind(aswPC,aswTr,aswEt)
    
    ruv3psallDE[[ct]] <- limma.matrix(y=ps.y.ct[,orig.s], design=design, gene.D=gene.D.ct)
    
    
  }
  
  # factors are the same bc computed by sample? no, other samples affect normalization but differences are minimal
  
  
  pval <- NULL 
  pval[['ind']] <- mapply(function(x,ct) paste0(ct,'_',row.names(x),'_',i),ruv3DE,names(ruv3DE),SIMPLIFY = F)
  pval[['UQ_Batch']] <- sapply(UQBDE, function(x) x$P.Value, simplify = F)
  pval[['RUV3']] <- sapply(ruv3DE, function(x) x$P.Value, simplify = F)  
  pval[['RUV3B']] <- sapply(ruv3BDE, function(x) x$P.Value, simplify = F)
  pval[['RUV3PS']] <- sapply(ruv3psDE, function(x) x$P.Value, simplify = F)
  pval[['RUV3PSall']] <- sapply(ruv3psallDE, function(x) x$P.Value, simplify = F)
  
  pval <- list_transpose(pval) %>% lapply(bind_cols)
  
  asw <- NULL
  asw[['UQ_Batch']] <- UQB_asw
  asw[['RUV3']] <- ruv3_asw
  asw[['RUV3B']] <- ruv3B_asw
  asw[['RUV3PS']] <- ruv3ps_asw
  asw[['RUV3PSall']] <- ruv3psall_asw
  
  asw2 <- list_transpose(asw) %>% lapply(bind_cols)
  
  pvalits[[i]] <- pval
  aswits[[i]] <- asw
  truthits[[i]] <- truth
  
  print(paste('it ',i))
  
}
saveRDS(truthits,paste0(path,'/truthpbps10_dgim_popE.rds'))
saveRDS(pvalits,paste0(path,'/pvalspbps10_dgim_popE.rds'))
saveRDS(aswits,paste0(path,'/aswpbps10_dgim_popE.rds'))

