rm(list=ls())
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
library(aricode)
library(cluster)
library(factoextra)
library(mvrsquared)

# Previously all cells with more than 100 genes, all cells with more than 400 counts, cellls with <15 percentage of mitochondrial genes


# Descriptives paper ------------------------------------------------------

# sc data -----------------------------------------------------------------

source("/mnt/auxf.R")
ds2 <- readRDS("/domino/datasets/local/RUV/sclupusds2.rds")
ds2 <- ds2[,ds2$Sex=='Female']
Idents(object = ds2) <- "cg_cov"
celltypes <- c('B', 'NK', 'T4', 'T8', 'cDC', 'cM', 'ncM', 'pDC')
ds2@meta.data <- droplevels(ds2@meta.data)
ds2s <- subset(x = ds2, idents = celltypes)
rm(list = 'ds2')

ds2s@meta.data <- droplevels(ds2s@meta.data)
pdata <- ds2s@meta.data
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


# NCG ---------------------------------------------------------------------

path <- ('/domino/datasets/local/RUV')
hkGagnon <- read.csv(paste0(path,"/Genes_Gagnon.txt"), sep="")
hkGagnon <- unlist(hkGagnon)
filter.vst<-FindVariableFeatures(PBC,selection.method = 'vst')
high_varg <- which(filter.vst[hkGagnon,]$vst.variance.standardized>1.8)


# Fake treatment ----------------------------------------------------------

subs <- unique(pbc.metaD$ind_cov)
set.seed(1)
treatment <- sample(c("A","B"),length(subs), replace=T)

pbc.metaD <- pbc.metaD %>% left_join(bind_cols(ind_cov = subs, fk.tr = treatment), by='ind_cov')
pdata <- pdata %>% left_join(bind_cols(ind_cov = subs, fk.tr = treatment), by='ind_cov')
bc.metaD <- unique(pbc.metaD[,!colnames(pbc.metaD) %in% c("cg_cov" ,"sample_cell" )])


# Global filter and DGEList -----------------------------------------------------------

PBC <- PBC[rowSums(PBC >= 10) >= 5,]
gene.D <- data.frame(gene=rownames(PBC))
rownames(gene.D) <- gene.D$gene
hk.ind <- rownames(PBC) %in% hkGagnon[-high_varg]

pbc.t3 <- split(pbc.metaD$sample_cell,pbc.metaD$ind_cov_batch_cov)
pbc.b <- sapply(pbc.t3, function(x) rowSums(PBC[,x])  )

# Graphs ------------------------------------------------------------------
path2 <- '/domino/datasets/local/RUV/paper/graphs/normalisation/'

group_by(pbc.metaD,Processing_Cohort,ind_cov,pop_cov,Age,fk.tr) %>%
  summarise(count=n()) %>% ungroup()%>% 
  mutate(ind_cov=str_extract(as.character(ind_cov),"[^_]+$")) %>%
  ggplot(aes(y=ind_cov,x=Processing_Cohort, group=ind_cov, color=pop_cov)) +
  geom_point(size=2) + geom_line() + scale_colour_brewer(palette='Set1') +
  labs(y = "Subject", x= "Processing cohort", color = "Ethnicity")+ theme_minimal()+
  ggtitle('Assays')

ggsave(paste0(path,'/paper/graphs/37reps.png'),height = 5, width = 7)


k=5
ct='T4'

orddata <- arrange(pbc.metaD,Processing_Cohort,cg_cov) 
orddata.bc <- arrange(bc.metaD,Processing_Cohort) 

# Type1 -------------------------------------------------------------------


pbc.all <- PBC
pbc.all <-  pbc.all[rowSums(pbc.all >= 10) >= 5,]
gene.D.all <- gene.D[ rownames(pbc.all),]

hk.ind <- rownames(pbc.all) %in% hkGagnon[-high_varg]
y.all <- DGEList(counts=pbc.all )
pmin.all <- find_p(pbc.all )
y.all <- calcNormFactors(y.all, method="upperquartile", p=pmin.all)
nf.all <- y.all$samples$norm.factors
logy.all <- voom(y.all)$E

median <- apply(logy.all, 1, median)
rle <- apply(logy.all, 2, function(x) abs(x - median))
splitsample <- split(rle, f=pbc.metaD$sample_cell)
sumsample <- sapply(splitsample ,sum)



pbc.metaD$ind_covT1 <- factor(paste(pbc.metaD$ind_cov, pbc.metaD$cg_cov,sep="."))

M <- replicate.matrix(pbc.metaD$ind_covT1)
rownames(M) <- pbc.metaD$ind_covT1

group <- factor(paste(pbc.metaD$fk.tr, pbc.metaD$cg_cov,sep="."))

ruv2T1 <- RUV2mod(Y = t(logy.all), X = group, ctl = hk.ind, k=k, include.intercept = F,Z = pbc.metaD$pop_cov)
ruv3T1 <- RUVIIIW(Y = t(logy.all), M = M, ctl = hk.ind,  k = k, return.info = T)
ruv4T1 <- RUV4mod(Y = t(logy.all), X = group, ctl = hk.ind, k=k, include.intercept = F, Z = pbc.metaD$pop_cov)

# Type3 -------------------------------------------------------------------


bc <- pbc.b
bc <-  bc[rowSums(bc >= 10) >= 5,]
gene.D.bc <- gene.D[ rownames(bc),]

hk.ind <- rownames(bc) %in% hkGagnon[-high_varg]
y.bc <- DGEList(counts=bc)
pmin.bc <- find_p(bc)
y.bc <- calcNormFactors(y.bc, method="upperquartile", p=pmin.bc)
nf.bc <- y.bc$samples$norm.factors
logy.bc <- voom(y.bc)$E

M.bc <- replicate.matrix(bc.metaD$ind_cov)
rownames(M.bc) <- bc.metaD$ind_cov

ruv2T3 <- RUV2mod(Y = t(logy.bc), X = bc.metaD$fk.tr, ctl = hk.ind, k=k, Z = bc.metaD$pop_cov)
ruv3T3 <- RUVIIIW(Y = t(logy.bc), M = M.bc, ctl = hk.ind,  k = k, return.info = T)
ruv4T3 <- RUV4mod(Y = t(logy.bc), X = bc.metaD$fk.tr, ctl = hk.ind, k=k, Z = bc.metaD$pop_cov)


median <- apply(ruv2T3$newY, 2, median)
rle <- apply(ruv2T3$newY, 1, function(x) x - median)

dataplot <- as.data.frame(rle) %>% 
  pivot_longer(everything(),names_to = 'ind_cov_batch_cov') %>% 
  left_join(bc.metaD, by='ind_cov_batch_cov')

dataplot$ind_cov_batch_cov <- factor(dataplot$ind_cov_batch_cov, levels=unique(orddata.bc$ind_cov_batch_cov))

p <- ggplot(dataplot,aes(x= ind_cov_batch_cov, y=value,colour = Processing_Cohort )) +
  geom_boxplot(outlier.shape = NA) + ylim(-2,2) + geom_hline(yintercept = 0,linetype="dashed") +
  labs( colour = "Processing cohort", y= "Relative log expression") + theme_minimal()+ scale_color_brewer(palette='Set1')+
  theme( legend.position = 'bottom',#strip.text =element_blank(),
         axis.title.x = element_blank(),
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank()) 

p 

ggsave(paste0(path2,ct,'RLEruv2T3bulk.png'),height = 5, width = 7)


ct <- 'T4'

pbc <- PBC[rowSums(PBC >= 10) >= 5,pbc.metaD$cg_cov==ct]
y <- DGEList(counts=pbc)
pmin <- find_p(pbc)
y <- calcNormFactors(y, method="upperquartile", p=pmin)
logy <- voom(y)$E
hk.ind <- rownames(pbc) %in% hkGagnon[-high_varg]


# UQ ----------------------------------------------------------------------


top_high_varg.ct <- FindVariableFeatures(logy)
top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]
pca <- calculatePCA(logy,ncomponents=10, subset_row=top_high_varg.ct)
pco <- as.integer(pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))

p <- fviz_silhouette(sil2)+ scale_fill_brewer(palette = 'Set1') + 
  scale_color_brewer(palette = 'Set1') + theme_minimal() + geom_hline(yintercept=0.41, color='black', linetype='dashed')+
  theme(axis.text.x=element_blank(),legend.position = "bottom") +
  labs(fill="Processing Cohort",color="Processing Cohort")

p

ggsave(paste0(path2,ct,'_raw_silhouette.png'),height = 5, width = 7)


avg_sw2 <- round(mean(sil2[,3]),2)
avg_swp <- round(mean(silp[,3]),2)
labelgg <- paste0(c('2 PCs: ','10 PCs: '), c(avg_sw2,avg_swp)) 
pos_x <- min(pca[,1])
pos_y <- min(pca[,2])

p <- data.frame(pca) %>%
  tibble::rownames_to_column("sample_cell") %>% # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(pbc.metaD, by = "sample_cell") %>% mutate(control= ind_cov=='ICC_control') %>%
  ggplot(aes(x = PC1,y = PC2,color= Processing_Cohort))+
  geom_point(aes(shape=pop_cov),size=3) + 
  geom_text(aes(label = ifelse(control, 'control', "")),color='black',vjust = 0, nudge_y = 0.5) +
  labs(color = "Processing cohort", shape = "Ethnicity") + 
  theme_minimal()+ theme(legend.position = "bottom") + scale_color_brewer(palette='Set1') +
  annotate("text", x = c(pos_x+3,pos_x+3), y = c(pos_y+3, pos_y), label = labelgg,   size=3,  fontface="bold")

p

ggsave(paste0(path2,ct,'_PCAraw.png'),height = 5, width = 7)


names <- orddata [orddata$cg_cov==ct,'sample_cell']

p <- ggRLE(logy, covs=pbc.metaD, names=names)+ scale_color_brewer(palette='Set1') +
  geom_hline(yintercept=0, linetype='dashed') + labs(y="Relative log expression", color='Processing cohort')

p

ggsave(paste0(path2,ct,'_RLEraw.png'),height = 5, width = 7)


design <- model.matrix(~ pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct] ) 
colnames(design) <- c('Intercept','fk.tr')
v <- voom(y, design, plot=F) 
vfit <- lmFit(v, design)
efit <- eBayes(vfit)
pvalsUQ <- topTable(efit, coef=2, number=dim(y)[1], adjust="BH")
# vp <- ggplot(pvalsUQ, aes(x=logFC,y=-log10(P.Value)))+ geom_point()+ scale_color_brewer(palette='Set1') + ggtitle("Volcano plot") + theme_minimal()

p <- ggplot(pvalsUQ, aes(x=P.Value)) + geom_histogram(bins=15) + scale_color_brewer(palette='Set1') + ggtitle("P-values histogram") + theme_minimal()
p

ggsave(paste0(path2,ct,'_histraw.png'),height = 5, width = 7)

# UQ+Batch ----------------------------------------------------------------

design <- model.matrix(~ pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct]+ pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct] + pbc.metaD$pop_cov[pbc.metaD$cg_cov==ct] ) 
colnames(design) <- c('(Intercept)','tr',paste0('Proc', 2:4),'pop')
v <- voom(y, design, plot=F) 
vfit <- lmFit(v, design)
efit <- eBayes(vfit)
pvalsUQB <- topTable(efit, coef=2, number=dim(y)[1], adjust="BH")

alpha <- vfit$coefficients[,3:5] # same coefficients as efit and returned by toptable
newY <- t(v$E) - design[,c(3:5)]%*%t(alpha)

top_high_varg.ct <- FindVariableFeatures(t(newY))
top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]
pca <- calculatePCA(t(newY),ncomponents=10, subset_row=top_high_varg.ct)
pco <- as.integer(pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))


avg_sw2 <- round(mean(sil2[,3]),2)
avg_swp <- round(mean(silp[,3]),2)
labelgg <- paste0(c('2 PCs: ','10 PCs: '), c(avg_sw2,avg_swp)) 
pos_x <- min(pca[,1])
pos_y <- min(pca[,2])

p <- data.frame(pca) %>%
  tibble::rownames_to_column("sample_cell") %>% # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(pbc.metaD, by = "sample_cell") %>% mutate(control= ind_cov=='ICC_control') %>%
  ggplot(aes(x = PC1,y = PC2,color= Processing_Cohort))+
  geom_point(aes(shape=pop_cov),size=3) + 
  geom_text(aes(label = ifelse(control, 'control', "")),color='black',vjust = 0, nudge_y = 0.5) +
  labs(color = "Processing cohort", shape = "Ethnicity") + 
  theme_minimal()+ theme(legend.position = "bottom") + scale_color_brewer(palette='Set1') +
  annotate("text", x = c(pos_x+3,pos_x+3), y = c(pos_y+3, pos_y), label = labelgg,   size=3,  fontface="bold")

p

ggsave(paste0(path2,ct,'_PCAbatch.png'),height = 5, width = 7)


p <- ggRLE(t(newY), covs=pbc.metaD, names=names) + scale_color_brewer(palette='Set1') +
  geom_hline(yintercept=0, linetype='dashed') + labs(y="Relative log expression", color='Processing cohort')

p
ggsave(paste0(path2,ct,'_RLEbatch.png'),height = 5, width = 7)

p <- ggplot(pvalsUQB, aes(x=P.Value)) + geom_histogram(bins=15) + scale_color_brewer(palette='Set1') + ggtitle("P-values histogram") + theme_minimal()
p

ggsave(paste0(path2,ct,'_histbatch.png'),height = 5, width = 7)


# T2 normalisation --------------------------------------------------------------------

## ruv2 --------------------------------------------------------------------


ruv2T2 <- RUV2mod(Y = t(logy), X = pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct], ctl = hk.ind, k=k, Z=pbc.metaD$pop_cov[pbc.metaD$cg_cov==ct])

top_high_varg.ct <- FindVariableFeatures(t(ruv2T2$newY))
top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]
pca <- calculatePCA(t(ruv2T2$newY),ncomponents=10, subset_row= top_high_varg.ct)
pco <- as.integer(pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))

avg_sw2 <- round(mean(sil2[,3]),2)
avg_swp <- round(mean(silp[,3]),2)
labelgg <- paste0(c('2 PCs: ','10 PCs: '), c(avg_sw2,avg_swp)) 
pos_x <- min(pca[,1])
pos_y <- min(pca[,2])

p <- data.frame(pca) %>%
  tibble::rownames_to_column("sample_cell") %>% # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(pbc.metaD, by = "sample_cell") %>% mutate(control= ind_cov=='ICC_control') %>%
  ggplot(aes(x = PC1,y = PC2,color= Processing_Cohort))+
  geom_point(aes(shape=pop_cov),size=3) + 
  geom_text(aes(label = ifelse(control, 'control', "")),color='black',vjust = 0, nudge_y = 0.5) +
  labs(color = "Processing cohort", shape = "Ethnicity") + 
  theme_minimal()+ theme(legend.position = "bottom") + scale_color_brewer(palette='Set1') +
  annotate("text", x = c(pos_x+3,pos_x+3), y = c(pos_y+3, pos_y), label = labelgg,   size=3,  fontface="bold")

p

ggsave(paste0(path2,ct,'_PCAruv2T2.png'),height = 5, width = 7)

p <- ggRLE(t(ruv2T2$newY), covs=pbc.metaD, names=names)+ scale_color_brewer(palette='Set1') +
  geom_hline(yintercept=0, linetype='dashed') + labs(y="Relative log expression", color='Processing cohort')

p 

ggsave(paste0(path2,ct,'_RLEruv2T2.png'),height = 5, width = 7)

design <- model.matrix(~ pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct] + ruv2T2$W + pbc.metaD$pop_cov[pbc.metaD$cg_cov==ct]) 
colnames(design) <- c('Intercept','fk.tr', paste0('W',1:k),'pop')
v <- voom(y, design, plot=F) 
vfit <- lmFit(v, design)
efit <- eBayes(vfit)
pvalsruv2t2 <- topTable(efit, coef=2, number=dim(y)[1], adjust="BH")
# vp <- ggplot(pvalsUQ, aes(x=logFC,y=-log10(P.Value)))+ geom_point()+ scale_color_brewer(palette='Set1') + ggtitle("Volcano plot") + theme_minimal()

p <- ggplot(pvalsruv2t2, aes(x=P.Value)) + geom_histogram(bins=15) + scale_color_brewer(palette='Set1') + ggtitle("P-values histogram") + theme_minimal()
p

ggsave(paste0(path2,ct,'_histruv2T2.png'),height = 5, width = 7)


## ruv3 --------------------------------------------------------------------

Mct <- replicate.matrix(pbc.metaD$ind_cov[pbc.metaD$cg_cov==ct])
rownames(Mct) <- pbc.metaD$ind_cov[pbc.metaD$cg_cov==ct]

ruv3T2 <- RUVIIIW(Y = t(logy), M = Mct, ctl=hk.ind,  k = k,
                  return.info = T)

top_high_varg.ct <- FindVariableFeatures(t(ruv3T2$newY))
top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]

pca <- calculatePCA(t(ruv3T2$newY),ncomponents=10, subset_row=top_high_varg.ct)
pco <- as.integer(pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))

p <- fviz_silhouette(sil2)+ scale_fill_brewer(palette = 'Set1') + 
  scale_color_brewer(palette = 'Set1') + theme_minimal() + geom_hline(yintercept=0.17, color='black', linetype='dashed')+
  theme(axis.text.x=element_blank(),legend.position = "bottom") +
  labs(fill="Processing Cohort",color="Processing Cohort")

p

ggsave(paste0(path2,ct,'_ruviii_silhouette.png'),height = 5, width = 7)



avg_sw2 <- round(mean(sil2[,3]),2)
avg_swp <- round(mean(silp[,3]),2)
labelgg <- paste0(c('2 PCs: ','10 PCs: '), c(avg_sw2,avg_swp)) 
pos_x <- min(pca[,1])
pos_y <- min(pca[,2])

p <- data.frame(pca) %>%
  tibble::rownames_to_column("sample_cell") %>% # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(pbc.metaD, by = "sample_cell") %>% mutate(control= ind_cov=='ICC_control') %>%
  ggplot(aes(x = PC1,y = PC2,color= Processing_Cohort))+
  geom_point(aes(shape=pop_cov),size=3) + 
  geom_text(aes(label = ifelse(control, 'control', "")),color='black',vjust = 0, nudge_y = 0.5) +
  labs(color = "Processing cohort", shape = "Ethnicity") + 
  theme_minimal()+ theme(legend.position = "bottom") + scale_color_brewer(palette='Set1') +
  annotate("text", x = c(pos_x+3,pos_x+3), y = c(pos_y+3, pos_y), label = labelgg,   size=3,  fontface="bold")
p

ggsave(paste0(path2,ct,'_PCAruv3T2.png'),height = 5, width = 7)


p <- ggRLE(t(ruv3T2$newY), covs=pbc.metaD, names=names) + scale_color_brewer(palette='Set1') +
  geom_hline(yintercept=0, linetype='dashed') + labs(y="Relative log expression", color='Processing cohort')


p

ggsave(paste0(path2,ct,'_RLEruv3T2.png'),height = 5, width = 7) 


design <- model.matrix(~ pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct] + ruv3T2$W + pbc.metaD$pop_cov[pbc.metaD$cg_cov==ct]) 
colnames(design) <- c('Intercept','fk.tr', paste0('W',1:k),'pop')
v <- voom(y, design, plot=F) 
vfit <- lmFit(v, design)
efit <- eBayes(vfit)
pvalsruv3t2 <- topTable(efit, coef=2, number=dim(y)[1], adjust="BH")
# vp <- ggplot(pvalsUQ, aes(x=logFC,y=-log10(P.Value)))+ geom_point()+ scale_color_brewer(palette='Set1') + ggtitle("Volcano plot") + theme_minimal()

p <- ggplot(pvalsruv3t2, aes(x=P.Value)) + geom_histogram(bins=15) + scale_color_brewer(palette='Set1') + ggtitle("P-values histogram") + theme_minimal()
p

ggsave(paste0(path2,ct,'_histruv3T2.png'),height = 5, width = 7)

# RUV3 + Batch --------------------------------------------------------------------


design <- model.matrix(~pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct]+ ruv3T2$W + 
                         pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct] + pbc.metaD$pop_cov[pbc.metaD$cg_cov==ct] ) 
colnames(design) <- c('(Intercept)','fk.tr',paste0('W', 1:k), paste0('PCohort', 2:4), 'pop')


fit <- lm( ruv3T2$W~pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])






yhat <- predict(fit, pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct], type = "raw")
r2_mv2 <- calc_rsquared(y = ruv3T2$W, yhat = yhat)





apply(ruv3T2$W,2, function (x) kruskal.test(x ~ pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct]))



summary(fit)


v <- voom(y, design, plot=F) 
vfit <- lmFit(v, design)
efit <- eBayes(vfit)
pvalsruv3B <- topTable(efit, coef=2, number=dim(y)[1], adjust="BH")

alpharuv3B <- vfit$coefficients[,3:10] # same coefficients as efit and returned by toptable
newY <- t(v$E) - design[,3:10]%*%t(alpharuv3B)

top_high_varg.ct <- FindVariableFeatures(t(newY))
top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]

pca <- calculatePCA(t(newY),ncomponents=10, subset_row=top_high_varg.ct)

aswPC <- ASW(pca,pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
aswpop <- ASW(pca,pbc.metaD$pop_cov[pbc.metaD$cg_cov==ct])

avg_sw2 <- aswPC[1]
avg_swp <- aswPC[2]
labelgg <- paste0(c('2 PCs: ','10 PCs: '), c(avg_sw2,avg_swp)) 
pos_x <- min(pca[,1])
pos_y <- min(pca[,2])

p <- data.frame(pca) %>%
  tibble::rownames_to_column("sample_cell") %>% # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(pbc.metaD, by = "sample_cell") %>% mutate(control= ind_cov=='ICC_control') %>%
  ggplot(aes(x = PC1,y = PC2,color= Processing_Cohort))+
  geom_point(aes(shape=pop_cov),size=3) + 
  geom_text(aes(label = ifelse(control, 'control', "")),color='black',vjust = 0, nudge_y = 0.5) +
  labs(color = "Processing cohort", shape = "Ethnicity") + 
  theme_minimal()+ theme(legend.position = "bottom") + scale_color_brewer(palette='Set1') +
  annotate("text", x = c(pos_x+3,pos_x+3), y = c(pos_y+3, pos_y), label = labelgg,   size=3,  fontface="bold")

p

ggsave(paste0(path2,ct,'_PCAruv3B.png'),height = 5, width = 7)


p <- ggRLE(t(newY), covs=pbc.metaD, names=names) + scale_color_brewer(palette='Set1') +
  geom_hline(yintercept=0, linetype='dashed') + labs(y="Relative log expression", color='Processing cohort')


p

ggsave(paste0(path2,ct,'_RLEruv3B.png'),height = 5, width = 7) 


p <- ggplot(pvalsruv3B, aes(x=P.Value)) + geom_histogram(bins=15) + scale_color_brewer(palette='Set1') + ggtitle("P-values histogram") + theme_minimal()
p

ggsave(paste0(path2,ct,'_histruv3B.png'),height = 5, width = 7)




## ruv3 pbps ---------------------------------------------------------------


pbps10rep_ds2 <- readRDS("/domino/datasets/local/RUV/pbps10rep_ds2.rds")
fullcount <- pbps10rep_ds2[[1]]
sifull <- pbps10rep_ds2[[2]]

ps.pbc <- fullcount[,sifull$cg_cov==ct]
ps.pbc <- ps.pbc[rowSums(bc >= 10) >= 5,]
ps.pbc <-  ps.pbc[rownames(pbc),]
ps.y <- DGEList(counts=ps.pbc )
pmin <- find_p(ps.pbc )
ps.y <- calcNormFactors(ps.y, method="upperquartile", p=pmin)
ps.nf <- ps.y$samples$norm.factors
ps.logy <- voom(ps.y,plot=F)$E

ps.M <- replicate.matrix(sifull$ind_cov[sifull$cg_cov==ct])
rownames(ps.M) <- sifull$ind_cov[sifull$cg_cov==ct]

ruv3ps <- RUVIIIW(Y = t(ps.logy), M = ps.M, ctl = hk.ind,  k = k,
                  return.info = T)

orig.s <- sifull$sample_cell[sifull$cg_cov==ct & sifull$pbps==0]

top_high_varg.ct <- FindVariableFeatures(t(ruv3ps$newY[orig.s,]))
top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]

pca <- calculatePCA(t(ruv3ps$newY[orig.s,]),ncomponents=10, subset_row=top_high_varg.ct)
aswPC <- ASW(pca,pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
avg_sw2 <- aswPC[1]
avg_swp <- aswPC[2]
labelgg <- paste0(c('2 PCs: ','10 PCs: '), c(avg_sw2,avg_swp)) 
pos_x <- min(pca[,1])
pos_y <- min(pca[,2])

p <- data.frame(pca) %>%
  tibble::rownames_to_column("sample_cell") %>% # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(pbc.metaD, by = "sample_cell") %>% mutate(control= ind_cov=='ICC_control') %>%
  ggplot(aes(x = PC1,y = PC2,color= Processing_Cohort))+
  geom_point(aes(shape=pop_cov),size=3) + 
  geom_text(aes(label = ifelse(control, 'control', "")),color='black',vjust = 0, nudge_y = 0.5) +
  labs(color = "Processing cohort", shape = "Ethnicity") + 
  theme_minimal()+ theme(legend.position = "bottom") + scale_color_brewer(palette='Set1') +
  annotate("text", x = c(pos_x+3,pos_x+3), y = c(pos_y+3, pos_y), label = labelgg,   size=3,  fontface="bold")

p

ggsave(paste0(path2,ct,'_PCAruv3pbps.png'),height = 5, width = 7)


p <- ggRLE(t(ruv3ps$newY[orig.s,]), covs=pbc.metaD, names=names) + scale_color_brewer(palette='Set1') +
  geom_hline(yintercept=0, linetype='dashed') + labs(y="Relative log expression", color='Processing cohort')


p

ggsave(paste0(path2,ct,'_RLEruv3pbps.png'),height = 5, width = 7) 


# T1 normalisation --------------------------------------------------------

## ruv2 --------------------------------------------------------------------

top_high_varg.ct <- FindVariableFeatures(t(ruv2T1$newY[pbc.metaD$cg_cov==ct,]))
top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]
pca <- calculatePCA(t(ruv2T1$newY[pbc.metaD$cg_cov==ct,]),ncomponents=10, subset_row=top_high_varg.ct)
pco <- as.integer(pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))

avg_sw2 <- round(mean(sil2[,3]),2)
avg_swp <- round(mean(silp[,3]),2)
labelgg <- paste0(c('2 PCs: ','10 PCs: '), c(avg_sw2,avg_swp)) 
pos_x <- min(pca[,1])
pos_y <- min(pca[,2])


p <- data.frame(pca)%>%
  tibble::rownames_to_column("sample_cell") %>% # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(pbc.metaD, by = "sample_cell") %>% mutate(control= ind_cov=='ICC_control') %>%
  ggplot(aes(x = PC1,y = PC2,color= Processing_Cohort))+
  geom_point(aes(shape=pop_cov),size=3) + 
  geom_text(aes(label = ifelse(control, 'control', "")),color='black',vjust = 0, nudge_y = 0.5) +
  labs(color = "Processing cohort", shape = "Ethnicity") + 
  theme_minimal()+ theme(legend.position = "bottom") + scale_color_brewer(palette='Set1') +
  annotate("text", x = c(pos_x+3,pos_x+3), y = c(pos_y+3, pos_y), label = labelgg,   size=3,  fontface="bold")

p  

ggsave(paste0(path2,ct,'_PCAruv2T1.png'),height = 5, width = 7)


p <- ggRLE(t(ruv2T1$newY[pbc.metaD$cg_cov==ct,]), covs=pbc.metaD, names=names) + scale_color_brewer(palette='Set1') +
  geom_hline(yintercept=0, linetype='dashed') + labs(y="Relative log expression", color='Processing cohort')


p

ggsave(paste0(path2,ct,'_RLEruv2T1.png'),height = 5, width = 7)


# T3 normalisation -------------------------------------------------------------

indsamp.ct <- bc.metaD$ind_cov_batch_cov %in% pbc.metaD$ind_cov_batch_cov[pbc.metaD$cg_cov==ct]
indg.ct <- rownames(bc) %in% rownames(logy)
indg.ct2 <- rownames(logy) %in% rownames(bc) 

newY <- t(logy[indg.ct2,]) - ruv2T3$W[indsamp.ct,] %*% ruv2T3$fullalpha[,indg.ct]

top_high_varg.ct <- FindVariableFeatures(t(newY))
top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]
pca <- calculatePCA(t(newY),ncomponents=10, subset_row=top_high_varg.ct)
pco <- as.integer(pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))

avg_sw2 <- round(mean(sil2[,3]),2)
avg_swp <- round(mean(silp[,3]),2)
labelgg <- paste0(c('2 PCs: ','10 PCs: '), c(avg_sw2,avg_swp)) 
pos_x <- min(pca[,1])
pos_y <- min(pca[,2])

p <- data.frame(pca) %>%
  tibble::rownames_to_column("sample_cell") %>% # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(pbc.metaD, by = "sample_cell") %>% mutate(control= ind_cov=='ICC_control') %>%
  ggplot(aes(x = PC1,y = PC2,color= Processing_Cohort))+
  geom_point(aes(shape=pop_cov),size=3) + 
  geom_text(aes(label = ifelse(control, 'control', "")),color='black',vjust = 0, nudge_y = 0.5) +
  labs(color = "Processing cohort", shape = "Ethnicity") + 
  theme_minimal()+ theme(legend.position = "bottom") + scale_color_brewer(palette='Set1') +
  annotate("text", x = c(pos_x+3,pos_x+3), y = c(pos_y+3, pos_y), label = labelgg,   size=3,  fontface="bold")

p

ggsave(paste0(path2,ct,'_PCAruv2T3.png'),height = 5, width = 7)

p <- ggRLE(t(newY), covs=pbc.metaD, names=names)  + scale_color_brewer(palette='Set1') +
  geom_hline(yintercept=0, linetype='dashed') + labs(y="Relative log expression", color='Processing cohort')


p
ggsave(paste0(path2,ct,'_RLEruv2T3.png'),height = 5, width = 7)


## ruv3 pbps samples and pseudosamples plotted ---------------------------------------------------------------


pbps10rep_ds2 <- readRDS("/domino/datasets/local/RUV/pbps10rep_ds2.rds")
fullcount <- pbps10rep_ds2[[1]]
sifull <- pbps10rep_ds2[[2]]


ps.pbc <- fullcount[,sifull$cg_cov==ct]
ps.pbc <- ps.pbc[rowSums(bc >= 10) >= 5,]
ps.pbc <-  ps.pbc[rownames(pbc),]
ps.y <- DGEList(counts=ps.pbc )
pmin <- find_p(ps.pbc )
ps.y <- calcNormFactors(ps.y, method="upperquartile", p=pmin)
ps.nf <- ps.y$samples$norm.factors
ps.logy <- voom(ps.y,plot=F)$E

ps.M <- replicate.matrix(sifull$ind_cov[sifull$cg_cov==ct])
rownames(ps.M) <- sifull$ind_cov[sifull$cg_cov==ct]

ruv3ps <- RUVIIIW(Y = t(ps.logy), M = ps.M, ctl = rep(T,nrow(ps.logy)),  k = k,
                  return.info = T)

orig.s <- sifull$sample_cell[sifull$cg_cov==ct & sifull$pbps==0]

top_high_varg.ct <- FindVariableFeatures(t(ruv3ps$newY[orig.s,]))
top_high_varg.ct <- rownames(arrange(top_high_varg.ct,-vst.variance.standardized))[1:500]

pca <- calculatePCA(t(ruv3ps$newY[orig.s,]),ncomponents=10, subset_row=top_high_varg.ct)
aswPC <- ASW(pca,pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
avg_sw2 <- aswPC[1]
avg_swp <- aswPC[2]
labelgg <- paste0(c('2 PCs: ','10 PCs: '), c(avg_sw2,avg_swp)) 
pos_x <- min(pca[,1])
pos_y <- min(pca[,2])

p <- data.frame(pca) %>%
  tibble::rownames_to_column("sample_cell") %>% # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(pbc.metaD, by = "sample_cell") %>% mutate(control= ind_cov=='ICC_control') %>%
  ggplot(aes(x = PC1,y = PC2,color= Processing_Cohort))+
  geom_point(aes(shape=pop_cov),size=3) + 
  geom_text(aes(label = ifelse(control, 'control', "")),color='black',vjust = 0, nudge_y = 0.5) +
  labs(color = "Processing cohort", shape = "Ethnicity") + 
  theme_minimal()+ theme(legend.position = "bottom") + scale_color_brewer(palette='Set1') +
  annotate("text", x = c(pos_x+3,pos_x+3), y = c(pos_y+3, pos_y), label = labelgg,   size=3,  fontface="bold")

p

ggsave(paste0(path2,ct,'_PCAruv3pbpsall.png'),height = 5, width = 7)


p <- ggRLE(t(ruv3ps$newY[orig.s,]), covs=pbc.metaD, names=names) + scale_color_brewer(palette='Set1') +
  geom_hline(yintercept=0, linetype='dashed') + labs(y="Relative log expression", color='Processing cohort')


p

ggsave(paste0(path2,ct,'_RLEruv3pbpsall.png'),height = 5, width = 7) 

