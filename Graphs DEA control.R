load('/domino/datasets/local/RUV//simsD1.RData')
#library(egg)
library(tidyverse)
library(rlist)
library(ggpubr)
library(iCOBRA)
source("/mnt/auxf.R")

path <- ('/domino/datasets/local/RUV')



# including pop cov E matrix balanced -------------------------------------------------------

pvals_ruvdg_tr_imbalance_pop_E <- readRDS("/domino/datasets/local/RUV/pvals_ruvdg_tr_imbalance_pop_E.rds")

pvalitsruv2 <- pvals_ruvdg_diffg_pop_E[[1]]
pvalitsruv3 <- pvals_ruvdg_diffg_pop_E[[2]]
pvalitsruv4 <- pvals_ruvdg_diffg_pop_E[[3]]
truthits <- pvals_ruvdg_diffg_pop_E[[4]]

pvalits2 <- list_transpose(pvalitsruv2)$T4

pvalits2 <- bind_rows(pvalits2) %>% remove_rownames()%>%column_to_rownames('ind') 
truthits2 <- list_transpose(truthits)$T4
truthits2 <- bind_rows(truthits2) %>% remove_rownames()%>%column_to_rownames('ind') 

COBRA <- COBRAData(pval=pvalits2,truth = truthits2)
COBRA <- calculate_adjp(COBRA)
COBRA <- calculate_performance(cobradata = COBRA, binary_truth = "trueDE_T4")
cobratoplot <- prepare_data_for_plot(COBRA, facetted = TRUE)

RUV2PLOT <- plot_fdrtprcurve(cobratoplot,plottype='points',pointsize=0.5) + 
  geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
  xlim(0,0.5) +
  ylim(0.82,0.95) +
  theme_minimal()+ ggtitle("RUV2") + 
  scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
  theme(strip.background = element_blank(),strip.text =element_blank())


pvalits2 <- list_transpose(pvalitsruv3)$T4
pvalits2 <- bind_rows(pvalits2) %>% remove_rownames()%>%column_to_rownames('ind') 

COBRA <- COBRAData(pval=pvalits2,truth = truthits2)
COBRA <- calculate_adjp(COBRA)
COBRA <- calculate_performance(cobradata = COBRA, binary_truth = "trueDE_T4")
cobratoplot <- prepare_data_for_plot(COBRA, facetted = TRUE)
RUV3PLOT <- plot_fdrtprcurve(cobratoplot,plottype='points',pointsize=0.5) + 
  geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
  xlim(0,0.5) +
  ylim(0.82,0.95) + 
  theme_minimal()+ ggtitle("RUVIII") +
  scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
  theme(strip.background = element_blank(),strip.text =element_blank())

pvalits2 <- list_transpose(pvalitsruv4)$T4
pvalits2 <- bind_rows(pvalits2) %>% remove_rownames()%>%column_to_rownames('ind') 


COBRA <- COBRAData(pval=pvalits2,truth = truthits2)
COBRA <- calculate_adjp(COBRA)
COBRA <- calculate_performance(cobradata = COBRA, binary_truth = "trueDE_T4")
cobratoplot <- prepare_data_for_plot(COBRA, facetted = TRUE)
RUV4PLOT <- plot_fdrtprcurve(cobratoplot,plottype='points',pointsize=0.5) + 
  geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
  xlim(0,0.5) +
  ylim(0.8,0.95) + 
  theme_minimal()+ ggtitle("RUV4") +
  scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
  theme(strip.background = element_blank(),strip.text =element_blank())



RUV2PLOT <- RUV2PLOT + xlim(0,0.5) +
  ylim(0.8,0.92) 

RUV3PLOT <- RUV3PLOT + xlim(0,0.5) +
  ylim(0.8,0.92)

RUV4PLOT <- RUV4PLOT + xlim(0,0.5) +
  ylim(0.8,0.92)

ggpubr::ggarrange(plotlist=list(RUV2PLOT,RUV3PLOT,RUV4PLOT),nrow=1,ncol=3,
                  common.legend = T, legend="bottom")
ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/FDRDGT4popE.png'),height = 4, width = 12)



pvalitsruv2 <- pvals_ruvdg_diffg_pop_E[[1]]
pvalitsruv3 <- pvals_ruvdg_diffg_pop_E[[2]]
pvalitsruv4 <- pvals_ruvdg_diffg_pop_E[[3]]
truthits <- pvals_ruvdg_diffg_pop_E[[4]]


pvalits2 <- list_transpose(pvalitsruv2)
pvalits2 <- lapply(pvalits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )

truthits2 <- list_transpose(truthits)
truthits2 <- lapply(truthits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )


COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
                        binary_truth = paste0("trueDE_",celltypes))
cobratoplot <- lapply(cobraperflist, prepare_data_for_plot, colorscheme = "Dark2",
                      facetted = TRUE)

fdrplots <- lapply(cobratoplot,function (x) plot_fdrtprcurve(x,plottype='points',pointsize=0.5) + geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                     xlim(0,0.6) +
                     ylim(0.75,0.92) + 
                     theme_minimal()+
                     scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
                     theme(strip.background = element_blank(),strip.text =element_blank())
)

ggpubr::ggarrange(plotlist=fdrplots,nrow=2,ncol=4,labels=names(fdrplots),
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/allcellsRUV2FDRG_popE.png'),height = 4, width = 12)

pvalits2 <- list_transpose(pvalitsruv3)
pvalits2 <- lapply(pvalits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )

COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
                        binary_truth = paste0("trueDE_",celltypes))
cobratoplot <- lapply(cobraperflist, prepare_data_for_plot, colorscheme = "Dark2",
                      facetted = TRUE)

fdrplots <- lapply(cobratoplot,function (x) plot_fdrtprcurve(x,plottype='points',pointsize=0.5) + geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                     xlim(0,0.6) +
                     ylim(0.75,0.92) + 
                     theme_minimal()+
                     scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
                     theme(strip.background = element_blank(),strip.text =element_blank())
)

ggpubr::ggarrange(plotlist=fdrplots,nrow=2,ncol=4,labels=names(fdrplots),
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/allcellsRUV3FDRDG_popE.png'),height = 4, width = 12)

pvalits2 <- list_transpose(pvalitsruv4)
pvalits2 <- lapply(pvalits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )

truthits2 <- list_transpose(truthits)
truthits2 <- lapply(truthits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )


COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
                        binary_truth = paste0("trueDE_",celltypes))
cobratoplot <- lapply(cobraperflist, prepare_data_for_plot, colorscheme = "Dark2",
                      facetted = TRUE)

fdrplots <- lapply(cobratoplot,function (x) plot_fdrtprcurve(x,plottype='points',pointsize=0.5) + geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                     xlim(0,0.6) +
                     ylim(0.75,0.92) + 
                     theme_minimal()+
                     scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
                     theme(strip.background = element_blank(),strip.text =element_blank())
)

ggpubr::ggarrange(plotlist=fdrplots,nrow=2,ncol=4,labels=names(fdrplots),
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/allcellsRUV4FDRDG_popE.png'),height = 4, width = 12)





# including pop cov E matrix imbalanced -------------------------------------------------------

pvals_ruvdg_tr_imbalance_pop_E <- readRDS("/domino/datasets/local/RUV/pvals_ruvdg_tr_imbalance_pop_E.rds")

pvalitsruv2 <- pvals_ruvdg_tr_imbalance_pop_E[[1]]
pvalitsruv3 <- pvals_ruvdg_tr_imbalance_pop_E[[2]]
pvalitsruv4 <- pvals_ruvdg_tr_imbalance_pop_E[[3]]
truthits <- pvals_ruvdg_tr_imbalance_pop_E[[4]]

pvalits2 <- list_transpose(pvalitsruv2)$T4

pvalits2 <- bind_rows(pvalits2) %>% remove_rownames()%>%column_to_rownames('ind') 
truthits2 <- list_transpose(truthits)$T4
truthits2 <- bind_rows(truthits2) %>% remove_rownames()%>%column_to_rownames('ind') 

COBRA <- COBRAData(pval=pvalits2,truth = truthits2)
COBRA <- calculate_adjp(COBRA)
COBRA <- calculate_performance(cobradata = COBRA, binary_truth = "trueDE_T4")
cobratoplot <- prepare_data_for_plot(COBRA, facetted = TRUE)

RUV2PLOT <- plot_fdrtprcurve(cobratoplot,plottype='points',pointsize=0.5) + 
  geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
  xlim(0,0.5) +
  ylim(0.82,0.95) +
  theme_minimal()+ ggtitle("RUV2") + 
  scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
  theme(strip.background = element_blank(),strip.text =element_blank())


pvalits2 <- list_transpose(pvalitsruv3)$T4
pvalits2 <- bind_rows(pvalits2) %>% remove_rownames()%>%column_to_rownames('ind') 

COBRA <- COBRAData(pval=pvalits2,truth = truthits2)
COBRA <- calculate_adjp(COBRA)
COBRA <- calculate_performance(cobradata = COBRA, binary_truth = "trueDE_T4")
cobratoplot <- prepare_data_for_plot(COBRA, facetted = TRUE)
RUV3PLOT <- plot_fdrtprcurve(cobratoplot,plottype='points',pointsize=0.5) + 
  geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
  xlim(0,0.5) +
  ylim(0.82,0.95) + 
  theme_minimal()+ ggtitle("RUVIII") +
  scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
  theme(strip.background = element_blank(),strip.text =element_blank())

pvalits2 <- list_transpose(pvalitsruv4)$T4
pvalits2 <- bind_rows(pvalits2) %>% remove_rownames()%>%column_to_rownames('ind') 


COBRA <- COBRAData(pval=pvalits2,truth = truthits2)
COBRA <- calculate_adjp(COBRA)
COBRA <- calculate_performance(cobradata = COBRA, binary_truth = "trueDE_T4")
cobratoplot <- prepare_data_for_plot(COBRA, facetted = TRUE)
RUV4PLOT <- plot_fdrtprcurve(cobratoplot,plottype='points',pointsize=0.5) + 
  geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
  xlim(0,0.5) +
  ylim(0.8,0.95) + 
  theme_minimal()+ ggtitle("RUV4") +
  scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
  theme(strip.background = element_blank(),strip.text =element_blank())



RUV2PLOT <- RUV2PLOT + xlim(0,0.6) +
  ylim(0.75,0.92) 

RUV3PLOT <- RUV3PLOT + xlim(0,0.6) +
  ylim(0.75,0.92)

RUV4PLOT <- RUV4PLOT + xlim(0,0.6) +
  ylim(0.75,0.92)

ggpubr::ggarrange(plotlist=list(RUV2PLOT,RUV3PLOT,RUV4PLOT),nrow=1,ncol=3,
                  common.legend = T, legend="bottom")
ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/FDRDimbT4popE.png'),height = 4, width = 12)



pvalitsruv2 <- pvals_ruvdg_tr_imbalance_pop_E[[1]]
pvalitsruv3 <- pvals_ruvdg_tr_imbalance_pop_E[[2]]
pvalitsruv4 <- pvals_ruvdg_tr_imbalance_pop_E[[3]]
truthits <- pvals_ruvdg_tr_imbalance_pop_E[[4]]


pvalits2 <- list_transpose(pvalitsruv2)
pvalits2 <- lapply(pvalits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )

truthits2 <- list_transpose(truthits)
truthits2 <- lapply(truthits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )


COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
                        binary_truth = paste0("trueDE_",celltypes))
cobratoplot <- lapply(cobraperflist, prepare_data_for_plot, colorscheme = "Dark2",
                      facetted = TRUE)

fdrplots <- lapply(cobratoplot,function (x) plot_fdrtprcurve(x,plottype='points',pointsize=0.5) + geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                     xlim(0,0.6) +
                     ylim(0.75,0.92) + 
                     theme_minimal()+
                     scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
                     theme(strip.background = element_blank(),strip.text =element_blank())
)

ggpubr::ggarrange(plotlist=fdrplots,nrow=2,ncol=4,labels=names(fdrplots),
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/allcellsRUV2FDRDimb_popE.png'),height = 4, width = 12)

pvalits2 <- list_transpose(pvalitsruv3)
pvalits2 <- lapply(pvalits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )

COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
                        binary_truth = paste0("trueDE_",celltypes))
cobratoplot <- lapply(cobraperflist, prepare_data_for_plot, colorscheme = "Dark2",
                      facetted = TRUE)

fdrplots <- lapply(cobratoplot,function (x) plot_fdrtprcurve(x,plottype='points',pointsize=0.5) + geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                     xlim(0,0.6) +
                     ylim(0.75,0.92) + 
                     theme_minimal()+
                     scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
                     theme(strip.background = element_blank(),strip.text =element_blank())
)

ggpubr::ggarrange(plotlist=fdrplots,nrow=2,ncol=4,labels=names(fdrplots),
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/allcellsRUV3FDRDimb_popE.png'),height = 4, width = 12)

pvalits2 <- list_transpose(pvalitsruv4)
pvalits2 <- lapply(pvalits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )

truthits2 <- list_transpose(truthits)
truthits2 <- lapply(truthits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )


COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
                        binary_truth = paste0("trueDE_",celltypes))
cobratoplot <- lapply(cobraperflist, prepare_data_for_plot, colorscheme = "Dark2",
                      facetted = TRUE)

fdrplots <- lapply(cobratoplot,function (x) plot_fdrtprcurve(x,plottype='points',pointsize=0.5) + geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                     xlim(0,0.6) +
                     ylim(0.75,0.92) + 
                     theme_minimal()+
                     scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
                     theme(strip.background = element_blank(),strip.text =element_blank())
)

ggpubr::ggarrange(plotlist=fdrplots,nrow=2,ncol=4,labels=names(fdrplots),
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/allcellsRUV4FDRDimb_popE.png'),height = 4, width = 12)


# including pop cov E with interaction matrix-------------------------------------------------------

pvals_ruvdg_tr_imbalance_pop_E <- readRDS("/domino/datasets/local/RUV/pvals_ruvdg_tr_imbalance_pop_Eint.rds")

pvalitsruv2 <- pvals_ruvdg_tr_imbalance_pop_E[[1]]
pvalitsruv3 <- pvals_ruvdg_tr_imbalance_pop_E[[2]]
pvalitsruv4 <- pvals_ruvdg_tr_imbalance_pop_E[[3]]
truthits <- pvals_ruvdg_tr_imbalance_pop_E[[4]]

pvalits2 <- list_transpose(pvalitsruv2)$T4

pvalits2 <- bind_rows(pvalits2) %>% remove_rownames()%>%column_to_rownames('ind') 
truthits2 <- list_transpose(truthits)$T4
truthits2 <- bind_rows(truthits2) %>% remove_rownames()%>%column_to_rownames('ind') 

COBRA <- COBRAData(pval=pvalits2,truth = truthits2)
COBRA <- calculate_adjp(COBRA)
COBRA <- calculate_performance(cobradata = COBRA, binary_truth = "trueDE_T4")
cobratoplot <- prepare_data_for_plot(COBRA, facetted = TRUE)

RUV2PLOT <- plot_fdrtprcurve(cobratoplot,plottype='points',pointsize=0.5) + 
  geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
  xlim(0,0.6) +
  ylim(0.7,0.95) +
  theme_minimal()+ ggtitle("RUV2") + 
  scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
  theme(strip.background = element_blank(),strip.text =element_blank())


pvalits2 <- list_transpose(pvalitsruv3)$T4
pvalits2 <- bind_rows(pvalits2) %>% remove_rownames()%>%column_to_rownames('ind') 

COBRA <- COBRAData(pval=pvalits2,truth = truthits2)
COBRA <- calculate_adjp(COBRA)
COBRA <- calculate_performance(cobradata = COBRA, binary_truth = "trueDE_T4")
cobratoplot <- prepare_data_for_plot(COBRA, facetted = TRUE)
RUV3PLOT <- plot_fdrtprcurve(cobratoplot,plottype='points',pointsize=0.5) + 
  geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
  xlim(0,0.6) +
  ylim(0.7,0.95) + 
  theme_minimal()+ ggtitle("RUVIII") +
  scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
  theme(strip.background = element_blank(),strip.text =element_blank())

pvalits2 <- list_transpose(pvalitsruv4)$T4
pvalits2 <- bind_rows(pvalits2) %>% remove_rownames()%>%column_to_rownames('ind') 


COBRA <- COBRAData(pval=pvalits2,truth = truthits2)
COBRA <- calculate_adjp(COBRA)
COBRA <- calculate_performance(cobradata = COBRA, binary_truth = "trueDE_T4")
cobratoplot <- prepare_data_for_plot(COBRA, facetted = TRUE)
RUV4PLOT <- plot_fdrtprcurve(cobratoplot,plottype='points',pointsize=0.5) + 
  geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
  xlim(0,0.6) +
  ylim(0.7,0.95) +
  theme_minimal()+ ggtitle("RUV4") +
  scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
  theme(strip.background = element_blank(),strip.text =element_blank())



RUV2PLOT <- RUV2PLOT + xlim(0,0.6) +
  ylim(0.6,0.92) 

RUV3PLOT <- RUV3PLOT + xlim(0,0.6) +
  ylim(0.6,0.92)

RUV4PLOT <- RUV4PLOT + xlim(0,0.6) +
  ylim(0.6,0.92)

ggpubr::ggarrange(plotlist=list(RUV2PLOT,RUV3PLOT,RUV4PLOT),nrow=1,ncol=3,
                  common.legend = T, legend="bottom")
ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/FDRDimbT4popEint.png'),height = 4, width = 12)



pvalitsruv2 <- pvals_ruvdg_tr_imbalance_pop_E[[1]]
pvalitsruv3 <- pvals_ruvdg_tr_imbalance_pop_E[[2]]
pvalitsruv4 <- pvals_ruvdg_tr_imbalance_pop_E[[3]]
truthits <- pvals_ruvdg_tr_imbalance_pop_E[[4]]


pvalits2 <- list_transpose(pvalitsruv2)
pvalits2 <- lapply(pvalits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )

truthits2 <- list_transpose(truthits)
truthits2 <- lapply(truthits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )


COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
                        binary_truth = paste0("trueDE_",celltypes))
cobratoplot <- lapply(cobraperflist, prepare_data_for_plot, colorscheme = "Dark2",
                      facetted = TRUE)

fdrplots <- lapply(cobratoplot,function (x) plot_fdrtprcurve(x,plottype='points',pointsize=0.5) + geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                     xlim(0,0.6) +
                     ylim(0.6,0.95) + 
                     theme_minimal()+
                     scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
                     theme(strip.background = element_blank(),strip.text =element_blank())
)

ggpubr::ggarrange(plotlist=fdrplots,nrow=2,ncol=4,labels=names(fdrplots),
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/allcellsRUV2FDRDimb_popEint.png'),height = 4, width = 12)

pvalits2 <- list_transpose(pvalitsruv3)
pvalits2 <- lapply(pvalits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )

COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
                        binary_truth = paste0("trueDE_",celltypes))
cobratoplot <- lapply(cobraperflist, prepare_data_for_plot, colorscheme = "Dark2",
                      facetted = TRUE)

fdrplots <- lapply(cobratoplot,function (x) plot_fdrtprcurve(x,plottype='points',pointsize=0.5) + geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                     xlim(0,0.6) +
                     ylim(0.6,0.95) + 
                     theme_minimal()+
                     scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
                     theme(strip.background = element_blank(),strip.text =element_blank())
)

ggpubr::ggarrange(plotlist=fdrplots,nrow=2,ncol=4,labels=names(fdrplots),
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/allcellsRUV3FDRDimb_popEint.png'),height = 4, width = 12)

pvalits2 <- list_transpose(pvalitsruv4)
pvalits2 <- lapply(pvalits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )

truthits2 <- list_transpose(truthits)
truthits2 <- lapply(truthits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )


COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
                        binary_truth = paste0("trueDE_",celltypes))
cobratoplot <- lapply(cobraperflist, prepare_data_for_plot, colorscheme = "Dark2",
                      facetted = TRUE)

fdrplots <- lapply(cobratoplot,function (x) plot_fdrtprcurve(x,plottype='points',pointsize=0.5) + geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                     xlim(0,0.6) +
                     ylim(0.6,0.95) + 
                     theme_minimal()+
                     scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
                     theme(strip.background = element_blank(),strip.text =element_blank())
)

ggpubr::ggarrange(plotlist=fdrplots,nrow=2,ncol=4,labels=names(fdrplots),
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/allcellsRUV4FDRDimb_popEint.png'),height = 4, width = 12)


# including pop cov Ewith interaction matrix balanced -------------------------------------------------------

pvals_ruvdg_diffg_pop_E <- readRDS("/domino/datasets/local/RUV/pvals_ruvdg_diffg_pop_Eint.rds")

pvalitsruv2 <- pvals_ruvdg_diffg_pop_E[[1]]
pvalitsruv3 <- pvals_ruvdg_diffg_pop_E[[2]]
pvalitsruv4 <- pvals_ruvdg_diffg_pop_E[[3]]
truthits <- pvals_ruvdg_diffg_pop_E[[4]]

pvalits2 <- list_transpose(pvalitsruv2)$T4

pvalits2 <- bind_rows(pvalits2) %>% remove_rownames()%>%column_to_rownames('ind') 
truthits2 <- list_transpose(truthits)$T4
truthits2 <- bind_rows(truthits2) %>% remove_rownames()%>%column_to_rownames('ind') 

COBRA <- COBRAData(pval=pvalits2,truth = truthits2)
COBRA <- calculate_adjp(COBRA)
COBRA <- calculate_performance(cobradata = COBRA, binary_truth = "trueDE_T4")
cobratoplot <- prepare_data_for_plot(COBRA, facetted = TRUE)

RUV2PLOT <- plot_fdrtprcurve(cobratoplot,plottype='points',pointsize=0.5) + 
  geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
  xlim(0,0.6) +
  ylim(0.6,0.95) +
  theme_minimal()+ ggtitle("RUV2") + 
  scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
  theme(strip.background = element_blank(),strip.text =element_blank())


pvalits2 <- list_transpose(pvalitsruv3)$T4
pvalits2 <- bind_rows(pvalits2) %>% remove_rownames()%>%column_to_rownames('ind') 

COBRA <- COBRAData(pval=pvalits2,truth = truthits2)
COBRA <- calculate_adjp(COBRA)
COBRA <- calculate_performance(cobradata = COBRA, binary_truth = "trueDE_T4")
cobratoplot <- prepare_data_for_plot(COBRA, facetted = TRUE)
RUV3PLOT <- plot_fdrtprcurve(cobratoplot,plottype='points',pointsize=0.5) + 
  geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
  xlim(0,0.6) +
  ylim(0.6,0.95) + 
  theme_minimal()+ ggtitle("RUVIII") +
  scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
  theme(strip.background = element_blank(),strip.text =element_blank())

pvalits2 <- list_transpose(pvalitsruv4)$T4
pvalits2 <- bind_rows(pvalits2) %>% remove_rownames()%>%column_to_rownames('ind') 


COBRA <- COBRAData(pval=pvalits2,truth = truthits2)
COBRA <- calculate_adjp(COBRA)
COBRA <- calculate_performance(cobradata = COBRA, binary_truth = "trueDE_T4")
cobratoplot <- prepare_data_for_plot(COBRA, facetted = TRUE)
RUV4PLOT <- plot_fdrtprcurve(cobratoplot,plottype='points',pointsize=0.5) + 
  geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
  xlim(0,0.6) +
  ylim(0.6,0.95) + 
  theme_minimal()+ ggtitle("RUV4") +
  scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
  theme(strip.background = element_blank(),strip.text =element_blank())



# RUV2PLOT <- RUV2PLOT + xlim(0,0.5) +
#   ylim(0.8,0.92) 
# 
# RUV3PLOT <- RUV3PLOT + xlim(0,0.5) +
#   ylim(0.8,0.92)
# 
# RUV4PLOT <- RUV4PLOT + xlim(0,0.5) +
#   ylim(0.8,0.92)

ggpubr::ggarrange(plotlist=list(RUV2PLOT,RUV3PLOT,RUV4PLOT),nrow=1,ncol=3,
                  common.legend = T, legend="bottom")
ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/FDRDGT4popEint.png'),height = 4, width = 12)



pvalitsruv2 <- pvals_ruvdg_diffg_pop_E[[1]]
pvalitsruv3 <- pvals_ruvdg_diffg_pop_E[[2]]
pvalitsruv4 <- pvals_ruvdg_diffg_pop_E[[3]]
truthits <- pvals_ruvdg_diffg_pop_E[[4]]


pvalits2 <- list_transpose(pvalitsruv2)
pvalits2 <- lapply(pvalits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )

truthits2 <- list_transpose(truthits)
truthits2 <- lapply(truthits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )


COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
                        binary_truth = paste0("trueDE_",celltypes))
cobratoplot <- lapply(cobraperflist, prepare_data_for_plot, colorscheme = "Dark2",
                      facetted = TRUE)

fdrplots <- lapply(cobratoplot,function (x) plot_fdrtprcurve(x,plottype='points',pointsize=0.5) + geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                     xlim(0,0.6) +
                     ylim(0.6,0.95) + 
                     theme_minimal()+
                     scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
                     theme(strip.background = element_blank(),strip.text =element_blank())
)

ggpubr::ggarrange(plotlist=fdrplots,nrow=2,ncol=4,labels=names(fdrplots),
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/allcellsRUV2FDRG_popEint.png'),height = 8, width = 12)

pvalits2 <- list_transpose(pvalitsruv3)
pvalits2 <- lapply(pvalits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )

COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
                        binary_truth = paste0("trueDE_",celltypes))
cobratoplot <- lapply(cobraperflist, prepare_data_for_plot, colorscheme = "Dark2",
                      facetted = TRUE)

fdrplots <- lapply(cobratoplot,function (x) plot_fdrtprcurve(x,plottype='points',pointsize=0.5) + geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                     xlim(0,0.6) +
                     ylim(0.6,0.95) + 
                     theme_minimal()+
                     scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
                     theme(strip.background = element_blank(),strip.text =element_blank())
)

ggpubr::ggarrange(plotlist=fdrplots,nrow=2,ncol=4,labels=names(fdrplots),
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/allcellsRUV3FDRDG_popEint.png'),height = 8, width = 12)

pvalits2 <- list_transpose(pvalitsruv4)
pvalits2 <- lapply(pvalits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )

truthits2 <- list_transpose(truthits)
truthits2 <- lapply(truthits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )


COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
                        binary_truth = paste0("trueDE_",celltypes))
cobratoplot <- lapply(cobraperflist, prepare_data_for_plot, colorscheme = "Dark2",
                      facetted = TRUE)

fdrplots <- lapply(cobratoplot,function (x) plot_fdrtprcurve(x,plottype='points',pointsize=0.5) + geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                     xlim(0,0.6) +
                     ylim(0.6,0.95) + 
                     theme_minimal()+
                     scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
                     theme(strip.background = element_blank(),strip.text =element_blank())
)

ggpubr::ggarrange(plotlist=fdrplots,nrow=2,ncol=4,labels=names(fdrplots),
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/allcellsRUV4FDRDG_popEint.png'),height = 8, width = 12)


# Misspecified negative control genes ------------------------------------------------------------


missp <- readRDS("/domino/datasets/local/RUV/pvals_misspecifiednc_tr_imbalance_popE.rds")

pvalsmiss <- missp[[1]] 
truthmiss <- missp[[2]] 

celltypes <- c('B', 'NK', 'T4', 'T8', 'cDC', 'cM', 'ncM', 'pDC')

pvalits2 <- list_transpose(pvalsmiss)
pvalits2 <- pvalits2[celltypes]
pvalits2  <- lapply(pvalits2 , function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )

truthits2 <- list_transpose(truthmiss )
truthits2 <- truthits2[celltypes]
truthits2 <- lapply(truthits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )


COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
                        binary_truth = paste0("trueDE_",celltypes))
cobratoplot <- lapply(cobraperflist, prepare_data_for_plot,
                      facetted = TRUE)

plotscobra <- lapply(cobratoplot,plot_fdrtprcurve,plottype='points',pointsize=0.5)


plotscobra1 <- mapply( function(x,y) x + 
                         geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                         xlim(0,0.6) +
                         ylim(0.3,0.93) +
                         theme_minimal()+ ggtitle(y) +
                         scale_color_brewer(palette='Set1',labels = c('RUV2', 'RUVIII', 'RUVIII no ncg', 'UQ')) +
                         theme(strip.background = element_blank(),strip.text =element_blank())
                       ,x=plotscobra,y=names(plotscobra),SIMPLIFY = F) 


ggpubr::ggarrange(plotlist=plotscobra1,nrow=2,ncol=4,
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/FDRruviiimsc.png'),height = 8, width = 12)

plotscobra[[3]]+ 
  geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
  xlim(0,0.6) +
  ylim(0.3,0.93) +
  theme_minimal()+
  scale_color_brewer(palette='Set1',labels = c('RUV2', 'RUVIII', 'RUVIII no ncg', 'UQ')) +
  theme(strip.background = element_blank(),strip.text =element_blank())

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/T4FDRruviiimsc.png'),height = 5, width = 7)

# PBPS 10 popE ------------------------------------------------------------


pvalspbps10_dgim <- readRDS("/domino/datasets/local/RUV/pvalspbps10_dgim_popE.rds")
truthpbps10_dgim <- readRDS("/domino/datasets/local/RUV/truthpbps10_dgim_popE.rds")

# celltypes <- c('B', 'NK', 'T4', 'T8', 'cDC', 'cM', 'ncM', 'pDC')
celltypes <- c('B', 'NK', 'T4', 'T8', 'cDC', 'cM', 'ncM', 'pDC')

pvalits2 <- list_transpose(pvalspbps10_dgim)
pvalits2 <- pvalits2[celltypes]
pvalits2 <- lapply(pvalits2 , function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )

pvalits2 <- lapply(pvalits2 , function (x) dplyr::select(x,-RUV3PSall) )


truthits2 <- list_transpose(truthpbps10_dgim )
truthits2 <- truthits2[celltypes]
truthits2 <- lapply(truthits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )


COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
                        binary_truth = paste0("trueDE_",celltypes))
cobratoplot <- lapply(cobraperflist, prepare_data_for_plot,
                      facetted = TRUE)

plotscobra <- lapply(cobratoplot,plot_fdrtprcurve,plottype='points',pointsize=0.5)


plotscobra1 <- mapply( function(x,y) x + 
                         geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                         xlim(0,0.4) +
                         ylim(0.81,0.93) +
                         theme_minimal()+ ggtitle(y) +
                         scale_color_brewer(palette='Set1',labels = c('RUVIII', 'RUVIII+Batch', 'RUVIII PBPS',  'UQ Batch')) +
                         theme(strip.background = element_blank(),strip.text =element_blank())
                       ,x=plotscobra,y=names(plotscobra),SIMPLIFY = F) 


ggpubr::ggarrange(plotlist=plotscobra1,nrow=2,ncol=4,
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/FDRruviii_popE.png'),height = 8, width = 12)

# PBPS 10 popE int------------------------------------------------------------


pvalspbps10_dgim <- readRDS("/domino/datasets/local/RUV/pvalspbps10_dgim_popEint.rds")
truthpbps10_dgim <- readRDS("/domino/datasets/local/RUV/truthpbps10_dgim_popEint.rds")

# celltypes <- c('B', 'NK', 'T4', 'T8', 'cDC', 'cM', 'ncM', 'pDC')

celltypes <-'T4'

pvalits2 <- list_transpose(pvalspbps10_dgim)
#pvalits2 <- pvalits2[celltypes]
pvalits2 <- pvalits2['T4']
pvalits2 <- lapply(pvalits2 , function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )

#pvalits2 <- lapply(pvalits2 , function (x) dplyr::select(x,-RUV3PSall) )


truthits2 <- list_transpose(truthpbps10_dgim )
# truthits2 <- truthits2[celltypes]
truthits2 <- truthits2['T4']
truthits2 <- lapply(truthits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )


COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
                        binary_truth = paste0("trueDE_",celltypes))
cobratoplot <- lapply(cobraperflist, prepare_data_for_plot,
                      facetted = TRUE)

plotscobra <- lapply(cobratoplot,plot_fdrtprcurve,plottype='points',pointsize=0.5)


plotscobra1 <- mapply( function(x,y) x + 
                         geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                         xlim(0,0.4) +
                         ylim(0.81,0.93) +
                         theme_minimal()+ ggtitle(y) +
                         scale_color_brewer(palette='Set1',labels = c('RUVIII', 'RUVIII+Batch', 'RUVIII PBPS',  'UQ Batch')) +
                         theme(strip.background = element_blank(),strip.text =element_blank())
                       ,x=plotscobra,y=names(plotscobra),SIMPLIFY = F) 

plotscobra1

ggpubr::ggarrange(plotlist=plotscobra1,nrow=2,ncol=4,
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/FDRruviii_popEint.png'),height = 8, width = 12)





# AWS ---------------------------------------------------------------------


asw_popE <- readRDS("/domino/datasets/local/RUV/aswp_ruvdg_diffg_pop_E.rds")
awsT4PC <- asw_plot_prep(asw_popE) |> aws_plot( ct='T4') 
awsT4PC <- awsT4PC+facet_wrap("RUV", ncol=4,labeller = as_labeller(c('no_ruv'='No RUV', 'ruv2'='RUV2', 'ruv3'='RUVIII', 'ruv4'='RUV4')))


ggsave('/domino/datasets/local/RUV/paper/graphs/normalisation/allT4_silhouettep.png',height = 5, width = 7)




