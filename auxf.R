gsea <- function(genes, background, geneSets, n=10, minSize=5, name=NULL){
  ### filter background to only include genes that we assessed.
  geneSets <- geneSets[geneSets$gene_symbol %in% background,]
  m_list <- geneSets %>% split(x = .$gene_symbol, f = .$gs_name)
  # gene set must have at least minSize genes in background.
  m_list <- m_list[unlist(lapply(m_list, length)) >= minSize]
  
  overlapPval <- unlist(lapply(m_list, function(gs){
    # genes in community and gene set
    inBoth <- sum(genes %in% gs)
    # genes in community and not in gene set
    inComOnly <- length(genes) - inBoth
    # genes in background and gene set
    inGsBack <- sum(background %in% gs)
    # genes in background and not in gene set
    outGsBack <- length(background) - inGsBack
    m <- matrix(c(inBoth, inComOnly,
                  inGsBack, outGsBack),
                nrow =2, ncol=2, byrow=TRUE,
                dimnames = list(c("in community", "out community"),
                                c("in gene set", "out gene set")))
    fis <- fisher.test(m, alternative = "greater")
    pval <- fis$p.value
    return(pval)
  }))
  padj <- p.adjust(overlapPval, "fdr")
  oo <- order(overlapPval, decreasing=FALSE)
  res <- data.frame(geneSet = names(m_list)[oo[1:n]],
                    pval = overlapPval[oo[1:n]],
                    padj = padj[oo[1:n]],
                    row.names = NULL)
  return(res)
}





find_p <- function(dge){
  p.selected <- NA
  for(i.p in c(seq(0.75,0.9,0.05),seq(0.905,0.995,0.005))){
    test <- tryCatch({dge <- calcNormFactors(dge, method="upperquartile",p=i.p)}, # p=0.975
                     error=function(e) e,
                     warning=function(w) w)
    if(inherits(test,"warning")){
      next
    }else{
      p.selected <- i.p
      break
    }
  }
  return(p.selected)
}


# RUVSeq ------------------------------------------------------------------


# x = "SeqExpressionSet"
betweenLaneNormalization2.SeqES <- function(x, which=c("median","upper","full"), offset=FALSE, round=TRUE, upper.p=0.75){
  
  if(all(is.na(normCounts(x)))) {
    counts <- counts(x)
  } else {
    counts <- normCounts(x)
  }
  if(offset) {
    o = offst(x) + betweenLaneNormalization2.matrix(counts, which, offset=TRUE, round, upper.p)
  }
  else {
    o = offst(x)
  }
  newSeqExpressionSet(counts=counts(x),
                      normalizedCounts=betweenLaneNormalization2.matrix(counts, which, offset=FALSE, round, upper.p),
                      offset=o,
                      phenoData=phenoData(x), featureData=featureData(x))
  
  
}

# x = "matrix"
betweenLaneNormalization2.matrix <- function(x, which=c("median", "upper", "full"), offset=FALSE, round=TRUE, upper.p=0.75){
  which <- match.arg(which)
  if(which=="full") {
    retval <- normalizeQuantileRank(as.matrix(x), robust=TRUE)
  } else {
    if(which=="upper") {
      sum <- apply(x, 2, quantile, upper.p)
    } else {
      sum <- apply(x, 2, median)
    }
    sum <- sum/mean(sum)
    retval <- scale(x, center=FALSE, scale=sum)
  }
  if(!offset) {
    if(round) {
      retval <- round(retval)
    }
    return(retval)
  } else {
    ret <- log(retval + 0.1) - log(x + 0.1)
    return(ret)
  }
}



simData <- function(subset,pDE=0.1){
  
  sim.a <- simulateDE(counts(subset)[!rownames(subset) %in% hkname,], which_cols = samples_to_swap, prop_DE =pDE)
  subset.a <- subset
  counts(subset.a)[!rownames(subset) %in% hkname,] <- assays(sim.a)$counts
  subset.a@featureData@data$trueDE <- FALSE
  subset.a@featureData@data$trueDE[!rownames(subset) %in% hkname] <- sim.a@elementMetadata@listData[["is_DE"]] 
}


ggRLE <- function (y, covs,names,is.log=T) 
{
  if(isFALSE(is.log)){
    y <- log(y + 1) 
  }
  
  median <- apply(y, 1, median)
  rle <- apply(y, 2, function(x) x - median)
  dataplot <- as.data.frame(rle) %>% 
    pivot_longer(everything(),names_to = 'sample_cell') %>% 
    left_join(covs, by='sample_cell')
  
  dataplot$sample_cell <- fct_relevel(dataplot$sample_cell,names)
  
  p <- ggplot(dataplot,aes(x=sample_cell, y=value,colour = Processing_Cohort )) +
    geom_boxplot(outlier.shape = NA) + ylim(-2,2) +
    facet_wrap('cg_cov')+#ggtitle(ct)+
    theme_minimal()+ theme( legend.position = 'bottom',strip.text =element_blank(),
                            axis.title.x = element_blank(),
                            axis.text.x = element_blank(),
                            axis.ticks.x = element_blank()
    ) 
  
  return(p)
}





DEA <- function(subset,trt='fake.tr',Z=NULL){
  factor <- pData(subset)[,trt]
  factors <- pData(subset)[,c(trt,Z)]
  if(!is.null(Z)){
    forml <- as.formula(paste(" ~ ", paste(names(factors),collapse="+")))
    design <- model.matrix(forml,data=factors) 
  } else {
    design <- model.matrix(~ factors) 
  }
  
  y <- DGEList(counts=counts(subset), group=factor)
  y <- calcNormFactors(y,  method="upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef=2)
  res <- residuals(fit, type="deviance")
  return(list(lrt=lrt,residuals=res))
  
}


DEA.matrix <- function(subset,trt,Z=NULL){
  factor <- trt
  factors <- bind_cols(trt=factor,Z)
  if(!is.null(Z)){
    forml <- as.formula(paste(" ~ ", paste(names(factors),collapse="+")))
    design <- model.matrix(forml,data=factors) 
  } else {
    design <- model.matrix(~ factor) 
  }
  
  y <- DGEList(counts=counts(subset), group=factor)
  y <- calcNormFactors(y,  method="upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef=2)
  res <- residuals(fit, type="deviance")
  return(list(lrt=lrt,residuals=res))
  
}




TPP<- function(pval,cut,DEG){
  
  apply(pval,2, function (x) sum(names(x[x<cut]) %in%  DEG)/length(DEG))
  
}


FDP<- function(pval,cut,DEG){
  
  apply(pval,2,function (x) sum(!names(x[x<cut]) %in%  DEG)/sum(x<cut))
  
}




TPFDcurve <- function(pval,truth){
  
  cuts <- c(seq(10^(-16),0.01,length.out=1000),seq(0.01,1,length.out=5000))
  DEG <- rownames(dplyr::filter(truth,status==T))
  TPC <- t(sapply(cuts,function (x) TPP(pval=pval,cut= x,DEG=DEG) ))
  TPC <- bind_cols(TPC,cuts=cuts)
  FDC <- t(sapply(cuts,function (x) FDP(pval=pval,cut= x,DEG=DEG) ))
  FDC <- bind_cols(FDC,cuts=cuts)
  return(list(TPC=TPC,FDC=FDC))
  
}

ret.pval<- function (x){
  y <- x$PValue
  names(y) <- row.names(x)
  return(y)
}

namepvals <- function (x,name){
  
  rownames(x) <- name
  return(x)
}



FDTPplot <- function(FDCits,TPCits){
  
  FDRC <-  FDCits %>% list_rbind(names_to='it') %>%
    group_by(cuts) %>% summarise(across(everything(), list(mean)))%>% 
    dplyr::select(-it_1)%>% pivot_longer(-cuts, names_to = 'type',values_to = 'FDR')%>%
    mutate(type=gsub("_1","",type))
  
  
  
  TPRC <-  TPCits %>% list_rbind(names_to='it') %>%
    group_by(cuts) %>% summarise(across(everything(), list(mean)))%>% 
    dplyr::select(-it_1)%>% pivot_longer(-cuts, names_to = 'type',values_to = 'TPR')%>%
    mutate(type=gsub("_1","",type))
  
  curves <- left_join(FDRC,TPRC,by=c('cuts','type'))
  
  p <- ggplot(curves, aes(x=FDR,y=TPR,color=type))+geom_line()#+theme_minimal()  
  
  return(p)
  
} 


gen_normCounts <- function(x,W,center=TRUE, round=TRUE, epsilon=1, tolerance=1e-8, isLog=FALSE){
  
  
  if(!isLog && !all(.isWholeNumber(x))) {
    warning(paste0("The expression matrix does not contain counts.\n",
                   "Please, pass a matrix of counts (not logged) or set isLog to TRUE to skip the log transformation"))
  }
  
  if(isLog) {
    Y <- t(x)
  } else {
    Y <- t(log(x+epsilon))
  }
  
  alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
  correctedY <- Y - W %*% alpha
  if(!isLog && all(.isWholeNumber(x))) {
    if(round) {
      correctedY <- round(exp(correctedY) - epsilon)
      correctedY[correctedY<0] <- 0
    } else {
      correctedY <- exp(correctedY) - epsilon
    }
  }
  colnames(W) <- paste("W", seq(1, ncol(W)), sep="_")
  return(list(W = W, normalizedCounts = t(correctedY)))
}



# ASW ---------------------------------------------------------------------

ASW <- function(pca, pco, plot=FALSE, plotlab=NULL){
  
  pco <- as.integer(as.factor(pco))
  sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
  silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
  
  avg_sw2 <- round(mean(sil2[,3]),2)
  avg_swp <- round(mean(silp[,3]),2)
  
  if(isTRUE(plot)){
    p2 <- fviz_silhouette(sil2)+theme_minimal() +
      theme(axis.text.x=element_blank(),legend.position = "bottom") +
      labs(fill=plotlab,color=plotlab)
    pp <- fviz_silhouette(silp)+theme_minimal() +
      theme(axis.text.x=element_blank(),legend.position = "bottom") +
      labs(fill=plotlab,color=plotlab)
    
    return(list(s2=avg_sw2,sp=avg_swp,p2=p2,pp=pp))
    
  }
  
  return(c(avg_sw2,avg_swp))
  
}

asw_plot_prep <- function(aswpDG){
  
  swpits<- list_transpose(aswpDG)
  
  con.silp <- purrr::imap(lapply(swpits,bind_rows), ~mutate(.x, cg_cov = .y)) %>% bind_rows()
  
  con.silp <- pivot_longer(con.silp,-cg_cov, names_to ='Method', values_to="ASW")
  
  con.silp <- con.silp %>%mutate(Approach =sub(".*ruv.", "", Method), RUV= str_extract(Method, "ruv.")) %>%
    mutate(Approach = trimws(Approach)) 
  
  con.silp$RUV[is.na(con.silp$RUV)] <- 'no_ruv'
  
  return(con.silp)
  
}


aws_plot <- function (con.silp, ct){
  p <- ggplot(filter(con.silp,cg_cov==ct), aes(y= ASW,  color=Approach, group=Approach)) + 
    geom_boxplot() + facet_wrap('RUV') + theme_minimal() + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")
  return(p)
}




#### extra RUVSeq #####

makeGroups <- function(xs) {
  xs <- factor(xs)
  groups <- matrix(-1, nrow = length(levels(xs)), ncol = max(table(xs)))
  for (i in 1:length(levels(xs))) {
    idxs <- which(xs == levels(xs)[i])
    groups[i,1:length(idxs)] <- idxs
  }
  groups
}

.isWholeNumber <- function(x, tol = .Machine$double.eps^0.5) {
  !is.na(x) & abs(x - round(x)) < tol
}



# ruv logy ----------------------------------------------------------------

limma.matrix <- function(y, design, gene.D){
  colnames(design)[2] <- 'tr'
  v <- voom(y, design, plot=F) 
  vfit <- lmFit(v, design)
  efit <- eBayes(vfit)
  
  etable <- topTable(efit,number=dim(y)[1] , coef="tr", adjust="BH")
  etable <- etable[gene.D$gene,]
  etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
  etable$trueDEG <- gene.D$trueDE
  
  return(etable)
  
}

limma.matrixint <- function(design, vfit, gene.D){
  contrast.matrix <- makeContrasts(
    SLEeuro =tr,
    SLEasian = interaction,
    levels = colnames(design))
  
  vfit <- contrasts.fit(vfit, contrast.matrix)
  efit <- eBayes(vfit)
  etable <- topTable(efit, number=dim(gene.D)[1], adjust="BH")
  etable <- etable[gene.D$gene,]
  etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
  etable$trueDEG <- gene.D$trueDE
  
  return(etable)
  
}





RUVIIIW <- function (Y, M, ctl, k = NULL, eta = NULL, include.intercept = TRUE, 
                     average = FALSE, fullalpha = NULL, return.info = FALSE, inputcheck = TRUE) 
{
  if (is.data.frame(Y)) 
    Y = data.matrix(Y)
  m = nrow(Y)
  n = ncol(Y)
  M = replicate.matrix(M)
  #ctl = tological(ctl, n)
  if (inputcheck) {
    if (m > n) 
      warning("m is greater than n!  This is not a problem itself, but may indicate that you need to transpose your data matrix.  Please ensure that rows correspond to observations (e.g. microarrays) and columns correspond to features (e.g. genes).")
    if (sum(is.na(Y)) > 0) 
      warning("Y contains missing values.  This is not supported.")
    if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) > 
        0) 
      warning("Y contains infinities.  This is not supported.")
  }
  Y = RUV1(Y, eta, ctl, include.intercept = include.intercept)
  if (ncol(M) >= m) 
    newY = Y
  else if (is.null(k)) {
    ycyctinv = solve(Y[, ctl] %*% t(Y[, ctl]))
    newY = (M %*% solve(t(M) %*% ycyctinv %*% M) %*% (t(M) %*% 
                                                        ycyctinv)) %*% Y
    fullalpha = NULL
  }
  else if (k == 0) {
    newY = Y
    fullalpha = NULL
  }
  else {
    if (is.null(fullalpha)) {
      Y0 = residop(Y, M)
      fullalpha = t(svd(Y0 %*% t(Y0))$u[, 1:min(m - ncol(M), 
                                                sum(ctl)), drop = FALSE]) %*% Y
    }
    alpha = fullalpha[1:min(k, nrow(fullalpha)), , drop = FALSE]
    ac = alpha[, ctl, drop = FALSE]
    W = Y[, ctl] %*% t(ac) %*% solve(ac %*% t(ac))
    newY = Y - W %*% alpha
  }
  if (average) 
    newY = ((1/apply(M, 2, sum)) * t(M)) %*% newY
  if (!return.info) 
    return(newY)
  else return(list(newY = newY, fullalpha = alpha, W = W))
}


RUV2mod <- function (Y, X, ctl, k, Z=1, include.intercept=T) 
{
  res <- RUV2(Y = Y, X = X, ctl = ctl, k = k, do_projectionplot = FALSE,Z=Z, include.intercept= include.intercept) 
  W <- res$W
  alpha <- res$alpha
  newY <- Y - W %*% alpha
  
  return (list(newY = newY, fullalpha = alpha, W = W))
}

RUV4mod <- function (Y, X, ctl, k, Z=1, include.intercept=T) 
{
  res <- RUV4(Y = Y, X = X, ctl = ctl, k = k, Z = Z, include.intercept= include.intercept) 
  W <- res$W
  alpha <- res$alpha
  newY <- Y - W %*% alpha
  
  return (list(newY = newY, fullalpha = alpha, W = W))
}


pseudosamp <- function(sc_samp, n=10, seed = 1){
  
  psscounts <- GetAssayData(object = sc_samp, assay = "originalexp", slot = "counts")
  cell_names <- colnames(sc_samp)
  cell_types <- sc_samp$cg_cov
  cell_names_list <- split(cell_names,  cell_types)
  cell_ab <- sapply(cell_names_list, length)
  cellt_filter <- cell_names_list[names(which(cell_ab>100))]
  cell_ab <- sapply(cellt_filter, length)
  #noise_ab <- rnorm() not needed since sum will also differ?
  set.seed(seed)
  pseudo_samp <- replicate(n,list(mapply(function (x,y) 
    base::sample(x=x,size=y,replace=T),cellt_filter, cell_ab, SIMPLIFY = F)))
  names(pseudo_samp) <-paste0('ps',1:n)
  pseudo_samp[['orig']] <- cellt_filter
  pbps <- map_depth(pseudo_samp,2, function(x) rowSums(psscounts [,x]))
  
  
  pbps.m <- as.matrix(suppressMessages(bind_cols(pbps)))
  rownames(pbps.m) <- rownames(psscounts)
  pdata <-  merge(names(cell_ab),names(pseudo_samp))
  colnames(pbps.m) <- paste0(pdata$x,pdata$y)
  colnames(pdata)<- c('cg_cov','ind_cov')
  rownames(pdata) <- colnames(pbps.m)
  
  set <- newSeqExpressionSet(pbps.m,
                             phenoData = pdata)
  
  return(set)
}


# PBPS --------------------------------------------------------------------

cellspbps <- function (subsample, pb_id='sample_cell',cell_id='cell_id', seed = 1,n=1){
  
  
  # compute the average number of cells per biological group
  total_cell <- subsample %>% group_by(across(all_of(pb_id))) %>% summarise(ncell=n())
  subsample <- left_join(subsample,total_cell)
  avg_cell <- group_by(subsample,bsg)%>%summarise(avg_ncell=mean(ncell))
  avg_cell <- dplyr::select(subsample,group,bsg)%>% unique()%>% left_join(avg_cell,by='bsg')
  rownames(avg_cell) <- avg_cell$group
  
  # select the cells from each group
  groups <- split(subsample[,cell_id],subsample$group)
  avg_cell <- avg_cell[names(groups),]
  set.seed(seed)
  
  #sampling with replacement
  pseudo_samp <- replicate(n,list(mapply(function (x,y) 
    base::sample(x=x,size=round(y),replace=T),groups, avg_cell$avg_ncell, SIMPLIFY = F)))
  
  names(pseudo_samp) <- paste0('pa',1:n)
  pseudo_samp <- unlist(pseudo_samp,recursive=FALSE)
  
  return(pseudo_samp)
  
}


pbps <- function(pseudo_samp,Gr, ds , cell_id='cell_id'){
  
  
  ds$orig.index <- ds@meta.data[[cell_id]] %in% pseudo_samp[[Gr]]
  Idents(object = ds) <- "orig.index"
  origgcounts <- GetAssayData(object = subset(x = ds, idents = TRUE), assay = "originalexp", slot = "counts")
  psscounts <- origgcounts[,pseudo_samp[[Gr]]]### cell names might be repeated, this does the trick
  pssamp <- rowSums(psscounts)
  return(pssamp)
}

gen_PBPS <- function(ds, PBC, ctype, BioVar, NVar, id_pb, id_sub, cell_id, n=2, Seed=2 ){
  
  groups_info <- group_by(PBC$samples,across(all_of(c(BioVar,NVar))))%>%
    summarise(nsample=n()) %>% mutate(bsg = paste(!!!syms(BioVar), sep = "_"),
                                      group=paste(bsg,!!!syms(NVar),sep='_'))
  
  bsg.keep <- names(table(groups_info$bsg)[table(groups_info$bsg)>1]) # remove bsg with only 1 group
  
  PBC$samples <- mutate(PBC$samples, bsg = paste(!!!syms(BioVar), sep = "_"),
                        group=paste(bsg,!!!syms(NVar),sep='_')) # add bsg and group info
  
  sc.ref <- ds@meta.data %>%mutate(bsg = paste(!!!syms(BioVar), sep = "_"),
                                   group=paste(bsg,!!!syms(NVar),sep='_')) %>% filter(bsg %in% bsg.keep)
  
  cells_to_pbps <- cellspbps(subsample = sc.ref, seed = Seed, n = n)
  
  pbpscounts <- sapply(names(cells_to_pbps),function(x)pbps(cells_to_pbps,Gr=x,ds=ds))
  
  pbps_info <- sapply(colnames(pbpscounts), function (x) regmatches(x, regexpr("\\.", x), invert = TRUE))
  pa <- sapply(pbps_info,function (x) x[1])
  pa_group <- sapply(pbps_info,function (x) x[2])
  
  pbps_info <- tibble(!!id_pb:= colnames(pbpscounts), pa = pa, group = pa_group) %>%
    left_join(groups_info, by='group') %>% select(-nsample)
  
  miss_var <- colnames(PBC$samples)[!colnames(PBC$samples) %in% colnames(pbps_info) ]
  otherv <- sapply(miss_var,function(x) pbps_info[,x ]<- rep(NA,dim(pbps_info)[1]))
  pbps_info <- cbind(pbps_info,otherv) %>% mutate(pbps=1)
  
  sifull <- mutate(PBC$samples,pbps=0, pa='orig') %>% rbind(pbps_info)
  sifull$ind_cov <- as.character(sifull$ind_cov)
  
  sifull[sifull$pbps==1,] <- sifull[sifull$pbps==1,] %>% mutate(!!id_sub:= paste(!!!syms(setdiff(BioVar, ctype)), sep = "_"))
  fullcount <- cbind(PBC$counts,pbpscounts[rownames(PBC),])
  
  PBPSC <- DGEList(counts = fullcount, samples = dplyr::select(sifull, -lib.size, -norm.factors))
  
  return(PBPSC)
  
}

# r2 ----------------------------------------------------------------------

r2 <- function (W, Y){
  fit <- lm(W ~ Y)
  res <-sapply(summary(fit), function (x) x$r.squared)
  return (res)
}


