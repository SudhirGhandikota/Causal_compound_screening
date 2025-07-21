#library(CWGCNA)
# https://github.com/yuabrahamliu/CWGCNA
setwd('CWGCNA/')
source('CWGCNA_methods_diff.R')
source('CWGCNA_methods_med.R')
library(doParallel)
library(limma)
getwd()
#Some parallel computing going on in the background that is not getting cleaned up
#fully between runs can cause Error in summary.connection(connection) : invalid
#connection. The following function is needed to be called to fix this error.
unregister_dopar <- function(){
  
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
  
}

anovaplot <- function(restab, 
                     plottype = 'pval', 
                     titlesufix = NULL, 
                     samplesize = NULL, 
                     titlesize = 15, 
                     textsize = 13, 
                     face = 'bold'){
  
  #library(ggplot2)
  
  if(!is.null(titlesufix))
    titlesufix <- paste0(' ', titlesufix)
  
  resmean <- colMeans(restab)
  resmean <- as.matrix(resmean)
  resmean <- t(resmean)
  row.names(resmean) <- 1:nrow(resmean)
  
  if(plottype == 'MSS'){
    resmean <- resmean[, grepl(pattern = '_MSS', colnames(resmean)), drop = FALSE]
    colnames(resmean) <- gsub(pattern = '_MSS', replacement = '', x = colnames(resmean))
    ytitle <- 'MSS'
  }else if(plottype == 'pval'){
    resmean <- resmean[, grepl(pattern = '_pval', colnames(resmean)), drop = FALSE]
    colnames(resmean) <- gsub(pattern = '_pval', replacement = '', x = colnames(resmean))
    resmean <- -log10(resmean)
    ytitle <- '-log10(pval)'
  }else{
    resmean <- resmean[, grepl(pattern = '_F', colnames(resmean)), drop = FALSE]
    colnames(resmean) <- gsub(pattern = '_F', replacement = '', x = colnames(resmean))
    ytitle <- 'F statistic'
  }
  
  resmean <- t(resmean)
  colnames(resmean) <- 'Statistic'
  resmean <- as.data.frame(resmean, stringsAsFactors = FALSE)
  
  resmean$Factor <- factor(row.names(resmean), levels = row.names(resmean), 
                           ordered = TRUE)
  row.names(resmean) <- 1:nrow(resmean)
  
  if(!is.null(samplesize))
    subtitle <- paste0('(sample size = ', samplesize, ')')
  
  p <- ggplot2::ggplot(data = resmean, 
                       mapping = ggplot2::aes(x = Factor, 
                                              y = Statistic, fill = Factor))
  p <- p + ggplot2::geom_bar(stat = 'identity') + 
    ggplot2::ggtitle(paste0('Type 3 ANOVA', titlesufix, ' ', ytitle), subtitle = subtitle) + 
    ggplot2::ylab(paste0('Averaged feature ', ytitle)) + 
    ggplot2::xlab('') + 
    ggplot2::scale_fill_discrete(guide = 'none') + 
    ggplot2::theme_bw() + 
    ggplot2::theme(axis.title.y = ggplot2::element_text(size = titlesize, face = face),
                   axis.text.x = ggplot2::element_text(face = face, size = textsize,
                                                       angle = 45, hjust = 1),
                   axis.text.y = ggplot2::element_text(face = face, size = textsize - 3),
                   plot.title = ggplot2::element_text(size = titlesize, face = face), 
                   plot.subtitle = ggplot2::element_text(size = textsize, face = face))
  
  if(plottype == 'pval')
    p <- p + ggplot2::geom_hline(yintercept = -log10(0.05), col = 'red')
  
  print(p)
  
}

featuresampling <- function(betas, 
                            topfeatures = 10000, 
                            pddat = NULL, 
                            variancetype = 'sd', 
                            anova = FALSE, 
                            responsevarname = NULL, 
                            filternovar = TRUE, 
                            threads = 1, 
                            plotannovares = TRUE, 
                            titlesize = 15, 
                            textsize = 13, 
                            face = 'bold', 
                            featuretype = NULL, 
                            plottitlesuffix = NULL){
  
  if(is.matrix(betas) == FALSE){
    
    rownum <- nrow(betas)
    colnum <- ncol(betas)
    
    betascolnames <- colnames(betas)
    betasrownames <- row.names(betas)
    
    betas <- as.matrix(betas)
    
    if(nrow(betas) == rownum & ncol(betas) == colnum){
      
      row.names(betas) <- betasrownames
      colnames(betas) <- betascolnames
      
      warning('The parameter `betas` needs a matrix, but the current one is not, \n
            so it has been converted to a matrix.\n')
      
    }else{
      
      warning('The parameter `betas` needs a matrix, but the current one is not, \n
            please convert it to a matrix.\n')
      
      return(NULL)
      
    }
    
  }
  
  if(topfeatures <= 0){
    topfeatures <- 1000
  }
  
  if(topfeatures > nrow(betas)){
    topfeatures <- nrow(betas)
  }
  
  if(variancetype == 'mad'){
    
    calvar <- function(line){
      
      linemedian <- median(line, na.rm = TRUE)
      res <- median(abs(line - linemedian), na.rm = TRUE)
      names(res) <- 'MAD'
      
      return(res)
      
    }
    
  }else{
    
    calvar <- function(line){
      
      res <- sd(line, na.rm = TRUE)
      names(res) <- 'SD'
      
      return(res)
      
    }
    
    
  }
  
  if(anova == TRUE & (!is.null(pddat))){
    betas <- betas[,pddat[,1], drop = FALSE]
    
    if(filternovar == TRUE){
      
      filteridx <- apply(X = betas, MARGIN = 1, FUN = sd, na.rm = TRUE)
      
      betas <- betas[filteridx != 0, , drop = FALSE]
      
    }
    
    calvar <- function(line, pds = pddat){
      
      #library(car)
      newdata <- pds[, -1, drop = FALSE]
      pdnames <- colnames(newdata)
      newdata$beta <- line
      formstr <- paste0(pdnames, collapse = ' + ')
      formstr <- paste0('beta ~ ', formstr)
      formstr <- as.formula(formstr)
      fit <- lm(formstr, data = newdata)
      
      aovfit <- car::Anova(fit, type = 3, singular.ok = TRUE)
      F <- aovfit$`F value`
      F <- F[2:(length(F)-1)]
      names(F) <- pdnames
      F <- as.data.frame(F, stringsAsFactors = FALSE)
      F <- as.data.frame(t(F))
      row.names(F) <- 1
      colnames(F) <- paste0(colnames(F), '_F')
      
      pvals <- aovfit$`Pr(>F)`
      pvals <- pvals[2:(length(pvals)-1)]
      names(pvals) <- pdnames
      pvals <- as.data.frame(pvals, stringsAsFactors = FALSE)
      pvals <- as.data.frame(t(pvals))
      row.names(pvals) <- 1
      colnames(pvals) <- paste0(colnames(pvals), '_pval')
      
      SS <- aovfit$`Sum Sq`
      DF <- aovfit$Df
      MSS <- SS/DF
      
      MSS <- MSS[2:(length(MSS)-1)]
      names(MSS) <- pdnames
      MSS <- as.data.frame(MSS, stringsAsFactors = FALSE)
      MSS <- as.data.frame(t(MSS))
      row.names(MSS) <- 1
      colnames(MSS) <- paste0(colnames(MSS), '_MSS')
      
      lineres <- cbind(F, pvals, MSS)
      row.names(lineres) <- row.names(line)
      
      return(lineres)
      
    }
    
  }
  
  if(threads == 1){
    varlist <- list()
    for(i in 1:nrow(betas)){
      varlist[[i]] <- calvar(line = betas[i,])
    }
    
  }else{
    
    #library(doParallel)
    
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(min(threads, cores))
    
    doParallel::registerDoParallel(cl)
    
    #date()
    `%dopar%` <- foreach::`%dopar%`
    varlist <- foreach::foreach(i = 1:nrow(betas), 
                                .export = NULL) %dopar% {
                                  calvar(line = betas[i,])
                                }
    
    parallel::stopCluster(cl)
    
    unregister_dopar()
    
    
  }
  
  varres <- do.call(rbind, varlist)
  
  if(is.null(dim(varres))){
    
    return(NULL)
    
  }
  
  row.names(varres) <- row.names(betas)
  
  if(anova == TRUE){
    
    if(is.null(responsevarname)){
      if(!is.null(pddat)){
        responsevarname <- names(pddat)[2]
      }
      
    }
    
    if(!is.null(responsevarname)){
      varres <- varres[order(varres[,paste0(responsevarname, '_pval')]), , drop = FALSE]
    }
    
  }else{
    if(variancetype == 'mad'){
      varres <- varres[order(-varres[,'MAD']), , drop = FALSE]
    }else{
      varres <- varres[order(-varres[,'SD']), , drop = FALSE]
    }
  }
  
  
  topvarres <- varres[1:topfeatures, , drop = FALSE]
  
  
  topbetas <- betas[(row.names(betas) %in% row.names(topvarres)), , drop = FALSE]
  
  #topbetas <- betas[1:topfeatures,]
  
  res <- list(betas = topbetas, 
              varicanceres = varres)
  
  if(anova == TRUE & (!is.null(pddat))){
    
    if(is.null(featuretype)){
      featuretype <- 'feature'
    }
    
    if(plotannovares == TRUE){
      
      if(!is.null(plottitlesuffix)){
        plottitlesuffix <- paste0('for ', plottitlesuffix, ' top ')
      }else{
        plottitlesuffix <- 'for top '
      }
      
      anovaplot(restab = topvarres, plottype = 'pval', 
                titlesufix = paste0(plottitlesuffix, topfeatures, ' ', featuretype, 's'), 
                titlesize = titlesize, textsize = textsize, face = face, 
                samplesize = ncol(betas))
      anovaplot(restab = topvarres, plottype = 'MSS', 
                titlesufix = paste0(plottitlesuffix, topfeatures, ' ', featuretype, 's'), 
                titlesize = titlesize, textsize = textsize, face = face, 
                samplesize = ncol(betas))
      anovaplot(restab = topvarres, plottype = 'F', 
                titlesufix = paste0(plottitlesuffix, topfeatures, ' ', featuretype, 's'), 
                titlesize = titlesize, textsize = textsize, face = face, 
                samplesize = ncol(betas))
      
    }
    
    
    
  }
  
  return(res)
  
  
}


maphuecolors <- function(numericvec){
  
  uninumericvec <- unique(numericvec)
  uninumericvec <- uninumericvec[order(uninumericvec)]
  uninumericvec <- unique(c(0, uninumericvec))
  unicolorvec <- c('#C0C0C0', scales::hue_pal()(length(uninumericvec) - 1))
  names(unicolorvec) <- as.character(uninumericvec)
  colorvec <- unicolorvec[as.character(numericvec)]
  colorvec <- as.vector(colorvec)
  return(colorvec)
  
}

runWGCNA = function(dat,
                    filterdat=FALSE,
                    filternum = 5000,
                    varicancetype = 'sd',
                    powers = seq(1, 20, 1),
                    rsqcutline = 0.8,
                    minclustersize = 50, 
                    mergecutheight = 0.2, 
                    deepSplit = 2,
                    TOMType = 'signed',
                    maxblocksize = 5000, 
                    returntom = FALSE, 
                    plot = FALSE, 
                    featuretype = 'gene', 
                    plottitleprefix = NULL, 
                    threads = 0, 
                    seed = 2022){
  
  if(filterdat){
    
    filternum <- min(filternum, nrow(dat))
    dat <- featuresampling(betas = dat,
                           topfeatures = filternum,
                           pddat = NULL,
                           variancetype = varicancetype,
                           anova = FALSE,
                           responsevarname = NULL,
                           threads = threads)
    dat <- dat$betas
    
  }
  
  print('****** Scale free analysis ******')
  sftpowers <- choosepower(dat = dat,
                           powers = powers,
                           rsqcutline = rsqcutline,
                           plot = FALSE)
  print(paste('*** Power Estimate', as.vector(unlist(sftpowers$powerEstimate)), '***'))
  if(!is.numeric(mergecutheight))
    mergecutheight <- WGCNA::dynamicMergeCut(n = ncol(dats[[1]]), mergeCor = 0.8)

  print('****** Running WGCNA ******')
  cor <- WGCNA::cor
  
  stamp <- Sys.time()
  stamp <- gsub(pattern = ' ', replacement = '_', x = stamp)
  stamp <- gsub(pattern = ':', replacement = '-', x = stamp)
  stamp <- paste0('blockwiseTOM.', stamp)
  
  wgcnares <-  WGCNA::blockwiseModules(datExpr = t(dat), 
                                       power = as.vector(unlist(sftpowers$powerEstimate)),
                                       TOMType = TOMType, 
                                       minModuleSize = minclustersize,
                                       deepSplit = deepSplit, 
                                       reassignThreshold = 0, 
                                       mergeCutHeight = mergecutheight,
                                       numericLabels = TRUE, 
                                       pamRespectsDendro = FALSE,
                                       saveTOMs = returntom, 
                                       saveTOMFileBase = stamp, 
                                       verbose = 3, 
                                       nThreads = threads, 
                                       maxBlockSize = maxblocksize, 
                                       randomSeed = seed)
  
  if(returntom == TRUE){
    
    if(length(wgcnares$TOMFiles) == 1){
      tomfile <- wgcnares$TOMFiles
      env2load <- environment()
      loadtom <- load(file = tomfile, envir = env2load)
      tom <- get(x = loadtom, envir = env2load)
      
      tom <- as.matrix(tom)
      colnames(tom) <- row.names(tom) <- 
        names(wgcnares$colors)[wgcnares$blockGenes[[1]]]
      
      unlink(tomfile)
      
    }else{
      tom <- list()
      i <- 1
      for(i in 1:length(wgcnares$TOMFiles)){
        tomfile <- wgcnares$TOMFiles[i]
        env2load <- environment()
        loadtom <- load(file = tomfile, envir = env2load)
        tom[[i]] <- get(x = loadtom, envir = env2load)
        
        tom[[i]] <- as.matrix(tom[[i]])
        colnames(tom[[i]]) <- row.names(tom[[i]]) <- 
          names(wgcnares$colors)[wgcnares$blockGenes[[i]]]
        
        unlink(tomfile)
      }
    }
    
  }
  MEs <- wgcnares$MEs
  numericlabels <- wgcnares$colors
  
  if(plot == TRUE){
    #colorlabels <- WGCNA::labels2colors(labels = numericlabels)
    colorlabels <- maphuecolors(numericvec = numericlabels)
    names(colorlabels) <- names(numericlabels)
    
    if(is.null(plottitleprefix)){
      
      captialstr <- paste0(toupper(substr(featuretype, 1, 1)), 
                           substr(featuretype, 2, nchar(featuretype)))
      
      plottitleprefix <- paste0(captialstr, ' dendrogram and module colors')
    }else{
      plottitleprefix <- paste0(plottitleprefix, ' ', featuretype, 
                                ' dendrogram and module colors')
    }
    
    if(length(wgcnares$dendrograms) == 1){
      
      print(
        WGCNA::plotDendroAndColors(wgcnares$dendrograms[[1]], 
                                   colorlabels, 
                                   groupLabels = 'Module colors', 
                                   dendroLabels = FALSE, 
                                   hang = 0.03,
                                   addGuide = TRUE, 
                                   guideHang = 0.05, 
                                   main = plottitleprefix)
        
      )
    }else{
      
      # Plot the dendrogram and the module colors underneath for each block
      for(j in 1:length(wgcnares$dendrograms)){
        
        print(
          WGCNA::plotDendroAndColors(wgcnares$dendrograms[[j]], 
                                     colorlabels[wgcnares$blockGenes[[j]]], 
                                     groupLabels = 'Module colors', 
                                     dendroLabels = FALSE, 
                                     hang = 0.03,
                                     addGuide = TRUE, 
                                     guideHang = 0.05,
                                     main = paste0(plottitleprefix, ' in block ', j))
        )
      }
    }
    
  }
  
  res <- list(mergedmelist = MEs, numericlabels = numericlabels)
  if(returntom == TRUE)
    res$tom = tom
  if(plot == TRUE)
    res$colorlabels = colorlabels
  if(filterdat)
    res$dat = dat
  
  return(res)
}

#data in datlist is a matrix using features as rows and samples as columns
choosepower = function(dat,
                        powers = seq(1, 20, 1), 
                        rsqcutline = 0.8, 
                        plot = FALSE, 
                        titlesize, 
                        textsize, 
                        face){
  
  dat <- t(dat)
  sft <- WGCNA::pickSoftThreshold(data = dat, 
                                  powerVector = powers, 
                                  RsquaredCut = rsqcutline, 
                                  verbose = 5)
  
  if(is.na(sft$powerEstimate)){
    
    if(max(sft$fitIndices$SFT.R.sq) > 0.5){
      sft$powerEstimate <- sft$fitIndices$Power[match(max(sft$fitIndices$SFT.R.sq), 
                                                      sft$fitIndices$SFT.R.sq)]
      warning(paste0('No power with a soft R square > ', rsqcutline, 
                     ', use soft R square = ', max(sft$fitIndices$SFT.R.sq), ' instead\n'))
    }else{
      return(NULL)
    }
    
  }
  
  if(plot == TRUE){
    sftdat <- data.frame(powers = rep(sft$fitIndices[,1], 2), 
                         vals = c(-sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
                                  sft$fitIndices[,5]), 
                         metrics = rep(c('Signed R Square', 'Mean Connectivity'), 
                                       each = nrow(sft$fitIndices)), 
                         stringsAsFactors = FALSE)
    
    sftdat$metrics <- factor(x = sftdat$metrics, levels = c('Signed R Square', 
                                                            'Mean Connectivity'), 
                             ordered = TRUE)
    p <- ggplot2::ggplot(sftdat, ggplot2::aes(x = powers, y = vals))
    
    print(
      p + ggplot2::geom_point(size = 2, color = scales::hue_pal()(1)) + 
        ggplot2::geom_line(size = 1, color = scales::hue_pal()(1)) + 
        ggplot2::xlab('Soft Threshold (Power)') + ggplot2::ylab('Value') + 
        ggplot2::ggtitle('Scale free topology for soft-thresholding') + 
        ggplot2::geom_vline(xintercept = sft$powerEstimate, color = 'blue', size = 1.5) + 
        ggplot2::facet_wrap(ggplot2::vars(metrics), scales = 'free_y', ncol = 2) + 
        ggplot2::theme_bw() + 
        ggplot2::theme(plot.title = ggplot2::element_text(size = titlesize, face = face), 
                       legend.title = ggplot2::element_text(size = textsize, face = face), 
                       legend.text = ggplot2::element_text(size = textsize, face = face), 
                       axis.text = ggplot2::element_text(size = textsize, face = face), 
                       axis.title = ggplot2::element_text(size = textsize, face = face), 
                       strip.text = ggplot2::element_text(size = textsize, face = face))
    )
  }
  return(sft)
}


balancesampling <- function(dat,
                            pddat, 
                            responsename, 
                            nround = 10,
                            seed = 2022, 
                            k = 5, 
                            threads = 1, 
                            samplingmethod = 'updn'){
  
  dat <- t(dat)
  
  if(is.null(responsename) & ncol(pddat) >= 2)
    responsename <- names(pddat)[2]
  
  if(is.null(responsename))
    return(NULL)
  
  sampleidmapping <- pddat[,c('sampleid', responsename)]
  names(sampleidmapping) <- c('sampleid', 'responsetype')
  dat <- dat[sampleidmapping$sampleid, , drop = FALSE]
  
  responsetypes <- unique(sampleidmapping$responsetype)
  print(responsetypes)

  samplingpds <- function(samplingsamples, 
                          pddat){
    samples <- gsub(pattern = '\\.[0-9]*$', replacement = '', x = samplingsamples)
    samplingpds <- pddat[match(samples, pddat$sampleid), , drop = FALSE]
    samplingpds$sampleid <- samplingsamples
    row.names(samplingpds) <- 1:nrow(samplingpds)
    return(samplingpds)
  }
  
  createnames <- function(rawnames){
    
    reversenames <- FALSE
    if(sum(grepl(pattern = '-', x = rawnames)) > 0){
      
      rawnames <- gsub(pattern = '-', replacement = '___', x = rawnames)
      reversenames <- TRUE
    }
    newnames <- make.names(names = rawnames, unique = TRUE)

    suffix <- gsub(pattern = '^.*\\.', replacement = '', x = newnames)
    suffix[suffix %in% rawnames] <- ''
    suffix <- paste0('.', suffix)
    suffix[suffix == '.'] <- ''
    newnames <- paste0(rawnames, suffix)
    
    if(reversenames == TRUE)
      newnames <- gsub(pattern = '___', replacement = '-', x = newnames)
    
    return(newnames)
  }
  
  singleroundsampling <- function(responsetypes,
                                  sampleidmapping,
                                  i = 2022,
                                  dat, 
                                  k = 5, 
                                  samplingmethod = 'updn'){
    if(samplingmethod == 'up')
      targetnum <- max(table(sampleidmapping$responsetype))
    else if(samplingmethod == 'dn')
      targetnum <- min(table(sampleidmapping$responsetype))
    else
      targetnum <- ceiling(nrow(sampleidmapping)/length(responsetypes))
    j = 1
    for(j in 1:length(responsetypes)){
      ResponseType <- responsetypes[j]
      sub <- subset(sampleidmapping, responsetype == ResponseType)
      subnum <- nrow(sub)
      
      if(subnum >= targetnum & subnum >= 1){
        sampleids <- sub$sampleid
        
        set.seed(i)
        samples <- sample(x = sampleids,
                          size = targetnum,
                          replace = FALSE)
        #Set replace to FALSE (random sampling), not TRUE (bootstrapping), 
        #because if using bootstrapping, the downstream limma diff feature 
        #analysis will return many false positive differential features 
        #due to the repeated features introduced by this bootstrapping step

        datsub <- dat[samples, , drop = FALSE]
        row.names(datsub) <- createnames(rawnames = row.names(datsub))
        
      }
      else if(subnum < targetnum & subnum > 1){
        subvar <- dat[sub$sampleid, , drop = FALSE]
        
        knnres <- dbscan::kNN(x = subvar, k = min(k, nrow(subvar) - 1))
        knnres <- knnres$id
        knnresidx <- seq(1, nrow(knnres)*ncol(knnres), by = 1)
        knnresidxmat <- matrix(data = knnresidx, 
                               nrow = nrow(knnres), 
                               byrow = FALSE)
        row.names(knnresidxmat) <- 1:nrow(knnresidxmat)
        colnames(knnresidxmat) <- 1:ncol(knnresidxmat)
        
        set.seed(i)
        sampledidx <- sample(x = knnresidx, size = targetnum - subnum, replace = TRUE)
        
        smoteidx <- do.call(rbind, 
                            lapply(X = sampledidx, 
                                   FUN = function(x){which(knnresidxmat == x, arr.ind = TRUE)}))
        row.names(smoteidx) <- 1:nrow(smoteidx)
        
        set.seed(i)
        zetas <- runif(nrow(smoteidx))
        
        synthesizeddat <- lapply(X = seq(1:nrow(smoteidx)), 
                                 FUN = function(x){subvar[row.names(knnres)[smoteidx[x, 1]],] + 
                                     zetas[x]*(subvar[knnres[smoteidx[x, 1], smoteidx[x, 2]],] - 
                                                 subvar[row.names(knnres)[smoteidx[x, 1]],])})
        names(synthesizeddat) <- row.names(knnres)[smoteidx[,1]]
        synthesizeddat <- do.call(rbind, synthesizeddat)
        datsub <- rbind(subvar, synthesizeddat)
        
        row.names(datsub) <- createnames(rawnames = row.names(datsub))
        
      }
      else if(subnum == 1){
        
        subvar <- dat
        knnres <- dbscan::kNN(x = subvar, k = min(k, nrow(subvar) - 1))
        
        mindist <- min(knnres$dist[knnres$dist > 0])
        mindistidx <- which(knnres$dist == mindist, arr.ind = TRUE)
        
        difference <- subvar[knnres$id[mindistidx[1, 1], mindistidx[1, 2]],] - 
          subvar[row.names(knnres$id)[mindistidx[1, 1]],]
        
        set.seed(i)
        zetas <- runif(targetnum - subnum)
        
        synthesizeddat <- lapply(X = 1:length(zetas), 
                                 FUN = function(x){subvar[sub$sampleid[1],] + zetas[x]*difference})
        
        names(synthesizeddat) <- rep(sub$sampleid[1], length(synthesizeddat))
        
        synthesizeddat <- do.call(rbind, synthesizeddat)
        datsub <- rbind(subvar[sub$sampleid[1], , drop = FALSE], synthesizeddat)
        row.names(datsub) <- createnames(rawnames = row.names(datsub))
      }
      if(j == 1)
        datsubs <- datsub
      else
        datsubs <- rbind(datsubs, datsub)
    }
    datsubs <- t(datsubs)
    return(datsubs)
    
  }
  
  
  if(threads == 1){
    
    sampleddats <- list()
    for(i in 1:nround){
      tmpdat <- singleroundsampling(responsetypes = responsetypes,
                                    sampleidmapping = sampleidmapping,
                                    i = seed + i - 1,
                                    dat = dat, 
                                    k = k, 
                                    samplingmethod = samplingmethod)
      tmppd <- samplingpds(samplingsamples = colnames(tmpdat), 
                           pddat = pddat)
      
      sampleddats[[i]] <- tmpdat
      
    }
  
  }
  else{
    iseqs <- 1:nround
    
    #library(doParallel)
    
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(min(threads, cores))
    
    doParallel::registerDoParallel(cl)
    
    #date()
    `%dopar%` <- foreach::`%dopar%`
    sampleddats <- foreach::foreach(i = iseqs,
                                    #.export = ls(name = globalenv())) %dopar% {
                                    .export = NULL) %dopar% {
                                      singleroundsampling(responsetypes = responsetypes,
                                                          sampleidmapping = sampleidmapping,
                                                          i = seed +i - 1,
                                                          dat = dat, 
                                                          k = k, 
                                                          samplingmethod = samplingmethod)
                                    }
    
    parallel::stopCluster(cl)
    
    unregister_dopar()
    
  }
  
  sampledpds <- list()
  for(i in 1:nround){
    
    sampledpds[[i]] <- samplingpds(samplingsamples = colnames(sampleddats[[i]]), 
                                   pddat = pddat)
    
  }
  
  
  names(sampleddats) <- names(sampledpds) <- paste0('round_', 1:nround)
  
  res <- list(dats = sampleddats, 
              pds = sampledpds)
  
  return(res)
  
}


diffmodules = function(dats, 
                       pddats,
                       wgcnares,
                       
                       responsevarname = NULL, 
                       responselevels = NULL, 
                       confoundings = NULL,
                       
                       balanceadj=FALSE,
                       samplingmethod = 'updn', 
                       nround = 10,
                       seed = 2022,
                       k = 5, 
                       threads = 1, 
                       method = 'Cauchy',
                       
                       plot=FALSE,
                       diffcutoff=0,
                       pvalcutoff=0.05,
                       pvalcolname='adj.P.Val',
                       
                       titleprefix = NULL,
                       titlesize = 14,
                       textsize = 12,
                       face = 'bold'){
  
  #limma on module level
  print('***** Limma analysis on WGCNA modules *****\n')
    
  # when multiple learner sets are used using balanced sampling
  if(balanceadj == TRUE){
    if(!is.list(dats)){
      cat('The ensemble mode needs multiple learner datasets \nHence, implement balanced sampling.\n')
      return()
    }
    print('**** Computing Module Eigengenes ****')
    MElist <- list()
    # computing module eigengenes for each learner dataset
    for(k in 1:length(dats)){
      MElist[[k]] <- WGCNA::moduleEigengenes(t(dats[[k]]), colors = wgcnares$numericlabels)
      MElist[[k]] <- MElist[[k]]$eigengenes
      names(MElist)[k] <- paste0('round_', k)
    }
    wgcnares$mergedmelist <- MElist
      
    mergedmelist <- list()
    for(i in 1:length(wgcnares$mergedmelist))
      mergedmelist[[i]] <- t(wgcnares$mergedmelist[[i]])
      
    names(mergedmelist) <- names(wgcnares$mergedmelist)
      
    # limma differential analysis (single threaded)
    if(threads == 1){
      print(paste('***** Differential Analysis (Module-level) *****'))
      limmareslist <- list()
      iseqs <- seq(1, nround, 1)
      for(i in iseqs){
        limmares <- diffsites(dats = mergedmelist,
                              pddats = pddats,
                              i = i,
                              responsevarname = responsevarname,
                              responselevels = responselevels,
                              confoundings = confoundings)
        limmareslist[[i]] <- limmares
      }
    }
    else{
      print(paste('***** Differential Analysis (Module-level) *****'))
      diffsites <- function(dats,
                            pddats,
                            i = 1,
                            responsevarname = NULL,
                            responselevels = NULL,
                            confoundings = NULL){

        dat <- dats[[i]]
        pddat <- pddats[[i]]

        if(is.null(responsevarname) & ncol(pddat) >= 2)
          responsevarname <- names(pddat)[2]

        if(is.null(responsevarname))
          return(NULL)

        #library(limma)
        simpddat <- pddat[c('sampleid', responsevarname, confoundings)]
        simpddat <- simpddat[complete.cases(simpddat),]

        if(!is.factor(pddat[[responsevarname]]) & !is.numeric(pddat[[responsevarname]])){

          if(is.null(responselevels)){
            responselevels <- simpddat[[responsevarname]]
            responselevels <- unique(responselevels)
            responselevels <- responselevels[order(responselevels)]
          }
          simpddat$Response <- factor(simpddat[[responsevarname]],
                                        levels = responselevels,
                                        ordered = TRUE)
        }
        else
          simpddat$Response <- simpddat[[responsevarname]]

        vars <- paste(c('Response', confoundings), collapse = ' + ')
        print(paste0('~ ', vars))
        design <- model.matrix(as.formula(paste0('~ ', vars)), data = simpddat)

        betasub <- dat[,simpddat$sampleid, drop = FALSE]
        fit1 <- limma::lmFit(betasub, design)
        fit1 <- limma::eBayes(fit1)
        allg.limma <- limma::topTable(fit1, coef = 2, n = dim(fit1)[1])

        if((nrow(allg.limma) == 1) & (nrow(betasub) == 1))
          row.names(allg.limma) <- row.names(betasub)

        return(allg.limma)
      }
      iseqs <- seq(1, nround, 1)

      #library(doParallel)
      cores <- parallel::detectCores()
      cl <- parallel::makeCluster(min(threads, cores))
      doParallel::registerDoParallel(cl)

      #date()
      `%dopar%` <- foreach::`%dopar%`
      limmareslist <- foreach::foreach(i = iseqs,
                                         #.export = ls(name = globalenv())) %dopar% {
                                         .export = NULL) %dopar% {
                                           diffsites(dats = mergedmelist,
                                                     pddats = pddats,
                                                     i,
                                                     responsevarname = responsevarname,
                                                     responselevels = responselevels,
                                                     confoundings = confoundings)
                                         }
      parallel::stopCluster(cl)
      unregister_dopar()
    }

    names(limmareslist) <- names(dats)
    combineres <- combinepval(limmareslist = limmareslist, method = method)
    metalimmares <- data.frame(logFC = combineres$logfccombs,
                                 P.Value = combineres$pcombs,
                                 adj.P.Val = combineres$padjs,
                                 stringsAsFactors = FALSE)
    limmares <- metalimmares
  }
  else{
    eigengenes <- WGCNA::moduleEigengenes(t(dats), 
                                           colors = wgcnares$numericlabels)
    eigengenes <- eigengenes$eigengenes
    wgcnares$mergedmelist = eigengenes
    mergedmelist <- list(t(eigengenes))

    limmares <- diffsites(dats = mergedmelist,
                          pddats = list(pddats),
                          i = 1,
                          responsevarname = responsevarname,
                          responselevels = responselevels,
                          confoundings = confoundings)
  }
  # plotting limma bar plots
  if(plot == TRUE){
    print('***** Plotting ME bar plot *****')
    if(tolower(pvalcolname) %in% c('pval', 'p-val', 'p.val',
                                     'pvalue', 'p-value', 'p.value'))
      pvalcolname <- 'P.Value'
    else
      pvalcolname <- 'adj.P.Val'

    limmabarplotsangel <- 0

    if(nrow(limmares) > 10)
      limmabarplotsangel <- 45
    plotdat <- limmabarplots(limmares = limmares,
                              diffcutoff = diffcutoff,
                              pvalcut = pvalcutoff,
                              pvalcolname = pvalcolname,
                              confoundings = confoundings,
                              responsevarname = responsevarname,
                              titleprefix = titleprefix,
                              titlesize = titlesize,
                              textsize = textsize,
                              face = face,
                              mesizes = table(wgcnares$numericlabels),
                              xtextangle = limmabarplotsangel)
  }

  sigmes <- limmares[(abs(limmares$logFC) > diffcutoff) &
                         (limmares[,pvalcolname] < pvalcutoff), , drop = FALSE]

  if(nrow(sigmes) > 0)
    sigmes <- row.names(sigmes)
  else
    sigmes <- NULL

  gc()
  res <- list(wgcnares = wgcnares, limmares = limmares)
  return(res)
  
}

diffgenes = function(dats,
                     pddats,
                     wgcnares,
                     
                     topn=FALSE,
                     
                     balanceadj=TRUE,
                     nround = 10,
                     seed = 2022,
                     k = 5, 
                     threads = 1, 
                     method = 'Cauchy',
                     
                     responsevarname = NULL, 
                     responselevels = NULL, 
                     confoundings = NULL,
                     
                     plot=FALSE,
                     featuretype = 'gene',
                     titleprefix='Response',
                     titlesize = 15,
                     textsize = 14,
                     face = 'bold',
                     labelnum = 5,
                     annotextsize = 4,
                     
                     diffcutoff=0,
                     pvalcutoff=0.05,
                     pvalcolname='adj.P.Val',
                     absxcutoff = 0,
                     saveplot=FALSE,
                     savedir=NULL){
  melimmares = NULL
  sigg.featurelimma = NULL

  print('***** limma analysis within WGCNA modules *****\n')
  
  melimmares = list()
  melimmasig = list()
  meplotdat = list()
  
  modules = names(table(wgcnares$numericlabels))
  l = 1
  print(paste('Number of modules:', length(modules)))
  # looping over each significant module 
  for(l in 1:length(modules)){
    print(paste('Module:', modules[l]))
    mename = modules[l]
    #mename = gsub(pattern = 'ME', replacement = '', x = mename)
    mename = as.numeric(mename)
    mefeatures <- names(wgcnares$numericlabels[wgcnares$numericlabels == mename])
    if(balanceadj == TRUE){
      if(!is.list(dats)){
        cat('The ensemble mode needs multiple learner datasets \nHence, implement balanced sampling.\n')
        return()
      }
      featurebalanceddats = dats
      for(m in 1:length(dats))
        featurebalanceddats[[m]] = featurebalanceddats[[m]][mefeatures, , drop = FALSE]
      
      featurelimmareslist = list()
      iseqs = seq(1, nround, 1)
      
      for(i in iseqs){
        featurelimmares = diffsites(dats = featurebalanceddats,
                                    pddats = pddats,
                                    i = i, 
                                    responsevarname = responsevarname, 
                                    responselevels = responselevels, 
                                    confoundings = confoundings)
        featurelimmareslist[[i]] = featurelimmares
      }
      names(featurelimmareslist) = names(dats$dats)
      combineres <- combinepval(limmareslist = featurelimmareslist, method = method)
      
      metafeaturelimmares = data.frame(logFC = combineres$logfccombs,
                                       P.Value = combineres$pcombs,
                                       adj.P.Val = combineres$padjs,
                                       stringsAsFactors = FALSE)
      
      featurelimmares <- metafeaturelimmares
    }
    else{
      featurebalanceddats <- list(dats = list(dats[mefeatures, , drop = FALSE]),
                                  pds = list(pddats))
      featurelimmares <- diffsites(dats = featurebalanceddats$dats,
                                   pddats = featurebalanceddats$pds, 
                                   i = 1,
                                   responsevarname = responsevarname,
                                   responselevels = responselevels,
                                   confoundings = confoundings)
    }
    
    mename = paste('ME', mename, sep='')
    featurelimmares$ME = mename
    
    if('remove' %in% names(featurelimmares))
      sigg.limma <- featurelimmares[(featurelimmares[,pvalcolname] < pvalcutoff) & 
                                      (abs(featurelimmares$logFC) > absxcutoff) & 
                               (featurelimmares$remove == FALSE),]
    else
      sigg.limma <- featurelimmares[(featurelimmares[,pvalcolname] < pvalcutoff) & 
                                      (abs(featurelimmares$logFC) > absxcutoff),]
    
    melimmares[[l]] = featurelimmares
    melimmasig[[l]] = sigg.limma
    
    # plotting feature scatter plots 
    if(plot == TRUE){
      if(saveplot)
        filename = paste(savedir, mename,'_differential.png', sep='')
      else
        filename = NULL
      
      if(tolower(pvalcolname) %in% c('pval', 'p-val', 'p.val', 'pvalue', 'p-value', 'p.value'))
        pvalcolname <- 'P.Value'
      else
        pvalcolname <- 'adj.P.Val'
      
      if(!is.null(titleprefix))
        featureplottitle <- paste0(titleprefix, ' ', mename)
      else
        featureplottitle <- mename
      
      # limma volcano plot
      featureplotdat <- singlevolcano(limmares = featurelimmares,
                                      absxcutoff = absxcutoff,
                                      pvalcut = pvalcutoff,
                                      pvalcolname = pvalcolname,
                                      confoundings = confoundings,
                                      responsevarname = responsevarname,
                                      titleprefix = featureplottitle,
                                      featuretype = featuretype,
                                      titlesize = titlesize,
                                      textsize = textsize,
                                      face = face,
                                      labelnum = labelnum,
                                      annotextsize = annotextsize,
                                      filename = filename)
      
      featureplotdat$ME <- mename
      meplotdat[[l]] <- featureplotdat
    }
  }
  melimmares <- do.call(rbind, melimmares)
  melimmasig <- do.call(rbind, melimmasig)
  
  res = list('melimmares'=melimmares, 'melimmasig'=melimmasig)
  return(res)
}

medgenes = function(dats,
                    pddats,
                    wgcnares,
                    limmares,
                    modules = c(),
                    topn=FALSE,
                     
                    balanceadj=TRUE,
                    nround = 10,
                    seed = 2022,
                    k = 5, 
                    threads = 1, 
                    method = 'Cauchy',
                     
                    responsevarname = NULL, 
                    responselevels = NULL, 
                    confoundings = NULL,
                     
                    plot=FALSE,
                    featuretype = 'gene',
                    titleprefix='Response',
                    titlesize = 15,
                    textsize = 14,
                    face = 'bold',
                    labelnum = 5,
                    labeltype = 'both',
                    annotextsize = 4,
                     
                    diffcutoff=0,
                    pvalcutoff=0.05,
                    pvalcolname='adj.P.Val',
                    absxcutoff = 0,
                    saveplot=FALSE,
                    savedir=NULL,
                    filterres=TRUE){
  melimmares = NULL
  sigg.featurelimma = NULL
  
  print('***** Mediation analysis within WGCNA modules *****\n')

  mediationdat = list()
  plotdat = list()
  mediatorprobedat = list()
  mediatorgenedat = list()
  
  if((is.null(modules)) | (length(modules) == 0))
    modules = as.vector(unique(limmares$ME))
  l = 1
  print(paste('Number of modules:', length(modules)))
  # looping over each significant module 
  for(l in 1:length(modules)){
    mename = modules[l]
    #menum = gsub(pattern = 'ME', replacement = '', x = mename)
    #menum = as.numeric(menum)
    mefeatures = row.names(limmares[limmares$ME == mename,])
    print(paste('***** Module:', mename, '*****', length(mefeatures)))
    probeannores = NULL
    geneannores = NULL
    
    if(balanceadj == TRUE){
      if(!is.list(dats)){
        cat('The ensemble mode needs multiple learner datasets \nHence, implement balanced sampling.\n')
        return()
      }
      for(n in 1:nround){
        print(paste('Round:', n))
        mediatordat <- dats[[n]][mefeatures, , drop = FALSE]
        mediatordat <- t(mediatordat)
        
        otherdat <- pddats[[n]][, c(responsevarname, confoundings), drop = FALSE]
        otherdat <- cbind(otherdat, wgcnares$mergedmelist[[n]][, mename, drop = FALSE])
        row.names(otherdat) <- row.names(mediatordat)

        otherdat = otherdat[complete.cases(otherdat[,responsevarname]),]
        mediatordat = mediatordat[row.names(otherdat),]
        
        mediationreslist <- mediationmod(mediatordat = mediatordat,
                                         predictdat = otherdat,
                                         source = mename,
                                         target = responsevarname,
                                         interaction = FALSE,
                                         bootnum = 100,
                                         IPW = TRUE,
                                         responselevels = responselevels,
                                         threads = threads)
        mediationres <- orgreslist(reslist = mediationreslist,
                                   targetname = responsevarname,
                                   sourcename = mename,
                                   anno = FALSE,
                                   propna = FALSE)

        if(n == 1)
          completeres <- mediationres$completeres
        else
          completeres <- rbind(completeres, mediationres$completeres)
        
        mediationrevreslist <- mediationmod(mediatordat = mediatordat,
                                            predictdat = otherdat,
                                            source = responsevarname,
                                            target = mename,
                                            interaction = FALSE,
                                            bootnum = 100,
                                            IPW = TRUE,
                                            sourcelevels = responselevels,
                                            threads = threads)
        
        mediationrevres <- orgreslist(reslist = mediationrevreslist,
                                      targetname = mename,
                                      sourcename = responsevarname,
                                      anno = FALSE,
                                      propna = FALSE)

        
        if(n == 1)
          completerevres <- mediationrevres$completeres
        else
          completerevres <- rbind(completerevres, mediationrevres$completeres)
      }
      completeres <- combineci(cires = completeres)
      completeres <- mediationanno(reslines = completeres, propna = TRUE)
      completeres = completeres[complete.cases(completeres),]
      if(filterres)
        subres <- subset(completeres, mediation == TRUE)
      else
        subres <- completeres
      
      completerevres <- combineci(cires = completerevres)
      completerevres <- mediationanno(reslines = completerevres, propna = TRUE)
      completerevres = completerevres[complete.cases(completerevres),]
      if(filterres)
        subrevres <- subset(completerevres, mediation == TRUE)
      else
        subrevres <- completerevres
      
      if(nrow(subres) > 0 & nrow(subrevres) > 0 & filterres){
        sharedmediators <- intersect(subres$mediator, subrevres$mediator)
        subres <- subset(subres, !(mediator %in% sharedmediators))
        subrevres <- subset(subrevres, !(mediator %in% sharedmediators))
      }

      qualifies <- rbind(subres, subrevres)
    }
    else{
      mediatordat <- dats[mefeatures, , drop = FALSE]
      mediatordat <- t(mediatordat)
      
      otherdat <- pddats[, c(responsevarname, confoundings), drop = FALSE]
      otherdat <- cbind(otherdat, wgcnares$mergedmelist[, mename, drop = FALSE])
      row.names(otherdat) <- row.names(mediatordat)
      
      # mediation analysis part
      mediationreslist <- mediationmod(mediatordat = mediatordat,
                                       predictdat = otherdat,
                                       source = mename,
                                       target = responsevarname,
                                       interaction = FALSE,
                                       bootnum = 100,
                                       IPW = TRUE,
                                       responselevels = responselevels,
                                       threads = threads)
      
      mediationres <- orgreslist(reslist = mediationreslist,
                                 targetname = responsevarname,
                                 sourcename = mename,
                                 propna = TRUE)
      
      subres <- mediationres$subres
      
      # mediation analysis part
      mediationrevreslist <- mediationmod(mediatordat = mediatordat,
                                          predictdat = otherdat,
                                          source = responsevarname,
                                          target = mename,
                                          interaction = FALSE,
                                          bootnum = 100,
                                          IPW = TRUE,
                                          sourcelevels = responselevels,
                                          threads = threads)
      mediationrevres <- orgreslist(reslist = mediationrevreslist,
                                    targetname = mename,
                                    sourcename = responsevarname,
                                    propna = TRUE)
      
      subrevres <- mediationrevres$subres
      
      if(!is.null(subres) & !is.null(subrevres)){
        sharedmediators <- intersect(subres$mediator, subrevres$mediator)
        subres <- subset(subres, !(mediator %in% sharedmediators))
        subrevres <- subset(subrevres, !(mediator %in% sharedmediators))
      }
      
      qualifies <- rbind(subres, subrevres)
    }
    
    # removing mediators with NA prop values 
    if(filterres)
      qualifies = qualifies[complete.cases(qualifies),]

    mediationdat[[l]] = qualifies
    
    if(plot==TRUE){
      if(saveplot)
        filename = paste(savedir, mename,'_mediation.png', sep='')
      else
        filename = NULL
      
      if(tolower(pvalcolname) %in% c('pval', 'p-val', 'p.val',
                                     'pvalue', 'p-value', 'p.value'))
        pvalcolname = 'P.Value'
      else
        pvalcolname = 'adj.P.Val'
      
      if(!is.null(qualifies)){
        if(nrow(qualifies)>0){
          qualifiesplotdat = limmares[qualifies$mediator, , drop = FALSE]
          othersplot = limmares[setdiff(row.names(limmares), 
                                        row.names(qualifiesplotdat)), , drop = FALSE]
          othersplot$color <- '#C0C0C0'
          qualifiesplotdat$color <- '#C0C0C0'
          qualifiesplotdat$color[qualifies$target == mename] <- '#0000FF'
          qualifiesplotdat$color[qualifies$source == mename] <- '#FF0000'
          qualifiesplotdat <- qualifiesplotdat[rev(seq(1, nrow(qualifiesplotdat))), , drop = FALSE]
          qualifiesplotdat <- rbind(othersplot, qualifiesplotdat)
          qualifiesplotdat$Type <- 'NS'
          qualifiesplotdat$Type[qualifiesplotdat$color == '#0000FF'] <- 'Rev'
          qualifiesplotdat$Type[qualifiesplotdat$color == '#FF0000'] <- 'Fwd'
          
          qualifiesplotdat$label <- FALSE
          
          if(!is.null(labelnum)){
            labelnum = as.numeric(labelnum)
            if(is.na(labelnum))
              labelnum = NULL
          }
          # labeling top mediator genes
          if(!is.null(labelnum)){
            top_fwds = head(qualifies[qualifies$source==mename,'mediator'],labelnum)
            top_revs = head(qualifies[qualifies$target==mename,'mediator'],labelnum)
            print(top_fwds)
            print(top_revs)
  
            if((labeltype == "Both") | (labeltype == 'Rev'))
              qualifiesplotdat[row.names(qualifiesplotdat) %in% top_revs,'label'] = TRUE
            if((labeltype == "Both") | (labeltype == 'Fwd'))
              qualifiesplotdat[row.names(qualifiesplotdat) %in% top_fwds, 'label'] = TRUE
          }
          if(!is.null(titleprefix))
            mediationplottitle <- paste0(titleprefix, ' mediation for ', mename)
          else
            mediationplottitle <- paste0('Mediation for ', mename)
          
          plotdat[[l]] = qualifiesplotdat
          
          # significant genes only 
          mediationplotdat <- singlevolcano(limmares = qualifiesplotdat,
                                            absxcutoff = absxcutoff,
                                            pvalcut = pvalcutoff,
                                            pvalcolname = pvalcolname,
                                            confoundings = confoundings,
                                            responsevarname = responsevarname,
                                            titleprefix = mediationplottitle,
                                            featuretype = featuretype,
                                            titlesize = titlesize,
                                            textsize = textsize,
                                            face = face,
                                            labelnum = labelnum,
                                            annotextsize = annotextsize,
                                            filename = filename)
        }
      }
    }
   
  }
  mediationdat <- do.call(rbind, mediationdat)
  return(list('medres' = mediationdat, 'plotdata'=plotdat))
}