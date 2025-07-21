#dat is a matrix using features as rows and samples as columns
#pds data.frame must contain a column named "sampleid"
# https://github.com/yuabrahamliu/CWGCNA
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
  design <- model.matrix(as.formula(paste0('~ ', vars)), 
                         data = simpddat)
  betasub <- dat[,simpddat$sampleid, drop = FALSE]

  fit1 <- limma::lmFit(betasub, design)
  fit1 <- limma::eBayes(fit1)
  
  allg.limma <- limma::topTable(fit1, coef = 2, n = dim(fit1)[1])
  
  if((nrow(allg.limma) == 1) & (nrow(betasub) == 1))
    row.names(allg.limma) <- row.names(betasub)
  
  return(allg.limma)
  
}

combinepval <- function(limmareslist, 
                        method = 'Cauchy'){
  
  probes <- row.names(limmareslist[[1]])
  pvallist <- list()
  logfclist <- list()
  i <- 1
  for(i in 1:length(limmareslist)){
    
    limmareslist[[i]] <- limmareslist[[i]][probes,]
    
    pvals <- limmareslist[[i]]$P.Value
    logfcs <- limmareslist[[i]]$logFC
    
    names(pvals) <- probes
    names(logfcs) <- probes
    
    pvallist[[i]] <- pvals
    logfclist[[i]] <- logfcs
    
  }
  
  names(pvallist) <- names(logfclist) <- names(limmareslist)
  pvaldat <- do.call(rbind, pvallist)
  logfcdat <- do.call(rbind, logfclist)
  
  if(method == 'Fisher'){
    
    pcombs <- apply(X = pvaldat, MARGIN = 2, FUN = function(x){metap::sumlog(x)$p})
    padjs <- p.adjust(pcombs, method = 'BH')
    
    
  }
  else if(method == 'harmonicmean'){
    
    pcombs <- apply(X = pvaldat, MARGIN = 2, 
                    FUN = harmonicmeanp::p.hmp, 
                    w = rep(1/nrow(pvaldat), nrow(pvaldat)), 
                    L = nrow(pvaldat))
    padjs <- p.adjust(pcombs, method = 'BH')
    
  }
  else if(method == 'Cauchy'){
    
    #pcombs <- apply(X = pvaldat, MARGIN = 2, FUN = ACAT::ACAT)
    pcombs <- apply(X = pvaldat, MARGIN = 2, FUN = ACAT)
    padjs <- p.adjust(pcombs, method = 'BH')
    
  }
  
  logfcsems <- apply(X = logfcdat, MARGIN = 1, FUN = function(x){sqrt(var(x)/length(x))})
  ivws <- 1/logfcsems/sum(1/logfcsems)
  logfccombs <- (t(logfcdat) %*% ivws)[,1]
  
  res <- list(pcombs = pcombs, 
              padjs = padjs, 
              logfccombs = logfccombs)
  
  return(res)
  
}

ACAT<-function(Pvals,weights=NULL,is.check=TRUE){
  Pvals<-as.matrix(Pvals)
  if (is.check){
    #### check if there is NA
    if (sum(is.na(Pvals))>0){
      stop("Cannot have NAs in the p-values!")
    }
    #### check if Pvals are between 0 and 1
    if ((sum(Pvals<0)+sum(Pvals>1))>0){
      stop("P-values must be between 0 and 1!")
    }
    #### check if there are pvals that are either exactly 0 or 1.
    is.zero<-(colSums(Pvals==0)>=1)
    is.one<-(colSums(Pvals==1)>=1)
    if (sum((is.zero+is.one)==2)>0)
      stop("Cannot have both 0 and 1 p-values in the same column!")
    
    if (sum(is.zero)>0)
      warning("There are p-values that are exactly 0!")
    
    if (sum(is.one)>0)
      warning("There are p-values that are exactly 1!")
    
  }
  #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
  if (is.null(weights))
    is.weights.null<-TRUE
  else{
    is.weights.null<-FALSE
    weights<-as.matrix(weights)
    if (sum(dim(weights)!=dim(Pvals))>0){
      stop("The dimensions of weights and Pvals must be the same!")
    }else if (is.check & (sum(weights<0)>0)){
      stop("All the weights must be nonnegative!")
    }else{
      w.sum<-colSums(weights)
      if (sum(w.sum<=0)>0){
        stop("At least one weight should be positive in each column!")
      }else{
        for (j in 1:ncol(weights)){
          weights[,j]<-weights[,j]/w.sum[j]
        }
      }
    }
    
  }
  
  #### check if there are very small non-zero p values and calcuate the cauchy statistics
  is.small<-(Pvals<1e-15)
  if (is.weights.null){
    Pvals[!is.small]<-tan((0.5-Pvals[!is.small])*pi)
    Pvals[is.small]<-1/Pvals[is.small]/pi
    cct.stat<-colMeans(Pvals)
  }else{
    Pvals[!is.small]<-weights[!is.small]*tan((0.5-Pvals[!is.small])*pi)
    Pvals[is.small]<-(weights[is.small]/Pvals[is.small])/pi
    cct.stat<-colSums(Pvals)
  }
  #### return the ACAT p value(s).
  pval<-pcauchy(cct.stat,lower.tail = F)
  return(pval)
}

# limma bar plot of modules
limmabarplots <- function(limmares, 
                          diffcutoff = 0, 
                          pvalcut = 0.05, 
                          pvalcolname = 'adj.P.Val', 
                          confoundings = NULL, 
                          responsevarname = NULL, 
                          titleprefix = NULL, 
                          titlesize = 15, 
                          textsize = 13, 
                          face = 'bold', 
                          mesizes = NULL, 
                          xtext = 'log2FC',
                          xtextangle = 0){
  
  facetname = xtext
  
  vars <- paste(c(responsevarname, confoundings), collapse = ' + ')
  sigg.limma <- limmares[(limmares[,pvalcolname] < pvalcut) & 
                           (abs(limmares$logFC) > diffcutoff),]
  nonsigg.limma <- limmares[(limmares[,pvalcolname] >= pvalcut) | 
                              (abs(limmares$logFC) <= diffcutoff),]
  
  if(nrow(sigg.limma) == 0){
    
    both_pos <- limmares[c('logFC', pvalcolname)]
    both_pos$Type <- 'NonSig'
    both_pos$ME <- row.names(both_pos)
    
    row.names(both_pos) <- 1:nrow(both_pos)
    
  }
  else{
    #Draw probe valcano
    both_pos <- sigg.limma[c('logFC', pvalcolname)]
    
    both_pos <- as.data.frame(both_pos, stringsAsFactors = FALSE)
    
    both_pos$Type <- 'NonSig'
    both_pos$Type[both_pos$logFC > diffcutoff] <- 'UP'
    both_pos$Type[both_pos$logFC < -diffcutoff] <- 'DN'
    
    hypernum <- nrow(subset(both_pos, Type == 'UP'))
    hyponum <- nrow(subset(both_pos, Type == 'DN'))
    
    both_pos_nonsig <- subset(both_pos, Type == 'NonSig')
    both_pos <- subset(both_pos, Type != 'NonSig')
    
    both_nonsig <- nonsigg.limma[c('logFC', pvalcolname)]
    if(dim(both_nonsig)[1]>0)
      both_nonsig$Type <- 'NonSig'
    both_nonsig <- rbind(both_nonsig, both_pos_nonsig)
    
    nonsignum <- nrow(both_nonsig)
    
    both_pos <- rbind(both_nonsig, both_pos)
    both_pos$ME <- row.names(both_pos)
    ncnum <- nrow(both_pos) - hypernum - hyponum
    row.names(both_pos) <- 1:nrow(both_pos)
  }
  
  names(both_pos)[names(both_pos) == pvalcolname] <- 'pval'
  if(pvalcolname == 'P.Value')
    pvalfacetname <- 'p-value'
  else
    pvalfacetname <- 'adjusted p-value'
  
  if(is.null(titleprefix))
    maintitle <- paste0('Differential WGCNA modules (', vars, ')')
  else
    maintitle <- paste0(titleprefix, ' differential WGCNA modules (', vars, ')')
  
  subtitle <- NULL
  
  if(!is.null(mesizes)){
    if(length(mesizes) <= 5){
      
      subtitle <- paste0('ME', seq(0, length(mesizes) - 1, 1), 
                         ' = ', as.vector(mesizes))
      subtitle <- paste0('(', paste(subtitle, collapse = '; '), ')')
      
    }
  }
  
  barplotdat <- reshape::melt(both_pos, c('Type', 'ME'))
  barplotdat$variable <- as.character(barplotdat$variable)
  barplotdat$variable[barplotdat$variable == 'pval'] <- pvalfacetname

  mes <- unique(as.numeric(gsub(pattern = 'ME', replacement = '', x = barplotdat$ME)))
  mes <- mes[order(mes)]
  mes <- paste0('ME', mes)
  barplotdat$ME <- factor(x = barplotdat$ME, levels = mes, ordered = TRUE)
  barplotdat$value[barplotdat$variable == pvalfacetname] <- 
    -log10(barplotdat$value[barplotdat$variable == pvalfacetname])
  barplotdat$variable <- factor(x = barplotdat$variable, 
                                levels = c(pvalfacetname, 'logFC'), 
                                ordered = TRUE)
  
  data_hline <- data.frame(variable = c(setdiff(barplotdat$variable, 'logFC'), 
                                        rep('logFC', 3)), 
                           hline = NA, 
                           stringsAsFactors = FALSE)
  data_hline$hline[data_hline$variable == pvalfacetname] <- -log10(pvalcut)
  data_hline$hline[data_hline$variable == 'logFC'] <- c(diffcutoff, -diffcutoff, 0)
  data_hline$color <- c('red', 'blue', 'blue', 'gray')
  
  facetnames <- list(pvalfacetname = paste0('-log10(', pvalfacetname, ')'), 
                     logFC = facetname)
  names(facetnames) <- c(pvalfacetname, 'logFC')
  
  facetlabeller <- function(variable, value)
    return(facetnames[value])
  
  p <- ggplot2::ggplot(barplotdat, ggplot2::aes(x = ME, y = value))
  
  print(
    p + ggplot2::geom_bar(ggplot2::aes(fill = ME),
                          stat = 'identity') + 
      ggplot2::xlab(label = '') +
      ggplot2::ylab(label = 'Value') +
      ggplot2::ggtitle(maintitle, subtitle = subtitle) + 
      ggplot2::scale_fill_manual(values = c('#C0C0C0', 
                                            scales::hue_pal()(length(unique(barplotdat$ME)) - 1))) + 
      ggplot2::facet_wrap(ggplot2::vars(variable), scales = 'free', ncol = 2, 
                          labeller = facetlabeller) + 
      
      ggplot2::geom_hline(data = data_hline, ggplot2::aes(yintercept = hline), 
                          color = data_hline$color) + 
      ggplot2::theme_bw() +
      
      ggplot2::theme(axis.text.y = ggplot2::element_text(face = face, size = textsize),
                     axis.title.y = ggplot2::element_text(face = face, size = textsize), 
                     axis.text.x = ggplot2::element_text(face = face, size = textsize, 
                                                         angle = xtextangle),
                     plot.title = ggplot2::element_text(size = titlesize, face = face),
                     plot.subtitle = ggplot2::element_text(size = textsize, face = face),
                     strip.text = ggplot2::element_text(size = textsize - 1, face = face), 
                     legend.position = 'none')
  )
  
  row.names(both_pos) <- both_pos$ME
  both_pos <- both_pos[,c('pval', 'logFC', 'Type')]
  names(both_pos)[1] <- pvalcolname
  
  return(both_pos)
  
}

#Get volcano plot for limma comparison result
singlevolcano <- function(limmares, 
                          absxcutoff = 0, 
                          pvalcut = 0.05, 
                          pvalcolname = 'adj.P.Val', 
                          confoundings = NULL, 
                          responsevarname = 'Response', 
                          titleprefix = NULL, 
                          featuretype = 'gene', 
                          titlesize = 15, 
                          textsize = 13, 
                          face = 'bold', 
                          labelnum = 5, 
                          annotextsize = 4, 
                          
                          responsecause = NULL, 
                          responseres = NULL,
                          
                          annotatetext=TRUE,
                          filename = NULL,
                          width=600,
                          height=800,
                          res=72){
  
  if(absxcutoff < 0)
    absxcutoff = -absxcutoff
  print(names(limmares))
  
  xaxis <- 'log2FC'
  vars <- paste(c(responsevarname, confoundings), collapse = ' + ')
  
  if('remove' %in% names(limmares)){
    sigg.limma <- limmares[(limmares[,pvalcolname] < pvalcut) & (abs(limmares$logFC) > absxcutoff) & 
                             (limmares$remove == FALSE),]
    nonsigg.limma <- limmares[(limmares[,pvalcolname] >= pvalcut) | (abs(limmares$logFC) <= absxcutoff) | 
                                (limmares$remove == TRUE),]
  }
  else{
    sigg.limma <- limmares[(limmares[,pvalcolname] < pvalcut) & (abs(limmares$logFC) > absxcutoff),]
    nonsigg.limma <- limmares[(limmares[,pvalcolname] >= pvalcut) | (abs(limmares$logFC) <= absxcutoff),]
  }
  
  if(('color' %in% names(limmares)) & ('Type' %in% names(limmares))){
    both_pos <- limmares
    names(both_pos)[names(both_pos) == 'color'] <- 'denscolor'
    both_pos$probe <- row.names(both_pos)
    types <- unique(both_pos$Type)
    subtitlename <- paste0(types, ' = ', as.vector(table(limmares$Type)[types]))
    subtitlename <- rev(subtitlename)
    subtitlename <- paste(subtitlename, collapse = ', ')
    subtitlename <- paste0('(', subtitlename, ')')
    
    both_pos$Sample_Group <- subtitlename
    row.names(both_pos) <- 1:nrow(both_pos)
    
    legendmapping <- both_pos[c('Type', 'denscolor')]
    legendmapping <- unique(legendmapping)
    legendmapping <- legendmapping[match(rev(types), legendmapping$Type),]
    row.names(legendmapping) <- 1:nrow(legendmapping)
  }
  else{
    if(nrow(sigg.limma) == 0){
      
      both_pos <- limmares[c('logFC', pvalcolname)]
      both_pos$Type <- 'NonSig'
      both_pos$denscolor <- '#C0C0C0'
      both_pos$red <- 192
      both_pos$green <- 192
      both_pos$blue <- 192
      both_pos$probe <- row.names(both_pos)
      
      subtitlename <- paste0('(NC = ', nrow(both_pos), ')')
      subtitlename = ''
      both_pos$Sample_Group <- subtitlename
      row.names(both_pos) <- 1:nrow(both_pos)
      
    }
    else{
      #Draw probe valcano
      both_pos <- sigg.limma[c('logFC', pvalcolname)]
      both_pos <- as.data.frame(both_pos, stringsAsFactors = FALSE)
      
      both_pos$Type <- 'NonSig'
      both_pos$Type[both_pos$logFC > absxcutoff] <- 'UP'
      both_pos$Type[both_pos$logFC < -absxcutoff] <- 'DN'
      
      hypernum <- nrow(subset(both_pos, Type == 'UP'))
      hyponum <- nrow(subset(both_pos, Type == 'DN'))
      
      both_pos_nonsig <- subset(both_pos, Type == 'NonSig')
      both_pos <- subset(both_pos, Type != 'NonSig')
      
      if(nrow(both_pos) > 50)
        myColor <- densCols(both_pos$logFC, -log10(both_pos[,pvalcolname]), 
                            colramp = colorRampPalette(rev(rainbow(10, end=4/6))))
      else
        myColor <- rep('blue', nrow(both_pos))
      
      both_pos$denscolor <- myColor
      rgbmat <- t(col2rgb(myColor))
      rgbmat <- as.data.frame(rgbmat, stringsAsFactors = FALSE)
      both_pos <- cbind(both_pos, rgbmat)
      both_pos <- both_pos[order(-both_pos$blue, both_pos$red, both_pos$green),]
      both_pos1 <- subset(both_pos, blue >= red)
      both_pos2 <- subset(both_pos, blue < red)
      both_pos2 <- both_pos2[order(-both_pos2$blue, both_pos2$red, -both_pos2$green),]
      both_pos <- rbind(both_pos1, both_pos2)
      both_nonsig <- nonsigg.limma[c('logFC', pvalcolname)]
      
      if(nrow(both_nonsig) > 0){
        both_nonsig$Type <- 'NonSig'
        both_nonsig <- rbind(both_nonsig, both_pos_nonsig)
        both_nonsig$denscolor <- '#C0C0C0'
        both_nonsig$red <- 192
        both_nonsig$green <- 192
        both_nonsig$blue <- 192
        both_pos <- rbind(both_nonsig, both_pos)
      }
      
      nonsignum <- nrow(both_nonsig)
      both_pos$probe <- row.names(both_pos)
      ncnum <- nrow(both_pos) - hypernum - hyponum
      
      # pieces <- c(paste0('UP = ', hypernum), paste0('DN = ', hyponum), 
      #             paste0('NC = ', ncnum))
      pieces <- c(paste0('UP = ', hypernum), paste0('DN = ', hyponum))
      pieces <- pieces[c(hypernum, hyponum, ncnum) != 0]
      subtitlename <- paste(pieces, collapse = ', ')
      subtitlename <- paste0('(', subtitlename, ')')
      
      both_pos$Sample_Group <- subtitlename
      row.names(both_pos) <- 1:nrow(both_pos)
    }
  }
  names(both_pos)[names(both_pos) == pvalcolname] <- 'pval'

  subtitlename = ''
  
  if(pvalcolname == 'P.Value')
    pvalaxisname <- 'p-value'
  else
    pvalaxisname <- 'adjusted p-value'
  
  maintitle = titleprefix
  if(!('label' %in% names(both_pos))){
    
    both_pos$label <- FALSE
    if(!is.null(labelnum)){
      labelnum <- as.numeric(labelnum)
      if(is.na(labelnum))
        labelnum <- NULL
    }
    
    if(!is.null(labelnum)){
      
      subup <- subset(both_pos, Type == 'UP')
      subdn <- subset(both_pos, Type == 'DN')
      if(nrow(subup) > 0){
        labelupnum <- min(labelnum, nrow(subup))
        labelupprobes <- subup[order(subup$pval, -abs(subup$logFC)),]$probe[seq(1, labelupnum, 1)]
        
      }
      else{
        labelupnum <- 0
        labelupprobes <- NULL
      }
      
      if(nrow(subdn) > 0){
        labeldnnum <- min(labelnum, nrow(subdn))
        labeldnprobes <- subdn[order(subdn$pval, -abs(subdn$logFC)),]$probe[seq(1, labeldnnum, 1)]
      }
      else{
        labeldnnum <- 0
        labeldnprobes <- NULL
      }
      both_pos$label[both_pos$probe %in% c(labelupprobes, labeldnprobes)] <- TRUE
    }
  }
  
  if(!is.null(filename)){
    if(grepl('.tiff', filename))
      tiff(filename, res = res, width = width, height = height, units = 'px')
    else
      png(filename = filename)
  }

  if(('color' %in% names(limmares)) & ('Type' %in% names(limmares))){
    p = ggplot2::ggplot(both_pos, ggplot2::aes(x = logFC, y = -log10(pval), 
                                                color = Type))
      
    p = p + ggplot2::geom_point(position='jitter') + 
        ggplot2::xlab(xaxis) + ggplot2::ylab(paste0('-log10(', pvalaxisname, ')')) + 
        ggplot2::ggtitle(maintitle, 
                         subtitle = subtitlename) + 
        ggplot2::geom_hline(yintercept = -log10(pvalcut), color = 'red') + 
        ggplot2::geom_vline(xintercept = absxcutoff, color = 'blue') + 
        ggplot2::geom_vline(xintercept = -absxcutoff, color = 'blue') + 
        ggplot2::geom_hline(yintercept = 0, color = 'gray', size = 0.5) + 
        ggplot2::geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
        ggplot2::theme_bw() + 
        ggplot2::scale_color_manual(name = 'Type',
                                    breaks = legendmapping$Type,
                                    values = legendmapping$denscolor, 
                                    guide = 'legend') + 
        ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                       plot.title = ggplot2::element_text(size = titlesize, face = face, hjust=0.5),
                       plot.subtitle = ggplot2::element_text(size = textsize, face = face), 
                       axis.text = ggplot2::element_text(size = textsize, face = face), 
                       axis.title = ggplot2::element_text(size = textsize, face = face), 
                       plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm"), 
                       legend.title = ggplot2::element_text(size = textsize, face = face), 
                       legend.text = ggplot2::element_text(size = textsize, face = face)) + 
        ggplot2::theme(legend.position = c(0.15 ,0.65))
        if(annotatetext)
          p = p + ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(label == TRUE, as.character(probe),'')), 
                                           size = annotextsize, fontface = face, max.overlaps = 1000, color = 'black')
    print(p)
    if(!is.null(filename))
      dev.off()
    
  }
  else{
    p = ggplot2::ggplot(both_pos, ggplot2::aes(x = logFC, y = -log10(pval)))
    
    p = p + ggplot2::geom_point(color = both_pos$denscolor, position='jitter') + 
        ggplot2::xlab(xaxis) + ggplot2::ylab(paste0('-log10(', pvalaxisname, ')')) + 
        ggplot2::ggtitle(maintitle, 
                         subtitle = subtitlename) + 
        ggplot2::geom_hline(yintercept = -log10(pvalcut), color = 'red') + 
        ggplot2::geom_vline(xintercept = absxcutoff, color = 'blue') + 
        ggplot2::geom_vline(xintercept = -absxcutoff, color = 'blue') + 
        ggplot2::geom_hline(yintercept = 0, color = 'gray', size = 0.5) + 
        ggplot2::geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
        ggplot2::theme_bw() + 
        ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                       plot.title = ggplot2::element_text(size = titlesize, face = face),
                       plot.subtitle = ggplot2::element_text(size = textsize, face = face), 
                       axis.text = ggplot2::element_text(size = textsize, face = face), 
                       axis.title = ggplot2::element_text(size = textsize, face = face), 
                       plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm")) + 
        
        ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(label == TRUE, as.character(probe),'')), 
                                 size = annotextsize, fontface = face, max.overlaps = 1000)
    print(p)
    if(!is.null(filename))
      dev.off()
  }
  
  row.names(both_pos) <- both_pos$probe
  both_pos <- both_pos[,c('pval', 'logFC', 'Type', 'denscolor', 'label')]
  names(both_pos)[1] <- pvalcolname
  names(both_pos)[4] <- 'color'
  names(both_pos)[5] <- 'labeled'
  return(both_pos)
}