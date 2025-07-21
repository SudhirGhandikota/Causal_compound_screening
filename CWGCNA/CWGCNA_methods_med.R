# https://github.com/yuabrahamliu/CWGCNA
judgevartype <- function(varvals){
  
  if(length(unique(varvals)) < 2)
    stop("The source variable must have at least two unique values.")
  else if(length(unique(varvals)) == 2)
    vartype <- "binary"
  else if(is.factor(varvals) | is.character(varvals))
    vartype <- "multinomial"
  else
    vartype <- "continuous"
  
  return(vartype)
  
}

mediationmod <- function(mediatordat, 
                         predictdat, 
                         source, 
                         target, 
                         interaction = FALSE, 
                         bootnum = 100, 
                         IPW = TRUE, 
                         responselevels = NULL, 
                         sourcelevels = NULL, 
                         threads = 1){
  
  colnames(predictdat)[colnames(predictdat) == source] <- 'source'
  
  colnames(predictdat)[colnames(predictdat) == target] <- 'target'
  
  confoundings <- colnames(predictdat)
  confoundings <- confoundings[-match('source', confoundings)]
  confoundings <- confoundings[-match('target', confoundings)]
  if(length(confoundings) == 0){
    confoundingstr <- NULL
  }else{
    confoundingstr <- paste(confoundings, collapse = '` + `')
    confoundingstr <- paste0('`', confoundings, '`')
    confoundingstr <- paste(confoundingstr, collapse = ' + ')
  }
  
  if(!is.factor(predictdat[['source']]) & !is.numeric(predictdat[['source']])){
    
    if(is.null(sourcelevels)){
      sourcelevels <- predictdat[['source']]
      sourcelevels <- unique(sourcelevels)
      sourcelevels <- sourcelevels[order(sourcelevels)]
    }
    predictdat[['source']] <- factor(predictdat[['source']], 
                                     levels = sourcelevels, 
                                     ordered = TRUE)
  }
  
  mod11str <- paste0('`', target, '` ~ `', source, '` + mediator', 
                     c(' + ', '')[c(!is.null(confoundingstr), is.null(confoundingstr))], 
                     confoundingstr)
  mod1str <- paste0('`', target, '` ~ `', source, '`*mediator', 
                    c(' + ', '')[c(!is.null(confoundingstr), is.null(confoundingstr))], 
                    confoundingstr)
  mod2str <- paste0('mediator ~ `', source, '`', 
                    c(' + ', '')[c(!is.null(confoundingstr), is.null(confoundingstr))], 
                    confoundingstr)
  
  mod11form <- as.formula(mod11str)
  mod1form <- as.formula(mod1str)
  mod2form <- as.formula(mod2str)
  
  #IPW model###
  
  if(IPW == TRUE & !is.null(confoundingstr)){
    
    ipwstr <- paste0('source ~ ', confoundingstr)
    ipwform <- as.formula(ipwstr)
    
    samplenames <- row.names(predictdat)
    `%>%` <- magrittr::`%>%`
    
    vartype <- judgevartype(varvals = predictdat$source)
    
    if(vartype == 'continuous'){
      
      suppressMessages(
        invwt_continuous_weightit <- tryCatch({
          WeightIt::weightit(ipwform, 
                             data = predictdat, 
                             stabilize = TRUE)
        }, warning = function(war){
          NULL
        })
      )
      
      if(is.null(invwt_continuous_weightit)){
        
        invwt_continuous_weightit <- WeightIt::weightit(ipwform, 
                                                        data = predictdat, 
                                                        stabilize = TRUE)
        
        invwt_continuous_weightit <- WeightIt::trim(invwt_continuous_weightit, 
                                                    at = 0.9)
        cat('Extreme weights have been trimed to avoid effective sample size reduction\n')
        
      }
      
      predictdat <- predictdat %>%
        plyr::mutate(invwt = invwt_continuous_weightit$weights)
    }else{
      suppressMessages(
        invwt_discrete_weightit <- tryCatch({
          WeightIt::weightit(ipwform, 
                             data = predictdat, 
                             estimand = "ATE", 
                             method = "ps")
        }, warning = function(war){
          NULL
        })
      )
      
      if(is.null(invwt_discrete_weightit)){
        
        invwt_discrete_weightit <- WeightIt::weightit(ipwform, 
                                                      data = predictdat, 
                                                      estimand = "ATE", 
                                                      method = "ps")
        
        invwt_discrete_weightit <- WeightIt::trim(invwt_discrete_weightit, 
                                                  at = 0.9)
        cat('Extreme weights have been trimed to avoid effective sample size reduction\n')
      }
      
      predictdat <- predictdat %>%
        plyr::mutate(invwt = invwt_discrete_weightit$weights)
    }
    
    row.names(predictdat) <- samplenames
    predictdat <- predictdat[complete.cases(predictdat),]
    predictdat <- predictdat[predictdat[,ncol(predictdat)] >= 0,]
    
  }
  
  colnames(predictdat)[colnames(predictdat) == 'source'] <- source
  colnames(predictdat)[colnames(predictdat) == 'target'] <- target
  
  ##Organize data
  mediatordat <- as.data.frame(mediatordat)
  mediatordat <- mediatordat[row.names(predictdat), , drop = FALSE]
  
  jseqs <- 1:ncol(mediatordat)
  
  if(threads == 1){
    reslist <- list()
    j <- 1
    
    for(j in jseqs){
      singleres <- singlemediation(j = j, 
                                   mediatordat = mediatordat, 
                                   predictdat = predictdat, 
                                   target = target, 
                                   source = source, 
                                   confoundings = confoundings, 
                                   interaction = interaction, 
                                   IPW = IPW, 
                                   responselevels = responselevels, 
                                   mod11form = mod11form, 
                                   mod1form = mod1form, 
                                   mod2form = mod2form, 
                                   bootnum = bootnum)
      reslist[[j]] <- singleres
      
    }
    
  }else{
    
    judgevartype <- function(varvals){
      
      if(length(unique(varvals)) < 2)
        stop("The source variable must have at least two unique values.")
      else if(length(unique(varvals)) == 2)
        vartype <- "binary"
      else if(is.factor(varvals) | is.character(varvals))
        vartype <- "multinomial"
      else
        vartype <- "continuous"
      
      return(vartype)
    }
    
    singlemediation <- function(j, 
                                mediatordat, 
                                predictdat, 
                                target, 
                                source, 
                                confoundings, 
                                interaction = FALSE, 
                                IPW, 
                                responselevels, 
                                mod11form, 
                                mod1form, 
                                mod2form, 
                                bootnum = 100){
      
      mediator <- colnames(mediatordat)[j]
      mediatorval <- mediatordat[, j, drop = FALSE]
      names(mediatorval) <- 'mediator'
      trainingdat <- cbind(predictdat, mediatorval)
      trainingdat <- trainingdat[complete.cases(trainingdat),]
      
      if(IPW == TRUE & length(confoundings) > 0)
        trainingdat <- trainingdat[c(target, source, 'mediator', confoundings, 'invwt')]
      else
        trainingdat <- trainingdat[c(target, source, 'mediator', confoundings)]
      
      if(!is.factor(trainingdat[[target]]) & !is.numeric(trainingdat[[target]])){
        
        if(is.null(responselevels)){
          responselevels <- trainingdat[[target]]
          responselevels <- unique(responselevels)
          responselevels <- responselevels[order(responselevels)]
        }
        trainingdat[[target]] <- factor(trainingdat[[target]], 
                                        levels = responselevels, 
                                        ordered = TRUE)
      }
      
      targettype <- judgevartype(varvals = trainingdat[[target]])
      
      if(targettype == 'binary'){
        #no interaction mod
        mod11 <- glm(mod11form, family = binomial(link = "logit"), data = trainingdat)
        #interaction mod
        mod1 <- glm(mod1form, family = binomial(link = "logit"), data = trainingdat)
        
      }else if(targettype == 'continuous'){
        #no interaction mod
        mod11 <- lm(mod11form, data = trainingdat)
        #interaction mod
        mod1 <- lm(mod1form, data = trainingdat)
        
      }
      
      if(sum(is.na(mod11$coefficients)) > 0)
        return(NULL)
      
      if(sum(is.na(mod1$coefficients)) > 0)
        return(NULL)
      
      #mediator prediction mod
      if(IPW == TRUE & length(confoundings) > 0)
        mod2 <- lm(mod2form, data = trainingdat, weights = invwt)
      else
        mod2 <- lm(mod2form, data = trainingdat)
      
      if(sum(is.na(mod2$coefficients)) > 0)
        return(NULL)
      
      #coef for source (interaction mod)
      theta1 <- mod1$coefficients[2]
      sd.theta1 <- summary(mod1)$coefficients[2, 2]
      pval.theta1 <- summary(mod1)$coefficients[2, 4]
      
      #coef for mediator (interaction mod)
      theta2 <- mod1$coefficients[3]
      sd.theta2 <- summary(mod1)$coefficients[3, 2]
      pval.theta2 <- summary(mod1)$coefficients[3, 4]
      
      #coef for interaction (interaction mod)
      interidx <- length(mod1$coefficients)
      theta3 <- mod1$coefficients[interidx]
      sd.theta3 <- summary(mod1)$coefficients[interidx, 2]
      pval.theta3 <- summary(mod1)$coefficients[interidx, 4]
      
      #coef for source (no interactio mod)
      theta11 <- mod11$coefficients[2]
      sd.theta11 <- summary(mod11)$coefficients[2, 2]
      pval.theta11 <- summary(mod11)$coefficients[2, 4]
      
      #coef for mediator (no interaction mod)
      theta21 <- mod11$coefficients[3]
      sd.theta21 <- summary(mod11)$coefficients[3, 2]
      pval.theta21 <- summary(mod11)$coefficients[3, 4]
      
      #intercept (mediator prediction mod)
      beta0 <- mod2$coefficients[1]
      
      #coef for source (mediator prediction mod)
      beta1 <- mod2$coefficients[2]
      sd.beta1 <- summary(mod2)$coefficients[2, 2]
      pval.beta1 <- summary(mod2)$coefficients[2, 4]
      
      #coef for other confoundings (mediator prediction mod)
      conftailidx <- length(mod2$coefficients)
      if(conftailidx >= 3)
        beta2 <- mod2$coefficients[3:conftailidx]
      else
        beta2 <- NULL
      
      #no interaction mod coefs
      mod11_spon_coef <- as.numeric(c(theta11, sd.theta11, pval.theta11, 
                                      theta21, sd.theta21, pval.theta21, 
                                      0, 0, NA, 
                                      beta1, sd.beta1, pval.beta1))
      
      #interaction mod coefs
      mod1_spon_coef <- as.numeric(c(theta1, sd.theta1, pval.theta1, 
                                     theta2, sd.theta2, pval.theta2, 
                                     theta3, sd.theta3, pval.theta3, 
                                     beta1, sd.beta1, pval.beta1))
      
      names(mod1_spon_coef) <- 
        names(mod11_spon_coef) <- 
        c('theta1', 'sd.theta1', 'pval.theta1', 
          'theta2', 'sd.theta2', 'pval.theta2', 
          'theta3', 'sd.theta3', 'pval.theta3', 
          'beta1', 'sd.beta1', 'pval.beta1')
      
      #(Residual standard error)^2 for mediator prediction model
      sigma2 <- (summary(mod2)$sigma)^2
      
      #Set all source value to 0 (a0) or 1 (a)
      a <- 1
      a0 <- 0
      
      continueconfs <- intersect(colnames(trainingdat), 
                                 gsub(pattern = '`', replacement = '', 
                                      names(beta2)))
      
      discreteconfs <- setdiff(names(beta2), continueconfs)
      
      #Set discrete confoundings to 0 and continuous confoundings to their medians
      if(length(continueconfs) > 0){
        continueconfs <- trainingdat[,continueconfs, drop = FALSE]
        continueconfs <- apply(X = continueconfs, MARGIN = 2, 
                               FUN = function(x){x <- rep(median(x), length(x))})
        row.names(continueconfs) <- row.names(trainingdat)
      }
      else
        continueconfs <- NULL
      
      if(length(discreteconfs) > 0){
        
        discreteconfs <- matrix(data = 0, nrow = nrow(trainingdat), 
                                ncol = length(discreteconfs))
        colnames(discreteconfs) <- setdiff(names(beta2), 
                                           intersect(colnames(trainingdat), 
                                                     gsub(pattern = '`', 
                                                          replacement = '', 
                                                          names(beta2))))
        row.names(discreteconfs) <- row.names(trainingdat)
        
      }
      else
        discreteconfs <- NULL
      
      confs <- cbind(discreteconfs, continueconfs)
      confs <- confs[,names(beta2), drop = FALSE]
      c.all <- confs
      
      #Set mediator value to its median value
      med <- median(trainingdat$mediator)
      
      #For interaction mod, set mod_spon_coef as mod1_spon_coef, 
      #for no interaction mod, set mod_spon_coef as mod11_spon_coef
      
      if(interaction == TRUE)
        mod_spon_coef <- mod1_spon_coef
      else
        mod_spon_coef <- mod11_spon_coef
      
      #Controlled Direct Effect (CDE) (interaction mod, 
      #need to consider interaction term, so theta3 is included)
      cde <- as.numeric((mod_spon_coef[['theta1']] + 
                           mod_spon_coef[['theta3']]*med)*(a - a0))
      
      #Natural Indirect Effect (NIE)
      #The change in the odds of target in assocaition with a defined 
      #change in source (a0 to a), while holding the meidator at the 
      #level it would have naturally been at when source is at a specific 
      #level (a) (interaction mod, need to consider interaction term, 
      #so theta3 is included)
      nie <- as.numeric((mod_spon_coef[['theta2']]*mod_spon_coef[['beta1']] + 
                           mod_spon_coef[['theta3']]*mod_spon_coef[['beta1']]*a)*(a - a0))
      
      #Natural Direct Effect (NDE)
      #The change in the odds of target in association with a defined 
      #change in source (a0 to a), while holding the mediator at the 
      #level it would have naturally been at when source is at the 
      #original level (a0) (interaction mod, need to consider interaction 
      #term, so theta3 is included)
      
      confterm <- tryCatch({
        c.all %*% beta2
      }, error = function(err){
        0
      })
      
      if(targettype == 'continuous'){
        nde <- as.numeric(
          (mod_spon_coef[['theta1']] + 
             mod_spon_coef[['theta3']]*(beta0 + 
                                          mod_spon_coef[['beta1']]*a0 + 
                                          confterm))*(a - a0)
        )
        
        if(unique(nie*nde) > 0)
          pct.med <- nie/(nde + nie)
        else
          pct.med <- abs(nie/nde)
        
      }
      else if(targettype == 'binary'){
        
        nde <- as.numeric(
          (mod_spon_coef[['theta1']] + 
             mod_spon_coef[['theta3']]*(beta0 + 
                                          mod_spon_coef[['beta1']]*a0 + 
                                          confterm + 
                                          mod_spon_coef[['theta2']]*sigma2))*(a - a0) + 
            0.5*mod_spon_coef[['theta3']]^2*sigma2*(a^2 - a0^2)
        )
        
        #pct.med <- exp(nde)*(exp(nie) - 1)/(exp(nde)*exp(nie) - 1)
        if(unique(nie*nde) > 0)
          pct.med <- nie/(nde + nie)
        else
          pct.med <- abs(nie/nde)
      }
      
      ## bootstrapping to estimate the standard errors
      nboot <- bootnum
      theta11.boot <- theta21.boot <- 
        theta1.boot <- theta2.boot <- theta3.boot <- 
        beta0.boot <- beta1.boot <- 
        sigma2.boot <- 
        a0.boot <- a.boot <- med.boot <- cde.boot <- nie.boot <- rep(0, nboot)
      if(!is.null(confs)){
        beta2.boot <- mat.or.vec(nboot, ncol(confs))
        
        if(is.null(dim(beta2.boot)))
          beta2.boot <- matrix(beta2.boot, ncol = 1)
        
      }
      else
        beta2.boot <- NULL
      
      #mat.or.vec create an nr by nc zero matrix if nc is greater than 1, 
      #and a zero vector of legnth nr if nc equals 1.
      nde.boot <- pct.med.boot <- mat.or.vec(nboot, nrow(trainingdat))
      #n is sample number
      
      B <- 1
      for(B in 1:nboot){
        set.seed(B)
        trainingdat.boot <- trainingdat[sample(1:nrow(trainingdat), 
                                               nrow(trainingdat), 
                                               replace = TRUE), ]
        
        if(IPW == TRUE & length(confoundings) > 0){
          fac.boot <- sum(trainingdat.boot$invwt)/sum(trainingdat$invwt)
          trainingdat.boot$invwt <- trainingdat.boot$invwt/fac.boot
          
        }
        
        if(targettype == 'binary'){
          
          #no interaction mod
          mod11.boot <- glm(mod11form, family = binomial(link = "logit"), data = trainingdat.boot)
          #interaction mod
          mod1.boot <- glm(mod1form, family = binomial(link = "logit"), data = trainingdat.boot)
          
        }else if(targettype == 'continuous'){
          
          #no interaction mod
          mod11.boot <- lm(mod11form, data = trainingdat.boot)
          #interaction mod
          mod1.boot <- lm(mod1form, data = trainingdat.boot)
          
        }
        
        #mediator prediction mod
        if(IPW == TRUE & length(confoundings) > 0)
          mod2.boot <- lm(mod2form, data = trainingdat.boot, weights = invwt)
        else
          mod2.boot <- lm(mod2form, data = trainingdat.boot)
        
        #coef for source (interaction mod)
        theta1.boot[B] <- mod1.boot$coefficients[2]
        #coef for mediator (interaction mod)
        theta2.boot[B] <- mod1.boot$coefficients[3]
        #coef for interaction (interaction mod)
        interidx <- length(mod1.boot$coefficients)
        theta3.boot[B] <- mod1$coefficients[interidx]
        
        #coef for source (no interactio mod)
        theta11.boot[B] <- mod11.boot$coefficients[2]
        #coef for mediator (no interaction mod)
        theta21.boot[B] <- mod11.boot$coefficients[3]
        
        #intercept (mediator prediction mod)
        beta0.boot[B] <- mod2.boot$coefficients[1]
        #coef for source (mediator prediction mod)
        beta1.boot[B] <- mod2.boot$coefficients[2]
        
        #coef for other confoundings (mediator prediction mod)
        conftailidx <- length(mod2.boot$coefficients)
        
        if(conftailidx >= 3){
          
          beta2boot <- mod2.boot$coefficients[3:conftailidx]
          remains <- setdiff(names(beta2), names(beta2boot))
          completes <- c(names(beta2boot), remains)
          beta2boot <- c(beta2boot, rep(0, length(remains)))
          names(beta2boot) <- completes
          beta2boot <- beta2boot[names(beta2)]
          beta2.boot[B,] <- beta2boot
          
        }
        
        #(Residual standard error)^2 for mediator prediction model
        sigma2.boot[B] <- (summary(mod2.boot)$sigma)^2
        
        a.boot[B] <- 1
        a0.boot[B] <- 0
        
        #Set discrete confoundings to 0 and continuous confoundings to their medians
        
        if(!is.null(c.all)){
          if(!is.null(continueconfs)){
            continueconfs.boot <- trainingdat.boot[,colnames(continueconfs), 
                                                   drop = FALSE]
            continueconfs.boot <- apply(X = continueconfs.boot, MARGIN = 2, 
                                        FUN = function(x){x <- rep(median(x), length(x))})
            row.names(continueconfs.boot) <- row.names(trainingdat.boot)
          }else
            continueconfs.boot <- NULL
          
          if(!is.null(discreteconfs)){
            discreteconfs.boot <- matrix(data = 0, nrow = nrow(trainingdat.boot), 
                                         ncol = ncol(discreteconfs))
            colnames(discreteconfs.boot) <- colnames(discreteconfs)
            row.names(discreteconfs.boot) <- row.names(trainingdat.boot)
            
          }else
            discreteconfs.boot <- NULL
          
          c.all.boot <- cbind(discreteconfs.boot, continueconfs.boot)
          c.all.boot <- c.all.boot[,colnames(c.all), drop = FALSE]
        }else
          c.all.boot <- NULL
        
        med.boot[B] <- median(trainingdat.boot$mediator)
        
        mod11_spon_coef.boot <- as.numeric(c(theta11.boot[B], 
                                             theta21.boot[B], 
                                             0, 
                                             beta1.boot[B]))
        mod1_spon_coef.boot <- as.numeric(c(theta1.boot[B], 
                                            theta2.boot[B], 
                                            theta3.boot[B], 
                                            beta1.boot[B]))
        
        names(mod11_spon_coef.boot) <- names(mod1_spon_coef.boot) <- 
          c('theta1.boot', 'theta2.boot', 'theta3.boot', 'beta1.boot')
        
        #For interaction mod, set mod_spon_coef as mod1_spon_coef, 
        #for no interaction mod, set mod_spon_coef as mo11_spon_coef
        
        if(interaction == TRUE)
          mod_spon_coef.boot <- mod1_spon_coef.boot
        else
          mod_spon_coef.boot <- mod11_spon_coef.boot
        
        #Controlled Direct Effect (CDE) (interaction mod, need to consider interaction term, 
        #so theta3 is included)
        cde.boot[B] <- as.numeric((mod_spon_coef.boot[['theta1.boot']] + 
                                     mod_spon_coef.boot[['theta3.boot']]*med.boot[B])*
                                    (a.boot[B] - a0.boot[B]))
        
        #Natural Indirect Effect (NIE)
        nie.boot[B] <- as.numeric((mod_spon_coef.boot[['theta2.boot']]*mod_spon_coef.boot[['beta1.boot']] + 
                                     mod_spon_coef.boot[['theta3.boot']]*mod_spon_coef.boot[['beta1.boot']]*a.boot[B])*
                                    (a.boot[B] - a0.boot[B]))
        
        #Natural Direct Effect (NDE)
        if(any(confterm != 0))
          confterm.boot <- c.all.boot %*% beta2.boot[B,]
        else
          confterm.boot <- 0
        
        if(targettype == 'continuous'){
          
          nde.boot[B,] <- as.numeric(
            (mod_spon_coef.boot[['theta1.boot']] + 
               mod_spon_coef.boot[['theta3.boot']]*(beta0.boot[B] + 
                                                      mod_spon_coef.boot[['beta1.boot']]*a0.boot[B] + 
                                                      confterm.boot))*(a.boot[B] - a0.boot[B])
          )
          
          if(unique(nie.boot[B]*nde.boot[B,]) > 0)
            pct.med.boot[B,] <- nie.boot[B]/(nde.boot[B,] + nie.boot[B])
          else
            pct.med.boot[B,] <- abs(nie.boot[B]/nde.boot[B,])
  
        }
        else if(targettype == 'binary'){
          nde.boot[B,] <- as.numeric(
            (mod_spon_coef.boot[['theta1.boot']] + 
               mod_spon_coef.boot[['theta3.boot']]*(beta0.boot[B] + 
                                                      mod_spon_coef.boot[['beta1.boot']]*a0.boot[B] + 
                                                      confterm.boot + 
                                                      mod_spon_coef.boot[['theta2.boot']]*sigma2.boot[B]))*
              (a.boot[B] - a0.boot[B]) + 
              0.5*mod_spon_coef.boot[['theta3.boot']]^2*sigma2.boot[B]*(a.boot[B]^2 - a0.boot[B]^2)
          )
          
          #pct.med.boot[B,] <- exp(nde.boot[B,])*(exp(nie.boot[B]) - 1)/(exp(nde.boot[B,])*exp(nie.boot[B]) - 1)
          if(unique(nie.boot[B]*nde.boot[B,]) > 0)
            pct.med.boot[B,] <- nie.boot[B]/(nde.boot[B,] + nie.boot[B])
          else
            pct.med.boot[B,] <- abs(nie.boot[B]/nde.boot[B,])
        }
        
      }
      
      summarycols <- c(setdiff(colnames(trainingdat)[4:max(4, 
                                                           ncol(trainingdat))], 
                               c('invwt')), 
                       #'CDE', 'CDE_lower', 'CDE_upper', 
                       'NDE', 'NDE_lower', 'NDE_upper', 
                       'NIE', 'NIE_lower', 'NIE_upper', 
                       'Prop', 'Prop_lower', 'Prop_upper')
      summarycols <- summarycols[!is.na(summarycols)]
      summary.tab <- mat.or.vec(nrow(trainingdat), length(summarycols))
      colnames(summary.tab) <- summarycols
      rownames(summary.tab) <- row.names(confs)
      
      if(!is.null(c.all)){
        for(B in 1:ncol(c.all)){
          
          summary.tab[,B] <- c.all[,B]
        }
      }
      
      summary.tab[,'NDE'] <- nde
      summary.tab[,'NDE_lower'] <- apply(nde.boot, 2, function(x){x[order(x)][nboot*0.025]})
      summary.tab[,'NDE_upper'] <- apply(nde.boot, 2, function(x){x[order(x)][nboot*0.975]})
      
      summary.tab[,'NIE'] <- nie
      summary.tab[,'NIE_lower'] <- nie.boot[order(nie.boot)][nboot*0.025]
      summary.tab[,'NIE_upper'] <- nie.boot[order(nie.boot)][nboot*0.975]
      
      summary.tab[,'Prop'] <- pct.med
      
      summary.tab[,'Prop_lower'] <- apply(pct.med.boot, 2, function(x){x[order(x)][nboot*0.025]})
      summary.tab[,'Prop_upper'] <- apply(pct.med.boot, 2, function(x){x[order(x)][nboot*0.975]})
      
      summary.tab <- as.data.frame(summary.tab, stringsAsFactors = FALSE)
      summary.tab$mediator <- mediator
      summary.tab <- unique(summary.tab)
      row.names(summary.tab) <- 1:nrow(summary.tab)
      
      cat(paste0('mediator number = ', j, '\n'))
      
      return(summary.tab)
      
    }
    
    #library(doParallel)
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(min(threads, cores))
    doParallel::registerDoParallel(cl)
    
    #date()
    `%dopar%` <- foreach::`%dopar%`
    reslist <- foreach::foreach(j = jseqs,
                                #.export = ls(name = globalenv())) %dopar% {
                                .export = NULL) %dopar% {
                                  
                                  singlemediation(j = j, 
                                                  mediatordat = mediatordat, 
                                                  predictdat = predictdat, 
                                                  target = target, 
                                                  source = source, 
                                                  confoundings = confoundings, 
                                                  interaction = interaction, 
                                                  IPW = IPW, 
                                                  responselevels = responselevels, 
                                                  mod11form = mod11form, 
                                                  mod1form = mod1form, 
                                                  mod2form = mod2form, 
                                                  bootnum = bootnum)
                                }
    
    parallel::stopCluster(cl)
    
    unregister_dopar()
  }
  
  if(length(reslist) == ncol(mediatordat))
    names(reslist) <- colnames(mediatordat)
  else if(length(reslist) > 0)
    names(reslist) <- colnames(mediatordat)[1:length(reslist)]
  
  return(reslist)
}

mediatorannotation <- function(mediators, 
                               platform = 450){
  
  mediatornames <- c(sub(pattern = '::.*$', replacement = '', mediators))
  mediatornames <- c(mediatornames, c(sub(pattern = '^.*::', replacement = '', mediators)))
  mediatornames <- unique(mediatornames)
  
  probeannores <- probeanno(platform = platform, probes = mediatornames)
  mediatorgeneannos <- getgeneanno(inputs = mediatornames)
  
  if(!is.null(probeannores)){
    
    probeannores <- probeannores[match(mediatornames, 
                                       probeannores$Probe)[!is.na(match(mediatornames, 
                                               probeannores$Probe))], , drop = FALSE]
    
    probeannores <- unique(probeannores)
    row.names(probeannores) <- 1:nrow(probeannores)
    probegeneannos <- getgeneanno(inputs = c(probeannores$UCSC_RefGene_Name, 
                                             probeannores$ENTREZID))
    
  }else
    probegeneannos <- NULL
  
  geneannores <- unique(rbind(mediatorgeneannos, probegeneannos))
  res <- list()
  res$geneannores <- geneannores
  
  if(!is.null(probeannores))
    res$probeannores <- probeannores
  
  return(res)
}

mediationanno <- function(reslines, propna = TRUE){
  
  reslines$mediation <- FALSE
  reslines$complete <- FALSE
  reslines$inconsistent <- FALSE
  
  reslines$mediation[reslines$NIE_lower*reslines$NIE_upper > 0] <- TRUE
  reslines$complete[(reslines$NDE_lower*reslines$NDE_upper < 0) & 
                      reslines$mediation] <- TRUE
  reslines$inconsistent[(reslines$NDE_low*reslines$NDE_upper > 0) & 
                          (reslines$NIE_lower*reslines$NIE_upper > 0) & 
                          (reslines$NDE*reslines$NIE < 0)] <- TRUE
  
  reslines <- reslines[order(-reslines$mediation, 
                             -reslines$complete, 
                             -reslines$inconsistent, 
                             -abs(reslines$Prop)),]
  
  if(propna == TRUE){
    reslines$Prop[(reslines$mediation == FALSE) | 
                    (reslines$mediation == TRUE & reslines$complete == TRUE) | 
                    (reslines$mediation == TRUE & reslines$inconsistent == TRUE)] <- 
      NA
    reslines$Prop_lower[is.na(reslines$Prop)] <- NA
    reslines$Prop_upper[is.na(reslines$Prop)] <- NA
  }
  
  row.names(reslines) <- 1:nrow(reslines)
  return(reslines)
}

orgreslist <- function(reslist, 
                       targetname, 
                       sourcename, 
                       anno = TRUE, 
                       propna = TRUE){
  
  i <- 1
  for(i in 1:length(reslist)){
    
    resline <- reslist[[i]]
    resline$source <- sourcename
    resline$target <- targetname
    
    resline <- resline[,c('source', 'target', 'mediator', 
                          colnames(resline)[1:(ncol(resline) - 3)])]
    
    if(i == 1)
      reslines <- resline
    else
      reslines <- rbind(reslines, resline)
  }
  
  row.names(reslines) <- 1:nrow(reslines)
  
  if(anno == TRUE){
    reslines <- mediationanno(reslines = reslines, propna = propna)
    qualifiedsubs <- subset(reslines, mediation == TRUE)
    
    if(nrow(qualifiedsubs) == 0){
      res <- list(completeres = reslines, 
                  subres = NULL)
    }else{
      row.names(qualifiedsubs) <- 1:nrow(qualifiedsubs)
      res <- list(completeres = reslines, 
                  subres = qualifiedsubs)
    }
    
  }else{
    res <- list(completeres = reslines, 
                subres = NULL)
  }
  return(res)
}

singlemediation <- function(j, 
                            mediatordat, 
                            predictdat, 
                            target, 
                            source, 
                            confoundings, 
                            interaction = FALSE, 
                            IPW, 
                            responselevels, 
                            mod11form, 
                            mod1form, 
                            mod2form, 
                            bootnum = 100){

  mediator <- colnames(mediatordat)[j]
  mediatorval <- mediatordat[, j, drop = FALSE]
  names(mediatorval) <- 'mediator'
  trainingdat <- cbind(predictdat, mediatorval)
  trainingdat <- trainingdat[complete.cases(trainingdat),]
  #print(trainingdat)
  
  if(IPW == TRUE & length(confoundings) > 0)
    trainingdat <- trainingdat[c(target, source, 'mediator', confoundings, 'invwt')]
  else
    trainingdat <- trainingdat[c(target, source, 'mediator', confoundings)]
  
  if(!is.factor(trainingdat[[target]]) & !is.numeric(trainingdat[[target]])){
    
    if(is.null(responselevels)){
      responselevels <- trainingdat[[target]]
      responselevels <- unique(responselevels)
      responselevels <- responselevels[order(responselevels)]
    }
    trainingdat[[target]] <- factor(trainingdat[[target]], 
                                    levels = responselevels, 
                                    ordered = TRUE)
    #trainingdat[[target]] <- as.numeric(trainingdat[[target]]) - 1
  }
  
  targettype <- judgevartype(varvals = trainingdat[[target]])
  if(targettype == 'binary'){
    #no interaction mod
    mod11 <- glm(mod11form, family = binomial(link = "logit"), data = trainingdat)
    #interaction mod
    mod1 <- glm(mod1form, family = binomial(link = "logit"), data = trainingdat)
  }
  else if(targettype == 'continuous'){
    #no interaction mod
    mod11 <- lm(mod11form, data = trainingdat)
    #interaction mod
    mod1 <- lm(mod1form, data = trainingdat)
  }
  
  if(sum(is.na(mod11$coefficients)) > 0)
    return(NULL)
  
  if(sum(is.na(mod1$coefficients)) > 0)
    return(NULL)
  #mediator prediction mod
  if(IPW == TRUE & length(confoundings) > 0)
    mod2 <- lm(mod2form, data = trainingdat, weights = invwt)
  else
    mod2 <- lm(mod2form, data = trainingdat)
  
  if(sum(is.na(mod2$coefficients)) > 0){
    return(NULL)
  }

  #coef for source (interaction mod)
  theta1 <- mod1$coefficients[2]
  sd.theta1 <- summary(mod1)$coefficients[2, 2]
  pval.theta1 <- summary(mod1)$coefficients[2, 4]
  
  #coef for mediator (interaction mod)
  theta2 <- mod1$coefficients[3]
  sd.theta2 <- summary(mod1)$coefficients[3, 2]
  pval.theta2 <- summary(mod1)$coefficients[3, 4]
  
  #coef for interaction (interaction mod)
  interidx <- length(mod1$coefficients)
  theta3 <- mod1$coefficients[interidx]
  sd.theta3 <- summary(mod1)$coefficients[interidx, 2]
  pval.theta3 <- summary(mod1)$coefficients[interidx, 4]
  
  #coef for source (no interactio mod)
  theta11 <- mod11$coefficients[2]
  sd.theta11 <- summary(mod11)$coefficients[2, 2]
  pval.theta11 <- summary(mod11)$coefficients[2, 4]
  
  #coef for mediator (no interaction mod)
  theta21 <- mod11$coefficients[3]
  sd.theta21 <- summary(mod11)$coefficients[3, 2]
  pval.theta21 <- summary(mod11)$coefficients[3, 4]
  
  #intercept (mediator prediction mod)
  beta0 <- mod2$coefficients[1]
  
  #coef for source (mediator prediction mod)
  beta1 <- mod2$coefficients[2]
  sd.beta1 <- summary(mod2)$coefficients[2, 2]
  pval.beta1 <- summary(mod2)$coefficients[2, 4]
  
  #coef for other confoundings (mediator prediction mod)
  conftailidx <- length(mod2$coefficients)
  if(conftailidx >= 3)
    beta2 <- mod2$coefficients[3:conftailidx]
  else
    beta2 <- NULL
  
  #no interaction mod coefs
  mod11_spon_coef <- as.numeric(c(theta11, sd.theta11, pval.theta11, 
                                  theta21, sd.theta21, pval.theta21, 
                                  0, 0, NA, 
                                  beta1, sd.beta1, pval.beta1))
  
  #interaction mod coefs
  mod1_spon_coef <- as.numeric(c(theta1, sd.theta1, pval.theta1, 
                                 theta2, sd.theta2, pval.theta2, 
                                 theta3, sd.theta3, pval.theta3, 
                                 beta1, sd.beta1, pval.beta1))
  
  names(mod1_spon_coef) <- 
    names(mod11_spon_coef) <- 
    c('theta1', 'sd.theta1', 'pval.theta1', 
      'theta2', 'sd.theta2', 'pval.theta2', 
      'theta3', 'sd.theta3', 'pval.theta3', 
      'beta1', 'sd.beta1', 'pval.beta1')
  
  #(Residual standard error)^2 for mediator prediction model
  sigma2 <- (summary(mod2)$sigma)^2
  
  #Set all source value to 0 (a0) or 1 (a)
  a <- 1
  a0 <- 0

  continueconfs <- intersect(colnames(trainingdat), 
                             gsub(pattern = '`', replacement = '', 
                                  names(beta2)))
  
  discreteconfs <- setdiff(names(beta2), continueconfs)
  
  #Set discrete confoundings to 0 and continuous confoundings to their medians
  if(length(continueconfs) > 0){
    continueconfs <- trainingdat[,continueconfs, drop = FALSE]
    continueconfs <- apply(X = continueconfs, MARGIN = 2, 
                           FUN = function(x){x <- rep(median(x), length(x))})
    row.names(continueconfs) <- row.names(trainingdat)
  }else{
    continueconfs <- NULL
  }
  
  if(length(discreteconfs) > 0){
    
    discreteconfs <- matrix(data = 0, nrow = nrow(trainingdat), 
                            ncol = length(discreteconfs))
    colnames(discreteconfs) <- setdiff(names(beta2), 
                                       intersect(colnames(trainingdat), 
                                                 gsub(pattern = '`', 
                                                      replacement = '', 
                                                      names(beta2))))
    row.names(discreteconfs) <- row.names(trainingdat)
    
  }else{
    discreteconfs <- NULL
  }
  
  confs <- cbind(discreteconfs, continueconfs)
  confs <- confs[,names(beta2), drop = FALSE]
  c.all <- confs
  
  #Set mediator value to its median value
  med <- median(trainingdat$mediator)
  
  #For interaction mod, set mod_spon_coef as mod1_spon_coef, 
  #for no interaction mod, set mod_spon_coef as mod11_spon_coef
  
  if(interaction == TRUE)
    mod_spon_coef <- mod1_spon_coef
  else
    mod_spon_coef <- mod11_spon_coef
  
  #Controlled Direct Effect (CDE) (interaction mod, 
  #need to consider interaction term, so theta3 is included)
  cde <- as.numeric((mod_spon_coef[['theta1']] + 
                       mod_spon_coef[['theta3']]*med)*(a - a0))
  
  #Natural Indirect Effect (NIE)
  #The change in the odds of target in assocaition with a defined 
  #change in source (a0 to a), while holding the meidator at the 
  #level it would have naturally been at when source is at a specific 
  #level (a) (interaction mod, need to consider interaction term, 
  #so theta3 is included)
  nie <- as.numeric((mod_spon_coef[['theta2']]*mod_spon_coef[['beta1']] + 
                       mod_spon_coef[['theta3']]*mod_spon_coef[['beta1']]*a)*(a - a0))
  
  #Natural Direct Effect (NDE)
  #The change in the odds of target in association with a defined 
  #change in source (a0 to a), while holding the mediator at the 
  #level it would have naturally been at when source is at the 
  #original level (a0) (interaction mod, need to consider interaction 
  #term, so theta3 is included)
  
  confterm <- tryCatch({
    c.all %*% beta2
  }, error = function(err){
    0
  })
  if(targettype == 'continuous'){
    nde <- as.numeric(
      (mod_spon_coef[['theta1']] + 
         mod_spon_coef[['theta3']]*(beta0 + 
                                      mod_spon_coef[['beta1']]*a0 + 
                                      #mod_spon_coef[['beta1']]*a + 
                                      confterm))*(a - a0)
    )
    
    if(unique(nie*nde) > 0)
      pct.med <- nie/(nde + nie)
    else
      pct.med <- abs(nie/nde)

  }
  else if(targettype == 'binary'){
    
    nde <- as.numeric(
      (mod_spon_coef[['theta1']] + 
         mod_spon_coef[['theta3']]*(beta0 + 
                                      mod_spon_coef[['beta1']]*a0 + 
                                      #mod_spon_coef[['beta1']]*a + 
                                      confterm + 
                                      mod_spon_coef[['theta2']]*sigma2))*(a - a0) + 
        0.5*mod_spon_coef[['theta3']]^2*sigma2*(a^2 - a0^2)
    )
    
    #pct.med <- exp(nde)*(exp(nie) - 1)/(exp(nde)*exp(nie) - 1)
    if(unique(nie*nde) > 0)
      pct.med <- nie/(nde + nie)
    else
      pct.med <- abs(nie/nde)
  }
  
  ## bootstrapping to estimate the standard errors
  nboot <- bootnum
  theta11.boot <- theta21.boot <- 
    theta1.boot <- theta2.boot <- theta3.boot <- 
    beta0.boot <- beta1.boot <- 
    sigma2.boot <- 
    a0.boot <- a.boot <- med.boot <- cde.boot <- nie.boot <- rep(0, nboot)
  if(!is.null(confs)){
    beta2.boot <- mat.or.vec(nboot, ncol(confs))
    
    if(is.null(dim(beta2.boot))){
      beta2.boot <- matrix(beta2.boot, ncol = 1)
    }
    
  }else
    beta2.boot <- NULL
  
  #mat.or.vec create an nr by nc zero matrix if nc is greater than 1, 
  #and a zero vector of legnth nr if nc equals 1.
  nde.boot <- pct.med.boot <- mat.or.vec(nboot, nrow(trainingdat))
  #n is sample number
  
  B <- 1
  for(B in 1:nboot){
    #print(paste('B',B))
    set.seed(B)
    trainingdat.boot <- trainingdat[sample(1:nrow(trainingdat), 
                                           nrow(trainingdat), 
                                           replace = TRUE), ]
    #print(sum(complete.cases(trainingdat.boot)))
    #print(names(trainingdat.boot))
    #print(table(trainingdat.boot$RACE))

    if(IPW == TRUE & length(confoundings) > 0){
      fac.boot <- sum(trainingdat.boot$invwt)/sum(trainingdat$invwt)
      trainingdat.boot$invwt <- trainingdat.boot$invwt/fac.boot
    }
    
    if(targettype == 'binary'){
      #no interaction mod
      mod11.boot <- glm(mod11form, family = binomial(link = "logit"), data = trainingdat.boot)
      #interaction mod
      mod1.boot <- glm(mod1form, family = binomial(link = "logit"), data = trainingdat.boot)
      
    }else if(targettype == 'continuous'){
      #no interaction mod
      mod11.boot <- lm(mod11form, data = trainingdat.boot)
      #interaction mod
      mod1.boot <- lm(mod1form, data = trainingdat.boot)
    }
    
    #mediator prediction mod
    if(IPW == TRUE & length(confoundings) > 0)
      mod2.boot <- lm(mod2form, data = trainingdat.boot, weights = invwt)
    else
      mod2.boot <- lm(mod2form, data = trainingdat.boot)
    
    #coef for source (interaction mod)
    theta1.boot[B] <- mod1.boot$coefficients[2]
    #coef for mediator (interaction mod)
    theta2.boot[B] <- mod1.boot$coefficients[3]
    #coef for interaction (interaction mod)
    interidx <- length(mod1.boot$coefficients)
    theta3.boot[B] <- mod1$coefficients[interidx]
    
    #coef for source (no interactio mod)
    theta11.boot[B] <- mod11.boot$coefficients[2]
    #coef for mediator (no interaction mod)
    theta21.boot[B] <- mod11.boot$coefficients[3]
    
    #intercept (mediator prediction mod)
    beta0.boot[B] <- mod2.boot$coefficients[1]
    #coef for source (mediator prediction mod)
    beta1.boot[B] <- mod2.boot$coefficients[2]
    
    #coef for other confoundings (mediator prediction mod)
    conftailidx <- length(mod2.boot$coefficients)
    
    if(conftailidx >= 3){
      beta2boot <- mod2.boot$coefficients[3:conftailidx]
      remains <- setdiff(names(beta2), names(beta2boot))
      completes <- c(names(beta2boot), remains)
      beta2boot <- c(beta2boot, rep(0, length(remains)))
      names(beta2boot) <- completes
      beta2boot <- beta2boot[names(beta2)]
      #print('names')
      #print(names(beta2boot))
      beta2.boot[B,] <- beta2boot
    }
    
    #(Residual standard error)^2 for mediator prediction model
    sigma2.boot[B] <- (summary(mod2.boot)$sigma)^2
    
    a.boot[B] <- 1
    a0.boot[B] <- 0
    
    #Set discrete confoundings to 0 and continuous confoundings to their medians
    
    if(!is.null(c.all)){
      if(!is.null(continueconfs)){
        continueconfs.boot <- trainingdat.boot[,colnames(continueconfs), 
                                               drop = FALSE]
        continueconfs.boot <- apply(X = continueconfs.boot, MARGIN = 2, 
                                    FUN = function(x){x <- rep(median(x), length(x))})
        row.names(continueconfs.boot) <- row.names(trainingdat.boot)
      }else{
        continueconfs.boot <- NULL
      }
      
      if(!is.null(discreteconfs)){
        discreteconfs.boot <- matrix(data = 0, nrow = nrow(trainingdat.boot), 
                                     ncol = ncol(discreteconfs))
        colnames(discreteconfs.boot) <- colnames(discreteconfs)
        row.names(discreteconfs.boot) <- row.names(trainingdat.boot)
        
      }else{
        discreteconfs.boot <- NULL
      }
      
      c.all.boot <- cbind(discreteconfs.boot, continueconfs.boot)
      c.all.boot <- c.all.boot[,colnames(c.all), drop = FALSE]
    }else{
      c.all.boot <- NULL
    }
    
    med.boot[B] <- median(trainingdat.boot$mediator)
    
    mod11_spon_coef.boot <- as.numeric(c(theta11.boot[B], 
                                         theta21.boot[B], 
                                         0, 
                                         beta1.boot[B]))
    mod1_spon_coef.boot <- as.numeric(c(theta1.boot[B], 
                                        theta2.boot[B], 
                                        theta3.boot[B], 
                                        beta1.boot[B]))
    
    names(mod11_spon_coef.boot) <- names(mod1_spon_coef.boot) <- 
      c('theta1.boot', 'theta2.boot', 'theta3.boot', 'beta1.boot')
    
    #For interaction mod, set mod_spon_coef as mod1_spon_coef, 
    #for no interaction mod, set mod_spon_coef as mo11_spon_coef
    
    if(interaction == TRUE)
      mod_spon_coef.boot <- mod1_spon_coef.boot
    else
      mod_spon_coef.boot <- mod11_spon_coef.boot
    
    #Controlled Direct Effect (CDE) (interaction mod, need to consider interaction term, 
    #so theta3 is included)
    cde.boot[B] <- as.numeric((mod_spon_coef.boot[['theta1.boot']] + 
                                 mod_spon_coef.boot[['theta3.boot']]*med.boot[B])*
                                (a.boot[B] - a0.boot[B]))
    
    #Natural Indirect Effect (NIE)
    nie.boot[B] <- as.numeric((mod_spon_coef.boot[['theta2.boot']]*mod_spon_coef.boot[['beta1.boot']] + 
                                 mod_spon_coef.boot[['theta3.boot']]*mod_spon_coef.boot[['beta1.boot']]*a.boot[B])*
                                (a.boot[B] - a0.boot[B]))

    #Natural Direct Effect (NDE)
    if(any(confterm != 0)){
      confterm.boot <- c.all.boot %*% beta2.boot[B,]
      #print('beta2.boot')
      #print(beta2.boot[B,])
    }
    else
      confterm.boot <- 0
    if(targettype == 'continuous'){
      #print(confterm.boot)
      nde.boot[B,] <- as.numeric(
        (mod_spon_coef.boot[['theta1.boot']] + 
           mod_spon_coef.boot[['theta3.boot']]*(beta0.boot[B] + 
                                                  mod_spon_coef.boot[['beta1.boot']]*a0.boot[B] + 
                                                  confterm.boot))*(a.boot[B] - a0.boot[B])
      )
      if(unique(nie.boot[B]*nde.boot[B,]) > 0)
        pct.med.boot[B,] <- nie.boot[B]/(nde.boot[B,] + nie.boot[B])
      else
        pct.med.boot[B,] <- abs(nie.boot[B]/nde.boot[B,])
      
    }
    else if(targettype == 'binary'){
      nde.boot[B,] <- as.numeric(
        (mod_spon_coef.boot[['theta1.boot']] + 
           mod_spon_coef.boot[['theta3.boot']]*(beta0.boot[B] + 
                                                  mod_spon_coef.boot[['beta1.boot']]*a0.boot[B] + 
                                                  confterm.boot + 
                                                  mod_spon_coef.boot[['theta2.boot']]*sigma2.boot[B]))*
          (a.boot[B] - a0.boot[B]) + 
          0.5*mod_spon_coef.boot[['theta3.boot']]^2*sigma2.boot[B]*(a.boot[B]^2 - a0.boot[B]^2)
      )
      #pct.med.boot[B,] <- exp(nde.boot[B,])*(exp(nie.boot[B]) - 1)/(exp(nde.boot[B,])*exp(nie.boot[B]) - 1)
      if(unique(nie.boot[B]*nde.boot[B,]) > 0)
        pct.med.boot[B,] <- nie.boot[B]/(nde.boot[B,] + nie.boot[B])
      else
        pct.med.boot[B,] <- abs(nie.boot[B]/nde.boot[B,])
    }
  }
  summarycols <- c(setdiff(colnames(trainingdat)[4:max(4, ncol(trainingdat))], 
                           c('invwt')), 
                   #'CDE', 'CDE_lower', 'CDE_upper', 
                   'NDE', 'NDE_lower', 'NDE_upper', 
                   'NIE', 'NIE_lower', 'NIE_upper', 
                   'Prop', 'Prop_lower', 'Prop_upper')
  summarycols <- summarycols[!is.na(summarycols)]
  summary.tab <- mat.or.vec(nrow(trainingdat), length(summarycols))
  colnames(summary.tab) <- summarycols
  rownames(summary.tab) <- row.names(confs)
  
  if(!is.null(c.all)){
    for(B in 1:ncol(c.all))
      summary.tab[,B] <- c.all[,B]
  }

  summary.tab[,'NDE'] <- nde
  summary.tab[,'NDE_lower'] <- apply(nde.boot, 2, function(x){x[order(x)][nboot*0.025]})
  summary.tab[,'NDE_upper'] <- apply(nde.boot, 2, function(x){x[order(x)][nboot*0.975]})
  
  summary.tab[,'NIE'] <- nie
  summary.tab[,'NIE_lower'] <- nie.boot[order(nie.boot)][nboot*0.025]
  summary.tab[,'NIE_upper'] <- nie.boot[order(nie.boot)][nboot*0.975]
  
  summary.tab[,'Prop'] <- pct.med
  summary.tab[,'Prop_lower'] <- apply(pct.med.boot, 2, function(x){x[order(x)][nboot*0.025]})
  summary.tab[,'Prop_upper'] <- apply(pct.med.boot, 2, function(x){x[order(x)][nboot*0.975]})
  
  summary.tab <- as.data.frame(summary.tab, stringsAsFactors = FALSE)
  summary.tab$mediator <- mediator
  summary.tab <- unique(summary.tab)
  row.names(summary.tab) <- 1:nrow(summary.tab)
  
  #cat(paste0('mediator number = ', j, '\n'))
  
  return(summary.tab)
  
}

singleci <- function(sub){
  
  n <- nrow(sub)
  #sub$CDE <- mean(sub$CDE)
  sub$NDE <- mean(sub$NDE)
  sub$NIE <- mean(sub$NIE)
  #sub$TT <- mean(sub$TT)
  sub$Prop <- mean(sub$Prop)
  
  #sub$CDE_lower <- sqrt(mean(sub$CDE_lower^2))
  #sub$CDE_upper <- sqrt(mean(sub$CDE_upper^2))
  sub$NDE_lower <- sqrt(mean(sub$NDE_lower^2))
  sub$NDE_upper <- sqrt(mean(sub$NDE_upper^2))
  sub$NIE_lower <- sqrt(mean(sub$NIE_lower^2))
  sub$NIE_upper <- sqrt(mean(sub$NIE_upper^2))
  #sub$TT_lower <- sqrt(mean(sub$TT_lower^2))
  #sub$TT_upper <- sqrt(mean(sub$TT_upper^2))
  sub$Prop_lower <- sqrt(mean(sub$Prop_lower^2))
  sub$Prop_upper <- sqrt(mean(sub$Prop_upper^2))
  
  #sub$CDE_lower <- sub$CDE - sub$CDE_lower
  #sub$CDE_upper <- sub$CDE + sub$CDE_upper
  sub$NDE_lower <- sub$NDE - sub$NDE_lower
  sub$NDE_upper <- sub$NDE + sub$NDE_upper
  sub$NIE_lower <- sub$NIE - sub$NIE_lower
  sub$NIE_upper <- sub$NIE + sub$NIE_upper
  #sub$TT_lower <- sub$TT - sub$TT_lower
  #sub$TT_upper <- sub$TT + sub$TT_upper
  sub$Prop_lower <- sub$Prop - sub$Prop_lower
  sub$Prop_upper <- sub$Prop + sub$Prop_upper
  
  sub <- sub[, c('source', 'target', 'mediator', 
                 'NDE', 'NDE_lower', 'NDE_upper', 
                 'NIE', 'NIE_lower', 'NIE_upper', 
                 'Prop', 'Prop_lower', 'Prop_upper'), drop = FALSE]
  
  sub <- unique(sub)
  
  return(sub)
}

combineci <- function(cires){
  
  #cires$CDE_lower <- cires$CDE - cires$CDE_lower
  #cires$CDE_upper <- cires$CDE_upper - cires$CDE
  
  cires$NDE_lower <- cires$NDE - cires$NDE_lower
  cires$NDE_upper <- cires$NDE_upper - cires$NDE
  
  cires$NIE_lower <- cires$NIE - cires$NIE_lower
  cires$NIE_upper <- cires$NIE_upper - cires$NIE
  
  #cires$TT_lower <- cires$TT - cires$TT_lower
  #cires$TT_upper <- cires$TT_upper - cires$TT
  
  cires$Prop_lower <- cires$Prop - cires$Prop_lower
  cires$Prop_upper <- cires$Prop_upper - cires$Prop
  
  cires <- cires[order(cires$mediator, cires$source, cires$target),]
  
  combineres <- plyr::ddply(.data = cires, 
                            .variables = c('mediator', 'source', 'target'), 
                            .fun = singleci)
  
  row.names(combineres) <- 1:nrow(combineres)
  
  return(combineres)
  
}