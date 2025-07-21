
library(MASS)
library(car)
library(pROC)
library(MLmetrics)
library(dplyr)
library(splitstackshape)
library(caret)
library(PRROC)
library(plyr)
library(caret)
library(glmnet)
library(splitTools)
library("variancePartition")
library(limma)
library(edgeR)

library(forestplot)
library(ggplot2)

# regularized glmnet models
compute_glmnet_reg = function(expr_df, labels, genes, 
                              clin_variables = c(), 
                              maxit = 20, nfolds=3, train=0.6, test=0.4, 
                              seed=123, plot_dir = plot_dir,
                              plotname=NA){
  # alpha = 0 => ridge regression
  # alpha = 1 => lasso regression
  set.seed(seed)
  X = as.data.frame(expr_df[,genes])
  if(dim(X)[2] == 1)
    # adding intercept
    X = cbind(rep(1, dim(X)[1]), X)
  y = labels
  y = as.numeric(y)

  X = as.data.frame(X[complete.cases(y),])
  y = y[complete.cases(y)]
  if(length(clin_variables)>0)
    X = cbind(X, trait_df[,clin_variables])
  
  if(dim(X)[2]>1)
    X = makeX(X, na.impute = TRUE)
  complete.cases(X)

  # stratified sampling
  inds = partition(y, p=c(train=train, test=test))
  #train_ind = createDataPartition(y, p=0.6, list = FALSE)
  X_train = as.data.frame(X[inds$train,])
  y_train = as.numeric(y[inds$train])
  X_test = as.data.frame(X[inds$test,])
  y_test = as.numeric(y[inds$test])
  
  if(nfolds=='LOOCV')
    nfolds=length(y_train)
  pos_samples = row.names(X_test[y_test==1,])
  neg_samples = row.names(X_test[y_test==0,])
  
  num_train = sum(X_train[,2]!=0)
  num_test = sum(X_test[,2]!=0)
  
  if((num_train<8)| (num_test<8))
    return(list("test_auprc"= NA, 'test_aucroc'=NA))
  
  alphas = seq(0, 1, 0.05)
  accs = c()
  models = list()
  for(i in 1:length(alphas)){
    alpha = alphas[i]
    # regularized regression
    lambdas = 10^seq(2, -3, by = -.1)
    model = cv.glmnet(as.matrix(X_train), y_train, family='binomial', alpha = alpha, 
                      lambda = lambdas, nfolds = nfolds, type.measure = 'class')
    models[[i]] = model
    lambda_best = model$lambda.min
    train_acc = 1-model$cvm[model$lambda==model$lambda.min]
    accs = c(accs, train_acc)
  }
  # best alpha parameter
  best_alpha = alphas[accs==max(accs)][1]
  best_model = models[[which(alpha %in% alphas)]]
  lambda_best = best_model$lambda.min
  
  coefs = predict(best_model, s = lambda_best, type='coef')
  sig_genes = names(coefs[order(as.vector(coefs), decreasing = T),])
  sig_genes = setdiff(sig_genes, "(Intercept)")
  preds = predict(best_model, newx = as.matrix(X_test), s = "lambda.min", type = "response")
  pr = pr.curve(as.vector(preds[pos_samples,]), 
                as.vector(preds[neg_samples,]), curve = T)
  auc_info = as.data.frame(cbind(pr$curve[,1], pr$curve[,2]))
  colnames(auc_info) = c("Recall", "Precision")
  
  if(!is.na(plotname)){
    # Plotting PR Curve
    filename = paste(plot_dir, "/PR_Curves/PR_", plotname, ".png", sep="")
    sizeGrWindow(9,8)
    png(filename = filename)
    # AUC=",round(pr$auc.integral,3)
    p = ggplot(auc_info,aes(x=Recall,y=Precision)) + 
      geom_line(size=0.5) + lims(y = c(0, 1.0)) +
      geom_smooth(se=FALSE, fullrange=FALSE, colour = "dark red") + 
      labs(x="Recall",y="Precision", title='Nonres_vs_Res', subtitle=paste("(AUC = ",round(pr$auc.integral,3), ")", sep="")) + 
      theme_bw(base_size=18) + theme(plot.title = element_text(hjust = 0.5)) + 
      theme(plot.subtitle = element_text(hjust = 0.5)) + 
      theme(legend.position = c(0.3, 0.3)) + 
      theme(plot.title = element_text(size=16, face="bold", hjust = 0.5),
            plot.subtitle = element_text(size=13, hjust = 0.5))
    print(p)
    dev.off()
  }
  
  roc = roc.curve(as.vector(preds[pos_samples,]), 
                  as.vector(preds[neg_samples,]), curve = T)
  auc_info = as.data.frame(cbind(roc$curve[,1], roc$curve[,2]))
  colnames(auc_info) = c("FP", "Recall")
  
  if(!is.na(plotname)){
    # Plotting ROC Curve
    filename = paste(plot_dir, "/ROC_Curves/ROC_", plotname, ".png", sep="")
    sizeGrWindow(9,8)
    png(filename = filename)
    # AUC=",round(pr$auc.integral,3)
    p = ggplot(auc_info,aes(x=FP,y=Recall)) + 
      geom_line(size=0.5) + lims(y = c(0, 1.0)) +
      geom_smooth(se=FALSE, fullrange=FALSE, colour = "dark red") + 
      labs(x="FP",y="Recall", title='Nonres_vs_Res', subtitle=paste("(AUC = ",round(roc$auc,3), ")", sep="")) + 
      theme_bw(base_size=18) + theme(plot.title = element_text(hjust = 0.5)) + 
      theme(plot.subtitle = element_text(hjust = 0.5)) + 
      theme(legend.position = c(0.3, 0.3)) + 
      theme(plot.title = element_text(size=16, face="bold", hjust = 0.5),
            plot.subtitle = element_text(size=13, hjust = 0.5))
    print(p)
    dev.off()
  }
  return(list("model" = best_model, "alpha"=best_alpha, "sig_genes" = sig_genes, "train_acc"=max(accs),
              "test_auprc"= pr$auc.integral, 'test_aucroc'=roc$auc))
}
clean_expr = function(expr_df){
  #removing samples and genes with too many missing values
  expr_df_good = goodSamplesGenes(expr_df,verbose=5)
  sum(!expr_df_good$goodGenes) # excluding 952 genes
  sum(!expr_df_good$goodSamples)
  if(!expr_df_good$allOK){
    #printing genes and samples that were removed
    if(sum(!expr_df_good$goodGenes)>0) #if any genes were "not" considered good
      printFlush(paste("Removing Genes:", paste(names(expr_df)[!expr_df_good$goodGenes], collapse=", ")))
    if(sum(!expr_df_good$goodSamples)>0) #if any samples were "not" considered good
      printFlush(paste("Removing Samples:", paste(names(expr_df)[!expr_df_good$goodSamples], collapse=", ")))
    
    #removing the "bad" samples and genes
    expr_df = expr_df[expr_df_good$goodSamples,expr_df_good$goodGenes]
    return(expr_df)
  }
  rm(expr_df_good)
}

############# GSE213001 #############
expr_df_213001 = read.csv('datasets/expr_data_GSE213001.txt', sep='\t')
expr_df_213001$X.NAME = NULL
row.names(expr_df_213001) = as.vector(expr_df_213001$Symbol)
expr_df_213001$Symbol = NULL

expr_df_213001 = as.data.frame(voom(expr_df_213001)$E)
expr_df_213001 = as.data.frame(t(expr_df_213001))

trait_info = read.csv('datasets/trait_info_GSE213001.txt', sep='\t')
cases = trait_info[trait_info$Diagnosis=="IPF",'Sample']
controls = trait_info[trait_info$Diagnosis=="NDC",'Sample']

base_samples = trait_info[(trait_info$Location=="Base") & 
                            (trait_info$Diagnosis=="IPF"),'Sample']
apex_samples = trait_info[(trait_info$Location=="Apex") & 
                            (trait_info$Diagnosis=="IPF"),'Sample']

severe_samples = trait_info[trait_info$Severity=="Severe",'Sample']
adv_samples = trait_info[trait_info$Severity=="Advanced",'Sample']

mediator_cands = read.csv('CWGCNA_Results/mediator_candidates.txt', sep = '\t')
mediator_cands = as.vector(mediator_cands$Gene.Symbol)
cands = intersect(mediator_cands, names(expr_df_213001))

plot_dir = 'GLMNet_Results/Plots/'

# IPF vs Controls
# new candidates (based on p-values)
labels = c(rep(0, length(controls)), rep(1, length(cases)))

all_results = c()
for(idx in 1:length(cands)){
  symbol = cands[idx]
  test_aurocs = c()
  test_auprcs = c()
  nruns = 100
  num_cases = sum(expr_df_213001[cases,symbol]>0)
  num_ctrls = sum(expr_df_213001[controls,symbol]>0)
  if((num_cases<8) | (num_ctrls<8))
    next
  for(run_idx in 1:nruns){
    seed = sample(1:100000, 1)
    # y_train # y_test
    # 0  1    # 0  1 
    # 21 12   # 14  7
    results_res = compute_glmnet_reg(expr_df_213001[c(controls, cases),], 
                                     labels, seed = seed,
                                     c(symbol), 
                                     maxit=500, nfolds = 'LOOCV', 
                                     train=0.6, test=0.4, 
                                     plot_dir = plot_dir,
                                     plotname = NA)
    test_auprcs = c(test_auprcs, results_res$test_auprc)
    test_aurocs = c(test_aurocs, results_res$test_aucroc)
  }
  print(paste(symbol,' ', num_cases, ' ', num_ctrls, ' AUPRC: ', round(mean(test_auprcs, na.rm=TRUE),3), '+/-', 
              round(sd(test_auprcs, na.rm=TRUE),3), ' AUROC: ', round(mean(test_aurocs, na.rm=TRUE),3), '+/-',
              round(sd(test_aurocs, na.rm = TRUE), 3), sep=''))
  all_results = rbind(all_results, c(symbol, num_cases, num_ctrls, round(mean(test_auprcs, na.rm=TRUE),3), 
                                     round(sd(test_auprcs, na.rm = TRUE), 3),
                                     round(mean(test_aurocs, na.rm = TRUE), 3), 
                                     round(sd(test_aurocs, na.rm = TRUE), 3)))
}
all_results = as.data.frame(all_results)
colnames(all_results) = c('Gene', 'Num_Cases', 'Num_Ctrls', 'AUPRC_avg', 'AUPRC_std', 'AUROC_avg', 'AUROC_std')
# Supplementary File 5
write.table(all_results, 'GLMNet_Results/Model_results.txt', sep='\t')

# Severe vs Controls
length(severe_samples)
labels = c(rep(0, length(controls)), rep(1, length(severe_samples)))

all_results = c()
for(idx in 1:length(cands)){
  symbol = cands[idx]
  test_aurocs = c()
  test_auprcs = c()
  nruns = 100
  num_cases = sum(expr_df_213001[severe_samples,symbol]>0)
  num_ctrls = sum(expr_df_213001[controls,symbol]>0)
  if((num_cases<8) | (num_ctrls<8))
    next
  for(run_idx in 1:nruns){
    seed = sample(1:100000, 1)
    # y_train # y_test
    # 0  1    # 0  1 
    # 21 12   # 14  7
    results_res = compute_glmnet_reg(expr_df_213001[c(controls, severe_samples),], 
                                     labels, seed = seed,
                                     c(symbol), 
                                     maxit=500, nfolds = 'LOOCV', 
                                     train=0.6, test=0.4, 
                                     plot_dir = plot_dir,
                                     plotname = NA)
    test_auprcs = c(test_auprcs, results_res$test_auprc)
    test_aurocs = c(test_aurocs, results_res$test_aucroc)
  }
  print(paste(symbol,' ', num_cases, ' ', num_ctrls, ' AUPRC: ', round(mean(test_auprcs, na.rm=TRUE),3), '+/-', 
              round(sd(test_auprcs, na.rm=TRUE),3), ' AUROC: ', round(mean(test_aurocs, na.rm=TRUE),3), '+/-',
              round(sd(test_aurocs, na.rm = TRUE), 3), sep=''))
  all_results = rbind(all_results, c(symbol, num_cases, num_ctrls, round(mean(test_auprcs, na.rm=TRUE),3), 
                                     round(sd(test_auprcs, na.rm = TRUE), 3),
                                     round(mean(test_aurocs, na.rm = TRUE), 3), 
                                     round(sd(test_aurocs, na.rm = TRUE), 3)))
}
all_results = as.data.frame(all_results)
colnames(all_results) = c('Gene', 'Num_Cases', 'Num_Ctrls', 'AUPRC_avg', 'AUPRC_std', 'AUROC_avg', 'AUROC_std')
# Supplementary File 5
write.table(all_results, 'GLMNet_Results/Model_results_Severe.txt', sep='\t')

# advanced vs Controls
length(adv_samples)
labels = c(rep(0, length(controls)), rep(1, length(adv_samples)))

all_results = c()
for(idx in 1:length(cands)){
  symbol = cands[idx]
  test_aurocs = c()
  test_auprcs = c()
  nruns = 100
  num_cases = sum(expr_df_213001[adv_samples,symbol]>0)
  num_ctrls = sum(expr_df_213001[controls,symbol]>0)
  if((num_cases<8) | (num_ctrls<8))
    next
  for(run_idx in 1:nruns){
    seed = sample(1:100000, 1)
    # y_train # y_test
    # 0  1    # 0  1 
    # 21 12   # 14  7
    results_res = compute_glmnet_reg(expr_df_213001[c(controls, adv_samples),], 
                                     labels, seed = seed,
                                     c(symbol), 
                                     maxit=500, nfolds = 'LOOCV', 
                                     train=0.6, test=0.4, 
                                     plot_dir = plot_dir,
                                     plotname = NA)
    test_auprcs = c(test_auprcs, results_res$test_auprc)
    test_aurocs = c(test_aurocs, results_res$test_aucroc)
  }
  print(paste(symbol,' ', num_cases, ' ', num_ctrls, ' AUPRC: ', round(mean(test_auprcs, na.rm=TRUE),3), '+/-', 
              round(sd(test_auprcs, na.rm=TRUE),3), ' AUROC: ', round(mean(test_aurocs, na.rm=TRUE),3), '+/-',
              round(sd(test_aurocs, na.rm = TRUE), 3), sep=''))
  all_results = rbind(all_results, c(symbol, num_cases, num_ctrls, round(mean(test_auprcs, na.rm=TRUE),3), 
                                     round(sd(test_auprcs, na.rm = TRUE), 3),
                                     round(mean(test_aurocs, na.rm = TRUE), 3), 
                                     round(sd(test_aurocs, na.rm = TRUE), 3)))
}
all_results = as.data.frame(all_results)
colnames(all_results) = c('Gene', 'Num_Cases', 'Num_Ctrls', 'AUPRC_avg', 'AUPRC_std', 'AUROC_avg', 'AUROC_std')
# Supplementary File 5
write.table(all_results, 'GLMNet_Results/Model_results_Adv.txt', sep='\t')

########## Forest Plots ###############
plot_forest2 = function(genes, avgs, stds, xlab='AUPRC', filename=NA, res=72){
  mean = avgs
  stds = stds
  lower = avgs - stds
  lower[lower<0] = 0.1
  upper = avgs + stds
  data <- tibble::tibble(mean  = avgs,
                         lower = lower,
                         upper = upper)
  
  data = as.data.frame(data)
  labeltext = genes
  if(!is.na(filename)){
    if(grepl('tiff', filename))
      tiff(filename, units = "in", width = 8, height = 9, res = res)
    else
      png(filename, units = "in", width = 8, height = 8, res = res)
  }
  xticks = seq(from = 0.85, to = 1.0, by = 0.025)
  ggfp = forestplot(labeltext,
                    data,
                    clip = c(0.3, 1),
                    xlog = TRUE,
                    xlab = xlab,
                    xticks = xticks)
  ggfp = fp_set_style(ggfp,
                      box = "royalblue",
                      line = "darkblue",
                      summary = "royalblue",
                      #hrz_lines = "#999999",
                      txt_gp = fpTxtGp(label = gpar(cex=1.7),
                                       ticks = gpar(cex=1.5),
                                       xlab = gpar(cex=1.5)))
  ggfp = fp_add_header(ggfp, "Gene")
  ggfp = fp_decorate_graph(ggfp, grid = structure(c(0.5, 0.7, 0.9),
                                                  gp = gpar(lty = 2, col = "#CCCCFF")))
  print(ggfp)
  if(!is.na(filename))
    dev.off()
}

model_results = read.csv('GLMNet_Results/Model_results.txt', sep = "\t")
model_results_fil = head(model_results, 10)
model_results_fil$AUPRC_avg
# Figure 4
plot_forest2(model_results_fil$Gene, 
             model_results_fil$AUPRC_avg,
             model_results_fil$AUPRC_std,
             xlab='AUPRC', 
             filename = 'GLMNet_Results/Plots//combined_AUPRC.png', res=300)
plot_forest2(model_results_fil$Gene, 
             model_results_fil$AUROC_avg,
             model_results_fil$AUROC_std,
             xlab='AUROC', 
             filename = 'GLMNet_Results/Plots//combined_AUROC.png', res=300)

model_results = read.csv('GLMNet_Results/Model_results_Severe.txt', sep = "\t")
model_results_fil = head(model_results, 10)
# Figure 4
plot_forest2(model_results_fil$Gene, 
             model_results_fil$AUPRC_avg,
             model_results_fil$AUPRC_std,
             xlab='AUPRC(Severe)', 
             filename = 'GLMNet_Results/Plots//severe_AUPRC.png', res=300)
plot_forest2(model_results_fil$Gene, 
             model_results_fil$AUROC_avg,
             model_results_fil$AUROC_std,
             xlab='AUROC(Severe)', 
             filename = 'GLMNet_Results/Plots//severe_AUROC.png', res=300)

model_results = model_results = read.csv('GLMNet_Results/Model_results_Adv.txt', sep = "\t")
model_results_fil = head(model_results, 10)
# Figure 4
plot_forest2(model_results_fil$Gene, 
             model_results_fil$AUPRC_avg,
             model_results_fil$AUPRC_std,
             xlab='AUPRC(Advanced)', 
             filename = 'GLMNet_Results/Plots//advanced_AUPRC.png', res=300)
plot_forest2(model_results_fil$Gene, 
             model_results_fil$AUROC_avg,
             model_results_fil$AUROC_std,
             xlab='AUROC(Advanced)', 
             filename = 'GLMNet_Results/Plots//advanced_AUROC.png', res=300)