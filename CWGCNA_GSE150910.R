#library(CWGCNA)
setwd('/home/aniljegga2/RScripts/CWGCNA/')
source('CWGCNA_methods_new.R')
source('CWGCNA_methods_diff.R')
library(DESeq2)
library(limma)

map_direction = function(x){
  if(grepl('ME', x[['source']]))
    return('Forward')
  else
    return('Reverse')
}
get_membership = function(x){
  return(signif(as.numeric(gene_membership[ x[['mediator']][1], x[['Module_num']][1] ]),2))
}
get_membership_pval = function(x){
  return(signif(as.numeric(mm_pval[ x[['mediator']][1], x[['Module_num']][1] ]),3))
}

expr_df = read.csv('datasets/expr_data.txt', sep='\t')
dim(expr_df)

row.names(expr_df) = as.vector(expr_df$Gene_Symbol)
expr_df$Gene_Symbol = NULL

# normalizing raw RNAseq counts
expr_df = as.data.frame(voom(expr_df)$E)

trait_df = read.csv('datasets/trait_data.txt', sep='\t')
row.names(trait_df) = as.vector(trait_df$sample)

trait_df = trait_df[names(expr_df),]
trait_df$sampleid = trait_df$sample
selected_traits = c('sampleid', 'sex', 'age', 'Diagnosis', 'smoking', 'FEV', 'FVC', 'DLCO')
trait_df_filtered = trait_df[,selected_traits]
print('done!')

############# With WGCNA ####################
# run WGCNA
wgcnares = runWGCNA(expr_df, 
                    filterdat=FALSE,
                    filternum = 5000,
                    varicancetype = 'sd',
                    powers = seq(1, 20, 1),
                    minclustersize = 50, 
                    mergecutheight = 0.2, 
                    deepSplit = 2,
                    TOMType = 'signed',
                    maxblocksize = 5000, 
                    returntom = TRUE, 
                    plot = TRUE, 
                    featuretype = 'gene', 
                    plottitleprefix = NULL, 
                    threads = 6, 
                    seed = 2022)
names(wgcnares)
eigen_genes = wgcnares$mergedmelist

wgcnares$colorlabels = WGCNA::labels2colors(wgcnares$numericlabels)
names(wgcnares$colorlabels) = names(wgcnares$numericlabels)

anovares = featuresampling(betas = expr_df,
                           pddat = trait_df_filtered,
                           topfeatures = nrow(expr_df),
                           anova = TRUE,
                           plotannovares = TRUE,
                           plottitlesuffix = "Diagnosis",
                           titlesize = 18, 
                           textsize = 16,
                           threads = 6)

balanceddats = balancesampling(dat = expr_df,
                               pddat = trait_df_filtered,
                               responsename = "Diagnosis",
                               nround = 10,
                               seed = 2022,
                               k = 5,
                               threads = 6,
                               samplingmethod = 'updn')
names(balanceddats)

diff_results = diffmodules(dats=balanceddats$dats,
                           pddats=balanceddats$pds,
                           wgcnares=wgcnares,
                           
                           responsevarname = 'Diagnosis',
                           confoundings = c("sex", "smoking"),
                           
                           balanceadj=TRUE,
                           samplingmethod = 'updn',
                           nround = 10,
                           seed = 2022,
                           k = 5,
                           threads = 6,
                           method = 'Cauchy',
                           plot=TRUE,
                           
                           diffcutoff=0,
                           pvalcutoff=0.05,
                           pvalcolname='adj.P.Val')

wgcnares = diff_results$wgcnares
length(wgcnares$mergedmelist)
dim(wgcnares$mergedmelist[[1]])

diff_gene_results = diffgenes(dats=balanceddats$dats,
                              pddats=balanceddats$pds, 
                              wgcnares=wgcnares,
                              
                              topn=FALSE,
                              balanceadj=TRUE,
                              nround = 10,
                              
                              responsevarname = 'Diagnosis',
                              confoundings = c("sex", "smoking"),
                              
                              plot=TRUE,
                              titleprefix='IPF',
                              featuretype = "gene",
                              
                              diffcutoff=0,
                              pvalcutoff=0.05,
                              pvalcolname='adj.P.Val',
                              saveplot = TRUE,
                              savedir = 'CWGCNA_Results/Plots/')
names(diff_gene_results)  

dim(diff_gene_results$melimmares)  
dim(diff_gene_results$melimmasig) 
table(diff_gene_results$melimmasig$ME)
table(wgcnares$numericlabels)

sig_mods = unique(diff_gene_results$melimmasig$ME)

combined_info = data.frame()
for(module in sig_mods){
  print(paste('******', module, '******'))
  filename = paste('CWGCNA_Results/Mediation_results_IPF_', module, '.txt', sep='')
  mediationres = medgenes(dats=balanceddats$dats,
                          pddats=balanceddats$pds, 
                          wgcnares=wgcnares,
                          limmares=diff_gene_results$melimmasig,
                          modules = c(module),
                          threads = 6,
                          topn=FALSE,
                          balanceadj=TRUE,
                          nround = 10,
                          
                          responsevarname = 'Diagnosis',
                          confoundings = c('sex', 'smoking'),
                          
                          plot=TRUE,
                          titleprefix='IPF',
                          featuretype = "gene",
                          labeltype = 'Both',
                          labelnum = 5,
                          annotextsize = 3,
                          
                          diffcutoff=0,
                          pvalcutoff=0.05,
                          pvalcolname='adj.P.Val',
                          saveplot = TRUE,
                          savedir = 'CWGCNA_Results/Plots/',
                          filterres = FALSE)
  mediationres = mediationres$medres
  print(dim(mediationres))
  if(dim(mediationres)[1]==0)
    next
  
  #mediationres$Direction = apply(mediationres, 1, map_direction)
  write.table(mediationres, filename, sep='\t')
  combined_info = rbind(combined_info, mediationres)
}
dim(combined_info)
write.table(combined_info, 'CWGCNA_Results/Mediation_results_IPF.txt', sep='\t')

sig_mods = unique(diff_gene_results$melimmasig$ME)
combined_info = data.frame()
for(module in sig_mods){
  filename = paste('CWGCNA_Results/Mediation_results_IPF_', module, '.txt', sep='')
  if(!file.exists(filename))
    next
  mediationres = read.csv(filename, sep='\t')
  mediationres$Direction = apply(mediationres, 1, map_direction)
  combined_info = rbind(combined_info, mediationres)
}
dim(combined_info)

all_fwd_meds = combined_info[combined_info$Direction=='Forward',]$mediator
all_rev_meds = combined_info[combined_info$Direction=='Reverse',]$mediator
common_meds = intersect(all_fwd_meds, all_rev_meds)
map_unique = function(x){
  symbol = x[['mediator']]
  if(symbol %in% common_meds)
    return("Both")
  else
    return("Unique")
}

combined_info$Unique = apply(combined_info, 1, map_unique)

#computing gene module membership as correlation between individual gene expression and module eigengene
gene_membership = as.data.frame(WGCNA::cor(as.data.frame(t(expr_df)), eigen_genes, use="p"))
mm_pval = as.data.frame(WGCNA::corPvalueStudent(as.matrix(gene_membership), ncol(expr_df)))
dim(gene_membership)

#computing correlation between individual gene expression and clinical traits
gs = as.data.frame(WGCNA::cor(t(expr_df[,rownames(trait_df_filtered)]),
                                            as.data.frame(trait_df_filtered),use="p"))
gs$sampleid = NULL
gs_pval = as.data.frame(WGCNA::corPvalueStudent(as.matrix(gs),
                                         length(rownames(trait_df_filtered))))

gs = signif(gs, 3)
gs_pval_adj = as.data.frame(apply(gs_pval, 2, p.adjust, method='BH'))
names(gs_pval_adj) = paste(names(gs_pval_adj), "_adj", sep="")
names(gs_pval) = paste(names(gs_pval), "_pval_ctrl", sep="")
gs_pval = signif(gs_pval, 4)
gs_pval_adj = signif(gs_pval_adj, 4)

combined_info$Module_num = paste("ME", as.vector(wgcnares$numericlabels[combined_info$mediator]), sep='')
combined_info$Module_color = as.vector(wgcnares$colorlabels[combined_info$mediator])

combined_info['Module_membership'] = apply(combined_info,1,get_membership)
combined_info['Membership_pval'] = apply(combined_info,1,get_membership_pval)

combined_info = cbind(combined_info, 
                      diff_gene_results$melimmares[combined_info$mediator,c('logFC', 'P.Value', 'adj.P.Val')],
                      gs[combined_info$mediator,], 
                      gs_pval[combined_info$mediator,], 
                      gs_pval_adj[combined_info$mediator,])
write.table(combined_info, 'CWGCNA_Results/Mediation_results_IPF.txt', sep='\t')

############ Without WGCNA ########################

load(file = "WGCNA_Results/RData/network_construction.RData")
table(merge_colors)

# converting module colors to numeric labels
named_list = seq(1:length(names(table(merge_colors))))
names(named_list) = names(table(merge_colors))

numeric_labels = c()
for(color in merge_colors){
  numeric_labels = c(numeric_labels, named_list[[color]])
}
names(numeric_labels) = row.names(expr_df)

wgcnares = list(numericlabels = numeric_labels)

anovares = featuresampling(betas = expr_df,
                           pddat = trait_df_filtered,
                           topfeatures = nrow(expr_df),
                           anova = TRUE,
                           plotannovares = TRUE,
                           plottitlesuffix = "Diagnosis",
                           titlesize = 18, 
                           textsize = 16,
                           threads = 6)

balanceddats = balancesampling(dat = expr_df,
                               pddat = trait_df_filtered,
                               responsename = "Diagnosis",
                               nround = 10,
                               seed = 2022,
                               k = 5,
                               threads = 6,
                               samplingmethod = 'updn')
names(balanceddats)
length(balanceddats$dats)

diff_results = diffmodules(dats=balanceddats$dats,
                           pddats=balanceddats$pds,
                           wgcnares=wgcnares,
                           
                           responsevarname = 'Diagnosis',
                           confoundings = c("sex", "smoking"),
                           
                           balanceadj=TRUE,
                           samplingmethod = 'updn',
                           nround = 10,
                           seed = 2022,
                           k = 5,
                           threads = 6,
                           method = 'Cauchy',
                           plot=TRUE,
                           
                           diffcutoff=0,
                           pvalcutoff=0.05,
                           pvalcolname='adj.P.Val')

names(diff_results)
names(diff_results$wgcnares)

wgcnares = diff_results$wgcnares

diff_gene_results = diffgenes(dats=balanceddats$dats,
                              pddats=balanceddats$pds, 
                              wgcnares=wgcnares,
                              
                              topn=FALSE,
                              balanceadj=TRUE,
                              nround = 10,
                              
                              responsevarname = 'Diagnosis',
                              confoundings = c("sex", "smoking"),
                              
                              #plot=TRUE,
                              plot = FALSE,
                              titleprefix='IPF',
                              featuretype = "gene",
                              
                              diffcutoff=0,
                              pvalcutoff=1.0,
                              pvalcolname='adj.P.Val',
                              #saveplot = TRUE,
                              saveplot = FALSE,
                              savedir = 'CWGCNA_Results/Plots/')
names(diff_gene_results)  


sig_mods = unique(diff_gene_results$melimmasig$ME)

# Figure 3
combined_info = data.frame()
for(module in sig_mods){
  print(paste('******', module, '******'))
  filename = paste('CWGCNA_Results/Mediation_results_IPF_', module, '.txt', sep='')
  print(filename)
  mediationres = medgenes(dats=balanceddats$dats,
                          pddats=balanceddats$pds, 
                          wgcnares=wgcnares,
                          limmares=diff_gene_results$melimmasig,
                          modules = c(module),
                          threads = 6,
                          topn=FALSE,
                          balanceadj=TRUE,
                          nround = 10,
                          
                          responsevarname = 'Diagnosis',
                          confoundings = c('sex', 'smoking'),
                          
                          plot=TRUE,
                          titleprefix='IPF',
                          featuretype = "gene",
                          labeltype = 'Both',
                          labelnum = 5,
                          annotextsize = 4,
                          
                          diffcutoff=0,
                          pvalcutoff=0.05,
                          pvalcolname='adj.P.Val',
                          saveplot = TRUE,
                          savedir = 'CWGCNA_Results/Plots/', # Figure 3
                          filterres = FALSE)
  mediationres = mediationres$medres
  if(dim(mediationres)[1]==0)
    next
  
  mediationres$Direction = apply(mediationres, 1, map_direction)
  write.table(mediationres, filename, sep='\t')
  combined_info = rbind(combined_info, mediationres)
}
write.table(combined_info, 'CWGCNA_Results/Mediation_results_IPF.txt', sep='\t')

all_fwd_meds = combined_info[combined_info$Direction=='Forward',]$mediator
all_rev_meds = combined_info[combined_info$Direction=='Reverse',]$mediator
common_meds = intersect(all_fwd_meds, all_rev_meds)
map_unique = function(x){
  symbol = x[['mediator']]
  if(symbol %in% common_meds)
    return("Both")
  else
    return("Unique")
}

combined_info$Unique = apply(combined_info, 1, map_unique)

modules_info = read.csv('WGCNA_Results/Modules_info.txt', sep='\t')
row.names(modules_info) = as.vector(modules_info$Gene)

combined_info = cbind(combined_info, modules_info[as.vector(combined_info$mediator),])
# Supplementary File 2
write.table(combined_info, 'CWGCNA_Results/Mediation_results_IPF.txt', sep='\t')
