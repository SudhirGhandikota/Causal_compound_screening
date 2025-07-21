setwd('/path/to/directory/')
library(dynamicTreeCut)
library(cluster)
library(flashClust)
library(Hmisc)
library(reshape)
library(foreach)
library(Rcpp)
library(stringi)
library(WGCNA)
library(impute)
library(doParallel)
library(MASS)
library(gplots)
library(limma)
library(vsn)
library(stringr)
registerDoParallel()
enableWGCNAThreads(nThreads = 20)
#to not treat strings as factors -> data.frame()
options(stringsAsFactors = FALSE)

# this function plots labeled heatmap for the given inputs
plot_mod_heatmap = function(mod_trait_cor,mod_trait_pval,cex.text = 0.6,cex.lab.y = 0.6,
                            cex.lab.x = 0.6,mar = c(6,6,2,3),xLabelsAngle = 45,xLabelsAdj = 1.0,
                            sep=' ', plotLegend=TRUE, res=300, width=800, height=600,
                            main='', filename=NA){
  sizeGrWindow(12,6)
  # saving the plot to a file if a name is provided
  if(!is.na(filename)){
    if(grepl('tiff', filename))
      tiff(filename, width = width, height = height, res=res, units = 'px')
    if(grepl('png', filename))
      png(filename, width = width, height = height, res=res, units = 'px')
  }
    
  #correlations and p-values
  text_mat = paste(signif(mod_trait_cor,2), sep, "(",signif(mod_trait_pval,1),")",sep="")
  dim(text_mat) = dim(mod_trait_cor)
  par(mar=mar) #c(bottom,left,top,right)
  
  #display the correlations within a heatmap
  labeledHeatmap(Matrix = mod_trait_cor,
                 xLabels = colnames(mod_trait_cor),
                 yLabels = row.names(mod_trait_cor),
                 ySymbols = row.names(mod_trait_cor),
                 colorLabels = FALSE,colors=blueWhiteRed(50),
                 textMatrix = text_mat,setStdMargins = FALSE,
                 cex.text = cex.text, zlim=c(-1,1),
                 cex.lab.y = cex.lab.y, cex.lab.x = cex.lab.x,
                 xLabelsAdj = xLabelsAdj,
                 xLabelsAngle = xLabelsAngle, main='', 
                 plotLegend = plotLegend)
  if(!is.na(filename))
    dev.off()
  
}
diff_test_limma = function(controls, cases){
  expr_df = t(expr_df_filtered[c(controls,cases),])
  labels = c(rep(1,length(controls)),rep(2,length(cases)))
  fit = lmFit(expr_df,design = model.matrix(~labels))
  rownames(expr_df)
  colnames(fit)
  fit = eBayes(fit)
  top_genes = topTable(fit,coef = 2,number = length(colnames(expr_df_filtered)))
  return(top_genes)
}

#this function performs a two sample t-test to identify if a given gene is differentially expressed or not
diff_test = function(x, y){
  # x -> expression vectors of cases
  # y -> expression vectors of controls
  # t.test requires a minimum of 2 observations from each group
  if((sum(!is.na(x)) < 2) | (sum(!is.na(y)) < 2))
    return(data.frame(FC = NA, Pval = 1.0)) 
  res = t.test(x, y, alternative = "two.sided")
  pval = res[['p.value']]
  # these are log2 fold changes since the VSN method uses log2 transformation on the counts data.
  fold_change = signif(res[['estimate']][1] - res[['estimate']][2],2)
  return(data.frame(FC = unlist(fold_change[[1]]),
                    Pval = unlist(pval[1])))
}

get_membership = function(x){
  module = x[['Module']][1]
  if(module=='grey')
    return(NA)
  return(signif(as.numeric(gene_membership[ x[['Gene']][1], x[['Module']][1] ]),2))
}
get_membership_pval = function(x){
  module = x[['Module']][1]
  if(module=='grey')
    return(NA)
  return(signif(as.numeric(mm_pval[ x[['Gene']][1], x[['Module']][1] ]),3))
}

expr_df = read.csv('datasets/expr_data.txt', sep='\t')
dim(expr_df)

row.names(expr_df) = as.vector(expr_df$Gene_Symbol)
expr_df$Gene_Symbol = NULL

# normalizing raw RNAseq counts
expr_df_filtered = as.data.frame(voom(expr_df)$E)
dim(expr_df_filtered)

expr_df_filtered = as.data.frame((t(expr_df_filtered)))
row.names(expr_df_filtered)
colnames(expr_df_filtered)

trait_df = read.csv('datasets/trait_data.txt', sep='\t')
row.names(trait_df) = as.vector(trait_df$sample)
trait_df = trait_df[row.names(expr_df_filtered),]
names(trait_df)
table(trait_df$Diagnosis)

trait_df$sample = NULL
selected_traits = c('sex', 'age', 'Diagnosis', 'smoking', 'FEV', 'FVC', 'DLCO')
trait_df_filtered = trait_df[,selected_traits]

#removing samples and proteins with too many missing values
expr_df_good = goodSamplesGenes(expr_df_filtered, verbose=10)
sum(!expr_df_good$goodGenes)
if(!expr_df_good$allOK){
  #removing the "bad" samples and genes
  expr_df_filtered = expr_df_filtered[,expr_df_good$goodGenes]
}
rm(expr_df_good)

#clustering the samples to identify any outliers
clusters = hclust(dist(expr_df_filtered), method="average") #average link clustering based on euclidean distance
sizeGrWindow(12,9) #setting dimensions for the graphic device window
par(cex=0.75) #magnification parameter
par(mar=c(0,5,2,0)) #setting margins

labels = ifelse(clusters$labels %in% ipf_cases, 'IPF', 'Ctrl')
clusters$labels = paste(clusters$labels, labels, sep="-")
plot(clusters, main="Clustering the sample to detect outliers", sub="",xlab=""
     ,cex.lab=1.5, cex.axis=1.3, cex = 1.0, cex.main=1.5)
rm(clusters)
dev.off()

# scale free topology analysis -> 
# plotting the k(number of neighbors) vs p(k) (frequency distribution of nodes having k neighbors)
powers = c(c(1:10), seq(from=12, to=30, by=2))
sft = pickSoftThreshold(expr_df_filtered, powerVector = powers, verbose=5)

sizeGrWindow(9,5)
par(mfrow=c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threhold(power)",
     ylab="Scale Free Topology Model Fit, signed R^2",type="n",main=paste("Scale independence-All genes"))
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=0.9,col="red")
#R^2 cutoff line
abline(h=0.9,col="red")
#mean connectivity vs soft-thresholding power
plot(sft$fitIndices[,1],sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity",type="n",
     main=paste("Mean connectivity-All genes"))
text(sft$fitIndices[,1],sft$fitIndices[,5], labels=powers, cex=0.9, col="red")
dev.off()

# Figure 2A
sizeGrWindow(6,4)
par(mfrow=c(1,1))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab="Soft Threhold(power)",
     ylab="Scale Free Topology Model Fit, signed R^2",type="n",
     main=paste("Scale independence-All genes"),
     cex.axis=1.2, cex.lab=1.2)
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1.2,col="red")
#R^2 cutoff line
abline(h=0.9,col="red")
dev.off()

rm(sft)

#Step-by-step network construction and module detection
beta = 8 #using Scale Free Topology Model analysis in the earlier section
#computes the adjacency matrix using power function(soft thresholding). Pearson correlation used for similarity
adj_mat = adjacency(expr_df_filtered, power = beta, type="signed")
TOM_sim = TOMsimilarity(adj_mat, TOMType = "signed")
genes = colnames(expr_df_filtered)

save(genes,TOM_sim,file="datasets/RData/TOM.RData")
dist_TOM = 1 - TOM_sim
rm(adj_mat)

load(file="datasets/RData/TOM.RData")
#hierarchical clustering of the TOM distance matrix
gene_tree = hclust(as.dist(dist_TOM), method="average")

#plot the resulting clustering dendrogram
sizeGrWindow(12,9)
plot(gene_tree, xlab="", sub="", main="Gene Clustering on TOM-based dissimilarity",
     labels=FALSE, hang=0.04)
dev.off()

dist_TOM = 1 - TOM_sim

#cross-validation with different module sizes i.e number of genes in each module
mod_size = 100
modules = cutreeDynamic(dendro = gene_tree, distM = dist_TOM,
                        deepSplit = 2, pamRespectsDendro = FALSE,minClusterSize = mod_size)

#deepsplit controls sensitivity to cluster splitting 
num_of_mods = length(unique(modules))
print(paste("Number of Modules: ",num_of_mods))

resultColors = labels2colors(modules)
table(modules) #module 0 is the set of genes unassigned to any modules
table(resultColors) #module "grey" is the set of genes unassigned to any modules

sizeGrWindow(8,6)
title = paste("Protein dendrogram and module colors(",num_of_mods,")")
plotDendroAndColors(gene_tree,resultColors,"Modules",dendroLabels = FALSE, 
                    hang=0.03, addGuide = TRUE, guideHang = 0.05,
                    cex.colorLabels = 1.5, cex.dendroLabels = 1.5, main = title)
dev.off()
rm(dist_TOM)

# merging modules whose expression profiles are very similar
#finding Module eigen-genes
MEs = moduleEigengenes(expr_df_filtered,colors = resultColors, impute = FALSE)
eigen_genes = MEs$eigengenes
# removing 'MEgrey' with all NA values
eigen_genes = eigen_genes[colnames(eigen_genes)!='MEgrey']

# #clustering module eigen-genes
dist_ME = 1 - cor(eigen_genes)
ME_tree = hclust(as.dist(dist_ME), method='average')
# 
# #plotting the module eigen-genes clusters
sizeGrWindow(7,6)
plot(ME_tree, main="Clustering of module eigengenes", xlab="",sub="")
threshold = 0.1 #distance = 0.25 => similarity = 0.75
# #add cut lines into the dendrogram using the set threshold = 0.25
abline(h=threshold, col="red")
dev.off()

#merging function
merge_mods = mergeCloseModules(expr_df_filtered,resultColors,
                               cutHeight = threshold,verbose=3)
merge_colors = merge_mods$colors
#eigen-genes of the new merged modules
merged_mes = merge_mods$newMEs
print(paste("Number of modules after merge: ",length(table(merge_colors)),sep=""))
table(merge_colors)

save(merged_mes,merge_colors,gene_tree, 
     file="WGCNA_Results/RData/network_construction.RData")

load(file="WGCNA_Results/RData/network_construction.RData")
table(merge_colors)

# Figure 2B
sizeGrWindow(8,6)
par(cex.lab = 1.5)
title = paste("Gene dendrogram and module colors(",
              length(table(merge_colors)),")")
plotDendroAndColors(gene_tree,merge_colors,"Modules",dendroLabels = FALSE, 
                    hang=0.03, addGuide = TRUE, guideHang = 0.05,
                    cex.colorLabels = 1.5, cex.dendroLabels = 1.5,
                    cex.rowText=1.5, 
                    main = title)
dev.off()

#computing correlations between eigen-genes and traits
mod_trait_cor = cor(merged_mes[rownames(trait_df_filtered),], trait_df_filtered, use="p")
#calculating the student asymptotic p-value
mod_trait_pval = corPvalueStudent(mod_trait_cor, length(row.names(trait_df_filtered)))

mod_names = substring(row.names(mod_trait_cor), 3)
#mod_names = paste(mod_names, '_base', sep='')
row.names(mod_trait_cor) = mod_names
row.names(mod_trait_pval) = mod_names
colnames(mod_trait_cor)

save(mod_trait_cor,mod_trait_pval, file="datasets/RData/mod_cor_info.RData")

mod_trait_cor = mod_trait_cor[order(-mod_trait_cor[,3]),]
mod_trait_pval = mod_trait_pval[row.names(mod_trait_cor),]
plot_mod_heatmap(mod_trait_cor, mod_trait_pval, 
                 cex.lab.x = 1.0, cex.lab.y = 1.0, cex.text = 1.0,
                 mar = c(3.5,5,1.5,1), #c(bottom,left,top,right)
                 sep='\n', main = "Module-Trait relations")

# Figure 3C
mod_trait_cor_fil = mod_trait_cor[,c('Diagnosis', 'FVC', 'DLCO')]
colnames(mod_trait_cor_fil) = c('IPF', 'FVC', 'DLCO')
mod_trait_pval_fil = mod_trait_pval[row.names(mod_trait_cor_fil),
                                    c('Diagnosis', 'FVC', 'DLCO')]
colnames(mod_trait_pval_fil) = c('IPF', 'FVC', 'DLCO')
plot_mod_heatmap(mod_trait_cor_fil, mod_trait_pval_fil, 
                 cex.lab.x = 1.0, cex.lab.y = 1.0, cex.text = 1.0,
                 mar = c(2.5,5,1.5,1), #c(bottom,left,top,right)
                 sep='\n', main = "")
plot_mod_heatmap(mod_trait_cor_fil, mod_trait_pval_fil, 
                 cex.lab.x = 1.0, cex.lab.y = 1.0, cex.text = 1.0,
                 mar = c(2.5,5,1.5,1), #c(bottom,left,top,right)
                 sep='\n', main = "", width = 1250, height = 2100,
                 filename = 'WGCNA_Results/Plots/corr_heatmap_fil.png')

#computing gene module membership as correlation between individual gene expression and module eigengene
gene_membership = as.data.frame(cor(expr_df_filtered, merged_mes, use="p"))
mm_pval = as.data.frame(corPvalueStudent(as.matrix(gene_membership), nrow(expr_df_filtered)))
names(gene_membership) = substring(names(gene_membership), 3)
names(mm_pval) = substring(names(mm_pval), 3)
dim(gene_membership)

#generating the excel file with all genes, module membership and correlation
symbols = names(expr_df_filtered)
#computing correlation between individual gene expression and clinical traits
gs = as.data.frame(cor(expr_df_filtered[rownames(trait_df_filtered),],
                                         as.data.frame(trait_df_filtered),use="p"))
gs_pval = as.data.frame(corPvalueStudent(as.matrix(gs), length(rownames(trait_df_filtered))))

gs = signif(gs, 3)
gs_pval_adj = as.data.frame(apply(gs_pval, 2, p.adjust, method='BH'))
names(gs_pval_adj) = paste(names(gs_pval_adj), "_adj", sep="")
names(gs_pval) = paste(names(gs_pval), "_pval", sep="")
gs_pval = signif(gs_pval, 4)
gs_pval_adj = signif(gs_pval_adj, 4)

#generating the excel file with all genes, module membership and correlation
symbols = names(expr_df_filtered)
all_genes_info = data.frame(Gene = symbols,
                            Module = as.vector(merge_colors))

all_genes_info['Module_membership'] = apply(all_genes_info,1,get_membership)
all_genes_info['Membership_pval'] = apply(all_genes_info,1,get_membership_pval)

#this function performs a two sample t-test to identify if a given gene is differentially expressed or not
get_fc = function(x, set){
  return(signif(as.numeric(top_genes[[set]][x[['Gene']],'logFC']),3))
}
get_adj_pval = function(x, set){
  return(as.numeric(top_genes[[set]][x[['Gene']],'adj.P.Val']))
}
get_pval = function(x, set){
  return(as.numeric(top_genes[[set]][x[['Gene']],'P.Value']))
}

top_genes = list()
top_genes[[1]] = diff_test_limma(controls, ipf_cases)
names(top_genes[[1]])

all_genes_info['logFC'] = apply(all_genes_info,1,get_fc, set = 1)
all_genes_info['adj_pval'] = apply(all_genes_info,1,get_adj_pval, set = 1)
all_genes_info['pval'] = apply(all_genes_info,1,get_pval, set = 1)

all_genes_info = cbind(all_genes_info, gs)
all_genes_info = cbind(all_genes_info, gs_pval)
all_genes_info = cbind(all_genes_info, gs_pval_adj)

all_genes_info = all_genes_info[,c(names(all_genes_info)[1:7], c(c(rbind(names(gs), names(gs_pval), names(gs_pval_adj)))))]
names(all_genes_info)

# Supplementary File 1
write.table(all_genes_info,file="WGCNA_Results/WGCNA_Modules_info.txt", row.names = F,sep="\t")
