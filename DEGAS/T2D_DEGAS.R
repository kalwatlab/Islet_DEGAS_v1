##############################################################
# The scRNAseq datasets are provided by A Mawla
# Setting up R environment
library(biomaRt)
library(dplyr)
library(Seurat)
library(DEGAS)
library(Rtsne)
library(ggplot2)
library(pROC)
library(tidyverse)
library(doParallel) # Added for parrallel marker selection TSJ 20230401
n_cores = 4 # Added for parallel marker selection TSJ 20230401
registerDoParallel(cores=n_cores) # Added for parrallel marker selection TSJ 20230401
work.dir = "/vol/data02/groy_data/Single_cell_RNA_seq/scRNAmerge_vs_Marselli/T2D/"


# Defining custom functions
selectFeats <- function(expression,features,selection_statistic)
{
  original_feats = row.names(expression)
  row.names(expression) = 1:dim(expression)[1]
  dup_feat_uniq = unique(features[duplicated(features)])
  message(length(dup_feat_uniq))
  if(length(dup_feat_uniq)>0){
    dup_feat = features[features %in% dup_feat_uniq]
    dup_feat_idx = which(features %in% dup_feat_uniq)
    rem_feat_idx = c()
    for(feat in dup_feat){
      feat_rowSums = apply(expression[dup_feat_idx[dup_feat==feat],],1,eval(parse(text=selection_statistic)))
      max_feat_idx = which(feat_rowSums==max(feat_rowSums))[1]
      rem_feat_idx = c(rem_feat_idx,as.numeric(names(feat_rowSums)[-max_feat_idx]))
    }
    expression = expression[-rem_feat_idx,]
    row.names(expression) = features[-rem_feat_idx]
    feat_df <- data.frame(new_feature=row.names(expression),original_feature=original_feats[-rem_feat_idx])
    row.names(feat_df) <- feat_df$new_feature
    return(list(expression,feat_df))
  }else{
    row.names(expression) = features
    feat_df <- data.frame(new_feature=row.names(expression),original_feature=original_feats[-rem_feat_idx])
    row.names(feat_df) <- feat_df$new_feature
    return(list(expression,feat_df))
  }
}

subcluster <- function(seurat.obj,ct){
  seurat.obj = subset(seurat.obj,subset=cell_type==ct)
  seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)
  seurat.obj <- ScaleData(seurat.obj,features=rownames(seurat.obj))
  seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj))
  seurat.obj <- RunUMAP(seurat.obj, dims=1:10)
  seurat.obj <- FindNeighbors(seurat.obj,reduction="umap",dims=1:2)
  seurat.obj <- FindClusters(seurat.obj, resolution = 0.8)
  return(seurat.obj)
}

getPval = function(x,g){
  tmp = kruskal.test(x~(g))
  return(tmp$p.value)
}

getLog2FC = function(x,g){
  out = log2(mean(x[g])/mean(x[!g]))
  return(out)
}

FindClusterMarkers = function(X,g,n){
  MarkerTable = list()
  for(grp in unique(g)){
    Pvals = apply(X,2,function(x) getPval(x,g==grp))
    Log2FCs = apply(X,2,function(x) getLog2FC(x,g==grp))
    FDRs = p.adjust(Pvals,method="BH")
    NegLog10_Pvals = -log10(Pvals)
    NegLog10_FDRs = -log10(FDRs)
    FDR_Ranks = rank(FDRs)
    Log2FC_Ranks = rank(-Log2FCs)
    FDR_Log2FC_Ranksums = FDR_Ranks+Log2FC_Ranks
    MarkerTable[[grp]] = data.frame(Gene = names(Pvals),Cluster=rep(grp,length(Pvals)),Log2FC=Log2FCs,Pval=Pvals,FDR=FDRs,NegLog10_Pval=NegLog10_Pvals,
                                    NegLog10_FDR=NegLog10_FDRs,Log2FC_Rank=Log2FC_Ranks,FDR_Rank=FDR_Ranks,FDR_Log2FC_Ranksum=FDR_Log2FC_Ranksums)
    MarkerTable[[grp]] = arrange(MarkerTable[[grp]],FDR_Log2FC_Ranksum) %>% top_n(-n)
  }
  MarkerTable = do.call(rbind,MarkerTable)
  return(MarkerTable)
}

calcAndSort <- function(X,g,n,grp){
  Pvals = apply(X,2,function(x) getPval(x,g==grp))
  Log2FCs = apply(X,2,function(x) getLog2FC(x,g==grp))
  FDRs = p.adjust(Pvals,method="BH")
  NegLog10_Pvals = -log10(Pvals)
  NegLog10_FDRs = -log10(FDRs)
  FDR_Ranks = rank(FDRs)
  Log2FC_Ranks = rank(-Log2FCs)
  FDR_Log2FC_Ranksums = FDR_Ranks+Log2FC_Ranks
  tmp = data.frame(Gene = names(Pvals),Cluster=rep(grp,length(Pvals)),Log2FC=Log2FCs,Pval=Pvals,FDR=FDRs,NegLog10_Pval=NegLog10_Pvals,
                   NegLog10_FDR=NegLog10_FDRs,Log2FC_Rank=Log2FC_Ranks,FDR_Rank=FDR_Ranks,FDR_Log2FC_Ranksum=FDR_Log2FC_Ranksums)
  rownames(tmp) <- NULL
  return(arrange(tmp,FDR_Log2FC_Ranksum) %>% top_n(-n))
}

FindClusterMarkersPar = function(X,g,n){
  MarkerTable = list()
  MarkerTable <- foreach(grp = unique(g)) %dopar% calcAndSort(X,g,n,grp)
  MarkerTable = do.call(rbind,MarkerTable)
  return(MarkerTable)
}


###########################################################################################################################
#                                          Loading Bulk data and preprocessing                                            #
###########################################################################################################################
bulk.expr = read.csv(paste0(work.dir,"bulk.expr.csv"), sep = ",", header = T, row.names = "symbol")
bulk.expr$X = NULL


###########################################################################################################################
#                                          Loading Bulk metadata                                                          #
###########################################################################################################################

bulk.meta <-  read.csv(paste0(work.dir,"marselli_GSE159984_metadata_clean.csv"), sep = ",", header = T, row.names = "SAMN_and_disease")

setdiff(rownames(bulk.meta),colnames(bulk.expr))
colnames(bulk.expr) = gsub("[.]","-",colnames(bulk.expr))
setdiff(rownames(bulk.meta),colnames(bulk.expr))
#bulk.expr = bulk.expr[,rownames(bulk.meta)]
identical(rownames(bulk.meta),colnames(bulk.expr))
all.equal(rownames(bulk.meta),colnames(bulk.expr))


###########################################################################################################################
#                                                       Loading SC data                                                   #
###########################################################################################################################

merged_seurat <- readRDS("~/data/Single_cell_RNA_seq/scRNAmerge_vs_Marselli/merged_seurat_AMawla_NEW.rds") 
sc.expr = as.data.frame(merged_seurat@assays$RNA@counts)


###Tabulate cells by cluster ID

library(data.table)
md <- merged_seurat@meta.data %>% as.data.table #ONE ROW PER CELL

## count the number of cells per unique combinations of "Sample" and "seurat_clusters"
md2 <- md[, .N, by = c("orig.ident", "seurat_clusters")]

## with additional casting after the counting
#md2 <- md[, .N, by = c("orig.ident", "seurat_clusters")] %>% dcast(., orig.ident ~ seurat_clusters, value.var = "N")
#write.csv(md2,file=paste0(work.dir,"beta cells per cluster.csv"),quote=FALSE)

md3 = md2 %>% group_by(seurat_clusters) %>% summarize(SumAvg = sum(N, na.rm = TRUE))



###########################################################################################################################
#                         Identifying final feature set (Changed this section TSJ 20230401)                               #
###########################################################################################################################

if(TRUE){
  # In this section we use the new parrallelized find cluster markers with
  # 800 top features for each cluster and patient group.
  int.feats = intersect(rownames(bulk.expr),rownames(sc.expr)) # how many is this? =
  sc.markers = FindClusterMarkersPar(t(sc.expr[int.feats,]),merged_seurat$seurat_clusters,800) 
  write.csv(sc.markers, file=paste0(work.dir,"sc.clustermarkers.csv"))
  tmp = t(bulk.expr[int.feats,])
  bulk.markers = FindClusterMarkersPar(tmp,bulk.meta$Disease,400) 
  write.csv(bulk.markers, file=paste0(work.dir,"bulk.clustermarkers.csv"))
  highvar.features = intersect(bulk.markers$Gene,sc.markers$Gene)
  highvar.features = intersect(unique(bulk.markers$Gene),unique(sc.markers$Gene))
  length(highvar.features)
  saveRDS(highvar.features,file=paste0(work.dir,"featureSet.rds"))
}



###########################################################################################################################
#                                Preprocessing data (log, normalization, scale)                                           #
###########################################################################################################################

scDat = preprocessCounts(sc.expr[highvar.features,]) #matrix of expression values from scRNAseq data
scLab = toOneHot(merged_seurat$seurat_clusters)
patDat = preprocessCounts(bulk.expr[highvar.features,]) #matrix of expression values from bulkRNAseq data, columns should be in same order as scDat
patLab = toOneHot(bulk.meta[colnames(bulk.expr),"Disease"])


###########################################################################################################################
#                                             Training DEGAS model                                                        #
###########################################################################################################################

# Intializing all of the hyperparameters such as batch size regularization weight, each loss term weight ...
initDEGAS()
# Setting a directory to write tmp files to
tmpDir = paste0(work.dir,'tmp/')
# Setting seed so that the model training is reproducible, to get the same answer in two different runs
set_seed_term(3)
# Default patient batch size is 50, since we have 115 patients, i.e. sampling all patients
set_patient_batch_size(50) 
set_single_cell_batch_size(500) 
# This writes the tensorflow model in python, runs the python script, then returns the model coefficients to R
DEGAS.model = runCCMTLBag(scDat,scLab,patDat,patLab,tmpDir,'ClassClass','DenseNet',3,5) #A trained bootstrap aggregated DEGAS model. 3 refers to the layers/complexity of neural network 
#and 5 is number of times to bootstrap aggregate the models. If increase the layers then bootstrap aggregation needs to be increased as well.

###########################################################################################################################
#                                               Predicting Diabetes associations                                          #
###########################################################################################################################

scpatPreds = predClassBag(DEGAS.model,scDat,"pat") #predClassBag(ccModel, Exp, scORpat)
colnames(scpatPreds) = colnames(patLab)
patpatPreds = predClassBag(DEGAS.model,patDat,"pat")
colnames(patpatPreds) = colnames(patLab)
scscPreds = predClassBag(DEGAS.model,scDat,"sc")
colnames(scscPreds) = colnames(scLab)
patscPreds = predClassBag(DEGAS.model,patDat,"sc")
colnames(patscPreds) = colnames(scLab)
# AUC in patients
roc(patLab[,1],patpatPreds[,1])
boxplot(patscPreds[,1]~patLab[,1])


###########################################################################################################################
#                                                   Post processing                                                       #
###########################################################################################################################

boxplot(scpatPreds[,"T2D"]~fromOneHot(scLab))
boxplot(scpatPreds[,"ND"]~fromOneHot(scLab))

T2D_assoc = knnSmooth(toCorrCoeff(scpatPreds[,"T2D"]),merged_seurat@reductions$umap@cell.embeddings)   # k-nearest neighbor smoothing (kNN-smoothing) algorithm, designed to reduce technical noise exhibited by UMI-filtered scRNA-Seq data by aggregating information from similar cells (neighbors) in a computationally efficient and statistically tractable manner. 
merged_seurat$T2D_assoc = T2D_assoc
merged_seurat$T2D_raw = toCorrCoeff(scpatPreds[,1])

A=FeaturePlot(merged_seurat, features = c("INS", "GCG", "SST", "PPY", "GHRL", "KRT19", "PRSS1", "REG1A", "CFTR", "COL6A1", "PDGFRA", "PLVAP","CD44", "KIT"), min.cutoff = "q9")
ggsave("A- Marker genes.svg", A, width = 12, height = 8, device = "svg")
ggsave("A- Marker genes.tiff", A, width = 12, height = 8, device = "tiff")

B= FeaturePlot(merged_seurat,features="T2D_assoc", pt.size = 1) + scale_color_gradient2(low = "black",mid="lavender",high="red")
ggsave("B- T2D_assoc.svg", B, width = 12, height = 8, device = "svg")
ggsave("B- T2D_assoc.tiff", B, width = 12, height = 8, device = "tiff")

FeaturePlot(merged_seurat,features="INS") + scale_color_gradient2(low = "black",mid="lavender",high="red")
FeaturePlot(merged_seurat,features="GCG") + scale_color_gradient2(low = "black",mid="lavender",high="red")
FeaturePlot(merged_seurat,features="SST") + scale_color_gradient2(low = "black",mid="lavender",high="red")
FeaturePlot(merged_seurat,features="PPY") + scale_color_gradient2(low = "black",mid="lavender",high="red")
FeaturePlot(merged_seurat,features="GHRL") + scale_color_gradient2(low = "black",mid="lavender",high="red")
boxplot(merged_seurat$T2D_assoc~merged_seurat$seurat_clusters)
DimPlot(merged_seurat, group.by = "seurat_clusters", label = T)+ NoLegend()

###########################################################################################################################
#                                        Diabetes associations in single cells                                            #
###########################################################################################################################

#set.seed(3)
scpatAssocSmooth = knnSmooth(toCorrCoeff(as.numeric(scpatPreds[,"T2D"])),as.matrix(merged_seurat@reductions$umap@cell.embeddings),k=20)
fig.df = data.frame(Diabetes = scpatAssocSmooth, UMAP1 = merged_seurat@reductions$umap@cell.embeddings[,"UMAP_1"],  UMAP2 = merged_seurat@reductions$umap@cell.embeddings[,"UMAP_2"], Cluster = fromOneHot(scLab),
                    GCG = log2(as.numeric(sc.expr["GCG",])+1), INS = log2(as.numeric(sc.expr["INS",])+1),
                    SST = log2(as.numeric(sc.expr["SST",])+1), MAFA = log2(as.numeric(sc.expr["MAFA",])+1),
                    PPY = log2(as.numeric(sc.expr["PPY",])+1), PRSS1 = log2(as.numeric(sc.expr["PRSS1",])+1),
                    KRT19 = log2(as.numeric(sc.expr["KRT19",])+1), MAFB = log2(as.numeric(sc.expr["MAFB",])+1),
                    CD73 = log2(as.numeric(sc.expr["NT5E",])+1), CD90 = log2(as.numeric(sc.expr["THY1",])+1))
ggplot(fig.df,aes(x=UMAP1,y=UMAP2,color=Cluster)) + geom_point()
ggplot(fig.df,aes(x=UMAP1,y=UMAP2,color=Diabetes)) + geom_point() + scale_color_gradient2(low = "black",mid="lavender",high="red") + ggtitle("T2D associated") + theme(plot.title = element_text(hjust = 0.5))
ggplot(fig.df,aes(x=UMAP1,y=UMAP2,color=GCG)) + geom_point() + scale_color_gradient(low = "black",high="red") + ggtitle("GCG") + theme(plot.title = element_text(hjust = 0.5))
ggplot(fig.df,aes(x=UMAP1,y=UMAP2,color=INS)) + geom_point() + scale_color_gradient(low = "black",high="red") + ggtitle("INS") + theme(plot.title = element_text(hjust = 0.5))
ggplot(fig.df,aes(x=UMAP1,y=UMAP2,color=SST)) + geom_point() + scale_color_gradient(low = "black",high="red") + ggtitle("SST") + theme(plot.title = element_text(hjust = 0.5))
ggplot(fig.df,aes(x=UMAP1,y=UMAP2,color=MAFA)) + geom_point() + scale_color_gradient(low = "black",high="red") + ggtitle("MAFA") + theme(plot.title = element_text(hjust = 0.5))
ggplot(fig.df,aes(x=UMAP1,y=UMAP2,color=PPY)) + geom_point() + scale_color_gradient(low = "black",high="red") + ggtitle("PPY") + theme(plot.title = element_text(hjust = 0.5))
ggplot(fig.df,aes(x=UMAP1,y=UMAP2,color=PRSS1)) + geom_point() + scale_color_gradient(low = "black",high="red") + ggtitle("PRSS1") + theme(plot.title = element_text(hjust = 0.5))
ggplot(fig.df,aes(x=UMAP1,y=UMAP2,color=KRT19)) + geom_point() + scale_color_gradient(low = "black",high="red") + ggtitle("KRT19") + theme(plot.title = element_text(hjust = 0.5))
ggplot(fig.df,aes(x=UMAP1,y=UMAP2,color=MAFB)) + geom_point() + scale_color_gradient(low = "black",high="red") + ggtitle("MAFB") + theme(plot.title = element_text(hjust = 0.5))
ggplot(fig.df,aes(x=UMAP1,y=UMAP2,color=CD73)) + geom_point() + scale_color_gradient(low = "black",high="red") + ggtitle("CD73") + theme(plot.title = element_text(hjust = 0.5))
ggplot(fig.df,aes(x=UMAP1,y=UMAP2,color=CD90)) + geom_point() + scale_color_gradient(low = "black",high="red") + ggtitle("CD90") + theme(plot.title = element_text(hjust = 0.5))


###########################################################################################################################
#                                                      Subsetting beta cells                                              #
###########################################################################################################################

FeaturePlot(merged_seurat, features = "INS")
boxplot(merged_seurat@assays$RNA@counts["INS",]~merged_seurat$seurat_clusters)

merged_seurat.beta <- subset(merged_seurat, subset = seurat_clusters %in% c(3,7,13,19,20,21))
DefaultAssay(merged_seurat.beta) = "integrated"    ## In integrated assay, you do not normalize the data
merged_seurat.beta <- ScaleData(merged_seurat.beta,features=rownames(merged_seurat.beta))
merged_seurat.beta <- RunPCA(merged_seurat.beta, features = VariableFeatures(object = merged_seurat.beta))
merged_seurat.beta <- RunUMAP(merged_seurat.beta, dims=1:10)
merged_seurat.beta <- FindNeighbors(merged_seurat.beta,reduction="umap",dims=1:2)
merged_seurat.beta <- FindClusters(merged_seurat.beta)
merged_seurat.beta$DiabetesBeta = knnSmooth(merged_seurat.beta$T2D_assoc,merged_seurat.beta@reductions$umap@cell.embeddings[,c("UMAP_1","UMAP_2")],20)


C= DimPlot(merged_seurat.beta, label = T, pt.size = 2, label.size = 5)+ ggtitle("Beta cell clusters")+ theme(plot.title = element_text(hjust = 0.5))
ggsave("C- Beta cell clusters.svg", C, width = 12, height = 8, device = "svg")
ggsave("C- Beta cell clusters.tiff", C, width = 12, height = 8, device = "tiff")

D1= FeaturePlot(merged_seurat.beta,features="T2D_assoc", label = FALSE, pt.size = 2) + scale_color_gradient2(low = "black",mid="lavender",high="red")+ ggtitle("T2D associated beta cells")
ggsave("D1- T2D associated beta cells.svg", D1, width = 12, height = 8, device = "svg")
ggsave("D1- T2D associated beta cells.tiff", D1, width = 12, height = 8, device = "tiff")
D2= FeaturePlot(merged_seurat.beta,features="DiabetesBeta", label = FALSE, pt.size = 2) + scale_color_gradient2(low = "black",mid="lavender",high="red")+ ggtitle("T2D associated beta cells after smoothing")
ggsave("D2- T2D associated beta cells after smoothing.svg", D2, width = 12, height = 8, device = "svg")
ggsave("D2- T2D associated beta cells after smoothing.tiff", D2, width = 12, height = 8, device = "tiff")

betacluster.markers = FindAllMarkers(merged_seurat.beta)
write.csv(betacluster.markers,file=paste0(work.dir,"betacellcluster_markers.csv"),quote=FALSE)

merged_seurat.beta.donor <- merged_seurat.beta
Idents(merged_seurat.beta.donor) <- 'orig.ident'
Idents(merged_seurat.beta.donor)
merged_seurat.beta.donor <- RenameIdents(merged_seurat.beta.donor, "human1"="Baron","human2"="Baron", "human3"="Baron","human4"="Baron","lawlor"="Lawlor","SRR4003787"="Muraro",  "SRR4003788"="Muraro",  "SRR4003789"="Muraro",
                                         "SRR4003790"="Muraro",  "SRR4003791"="Muraro",  "SRR4003792"="Muraro",  "SRR4003793"="Muraro",  "SRR4003794"="Muraro", "SRR4003795"="Muraro",  "SRR4003796"="Muraro",  "SRR4003797"="Muraro",
                                         "SRR4003798"="Muraro",  "SRR4003799"="Muraro",  "SRR4003800"="Muraro",  "SRR4003801"="Muraro",  "SRR4003802"="Muraro",  "SRR4003803"="Muraro",  "SRR4003804"="Muraro",  "SRR4003805"="Muraro",
                                         "SRR4003806"="Muraro",  "SRR4003807"="Muraro",  "SRR4003808"="Muraro", "SRR4003809"="Muraro",  "SRR4003810"="Muraro",  "SRR4003811"="Muraro",  "SRR4003812"="Muraro",  "SRR4003813"="Muraro",
                                         "SRR4003814"="Muraro",  "SRR4003815"="Muraro",  "SRR4003816"="Muraro",  "SRR4003817"="Muraro",  "SRR4003818"="Muraro",  "segerstolpe"="Segerstolpe", "xin"="Xin") 
E1= DimPlot(merged_seurat.beta.donor, cols = c("Baron"="#F8766D", "Lawlor"="#ABA300", "Muraro"="#00A9FF", "Segerstolpe"="#C77CFF", "Xin"="#E68613"), pt.size = 2)+ ggtitle("Identity of beta cells based on source dataset")
Idents(merged_seurat.beta.donor) <- 'Donor'
levels(merged_seurat.beta.donor)
E2= DimPlot(merged_seurat.beta.donor, pt.size = 2)+ ggtitle("Donor identity of beta cells")
E= E1+E2
print(E)
ggsave("E- Donor identity of beta cells.svg", E, width = 20, height = 8, device = "svg")
ggsave("E- Donor identity of beta cells.tiff", E, width = 20, height = 8, device = "tiff")

###Tabulate cells by cluster ID
library(data.table)
md <- merged_seurat.beta.donor@meta.data %>% as.data.table #ONE ROW PER CELL

## count the number of cells per unique combinations of "Sample" and "seurat_clusters"
md2 <- md[, .N, by = c("orig.ident", "seurat_clusters")]

## with additional casting after the counting
#md2 <- md[, .N, by = c("orig.ident", "seurat_clusters")] %>% dcast(., orig.ident ~ seurat_clusters, value.var = "N")
#write.csv(md2,file=paste0(work.dir,"beta cells per cluster.csv"),quote=FALSE)

md3 = md2 %>% group_by(seurat_clusters) %>% summarize(SumAvg = sum(N, na.rm = TRUE))

###########################################################################################################################
#                               Taking the quantiles instead of median                                                    #
###########################################################################################################################


quants = quantile(merged_seurat.beta$DiabetesBeta,c(0.1,0.15,0.2,0.25,0.33,0.5,0.66,0.75,0.8,0.85,0.9))
names(quants) = c("10","15","20","25","33","50","66","75","80","85","90")

###########################################################################################################################
#                                    Defining high_mid_low groups                                                         #
###########################################################################################################################

highmidlow = ifelse(merged_seurat.beta$DiabetesBeta<quants["20"],"low","mid")
highmidlow[highmidlow=="mid"] = ifelse(merged_seurat.beta$DiabetesBeta[highmidlow=="mid"]>quants["80"],"high","mid")
merged_seurat.beta$highDiabetes3 = highmidlow
unique(highmidlow)

merged_seurat.beta.highmidlow <- merged_seurat.beta
Idents(merged_seurat.beta.highmidlow) = merged_seurat.beta.highmidlow$highDiabetes3

F= DimPlot(merged_seurat.beta.highmidlow, label = FALSE, pt.size = 2, label.size = 5)+ ggtitle("Beta cell clusters-high_mid_low")+ theme(plot.title = element_text(hjust = 0.5))
ggsave("F- Beta cell clusters-high_mid_low.svg", F, width = 12, height = 8, device = "svg")
ggsave("F- Beta cell clusters-high_mid_low.tiff", F, width = 12, height = 8, device = "tiff")

markers.to.plot <- c("RPS4Y1","SDF2L1","SIX3","SPP1","PTPRN2","BTG1","SYT7","PSAP")
Idents(merged_seurat.beta.highmidlow)
DotPlot(merged_seurat.beta.highmidlow, features = markers.to.plot, cols = c("orange","brown"), idents = c("high","low"), dot.scale = 8) +RotatedAxis()+coord_flip()+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=12,face="bold"))


### Calculating differential gene expression with BH FDR (comparing between high and low groups, ignoring the mid)

DefaultAssay(merged_seurat.beta) = "RNA"
merged_seurat.beta = NormalizeData(merged_seurat.beta,normalization.method = "LogNormalize")
highDiabetes.markers2 = FindMarkers(merged_seurat.beta,ident.1=colnames(merged_seurat.beta)[merged_seurat.beta$highDiabetes3=="high"],
                                    ident.2=colnames(merged_seurat.beta)[merged_seurat.beta$highDiabetes3=="low"])
highDiabetes.markers2[c("SDF2L1","QPCT","HADH","IAPP","NPTX2","CNIH2", "PCDH7"),]
write.csv(highDiabetes.markers2,file=paste0(work.dir,"highDiabetes.markers >80 high vs <20low.csv"))


### Calculating differential gene expression with BH FDR (comparing between high and mid+low groups combined)

highDiabetes.markers3 = FindMarkers(merged_seurat.beta,ident.1=colnames(merged_seurat.beta)[merged_seurat.beta$highDiabetes3=="high"],
                                    ident.2=colnames(merged_seurat.beta)[merged_seurat.beta$highDiabetes3==c("low", "mid")])
highDiabetes.markers3[c("SDF2L1","QPCT","HADH","IAPP","NPTX2","CNIH2", "PCDH7"),]
write.csv(highDiabetes.markers3,file=paste0(work.dir,"highDiabetes.markers >80 high vs mid+low.csv"))

###########################################################################################################################
#                                              Defining high_low groups                                                   #
###########################################################################################################################

highlow = ifelse(merged_seurat.beta$T2D_assoc>quants["50"],"high","low")
merged_seurat.beta$highDiabetes = highlow
unique(highlow)

merged_seurat.beta.highlow <- merged_seurat.beta
Idents(merged_seurat.beta.highlow) = merged_seurat.beta.highlow$highDiabetes

G= DimPlot(merged_seurat.beta.highlow, label = FALSE, pt.size = 2, label.size = 5)+ ggtitle("Beta cell clusters-high_low")+ theme(plot.title = element_text(hjust = 0.5))
ggsave("G- Beta cell clusters-high_low.svg", G, width = 12, height = 8, device = "svg")
ggsave("G- Beta cell clusters-high_low.tiff", G, width = 12, height = 8, device = "tiff")


### Calculating differential gene expression with BH FDR
DefaultAssay(merged_seurat.beta.highlow) = "RNA"
highDiabetes.markers = FindMarkers(merged_seurat.beta.highlow,ident.1=colnames(merged_seurat.beta.highlow)[merged_seurat.beta.highlow$highDiabetes=="high"],
                                   ident.2=colnames(merged_seurat.beta.highlow)[merged_seurat.beta.highlow$highDiabetes=="low"])
highDiabetes.markers[c("SDF2L1","QPCT","HADH","IAPP","NPTX2","CNIH2", "PCDH7"),]
write.csv(highDiabetes.markers,file=paste0(work.dir,"highDiabetes.markers >50 high.csv"))


###########################################################################################################################
#                                      Overlay genes on T2D risk                                                          #
###########################################################################################################################

T= FeaturePlot(merged_seurat.beta,features=c("DiabetesBeta","DLK1"),blend=TRUE, combine=FALSE, pt.size = 2)
ggsave("T- DLK1 overlap with disease score.svg", T, width = 12, height = 8, device = "svg")
ggsave("T- DLK1 overlap with disease score.tiff", T, width = 12, height = 8, device = "tiff")

FeaturePlot(merged_seurat.beta,features=c("DiabetesBeta","QPCT"),blend=TRUE)
FeaturePlot(merged_seurat.beta,features=c("DiabetesBeta","IAPP"),blend=TRUE)

U= FeaturePlot(merged_seurat.beta,features=c("DiabetesBeta","SDF2L1"),blend=TRUE, combine=FALSE, pt.size = 2)
ggsave("U- SDF2L1 overlap with disease score.svg", U, width = 12, height = 8, device = "svg")
ggsave("U- SDF2L1 overlap with disease score.tiff", U, width = 12, height = 8, device = "tiff")

FeaturePlot(merged_seurat.beta,features=c("DiabetesBeta","CNIH2"),blend=TRUE)
FeaturePlot(merged_seurat.beta,features=c("DiabetesBeta","NPTX2"),blend=TRUE)
FeaturePlot(merged_seurat.beta,features=c("DiabetesBeta","RGS16"),blend=TRUE)
FeaturePlot(merged_seurat.beta,features=c("DiabetesBeta","PPP1R1A"),blend=TRUE)
FeaturePlot(merged_seurat.beta,features=c("DiabetesBeta","CDKN1C"),blend=TRUE)
FeaturePlot(merged_seurat.beta,features=c("CDKN1C","BTG1"),blend=TRUE)
FeaturePlot(merged_seurat.beta,features=c("CDKN1C","KCNQ1OT1"),blend=TRUE)

L= FeaturePlot(merged_seurat.beta,features=c("SYT7","SPP1","SDF2L1", "RGS16", "PPP1R1A", "DLK1", "BTG1","CDKN1C","RPS4Y1","PTPRN2","PSAP","SIX3")) + NoLegend()
ggsave("L- markers on beta cell UMAP.svg", L, width = 12, height = 8, device = "svg")
ggsave("L- markers on beta cell UMAP.tiff", L, width = 12, height = 8, device = "tiff")

BTG1 = FeaturePlot(merged_seurat.beta,features="BTG1", pt.size = 2)+NoLegend()
ggsave("BTG1.tiff", BTG1, width = 12, height = 8, device = "tiff")

SIX3 = FeaturePlot(merged_seurat.beta,features="SIX3", pt.size = 2)+NoLegend()
ggsave("SIX3.tiff",SIX3, width = 12, height = 8, device = "tiff")

SPP1 = FeaturePlot(merged_seurat.beta,features="SPP1", pt.size = 2)+NoLegend()
ggsave("SPP1.tiff", SPP1, width = 12, height = 8, device = "tiff")

SYT7 = FeaturePlot(merged_seurat.beta,features="SYT7", pt.size = 2)+NoLegend()
ggsave("SYT7.tiff", SYT7, width = 12, height = 8, device = "tiff")

SDF2L1 = FeaturePlot(merged_seurat.beta,features="SDF2L1", pt.size = 2)+NoLegend()
ggsave("SDF2L1.tiff", SDF2L1, width = 12, height = 8, device = "tiff")

RGS16 = FeaturePlot(merged_seurat.beta,features="RGS16", pt.size = 2)+NoLegend()
ggsave("RGS16.tiff", RGS16, width = 12, height = 8, device = "tiff")

PPP1R1A = FeaturePlot(merged_seurat.beta,features="PPP1R1A", pt.size = 2)+NoLegend()
ggsave("PPP1R1A.tiff", PPP1R1A, width = 12, height = 8, device = "tiff")

CDKN1C = FeaturePlot(merged_seurat.beta,features="CDKN1C", pt.size = 2)+NoLegend()
ggsave("CDKN1C.tiff",CDKN1C, width = 12, height = 8, device = "tiff")

DLK1 = FeaturePlot(merged_seurat.beta,features="DLK1", pt.size = 2)+NoLegend()
ggsave("DLK1.tiff", DLK1, width = 12, height = 8, device = "tiff")

RPS4Y1 = FeaturePlot(merged_seurat.beta,features="RPS4Y1", pt.size = 2)+NoLegend()
ggsave("RPS4Y1.tiff", RPS4Y1, width = 12, height = 8, device = "tiff")


H= DimPlot(merged_seurat.beta, reduction = "umap", group.by = c("highDiabetes3","highDiabetes", "Sex", "Disease", "Ethnicity", "BMI", "Age"))+ NoLegend()
ggsave("H- aggregated plots.svg", H, width = 12, height = 8, device = "svg")
ggsave("H- aggregated plots.tiff", H, width = 12, height = 8, device = "tiff")

V1= DimPlot(merged_seurat.beta, reduction = "umap", pt.size = 2, group.by = "Disease")
ggsave("V1- Disease overlay on beta cells.tiff", V1, width = 12, height = 8, device = "tiff")

V2= DimPlot(merged_seurat.beta, reduction = "umap", pt.size = 2, group.by = "BMI")
ggsave("V2- BMI overlay on beta cells.tiff", V2, width = 12, height = 8, device = "tiff")

DimPlot(merged_seurat.beta, reduction = "umap", pt.size = 2, group.by = "Ethnicity")
DimPlot(merged_seurat.beta, reduction = "umap", pt.size = 2, group.by = "T2D_assoc", label = F)




save.image(file="/vol/data02/groy_data/Single_cell_RNA_seq/scRNAmerge_vs_Marselli/T2D/DEGAS_T2D_AMawla&Marselli.RData")

###########################################################################################################################
#                                              Volcano plots                                                              #
###########################################################################################################################

library(devtools)
devtools::install_github('kevinblighe/EnhancedVolcano')
library(EnhancedVolcano)
I= EnhancedVolcano(highDiabetes.markers,
                   lab = rownames(highDiabetes.markers),
                   x = 'avg_log2FC',FCcutoff = 0.5,#xlim = c(-2.5, 2.5),
                   y = 'p_val_adj',pCutoff = 5.00E-02,
                   title = 'high diabetes markers (>50 quantile=high, <50 quantile=low)',max.overlaps=15,maxoverlapsConnectors = NULL,
                   pointSize = 2.0, labSize = 4, legendPosition = 'right',
                   legendLabSize = 10,legendIconSize = 4.0,
                   drawConnectors = TRUE,typeConnectors = "closed",
                   widthConnectors = 0.5, gridlines.major = FALSE,
                   gridlines.minor = FALSE, border = 'full')
ggsave("I- Volcano: 50 by 50.svg", I, width = 12, height = 8, device = "svg")
ggsave("I- Volcano: 50 by 50.tiff", I, width = 12, height = 8, device = "tiff")

J= EnhancedVolcano(highDiabetes.markers2,
                   lab = rownames(highDiabetes.markers2),
                   x = 'avg_log2FC',FCcutoff = 0.5,xlim = c(-3.0, 3.0),
                   y = 'p_val_adj',pCutoff = 5.00E-02,
                   title = 'high diabetes markers (>80 quantile=high, <20 quantile=low)',max.overlaps=20,maxoverlapsConnectors = NULL,
                   pointSize = 2.0, labSize = 4, legendPosition = 'right',
                   legendLabSize = 10,legendIconSize = 4.0,
                   drawConnectors = TRUE,typeConnectors = "closed",
                   widthConnectors = 0.5, gridlines.major = FALSE,
                   gridlines.minor = FALSE, border = 'full')
ggsave("J- Volcano: high vs low.svg", J, width = 12, height = 8, device = "svg")
ggsave("J- Volcano: high vs low.tiff", J, width = 12, height = 8, device = "tiff")


K= EnhancedVolcano(highDiabetes.markers3,
                   lab = rownames(highDiabetes.markers3),
                   x = 'avg_log2FC',FCcutoff = 0.5,xlim = c(-3.0, 3.0),
                   y = 'p_val_adj',pCutoff = 5.00E-02,
                   title = 'high diabetes markers (>80 quantile=high, <80 quantile=mid+low)',max.overlaps=25,maxoverlapsConnectors = 25,
                   pointSize = 2.0, labSize = 4, legendPosition = 'right',
                   legendLabSize = 10,legendIconSize = 4.0,
                   drawConnectors = TRUE,typeConnectors = "closed",#lengthConnectors = unit(0.005, "npc"),
                   widthConnectors = 0.5, gridlines.major = FALSE,
                   gridlines.minor = FALSE, border = 'full')
ggsave("K- Volcano: high vs mid+low.svg", K, width = 12, height = 8, device = "svg")
ggsave("K- Volcano: high vs mid+low.tiff", K, width = 12, height = 8, device = "tiff")




###########################################################################################################################
#                                Filtering out significant genes found by DEGAS only                                      #
###########################################################################################################################

T2D_markers_betacells <- read.csv(paste0("~/data/Single_cell_RNA_seq/scRNAmerge_vs_Marselli/T2D/T2D_markers_betacells_defaultpct.csv"), sep = ",", row.names = 1)

T2D_markers_betacells = T2D_markers_betacells[T2D_markers_betacells$p_val_adj<0.05 & abs(T2D_markers_betacells$avg_log2FC)>=0.58,]

length(rownames(highDiabetes.markers2)) ## 80/20 cutoff
length(intersect(rownames(highDiabetes.markers2), rownames(T2D_markers_betacells)))
length(setdiff(rownames(highDiabetes.markers2),rownames(T2D_markers_betacells))) ## 80/20 cutoff
length(rownames(highDiabetes.markers3))
length(intersect(rownames(highDiabetes.markers3), rownames(T2D_markers_betacells)))
length(setdiff(rownames(highDiabetes.markers3),rownames(T2D_markers_betacells))) ## 80/ALL cutoff

uniqDEGASgenes_mediancutoff <- highDiabetes.markers[!rownames(highDiabetes.markers) %in% rownames(T2D_markers_betacells),]
write.csv(uniqDEGASgenes_mediancutoff,file=paste0("~/data/Single_cell_RNA_seq/scRNAmerge_vs_Marselli/T2D/uniqDEGASgenes_mediancutoff.csv"))
M= EnhancedVolcano(uniqDEGASgenes_mediancutoff,
                   lab = rownames(uniqDEGASgenes_mediancutoff),
                   x = 'avg_log2FC',FCcutoff = 0.5,xlim = c(-3.0, 3.0),
                   y = 'p_val_adj',pCutoff = 5.00E-02,
                   title = 'Unique genes found by DEGAS only (median cutoff)',max.overlaps=15,maxoverlapsConnectors = NULL,
                   pointSize = 2.0, labSize = 4, legendPosition = 'right',
                   legendLabSize = 10,legendIconSize = 4.0,
                   drawConnectors = TRUE,typeConnectors = "closed",#lengthConnectors = unit(0.005, "npc"),
                   widthConnectors = 0.5, gridlines.major = FALSE,
                   gridlines.minor = FALSE, border = 'full')
ggsave("M- Unique genes found by DEGAS only (median cutoff).svg", M, width = 12, height = 8, device = "svg")
ggsave("M- Unique genes found by DEGAS only (median cutoff).tiff", M, width = 12, height = 8, device = "tiff")

uniqDEGASgenes_80by20cutoff <- highDiabetes.markers2[!rownames(highDiabetes.markers2) %in% rownames(T2D_markers_betacells),]
write.csv(uniqDEGASgenes_80by20cutoff,file=paste0("~/data/Single_cell_RNA_seq/scRNAmerge_vs_Marselli/T2D/uniqDEGASgenes_80by20cutoff.csv"))
N= EnhancedVolcano(uniqDEGASgenes_80by20cutoff,
                   lab = rownames(uniqDEGASgenes_80by20cutoff),
                   x = 'avg_log2FC',FCcutoff = 0.58,xlim = c(-3.0, 3.0),
                   y = 'p_val_adj',pCutoff = 5.00E-02,
                   title = 'Unique genes found by DEGAS only (80/20 cutoff)',max.overlaps=15,maxoverlapsConnectors = NULL,
                   pointSize = 2.0, labSize = 4, legendPosition = 'right',
                   legendLabSize = 10,legendIconSize = 4.0, #selectLab = c("SDF2L1", "MANF"),
                   drawConnectors = TRUE,typeConnectors = "closed", #max.overlaps=20, maxoverlapsConnectors= 20,#lengthConnectors = unit(0.005, "npc")
                   widthConnectors = 0.5, gridlines.major = FALSE,
                   gridlines.minor = FALSE, border = 'full')
ggsave("N- Unique genes found by DEGAS only (80by20 cutoff).svg", N, width = 12, height = 8, device = "svg")
ggsave("N- Unique genes found by DEGAS only (80by20 cutoff).tiff", N, width = 12, height = 8, device = "tiff")
ggsave("N- Unique genes found by DEGAS only (80by20 cutoff).pdf", N, width = 12, height = 8, device = "pdf")

uniqDEGASgenes_80byallcutoff <- highDiabetes.markers3[!rownames(highDiabetes.markers3) %in% rownames(T2D_markers_betacells),]
write.csv(uniqDEGASgenes_80byallcutoff,file=paste0("~/data/Single_cell_RNA_seq/scRNAmerge_vs_Marselli/T2D/uniqDEGASgenes_80byallcutoff.csv"))
O= EnhancedVolcano(uniqDEGASgenes_80byallcutoff,
                   lab = rownames(uniqDEGASgenes_80byallcutoff),
                   x = 'avg_log2FC',FCcutoff = 0.58,xlim = c(-3.0, 3.0),
                   y = 'p_val_adj',pCutoff = 5.00E-02,
                   title = 'Unique genes found by DEGAS only (80/ALL cutoff)',max.overlaps=15,maxoverlapsConnectors = NULL,
                   pointSize = 2.0, labSize = 4, legendPosition = 'right',
                   legendLabSize = 10,legendIconSize = 4.0,
                   drawConnectors = TRUE,typeConnectors = "closed",#lengthConnectors = unit(0.005, "npc"),
                   widthConnectors = 0.5, gridlines.major = FALSE,
                   gridlines.minor = FALSE, border = 'full')
ggsave("O- Unique genes found by DEGAS only (80byALL cutoff).svg", O, width = 12, height = 8, device = "svg")
ggsave("O- Unique genes found by DEGAS only (80byALL cutoff).tiff", O, width = 12, height = 8, device = "tiff")

###########################################################################################################################
#                                         Annotate clusters                                                               #
###########################################################################################################################
# DefaultAssay(merged_seurat) = "RNA"
FeaturePlot(merged_seurat, features = c("INS", "GCG", "SST", "PPY", "GHRL", "KRT19", "SOX9", "PRSS1", "COL6A1", "CD44", "INHBA", "KIT", "CD69"), min.cutoff = "q9")
FeaturePlot(merged_seurat, features = "PLVAP")

merged.fivedatasets <- merged_seurat
Idents(merged.fivedatasets) <- "seurat_clusters"
merged.fivedatasets <- RenameIdents(merged.fivedatasets, `0` = "Alpha cells", `1` = "Alpha cells", `2` = "Alpha cells",
                                    `3` = "Beta cells", `4` = "Acinar cells", `5` = "Delta cells", `6` = "Ductal cells", `7` = "Beta cells", `8` = "Gamma+Ghrelin cells", `9` = "Stellate+Mesenchymal cells",
                                    `10` = "Mast+Endothelial cells", `11` = "Ductal+Acinar cells", `12` = "Ductal cells", `13` = "Beta cells", `14` = "Alpha cells", `15` = "Stellate+Mesenchymal cells", `16` = "Mast cells",
                                    `17` = "Stellate+Mesenchymal cells", `18` = "Alpha cells", `19` = "Beta cells", `20` = "Beta cells", `21`="Beta cells", `22`=" ")

DimPlot(merged.fivedatasets, label = TRUE)

S= VlnPlot(merged.fivedatasets, features = c("T2D_assoc"))
ggsave("S- Violin plot of T2D association.svg", S, width = 12, height = 8, device = "svg")
ggsave("S- Violin plot of T2D association.tiff", S, width = 12, height = 8, device = "tiff")


###########################################################################################################################
#                                                 Gene Ontology/Pathway Analysis                                          #
###########################################################################################################################

library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(cowplot)
library(stringr)

# highdiabetesgenes <- highDiabetes.markers2
# highdiabetesgenes$entrez <- mapIds(org.Hs.eg.db, rownames(highDiabetes.markers), 'ENTREZID', 'SYMBOL', multiVals="first")
up_cutoff <- 0.58
down_cutoff <- -0.58

uniqDEGASgenes_80by20cutoff_up <- uniqDEGASgenes_80by20cutoff[uniqDEGASgenes_80by20cutoff$avg_log2FC>up_cutoff & uniqDEGASgenes_80by20cutoff$p_val_adj<0.05,]
uniqDEGASgenes_80by20cutoff_down <- uniqDEGASgenes_80by20cutoff[uniqDEGASgenes_80by20cutoff$avg_log2FC<down_cutoff & uniqDEGASgenes_80by20cutoff$p_val_adj<0.05,]

q1 <- enrichGO(rownames(uniqDEGASgenes_80by20cutoff_up), OrgDb = org.Hs.eg.db, pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05,
               keyType = "SYMBOL", ont="BP") ##if using SYMBOL do not use setreadable()
q1 <- pairwise_termsim(q1)
q1 <- simplify(q1, cutoff=0.5, by="p.adjust", select_fun=min)
q1 <- dotplot(q1, label_format = 10, showCategory=20, orderBy="GeneRatio") + ggtitle("Upregulated-High Diabetes markers (high/low)") + theme(panel.grid=element_blank())

q2 <- enrichGO(rownames(uniqDEGASgenes_80by20cutoff_down), OrgDb = org.Hs.eg.db, pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05,
               keyType = "SYMBOL", ont="BP") ##if using SYMBOL do not use setreadable()
q2 <- pairwise_termsim(q2)
q2 <- simplify(q2, cutoff=0.5, by="p.adjust", select_fun=min)
q2 <- dotplot(q2, label_format = 10, showCategory=20, orderBy="GeneRatio") +ggtitle("Downregulated-High Diabetes markers (high/low)")+ theme(panel.grid=element_blank(),axis.text=element_text(size=16))
theme(axis.text=element_text(size=16),axis.title=element_text(size=12,face="bold"))
#Q= plot_grid(q1, q2, ncol=2, labels=LETTERS[1:2])
#ggsave("Q- GO Biological Process-Unique DEGAS genes (highbylow).svg", Q, width = 30, height = 8, device = "svg")
#ggsave("Q- GO Biological Process-Unique DEGAS genes (highbylow).tiff", Q, width = 30, height = 8, device = "tiff")

ggsave("q1- GO Biological Process-Unique DEGAS genes (highbylow).pdf", q1, width = 12, height = 8, device = "pdf")
write_csv(q1$data, "GO-BP-Upregulated-High Diabetes markers (highbylow).csv")

ggsave("q2- GO Biological Process-Unique DEGAS genes (highbylow).pdf", q2, width = 12, height = 8, device = "pdf")
write.csv(q2$data, "GO-BP-Downregulated-High Diabetes markers (highbylow).csv")



uniqDEGASgenes_80byallcutoff_up <- uniqDEGASgenes_80byallcutoff[uniqDEGASgenes_80byallcutoff$avg_log2FC>up_cutoff & uniqDEGASgenes_80byallcutoff$p_val_adj<0.05,]
uniqDEGASgenes_80byallcutoff_down <- uniqDEGASgenes_80byallcutoff[uniqDEGASgenes_80byallcutoff$avg_log2FC>down_cutoff & uniqDEGASgenes_80byallcutoff$p_val_adj<0.05,]

r1 <- enrichGO(rownames(uniqDEGASgenes_80byallcutoff_up), OrgDb = org.Hs.eg.db, pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05,
               keyType = "SYMBOL", ont="BP") ##if using SYMBOL do not use setreadable()
r1 <- pairwise_termsim(r1)
r1 <- simplify(r1, cutoff=0.5, by="p.adjust", select_fun=min)
r1 <- dotplot(r1, label_format = 10, showCategory=20, orderBy="GeneRatio") +ggtitle("Upregulated-High Diabetes markers (high/mid+low)")

r2 <- enrichGO(rownames(uniqDEGASgenes_80byallcutoff_down), OrgDb = org.Hs.eg.db, pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05,
               keyType = "SYMBOL", ont="BP") ##if using SYMBOL do not use setreadable()
r2 <- pairwise_termsim(r2)
r2 <- simplify(r2, cutoff=0.5, by="p.adjust", select_fun=min)
r2 <- dotplot(r2, label_format = 10, showCategory=20, orderBy="GeneRatio") +ggtitle("Downregulated-High Diabetes markers (high/mid+low)")

#R= plot_grid(r1, r2, ncol=2, labels=LETTERS[1:2])
ggsave("r1- GO Biological Process-Unique DEGAS genes (highbyALL).pdf", r1, width = 12, height = 8, device = "pdf")
ggsave("r2- GO Biological Process-Unique DEGAS genes (highbyALL).pdf", r2, width = 12, height = 8, device = "pdf")


###########################################################################################################################
#                                                     GSEA                                                                #
###########################################################################################################################
highDiabetes.markers2.noFCthreshold = FindMarkers(merged_seurat.beta,ident.1=colnames(merged_seurat.beta)[merged_seurat.beta$highDiabetes3=="high"],
                                                  ident.2=colnames(merged_seurat.beta)[merged_seurat.beta$highDiabetes3=="low"], logfc.threshold = 0.0)
P= EnhancedVolcano(highDiabetes.markers2.noFCthreshold,
                   lab = rownames(highDiabetes.markers2.noFCthreshold),
                   x = 'avg_log2FC',FCcutoff = 0.5,xlim = c(-3.0, 3.0),
                   y = 'p_val_adj',pCutoff = 5.00E-02,
                   title = 'high diabetes markers (>80 quantile=high, <20 quantile=low) w/o FC threshold',max.overlaps=10,maxoverlapsConnectors = NULL,
                   pointSize = 2.0, labSize = 4, legendPosition = 'right',
                   legendLabSize = 10,legendIconSize = 4.0,
                   drawConnectors = TRUE,typeConnectors = "closed",
                   widthConnectors = 0.5, gridlines.major = FALSE,
                   gridlines.minor = FALSE, border = 'full')
ggsave("P- Volcano-high diabetes markers (>80 quantile=high, <20 quantile=low) w-o FC threshold.svg", P, width = 12, height = 8, device = "svg")
ggsave("P- Volcano-high diabetes markers (>80 quantile=high, <20 quantile=low) w-o FC threshold.tiff", P, width = 12, height = 8, device = "tiff")

library(ggplot2)
library(dplyr)
library(hrbrthemes)

diabetic_GSEA <- read.csv(paste0("/vol/data02/groy_data/Single_cell_RNA_seq/scRNAmerge_vs_Marselli/T2D/diabetic_GSEA.csv"), sep = ",")
diabetic_GSEA <- diabetic_GSEA[, c(1,5,7)]

# Reorder the data
diabetic_GSEA_sort <- diabetic_GSEA %>%
  arrange(NES) %>%
  mutate(GS=factor(GS,GS))

FDR_cutoff = 0.1

ggplot(diabetic_GSEA_sort, aes(x=GS, y=NES)) +
  geom_segment(
    aes(x=GS, xend=GS, y=0, yend=NES), 
    color=ifelse(diabetic_GSEA_sort$FDR.q.val<0.1, "orange", "grey"), 
    size=1.5
  ) +coord_flip()+
  geom_point(
    color=ifelse(diabetic_GSEA_sort$FDR.q.val<0.1, "orange", "grey"), 
    size=ifelse(diabetic_GSEA_sort$GS %in% FDR_cutoff, 5, 2)
  ) +
  theme_ipsum() +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    legend.position="none"
  ) +
  xlab("") +
  ylab("NES") +
  ggtitle("GSEA high vs low T2D")+
  geom_hline(yintercept = -4, linetype = "solid", color = "black") +  
  geom_vline(xintercept = 0.4, linetype = "solid", color = "black")


save.image(file="/vol/data02/groy_data/Single_cell_RNA_seq/scRNAmerge_vs_Marselli/T2D/DEGAS_T2D_AMawla&Marselli.RData")








