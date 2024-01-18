##############################################################
# The scRNAseq datasets are provided by A Mawla
# Setting up R environment
library(biomaRt)
library(dplyr)
library(EnhancedVolcano)
library(Seurat)
library(DEGAS)
library(Rtsne)
library(ggplot2)
library(pROC)
library(tidyverse)
library(doParallel) # Added for parrallel marker selection TSJ 20230401
n_cores = 4 # Added for parallel marker selection TSJ 20230401
registerDoParallel(cores=n_cores) # Added for parrallel marker selection TSJ 20230401
work.dir = "/vol/data02/groy_data/Single_cell_RNA_seq/scRNAmerge_vs_Marselli/BMI/"


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
bulk.expr = bulk.expr[, -83]

###########################################################################################################################
#                                          Loading Bulk metadata                                                          #
###########################################################################################################################

bulk.meta <-  read.csv(paste0("/vol/data02/groy_data/Single_cell_RNA_seq/finaltest/marselli_GSE159984_metadata_clean.csv"), sep = ",", header = T, row.names = "SAMN_and_disease")
bulk.meta <- bulk.meta[-83,]

#Check
setdiff(rownames(bulk.meta),colnames(bulk.expr))
identical(rownames(bulk.meta),colnames(bulk.expr))
all.equal(rownames(bulk.meta),colnames(bulk.expr))


###########################################################################################################################
#                                                       Loading SC data                                                   #
###########################################################################################################################

merged_seurat <- readRDS("~/data/Single_cell_RNA_seq/finaltest/merged_seurat_AMawla_NEW.rds") 
sc.expr = as.data.frame(merged_seurat@assays$RNA@counts)


###Tabulate cells by cluster ID
library(data.table)
md <- merged_seurat@meta.data %>% as.data.table #ONE ROW PER CELL

## count the number of cells per unique combinations of "Sample" and "seurat_clusters"
md2 <- md[, .N, by = c("orig.ident", "seurat_clusters")]

## with additional casting after the counting
#md2 <- md[, .N, by = c("orig.ident", "seurat_clusters")] %>% dcast(., orig.ident ~ seurat_clusters, value.var = "N")
write.csv(md2,file=paste0(work.dir,"beta cells per cluster.csv"),quote=FALSE)

md3 = md2 %>% group_by(seurat_clusters) %>% summarize(SumAvg = sum(N, na.rm = TRUE))


###########################################################################################################################
#                         Identifying final feature set (Changed this section TSJ 20230401)                               #
###########################################################################################################################
if(TRUE){
  # In this section we use the new parrallelized find cluster markers with
  # 800 top features for each cluster and patient group.
  int.feats = intersect(rownames(bulk.expr),rownames(sc.expr)) # how many is this?
  sc.markers = FindClusterMarkersPar(t(sc.expr[int.feats,]),merged_seurat$seurat_clusters,800) 
  write.csv(sc.markers, file=paste0(work.dir,"sc.clustermarkers.csv"))
  bulk.markers = FindClusterMarkersPar(t(bulk.expr[int.feats,]),bulk.meta$BMI,400) 
  write.csv(bulk.markers, file=paste0(work.dir,"bulk.clustermarkers.csv"))
  highvar.features = intersect(bulk.markers$Gene,sc.markers$Gene)
  highvar.features = intersect(unique(bulk.markers$Gene),unique(sc.markers$Gene))
  length(highvar.features)
  saveRDS(highvar.features,file=paste0(work.dir,"featureSet_bmi.rds"))
}

###########################################################################################################################
#                                Preprocessing data (log, normalization, scale)                                           #
###########################################################################################################################

scDat = preprocessCounts(sc.expr[highvar.features,]) #matrix of expression values from scRNAseq data
scLab = toOneHot(merged_seurat$seurat_clusters)
patDat = preprocessCounts(bulk.expr[highvar.features,]) #matrix of expression values from bulkRNAseq data, columns should be in same order as scDat
#patLab = toOneHot(bulk.meta[colnames(bulk.expr),"BMI"])
bmi_tmp = bulk.meta[colnames(bulk.expr),"BMI"]
obese_status = ifelse(bmi_tmp>30,"obese","other")
obese_status[obese_status=="other"] = ifelse(bmi_tmp[obese_status=="other"]>25,"overweight","healthyweight")
patLab = toOneHot(obese_status)



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
# AUC in cells

###########################################################################################################################
#                                                   Post processing                                                       #
###########################################################################################################################

boxplot(scpatPreds[,"overweight"]~fromOneHot(scLab))
boxplot(scpatPreds[,"obese"]~fromOneHot(scLab))

BMIabove30_assoc = knnSmooth(toCorrCoeff(scpatPreds[,"obese"]),merged_seurat@reductions$umap@cell.embeddings) 
merged_seurat$BMIabove30_assoc = BMIabove30_assoc

BMIbelow25_assoc = knnSmooth(toCorrCoeff(scpatPreds[,"healthyweight"]),merged_seurat@reductions$umap@cell.embeddings)
merged_seurat$BMIbelow25_assoc = BMIbelow25_assoc

BMI25to30_assoc = knnSmooth(toCorrCoeff(scpatPreds[,"overweight"]),merged_seurat@reductions$umap@cell.embeddings)
merged_seurat$BMI25to30_assoc = BMI25to30_assoc

FeaturePlot(merged_seurat, features = c("BMIabove30_assoc", "BMIbelow25_assoc", "BMI25to30_assoc", "INS", "GCG", "SST", "PPY", "GHRL"), min.cutoff = "q9")
A= FeaturePlot(merged_seurat,features="BMIabove30_assoc", pt.size = 2) + scale_color_gradient2(low = "black",mid="lavender",high="red")
ggsave("A- BMIabove30_assoc.svg", A, width = 12, height = 8, device = "svg")
ggsave("A- BMIabove30_assoc.tiff", A, width = 12, height = 8, device = "tiff")
ggsave("A- BMIabove30_assoc.pdf", A, width = 12, height = 8, device = "pdf")

B= FeaturePlot(merged_seurat,features="BMIbelow25_assoc", pt.size = 2) + scale_color_gradient2(low = "black",mid="lavender",high="red")
ggsave("B- BMIbelow25_assoc.svg", B, width = 12, height = 8, device = "svg")
ggsave("B- BMIbelow25_assoc.tiff", B, width = 12, height = 8, device = "tiff")
ggsave("B- BMIbelow25_assoc.pdf", B, width = 12, height = 8, device = "pdf")

C= FeaturePlot(merged_seurat,features="BMI25to30_assoc", pt.size = 2) + scale_color_gradient2(low = "black",mid="lavender",high="red")
ggsave("C- BMI25to30_assoc.svg", C, width = 12, height = 8, device = "svg")
ggsave("C- BMI25to30_assoc.tiff", C, width = 12, height = 8, device = "tiff")
ggsave("C- BMI25to30_assoc.pdf", C, width = 12, height = 8, device = "pdf")

FeaturePlot(merged_seurat,features="INS", label = T) + scale_color_gradient2(low = "black",mid="lavender",high="red")
FeaturePlot(merged_seurat,features="GCG", label = T) + scale_color_gradient2(low = "black",mid="lavender",high="red")
FeaturePlot(merged_seurat,features="SST", label = T) + scale_color_gradient2(low = "black",mid="lavender",high="red")
FeaturePlot(merged_seurat,features="PPY", label = T) + scale_color_gradient2(low = "black",mid="lavender",high="red")
FeaturePlot(merged_seurat,features="GHRL", label = T) + scale_color_gradient2(low = "black",mid="lavender",high="red")
boxplot(merged_seurat$T2D_assoc~merged_seurat$seurat_clusters)
DimPlot(merged_seurat, group.by = "seurat_clusters", label=T)


###########################################################################################################################
#                                                      Subsetting beta cells                                              #
###########################################################################################################################

FeaturePlot(merged_seurat, features = "INS")
DimPlot(merged_seurat, group.by = "seurat_clusters", label=T)+ theme(plot.title = element_text(hjust = 0.5))
boxplot(merged_seurat@assays$RNA@counts["INS",]~merged_seurat$seurat_clusters)
merged_seurat.beta <- subset(merged_seurat, subset = seurat_clusters %in% c(3,7,13,19,20,21))
DefaultAssay(merged_seurat.beta) = "integrated"
merged_seurat.beta <- ScaleData(merged_seurat.beta,features=rownames(merged_seurat.beta))
merged_seurat.beta <- RunPCA(merged_seurat.beta, features = VariableFeatures(object = merged_seurat.beta))
merged_seurat.beta <- RunUMAP(merged_seurat.beta, dims=1:10)
merged_seurat.beta <- FindNeighbors(merged_seurat.beta,reduction="umap",dims=1:2)
merged_seurat.beta <- FindClusters(merged_seurat.beta)
merged_seurat.beta$obeseBeta = knnSmooth(merged_seurat.beta$BMIabove30_assoc,merged_seurat.beta@reductions$umap@cell.embeddings[,c("UMAP_1","UMAP_2")],20)
merged_seurat.beta$overweightBeta = knnSmooth(merged_seurat.beta$BMI25to30_assoc,merged_seurat.beta@reductions$umap@cell.embeddings[,c("UMAP_1","UMAP_2")],20)
merged_seurat.beta$healthyBeta = knnSmooth(merged_seurat.beta$BMIbelow25_assoc,merged_seurat.beta@reductions$umap@cell.embeddings[,c("UMAP_1","UMAP_2")],20)



DimPlot(merged_seurat.beta, label = F, pt.size = 1, label.size = 5)+ ggtitle("Beta cell clusters")+ theme(plot.title = element_text(hjust = 0.5))
DimPlot(merged_seurat.beta, reduction = "umap", group.by = "BMI", pt.size = 2)

D= FeaturePlot(merged_seurat.beta,features="obeseBeta", pt.size = 2, label = T)+ ggtitle("Beta cells associated with >30 BMI (Obese)")+ theme(plot.title = element_text(hjust = 0.5))+ scale_color_gradient2(low = "black",mid="lavender",high="red")
ggsave("D- Beta cells associated with >30 BMI (Obese).svg", D, width = 12, height = 8, device = "svg")
ggsave("D- Beta cells associated with >30 BMI (Obese).tiff", D, width = 12, height = 8, device = "tiff")
ggsave("D- Beta cells associated with >30 BMI (Obese).pdf", D, width = 12, height = 8, device = "pdf")

E= FeaturePlot(merged_seurat.beta,features="overweightBeta", pt.size = 2, label = T)+ ggtitle("Beta cells associated with 25-30 BMI (Overweight)")+ theme(plot.title = element_text(hjust = 0.5))+ scale_color_gradient2(low = "black",mid="lavender",high="red")
ggsave("E- Beta cells associated with 25-30 BMI (Overweight).svg", E, width = 12, height = 8, device = "svg")
ggsave("E- Beta cells associated with 25-30 BMI (Overweight).tiff", E, width = 12, height = 8, device = "tiff")
ggsave("E- Beta cells associated with 25-30 BMI (Overweight).pdf", E, width = 12, height = 8, device = "pdf")

F= FeaturePlot(merged_seurat.beta,features="healthyBeta", pt.size = 2, label = T)+ ggtitle("Beta cells associated with <25 BMI (Healthyweight)")+ theme(plot.title = element_text(hjust = 0.5))+ scale_color_gradient2(low = "black",mid="lavender",high="red")
ggsave("F- Beta cells associated with <25 BMI (Healthyweight).svg", F, width = 12, height = 8, device = "svg")
ggsave("F- Beta cells associated with <25 BMI (Healthyweight).tiff", F, width = 12, height = 8, device = "tiff")
ggsave("F- Beta cells associated with <25 BMI (Healthyweight).pdf", F, width = 12, height = 8, device = "pdf")

FeaturePlot(merged_seurat.beta,features=c("healthyBeta","overweightBeta", "obeseBeta"), label = T)

betacluster.markers = FindAllMarkers(merged_seurat.beta)
write.csv(betacluster.markers,file=paste0(work.dir,"betacellcluster_markers.csv"),quote=FALSE)


###Tabulate cells by cluster ID
library(data.table)
md <- merged_seurat.beta@meta.data %>% as.data.table #ONE ROW PER CELL

## count the number of cells per unique combinations of "Sample" and "seurat_clusters"
md2 <- md[, .N, by = c("orig.ident", "seurat_clusters")]

## with additional casting after the counting
#md2 <- md[, .N, by = c("orig.ident", "seurat_clusters")] %>% dcast(., orig.ident ~ seurat_clusters, value.var = "N")
#write.csv(md2,file=paste0(work.dir,"beta cells per cluster.csv"),quote=FALSE)

md3 = md2 %>% group_by(seurat_clusters) %>% summarize(SumAvg = sum(N, na.rm = TRUE))




###########################################################################################################################
#                                              DEGs from Obese group                                                      #
###########################################################################################################################

merged_seurat.betaOBESE <- merged_seurat.beta
merged_seurat.betaOBESE$obesegroup_BMI = ifelse(merged_seurat.betaOBESE$obeseBeta>median(merged_seurat.betaOBESE$obeseBeta),"obesegroup_higherBMI","obesegroup_lowerBMI")

G= DimPlot(merged_seurat.betaOBESE, group.by = "obesegroup_BMI", pt.size = 2) + ggtitle("Applied median threshold on beta cells from Obese group")
ggsave("G- Applied median threshold on beta cells from Obese group.svg", G, width = 12, height = 8, device = "svg")
ggsave("G- Applied median threshold on beta cells from Obese group.tiff", G, width = 12, height = 8, device = "tiff")

# Calculating differential gene expression with BH FDR
DefaultAssay(merged_seurat.betaOBESE) = "RNA"
merged_seurat.betaOBESE = NormalizeData(merged_seurat.betaOBESE,normalization.method = "LogNormalize")

higherBMImarkersinObesity = FindMarkers(merged_seurat.betaOBESE,ident.1=colnames(merged_seurat.betaOBESE)[merged_seurat.betaOBESE$obesegroup_BMI=="obesegroup_higherBMI"],
                                        ident.2=colnames(merged_seurat.betaOBESE)[merged_seurat.betaOBESE$obesegroup_BMI=="obesegroup_lowerBMI"])

higherBMImarkersinObesity[c("SDF2L1","QPCT","HADH","IAPP","NPTX2","CNIH2", "PCDH7"),]
H= EnhancedVolcano(higherBMImarkersinObesity,
                   lab = rownames(higherBMImarkersinObesity),
                   x = 'avg_log2FC',FCcutoff = 0.5,xlim = c(-5.0, 5.0),
                   y = 'p_val_adj',pCutoff = 5.00E-02,
                   title = 'Markers of Higher BMI in Obese group',
                   pointSize = 3.0, labSize = 6, legendPosition = 'right',
                   legendLabSize = 10,legendIconSize = 4.0,
                   drawConnectors = TRUE,typeConnectors = "closed",
                   widthConnectors = 0.5, gridlines.major = FALSE,
                   gridlines.minor = FALSE, border = 'full')
ggsave("H- Markers of Higher BMI in Obese group.svg", H, width = 12, height = 8, device = "svg")
ggsave("H- Markers of Higher BMI in Obese group.tiff", H, width = 12, height = 8, device = "tiff")

write.csv(deg_table_obese,file=paste0("/vol/data02/groy_data/Single_cell_RNA_seq/scRNAmerge_vs_Marselli/BMI/deg_table_obese.csv"))
write.csv(higherBMImarkersinObesity,file=paste0("/vol/data02/groy_data/Single_cell_RNA_seq/scRNAmerge_vs_Marselli/BMI/higherBMImarkersinObesity.csv"))


###########################################################################################################################
#                                           DEGs from ND obese vs T2D obese                                               #
###########################################################################################################################

## Comparing cluster 1,8,12 to clusters 7 & 10
cluster1812markersvs710 <- FindMarkers(merged_seurat.betaOBESE, ident.1 = c(1,8,12), ident.2 = c(7,10), min.pct = 0.25) ## compare cluster 1,8,12 against 7,10
I= EnhancedVolcano(cluster1812markersvs710,
                   lab = rownames(cluster1812markersvs710),
                   x = 'avg_log2FC',FCcutoff = 0.5,xlim = c(-5.0, 5.0),
                   y = 'p_val_adj',pCutoff = 5.00E-02,
                   title = 'Markers of Higher BMI in Obese group clusters 1,8,12 vs 7,10',
                   pointSize = 3.0, labSize = 6, legendPosition = 'right',
                   legendLabSize = 10,legendIconSize = 4.0,
                   drawConnectors = TRUE,typeConnectors = "closed",
                   widthConnectors = 0.5, gridlines.major = FALSE,
                   gridlines.minor = FALSE, border = 'full')
ggsave("I- Volcano- Markers of Higher BMI in Obese group clusters 1,8,12 vs 7,10.svg", I, width = 12, height = 8, device = "svg")
ggsave("I- Volcano- Markers of Higher BMI in Obese group clusters 1,8,12 vs 7,10.tiff", I, width = 12, height = 8, device = "tiff")

# write.csv(cluster8markersvs23,file=paste0("/vol/data02/groy_data/Single_cell_RNA_seq/realigned_bulk_microarray_dropped/BMI/cluster8markersvs23.csv"))

merged_seurat.obeseBMI <- subset(merged_seurat.betaOBESE, subset = obeseBeta>0.1)
DimPlot(merged_seurat.obeseBMI, label = FALSE, pt.size = 2)+ggtitle("Subset of β-cells with high obesity-association scores")

md4 <- merged_seurat.obeseBMI@meta.data %>% as.data.table #ONE ROW PER CELL
## count the number of cells per unique combinations of "Sample" and "seurat_clusters"
md5 <- md4[, .N, by = c("orig.ident", "seurat_clusters")]
md6 = md5 %>% group_by(seurat_clusters) %>% summarize(SumAvg = sum(N, na.rm = TRUE))



###########################################################################################################################
#                                         Annotate obese clusters                                                         #
###########################################################################################################################

DefaultAssay(merged_seurat.obeseBMI) = "RNA"
merged_seurat.obeseBMI <- RenameIdents(merged_seurat.obeseBMI, `1`="T2D obese", `8`="T2D obese", `12`="T2D obese", `4`="ND obese",`7`="ND obese", `10`="ND obese" )

J= DimPlot(merged_seurat.obeseBMI, label = FALSE, pt.size = 2)
ggsave("J- ND obese and T2D obese.svg", J, width = 12, height = 8, device = "svg")
ggsave("J- ND obese and T2D obese.tiff", J, width = 12, height = 8, device = "tiff")
ggsave("J- ND obese and T2D obese.pdf", J, width = 12, height = 8, device = "pdf")

NDobese_vs_T2Dobese <- FindMarkers(merged_seurat.obeseBMI, ident.1 = "ND obese", ident.2 = "T2D obese", min.pct = 0.25)
K= EnhancedVolcano(NDobese_vs_T2Dobese,
                   lab = rownames(NDobese_vs_T2Dobese),
                   x = 'avg_log2FC',FCcutoff = 0.5,xlim = c(-5.0, 5.0),
                   y = 'p_val_adj',pCutoff = 5.00E-02,
                   title = 'Markers- ND obese vs T2D obese',
                   pointSize = 3.0, labSize = 6, legendPosition = 'right',
                   legendLabSize = 10,legendIconSize = 4.0,
                   drawConnectors = TRUE,typeConnectors = "closed",
                   widthConnectors = 0.5, gridlines.major = FALSE,
                   gridlines.minor = FALSE, border = 'full')
ggsave("K- Volcano-Markers- ND obese vs T2D obese.svg", K, width = 20, height = 20, device = "svg")
ggsave("K- Volcano- Markers- ND obese vs T2D obese.tiff", K, width = 18, height = 15, device = "tiff")

NDobese_vs_T2Dobese[c("SDF2L1","QPCT","HADH","IAPP","NPTX2","CNIH2", "PCDH7", "MANF","ATF4","DDIT3","HSPA5"),]
write.csv(NDobese_vs_T2Dobese,file=paste0("/vol/data02/groy_data/Single_cell_RNA_seq/finaltest/NDobese_vs_T2Dobese.csv"))


##############################################################
# Annotate clusters
merged.fivedatasets_for_annotation <- merged_seurat
Idents(merged.fivedatasets_for_annotation) <- "seurat_clusters"
merged.fivedatasets_for_annotation <- RenameIdents(merged.fivedatasets_for_annotation, `0` = "Alpha cells", `1` = "Alpha cells", `2` = "Alpha cells",
                                    `3` = "Beta cells", `4` = "Acinar cells", `5` = "Delta cells", `6` = "Ductal cells", `7` = "Beta cells", `8` = "Gamma+Ghrelin cells", `9` = "Stellate+Mesenchymal cells",
                                    `10` = "Mast+Endothelial cells", `11` = "Ductal+Acinar cells", `12` = "Ductal cells", `13` = "Beta cells", `14` = "Alpha cells", `15` = "Stellate+Mesenchymal cells", `16` = "Mast cells",
                                    `17` = "Stellate+Mesenchymal cells", `18` = "Alpha cells", `19` = "Beta cells", `20` = "Beta cells", `21`="Beta cells", `22`=" ")

DimPlot(merged.fivedatasets_for_annotation, label = TRUE)

M= VlnPlot(merged.fivedatasets_for_annotation, features = c("BMIabove30_assoc"))
ggsave("M- Violin plot of BMIabove30 association.svg", M, width = 12, height = 8, device = "svg")
ggsave("M- Violin plot of BMIabove30 association.tiff", M, width = 12, height = 8, device = "tiff")

###########################################################################################################################
#                                                 Gene Ontology/Pathway Analysis                                          #
###########################################################################################################################

library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(org.Hs.eg.db)
library(AnnotationDbi)

up_cutoff <- 0.58
down_cutoff <- -0.58

NDobese_vs_T2Dobese_up <- NDobese_vs_T2Dobese[NDobese_vs_T2Dobese$avg_log2FC>up_cutoff & NDobese_vs_T2Dobese$p_val_adj<0.05,]
NDobese_vs_T2Dobese_down <- NDobese_vs_T2Dobese[NDobese_vs_T2Dobese$avg_log2FC<down_cutoff & NDobese_vs_T2Dobese$p_val_adj<0.05,]

l1 <- enrichGO(rownames(NDobese_vs_T2Dobese_up), OrgDb = org.Hs.eg.db,
               keyType = "SYMBOL",ont="BP",  pAdjustMethod = "BH",pvalueCutoff=0.01, qvalueCutoff=0.05) ##if using SYMBOL do not use setreadable()
l1 <- pairwise_termsim(l1)
l1 <- simplify(l1, cutoff=0.5, by="p.adjust", select_fun=min)
l1 <- dotplot(l1, label_format = 10, showCategory=20, orderBy="GeneRatio") +ggtitle("Upregulated-ND Beta cells with high BMI-association score")+ theme(panel.grid=element_blank())

l2 <- enrichGO(rownames(NDobese_vs_T2Dobese_down), OrgDb = org.Hs.eg.db,
               keyType = "SYMBOL",ont="BP",  pAdjustMethod = "BH",pvalueCutoff=0.01, qvalueCutoff=0.05) ##if using SYMBOL do not use setreadable()
l2 <- pairwise_termsim(l2)
l2 <- simplify(l2, cutoff=0.5, by="p.adjust", select_fun=min)
l2 <- dotplot(l2, label_format = 10, showCategory=20, orderBy="GeneRatio") +ggtitle("Upregulated-T2D Beta cells with high BMI-association score")+ theme(panel.grid=element_blank())

ggsave("l1- GO BP-NDobese_vs_T2Dobese.pdf", l1, width = 12, height = 8, device = "pdf")
write.csv(l1$data, "GO-BP-Upregulated-ND Beta cells with high BMI-association score.csv")

ggsave("l2- GO BP-NDobese_vs_T2Dobese.pdf", l2, width = 12, height = 8, device = "pdf")
write.csv(l2$data, "GO-BP-Upregulated-T2D Beta cells with high BMI-association score.csv")

###########################################################################################################################
#                                                 GSEA                                                                    #
###########################################################################################################################

library(ggplot2)
library(dplyr)
library(hrbrthemes)

obesity_GSEA <- read.csv(paste0("/vol/data02/groy_data/Single_cell_RNA_seq/scRNAmerge_vs_Marselli/BMI/obesity-GSEA.csv"), sep = ",", header = T)
obesity_GSEA <- obesity_GSEA[, c(1,5,7)]

# Reorder the data
obesity_GSEA_sort <- obesity_GSEA %>%
  arrange(NES) %>%
  mutate(GS=factor(GS,GS))

ggplot(obesity_GSEA_sort, aes(x=GS, y=NES)) +
  geom_segment(
    aes(x=GS, xend=GS, y=0, yend=NES), 
    color=ifelse(obesity_GSEA_sort$FDR.q.val<0.09, "orange", "grey"), 
    size=1.5
  ) +coord_flip()+
  geom_point(
    color=ifelse(obesity_GSEA_sort$FDR.q.val<0.1, "orange", "grey"), 
    size=ifelse(obesity_GSEA_sort$GS %in% FDR_cutoff, 5, 2)
  ) +
  theme_ipsum() +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    legend.position="none"
  ) +
  xlab("") +
  ylab("NES") +
  ggtitle("GSEA ND vs T2D")+
  geom_hline(yintercept = -3.5, linetype = "solid", color = "black") +  
  geom_vline(xintercept = 0.4, linetype = "solid", color = "black")


###########################################################################################################################
#                                      Overlay genes on T2D risk                                                          #
###########################################################################################################################

N= FeaturePlot(merged_seurat.beta,features=c("obeseBeta","DLK1"),blend=TRUE, combine=FALSE, pt.size = 2)
ggsave("N- DLK1 overlap with OBESITY score.svg", N, width = 12, height = 8, device = "svg")
ggsave("N- DLK1 overlap with OBESITY score.tiff", N, width = 12, height = 8, device = "tiff")
ggsave("N- DLK1 overlap with OBESITY score.pdf", N, width = 12, height = 8, device = "pdf")

FeaturePlot(merged_seurat.beta,features=c("obeseBeta","QPCT"),blend=TRUE)
FeaturePlot(merged_seurat.beta,features=c("obeseBeta","IAPP"),blend=TRUE)

FeaturePlot(merged_seurat.obeseBMI,features=c("obeseBeta","SDF2L1"),blend=TRUE, combine=FALSE, pt.size = 2)
ggsave("O- SDF2L1 overlap with OBESITY score.svg", O, width = 12, height = 8, device = "svg")
ggsave("O- SDF2L1 overlap with OBESITY score.tiff", O, width = 12, height = 8, device = "tiff")

FeaturePlot(merged_seurat.obeseBMI,features=c("obeseBeta","MANF"),blend=TRUE, combine=FALSE, pt.size = 2)

FeaturePlot(merged_seurat.beta,features=c("obeseBeta","DDIT3"),blend=TRUE)
FeaturePlot(merged_seurat.beta,features=c("obeseBeta","NPTX2"),blend=TRUE)
FeaturePlot(merged_seurat.beta,features=c("obeseBeta","RGS16"),blend=TRUE)
FeaturePlot(merged_seurat.beta,features=c("obeseBeta","PPP1R15A"),blend=TRUE, combine=FALSE, pt.size = 2)
UCN3= FeaturePlot(merged_seurat.beta,features="UCN3", pt.size = 2) + NoLegend()
ggsave("UCN3.tiff", UCN3, width = 12, height = 8, device = "tiff")

NEUROD1= FeaturePlot(merged_seurat.beta,features="NEUROD1", pt.size = 2) + NoLegend()
ggsave("NEUROD1.tiff", NEUROD1, width = 12, height = 8, device = "tiff")

G6PC2= FeaturePlot(merged_seurat.beta,features="G6PC2", pt.size = 2) + NoLegend()
ggsave("G6PC2.tiff", G6PC2, width = 12, height = 8, device = "tiff")

DDIT3= FeaturePlot(merged_seurat.betatest,features="DDIT3", pt.size = 2) + NoLegend()
ggsave("DDIT3.tiff", DDIT3, width = 12, height = 8, device = "tiff")

HSPA5= FeaturePlot(merged_seurat.beta,features="HSPA5", pt.size = 2) + NoLegend()
ggsave("HSPA5.tiff", HSPA5, width = 12, height = 8, device = "tiff")

SDF2L1= FeaturePlot(merged_seurat.betatest,features="SDF2L1", pt.size = 2) + NoLegend()
ggsave("SDF2L1.tiff", SDF2L1, width = 12, height = 8, device = "tiff")

MANF= FeaturePlot(merged_seurat.betatest,features="MANF", pt.size = 2) + NoLegend()
ggsave("MANF.tiff", MANF, width = 12, height = 8, device = "tiff")

SEL1L= FeaturePlot(merged_seurat.betatest,features="SEL1L",pt.size = 2) + NoLegend()
ggsave("SEL1L.tiff", SEL1L, width = 12, height = 8, device = "tiff")

PPP1R15A= FeaturePlot(merged_seurat.beta,features="PPP1R15A", pt.size = 2)
ggsave("PPP1R15A.tiff", PPP1R15A, width = 12, height = 8, device = "tiff")


Q= FeaturePlot(merged_seurat.beta,features="DLK1",pt.size = 2) + ggtitle("DLK1 expression in beta cells")
ggsave("Q- DLK1 expression in beta cells.tiff", Q, width = 12, height = 8, device = "tiff")
ggsave("Q- DLK1 expression in beta cells.pdf", Q, width = 12, height = 8, device = "pdf")

save.image(file="/vol/data02/groy_data/Single_cell_RNA_seq/scRNAmerge_vs_Marselli/BMI/DEGASBMI_AMawla_data.RData")


