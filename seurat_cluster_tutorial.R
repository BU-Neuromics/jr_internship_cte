library(dplyr)
library(Seurat)
library(patchwork)

#loading datasets
counts_data <- read.csv("/restricted/projectnb/cteseq/projects/ALS_rnaseq/data/raw/exp_counts.csv", row.names = 1)
meta_data <- read.csv("/restricted/projectnb/cteseq/projects/ALS_rnaseq/data/raw/exp_meta.csv")
#filter metadata for ALS only under Primary_Dx
meta_data_filtered <- filter(meta_data, Primary_Dx == "ALS")
meta_data_filtered <- as.data.frame(meta_data_filtered)
rownames(meta_data_filtered) <- meta_data_filtered$Core_ID

#Gets each unique value from the Core_ID column in metadata and assigns it to "sample"
sample <- unique(meta_data_filtered$Core_ID)
#subsets the values from "sample" from the counts data so only ALS data remains
counts_data_filtered <- counts_data[,sample]
#passing counts and meta data to a seurat objext
obj <- CreateSeuratObject(counts = counts_data_filtered, meta.data = meta_data_filtered)
obj

#QC section of tutorial
#example plots
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
ftplot <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ftplot
#Filtering step. This deviates from the tutorial in that mitochondrial genome is omitted
#and there is no upper limit on features to allow PCA to work properly
obj <- subset(obj, subset = nFeature_RNA > 200)

#Normalizing data
obj <- NormalizeData(obj)
#Identifying most high variable features
obj <- FindVariableFeatures(obj, selection.method = "vst")
#Plot variable features, displaying without top 10 and with. Produces viewport error
top10 <- head(VariableFeatures(obj), 10)
fsplot1 <- VariableFeaturePlot(obj)
fsplot2 <- LabelPoints(plot = fsplot1, points = top10, repel = TRUE)
fsplot1 + fsplot2

#Data scaling, produces viewport error
all.genes <- rownames(obj)
obj <- ScaleData(obj, features = all.genes)

#Linear dimensional reduction
#runs pca
obj <- RunPCA(obj, features = VariableFeatures(object = obj))
print(obj[["pca"]], dims = 1.5, nFeatures = 5)
#plots
VizDimLoadings(obj, dims = 1:2, reduction = "pca")
DimPlot(obj, reduction = "pca")
#heatmaps
DimHeatmap(obj, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(obj, dims = 1:15, cells = 500, balanced = TRUE)
#dataset dimensionality
obj <- JackStraw(obj, num.replicate = 100)
obj <- ScoreJackStraw(obj, dims = 1:20)
JackStrawPlot(obj, dims = 1:15)
ElbowPlot(obj)
#Clustering cells
obj <- FindNeighbors(obj, dims = 1:10)
obj <- FindClusters(obj, resolution = 0.5)
head(Idents(obj),5)
#UMAP/tSNE
obj <- RunUMAP(obj, dims = 1:10)
DimPlot(obj, reduction = "umap")
