# Preparing for analysis
library(Seurat)
library(dplyr)
genome = 'hg19'

setwd('/Users/arrenliu/Library/Mobile Documents/com~apple~CloudDocs/Important/Research/Plaisier/Single cell data/')
load('rf1.RData')

####################
### 827 KAT5 K/O ###
####################
# Setup Seurat object
# setwd('C:/Users/plais/Dropbox (ASU)/Single cell data')
exp_mat = Read10X(data.dir='827_kat5ko')
Kat5KO827 = CreateSeuratObject(counts = exp_mat, project = "Kat5KO827")
Kat5KO827@meta.data$pert = "Kat5KO827"
mito.genes = grep(pattern = "^MT-", x = rownames(x = Kat5KO827@meta.data), value = TRUE)
Kat5KO827[["percent.mt"]] = PercentageFeatureSet(object = Kat5KO827, pattern = "^MT-")

pdf("Kat5KO827 Vlnplot.pdf")
VlnPlot(object = Kat5KO827, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# Kat5KO827 = AddMetaData(object = Kat5KO827, metadata = percent.mito, col.name = "percent.mito")
dev.off()

plot1 <- FeatureScatter(object = Kat5KO827, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = Kat5KO827, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("Kat5KO827 FeatureScatter.pdf")
CombinePlots(plots = list(plot1, plot2))
dev.off()

# Kat5KO827 <- subset(x = Kat5KO827, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

Kat5KO827 <- NormalizeData(object = Kat5KO827, normalization.method = "LogNormalize", scale.factor = 10000)
Kat5KO827 <- FindVariableFeatures(object = Kat5KO827, selection.method = "mean.var.plot", nfeatures = 5000)
top10 <- head(x = VariableFeatures(object = Kat5KO827), 10)

plot1 <- VariableFeaturePlot(object = Kat5KO827)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
# CombinePlots(plots = list(plot1, plot2))
pdf("Kat5KO827 VariableFeaturePlot.pdf")
plot1
plot2
dev.off()

# Kat5KO827 = ScaleData(object = Kat5KO827, vars.to.regress = c('percent.mito','nUMI'), genes.use = Kat5KO827@var.genes, model.use = 'negbinom')
all.genes <- rownames(x = Kat5KO827)
Kat5KO827 <- ScaleData(object = Kat5KO827, features = all.genes)
Kat5KO827 <- ScaleData(object = Kat5KO827, vars.to.regress = "percent.mt")

# Kat5KO827 = RunPCA(object = Kat5KO827, features = VariableFeatures(object = Kat5KO827), verbose = FALSE)
Kat5KO827 <- RunPCA(object = Kat5KO827, features = VariableFeatures(object = Kat5KO827))
print(x = Kat5KO827[["pca"]], dims = 1:5, nfeatures = 5)

pdf("Kat5KO827DimPlot.pdf")
VizDimLoadings(object = Kat5KO827, dims = 1:5, reduction = "pca")
DimPlot(object = Kat5KO827, reduction = "pca")
dev.off()

pdf("Kat5KO827 1000 Headmap.pdf")
DimHeatmap(object = Kat5KO827, dims = 1:6, cells = 1000, balanced = TRUE)
DimHeatmap(object = Kat5KO827, dims = 7:12, cells = 1000, balanced = TRUE)
DimHeatmap(object = Kat5KO827, dims = 13:15, cells = 1000, balanced = TRUE)
dev.off()

Kat5KO827 <- JackStraw(object = Kat5KO827, num.replicate = 100)
Kat5KO827 <- ScoreJackStraw(object = Kat5KO827, dims = 1:20)

pdf("Kat5KO827 JackElbowPlot.pdf")
JackStrawPlot(object = Kat5KO827, dims = 1:15)
ElbowPlot(object = Kat5KO827)
dev.off()

# Clustering
# Kat5KO827 = FindClusters(Kat5KO827, reduction.type = "pca", resolution = 0.6, dims.use = 1:20)
# Kat5KO827 = SetAllIdent(object = Kat5KO827, id='res.0.6')
# Kat5KO827 = RunTSNE(object = Kat5KO827, dims.use = 1:10, do.fast = TRUE)
# TSNEPlot(object = Kat5KO827,do.label=T)

Kat5KO827 = RunTSNE(object = Kat5KO827, dims.use = 1:4, do.fast = TRUE, perplexity = 100)
TSNEPlot(Kat5KO827)

# Find clusters
Kat5KO827 = FindClusters(object = Kat5KO827, reduction.type = "pca", dims.use = 1:10, resolution = 0.3, print.output = 0, save.SNN = TRUE)
TSNEPlot(Kat5KO827)

# Find clusters at varying resolutions
Kat5KO827 = FindClusters(object = Kat5KO827, reduction.type = "pca", resolution = c(0.02, 0.1, .5, 1, 3), print.output = 0, reuse.SNN  = TRUE)
TSNEPlot(Kat5KO827)

Kat5KO827 <- FindNeighbors(object = Kat5KO827, dims = 1:15)
Kat5KO827 <- FindClusters(object = Kat5KO827, resolution = 1)
# Look at cluster IDs of the first 5 cells
head(x = Idents(object = Kat5KO827), 5)

Kat5KO827 <- RunUMAP(object = Kat5KO827, dims = 1:15)
plot1 = DimPlot(object = Kat5KO827, reduction = "umap")
plot2 = TSNEPlot(Kat5KO827)
pdf("Kat5KO827 DimTSNE.pdf")
# CombinePlots(plots = list(plot1, plot2))
plot1
plot2
dev.off()

# find all markers of cluster 0
cluster0.markers <- FindMarkers(object = Kat5KO827, ident.1 = 0, min.pct = 0.25)
head(x = cluster0.markers, n = 5)

# find all markers of cluster 1
cluster1.markers <- FindMarkers(object = Kat5KO827, ident.1 = 1, min.pct = 0.25)
head(x = cluster1.markers, n = 5)

# find all markers of cluster 2
cluster2.markers <- FindMarkers(object = Kat5KO827, ident.1 = 2, min.pct = 0.25)
head(x = cluster2.markers, n = 5)

# find all markers of cluster 3
cluster3.markers <- FindMarkers(object = Kat5KO827, ident.1 = 3, min.pct = 0.25)
head(x = cluster3.markers, n = 5)

# # find all markers of cluster 4
# cluster4.markers <- FindMarkers(object = Kat5KO827, ident.1 = 4, min.pct = 0.25)
# head(x = cluster4.markers, n = 5)

# find all markers distinguishing cluster 0 from clusters 1, 2, 3, and 4
cluster0.markers <- FindMarkers(object = Kat5KO827, ident.1 = 0, ident.2 = c(1,2,3), min.pct = 0.25)
head(x = cluster0.markers, n = 5)

# find all markers distinguishing cluster 1 from clusters 0, 2, 3, and 4
cluster1.markers <- FindMarkers(object = Kat5KO827, ident.1 = 1, ident.2 = c(0,2,3), min.pct = 0.25)
head(x = cluster1.markers, n = 5)

# find all markers distinguishing cluster 2 from clusters 0, 1, 3, and 4
cluster2.markers <- FindMarkers(object = Kat5KO827, ident.1 = 2, ident.2 = c(0,1,3), min.pct = 0.25)
head(x = cluster2.markers, n = 5)

# find all markers distinguishing cluster 3 from clusters 0, 1, 2, and 4
cluster3.markers <- FindMarkers(object = Kat5KO827, ident.1 = 3, ident.2 = c(0,1,2), min.pct = 0.25)
head(x = cluster3.markers, n = 5)

# find all markers distinguishing cluster 4 from clusters 0, 1, 2, and 3
cluster4.markers <- FindMarkers(object = Kat5KO827, ident.1 = 4, ident.2 = c(0,1,2,3), min.pct = 0.25)
head(x = cluster4.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
Kat5KO827.markers <- FindAllMarkers(object = Kat5KO827, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Kat5KO827.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

cluster1.markers <- FindMarkers(object = Kat5KO827, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", 
                                only.pos = TRUE)

# pdf("Kat5KO827 VlnFeature.pdf")
# VlnPlot(object = Kat5KO827, features = c("TOP2A", "CKS2"))
# VlnPlot(object = Kat5KO827, features = c("HIST1H4C", "HIST1H1D"))
# VlnPlot(object = Kat5KO827, features = c("MALAT1", "SNHG10"))
# VlnPlot(object = Kat5KO827, features = c("BIRC5", "HMGB2"))
# VlnPlot(object = Kat5KO827, features = c("HSP90AA1", "ACTB"))
# FeaturePlot(object = Kat5KO827, features = c("TOP2A", "CKS2"))
# FeaturePlot(object = Kat5KO827, features = c("HIST1H4C", "HIST1H1D"))
# FeaturePlot(object = Kat5KO827, features = c("MALAT1", "SNHG10"))
# FeaturePlot(object = Kat5KO827, features = c("BIRC5", "HMGB2"))
# FeaturePlot(object = Kat5KO827, features = c("HSP90AA1", "ACTB"))
# dev.off()

pdf("2Kat5KO827 VlnFeature.pdf")
TSNEPlot(Kat5KO827)
VlnPlot(object = Kat5KO827, features = c("NDFIP1", "SPARC"))
VlnPlot(object = Kat5KO827, features = c("ACTB", "PKM"))
VlnPlot(object = Kat5KO827, features = c("CPA4", "S100B"))
VlnPlot(object = Kat5KO827, features = c("PBK", "TOP2A"))
VlnPlot(object = Kat5KO827, features = c("TUBA1A", "HSP90AA1"))
FeaturePlot(object = Kat5KO827, features = c("NDFIP1", "SPARC"))
FeaturePlot(object = Kat5KO827, features = c("ACTB", "PKM"))
FeaturePlot(object = Kat5KO827, features = c("CPA4", "S100B"))
FeaturePlot(object = Kat5KO827, features = c("PBK", "TOP2A"))
FeaturePlot(object = Kat5KO827, features = c("TUBA1A", "HSP90AA1"))
dev.off()

# Kat5KO827.markers = FindAllMarkers(object = Kat5KO827, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
# Kat5KO827.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)

top10 <- Kat5KO827.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(object = Kat5KO827, features = top10$gene) + NoLegend()