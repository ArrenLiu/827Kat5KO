# Preparing for analysis
library(Seurat)
library(dplyr)
genome = 'hg19'

setwd('/Users/arrenliu/Library/Mobile Documents/com~apple~CloudDocs/Important/Research/Plaisier/Single cell data/')
load('rf1.RData')

###########
### 827 ###
###########
# Setup Seurat object
exp_mat = Read10X(data.dir='827_ctrl data')
ctrl827 = CreateSeuratObject(counts = exp_mat, project = "ctrl827")

####################
### 827 KAT5 K/O ###
####################
# Setup Seurat object
exp_mat = Read10X(data.dir='827_kat5ko')
KAT5ko827 = CreateSeuratObject(counts = exp_mat, project = "KAT5ko827")

###########
### 131 ###
###########
exp_mat = Read10X(data.dir='131')
ctrl131 = CreateSeuratObject(counts = exp_mat, project = "131")

##################
### 131_kat5ko ###
##################
# Setup Seurat object
exp_mat = Read10X(data.dir='131_kat5ko')
ctrl131_kat5ko = CreateSeuratObject(counts = exp_mat, project = "131_kat5ko")

#################
### U5_kat5ko ###
#################
# Setup Seurat object
exp_mat = Read10X(data.dir='U5_kat5ko')
ctrlU5_kat5ko = CreateSeuratObject(counts = exp_mat, project = "ctrlU5_kat5ko")

#############
### U5 WT ###
#############
# Setup Seurat object
exp_mat = Read10X(data.dir='U5_WT')
WT = CreateSeuratObject(counts = exp_mat, project = "WT")


################
### Combining ###
################
single_cell.combined = merge(x = ctrl827, y = KAT5ko827, 
                                            add.cell.ids = c("ctrl827", "KAT5ko827"), project = "combined")

single_cell.combined = merge(x = single_cell.combined, y = ctrl131, 
                             add.cell.ids = c("single_cell.combined", "ctrl131"), project = "combined")

single_cell.combined = merge(x = single_cell.combined, y = ctrl131_kat5ko, 
                             add.cell.ids = c("single_cell.combined", "ctrl131_kat5ko"), project = "combined")

single_cell.combined = merge(x = single_cell.combined, y = ctrlU5_kat5ko, 
                             add.cell.ids = c("single_cell.combined", "ctrlU5_kat5ko"), project = "combined")

single_cell.combined = merge(x = single_cell.combined, y = WT, 
                             add.cell.ids = c("single_cell.combined", "WT"), project = "combined")

#################
### Filtering ###
#################
single_cell.combined@meta.data$pert = "single_cell.combined"
mito.genes = grep(pattern = "^MT-", x = rownames(x = single_cell.combined@meta.data), value = TRUE)
# mito.genes = grep(pattern = "^MT-", x = rownames(single_cell.combined), value = TRUE)
single_cell.combined[["percent.mt"]] = PercentageFeatureSet(object = single_cell.combined, pattern = "^MT-")
# percent.mito = colSums(single_cell.combined@meta.data[mito.genes, ]) / colSums(single_cell.combined@meta.data)
VlnPlot(object = single_cell.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# single_cell.combined = AddMetaData(object = single_cell.combined, metadata = percent.mito, col.name = "percent.mito")

plot1 <- FeatureScatter(object = single_cell.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = single_cell.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# single_cell.combined <- subset(x = single_cell.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
single_cell.combined <- subset(x = single_cell.combined, subset = percent.mt < 5)

single_cell.combined <- NormalizeData(object = single_cell.combined, normalization.method = "LogNormalize", scale.factor = 10000)
# single_cell.combined = NormalizeData(object = single_cell.combined)

# single_cell.combined = FindVariableGenes(object = single_cell.combined)
single_cell.combined <- FindVariableFeatures(object = single_cell.combined, selection.method = "mean.var.plot", nfeatures = 5000)
top10 <- head(x = VariableFeatures(object = single_cell.combined), 10)

plot1 <- VariableFeaturePlot(object = single_cell.combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
# CombinePlots(plots = list(plot1, plot2))
plot1
plot2

# single_cell.combined = ScaleData(object = single_cell.combined, vars.to.regress = c('percent.mito','nUMI'), genes.use = single_cell.combined@var.genes, model.use = 'negbinom')
all.genes <- rownames(x = single_cell.combined)
single_cell.combined <- ScaleData(object = single_cell.combined, features = all.genes)
single_cell.combined <- ScaleData(object = single_cell.combined, vars.to.regress = "percent.mt")

# single_cell.combined = RunPCA(object = single_cell.combined, features = VariableFeatures(object = single_cell.combined), verbose = FALSE)
single_cell.combined <- RunPCA(object = single_cell.combined, features = VariableFeatures(object = single_cell.combined))
print(x = single_cell.combined[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(object = single_cell.combined, dims = 1:2, reduction = "pca")
DimPlot(object = single_cell.combined, reduction = "pca")

# single_cell.combined = JackStraw(object = single_cell.combined, num.replicate = 100)
# ScoreJackStraw(object = single_cell.combined)
# JackStrawPlot(object = single_cell.combined)
# PCElbowPlot(object = single_cell.combined)

pdf("single_cell.combined 1000 Headmap.pdf")
DimHeatmap(object = single_cell.combined, dims = 1:6, cells = 1000, balanced = TRUE)
DimHeatmap(object = single_cell.combined, dims = 7:12, cells = 1000, balanced = TRUE)
DimHeatmap(object = single_cell.combined, dims = 13:15, cells = 1000, balanced = TRUE)
dev.off()

# single_cell.combined <- JackStraw(object = single_cell.combined, num.replicate = 100)
# single_cell.combined <- ScoreJackStraw(object = single_cell.combined, dims = 1:20)
# JackStrawPlot(object = single_cell.combined, dims = 1:15)

# ElbowPlot(object = single_cell.combined)

# Clustering
# single_cell.combined = FindClusters(single_cell.combined, reduction.type = "pca", resolution = 0.6, dims.use = 1:20)
single_cell.combined <- FindNeighbors(object = single_cell.combined, dims = 1:20)
single_cell.combined <- FindClusters(object = single_cell.combined, resolution = 0.6)

# single_cell.combined = SetAllIdent(object = single_cell.combined, id='res.0.6')
object = single_cell.combined
# id='res.1'
id='res.0.6'
if (id %in% colnames(object@meta.data)) {
  cells.use <- rownames(object@meta.data)
  ident.use <- object@meta.data[, id]
  object <- SetIdent(object,cells.use, ident.use)
}
single_cell.combined = object

single_cell.combined = RunTSNE(object = single_cell.combined, dims.use = 1:20, do.fast = TRUE)
pdf("single_cell.combined_tSNE_4_29_19.pdf")
TSNEPlot(object = single_cell.combined,do.label=T)
dev.off()

single_cell.combined.markers = FindAllMarkers(object = single_cell.combined, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
single_cell.combined.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
top_markers = single_cell.combined.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
top_markers

# rf2 = BuildRFClassifier(single_cell.combined, training.genes=single_cell.combined@var.genes, training.classes=single_cell.combined@ident)
# object = single_cell.combined
# training.genes = single_cell.combined@assays$RNA@var.features
# training.classes = single_cell.combined@active.ident
# training.classes <- as.vector(training.classes)
# # training.genes <- SetIfNull(training.genes, default = rownames(object@assays$RNA@data))
# # training.data <- as.data.frame(as.matrix(t(object@assays$RNA@data[training.genes, ])))
# training.data <- object@assays$RNA@data[training.genes, ]
# training.data <- as.data.frame(t(as.matrix(object@assays$RNA@data[training.genes, ])))
# training.data$class <- factor(training.classes)
# classifier <- ranger::ranger(data = training.data, dependent.variable.name = "class", classification = TRUE, write.forest = TRUE)
# rf2 = classifier

# Classifying single_cell.combined
# clusts.WT.single_cell.combined = ClassifyCells(single_cell.combined,classifier=rf1,new.data=single_cell.combined@data)
object = single_cell.combined
classifier = rf1
new.data = single_cell.combined@assays$RNA@data
features <- classifier$forest$independent.variable.names
genes.to.add <- setdiff(features, rownames(new.data))
data.to.add <- matrix(data = 0, nrow = length(genes.to.add), ncol = ncol(new.data))
rownames(data.to.add) <- genes.to.add
new.data <- rbind(new.data, data.to.add)
new.data <- new.data[features, ]
new.data <- t(as.matrix((new.data)))
prediction <- predict(classifier, new.data)
new.classes <- prediction$predictions
clusts.WT.single_cell.combined = new.classes

names(clusts.WT.single_cell.combined) = colnames(single_cell.combined@assays$RNA@data)
levels(clusts.WT.single_cell.combined)[2] = 'Neural/G0'
single_cell.combined = AddMetaData(object = single_cell.combined, metadata = clusts.WT.single_cell.combined, col.name = "clusts_WT")
table(single_cell.combined@meta.data$clusts_WT)/ncol(single_cell.combined@assays$RNA@data)

# single_cell.combined = SetAllIdent(object = single_cell.combined, id='clusts_WT')
object = single_cell.combined
id='clusts_WT'
if (id %in% colnames(object@meta.data)) {
  cells.use <- rownames(object@meta.data)
  ident.use <- object@meta.data[, id]
  object <- SetIdent(object,cells.use, ident.use)
}
single_cell.combined = object

# table(single_cell.combined@meta.data$clusts_WT,single_cell.combined@meta.data$RNA_snn_res.1)
table(single_cell.combined@meta.data$clusts_WT,single_cell.combined@meta.data$RNA_snn_res.0.6)

# Clustering
single_cell.combined = FindClusters(single_cell.combined, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:20)
# single_cell.combined = FindClusters(single_cell.combined, reduction.type = "cca.aligned", resolution = 1.2, dims.use = 1:20)

# single_cell.combined = SetAllIdent(object = single_cell.combined, id='res.0.6')
object = single_cell.combined
# id='res.1'
id='res.0.6'
if (id %in% colnames(object@meta.data)) {
  cells.use <- rownames(object@meta.data)
  ident.use <- object@meta.data[, id]
  object <- SetIdent(object,cells.use, ident.use)
}
single_cell.combined = object

current.cluster.ids = c(0,1,2,3,4,5,6,7,8,9)
# new.cluster.ids = c('G1','G1/Differentiation','G2/M','S/G2','M/Early G1','S','Late G1','G1/other', '8', '9', '10', '11')
new.cluster.ids = c('G1','G1/Differentiation','G2/M','S/G2','M/Early G1','S','Late G1','G1/other', '8', '9')
# current.cluster.ids = c(0,1,2,3,4)
# new.cluster.ids = c('G1','G1/Differentiation','G2/M','S/G2','M/Early G1')
tmp = plyr::mapvalues(single_cell.combined@active.ident, from=current.cluster.ids, to=new.cluster.ids)
single_cell.combined = AddMetaData(object = single_cell.combined, metadata = tmp, col.name = "clusts_named")
# single_cell.combined = SetAllIdent(object = single_cell.combined, id='clusts_named')
object = single_cell.combined
id='clusts_named'
if (id %in% colnames(object@meta.data)) {
  cells.use <- rownames(object@meta.data)
  ident.use <- object@meta.data[, id]
  object <- SetIdent(object,cells.use, ident.use)
}
single_cell.combined = object
pdf("single_cell.combined_tSNE_Cycle.pdf")
TSNEPlot(single_cell.combined)
dev.off()

names(x = new.cluster.ids) <- levels(x = single_cell.combined)
single_cell.combined <- RenameIdents(object = single_cell.combined, new.cluster.ids)

single_cell.combined <- RunUMAP(object = single_cell.combined, dims = 1:20)
pdf("single_cell.combined_uMAP_Cycle.pdf")
FeaturePlot(object = single_cell.combined, features = c("IGFBP5", "TAGLN", "SPP1", "MYL9"))
FeaturePlot(object = single_cell.combined, features = c("FABP7", "GPM6B", "SOX4", "PRCP"))
FeaturePlot(object = single_cell.combined, features = c("S100A6", "SERPINE2", "FTH1", "ATP5F1E"))
FeaturePlot(object = single_cell.combined, features = c("CDK1", "KIAA0101", "RRM2", "HIST1H4C"))
FeaturePlot(object = single_cell.combined, features = c("CCNB1", "CENPF", "PTTG1", "HMGB2"))
FeaturePlot(object = single_cell.combined, features = c("ATP5MG", "ATP5F1E", "ATP5MF", "FLG"))
FeaturePlot(object = single_cell.combined, features = c("BST2", "EIF5A", "MGST1", "THY1"))
FeaturePlot(object = single_cell.combined, features = c("RP11-745C15.2", "TPM2", "RPS29", "SEC61G"))
FeaturePlot(object = single_cell.combined, features = c("AC090498.1", "NEFL", "RPS17", "UQCR11"))
FeaturePlot(object = single_cell.combined, features = c("CD24", "NEFL", "C1QL1", "RPS17"))
DimPlot(object = single_cell.combined, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

pdf("single_cell.combined_barplot.pdf")
VlnPlot(single_cell.combined, c('TOP2A','S100B','UBE2S','HIST1H4C','UBE2C'))
VlnPlot(single_cell.combined, 'S100B')
# boxplot(single_cell.combined@assays$RNA@scale.data['S100B',]~ctrl827@meta.data$clusts_WT,las=2)
hist(single_cell.combined@meta.data$nCount_RNA,breaks=50,col=rgb(1,0,0,0.5),freq=T,ylim=c(0,125))
dev.off()
