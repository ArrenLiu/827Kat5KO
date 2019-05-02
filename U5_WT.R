# Preparing for analysis
library(Seurat)
library(dplyr)
genome = 'hg19'

setwd('/Users/arrenliu/Library/Mobile Documents/com~apple~CloudDocs/Important/Research/Plaisier/Single cell data/')
load('rf1.RData')

#############
### U5 WT ###
#############
# Setup Seurat object
exp_mat = Read10X(data.dir='U5_WT')
U5_WT = CreateSeuratObject(counts = exp_mat, project = "WT")

#################
### Filtering ###
#################
U5_WT@meta.data$pert = "U5_WT"
mito.genes = grep(pattern = "^MT-", x = rownames(x = U5_WT@meta.data), value = TRUE)
# mito.genes = grep(pattern = "^MT-", x = rownames(U5_WT), value = TRUE)
U5_WT[["percent.mt"]] = PercentageFeatureSet(object = U5_WT, pattern = "^MT-")
# percent.mito = colSums(U5_WT@meta.data[mito.genes, ]) / colSums(U5_WT@meta.data)
VlnPlot(object = U5_WT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# U5_WT = AddMetaData(object = U5_WT, metadata = percent.mito, col.name = "percent.mito")

plot1 <- FeatureScatter(object = U5_WT, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = U5_WT, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# U5_WT <- subset(x = U5_WT, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
U5_WT <- subset(x = U5_WT, subset = percent.mt < 25)

U5_WT <- NormalizeData(object = U5_WT, normalization.method = "LogNormalize", scale.factor = 10000)
# U5_WT = NormalizeData(object = U5_WT)

# U5_WT = FindVariableGenes(object = U5_WT)
U5_WT <- FindVariableFeatures(object = U5_WT, selection.method = "mean.var.plot", nfeatures = 5000)
top10 <- head(x = VariableFeatures(object = U5_WT), 10)

plot1 <- VariableFeaturePlot(object = U5_WT)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
# CombinePlots(plots = list(plot1, plot2))
plot1
plot2

# U5_WT = ScaleData(object = U5_WT, vars.to.regress = c('percent.mito','nUMI'), genes.use = U5_WT@var.genes, model.use = 'negbinom')
all.genes <- rownames(x = U5_WT)
U5_WT <- ScaleData(object = U5_WT, features = all.genes)
U5_WT <- ScaleData(object = U5_WT, vars.to.regress = "percent.mt")

# U5_WT = RunPCA(object = U5_WT, features = VariableFeatures(object = U5_WT), verbose = FALSE)
U5_WT <- RunPCA(object = U5_WT, features = VariableFeatures(object = U5_WT))
print(x = U5_WT[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(object = U5_WT, dims = 1:2, reduction = "pca")
DimPlot(object = U5_WT, reduction = "pca")

# U5_WT = JackStraw(object = U5_WT, num.replicate = 100)
# ScoreJackStraw(object = U5_WT)
# JackStrawPlot(object = U5_WT)
# PCElbowPlot(object = U5_WT)

pdf("U5_WT 1000 Headmap.pdf")
DimHeatmap(object = U5_WT, dims = 1:6, cells = 1000, balanced = TRUE)
DimHeatmap(object = U5_WT, dims = 7:12, cells = 1000, balanced = TRUE)
DimHeatmap(object = U5_WT, dims = 13:15, cells = 1000, balanced = TRUE)
dev.off()

# U5_WT <- JackStraw(object = U5_WT, num.replicate = 100)
# U5_WT <- ScoreJackStraw(object = U5_WT, dims = 1:20)
# JackStrawPlot(object = U5_WT, dims = 1:15)

# ElbowPlot(object = U5_WT)

# Clustering
# U5_WT = FindClusters(U5_WT, reduction.type = "pca", resolution = 0.6, dims.use = 1:20)
U5_WT <- FindNeighbors(object = U5_WT, dims = 1:20)
U5_WT <- FindClusters(object = U5_WT, resolution = 0.6)

# U5_WT = SetAllIdent(object = U5_WT, id='res.0.6')
object = U5_WT
# id='res.1'
id='res.0.6'
if (id %in% colnames(object@meta.data)) {
  cells.use <- rownames(object@meta.data)
  ident.use <- object@meta.data[, id]
  object <- SetIdent(object,cells.use, ident.use)
}
U5_WT = object

U5_WT = RunTSNE(object = U5_WT, dims.use = 1:20, do.fast = TRUE)
pdf("U5_WT_tSNE_4_29_19.pdf")
TSNEPlot(object = U5_WT,do.label=T)
dev.off()

U5_WT.markers = FindAllMarkers(object = U5_WT, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
U5_WT.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
top_markers = U5_WT.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
top_markers



################################################################################################
################################################################################################
################################ Create random forest classifier################################
################################################################################################
################################################################################################

##rf2 = BuildRFClassifier(U5_WT, training.genes=U5_WT@var.genes, training.classes=U5_WT@ident)
# Make data object for classification
object = U5_WT
training.genes = U5_WT@assays$RNA@var.features
training.classes = U5_WT@active.ident
training.classes <- as.vector(training.classes)
training.data <- object@assays$RNA@data[training.genes, ]
training.data <- as.data.frame(t(as.matrix(object@assays$RNA@data[training.genes, ])))
# training.data2 <- as.data.frame(t(as.matrix(object@assays$RNA@data[training.genes, ])))
training.data$class <- factor(training.classes)

U5_WT@assays$RNA@scale.data[1:5,1:5]
dim(U5_WT@assays$RNA@scale.data)
training.data$class[1:5]
length(training.data$class)
# tmp = cbind(t(U5_WT@assays$RNA@scale.data),status=as.factor(as.numeric(training.data$class)))

# mads = apply(U5_WT@assays$RNA@scale.data,1,mad)
# tmp = U5_WT@assays$RNA@scale.data[order(mads,decreasing=T),]
# tmp = cbind(t(tmp),status=as.numeric(as.factor(training.data$class)))
tmp2 = U5_WT@assays$RNA@scale.data
rownames(tmp2) <- make.names(rownames(tmp2))
tmp = cbind(t(tmp2),status = as.numeric(as.factor(training.data$class)))

library(randomForest)
U5.rf = randomForest(as.factor(status) ~ .,data=tmp,importance=T,proximity=T)
U5.rf$importance[order(U5.rf$importance[,4],decreasing=T),]
votes = U5.rf$votes
# Get gene names for most important features
# featureData(gset)@data[rownames(U5.rf$importance[order(U5.rf$importance[,4],decreasing=T)[1:20],]),'Symbol']
# U5_WT@meta.data[rownames(U5.rf$importance[order(U5.rf$importance[,4],decreasing=T)[1:20],]),'Symbol']
U5.rf

# # Calculate ROC AUC
# library(ROCR)
# predictions = as.vector(U5.rf$votes[,2])
# pred1 = prediction(predictions, tmp[,'status'])
# auc1 = performance(pred1, 'auc')
# auc1@y.values[[1]]

# Calculate predictive value of classifier
confusion1 = U5.rf$confusion[c(2,1),c(2,1)]
ppv1 = ((confusion1[1,1])/(confusion1[1,1]+confusion1[2,1]))
ppv1
npv1 = ((confusion1[2,2])/(confusion1[1,2]+confusion1[2,2]))
npv1

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################

# # Classifying U5_WT
# # clusts.WT.U5_WT = ClassifyCells(U5_WT,classifier=rf1,new.data=U5_WT@data)
# object = U5_WT
# classifier = rf1
# new.data = U5_WT@assays$RNA@data
# features <- classifier$forest$independent.variable.names
# genes.to.add <- setdiff(features, rownames(new.data))
# data.to.add <- matrix(data = 0, nrow = length(genes.to.add), ncol = ncol(new.data))
# rownames(data.to.add) <- genes.to.add
# new.data <- rbind(new.data, data.to.add)
# new.data <- new.data[features, ]
# new.data <- t(as.matrix((new.data)))
# prediction <- predict(classifier, new.data)
# new.classes <- prediction$predictions
# clusts.WT.U5_WT = new.classes

names(clusts.WT.U5_WT) = colnames(U5_WT@assays$RNA@data)
levels(clusts.WT.U5_WT)[2] = 'Neural/G0'
U5_WT = AddMetaData(object = U5_WT, metadata = clusts.WT.U5_WT, col.name = "clusts_WT")
table(U5_WT@meta.data$clusts_WT)/ncol(U5_WT@assays$RNA@data)

# U5_WT = SetAllIdent(object = U5_WT, id='clusts_WT')
object = U5_WT
id='clusts_WT'
if (id %in% colnames(object@meta.data)) {
  cells.use <- rownames(object@meta.data)
  ident.use <- object@meta.data[, id]
  object <- SetIdent(object,cells.use, ident.use)
}
U5_WT = object

# table(U5_WT@meta.data$clusts_WT,U5_WT@meta.data$RNA_snn_res.1)
table(U5_WT@meta.data$clusts_WT,U5_WT@meta.data$RNA_snn_res.0.6)

# Clustering
U5_WT = FindClusters(U5_WT, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:20)
# U5_WT = FindClusters(U5_WT, reduction.type = "cca.aligned", resolution = 1.2, dims.use = 1:20)

# U5_WT = SetAllIdent(object = U5_WT, id='res.0.6')
object = U5_WT
# id='res.1'
id='res.0.6'
if (id %in% colnames(object@meta.data)) {
  cells.use <- rownames(object@meta.data)
  ident.use <- object@meta.data[, id]
  object <- SetIdent(object,cells.use, ident.use)
}
U5_WT = object

current.cluster.ids = c(0,1,2,3,4,5,6,7,8,9)
# new.cluster.ids = c('G1','G1/Differentiation','G2/M','S/G2','M/Early G1','S','Late G1','G1/other', '8', '9', '10', '11')
new.cluster.ids = c('G1','G1/Differentiation','G2/M','S/G2','M/Early G1','S','Late G1','G1/other', '8', '9')
# current.cluster.ids = c(0,1,2,3,4)
# new.cluster.ids = c('G1','G1/Differentiation','G2/M','S/G2','M/Early G1')
tmp = plyr::mapvalues(U5_WT@active.ident, from=current.cluster.ids, to=new.cluster.ids)
U5_WT = AddMetaData(object = U5_WT, metadata = tmp, col.name = "clusts_named")
# U5_WT = SetAllIdent(object = U5_WT, id='clusts_named')
object = U5_WT
id='clusts_named'
if (id %in% colnames(object@meta.data)) {
  cells.use <- rownames(object@meta.data)
  ident.use <- object@meta.data[, id]
  object <- SetIdent(object,cells.use, ident.use)
}
U5_WT = object
pdf("U5_WT_tSNE_Cycle.pdf")
TSNEPlot(U5_WT)
dev.off()

names(x = new.cluster.ids) <- levels(x = U5_WT)
U5_WT <- RenameIdents(object = U5_WT, new.cluster.ids)

U5_WT <- RunUMAP(object = U5_WT, dims = 1:20)
pdf("U5_WT_U5_WT_U5_kat5ko.pdf")
FeaturePlot(object = U5_WT, features = top_markers$gene[1:5])    #0
FeaturePlot(object = U5_WT, features = top_markers$gene[6:10])   #1
FeaturePlot(object = U5_WT, features = top_markers$gene[11:15])  #2
FeaturePlot(object = U5_WT, features = top_markers$gene[16:20])  #3
FeaturePlot(object = U5_WT, features = top_markers$gene[21:25])  #4
FeaturePlot(object = U5_WT, features = top_markers$gene[26:30])  #5
FeaturePlot(object = U5_WT, features = top_markers$gene[31:35])  #6
FeaturePlot(object = U5_WT, features = top_markers$gene[36:40])  #7
FeaturePlot(object = U5_WT, features = top_markers$gene[41:44])  #8
FeaturePlot(object = U5_WT, features = top_markers$gene[46:50])  #9
FeaturePlot(object = U5_WT, features = top_markers$gene[51:55])  #10
FeaturePlot(object = U5_WT, features = top_markers$gene[56:60])  #11
DimPlot(object = U5_WT, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

# pdf("U5_WT_U5_WT_U5_kat5ko.pdf")
VlnPlot(U5_WT_U5_WT_U5_kat5ko, 'S100B')
# boxplot(U5_WT@assays$RNA@scale.data['S100B',]~ctrl827@meta.data$clusts_WT,las=2)
# hist(U5_WT@meta.data$nCount_RNA,breaks=50,col=rgb(1,0,0,0.5),freq=T,ylim=c(0,125))
# dev.off()
