#### Add necessary tools to library ####
library(Seurat)
library(devtools)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
library(reticulate)
library(monocle)
library(RColorBrewer)
library(ggpubr)
library(GGally)
library(slingshot)
library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(scales)
library(viridis)
library(scater)
library(PseudotimeDE)
library(SingleCellExperiment)
library(tibble)
library(irlba)

####Loading data####
##HiMYC
Ctrlunfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Hi-myc-ARKO/2nd Run 6-4-2020/Sun Lab Analysis/Project_COHP_36595_1_XCTC1_count/outs/filtered_feature_bc_matrix")
Ctrlunfiltered <- CreateSeuratObject(counts = Ctrlunfiltered.data,  min.cells = 3, min.features = 500, project = "Ctrlunfiltered")
Ctrlunfiltered <- NormalizeData(Ctrlunfiltered)
##HiMYC-ARKO
MycARKOunfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Hi-myc-ARKO/190408_190425_combined/28345_count/outs/filtered_feature_bc_matrix")
MycARKOunfiltered <- CreateSeuratObject(counts = MycARKOunfiltered.data, project = "MycARKOSC", min.cells = 3, min.features = 200)
MycARKOunfiltered <- NormalizeData(MycARKOunfiltered)

####Initial processing, Filtering and Clustering####
##HiMYC
Ctrlunfiltered[["percent.mt"]] <- PercentageFeatureSet(Ctrlunfiltered, pattern = "^mt-")
VlnPlot(Ctrlunfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(Ctrlunfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "Ctrl Pre-filteration")
hist(Ctrlunfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "Ctrl Pre-filteration")
Ctrl <- subset(Ctrlunfiltered, subset = nFeature_RNA > 600 & nFeature_RNA < 6500 & percent.mt < 15)
VlnPlot(Ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(Ctrl@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "Ctrl Post-filteration")
hist(Ctrl@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "Ctrl Post-filteration")
Ctrl <- FindVariableFeatures(Ctrl, selection.method = "vst", nfeatures = 5000)
VariableFeaturePlot(Ctrl)
Ctrl <- ScaleData(Ctrl, verbose = FALSE)
Ctrl <- RunPCA(Ctrl, npcs = 30, verbose = FALSE)
ElbowPlot(Ctrl, ndims = 50)
Ctrl <- FindNeighbors(Ctrl, reduction = "pca", dims = 1:20)
Ctrl <- FindClusters(Ctrl, resolution = 0.5)
Ctrl <- RunTSNE(Ctrl, reduction = "pca", dims = 1:20)
Ctrl <- RunUMAP(Ctrl, reduction = "pca", dims = 1:20)
DimPlot(Ctrl, reduction = "umap", pt.size = 0.3)
##HiMYC-ARKO
MycARKOunfiltered[["percent.mt"]] <- PercentageFeatureSet(MycARKOunfiltered, pattern = "^mt-")
VlnPlot(MycARKOunfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(MycARKOunfiltered@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "MycARKO Pre-filteration")
hist(MycARKOunfiltered@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "MycARKO Pre-filteration")
MycARKO <- subset(MycARKOunfiltered, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
VlnPlot(MycARKO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(MycARKO@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "MycARKO Post-filteration")
hist(MycARKO@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "MycARKO Post-filteration")
MycARKO <- FindVariableFeatures(MycARKO, selection.method = "vst", nfeatures = 5000)
VariableFeaturePlot(MycARKO)
MycARKO <- ScaleData(MycARKO, verbose = FALSE)
MycARKO <- RunPCA(MycARKO, npcs = 30, verbose = FALSE)
ElbowPlot(MycARKO, ndims = 50)
MycARKO <- FindNeighbors(MycARKO, reduction = "pca", dims = 1:20)
MycARKO <- FindClusters(MycARKO, resolution = 0.5)
MycARKO <- RunTSNE(MycARKO, reduction = "pca", dims = 1:20)
MycARKO <- RunUMAP(MycARKO, reduction = "pca", dims = 1:20)
DimPlot(MycARKO, reduction = "umap", pt.size = 0.3)

####Merging Datasets####
Ctrl[["orig.clusters"]] <- Idents(object = Ctrl)
MycARKO[["orig.clusters"]] <- Idents(object = MycARKO)
Idents(object = Ctrl) <- "seurat_clusters"
Idents(object = MycARKO) <- "seurat_clusters"
Ctrl$stim <- "Ctrl"
MycARKO$stim <- "ARKO"
CtrlvARKO.anchors <- FindIntegrationAnchors(object.list = list(Ctrl, MycARKO), dims = 1:20)
CtrlvARKO.combined <- IntegrateData(anchorset = CtrlvARKO.anchors, dims = 1:20)
DefaultAssay(CtrlvARKO.combined) <- "integrated"
#Run the standard workflow for visualization and clustering
CtrlvARKO.combined <- ScaleData(CtrlvARKO.combined, verbose = FALSE)
CtrlvARKO.combined <- RunPCA(CtrlvARKO.combined, npcs = 30, verbose = FALSE)
CtrlvARKO.combined <- FindNeighbors(CtrlvARKO.combined, reduction = "pca", dims = 1:20)
CtrlvARKO.combined <- FindClusters(CtrlvARKO.combined, resolution = 0.5)
CtrlvARKO.combined <- RunUMAP(CtrlvARKO.combined, reduction = "pca", dims = 1:20)
CtrlvARKO.combined <- RunTSNE(CtrlvARKO.combined, reduction = "pca", dims = 1:20)
Idents(object = CtrlvARKO.combined) <- "stim"
DimPlot(CtrlvARKO.combined, reduction = "umap", pt.size = 0.3, cols = c("light grey", "darkblue")) 
#DEGS
DefaultAssay(CtrlvARKO.combined) <- "RNA"
CtrlvARKO.combined <- ScaleData(CtrlvARKO.combined, features = rownames(CtrlvARKO.combined))
CtrlvARKO.combined.all.Markers <- FindAllMarkers(CtrlvARKO.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CtrlvARKO.combined.all.Markers, "CtrlvARKO.combined.all.seurat.Markers.csv")
#New Labeling
CtrlvARKO.combined <- RenameIdents(object = CtrlvARKO.combined, '0' = "BE", '1' = "FB",
                                   '2' = "FB", '3' = "BE", '4' = "LE", '5' = "LE", 
                                   '6' = "FB", '7' = "Leu", '8' = "LE", '9' = "Endo",
                                   '10' = "Leu", '11' = "SM", '12' = "Leu", '13' = "Leu",
                                   '14' = "LE", '15' = "Peri", '16' = "FB", '17' = "Glia",
                                   '18' = "Leu")
CtrlvARKO.combined[["celltype"]] <- Idents(object = CtrlvARKO.combined)
#Dotplot
Idents(object = CtrlvARKO.combined) <- "celltype"
DotPlot(CtrlvARKO.combined, features = c("Ar", "Pbsn", "MYC-transgene", "EGFP", "Gli1", 
                                         "Krt5", "Krt14", "Krt15", "Krt17", "Col17a1", 
                                         "Tspan1", "Cldn3", "Agr2", "Bex1", "Tspan8",
                                         "Fbln1", "Apod", "Lum", "Igf1", "Col3a1",
                                         "Acta1", "Acta2", "Myh11", "Cnn1", "Actg2",
                                         "Rgs4", "Rgs5", "Gja4", "S1pr3", "Ndufa4l2",
                                         "Rgs1", "Ccl5", "Tyrobp", "Fcer1g", "Cd52",
                                         "Plvap", "Pecam1", "Aqp1", "Cd93", "Cdh5",
                                         "Kcna1", "Cd59a", "Mpz", "Plp1", "Fabp7"), cols = c("light grey", "red")) + RotatedAxis()
#Umap plots
Idents(object = CtrlvARKO.combined) <- "stim"
CtrlvARKO.combined <- RenameIdents(object = CtrlvARKO.combined, 'Ctrl' = "Ctrl", 'ARKO' = "ARKO")
CtrlvARKO.combined[["stim"]] <- Idents(object = CtrlvARKO.combined)
Idents(object = CtrlvARKO.combined) <- "celltype"
DimPlot(CtrlvARKO.combined, reduction = "umap", pt.size = 0.3, cols = c("Red", "Purple", "Blue", "Green", "Orange", "Brown", "Gold", "Black"))
DimPlot(CtrlvARKO.combined, reduction = "umap", split.by = "stim", pt.size = 0.3, cols = c("Red", "Purple", "Blue", "Green", "Orange", "Brown", "Gold", "Black"))
#Feature plots
DefaultAssay(CtrlvARKO.combined) <- "RNA"
FeaturePlot(CtrlvARKO.combined, reduction = "umap", split.by = "stim", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", split.by = "stim", features = c('MYC-transgene'), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = 2.5)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", split.by = "stim", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
#Numbers of cells
Idents(object = CtrlvARKO.combined) <- "celltype"
CtrlvARKO.combined$celltype.stim <- paste(Idents(CtrlvARKO.combined), CtrlvARKO.combined$stim, sep = "_")
Idents(object = CtrlvARKO.combined) <- "celltype.stim"
table(Idents(CtrlvARKO.combined))

#### Re-clustering Epithelial cells ####
Idents(object = CtrlvARKO.combined) <- "celltype"
CtrlvARKO.combined.Epi <- subset(CtrlvARKO.combined, idents = c("BE", "LE"))
DefaultAssay(CtrlvARKO.combined.Epi) <- "integrated"
#Run the standard workflow for visualization and clustering
CtrlvARKO.combined.Epi <- ScaleData(CtrlvARKO.combined.Epi, verbose = FALSE)
CtrlvARKO.combined.Epi <- RunPCA(CtrlvARKO.combined.Epi, npcs = 30, verbose = FALSE)
ElbowPlot(CtrlvARKO.combined.Epi, ndims = 50)
#Umap, Tsne and Clustering
CtrlvARKO.combined.Epi <- FindNeighbors(CtrlvARKO.combined.Epi, reduction = "pca", dims = 1:16)
CtrlvARKO.combined.Epi <- FindClusters(CtrlvARKO.combined.Epi, resolution = 0.5)
CtrlvARKO.combined.Epi <- RunTSNE(CtrlvARKO.combined.Epi, reduction = "pca", dims = 1:16)
CtrlvARKO.combined.Epi <- RunUMAP(CtrlvARKO.combined.Epi, reduction = "pca", dims = 1:16)
DimPlot(CtrlvARKO.combined.Epi, reduction = "umap", label = TRUE)
#Cell cycle regression
mouse_cell_cycle_genes <- readRDS(file = "//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARKOvCtrl_3timepoint/E18.5/mouse_cell_cycle_genes.rds")
s.genes <- mouse_cell_cycle_genes$s.genes
g2m.genes <- mouse_cell_cycle_genes$g2m.genes
DefaultAssay(CtrlvARKO.combined.Epi) <- "RNA"
all.genes <- rownames(CtrlvARKO.combined.Epi)
CtrlvARKO.combined.Epi <- ScaleData(CtrlvARKO.combined.Epi, features = all.genes)
CtrlvARKO.combined.Epi <- CellCycleScoring(CtrlvARKO.combined.Epi, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Idents(object = CtrlvARKO.combined.Epi) <- "Phase"
CtrlvARKO.combined.Epi <- RenameIdents(object = CtrlvARKO.combined.Epi, 'G1' = "G1", 'G2M' = "G2M", 'S' = "S")
CtrlvARKO.combined.Epi[["Phase"]] <- Idents(object = CtrlvARKO.combined.Epi)
DimPlot(CtrlvARKO.combined.Epi, reduction = "umap", pt.size = 0.3, cols = c("#FFD966", "purple", "#1D762E"))
#Take Cell cycle out 
CtrlvARKO.combined.Epi1 <- CtrlvARKO.combined.Epi
DefaultAssay(CtrlvARKO.combined.Epi1) <- "integrated"
CtrlvARKO.combined.Epi1 <- ScaleData(CtrlvARKO.combined.Epi1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(CtrlvARKO.combined.Epi1))
CtrlvARKO.combined.Epi1 <- RunPCA(CtrlvARKO.combined.Epi1, features = VariableFeatures(CtrlvARKO.combined.Epi1))
ElbowPlot(CtrlvARKO.combined.Epi1, ndims = 50)
CtrlvARKO.combined.Epi1 <- FindNeighbors(CtrlvARKO.combined.Epi1, reduction = "pca", dims = 1:15)
CtrlvARKO.combined.Epi1 <- FindClusters(CtrlvARKO.combined.Epi1, resolution = 0.4)
CtrlvARKO.combined.Epi1 <- RunUMAP(CtrlvARKO.combined.Epi1, reduction = "pca", dims = 1:15)
CtrlvARKO.combined.Epi1 <- RunTSNE(CtrlvARKO.combined.Epi1, reduction = "pca", dims = 1:15)
DimPlot(CtrlvARKO.combined.Epi1, reduction = "umap", pt.size = 0.3, label = TRUE)
Idents(object = CtrlvARKO.combined.Epi1) <- "Phase"
CtrlvARKO.combined.Epi1 <- RenameIdents(object = CtrlvARKO.combined.Epi1, 'G1' = "G1", 'G2M' = "G2M", 'S' = "S")
CtrlvARKO.combined.Epi1[["Phase"]] <- Idents(object = CtrlvARKO.combined.Epi1)
DimPlot(CtrlvARKO.combined.Epi1, reduction = "umap", pt.size = 0.3, cols = c("#FFD966", "purple", "#1D762E"))
#Rename
Idents(object = CtrlvARKO.combined.Epi1) <- "seurat_clusters"
CtrlvARKO.combined.Epi1 <- RenameIdents(object = CtrlvARKO.combined.Epi1, '2' = "BE1", '0' = "BE2", '1' = "BE3", '4' = "LE1", '3' = "LE2", '5' = "LE3", '7' = "LE4", '8' = "LE5", '6' = "LE6", '9' = "OE")
DimPlot(CtrlvARKO.combined.Epi1, reduction = "umap", pt.size = 1, label = TRUE)
CtrlvARKO.combined.Epi1[["Epicelltype"]] <- Idents(object = CtrlvARKO.combined.Epi1)
#Dimplot
Idents(object = CtrlvARKO.combined.Epi1) <- "Epicelltype"
DimPlot(CtrlvARKO.combined.Epi1, reduction = "umap", pt.size = 0.3, cols = c("darkcyan", "saddlebrown", "red1", "Purple2", "Gold", "Green3", "hotpink", "#8494FF", "Blue", "Light grey")) 
#Featureplot
DefaultAssay(CtrlvARKO.combined.Epi1) <- "RNA"
FeaturePlot(CtrlvARKO.combined.Epi1, reduction = "umap", features = c("Krt8"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.7, min.cutoff = 0, max.cutoff = 2.5)
FeaturePlot(CtrlvARKO.combined.Epi1, reduction = "umap", features = c("Krt14"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.7, min.cutoff = 0, max.cutoff = 2.5)
FeaturePlot(CtrlvARKO.combined.Epi1, reduction = "umap", features = c("MYC-transgene"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.7, min.cutoff = 0, max.cutoff = 2.5)
#DEGs_Epicelltype
DefaultAssay(CtrlvARKO.combined.Epi1) <- "RNA"
Idents(object = CtrlvARKO.combined.Epi1) <- "Epicelltype"
CtrlvARKO.combined.Epi1 <- ScaleData(CtrlvARKO.combined.Epi1, features = rownames(CtrlvARKO.combined.Epi1))
CtrlvARKO.combined.Epi1.celltype.marker <- FindAllMarkers(CtrlvARKO.combined.Epi1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CtrlvARKO.combined.Epi1.celltype.marker, "CtrlvARKO.combined.Epi1.celltype.marker.csv")
#Dotplot
Idents(object = CtrlvARKO.combined.Epi1) <- "Epicelltype"
DotPlot(CtrlvARKO.combined.Epi1, features = c("Ar", "Pbsn", "MYC-transgene", 
                                              "Pim1", "Cyr61", "Ier3", "Tead1", "Slco2a1", 
                                              "Cp", "Lbp", "Col4a1", "Pltp", "Smoc2",
                                              "Calcb", "Anxa8", "Krt16", "Lypd3", "Adam8",
                                              "Pigr", "Cxcl17", "Pglyrp1", "Cldn10", "Spns2", 
                                              "Crabp1", "Asrgl1", "Defb29", "Gdpd1", "Sbpl",
                                              "Wfdc15b", "Serpinb11", "Fmo1",  "Fmo3", "Folh1",
                                              "Glb1l2", "Fuca1", "Mmp7", "Qsox1", "Tcaf2", 
                                              "Pnliprp1", "Gm37223", "Chn2", "Mt3", "Gsdma",
                                              "Gm26917", "Lars2", "mt-Co3", "AY036118", "Erdr1", 
                                              "Col6a1", "Lox", "Tcf21", "Col5a2", "Col1a1"
), cols = c("light grey", "red")) + RotatedAxis()
#Number of Cells
Idents(object = CtrlvARKO.combined.Epi1) <- "Epicelltype"
CtrlvARKO.combined.Epi1$Epicelltype.stim <- paste(Idents(CtrlvARKO.combined.Epi1), CtrlvARKO.combined.Epi1$stim, sep = "_")
Idents(object = CtrlvARKO.combined.Epi1) <- "Epicelltype.stim"
table(Idents(CtrlvARKO.combined.Epi1))
#Add hMYCtg information
DefaultAssay(CtrlvARKO.combined.Epi1) <- "RNA"
CtrlvARKO.combined.Epi1MYCPos <- subset(x=CtrlvARKO.combined.Epi1,  subset = `MYC-transgene` > 0)
CtrlvARKO.combined.Epi1MYCNeg <- subset(x=CtrlvARKO.combined.Epi1,  subset = `MYC-transgene` == 0)
Idents(object = CtrlvARKO.combined.Epi1MYCPos) <- "MYCPos"
Idents(object = CtrlvARKO.combined.Epi1MYCNeg) <- "MYCNeg"
CtrlvARKO.combined.Epi1MYCPos[["MYCExp"]] <- Idents(object = CtrlvARKO.combined.Epi1MYCPos)
CtrlvARKO.combined.Epi1MYCNeg[["MYCExp"]] <- Idents(object = CtrlvARKO.combined.Epi1MYCNeg)
CtrlvARKO.combined.Epi1MYC <- merge(x = CtrlvARKO.combined.Epi1MYCPos, y = CtrlvARKO.combined.Epi1MYCNeg)
Idents(object = CtrlvARKO.combined.Epi1MYC) <- "MYCExp"
CtrlvARKO.combined.Epi1$MYCExp <- Idents(object = CtrlvARKO.combined.Epi1MYC)
Idents(object = CtrlvARKO.combined.Epi1) <- "MYCExp"
CtrlvARKO.combined.Epi1.MYC <- subset(CtrlvARKO.combined.Epi1, idents = c("MYCPos", "MYCNeg"))
#hMYCtg+ cells in total BE or LE cells
Idents(object = CtrlvARKO.combined.Epi1) <- "Epicelltype.stim"
CtrlvARKO.combined.Epi1$Epicelltype.stim.MYCExp <- paste(Idents(CtrlvARKO.combined.Epi1), CtrlvARKO.combined.Epi1$MYCExp, sep = "_")
Idents(object = CtrlvARKO.combined.Epi1) <- "Epicelltype.stim.MYCExp"
table(Idents(CtrlvARKO.combined.Epi1))
#Dimplot
Idents(object = CtrlvARKO.combined.Epi1) <- "Epicelltype"
DimPlot(CtrlvARKO.combined.Epi1, reduction = "umap", split.by = "stim", pt.size = 0.3, cols = c("darkcyan", "saddlebrown", "red1", "Purple2", "Gold", "Green3", "hotpink", "#8494FF", "Blue", "Light grey")) 
#Featureplot
DefaultAssay(CtrlvARKO.combined.Epi1) <- "RNA"
FeaturePlot(CtrlvARKO.combined.Epi1, reduction = "umap", features = c("Igf1r"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.7, min.cutoff = 0, max.cutoff = "q90")
FeaturePlot(CtrlvARKO.combined.Epi1, reduction = "umap", features = c("Ctnnb1"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.7, min.cutoff = 0, max.cutoff = 2.5)
FeaturePlot(CtrlvARKO.combined.Epi1, reduction = "umap", features = c("Ccnd1"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.7, min.cutoff = 0, max.cutoff = 2.5)
FeaturePlot(CtrlvARKO.combined.Epi1, reduction = "umap", features = c("Cd44"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.7, min.cutoff = 0, max.cutoff = 2.5)
####subclustering MYC+ BE####
Idents(object = CtrlvARKO.combined.Epi) <- "MYCExp"
CtrlvARKO.combined.Epi$MYCExp.Epicelltype <- paste(Idents(CtrlvARKO.combined.Epi), CtrlvARKO.combined.Epi$Epicelltype, sep = "_")
Idents(object = CtrlvARKO.combined.Epi) <- "MYCExp.Epicelltype"
CtrlvARKO.combined.MYCBE <- subset(CtrlvARKO.combined.Epi, idents = c("MYCPos_BE1", "MYCPos_BE2", "MYCPos_BE3"))
#DEGs_hMYCtg+BE in ctrl vs hMYCtg+BE in ARKO
DefaultAssay(CtrlvARKO.combined.MYCBE) <- "RNA"
Idents(object = CtrlvARKO.combined.MYCBE) <- "stim"
CtrlvARKO.combined.MYCBE <- ScaleData(CtrlvARKO.combined.MYCBE, features = rownames(CtrlvARKO.combined.MYCBE))
CtrlvARKO.combined.MYCBE.GSEA.Markers <- FindMarkers(CtrlvARKO.combined.MYCBE, ident.1 = "Ctrl", ident.2 = "ARKO", min.pct = 0.01, logfc.threshold = 0.01)
write.csv(CtrlvARKO.combined.MYCBE.GSEA.Markers, "CtrlvARKO.combined.MYCBE.GSEA.Markers.csv")
#Boxplot Generation
boxdata = FetchData(CtrlvARKO.combined.MYCBE, c("stim", "Igf1r", "Hras", "Irs2", "Akt1", "Ctnnb1", "Ccnd1", "Cd44", "Tcf7l2"))
tail(boxdata,9)
ggplot(boxdata, aes(x=stim, y=Igf1r, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
ggplot(boxdata, aes(x=stim, y=Hras, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
ggplot(boxdata, aes(x=stim, y=Irs2, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
ggplot(boxdata, aes(x=stim, y=Akt1, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
ggplot(boxdata, aes(x=stim, y=Ctnnb1, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
ggplot(boxdata, aes(x=stim, y=Ccnd1, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
ggplot(boxdata, aes(x=stim, y=Cd44, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
ggplot(boxdata, aes(x=stim, y=Tcf7l2, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
#hMYCtg+BE_in_Ctrl
Idents(object = CtrlvARKO.combined.Epi1) <- "Epicelltype.stim"
DimPlot(CtrlvARKO.combined.Epi1, reduction = "umap", pt.size = 0.3)
Ctrl.combined.BE <- subset(CtrlvARKO.combined.Epi1, idents = c("BE1_Ctrl", "BE2_Ctrl", "BE3_Ctrl"))
Idents(object = Ctrl.combined.BE) <- "MYCExp"
GOI <- c('Igf1r','Hras','Irs2', 'Akt1', 'Ctnnb1', 'Ccnd1', 'Cd44', 'Tcf7l2')  
GOI_index <- is.element(rownames(Ctrl.combined.BE),GOI)
Cell_index <- is.element(Idents(Ctrl.combined.BE),c('MYCPos'))
expr_GOI <- Ctrl.combined.BE@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- Ctrl.combined.BE@assays$RNA@counts[GOI_index,Cell_index]
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE)
#hMYCtg+BE_in_ARKO
Idents(object = CtrlvARKO.combined.Epi1) <- "Epicelltype.stim"
DimPlot(CtrlvARKO.combined.Epi1, reduction = "umap", pt.size = 0.3)
ARKO.combined.BE <- subset(CtrlvARKO.combined.Epi1, idents = c("BE1_ARKO", "BE2_ARKO", "BE3_ARKO"))
Idents(object = ARKO.combined.BE) <- "MYCExp"
GOI <- c('Igf1r','Hras','Irs2', 'Akt1', 'Ctnnb1', 'Ccnd1', 'Cd44', 'Tcf7l2')  
GOI_index <- is.element(rownames(ARKO.combined.BE),GOI)
Cell_index <- is.element(Idents(ARKO.combined.BE),c('MYCPos'))
expr_GOI <- ARKO.combined.BE@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- ARKO.combined.BE@assays$RNA@counts[GOI_index,Cell_index]
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE)
####Subclustering Stro ####
Idents(object = CtrlvARKO.combined) <- "celltype"
CtrlvARKO.combined.Stro <- subset(CtrlvARKO.combined, idents = c("FB", "SM", "Pericyte", "Leukocyte", "Endothelia", "Glia"))
DefaultAssay(CtrlvARKO.combined.Stro) <- "integrated"
#Run the standard workflow for visualization and clustering
CtrlvARKO.combined.Stro <- ScaleData(CtrlvARKO.combined.Stro, verbose = FALSE)
CtrlvARKO.combined.Stro <- RunPCA(CtrlvARKO.combined.Stro, npcs = 30, verbose = FALSE)
ElbowPlot(CtrlvARKO.combined.Stro, ndims = 50)
#Umap and Clustering
CtrlvARKO.combined.Stro <- FindNeighbors(CtrlvARKO.combined.Stro, reduction = "pca", dims = 1:15)
CtrlvARKO.combined.Stro <- FindClusters(CtrlvARKO.combined.Stro, resolution = 0.5)
CtrlvARKO.combined.Stro <- RunUMAP(CtrlvARKO.combined.Stro, reduction = "pca", dims = 1:15)
DimPlot(CtrlvARKO.combined.Stro, reduction = "umap", label = TRUE)
#DEGs_celltype
DefaultAssay(CtrlvARKO.combined.Stro) <- "RNA"
Idents(object = CtrlvARKO.combined.Stro) <- "seurat_clusters"
CtrlvARKO.combined.Stro <- ScaleData(CtrlvARKO.combined.Stro, features = rownames(CtrlvARKO.combined.Stro))
CtrlvARKO.combined.Stro.marker <- FindAllMarkers(CtrlvARKO.combined.Stro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CtrlvARKO.combined.Stro.marker, "CtrlvARKO.combined.Stro.marker.csv")
#Rename
CtrlvARKO.combined.Stro <- RenameIdents(object = CtrlvARKO.combined.Stro, '3' = "FB1", '2' = "FB2",
                                        '0' = "FB3", '1' = "FB4", '4' = "FB5", '8' = "SM", 
                                        '12' = "Peri", '6' = "Endo", '13' = "Glia", '5' = "Leu",
                                        '11' = "Leu", '10' = "Leu", '7' = "Leu", '9' = "OS")
CtrlvARKO.combined.Stro[["Strocelltype"]] <- Idents(object = CtrlvARKO.combined.Stro)
####Subclustering FBSM####
Idents(object = CtrlvARKO.combined.Stro) <- "Strocelltype"
DimPlot(CtrlvARKO.combined.Stro, reduction = "umap", pt.size = 0.3, cols = c("red", "purple", "blue", "green", "orange", "brown", "light grey", "light grey", "light grey", "light grey", "light grey"))
CtrlvARKO.combined.FBSM <- subset(CtrlvARKO.combined.Stro, idents = c("FB1", "FB2", "FB3", "FB4", "FB5", "SM"))
CtrlvARKO.combined.FBSM <- RunUMAP(CtrlvARKO.combined.FBSM, reduction = "pca", dims = 1:15)
DimPlot(CtrlvARKO.combined.FBSM, reduction = "umap", label = TRUE)
#Rename
Idents(object = CtrlvARKO.combined.FBSM) <- "Strocelltype"
DimPlot(CtrlvARKO.combined.FBSM, reduction = "umap", label = TRUE)
CtrlvARKO.combined.FBSM <- RenameIdents(object = CtrlvARKO.combined.FBSM, 'FB1' = "FB1", 'FB2' = "FB2",
                                        'FB3' = "FB3", 'FB4' = "FB4", 'FB5' = "FB5", 'SM' = "SM")
CtrlvARKO.combined.FBSM[["Strocelltype"]] <- Idents(object = CtrlvARKO.combined.FBSM)
#Umap
Idents(object = CtrlvARKO.combined.FBSM) <- "Strocelltype"
DimPlot(CtrlvARKO.combined.FBSM, reduction = "umap", pt.size = 0.3, cols = c("orange", "purple", "blue", "red", "green", "brown"))
DimPlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", pt.size = 0.3, cols = c("orange", "purple", "blue", "red", "green", "brown"))
#Featureplots
DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", features = c("EGFP"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = 0, max.cutoff = 2.5)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", features = c("Ar"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = 0, max.cutoff = 2.5)
#DEGs_Strocelltype
DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
Idents(object = CtrlvARKO.combined.FBSM) <- "Strocelltype"
CtrlvARKO.combined.FBSM <- ScaleData(CtrlvARKO.combined.FBSM, features = rownames(CtrlvARKO.combined.FBSM))
CtrlvARKO.combined.FBSM.celltype.marker <- FindAllMarkers(CtrlvARKO.combined.FBSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CtrlvARKO.combined.FBSM.celltype.marker, "CtrlvARKO.combined.FBSM.celltype.marker.csv")
#Dotplot
Idents(object = CtrlvARKO.combined.FBSM) <- "Strocelltype"
tiff(file = "CtrlvARKO.combined.FBSM Celltype dotplot.tiff", width = 12, height = 3, units = "in", compression = "lzw", res = 800)
DotPlot(CtrlvARKO.combined.FBSM, features = c("Ar", "EGFP", "Gli1", 
                                              "Cxcl14", "F2r", "Itgbl1", "Hsd11b1", "Isg15",
                                              "Pcolce2", "Cd248", "Sfrp4", "Tmem100", "Mest",
                                              "Jun", "Klf2", "Hspa1a", "Nr1d1", "Tgfbi",
                                              "Ifi205", "Thbd", "Ccl11", "Ptx3", "Sult1e1",
                                              "Kdm6b", "Maff", "Pim1", "Egr1", "Nr4a1",
                                              "Acta2", "Tagln", "Myl9", "Actg2", "Myh11"
), cols = c("light grey", "red")) + RotatedAxis()
dev.off()
#Number of Cells
Idents(object = CtrlvARKO.combined.FBSM) <- "Strocelltype"
CtrlvARKO.combined.FBSM$Strocelltype.stim <- paste(Idents(CtrlvARKO.combined.FBSM), CtrlvARKO.combined.FBSM$stim, sep = "_")
Idents(object = CtrlvARKO.combined.FBSM) <- "Strocelltype.stim"
table(Idents(CtrlvARKO.combined.FBSM))
#Violin plots
Idents(object = CtrlvARKO.combined.FBSM) <- "Strocelltype"
VlnPlot(CtrlvARKO.combined.FBSM, features = "EGFP", pt.size = 0, split.by = "stim", split.plot = TRUE, cols = c("#3399FF", "#E06666"))
VlnPlot(CtrlvARKO.combined.FBSM, features = "Ar", pt.size = 0, split.by = "stim", split.plot = TRUE, cols = c("#3399FF", "#E06666"))
VlnPlot(CtrlvARKO.combined.FBSM, features = "Igfbp3", pt.size = 0, split.by = "stim", split.plot = TRUE, cols = c("#3399FF", "#E06666"))
#Featureplots
DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", features = c("Igfbp3"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.3, min.cutoff = 0, max.cutoff = 2.5)
#Coexpression plots
Idents(object = CtrlvARKO.combined.FBSM) <- "stim"
CtrlvARKO.combined.FBSM.Ctrl <- subset(CtrlvARKO.combined.FBSM, idents = c("Ctrl"))
DefaultAssay(CtrlvARKO.combined.FBSM.Ctrl) <- "RNA"
FeaturePlot(CtrlvARKO.combined.FBSM.Ctrl, reduction = "umap", features = c("EGFP", "Igfbp3"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
FeaturePlot(CtrlvARKO.combined.FBSM.Ctrl, reduction = "umap", features = c("Ar", "Igfbp3"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
Idents(object = CtrlvARKO.combined.FBSM) <- "stim"
CtrlvARKO.combined.FBSM.ARKO <- subset(CtrlvARKO.combined.FBSM, idents = c("ARKO"))
DefaultAssay(CtrlvARKO.combined.FBSM.ARKO) <- "RNA"
FeaturePlot(CtrlvARKO.combined.FBSM.ARKO, reduction = "umap", features = c("EGFP", "Igfbp3"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
FeaturePlot(CtrlvARKO.combined.FBSM.ARKO, reduction = "umap", features = c("Ar", "Igfbp3"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
####Subclustering EGFP+ FB####
#Add EGFP info
DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
CtrlvARKO.combined.FBSM.EGFPPos <- subset(x=CtrlvARKO.combined.FBSM,  subset = `EGFP` > 0)
CtrlvARKO.combined.FBSM.EGFPNeg <- subset(x=CtrlvARKO.combined.FBSM,  subset = `EGFP` == 0)
Idents(object = CtrlvARKO.combined.FBSM.EGFPPos) <- "EGFPPos"
Idents(object = CtrlvARKO.combined.FBSM.EGFPNeg) <- "EGFPNeg"
CtrlvARKO.combined.FBSM.EGFPPos[["EGFPExp"]] <- Idents(object = CtrlvARKO.combined.FBSM.EGFPPos)
CtrlvARKO.combined.FBSM.EGFPNeg[["EGFPExp"]] <- Idents(object = CtrlvARKO.combined.FBSM.EGFPNeg)
CtrlvARKO.combined.FBSM.EGFP <- merge(x = CtrlvARKO.combined.FBSM.EGFPPos, y = CtrlvARKO.combined.FBSM.EGFPNeg)
Idents(object = CtrlvARKO.combined.FBSM.EGFP) <- "EGFPExp"
CtrlvARKO.combined.FBSM$EGFPExp <- Idents(object = CtrlvARKO.combined.FBSM.EGFP)
#Subset EGFP+FB
Idents(object = CtrlvARKO.combined.FBSM) <- "EGFPExp"
combined.FBSM.EGFP <- subset(CtrlvARKO.combined.FBSM, idents = c("EGFPPos"))
Idents(object = combined.FBSM.EGFP) <- "celltype"
combined.FB.EGFP <- subset(combined.FBSM.EGFP, idents = c("FB"))
#Violinplot
DefaultAssay(combined.FB.EGFP) <- "RNA"
VlnPlot(combined.FB.EGFP, features = "Igfbp3", group.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"))
####CAF analysis####
#Heatmap_HiMYC_FBvHiMYC-ARKO_FB
Idents(object = CtrlvARKO.combined.FBSM) <- "Strocelltype"
CtrlvARKO.combined.FB <- subset(CtrlvARKO.combined.FBSM, idents = c("FB1", "FB2", "FB3", "FB4", "FB5"))
Idents(object = CtrlvARKO.combined.FB) <- "stim"
CtrlvARKO.combined.FB <- ScaleData(CtrlvARKO.combined.FB, features = rownames(CtrlvARKO.combined.FB))
DoHeatmap(CtrlvARKO.combined.FB, features = c("Ar", "Vim", "Pdgfrb", "Fap", "S100a4", "Twist1", "Foxf1", "Il11",  "Sox9", "Cxcl10")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
#Heatmap_Strocelltype_Both
Idents(object = CtrlvARKO.combined.FB) <- "Strocelltype"
CtrlvARKO.combined.FB <- ScaleData(CtrlvARKO.combined.FB, features = rownames(CtrlvARKO.combined.FB))
DoHeatmap(CtrlvARKO.combined.FB, features = c("Ar", "Vim", "Pdgfrb", "Fap", "S100a4", "Twist1", "Foxf1", "Il11",  "Sox9", "Cxcl10")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
#Heatmap_Strocelltype_Ctrl
Idents(object = CtrlvARKO.combined.FB) <- "stim"
Ctrl.combined.FB <- subset(CtrlvARKO.combined.FB, idents = c("Ctrl"))
Idents(object = Ctrl.combined.FB) <- "Strocelltype"
Ctrl.combined.FB <- ScaleData(Ctrl.combined.FB, features = rownames(Ctrl.combined.FB))
DoHeatmap(Ctrl.combined.FB, features = c("Ar", "Vim", "Pdgfrb", "Fap", "S100a4", "Twist1", "Foxf1", "Il11",  "Sox9", "Cxcl10")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
#Heatmap_Himyc_mGFP+Ar+FB1_vs_HiMYC-ARKO_mGFP+Ar-FB1
Idents(object = CtrlvARKO.combined.FBSM) <- "ArExp"
CtrlvARKO.combined.FBSM$ArExp.EGFPExp <- paste(Idents(CtrlvARKO.combined.FBSM), CtrlvARKO.combined.FBSM$EGFPExp, sep = "_")
Idents(object = CtrlvARKO.combined.FBSM) <- "ArExp.EGFPExp"
CtrlvARKO.combined.FBSM$ArExp.EGFPExp.Strocelltype <- paste(Idents(CtrlvARKO.combined.FBSM), CtrlvARKO.combined.FBSM$Strocelltype, sep = "_")
Idents(object = CtrlvARKO.combined.FBSM) <- "ArExp.EGFPExp.Strocelltype"
CtrlvARKO.combined.FBSM$ArExp.EGFPExp.Strocelltype.stim <- paste(Idents(CtrlvARKO.combined.FBSM), CtrlvARKO.combined.FBSM$stim, sep = "_")
Idents(object = CtrlvARKO.combined.FBSM) <- "ArExp.EGFPExp.Strocelltype.stim"
DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
Idents(object = CtrlvARKO.combined.FBSM) <- "ArExp.EGFPExp.Strocelltype.stim"
CtrlvARKO.combined.FB1 <- subset(CtrlvARKO.combined.FBSM, idents = c("ArPos_EGFPPos_FB1_Ctrl", "ArNeg_EGFPPos_FB1_ARKO"))
DefaultAssay(CtrlvARKO.combined.FB1) <- "RNA"
CtrlvARKO.combined.FB1 <- ScaleData(CtrlvARKO.combined.FB1, features = rownames(CtrlvARKO.combined.FB1))
DoHeatmap(CtrlvARKO.combined.FB1, features = c("Ar", "Fkbp5", "Aldh1a1", "Ly6a", "Pbsn", "Wnt2", "Wnt6", "Fos", "Fosb", "Jun", "Junb", "Egr1", "Serpine1", "Cdkn1a", "Dnajb1", "Bcl2", "Igfbp3", "Igfbp7")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
#DEGs_Himyc_mGFP+Ar+FB1_vs_HiMYC-ARKO_mGFP+Ar-FB1
DefaultAssay(CtrlvARKO.combined.FB1) <- "RNA"
Idents(object = CtrlvARKO.combined.FB1) <- "stim"
CtrlvARKO.combined.FB1 <- ScaleData(CtrlvARKO.combined.FB1, features = rownames(CtrlvARKO.combined.FB1))
CtrlvARKO.combined.FB1.GSEA.Markers <- FindMarkers(CtrlvARKO.combined.FB1, ident.1 = "Ctrl", ident.2 = "ARKO", min.pct = 0.01, logfc.threshold = 0.01)
write.csv(CtrlvARKO.combined.FB1.GSEA.Markers, "CtrlvARKO.combined.FB1.GSEA.Markers.csv")
#Violinplots
Idents(object = CtrlvARKO.combined.FB) <- "Strocelltype"
VlnPlot(CtrlvARKO.combined.FB, features = "Pdgfrb", group.by = "Strocelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
VlnPlot(CtrlvARKO.combined.FB, features = "Twist1", group.by = "Strocelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
VlnPlot(CtrlvARKO.combined.FB, features = "Foxf1", group.by = "Strocelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
VlnPlot(CtrlvARKO.combined.FB, features = "Il11", group.by = "Strocelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
VlnPlot(CtrlvARKO.combined.FB, features = "Cxcl10", group.by = "Strocelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
VlnPlot(CtrlvARKO.combined.FB, features = "Sox9", group.by = "Strocelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)