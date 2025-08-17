# NSCLC SINGLE-CELL RNA-SEQ ANALYSIS PIPELINE

# -------------------------------------
# STEP 1: LOAD LIBRARIES
# -------------------------------------
library(Seurat)
library(tidyverse)

# -------------------------------------
# STEP 2: LOAD 10X DATA
# -------------------------------------
h5_path <- file.choose()  
data <- Read10X_h5(h5_path)
cts <- data$`Gene Expression`


# -------------------------------------
# STEP 3: CREATE SEURAT OBJECT
# -------------------------------------
seurat_obj <- CreateSeuratObject(counts = cts, project = "NSCLC", min.cells = 3, min.features = 200)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# -------------------------------------
# STEP 4: QUALITY CONTROL PLOTS
# -------------------------------------
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")

# -------------------------------------
# STEP 5: FILTERING
# -------------------------------------
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# -------------------------------------
# STEP 6: NORMALIZATION
# -------------------------------------
seurat_obj <- NormalizeData(seurat_obj)

# -------------------------------------
# STEP 7: HIGHLY VARIABLE GENES
# -------------------------------------
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seurat_obj), 10)
VariableFeaturePlot(seurat_obj) %>% LabelPoints(points = top10)

# -------------------------------------
# STEP 8: SCALING
# -------------------------------------
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

# -------------------------------------
# STEP 9: PCA
# -------------------------------------
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
ElbowPlot(seurat_obj)

# -------------------------------------
# STEP 10: CLUSTERING + UMAP
# -------------------------------------
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:14)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:14)
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

# -------------------------------------
# STEP 11: MARKER GENES
# -------------------------------------
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()

# -------------------------------------
# STEP 12: SINGLE R ANNOTATION
# -------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("SingleR", "celldex"), ask = FALSE, update = FALSE)

library(SingleR)
library(celldex)
ref <- celldex::HumanPrimaryCellAtlasData()
res <- SingleR(test = GetAssayData(seurat_obj, slot = "data"), ref = ref, labels = ref$label.main)
seurat_obj$SingleR.labels <- res$labels[match(colnames(seurat_obj), rownames(res))]
DimPlot(seurat_obj, group.by = "SingleR.labels", label = TRUE, repel = TRUE) + NoLegend()

# -------------------------------------
# STEP 13: ENRICHMENT (EXAMPLE: T CELLS)
# -------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)

Idents(seurat_obj) <- "SingleR.labels"
t_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  filter(cluster == "T_cells" & p_val_adj < 0.05)

gene.df <- bitr(t_markers$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene = gene.df$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                ont = "BP", pvalueCutoff = 0.05)
barplot(ego, showCategory = 10, title = "GO: Biological Processes (T cells)")

# -------------------------------------
# STEP 14: MONOCLE3 PSEUDOTIME
# -------------------------------------
library(SeuratWrappers)
library(monocle3)

cds <- as.cell_data_set(seurat_obj)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)

plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = FALSE,
           label_leaves = TRUE, label_branch_points = TRUE)

# -------------------------------------
# STEP 15: SAVE RESULTS
# -------------------------------------
saveRDS(seurat_obj, file = "results/NSCLC_seurat_obj.rds")
saveRDS(cds, file = "results/NSCLC_monocle3_cds.rds")
write.csv(seurat_obj@meta.data, "results/NSCLC_metadata_with_SingleR.csv")

