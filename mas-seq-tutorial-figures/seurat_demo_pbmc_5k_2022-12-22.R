library(Seurat)
library(ggplot2)
library(dplyr)
theme_set(theme_bw(base_size=18))

matrix <- ReadMtx("genes_seurat/matrix.mtx",
                  features = "genes_seurat/genes.tsv",
                  cells = "genes_seurat/barcodes.tsv",
                  feature.column = 2)

seurat_df <- CreateSeuratObject(counts = matrix, project="PBMC5k", 
                                min.cells = 3, min.features = 200)

# Look at mitochondrial reads and ribosomal reads. They should already be filtered
# in SMRT Link
seurat_df[["percent.mt"]] <- PercentageFeatureSet(seurat_df, pattern = "^MT-")
seurat_df[["percent.ribo"]] <- PercentageFeatureSet(seurat_df, pattern = "^RP[SL]")

VlnPlot(seurat_df, 
        features = c("nFeature_RNA", "nCount_RNA", 
                     "percent.mt", "percent.ribo"), 
        ncol = 2)

# Correlation between number of counts to number of features 
FeatureScatter(seurat_df, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# Normalize data
seurat_df <- NormalizeData(seurat_df, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_df <- FindVariableFeatures(seurat_df, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_df), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_df)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
# PCA
all.genes <- rownames(seurat_df)
# Scaling with all genes to make visualization easier
seurat_df <- ScaleData(seurat_df, features = all.genes)
# PCA with 2000 variable features
seurat_df <- RunPCA(seurat_df, features = VariableFeatures(object = seurat_df))
# Look at important PC:
ElbowPlot(seurat_df)

# Choose cut-off based on elbow plot
seurat_df <- FindNeighbors(seurat_df, dims = 1:13)
seurat_df <- FindClusters(seurat_df, resolution = 0.5)
# UMAP
seurat_df <- RunUMAP(seurat_df, dims = 1:13)
plot_umap <- DimPlot(seurat_df, reduction = "umap")
ggsave("UMAP_5kCells.pdf", plot_umap, width=8, height=6, useDingbats=FALSE)

# Find marker genes for each cluster
seurat.marker <-FindAllMarkers(seurat_df, only.pos = TRUE,
                               min.pct = 0.25, logfc.threshold = 0.25)

# Top 5 marker genes in each cluster
ordered_markers <- seurat.marker |>
  group_by(cluster) |>
  # Top 3 markers in each cluster
  slice_max(n = 5, order_by = avg_log2FC)

# Plot marker genes with heatmap or dotplot
DoHeatmap(seurat_df, features = ordered_markers$gene) + 
  NoLegend()
plot_dot <- DotPlot(seurat_df, features = unique(ordered_markers$gene)) +
  RotatedAxis()
ggsave("Dotplot_5kCells.pdf", plot_dot, width=20, height=9, useDingbats=FALSE)

# PBMC marker genes from Seurat tutorial
pbmc_markers <- c("IL7R", "CCR7", "CD14", "LYZ", "S100A4", "MS4A1",
                  "CD8A", "FCGR3A", "MS4A7", "GNLY", "NKG7",
                  "FCER1A", "CST3", "PPBP", "CD3E")
DoHeatmap(seurat_df, features = pbmc_markers)
# Rename identity based on the heatmap above
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD8 T", "CD14+ Mono", "B", "FCGR3A+ Mono",
                     "NK", "", "Platelet", "DC", " ")
names(new.cluster.ids) <- levels(seurat_df)
seurat_df <- RenameIdents(seurat_df, new.cluster.ids)

# UMAP plot with the cell types
plot_with_cell <- DimPlot(seurat_df, reduction = "umap", 
                          label = TRUE, pt.size = 0.5) +
  NoLegend()
ggsave("UMAP_5kCells_withCells.pdf", plot_with_cell, width=8, height=6, useDingbats=FALSE)

