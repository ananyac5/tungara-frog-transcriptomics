#--------------------------------------------------------------
# This is the code for the scType with just the default dataset 
# using brain cell markers, not that important
#--------------------------------------------------------------

# Download default ScType marker database
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"

gs_list <- gene_sets_prepare(
  db_,
  "Brain"
)

# Run ScType
es.max <- sctype_score(
  scRNAseqData = combined[["RNA"]]$scale.data,
  scaled = TRUE,
  gs = gs_list$gs_positive,
  gs2 = gs_list$gs_negative
)

# Assign best label per cluster
sctype_top <- do.call("rbind", lapply(unique(combined$seurat_clusters), function(cl) {
  
  cells <- rownames(combined@meta.data)[combined$seurat_clusters == cl]
  
  scores <- sort(
    rowSums(es.max[, cells, drop = FALSE]),
    decreasing = TRUE
  )
  
  data.frame(
    cluster = cl,
    cell_type = names(scores)[1],
    score = scores[1]
  )
}))

# Add labels to Seurat object
combined$sctype_label <- sctype_top$cell_type[
  match(combined$seurat_clusters, sctype_top$cluster)
]

# Plot
DimPlot(
  combined,
  group.by = "sctype_label",
  label = TRUE,
  repel = TRUE
)

# View assignments
sctype_top
table(combined$sctype_label)


#------------------------------------------------------------------------
# Integration approach
#------------------------------------------------------------------------
# Make a seperate combined object for integration
combined_int <- merge(
  x = control,
  y = list(condW, condWC),
  add.cell.ids = c("GN", "W", "WC"),
  project = "frog_conditions_integrated"
)

# normal preprocessing before integration
combined_int <- NormalizeData(combined_int)
combined_int <- FindVariableFeatures(combined_int)
combined_int <- ScaleData(combined_int)
combined_int <- RunPCA(combined_int)

# Seurat v5 integration
combined_int <- IntegrateLayers(
  object = combined_int,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)

# re-join RNA layers after integration
combined_int[["RNA"]] <- JoinLayers(combined_int[["RNA"]])

# clustering using integrated reduction
combined_int <- FindNeighbors(
  combined_int,
  reduction = "integrated.cca",
  dims = 1:20
)

combined_int <- FindClusters(
  combined_int,
  resolution = 0.6
)

combined_int <- RunUMAP(
  combined_int,
  reduction = "integrated.cca",
  dims = 1:20
)

#adding back the conditions labels
combined_int$condition <- case_when(
  grepl("^GN_", colnames(combined_int)) ~ "GN",
  grepl("^W_", colnames(combined_int)) ~ "W",
  grepl("^WC_", colnames(combined_int)) ~ "WC"
)

# plots
DimPlot(combined_int, reduction = "umap", label = TRUE)
DimPlot(combined_int, reduction = "umap", group.by = "condition", alpha = 0.5)

#------------------------------------------
# Finding markers from the integrated data
#------------------------------------------
DefaultAssay(combined_int) <- "RNA"
Idents(combined_int) <- "seurat_clusters"

markers_int <- FindAllMarkers(
  combined_int,
  only.pos = TRUE,
  logfc.threshold = 0.693
)

top_markers_int <- markers_int %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 8)

qsave(combined_int, "combined_int.qs")
qsave(markers_int, "markers_int.qs")

combined_int <- qread("combined_int.qs")
markers_int <- qread("markers_int.rds")

#-------------------------------------------
#Testing out different clustering strategies
#-------------------------------------------
combined_int_test <- combined_int

combined_int_test <- FindNeighbors(
  combined_int_test,
  reduction = "integrated.cca",
  dims = 1:45,
  graph.name = c("test_nn", "test_snn")
)

combined_int_test <- FindClusters(
  combined_int_test,
  graph.name = "test_snn",
  resolution = 0.6,
  cluster.name = "clusters_30PC_0.6"
)

combined_int_test <- RunUMAP(
  combined_int_test,
  reduction = "integrated.cca",
  dims = 1:45,
  reduction.name = "umap_30PC_0.6"
)

DimPlot(combined_int_test, reduction = "umap_30PC_0.6", group.by = "clusters_30PC_0.6", label = TRUE)

table(combined_int_test$clusters_30PC_0.6)

qsave(combined_int_test, "combined_int_test.qs")
combined_int_test <- qread("combined_int_test.qs")
#----------------------------------------------
#scType & checking FindAllMarkers() for this dimension 0.6 
es.max <- sctype_score(
  scRNAseqData = combined_int_test[["RNA"]]$scale.data,
  scaled = TRUE,
  gs = gs_list$gs_positive,
  gs2 = NULL
)

sctype_top <- do.call("rbind", lapply(unique(combined_int_test$clusters_30PC_0.6), function(cl) {
  cells <- rownames(combined_int_test@meta.data)[combined_int_test$clusters_30PC_0.6 == cl]
  scores <- sort(rowSums(es.max[, cells, drop = FALSE]), decreasing = TRUE)
  data.frame(cluster = cl, cell_type = unname(names(scores)[1]), score = unname(scores[1]))
}))

rownames(sctype_top) <- NULL

combined_int_test$sctype_label_30PC_0.6 <- sctype_top$cell_type[
  match(combined_int_test$clusters_30PC_0.6, sctype_top$cluster)
]

DimPlot(combined_int_test, reduction = "umap_30PC_0.6", group.by = "sctype_label_30PC_0.6", label = F)
#------------------------------------------------
DefaultAssay(combined_int_test) <- "RNA"
Idents(combined_int_test) <- "clusters_30PC_0.6"

markers_30PC_0.6 <- FindAllMarkers(
  combined_int_test,
  only.pos = TRUE,
  logfc.threshold = 0.693 #log2Fc
)

#----------------------------------------
#Running scType
#----------------------------------------
es.max <- sctype_score(
  scRNAseqData = combined_int[["RNA"]]$scale.data,
  scaled = TRUE,
  gs = gs_list$gs_positive,
  gs2 = NULL
)

sctype_top = do.call("rbind", lapply(unique(combined_int$seurat_clusters), function(cl) {
  cells = rownames(combined_int@meta.data)[combined_int$seurat_clusters == cl]
  scores = sort(rowSums(es.max[, cells, drop = FALSE]), decreasing = TRUE)
  
  data.frame(
    cluster = cl,
    cell_type = names(scores)[1],
    score = scores[1]
  )
}))

combined_int$sctype_label <- sctype_top$cell_type[
  match(combined_int$seurat_clusters, sctype_top$cluster)
]

DimPlot(combined_int, group.by = "sctype_label", label = T)
table(sctype_top$cell_type)

#=========================================================================
# MAPMYCELL convertion
#=========================================================================
install.packages("anndata")
library(anndata)
library(Matrix)

# Get raw RNA counts: genes x cells
counts <- GetAssayData(
  combined_int,
  assay = "RNA",
  layer = "counts"
)

# MapMyCells needs cells x genes
counts_t <- Matrix::t(counts)

# Create AnnData object
adata <- AnnData(
  X = counts_t,
  obs = data.frame(
    row.names = rownames(counts_t)
  ),
  var = data.frame(
    row.names = colnames(counts_t)
  )
)

# Save as h5ad
adata$write_h5ad(
  "tungara_brain_raw_counts.h5ad",
  compression = "gzip"
)

#------------------------------------------------------
# Looking at the MapMyCell results for exploration
#------------------------------------------------------
#reading in the mapmycell csv result
mapmycells <- read.csv(
  "tungara_brain_10xWholeMouse_Correlation.csv",
  skip = 4,
  check.names = FALSE
)

#creating a seperate obj to keep the original annotations
mapmycell_cor <- combined_int

#matching mapmycell ids
match_index <- match(
  colnames(mapmycell_cor),
  mapmycells$cell_id
)

#taxonomy labels
mapmycell_cor$mmc_class <- mapmycells$class_name[match_index]
mapmycell_cor$mmc_subclass <- mapmycells$subclass_name[match_index]
mapmycell_cor$mmc_supertype <- mapmycells$supertype_name[match_index]
mapmycell_cor$mmc_cluster <- mapmycells$cluster_name[match_index]

#correlation scores
mapmycell_cor$mmc_class_cor <- mapmycells$class_correlation_coefficient[match_index]
mapmycell_cor$mmc_subclass_cor <- mapmycells$subclass_correlation_coefficient[match_index]
mapmycell_cor$mmc_supertype_cor <- mapmycells$supertype_correlation_coefficient[match_index]
mapmycell_cor$mmc_cluster_cor <- mapmycells$cluster_correlation_coefficient[match_index]

#dim plot broader class
DimPlot(mapmycell_cor, group.by = "mmc_class", label = FALSE)

#comparison against my scType annotation
table(ScType = mapmycell_cor$sctype_label, MapMyCells = mapmycell_cor$mmc_class)

#collapsing the mouse classes
mapmycell_cor$mmc_broad <- dplyr::case_when(
  grepl("Glut", mapmycell_cor$mmc_class) ~ "Glutamatergic neuron",
  grepl("GABA", mapmycell_cor$mmc_class) ~ "GABAergic neuron",
  grepl("Astro-Epen", mapmycell_cor$mmc_class) ~ "Astrocyte/Ependymal",
  grepl("OPC-Oligo", mapmycell_cor$mmc_class) ~ "Oligodendrocyte lineage",
  grepl("Vascular", mapmycell_cor$mmc_class) ~ "Vascular",
  grepl("Immune", mapmycell_cor$mmc_class) ~ "Immune",
  TRUE ~ "Other")

DimPlot(mapmycell_cor, group.by = "mmc_broad", label = FALSE)

#looking at scType simalarities by percentage
comparison_broad <- table(
  ScType = mapmycell_cor$sctype_label,
  MapMyCells = mapmycell_cor$mmc_broad)

comparison_percent <- prop.table(comparison_broad, margin = 1) * 100
round(comparison_percent, 1)
#----------------------------------------------------------
# Looking at hierarchical clustering instead or correlation
#----------------------------------------------------------
# reading in the MapMyCells hierarchical CSV result
mapmycells_h <- read.csv(
  "tungara_brain_10xWholeMouseBrain_Hierarchical.csv",
  skip = 4,
  check.names = FALSE
)

# creating a separate object to keep the original annotations
mapmycell_hier <- combined_int

# matching MapMyCells cell IDs to Seurat cell barcodes
match_index <- match(colnames(mapmycell_hier), mapmycells_h$cell_id)

# taxonomy labels
mapmycell_hier$mmc_class <- mapmycells_h$class_name[match_index]
mapmycell_hier$mmc_subclass <- mapmycells_h$subclass_name[match_index]
mapmycell_hier$mmc_supertype <- mapmycells_h$supertype_name[match_index]
mapmycell_hier$mmc_cluster <- mapmycells_h$cluster_name[match_index]

# DimPlot for broader class
DimPlot(mapmycell_hier, group.by = "mmc_class", label = FALSE)

# comparison against ScType annotation
table(
  ScType = mapmycell_hier$sctype_label,
  MapMyCells = mapmycell_hier$mmc_class)

# collapsing the mouse classes into broad categories
mapmycell_hier$mmc_broad <- dplyr::case_when(
  grepl("Glut", mapmycell_hier$mmc_class) ~ "Glutamatergic neuron",
  grepl("GABA", mapmycell_hier$mmc_class) ~ "GABAergic neuron",
  grepl("Astro-Epen", mapmycell_hier$mmc_class) ~ "Astrocyte/Ependymal",
  grepl("OPC-Oligo", mapmycell_hier$mmc_class) ~ "Oligodendrocyte lineage",
  grepl("Vascular", mapmycell_hier$mmc_class) ~ "Vascular",
  grepl("Immune", mapmycell_hier$mmc_class) ~ "Immune",
  TRUE ~ "Other")

DimPlot(mapmycell_hier, group.by = "mmc_broad", label = FALSE)

# looking at ScType similarities by percentage
comparison_broad_h <- table(
  ScType = mapmycell_hier$sctype_label,
  MapMyCells = mapmycell_hier$mmc_broad
)

comparison_percent_h <- prop.table(comparison_broad_h, margin = 1) * 100
round(comparison_percent_h, 1)

#-----------------------------
#Saving both objects
qsave(mapmycell_cor, "mapmycell_cor.qs")
qsave(mapmycell_hier, "mapmycell_hier.qs")

