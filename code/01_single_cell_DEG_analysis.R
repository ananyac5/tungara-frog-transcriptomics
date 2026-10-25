
# ------------------------------------------------------------
# 0. THE LOADING DOCK (ALL LIBRARIES)
# ------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(scales) 
library(reshape2)
library(ggrepel)
library(openxlsx) 
library(scXCL)
library(qs)
library(pheatmap)
library(gridExtra) # Needed to save your tables as PDFs!

# ------------------------------------------------------------
# PART 1: DATA LOADING AND SETUP
# ------------------------------------------------------------
# Define the base directory where your data lives
base_dir <- "/stor/scratch/FRI-BigDataBio/FRI_summer_2026/frog_data"

# Load the 10X data matrices using the .h5 files
gn_data <- Read10X_h5(filename = file.path(base_dir, "TF_GN_outs/outs/filtered_feature_bc_matrix.h5"))
w_data  <- Read10X_h5(filename = file.path(base_dir, "TF_W_outs/outs/filtered_feature_bc_matrix.h5"))
wc_data <- Read10X_h5(filename = file.path(base_dir, "TF_WC_outs/outs/filtered_feature_bc_matrix.h5"))

# Create individual Seurat objects
gn_obj <- CreateSeuratObject(counts = gn_data, project = "Green_Noise")
w_obj  <- CreateSeuratObject(counts = w_data, project = "Condition_W")
wc_obj <- CreateSeuratObject(counts = wc_data, project = "Condition_WC")

# Add a metadata column for the stimulation condition
gn_obj$mating_call <- "Green_Noise"
w_obj$mating_call  <- "Condition_W"
wc_obj$mating_call <- "Condition_WC"

# Merge the three objects into one master object
frog_combined <- merge(x = gn_obj, y = c(w_obj, wc_obj), add.cell.ids = c("GN", "W", "WC"))

# Smash the separate treatment layers into one unified layer so the math doesn't crash later!
frog_combined <- JoinLayers(frog_combined)

# ------------------------------------------------------------
# PART 2: STANDARD SEURAT PROCESSING
# ------------------------------------------------------------
# Run the standard analysis workflow directly on the layered object
frog_combined <- NormalizeData(frog_combined)
frog_combined <- FindVariableFeatures(frog_combined)
frog_combined <- ScaleData(frog_combined)
frog_combined <- RunPCA(frog_combined)

# Generate the Elbow Plot to look at the top 50 Principal Components
my_elbow <- ElbowPlot(frog_combined, ndims = 50)
print(my_elbow)

# Build the neighborhood graph and cluster (resolution 0.2)
frog_combined <- FindNeighbors(frog_combined, dims = 1:20, reduction = "pca")
frog_combined <- FindClusters(frog_combined, resolution = 0.2) 

# Run UMAP and generate plots
frog_combined <- RunUMAP(frog_combined, dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated")
my_umap <- DimPlot(frog_combined, reduction = "umap.unintegrated", group.by = "seurat_clusters")
print(my_umap)

# ------------------------------------------------------------
# PART 3: AUTOMATED SCTYPE ANNOTATION
# ------------------------------------------------------------
# Download the scType functions directly from their GitHub
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# Load the official database and tell it we are looking at "Brain" tissue
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
tissue <- "Brain"
gs_list <- gene_sets_prepare(db_, tissue)

# Grab scaled data and run scoring math
scRNAseqData <- frog_combined[["RNA"]]$scale.data
es.max <- sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# Figure out which cell type scored highest
cL_results <- do.call("rbind", lapply(unique(frog_combined@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(frog_combined@meta.data[frog_combined@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(frog_combined@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores <- cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

# Add these newly discovered labels
frog_combined@meta.data$sctype_classification = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]
  frog_combined@meta.data$sctype_classification[frog_combined@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

# Generate the annotated UMAP plots
annotated_umap <- DimPlot(frog_combined, reduction = "umap.unintegrated", label = TRUE, repel = TRUE, group.by = 'sctype_classification') + ggtitle("sctype_classification")
print(annotated_umap)

split_umap <- DimPlot(frog_combined, reduction = "umap.unintegrated", split.by = "orig.ident")
print(split_umap)

# ------------------------------------------------------------
# PART 4: AUTOMATED XENOPUS scXCL MAPPING
# ------------------------------------------------------------
# Extract counts and force Túngara genes to lowercase
counts = LayerData(frog_combined, assay = "RNA", layer = "counts")
tung_genes = tolower(rownames(counts))
rownames(counts) = tung_genes

# Extract the Xenopus dictionary and match genes
xeno_genes = rownames(scXCL:::ref.expr)
xeno_base = sub("\\..*", "", xeno_genes)
shared_bases = intersect(tung_genes, xeno_base)
print(paste("Successfully mapped", length(shared_bases), "genes!"))

# Build Translation Dictionary and run prediction
xeno_full_names = xeno_genes[match(shared_bases, xeno_base)]
counts_matched = counts[shared_bases, ]
rownames(counts_matched) = xeno_full_names

scXCL_results = scXCL(scdata = counts_matched, numbers_plot = 1)
predictions = scXCL_results$scXCL

# Add predictions to Seurat object
if (length(predictions) > 0) {
  names(predictions) = colnames(frog_combined)
  frog_combined = AddMetaData(frog_combined, metadata = predictions, col.name = "Xenopus_Prediction")
  final_xenopus_umap <- DimPlot(frog_combined, reduction = "umap.unintegrated", group.by = "Xenopus_Prediction", label = TRUE, repel = TRUE) + ggtitle("Official Xenopus Validation")
  print(final_xenopus_umap)
} else {
  print("Error: Automated prediction failed.")
}

# ------------------------------------------------------------
# PART 5: STANDARD DIFFERENTIAL EXPRESSION (WC vs GN)
# ------------------------------------------------------------
frog_combined$cluster_treatment <- paste(frog_combined$seurat_clusters, frog_combined$orig.ident, sep = "_")
Idents(frog_combined) <- "cluster_treatment"

all_deg_results <- data.frame() 
cluster_numbers <- levels(frog_combined$seurat_clusters)

for (cluster in cluster_numbers) {
  print(paste("Calculating standard Seurat math for Cluster:", cluster))
  experimental_group <- paste(cluster, "Condition_WC", sep = "_")
  control_group <- paste(cluster, "Green_Noise", sep = "_")
  
  try({
    cluster_degs <- FindMarkers(frog_combined, ident.1 = experimental_group, ident.2 = control_group, min.pct = 0.1, logfc.threshold = 0.1)
    cluster_degs$Cluster <- cluster
    cluster_degs$Gene <- rownames(cluster_degs)
    all_deg_results <- rbind(all_deg_results, cluster_degs)
  }, silent = TRUE)
}

significant_degs <- subset(all_deg_results, p_val_adj < 0.05)
write.csv(significant_degs, file = "All_Clusters_WC_vs_GN_Standard.csv", row.names = FALSE)
print("Standard Pipeline Complete! File saved.")

# ============================================================
# PART 6: CELL TYPE INVENTORY (30PC 0.6 RESOLUTION)
# ============================================================

cell_roster <- frog@meta.data

# Group using the 30PC 0.6 resolution labels
inventory <- cell_roster %>%
  group_by(orig.ident, sctype_label_30PC_0.6) %>%
  summarise(Cell_Count = n(), .groups = "drop")

# 1. Raw Counts Plot
counts_plot <- ggplot(inventory, aes(x = orig.ident, y = Cell_Count, fill = sctype_label_30PC_0.6)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "Cell Type Inventory: Raw Counts (30PC 0.6)", 
    x = "Acoustic Treatment", 
    y = "Total Number of Cells", 
    fill = "Cell Type"
  )

pdf("Cell_Type_Inventory_Raw_Counts_30PC.pdf", width = 10, height = 7)
print(counts_plot)
dev.off()

# 2. Relative Proportions Plot
proportions_plot <- ggplot(inventory, aes(x = orig.ident, y = Cell_Count, fill = sctype_label_30PC_0.6)) +
  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous(labels = scales::percent) + 
  theme_minimal() +
  labs(
    title = "Cell Type Inventory: Relative Proportions (30PC 0.6)", 
    x = "Acoustic Treatment", 
    y = "Percent of Total Brain (%)", 
    fill = "Cell Type"
  )

pdf("Cell_Type_Inventory_Proportions_30PC.pdf", width = 10, height = 7)
print(proportions_plot)
dev.off()

print("Inventory graphs successfully refreshed using sctype_label_30PC_0.6!")

# ------------------------------------------------------------
# PART 7: WHINE-CHUCK vs WHINE (NEURONS ONLY)
# ------------------------------------------------------------
DefaultAssay(frog) <- "RNA"
frog$cell_treatment <- paste(frog$sctype_classification, frog$orig.ident, sep = "_")
Idents(frog) <- "cell_treatment"

all_deg_results <- data.frame()
cell_types <- unique(frog$sctype_classification)
cell_types <- cell_types[grepl("neuron", cell_types, ignore.case = TRUE)]

print(paste("Found", length(cell_types), "neuron types! Starting WC vs W math..."))

for (cell in cell_types) {
  experimental_group <- paste(cell, "WC", sep = "_")
  control_group <- paste(cell, "W", sep = "_")
  
  try({
    cluster_degs <- FindMarkers(frog, ident.1 = experimental_group, ident.2 = control_group, min.pct = 0.1, logfc.threshold = 0.1)
    cluster_degs$Cell_Type <- cell
    cluster_degs$Gene <- rownames(cluster_degs)
    all_deg_results <- rbind(all_deg_results, cluster_degs)
  }, silent = TRUE) 
}

significant_degs <- subset(all_deg_results, p_val_adj < 0.05)
write.csv(significant_degs, file = "/stor/scratch/FRI-BigDataBio/FRI_summer_2026/frog_data/Significant_DEGs_WC_vs_W_Neurons.csv", row.names = FALSE)

# Chuck Neuron Volcano Gallery
my_degs <- read.csv("/stor/scratch/FRI-BigDataBio/FRI_summer_2026/frog_data/Significant_DEGs_WC_vs_W_Neurons.csv")
neuron_types <- unique(my_degs$Cell_Type)

pdf("/stor/scratch/FRI-BigDataBio/FRI_summer_2026/frog_data/Chuck_Neuron_Volcanos.pdf", width=8, height=6)
for (target_cell in neuron_types) {
  volcano_data <- my_degs %>% filter(Cell_Type == target_cell)
  volcano_data$Expression <- "Not Significant"
  volcano_data$Expression[volcano_data$avg_log2FC > 0.5 & volcano_data$p_val_adj < 0.05] <- "Up in WC"
  volcano_data$Expression[volcano_data$avg_log2FC < -0.5 & volcano_data$p_val_adj < 0.05] <- "Down in WC"
  
  volcano_plot <- ggplot(volcano_data, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Expression)) +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_manual(values = c("Down in WC" = "blue", "Up in WC" = "red", "Not Significant" = "grey")) +
    theme_minimal() +
    labs(title = paste("Volcano Plot:", target_cell), subtitle = "Whine-Chuck vs. Whine", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
    geom_vline(xintercept = c(-0.5, 0.5), col = "black", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed")
  
  print(volcano_plot)
}
dev.off()


# ============================================================
# PART 8: THE ZOOM-IN ANALYSIS (GABAergic vs EXCITATORY)
# ============================================================

# ------------------------------------------------------------
# SUB-SECTION A: THE BRAKES (GABAergic NEURONS)
# ------------------------------------------------------------
# 1. Isolate strictly the GABAergic clusters
gaba_neurons <- subset(frog, idents = c(
  "GABAergic neuron_tac1 high",
  "GABAergic neuron_sst high",
  "GABAergic neuron_scgn high"
))
Idents(gaba_neurons) <- "orig.ident"

# 2. Run the DEG Math 
gaba_degs <- FindMarkers(gaba_neurons, ident.1 = "WC", ident.2 = "W", logfc.threshold = 0.25)
gaba_degs$gene <- rownames(gaba_degs)
write.csv(gaba_degs, "GABAergic_DEGs_WC_vs_W.csv", row.names = FALSE)

# 3. Generate the GABA Volcano Plot
gaba_degs$Significance <- "Not Significant"
gaba_degs$Significance[gaba_degs$p_val_adj < 0.05 & gaba_degs$avg_log2FC > 0.5] <- "Upregulated in WC"
gaba_degs$Significance[gaba_degs$p_val_adj < 0.05 & gaba_degs$avg_log2FC < -0.5] <- "Downregulated in WC"

gaba_volcano <- ggplot(gaba_degs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("blue", "gray", "red")) +
  theme_minimal() +
  labs(title = "GABAergic Neurons: Whine-Chuck vs Whine", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value")
ggsave("GABAergic_Volcano.pdf", plot = gaba_volcano, width = 7, height = 5)

# 4. Generate the GABA Clustered Heatmap (Top 25%)
gaba_sig <- gaba_degs %>% filter(p_val_adj < 0.05) %>% arrange(p_val_adj, desc(abs(avg_log2FC)))
n_top_gaba <- min(50, ceiling(nrow(gaba_sig) * 0.25))

if(n_top_gaba > 0) {
  top_genes_gaba <- gaba_sig %>% slice(1:n_top_gaba) %>% pull(gene)
  avg_data_gaba <- AverageExpression(gaba_neurons, features = top_genes_gaba, group.by = "orig.ident")$RNA
  scaled_gaba <- t(scale(t(avg_data_gaba)))
  scaled_gaba <- scaled_gaba[complete.cases(scaled_gaba), ]
  
  pdf("GABAergic_Clustered_Heatmap.pdf", width = 6, height = 8)
  pheatmap(
    mat = scaled_gaba,           
    cluster_rows = TRUE,         
    cluster_cols = TRUE,         
    main = "Top 25% GABAergic DEGs",
    color = colorRampPalette(c("blue", "black", "yellow"))(50), 
    angle_col = "45",
    border_color = NA 
  )
  dev.off()
}

# 5. Table 2: GABAergic Top 10 VIP Genes 
top_10_gaba <- gaba_degs %>%
  arrange(p_val_adj) %>%
  select(gene, avg_log2FC, p_val_adj, Significance) %>%
  head(10)

pdf("GABAergic_Top10_VIP_Genes_R.pdf", height = 4, width = 8)
grid.table(top_10_gaba)
dev.off()
print("GABAergic analysis complete!")


# ------------------------------------------------------------
# SUB-SECTION B: THE GAS (EXCITATORY NEURONS)
# ------------------------------------------------------------
# 1. Isolate strictly the Excitatory clusters
ex_neurons <- subset(frog, idents = c(
  "Excitatory neuron_zbtb18 high",
  "Excitatory neuron_tac1 high",
  "Excitatory neuron_mex3a high"
))
Idents(ex_neurons) <- "orig.ident"

# 2. Run the DEG Math 
ex_degs <- FindMarkers(ex_neurons, ident.1 = "WC", ident.2 = "W", logfc.threshold = 0.25)
ex_degs$gene <- rownames(ex_degs)
write.csv(ex_degs, "Excitatory_DEGs_WC_vs_W.csv", row.names = FALSE)

# 3. Generate the Excitatory Volcano Plot
ex_degs$Significance <- "Not Significant"
ex_degs$Significance[ex_degs$p_val_adj < 0.05 & ex_degs$avg_log2FC > 0.5] <- "Upregulated in WC"
ex_degs$Significance[ex_degs$p_val_adj < 0.05 & ex_degs$avg_log2FC < -0.5] <- "Downregulated in WC"

ex_volcano <- ggplot(ex_degs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("blue", "gray", "red")) +
  theme_minimal() +
  labs(title = "Excitatory Neurons: Whine-Chuck vs Whine", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value")
ggsave("Excitatory_Volcano.pdf", plot = ex_volcano, width = 7, height = 5)

# 4. Generate the Excitatory Clustered Heatmap (Top 25%) - FIXED!
ex_sig <- ex_degs %>% filter(p_val_adj < 0.05) %>% arrange(p_val_adj, desc(abs(avg_log2FC)))
n_top_ex <- min(50, ceiling(nrow(ex_sig) * 0.25))

if(n_top_ex > 0) {
  top_genes_ex <- ex_sig %>% slice(1:n_top_ex) %>% pull(gene)
  avg_data_ex <- AverageExpression(ex_neurons, features = top_genes_ex, group.by = "orig.ident")$RNA
  
  scaled_ex <- t(scale(t(avg_data_ex)))
  scaled_ex <- scaled_ex[complete.cases(scaled_ex), ]
  
  pdf("Excitatory_Clustered_Heatmap.pdf", width = 6, height = 8)
  pheatmap(
    mat = scaled_ex,             
    cluster_rows = TRUE,         
    cluster_cols = TRUE,         
    main = "Top 25% Excitatory DEGs",
    color = colorRampPalette(c("blue", "black", "yellow"))(50), 
    angle_col = "45",
    border_color = NA 
  )
  dev.off()
}

# 5. Table 1: Excitatory Top 10 VIP Genes - FIXED TO PDF!
top_10_ex <- ex_degs %>%
  arrange(p_val_adj) %>%
  select(gene, avg_log2FC, p_val_adj, Significance) %>%
  head(10)

pdf("Excitatory_Top10_VIP_Genes_R.pdf", height = 4, width = 8)
grid.table(top_10_ex)
dev.off()

print("Clustered Excitatory heatmap and table complete! Both pipelines are now fully upgraded.")


# ============================================================
# VENN DIAGRAM: GLUTAMATERGIC VS GABAERGIC OVERLAP
# ============================================================

# 1. Load the required packages
library(ggVennDiagram)
library(ggplot2)
library(dplyr)

# 2. Extract ONLY the statistically significant genes (p_adj < 0.05)
gaba_sig_genes <- gaba_degs %>% 
  filter(p_val_adj < 0.05) %>% 
  pull(gene)

ex_sig_genes <- ex_degs %>% 
  filter(p_val_adj < 0.05) %>% 
  pull(gene)

# 3. Group them into a list 
# (Doubled the spaces on the RIGHT to push them further left!)
venn_list <- list(
  "GABAergic                    " = gaba_sig_genes,
  "Glutamatergic                    " = ex_sig_genes
)

# 4. Draw the Venn Diagram with polished formatting
venn_plot <- ggVennDiagram(
  venn_list, 
  color = "black", 
  lwd = 0.8, 
  lty = 1,
  set_name_size = 4.5 # Slightly shrinks the labels so they look cleaner
) + 
  scale_fill_gradient(low = "white", high = "dodgerblue") +
  # Increased the expansion from 0.4 to 0.5 to make room for the extra push
  scale_x_continuous(expand = expansion(mult = 0.5)) + 
  theme_void() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) + 
  labs(
    title = "DEG Overlap: Glutamatergic vs GABAergic Neurons",
    subtitle = "Significant Genes (Adjusted P-Value < 0.05) in Whine-Chuck vs Whine"
  )

# 5. Save it directly to your folder
ggsave("Venn_Diagram_Gluta_vs_GABA.pdf", plot = venn_plot, width = 8, height = 5)

print("Final Venn diagram completely saved! Check your folder.")

# ============================================================
# PART 9: CANDIDATE GENE MODULE MAPPING & GO ANALYSIS
# ============================================================

# 1. Load the required packages
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

print("Starting Candidate Gene Analysis...")

# 2. Define your exact candidate list (The Social Behavior MVPs!)
candidate_list <- c("AVP", "LOC140113373", "OXTR", "AVPR1A", "AVPR1B", "DRD1", "DRD2", "DRD4", "DRD5")

# 3. Check which genes actually survived the filtering in your blast object
# (Ensuring R only looks for genes that actually exist in your dataset)
valid_candidates <- candidate_list[candidate_list %in% rownames(tunga_frog_merged_blast)]
print(paste("Found", length(valid_candidates), "valid candidates in the dataset!"))

# 4. Find exactly which modules/cell types these genes are active in!
# We use the 'significant_degs' master spreadsheet we built earlier
candidate_modules <- significant_degs %>%
  filter(Gene %in% valid_candidates) %>%
  select(Gene, Cell_Type, avg_log2FC, p_val_adj) %>%
  arrange(Gene, p_val_adj)

# Print the results so you can see exactly where they are working in the console
print("Here is where your candidate genes are active:")
print(candidate_modules)

# Save this specific mini-table as a CSV for your PI
write.csv(candidate_modules, "Candidate_Genes_Modules.csv", row.names = FALSE)
print("Candidate_Genes_Modules.csv saved successfully!")

# ------------------------------------------------------------
# GO PATHWAY ENRICHMENT (CANDIDATE GENES)
# ------------------------------------------------------------
print("Running GO Pathway Analysis on candidate genes...")

# 5. Run the GO Enrichment Math on your valid candidates
go_results <- enrichGO(
  gene          = valid_candidates,  # Feeding your specific valid list here
  OrgDb         = org.Hs.eg.db,      # Using the comprehensive human database
  keyType       = "SYMBOL",          
  ont           = "BP",              # Searching for "Biological Processes"
  pAdjustMethod = "BH",              
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

# 6. Check if the math found significant pathways before drawing the plot
if (!is.null(go_results) && nrow(go_results) > 0) {
  
  # Draw the Dot Plot
  go_plot <- dotplot(go_results, showCategory = 10) +
    ggtitle("GO Pathway Analysis: Social Candidate Genes",
            subtitle = "Top Biological Processes") +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  
  # Save the final pathway plot to your folder
  ggsave("GO_Pathway_Candidate_Genes.pdf", plot = go_plot, width = 8, height = 6)
  print("GO Analysis completely finished! Check your folder for GO_Pathway_Candidate_Genes.pdf")
  
} else {
  print("Note: The candidate list was too small to hit the statistical threshold for broad GO pathways. This is completely normal for hyper-specific targeted gene lists!")
}

# ============================================================
# PART 12: INDIVIDUAL GO PATHWAY ENRICHMENT (GABA vs EXCITATORY)
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# ------------------------------------------------------------
# SUB-SECTION A: GABAergic INDIVIDUAL GO ENRICHMENT
# ------------------------------------------------------------
print("Starting Individual GABAergic GO Pathway Analysis...")

# 1. Grab the significant GABAergic genes vector (from your Venn setup)
print(paste("Running GO Enrichment on", length(gaba_sig_genes), "GABAergic significant genes..."))

# 2. Run the GO Enrichment Math
go_gaba <- enrichGO(
  gene          = gaba_sig_genes,
  OrgDb         = org.Hs.eg.db,      
  keyType       = "SYMBOL",          
  ont           = "BP",              # Biological Processes
  pAdjustMethod = "BH",              
  pvalueCutoff  = 0.15,
  qvalueCutoff  = 0.15 
)

# 3. Check for results, plot, and save text tables
if (!is.null(go_gaba) && nrow(go_gaba) > 0) {
  
  # Generate Dot Plot
  gaba_go_plot <- dotplot(go_gaba, showCategory = 15) +
    ggtitle("GO Enrichment: GABAergic Neuronal Responses",
            subtitle = "Biological Processes Altered in Whine-Chuck vs Whine") +
    theme_minimal() +
    theme(plot.title = element_text(size = 13, face = "bold"))
  
  ggsave("GABAergic_Individual_GO_Plot.pdf", plot = gaba_go_plot, width = 10, height = 7)
  write.csv(as.data.frame(go_gaba), "GABAergic_Individual_GO_Table.csv", row.names = FALSE)
  print("GABAergic individual GO pathway analysis saved successfully!")
  
} else {
  print("Note: No hyper-specific pathways passed the 0.05 cutoff for the individual GABA list.")
}


# ------------------------------------------------------------
# SUB-SECTION B: EXCITATORY INDIVIDUAL GO ENRICHMENT
# ------------------------------------------------------------
print("Starting Individual Excitatory GO Pathway Analysis...")

# 1. Grab the significant Excitatory genes vector (from your Venn setup)
print(paste("Running GO Enrichment on", length(ex_sig_genes), "Excitatory significant genes..."))

# 2. Run the GO Enrichment Math
go_ex <- enrichGO(
  gene          = ex_sig_genes,
  OrgDb         = org.Hs.eg.db,      
  keyType       = "SYMBOL",          
  ont           = "BP",              # Biological Processes
  pAdjustMethod = "BH",              
  pvalueCutoff  = 0.15,
  qvalueCutoff  = 0.15 
)

# 3. Check for results, plot, and save text tables
if (!is.null(go_ex) && nrow(go_ex) > 0) {
  
  # Generate Dot Plot
  ex_go_plot <- dotplot(go_ex, showCategory = 15) +
    ggtitle("GO Enrichment: Excitatory Neuronal Responses",
            subtitle = "Biological Processes Altered in Whine-Chuck vs Whine") +
    theme_minimal() +
    theme(plot.title = element_text(size = 13, face = "bold"))
  
  ggsave("Excitatory_Individual_GO_Plot.pdf", plot = ex_go_plot, width = 10, height = 7)
  write.csv(as.data.frame(go_ex), "Excitatory_Individual_GO_Table.csv", row.names = FALSE)
  print("Excitatory individual GO pathway analysis saved successfully!")
  
} else {
  print("Note: No hyper-specific pathways passed the 0.05 cutoff for the individual Excitatory list.")
}

# ============================================================
# PART 14: CALCULATE DEGS & DRAW DENDROGRAM HEATMAP (FIXED)
# ============================================================
library(dplyr)
library(pheatmap)
library(Seurat)

DefaultAssay(frog) <- "RNA"

# 1. Set the exact identity to Jace's updated 45 PC column!
Idents(frog) <- "sctype_label_45PC_0.6_updated"

# 2. Calculate the markers from scratch
print("Calculating all 45 PC markers... this might take a minute!")
fresh_markers <- FindAllMarkers(frog, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save the real CSV so you have it safely backed up!
write.csv(fresh_markers, "Jace_45PC_Cluster_Markers.csv", row.names = FALSE)

# 3. Extract the Top 10 most significant genes per cluster
# Note: Seurat usually outputs lowercase 'cluster', but we will use the exact column name
top10_markers <- fresh_markers %>%
  group_by(cluster) %>%  
  top_n(n = 10, wt = avg_log2FC)

# Pull out just the unique gene names for the plot
top_genes <- unique(top10_markers$gene)

# 4. Calculate average expression for the heatmap
print("Calculating average expression per cluster...")
avg_expression <- AverageExpression(frog, features = top_genes, return.seurat = FALSE)
expr_matrix <- avg_expression$RNA

# 5. Generate the heatmap with dendrograms
print("Drawing the heatmap...")
pdf("Global_Dendrogram_Heatmap_Jace_Pretty.pdf", width = 12, height = 10)

pheatmap(
  mat = expr_matrix,
  scale = "row",              
  cluster_rows = TRUE,        
  cluster_cols = TRUE,        
  show_rownames = FALSE,      
  color = colorRampPalette(c("skyblue", "white", "lightcoral"))(100), # <--- ADDED THE COMMA RIGHT HERE!
  main = "Top DEGs Across 45 PC Clusters"
)

dev.off()

print("Success! Check your folder for the new pretty heatmap!")