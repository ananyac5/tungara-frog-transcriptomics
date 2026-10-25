# Brain gene-expression graphs
# ----------------------------
library(Seurat)
library(patchwork)
library(ggplot2)

DefaultAssay(combined_int_test) <- "RNA"

genes <- c(
  "CHGA", "CHGB",
  "FSHB", "LOC140066080", "PRL", "TSHB",
  "SLC32A1", "GAD2",
  "MBP", "OLIG2", "PDGFRA", "PMP22"
)

genes_present <- genes[genes %in% rownames(combined_int_test)]
genes_missing <- setdiff(genes, genes_present)

genes_present
genes_missing

plots_all <- FeaturePlot(
  combined_int_test,
  features = genes_present,
  reduction = "umap_30PC_0.6",
  order = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q95",
  combine = FALSE
)

plots_all <- lapply(plots_all, function(p) {
  p +
    NoAxes() +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = 8)
    )
})

p_all <- wrap_plots(plots_all, ncol = 4)
p_all
#=========================================
#Graphs of genes grouped by cell type
#=========================================

# Endocrine genes
endocrine_genes <- intersect(
  c("CHGA", "CHGB", "FSHB", "LOC140066080", "PRL", "TSHB"),
  rownames(combined_int_test)
)

p_endocrine <- FeaturePlot(
  combined_int_test,
  features = endocrine_genes,
  reduction = "umap_30PC_0.6",
  order = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q95",
  ncol = 2
)

p_endocrine & scale_color_gradient(low = "grey90", high = "blue")

# GABAergic genes
neuron_genes <- intersect(
  c("SLC32A1", "GAD2"),
  rownames(combined_int_test)
)

p_neurons <- FeaturePlot(
  combined_int_test,
  features = neuron_genes,
  reduction = "umap_30PC_0.6",
  order = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q95",
  ncol = 2
)

p_neurons & scale_color_gradient(low = "grey90", high = "blue")

# Glial genes
glia_genes <- intersect(
  c("MBP", "OLIG2", "PDGFRA", "PMP22"),
  rownames(combined_int_test)
)

p_glia <- FeaturePlot(
  combined_int_test,
  features = glia_genes,
  reduction = "umap_30PC_0.6",
  order = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q95",
  ncol = 2
)

p_glia & scale_color_gradient(low = "grey90", high = "blue")

#-----------------------------------------------
#comparison of blast and ortholog
#-----------------------------------------------
library(dplyr)
library(readxl)

setwd("/stor/scratch/FRI-BigDataBio/FRI_summer_2026/frog_data/ortho_test")

blast_table <- read_excel(
  "translation_table_xenopus_to_tungara.xlsx"
)

ortholog_table <- read.csv(
  "xenopus_to_tungara_ortholog_table.csv",
  stringsAsFactors = FALSE
)

# BLAST: X. tropicalis cleaned gene symbol
blast_clean <- blast_table %>%
  transmute(
    xenopus_key = tolower(trimws(xenopus_clean)),
    blast_xenopus_symbol = trimws(xenopus_symbol),
    blast_tungara_protein = trimws(tungara_protein),
    blast_tungara_symbol = toupper(trimws(tungara_symbol))
  ) %>%
  filter(
    !is.na(xenopus_key),
    xenopus_key != ""
  ) %>%
  distinct()

# OrthoFinder: X. laevis gene name is stored in xenopus_gene
ortholog_clean <- ortholog_table %>%
  transmute(
    Orthogroup,
    xenopus_key = tolower(trimws(xenopus_gene)),
    ortholog_xenopus_gene = trimws(xenopus_gene),
    ortholog_tungara_protein = trimws(tungara_protein),
    ortholog_tungara_gene = trimws(tungara_gene),
    ortholog_tungara_symbol = toupper(trimws(tungara_gene_name))
  ) %>%
  mutate(
    # Remove X. laevis homeolog suffixes
    xenopus_key = sub("\\.(l|s)$", "", xenopus_key)
  ) %>%
  filter(
    !is.na(xenopus_key),
    xenopus_key != ""
  ) %>%
  distinct()

#unqiue genes
blast_unique_genes <- blast_table %>%
  mutate(
    xenopus_key = tolower(trimws(xenopus_clean))
  ) %>%
  filter(
    !is.na(xenopus_key),
    xenopus_key != ""
  ) %>%
  distinct(xenopus_key)

ortholog_unique_genes <- ortholog_table %>%
  mutate(
    xenopus_key = tolower(trimws(xenopus_gene)),
    xenopus_key = sub("\\.(l|s)$", "", xenopus_key)
  ) %>%
  filter(
    !is.na(xenopus_key),
    xenopus_key != ""
  ) %>%
  distinct(xenopus_key)

nrow(blast_unique_genes)
nrow(ortholog_unique_genes)

