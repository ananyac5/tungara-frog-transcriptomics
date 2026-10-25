setwd("/stor/scratch/FRI-BigDataBio/FRI_summer_2026/frog_data/Hedi_R")

library(Seurat)
library(dplyr)
library(readxl)
library(tidyr)
library(HGNChelper)
library(ggplot2)
library(future)
library(devtools)
library(hdWGCNA)
library(WGCNA)
library(writexl)
library(pheatmap)

options(future.globals.maxSize = 64000 * 1024^2)
Sys.setenv(R_MAX_VSIZE = "64Gb")





tunga_frog_merged_blast = readRDS("tunga_frog_annotated_blast_wgcna.rds")





PlotDendrogram(tunga_frog_merged_blast, main = 'Tungara Frog Brain Gene Network')

ModuleFeaturePlot(
  tunga_frog_merged_blast,
  features = 'hMEs', # Module Eigengenes
  order = TRUE,       # cells with the highest expression to the front
  reduction = 'umap_by_pca'
)

head(GetModules(tunga_frog_merged_blast))


#graphics.off()

HubGeneNetworkPlot(
  tunga_frog_merged_blast,
  n_hubs = 10,       # 25 most connected genes
#  n_other = 2,       # 5 loosely connected genes for context
  edge_prop = 0.75,  # thickness of the lines
  mods = "green" # Specific module color
)


# VlnPlot(
#   tunga_frog_merged_blast,
#   features = "blue", 
#   group.by = "orig.ident", 
#   pt.size = 0
# ) + ggtitle("Blue: Gaba Only")

VlnPlot(
  tunga_frog_merged_blast,
  features = c("brown", "yellow", "blue", "green"), 
  group.by = "orig.ident", #
  pt.size = 0,
  ncol = 2
)

modules = GetModules(tunga_frog_merged_blast)

# Count how many genes are in each color module
round(table(modules$module)/sum(table(modules$module)) * 100, 2)

# Calculate the percentage of genes in the grey module
grey_genes = sum(modules$module == "grey")
total = nrow(modules)
grey_percentage = (grey_genes / total) * 100
print(round(grey_percentage, 2))



hub_green_genes = modules %>%
  dplyr::filter(module == "green") %>%
  dplyr::arrange(desc(kME_green)) %>%
  head(10)


# LOAD BLAST MATCHES
orthos = read.csv("xenopus_to_tungara_RBH.csv")
orthos$xenopus_protein = trimws(gsub("\\..*", "", orthos$xenopus_protein))
orthos$tungara_protein = trimws(gsub("\\..*", "", orthos$tungara_protein))
print(paste(nrow(orthos), "protein matches"))

# LOAD THE GOLDEN TROPICALIS DICT FOUND
xeno_id_to_symbol = read.csv("xenopus_id_to_symbol_tropicalis.csv") %>% 
  dplyr::rename(xenopus_protein = protein_id, xenopus_symbol = gene_symbol) %>%
  dplyr::filter(grepl("^[XN]P_", xenopus_protein)) 

xeno_id_to_symbol$xenopus_protein = trimws(gsub("\\..*", "", xeno_id_to_symbol$xenopus_protein))

# LOAD THE REFSEQ DICTIONARY
tunga_id_to_symbol = read.csv("tungara_id_to_symbol_FIXED.csv") %>% 
  dplyr::rename(tungara_protein = protein_id, tungara_symbol = gene_symbol)  

tunga_id_to_symbol$tungara_protein = trimws(gsub("\\..*", "", tunga_id_to_symbol$tungara_protein))


translation_table = orthos %>%
  dplyr::inner_join(xeno_id_to_symbol, by = "xenopus_protein") %>%
  dplyr::inner_join(tunga_id_to_symbol, by = "tungara_protein") %>%
  dplyr::mutate(
    xenopus_clean = toupper(xenopus_symbol),
    tungara_symbol = toupper(tungara_symbol)
  ) %>%
  dplyr::distinct(xenopus_clean, tungara_symbol, .keep_all = TRUE)



green_gene_names = modules %>% 
  dplyr::filter(module == "green") %>%
  dplyr::pull(gene_name) # turns column into a text list


green_module_list = data.frame(tungara_symbol = green_gene_names)

# Join with table to find symbols for LOCs
translated_list = green_module_list %>%
  dplyr::left_join(tunga_id_to_symbol, by = "tungara_symbol") %>%
  dplyr::mutate(final_id = ifelse(is.na(tungara_symbol), tungara_symbol, 
                                  ifelse(grepl("^LOC", tungara_symbol) & !is.na(tungara_symbol), 
                                         tungara_symbol, tungara_symbol))) 

# We filter out the 'LOC' entries because GO databases don't recognize them
shinygo_input = green_gene_names[!grepl("^LOC", green_gene_names)]

print(paste(unique(shinygo_input), collapse = "\n"))
length(unique(shinygo_input))


green_genes_tt = data.frame(tungara_symbol = green_gene_names)

# Join against translation table to find the Xenopus Protein ID
full_protein_list = green_genes_tt %>%
  dplyr::left_join(translation_table, by = "tungara_symbol")

# Xenopus Protein IDs
xenopus_protein_ids = unique(full_protein_list$xenopus_protein)

# Filter out any NAs, genes that had no BLAST hit
xenopus_protein_ids_clean = xenopus_protein_ids[!is.na(xenopus_protein_ids)]

print(paste(xenopus_protein_ids_clean, collapse = "\n"))
length(xenopus_protein_ids_clean)


# From csv table
# XP_004910415,XP_012821611, XP_004912120

jnk_genes = c("MAPK4", "MAPK8", "TNFRSF19")

DotPlot(
  tunga_frog_merged_blast, 
  features = jnk_genes, 
  group.by = "orig.ident"
) + 
  ggtitle("JNK expression by treatment")


Idents(tunga_frog_merged_blast) = "orig.ident"

# Targeted DE test for JNK genes
jnk_stats = FindMarkers(
  tunga_frog_merged_blast, 
  ident.1 = "Whine-chuck", 
  ident.2 = "Whine", 
  features = c("MAPK4", "MAPK8", "TNFRSF19"), 
  logfc.threshold = 0 # output results even if fold change is low
)

print(jnk_stats)


DotPlot(
  tunga_frog_merged_blast, 
  features = c("MAPK8", "MAPK4", "TNFRSF19"), 
  group.by = "refined_annotation", 
  cols = c("lightgrey", "blue", "red")
)







tunga_frog_merged_blast$Celltype_treatment = paste(
  tunga_frog_merged_blast$refined_annotation,
  tunga_frog_merged_blast$orig.ident,       
  sep = "_"
)

cells_to_plot = c("Astrocyte_slc1a3 high", 
                  "Myelinating oligodendrocyte progenitor cell", 
                  "Capillary endothelial cell", 
                  "Antigen-presenting cell",
                  "GABAergic neuron_sst high",
                  "Excitatory neuron_tac1 high")

target_subset = subset(tunga_frog_merged_blast, subset = refined_annotation %in% cells_to_plot)

DotPlot(
  target_subset, 
  features = c("MAPK8"), 
  group.by = "Celltype_treatment"
) 










# DE test on 142 genes
green_genes_stats = FindMarkers(
  tunga_frog_merged_blast, 
  ident.1 = "Whine-chuck", 
  ident.2 = "Whine", 
  features = green_gene_names,
  logfc.threshold = 0 
)

# ONLY genes that are statistically significant
significant_green_genes = green_genes_stats %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::arrange(p_val_adj)

head(significant_green_genes, n = 15)


VlnPlot(
  tunga_frog_merged_blast, 
  features = c("NRTN", "CAPRIN2"), 
  group.by = "orig.ident", 
  pt.size = 0
)

# Vasopressin and Oxytocin related symbols
translation_table %>% 
  dplyr::filter(grepl("AVP|OXT|TOCIN|AVT", xenopus_symbol, ignore.case = TRUE) | 
                  grepl("AVP|OXT|TOCIN|AVT", tungara_symbol, ignore.case = TRUE)) %>% 
  dplyr::select(xenopus_symbol, tungara_symbol) %>% 
  distinct()

# Dopamine receptors
translation_table %>% 
  dplyr::filter(grepl("DRD|DOPAMINE", xenopus_symbol, ignore.case = TRUE) | 
                  grepl("DRD|DOPAMINE", tungara_symbol, ignore.case = TRUE)) %>% 
  dplyr::select(xenopus_symbol, tungara_symbol) %>% 
  distinct()


candidate_list = c("AVP", "LOC140113373", "OXTR", "AVPR1A", "AVPR1B", "DRD1", "DRD2", "DRD4", "DRD5")
candidate_list[candidate_list %in% rownames(tunga_frog_merged_blast)]

# Check which modules these genes belong
modules %>% 
  dplyr::filter(gene_name %in% candidate_list) %>%
  dplyr::select(gene_name, module, matches("kME_"))



DotPlot(
  tunga_frog_merged_blast, 
  features = candidate_list, 
  group.by = "Celltype_treatment"
) + RotatedAxis() +theme(axis.text = element_text(size = 8)) 



DotPlot(
  tunga_frog_merged_blast, 
  features = candidate_list, 
  group.by = "refined_annotation"
) + RotatedAxis() + theme(axis.text = element_text(size = 8)) 



DotPlot(
  target_subset, 
  features = c("NRTN", "CAPRIN2"), 
  group.by = "Celltype_treatment"
) + RotatedAxis() 


unique(tunga_frog_merged_blast$broad_annotation)

cells_of_interest = c(
  "GABAergic Neuron", 
  "Gonadotroph cell_cga high",
  "Excitatory Neuron"
)

key_features = c("DRD1", "DRD2", "OXTR", "MAPK8", "NRTN", "CAPRIN2")

target_subset = subset(tunga_frog_merged_blast, subset = broad_annotation %in% cells_of_interest)

DotPlot(
  target_subset, 
  features = key_features, 
  group.by = "Celltype_treatment"
) + RotatedAxis() 








# Do cells with higher dopamine sensitivity, more DRD2 baseline, have a stronger growth response to the mating call?

# Subset to only Whine-chuck
wc_cells = subset(tunga_frog_merged_blast, subset = orig.ident == "Whine-chuck")

# DRD2 expression of cells
drd2_expr = GetAssayData(wc_cells, assay = "RNA", slot = "data")["DRD2", ]

# Categorize cells, if expressing (> 0), if not expressing (== 0)
wc_cells$DRD2_status = ifelse(drd2_expr > 0, "DRD2_positive", "DRD2_negative")

Idents(wc_cells) = "DRD2_status"

# Statistical test for growth and JNK
sensitivity_stats = FindMarkers(
  wc_cells,
  ident.1 = "DRD2_positive",
  ident.2 = "DRD2_negative",
  features = c("NRTN", "MAPK8", "CAPRIN2"),
  logfc.threshold = 0
)

print(sensitivity_stats)

VlnPlot(
  wc_cells, 
  features = c("NRTN", "MAPK8", "CAPRIN2"), 
  group.by = "DRD2_status", 
  pt.size = 0,
  cols = c("red", "blue")
)


# Do cells with higher dopamine sensitivity, more DRD2 baseline, have a stronger growth response to the mating call?
#-- Answer => No
# favors neurons that don't have DRD2







# Checking if DRD2 positive cells are producing AVP or OXT

FeaturePlot(
  wc_cells,
  features = c("DRD2", "CAPRIN2", "AVP", "OXTR"),
  reduction = "umap_by_pca",
  ncol = 2,
  order = TRUE,
  pt.size = 1
)

DotPlot(
  wc_cells,
  features = c("DRD2", "CAPRIN2", "AVP", "OXTR", "LOC140113373"), 
  group.by = "refined_annotation"
) + RotatedAxis() 

# Checking if DRD2 positive cells are producing AVP or OXT
#-- Answer => Yes
# DRD2 positive Gonadotroph cells are producing Vasopressin and are highly expressed in OXT



# How AVP OXTR and CAPRIN2 react in Gonadotrophs cells when the frog hears the Whine-chuck?

gonadotroph_cells = subset(tunga_frog_merged_blast, subset = refined_annotation == "Gonadotroph cell_cyp21a2.1 high")

Idents(gonadotroph_cells) = "orig.ident"

# Whine-chuck against Whine
hormone_response_stats = FindMarkers(
  gonadotroph_cells,
  ident.1 = "Whine-chuck",
  ident.2 = "Whine",
  features = c("AVP", "OXTR", "DRD2", "CAPRIN2"),
  logfc.threshold = 0 
)

print(hormone_response_stats)

# How AVP OXTR and CAPRIN2 react in Gonadotrophs cells when the frog hears the Whine-chuck?
#-- Answer => When the frog hears the Whine-chuck,in Gonadotroph cells, AVP and CAPRIN2 go from zero expression to somewhat expressed
# OXTR is stable and unaffected.




# Dopamine isn't the chemical signal telling these neurons to grow and remodel when they hear chuck

# We need to figure out the neural remodeling track.

# Search for glutamate receptors 
translation_table %>% 
  dplyr::filter(grepl("GRIN|GRIA|GLUTAMATE RECEPTOR", xenopus_symbol, ignore.case = TRUE) | 
                  grepl("GRIN|GRIA|GLUTAMATE RECEPTOR", tungara_symbol, ignore.case = TRUE)) %>% 
  dplyr::select(xenopus_symbol, tungara_symbol) %>% 
  distinct()

glutamate_candidates = c(
  "GRIA1", "GRIA2", "GRIA3", "GRIA4", 
  "GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D", "GRIN3A", "GRIN3B"
)

# Keep only the ones present
glutamate_candidates = glutamate_candidates[glutamate_candidates %in% rownames(tunga_frog_merged_blast)]
print(glutamate_candidates)

# Neurons
neurons_clusters = c(
  "GABAergic neuron_scgn high", 
  "GABAergic neuron_sst high",     
  "Excitatory neuron_tac1 high",
  "Excitatory neuron_zbtb18 high"
)

neural_subset = subset(tunga_frog_merged_blast, subset = refined_annotation %in% neurons_clusters)

DotPlot(
  neural_subset, 
  features = c(glutamate_candidates, "MAPK8", "NRTN"), 
  group.by = "Celltype_treatment"
) + RotatedAxis()

# These neurons have a big baseline of glutamate receptors but they only push the JNK pathway MAPK8 and NRTN when the Whine-chuck mating call




# Co-expression matrix to check the remaining genes in the green module

green_expr = as.matrix(GetAssayData(wc_cells, assay = "RNA", slot = "data")[green_gene_names, ])

# Pearson correlation matrix
gene_cor_matrix_green = cor(t(green_expr), method = "pearson")


pheatmap(
  gene_cor_matrix_green,
  show_rownames = FALSE, 
  show_colnames = FALSE, 
  main = "Green module networks for Whine-Chuck call",
  color = colorRampPalette(c("blue", "grey", "red"))(100),
  clustering_method = "ward.D2"
)
# red square


# Extract the gene names driving it
# Top 10 genes correlated with CAPRIN2
caprin2_correlations = sort(gene_cor_matrix_green["CAPRIN2", ], decreasing = TRUE)
print(head(caprin2_correlations, 11)) # 11 because CAPRIN2


# Top 10 genes correlated with MAPK8
mapk8_correlations = sort(gene_cor_matrix_green["MAPK8", ], decreasing = TRUE)
print(head(mapk8_correlations, 11))
# NRTN as the most correlated gene for mapk8
# NUMB and WIF1: neural rewiring genes

# gprofiler on entire 142 gene green module 
# foud the JNK cascade but results were not great because half the module was on endocrine work.
# rerun gprofiler again

# gprofiler on neural remodeling

# top 30 gene symbols correlated with MAPK8
top_neural_genes = names(head(mapk8_correlations, 31))

# Filter with translation table to get protein ids
neural_translation = translation_table %>%
  dplyr::filter(tungara_symbol %in% top_neural_genes)

# Extract the clean Xenopus Protein IDs
neural_xenopus_protein_ids = unique(neural_translation$xenopus_protein)

# Filter out NAs
neural_xenopus_protein_ids = neural_xenopus_protein_ids[!is.na(neural_xenopus_protein_ids)]

print(length(neural_xenopus_protein_ids))
print(neural_xenopus_protein_ids, collapse = "\n")


# Gene names to map on human instead

print(paste(top_neural_genes, collapse = "\n"))



# gprofiler results

gprofiler_results = read.csv("gProfiler_hsapiens_2026-06-23_21-29-12__intersections.csv")

hits = gprofiler_results %>%
  # only biological processes: what the cell is trying to achieve
  dplyr::filter(source == "GO:BP") %>%
  # only statistically significant results
  dplyr::filter(adjusted_p_value < 0.05) %>%
  # remove broad
  dplyr::filter(term_size < 2000) %>%
  # sort by the highest number of genes working together
  dplyr::arrange(desc(intersection_size), adjusted_p_value)

print(head(hits[, c("term_name", "adjusted_p_value", "term_size", "intersection_size")], 5))


# genes cell projection organization
projection_genes = c("NRTN", "RPGRIP1L", "NUMB", "CFAP74", "FGD6", "ARHGEF4", "DDR2", "AK9")

projection_genes = projection_genes[projection_genes %in% rownames(tunga_frog_merged_blast)]
print(projection_genes)

DotPlot(
  neural_subset, 
  features = projection_genes, 
  group.by = "Celltype_treatment"
) + RotatedAxis() 



# important baseline expression across treatment then more expressed during the Whine-chuck for NUMB and RPGRIP1L


# TF: Transcription factor: protein that grab onto DNA and turns genes on or off
math_table = gprofiler_results %>%
  dplyr::filter(source == "TF") %>%
  dplyr::mutate(
    Background_rate = term_size / effective_domain_size,
    Expected_hits = Background_rate * query_size,
    Enrichment_ratio = intersection_size / Expected_hits
  ) %>%
  # highest to lowest
  dplyr::arrange(desc(Enrichment_ratio)) %>%
  dplyr::select(term_name, term_size, Expected_hits, intersection_size, Enrichment_ratio)

print(math_table)



tf_candidates = c("SP4", "E2F1", "E2F2", "E2F3", "SMAD1", "SMAD2", "SMAD3", "SMAD4")

tf_candidates = tf_candidates[tf_candidates %in% rownames(tunga_frog_merged_blast)]
print(tf_candidates)

DotPlot(
  neural_subset, 
  features = tf_candidates, 
  group.by = "Celltype_treatment"
) + RotatedAxis()



# E2Fs and SMAD1, SMAD4: low expression

# SMAD2, SMAD3, and SP4: imortant baseline expression, even more on whine chuck








all_neurons = unique(tunga_frog_merged_blast$refined_annotation[
  grepl("neuron", tunga_frog_merged_blast$refined_annotation, ignore.case = TRUE)
])

neural_obj = subset(tunga_frog_merged_blast, subset = refined_annotation %in% all_neurons)



dopamine_synthesis = c("TH", "DDC")
dopamine_receptors = c("DRD1", "DRD2", "DRD3", "DRD4", "DRD5")
growth_factors = c("NRTN", "CAPRIN2") 

dopamine_candidates = c(dopamine_synthesis, dopamine_receptors, growth_factors)

dopamine_candidates = dopamine_candidates[dopamine_candidates %in% rownames(neural_obj)]


Idents(neural_obj) = "Celltype_treatment"

DotPlot(
  neural_obj, 
  features = dopamine_candidates,
  group.by = "Celltype_treatment"
) + RotatedAxis()

DotPlot(
  neural_obj, 
  features = dopamine_synthesis,
  group.by = "Celltype_treatment"
) + RotatedAxis()





# CHECKPOINT 12: WGCNA for neurons

if (file.exists("neural_wgcna_obj.rds")) {
  
  neural_obj = readRDS("neural_wgcna_obj.rds")
  
} else { 
  
  neural_obj = SetupForWGCNA(
    neural_obj,
    gene_select = "fraction",
    fraction = 0.05,
    wgcna_name = "neural_WGCNA"
  )
  
  neural_obj = MetacellsByGroups(
    neural_obj,
    group.by = c("refined_annotation", "orig.ident"),
    k = 25, 
    max_shared = 10,
    ident.group = "refined_annotation"
  )
  
  neural_obj = NormalizeMetacells(neural_obj)
  
  metacell_obj = GetMetacellObject(neural_obj)
  
  valid_neurons = as.character(unique(metacell_obj$refined_annotation))
  
  neural_obj = SetDatExpr(
    neural_obj,
    group_name = valid_neurons,
    group.by = "refined_annotation",
    assay = "RNA"
  )
  
  neural_obj = TestSoftPowers(
    neural_obj,
    networkType = 'signed'
  )
  
  neural_obj = ConstructNetwork(
    neural_obj, 
    tom_name = 'neural_WGCNA',
    soft_power = 6, # Power from TestSoftPowers
    setDatExpr = FALSE,
    #  overwrite_tom = TRUE, # replacing old
    
    deepSplit = 4,         # split into smaller modules
    minModuleSize = 50,    # Prevents  from being too small
    mergeCutHeight = 0.05   # Merges modules that are >80% correlated
  )
  
  neural_obj = ModuleEigengenes(neural_obj)
  neural_obj = ModuleConnectivity(neural_obj)
  
  eigengenes_n = GetMEs(neural_obj)
  neural_obj = AddMetaData(neural_obj, metadata = eigengenes_n)
  
  saveRDS(neural_obj, "neural_wgcna_obj.rds")
}

PlotDendrogram(neural_obj, main = 'Neurons Brain Gene Network')

ModuleFeaturePlot(
  neural_obj,
  features = 'hMEs', # Module Eigengenes
  order = TRUE,       # cells with the highest expression to the front
  reduction = 'umap_by_pca'
)

head(GetModules(tunga_frog_merged_blast))


VlnPlot(
  neural_obj,
  features = c("turquoise", "yellow", "red", "green"), 
  group.by = "orig.ident", #
  pt.size = 0,
  ncol = 2
)


VlnPlot(
  neural_obj,
  features = c("greenyellow", "blue", "magenta", "tan"), 
  group.by = "orig.ident", #
  pt.size = 0,
  ncol = 2
)


VlnPlot(
  neural_obj,
  features = c("pink", "purple", "black", "brown"), 
  group.by = "orig.ident", #
  pt.size = 0,
  ncol = 2
)



top_yellow_hubs = modules %>%
  dplyr::filter(module == "yellow") %>%
  dplyr::arrange(desc(kME_yellow)) %>%
  head(30) %>%
  dplyr::pull(gene_name)

# 30 highest kme score
print(top_yellow_hubs)



dopamine_makers = c("TH", "DDC")

# Combine all candidate lists
full_candidate_list = unique(c(
  dopamine_makers,            # Dopamine makers
  dopamine_receptors,         # Dopamine receptors
  growth_factors,             # growth factors
  glutamate_candidates,       # glutamate candidates
  projection_genes,           # projection_genes: neural remodeling 
  tf_candidates              # Transcription factors genes
))


DotPlot(
  neural_obj, 
  features = c(full_candidate_list,   top_yellow_hubs),
  group.by = "Celltype_treatment"
) + 
  RotatedAxis() + 
  theme(
    axis.text.x = element_text(size = 3),
    axis.text.y = element_text(size = 9)
  ) +
  ggtitle("All candidate genes expression Neurons only")



DotPlot(
  neural_obj, 
  features = top_yellow_hubs 
#  group.by = "Celltype_treatment"
) + 
  RotatedAxis() + 
  theme(
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 9)
  ) +
  ggtitle("top_yellow_hubs expression Neurons only")


# Checking if my candidate genes are in the module
neural_modules = GetModules(neural_obj)

check = neural_modules %>%
  dplyr::filter(gene_name %in%  
                  c(dopamine_makers, dopamine_receptors, growth_factors, glutamate_candidates, projection_genes, tf_candidates)) %>%
  dplyr::select(gene_name, module, kME_yellow) %>%
  arrange(module)

print("Are validated targets in yellow module?")
print(check)


# Check change in gene expression from whine to whine chuck
yellow_gene_names = neural_modules %>% 
  dplyr::filter(module == "yellow") %>%
  dplyr::pull(gene_name) # turns column into a text list

yellow_stats = FindMarkers(
  tunga_frog_merged_blast, 
  ident.1 = "Whine-chuck", 
  ident.2 = "Whine", 
  features = yellow_gene_names,
  logfc.threshold = 0 # output results even if fold change is low
)

print(yellow_stats)
nrow(yellow_stats)
# Comparing the two modules. green only on gabas, yellow on all neurons
print(length(intersect(yellow_gene_names, green_gene_names)))


# Extracting top 50 most shifting Yellow genes but excluding LOCs
top_yellow = rownames(yellow_stats)[!grepl("^LOC", rownames(yellow_stats))]
top_yellow = head(top_yellow, 65)

DotPlot(
  neural_obj, 
  features = top_yellow, 
  group.by = "orig.ident"
) + 
  RotatedAxis() + 
  theme(axis.text.x = element_text(size = 9)) +
  ggtitle("top 65 yellow de genes. pvaljue adj less than 0.05")


# G profiler on those genes

# Extracting most significant shifting Yellow genes 
#(but excluding LOCs)
#significant_yellow = rownames(yellow_stats)[yellow_stats$p_val_adj < 0.05 & !grepl("^LOC", rownames(yellow_stats))]
significant_yellow = rownames(yellow_stats)[yellow_stats$p_val_adj < 0.05]

print(length(significant_yellow))
# 102

# Map to Xenopus protein
yellow_translation = translation_table %>%
  dplyr::filter(tungara_symbol %in% significant_yellow)

yellow_xenopus_ids = unique(yellow_translation$xenopus_protein)
yellow_xenopus_ids = yellow_xenopus_ids[!is.na(yellow_xenopus_ids)]

print(yellow_xenopus_ids, collapse = "\n")
print(length(yellow_xenopus_ids))
# 73


# Using all genes in module
# Map to Xenopus protein
all_yellow_translation = translation_table %>%
  dplyr::filter(tungara_symbol %in% yellow_gene_names)

all_yellow_xenopus_ids = unique(all_yellow_translation$xenopus_protein)
all_yellow_xenopus_ids = all_yellow_xenopus_ids[!is.na(all_yellow_xenopus_ids)]

cat(paste(all_yellow_xenopus_ids, collapse = "\n"))
print(length(all_yellow_xenopus_ids))





wc_neurons = subset(neural_obj, subset = orig.ident == "Whine-chuck")

global_mobilizers = c("SLC22A2", "ATP7B", "TRAF2")  # Yellow significiant

DotPlot(
  wc_neurons, 
  features = c(full_candidate_list, global_mobilizers), 
  group.by = "broad_annotation"
) + 
  RotatedAxis() + 
  theme(
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 12)
  ) +
  ggtitle("Candidate genes for all neurons: Dopamine makers Dopamine receptors Growth factors Glutamate candidates Projection_genes: neural remodeling Transcription factors genes
")

DotPlot(
  neural_obj, 
  features = c(full_candidate_list, global_mobilizers), 
  group.by = "Celltype_treatment"
 # group.by = "broad_annotation"
) + 
  RotatedAxis() + 
  theme(
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 12)
  ) +
  ggtitle("Candidate genes for all neurons: Dopamine makers Dopamine receptors Growth factors Glutamate candidates Projection_genes: neural remodeling Transcription factors genes
")

# Plot Yellow module activation by brain map
FeaturePlot(
  neural_obj,
  features = "yellow", # Yellow WGCNA score
  split.by = "orig.ident",
  reduction = "umap_by_pca",
  pt.size = 1,
  order = TRUE
)




yellow_response_genes = c("NRTN", "CAPRIN2", "NUMB", "RPGRIP1L", "SLC22A2", "ATP7B", "TRAF2")

DotPlot(
  neural_obj, 
  features = yellow_response_genes, 
  group.by = "Celltype_treatment"
) + 
  RotatedAxis() + 
  theme(axis.text.x = element_text(size = 10)) + coord_flip() +
  ggtitle("ONLY yellow top genes")


all_response_genes = c(yellow_response_genes, "SP4", "SMAD2", "E2F3", "MAPK8", "AVP")

DotPlot(
  neural_obj, 
  features = all_response_genes, 
  group.by = "Celltype_treatment"
) + 
  RotatedAxis() + 
  theme(axis.text.x = element_text(size = 10)) + coord_flip() +
  ggtitle("ALL top genes")










all_neurons_stats = FindMarkers(
  neural_obj,
  ident.1 = "Whine-chuck",
  ident.2 = "Whine",
  logfc.threshold = 0.25, # Strict threshold
  only.pos = TRUE         # Focus on turned ON by sound
)

# Top genes without WGCNA
all_neurons_stats$gene = rownames(all_neurons_stats)
top_20_neurons_genes = head(all_neurons_stats[order(-all_neurons_stats$avg_log2FC), ], 20)$gene
print(top_20_neurons_genes)

DotPlot(
  neural_obj,
  features = top_20_neurons_genes,
  group.by = "Celltype_treatment"
) + RotatedAxis() + ggtitle("Top 20 response genes all neurons") + coord_flip()



final_gene_list = translation_table %>%
  dplyr::filter(tungara_symbol %in% top_20_neurons_genes) %>%
  dplyr::filter(!is.na(xenopus_clean) & !grepl("^LOC", xenopus_clean)) %>%
  dplyr::pull(tungara_symbol) %>%
  unique()
print('SPTLC3"  "DSCC1"   "SMYD1"   "CAPRIN2" "GJD3"    "ARL3L2" ')

DotPlot(
  neural_obj,
  features = final_gene_list,
  group.by = "Celltype_treatment"
) + RotatedAxis() + ggtitle("Top 20 response genes all neurons with a xenopus correspondent") + coord_flip()


response_matrix = GetAssayData(neural_obj, assay = "RNA", slot = "data")[all_response_genes, ]

# Transpose and cor
gene_cor = cor(t(as.matrix(response_matrix)))

# Plot a heatmap of the genes that respond to the call
pheatmap(
  gene_cor,
  show_rownames = T,
  show_colnames = T,
  main = "Correlation matrix on response genes",
  clustering_method = "ward.D2"
)


wc_expr = GetAssayData(subset(neural_obj, orig.ident == "Whine-chuck"), assay = "RNA", slot = "data")[all_response_genes, ]
w_expr = GetAssayData(subset(neural_obj, orig.ident == "Whine"), assay = "RNA", slot = "data")[all_response_genes, ]

#Calculate the difference in correlation
pheatmap(cor(t(as.matrix(wc_expr))), main = "Whine-Chuck correlation ")
pheatmap(cor(t(as.matrix(w_expr))), main = "Whine correlation")



response_intensity = colMeans(GetAssayData(neural_obj, slot="data")[all_response_genes, ])
priming_score = GetAssayData(neural_obj, slot="data")["SMAD2", ]

plot(priming_score, response_intensity, 
     main="Does priming predict the response intensity?",
     xlab="Priming of SMAD2 level", ylab="Response intensity")











# Calculating grey percentage in neural object
grey_genes = sum(neural_modules$module == "grey")
total = nrow(neural_modules)
grey_percentage = (grey_genes / total) * 100
print(round(grey_percentage, 2))
