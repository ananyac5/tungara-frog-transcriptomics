# Neural Gene Expression and Co-expression Networks in Túngara Frogs

### Single-Nucleus RNA Sequencing Analysis | University of Texas at Austin

## Project Overview

This undergraduate computational biology research project investigates how different acoustic mating calls influence gene expression in the brains of female túngara frogs (*Engystomops pustulosus*).

Using single-nucleus RNA sequencing (snRNA-seq), our research team analyzed brain cell populations, identified differentially expressed genes, and explored gene co-expression networks associated with responses to mating calls.

The project was conducted as part of the University of Texas at Austin's Big Data in Biology research program in collaboration with the Hofmann Lab.

## Research Objectives

The primary research question was:

**How do different acoustic mating calls influence gene expression and neural activity-related molecular pathways in female túngara frogs?**

We investigated three experimental conditions:

- **Whine-Chuck (WC):** More attractive mating call.
- **Whine (W):** Less attractive mating call.
- **Green Noise (GN):** Control acoustic stimulus.

Our objectives included:

1. Identifying and annotating distinct brain cell populations.
2. Comparing cell-type distributions across experimental conditions.
3. Identifying differentially expressed genes associated with acoustic stimuli.
4. Investigating transcriptional responses in GABAergic and excitatory neurons.
5. Identifying co-expressed gene networks and potential hub genes.
6. Exploring biological pathways associated with neural responses and mate-choice behavior.

## Technologies and Tools

**Programming and Analysis**
- R
- Seurat
- dplyr
- ggplot2
- patchwork
- pheatmap

**Bioinformatics and Statistical Methods**
- Single-nucleus RNA sequencing
- Principal Component Analysis (PCA)
- Uniform Manifold Approximation and Projection (UMAP)
- SCtype cell-type annotation
- Cross-species ortholog mapping using BLAST and OrthoFinder results
- Differential gene expression analysis
- Gene Ontology (GO) enrichment analysis
- hdWGCNA gene co-expression network analysis
- Hierarchical clustering
- MapMyCells annotation comparison

## Methodology

### 1. Data Preprocessing and Integration

Single-nucleus RNA-seq datasets from three acoustic treatment conditions were processed using Seurat.

The workflow included normalization, variable feature selection, scaling, dimensionality reduction, and clustering.

Seurat's integration framework was used to examine cell populations across experimental conditions.

### 2. Cell-Type Annotation

Brain cell populations were annotated using marker gene expression and computational annotation methods.

SCtype and reference-based approaches were used to investigate cell identities. Cross-species gene mapping helped connect túngara frog genes with reference annotations from *Xenopus* species.

UMAP visualizations were used to examine cell clusters and their predicted identities.

### 3. Differential Gene Expression Analysis

Differential gene expression analyses compared transcriptional responses across acoustic treatments, with particular attention to Whine-Chuck versus Whine.

Analyses focused on GABAergic and excitatory neuronal populations.

Adjusted p-values and log2 fold changes were used to investigate genes showing differences in expression.

Heatmaps and additional visualizations were generated to examine gene expression patterns.

### 4. Functional Enrichment Analysis

Gene Ontology enrichment analysis was used to investigate biological processes associated with selected differentially expressed genes and candidate gene sets.

The analyses explored potential relationships between transcriptional responses, cellular signaling, and neural function.

### 5. Gene Co-expression Network Analysis

High-dimensional Weighted Gene Co-expression Network Analysis (hdWGCNA) was used to identify groups of co-expressed genes.

The research examined gene modules, module eigengenes, and highly connected candidate hub genes.

Network visualizations and correlation analyses were used to explore relationships among genes within neuronal populations.

## Key Findings

Our analyses identified distinct brain cell populations and transcriptional differences associated with acoustic treatment conditions.

The group's findings included:

- Identification and annotation of multiple neuronal and non-neuronal cell populations.
- Differential gene expression patterns in GABAergic and excitatory neurons.
- Candidate genes including **CAPRIN2, NRTN, and SLC22A2**.
- Gene co-expression modules associated with biological processes including cellular signaling and neural development.
- Potential molecular pathways relevant to neural responses during mate-choice behavior.

These results identify candidate genes and pathways for further investigation rather than establishing direct causal mechanisms.

## Research Visualizations

Selected visualizations are available in the `figures/` directory.

They include:

- Cell-type annotation UMAP
- GABAergic differential expression heatmap
- Gene Ontology enrichment visualization
- hdWGCNA module dendrogram
- Green module hub gene network

The complete research presentation is available in the `presentation/` directory.

## Repository Structure

```text
Tungara-Frog-Research/
├── README.md
├── code/
│   ├── 01_single_cell_DEG_analysis.R
│   ├── 02_gene_expression_ortholog_comparison.R
│   ├── 03_gene_coexpression_network_analysis.R
│   └── 04_cell_integration_annotation.R
├── figures/
├── presentation/
│   └── Tungara_Frog_Research_Presentation.pdf
└── documentation/
```

### Code Descriptions

| Script | Description |
|---|---|
| `01_single_cell_DEG_analysis.R` | Single-cell preprocessing, clustering, annotation, differential expression, heatmaps, and GO enrichment. |
| `02_gene_expression_ortholog_comparison.R` | Gene expression feature plots and comparison of BLAST and OrthoFinder mapping results. |
| `03_gene_coexpression_network_analysis.R` | hdWGCNA analysis, co-expression modules, hub gene networks, and candidate gene exploration. |
| `04_cell_integration_annotation.R` | Seurat integration, clustering parameter exploration, SCtype annotation, and MapMyCells comparison. |

The scripts are research workflows and exploratory analyses. Some depend on preprocessed Seurat objects, reference files, and external data not included in this repository.

## Research Team

**University of Texas at Austin — Big Data in Biology, Summer 2026**

Research team:
- Ananya Chintalapati
- Jace Ramo
- Hedi El Matri
- Sakthi Priya Vijaya Basker

Research conducted in collaboration with the Hofmann Lab at UT Austin.

### Individual Contributions

*To be completed with a description of individual responsibilities, analyses performed, scripts developed, and contributions to the research presentation.*

This repository includes collaborative research work. Authorship and contributions should be attributed appropriately.

## Limitations

The study included relatively small numbers of biological samples for certain treatment conditions.

Cross-species gene annotation presented challenges due to incomplete mapping, gene naming differences, and limited reference resources.

Some scripts rely on research computing infrastructure, previously generated objects, or datasets that are not distributed with this repository.

## Data Availability and Reproducibility

Raw sequencing data and large intermediate analysis objects are not included in this repository.

The code is provided to document research methods and analytical workflows. It has not been packaged as a fully self-contained, independently reproducible pipeline.

The availability of research data, code, and results is subject to laboratory permissions and applicable research-sharing policies.

## Acknowledgments

We acknowledge the Hofmann Lab and the Big Data in Biology research program at the University of Texas at Austin for providing research guidance, computational resources, and opportunities to participate in computational biology research.

---

**Research Area:** Computational Biology | Bioinformatics | Single-Cell Genomics | Data Science
