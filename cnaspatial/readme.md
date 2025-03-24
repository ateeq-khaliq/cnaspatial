# CNASpatial

## Overview

CNASpatial is an R package for analyzing Copy Number Alterations (CNAs) in spatial transcriptomics data. The package provides tools to categorize cells based on functional signatures like hypoxia and proliferation, detect CNAs using CopyKAT, and visualize CNAs in relation to these functional categories.

## Features

- Calculate hypoxia and proliferation scores for cells in spatial data
- Categorize cells based on their functional signatures
- Create subsets for CopyKAT analysis
- Run CopyKAT to infer CNAs from gene expression data
- Visualize CNAs in relation to functional categories using:
  - Heatmaps
  - Chromosome ideograms
  - Circos plots
  - Barplots
- Generate comprehensive Excel reports summarizing CNA findings
- Compare CNA profiles between different treatment conditions

## Installation

You can install the development version of CNASpatial from GitHub:

```r
# Install dependencies first
install.packages(c("Seurat", "dplyr", "ggplot2", "pheatmap", "RColorBrewer", 
                  "patchwork", "gridExtra", "grid", "writexl", "tidyr", 
                  "circlize", "ComplexHeatmap"))

# Install msigdbr for gene signatures
install.packages("msigdbr")

# Install copykat
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("copykat")

# Install CNASpatial from GitHub
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("yourusername/cnaspatial")
```

## Quick Start

```r
library(cnaspatial)

# Load a Seurat object with spatial data
# seurat_obj <- ...

# Run the complete pipeline
results <- runCNASpatialPipeline(
    seurat_obj, 
    output_dir = "results", 
    sample_name = "my_sample"
)

# Access the results
seurat_obj_analyzed <- results$seurat_obj
copykat_results <- results$copykat_results
cna_profiles <- results$cna_profiles
```

## Workflow

The typical workflow consists of:

1. Calculate functional scores (hypoxia and proliferation)
2. Categorize cells based on these scores
3. Create a subset for CopyKAT analysis
4. Run CopyKAT to detect CNAs
5. Add CNA data to the Seurat object
6. Create simplified functional classes
7. Extract and analyze CNA profiles by functional category
8. Visualize results using various plotting functions
9. Generate a comprehensive Excel report

## Example Usage

### Basic Pipeline

```r
library(cnaspatial)
library(Seurat)

# Load your Seurat object with spatial data
# For example:
# seurat_obj <- readRDS("my_spatial_data.rds")

# Run the complete pipeline
results <- runCNASpatialPipeline(
    seurat_obj,
    output_dir = "results/my_sample",
    sample_name = "sample1",
    run_copykat = TRUE,   # Set to FALSE if you want to skip CopyKAT analysis
    max_cells = 20000     # Maximum number of cells for CopyKAT
)

# Access the resulting Seurat object with all annotations
seurat_obj_analyzed <- results$seurat_obj

# Visualize spatial distribution of functional categories
SpatialDimPlot(seurat_obj_analyzed, group.by = "functional_class")

# Visualize CopyKAT predictions
SpatialDimPlot(seurat_obj_analyzed, group.by = "copykat_prediction")
```

### Step-by-Step Approach

If you prefer more control over each step:

```r
library(cnaspatial)
library(Seurat)

# 1. Calculate functional scores
seurat_obj <- calculateFunctionalScores(seurat_obj)

# 2. Categorize cells by function
seurat_obj <- categorizeCellsByFunction(seurat_obj)

# 3. Create expanded subset for CopyKAT
seurat_obj <- createCopyKATSubset(seurat_obj)

# 4. Prepare count matrix for CopyKAT
copykat_matrix <- prepareCopyKATMatrix(
    seurat_obj,
    output_file = "copykat_input.txt"
)

# 5. Run CopyKAT
copykat_results <- runCopyKAT(
    copykat_matrix,
    output_dir = "copykat_results",
    sample_name = "my_sample"
)

# 6. Add CopyKAT results to Seurat object
seurat_obj <- addCopyKATResults(seurat_obj, copykat_results)

# 7. Create simplified functional classes
seurat_obj <- createSimplifiedClasses(seurat_obj)

# 8. Extract and analyze CNA data
cna_data <- extractCNAData(seurat_obj)
cna_profiles <- calculateCNAProfiles(cna_data)
arm_summary <- aggregateByChromosomeArm(cna_profiles)

# 9. Visualize results
plotCNAHeatmap(arm_summary)
plotChromosomeIdeogram(cna_profiles)
plotCNABarplot(arm_summary)
plotCNACircos(arm_summary)

# 10. Create Excel summary
createCNASummaryExcel(arm_summary, output_file = "CNA_summary.xlsx")
```

### Treatment Comparison

If you have treated and untreated samples:

```r
# Assuming you have a 'treated' column in your Seurat metadata
# with values "Yes" or "No"
plotTreatmentComparison(
    seurat_obj,
    treatment_col = "treated",
    output_dir = "treatment_results"
)
```

## Visualization Examples

### Functional Categories

![Functional Categories](docs/images/functional_classes.png)

### CNA Heatmap

![CNA Heatmap](docs/images/cna_heatmap.png)

### Chromosome Ideogram

![Chromosome Ideogram](docs/images/chromosome_ideogram.png)

### Treatment Comparison

![Treatment Comparison](docs/images/treatment_comparison.png)

## Citation

If you use CNASpatial in your research, please cite:

```
[Your citation information will go here]
```

## License

This package is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a pull request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request
