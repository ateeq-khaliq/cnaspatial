#' Calculate Hypoxia and Proliferation Scores
#'
#' This function calculates hypoxia and proliferation scores for cells in a
#' Seurat object based on predefined gene signatures.
#'
#' @param seurat_obj A Seurat object containing gene expression data
#' @param assay The assay to use for calculations. Default is "Spatial"
#' @param species Species for the gene set. Default is "Homo sapiens"
#' @param hypoxia_geneset_name Name of the MSigDB hypoxia gene set. Default is "HALLMARK_HYPOXIA"
#' @param prolif_geneset_names Character vector of MSigDB proliferation gene set names.
#'        Default is c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2")
#' @param hypoxia_score_name Name to use for the hypoxia score. Default is "HypoxiaScore"
#' @param prolif_score_name Name to use for the proliferation score. Default is "ProliferationScore"
#'
#' @return A Seurat object with hypoxia and proliferation scores added
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_obj <- calculateFunctionalScores(seurat_obj)
#' }
calculateFunctionalScores <- function(seurat_obj, 
                                    assay = "Spatial",
                                    species = "Homo sapiens",
                                    hypoxia_geneset_name = "HALLMARK_HYPOXIA",
                                    prolif_geneset_names = c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"),
                                    hypoxia_score_name = "HypoxiaScore",
                                    prolif_score_name = "ProliferationScore") {
  
  # Check if required packages are available
  if (!requireNamespace("msigdbr", quietly = TRUE) || 
      !requireNamespace("Seurat", quietly = TRUE)) {
    stop("Required packages msigdbr and Seurat must be installed")
  }
  
  # Get hypoxia gene set
  hallmark_hypoxia <- msigdbr::msigdbr(species = species, category = "H") %>%
    dplyr::filter(gs_name == hypoxia_geneset_name) %>%
    dplyr::pull(gene_symbol)
  
  # Get proliferation gene sets (MYC targets)
  myc_targets <- msigdbr::msigdbr(species = species, category = "H") %>%
    dplyr::filter(gs_name %in% prolif_geneset_names) %>%
    dplyr::pull(gene_symbol)
  
  # Check that genes exist in the dataset
  genes_in_data <- rownames(seurat_obj)
  hypoxia_genes_present <- hallmark_hypoxia[hallmark_hypoxia %in% genes_in_data]
  prolif_genes_present <- myc_targets[myc_targets %in% genes_in_data]
  
  if (length(hypoxia_genes_present) < 5) {
    warning("Few hypoxia signature genes found in the dataset. Scores may be unreliable.")
  }
  if (length(prolif_genes_present) < 5) {
    warning("Few proliferation signature genes found in the dataset. Scores may be unreliable.")
  }
  
  # Set assay
  Seurat::DefaultAssay(seurat_obj) <- assay
  
  # Calculate scores
  seurat_obj <- Seurat::AddModuleScore(seurat_obj, 
                                     features = list(hypoxia_genes_present), 
                                     name = hypoxia_score_name)
  
  seurat_obj <- Seurat::AddModuleScore(seurat_obj, 
                                     features = list(prolif_genes_present), 
                                     name = prolif_score_name)
  
  return(seurat_obj)
}

#' Categorize Cells Based on Functional Scores
#'
#' This function categorizes cells in a Seurat object based on their hypoxia and
#' proliferation scores.
#'
#' @param seurat_obj A Seurat object with hypoxia and proliferation scores
#' @param hypoxia_score_col Name of column containing hypoxia scores. Default is "HypoxiaScore1"
#' @param prolif_score_col Name of column containing proliferation scores. Default is "ProliferationScore1"
#' @param high_quantile Quantile threshold for "high" designation. Default is 0.7
#' @param low_quantile Quantile threshold for "low" designation. Default is 0.3
#' @param category_col Name for the category column. Default is "HypoxProlifCategory"
#' @param refined_category_col Name for the refined category column. Default is "HypoxProlifCategory_Refined"
#'
#' @return A Seurat object with category designations added to metadata
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_obj <- categorizeCellsByFunction(seurat_obj)
#' }
categorizeCellsByFunction <- function(seurat_obj,
                                     hypoxia_score_col = "HypoxiaScore1",
                                     prolif_score_col = "ProliferationScore1",
                                     high_quantile = 0.7,
                                     low_quantile = 0.3,
                                     category_col = "HypoxProlifCategory",
                                     refined_category_col = "HypoxProlifCategory_Refined") {
  
  if (!all(c(hypoxia_score_col, prolif_score_col) %in% colnames(seurat_obj@meta.data))) {
    stop("Hypoxia and/or proliferation score columns not found in Seurat metadata")
  }
  
  # Extract scores
  hypoxia_scores <- seurat_obj@meta.data[[hypoxia_score_col]]
  prolif_scores <- seurat_obj@meta.data[[prolif_score_col]]
  
  # Define thresholds for high and low scores
  hypoxia_high_threshold <- quantile(hypoxia_scores, high_quantile, na.rm = TRUE)
  prolif_high_threshold <- quantile(prolif_scores, high_quantile, na.rm = TRUE)
  hypoxia_low_threshold <- quantile(hypoxia_scores, low_quantile, na.rm = TRUE)
  prolif_low_threshold <- quantile(prolif_scores, low_quantile, na.rm = TRUE)
  
  # Print thresholds for reference
  message(paste("Hypoxia high threshold:", round(hypoxia_high_threshold, 3)))
  message(paste("Hypoxia low threshold:", round(hypoxia_low_threshold, 3)))
  message(paste("Proliferation high threshold:", round(prolif_high_threshold, 3)))
  message(paste("Proliferation low threshold:", round(prolif_low_threshold, 3)))
  
  # Create a new column for categorization
  seurat_obj@meta.data[[category_col]] <- NA
  
  # Categorize cells
  seurat_obj@meta.data[[category_col]][
    seurat_obj@meta.data[[hypoxia_score_col]] > hypoxia_high_threshold & 
    seurat_obj@meta.data[[prolif_score_col]] > prolif_high_threshold
    ] <- "High_Both"
  
  seurat_obj@meta.data[[category_col]][
    seurat_obj@meta.data[[hypoxia_score_col]] > hypoxia_high_threshold & 
    seurat_obj@meta.data[[prolif_score_col]] <= prolif_high_threshold & 
    seurat_obj@meta.data[[prolif_score_col]] > prolif_low_threshold
    ] <- "High_Hypoxia"
  
  seurat_obj@meta.data[[category_col]][
    seurat_obj@meta.data[[prolif_score_col]] > prolif_high_threshold & 
    seurat_obj@meta.data[[hypoxia_score_col]] <= hypoxia_high_threshold &
    seurat_obj@meta.data[[hypoxia_score_col]] > hypoxia_low_threshold
    ] <- "High_Proliferation"
  
  seurat_obj@meta.data[[category_col]][
    seurat_obj@meta.data[[hypoxia_score_col]] <= hypoxia_low_threshold & 
    seurat_obj@meta.data[[prolif_score_col]] <= prolif_low_threshold
    ] <- "Undefined"
  
  # Assign "Medium" to any remaining uncategorized cells
  seurat_obj@meta.data[[category_col]][is.na(seurat_obj@meta.data[[category_col]])] <- "Medium"
  
  # Convert to factor with specific level order for better visualization
  seurat_obj@meta.data[[category_col]] <- factor(
    seurat_obj@meta.data[[category_col]], 
    levels = c("High_Both", "High_Hypoxia", "High_Proliferation", "Medium", "Undefined")
  )
  
  # Save the current categorization
  seurat_obj@meta.data[[paste0(category_col, "_Original")]] <- seurat_obj@meta.data[[category_col]]
  
  # Create refined category based on dominant signature
  seurat_obj@meta.data[[refined_category_col]] <- as.character(seurat_obj@meta.data[[category_col]])
  
  # For High_Both cells, compare the relative strength of each signal
  high_both_cells <- which(seurat_obj@meta.data[[category_col]] == "High_Both")
  
  if (length(high_both_cells) > 0) {
    for (i in high_both_cells) {
      # Calculate each score's percentile rank within its own distribution
      hypoxia_percentile <- ecdf(hypoxia_scores)(hypoxia_scores[i])
      prolif_percentile <- ecdf(prolif_scores)(prolif_scores[i])
      
      # Assign to the category with the higher relative signal
      if (hypoxia_percentile > prolif_percentile) {
        seurat_obj@meta.data[[refined_category_col]][i] <- "High_Hypoxia_Dominant"
      } else {
        seurat_obj@meta.data[[refined_category_col]][i] <- "High_Proliferation_Dominant"
      }
    }
  }
  
  # Convert to factor with specific level order
  seurat_obj@meta.data[[refined_category_col]] <- factor(
    seurat_obj@meta.data[[refined_category_col]],
    levels = c(
      "High_Hypoxia_Dominant", "High_Proliferation_Dominant",
      "High_Hypoxia", "High_Proliferation",
      "Medium", "Undefined"
    )
  )
  
  # Print the distribution of the refined categories
  refined_cat_table <- table(seurat_obj@meta.data[[refined_category_col]])
  message("Refined category counts:")
  print(refined_cat_table)
  
  return(seurat_obj)
}

#' Create CopyKAT Ready Subset
#'
#' Subsets a Seurat object to include high-interest cells for CopyKAT analysis
#' and optionally includes additional cells to reach a target count.
#'
#' @param seurat_obj A Seurat object with functional categories
#' @param category_col Column name with refined functional categories. Default is "HypoxProlifCategory_Refined"
#' @param target_cells Target number of cells for CopyKAT analysis. Default is 20000
#' @param high_interest_categories Vector of category names to include. Default is 
#'        c("High_Hypoxia_Dominant", "High_Hypoxia", "High_Proliferation_Dominant", "High_Proliferation")
#' @param output_col Name for the expanded subset column. Default is "ExpandedHighSubset"
#' @param hypoxia_score_col Name of column containing hypoxia scores. Default is "HypoxiaScore1"
#' @param prolif_score_col Name of column containing proliferation scores. Default is "ProliferationScore1"
#'
#' @return A Seurat object with expanded subset column added to metadata
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_obj <- createCopyKATSubset(seurat_obj)
#' subset_obj <- subset(seurat_obj, subset = ExpandedHighSubset %in% 
#'                      c("High_Hypoxia_Dominant", "High_Hypoxia", 
#'                        "High_Proliferation_Dominant", "High_Proliferation", "Medium_High"))
#' }
createCopyKATSubset <- function(seurat_obj,
                               category_col = "HypoxProlifCategory_Refined",
                               target_cells = 20000,
                               high_interest_categories = c("High_Hypoxia_Dominant", "High_Hypoxia", 
                                                         "High_Proliferation_Dominant", "High_Proliferation"),
                               output_col = "ExpandedHighSubset",
                               hypoxia_score_col = "HypoxiaScore1",
                               prolif_score_col = "ProliferationScore1") {
  
  if (!category_col %in% colnames(seurat_obj@meta.data)) {
    stop("Category column not found in Seurat metadata")
  }
  
  # Get current counts from categories
  current_counts <- table(seurat_obj@meta.data[[category_col]])
  message("Current counts in the categories:")
  print(current_counts)
  
  # Calculate how many cells are in high-interest categories
  high_category_counts <- current_counts[names(current_counts) %in% high_interest_categories]
  total_high_count <- sum(high_category_counts)
  message(paste("Total count for all high categories:", total_high_count))
  
  # Calculate how many additional cells needed to reach target
  spots_to_add <- min(target_cells - total_high_count, 
                     ncol(seurat_obj) - total_high_count)  # Can't add more than available
  message(paste("Need to add approximately", spots_to_add, "spots"))
  
  # Create a new column for the expanded subset
  seurat_obj@meta.data[[output_col]] <- as.character(seurat_obj@meta.data[[category_col]])
  
  # Mark non-high-interest categories as "Other"
  seurat_obj@meta.data[[output_col]][
    !(seurat_obj@meta.data[[category_col]] %in% high_interest_categories)
  ] <- "Other"
  
  if (spots_to_add > 0) {
    # Identify "Other" cells
    other_indices <- which(seurat_obj@meta.data[[output_col]] == "Other")
    
    if (length(other_indices) > 0) {
      # Get scores for these cells
      other_hypoxia_scores <- seurat_obj@meta.data[[hypoxia_score_col]][other_indices]
      other_prolif_scores <- seurat_obj@meta.data[[prolif_score_col]][other_indices]
      
      # Calculate a combined score (weighted average)
      combined_scores <- (other_hypoxia_scores + other_prolif_scores) / 2
      
      # Create a data frame to help with selection
      selection_df <- data.frame(
        index = other_indices,
        hypoxia = other_hypoxia_scores,
        proliferation = other_prolif_scores,
        combined = combined_scores,
        stringsAsFactors = FALSE
      )
      
      # Sort by combined score (descending)
      selection_df <- selection_df[order(-selection_df$combined), ]
      
      # Select top cells to reach desired total
      top_count <- min(spots_to_add, nrow(selection_df))
      top_indices <- selection_df$index[1:top_count]
      
      # Update the expanded subset column with selected cells
      seurat_obj@meta.data[[output_col]][top_indices] <- "Medium_High"
    }
  }
  
  # Convert to factor with ordered levels
  seurat_obj@meta.data[[output_col]] <- factor(
    seurat_obj@meta.data[[output_col]], 
    levels = c(high_interest_categories, "Medium_High", "Other")
  )
  
  # Check the final counts
  expanded_counts <- table(seurat_obj@meta.data[[output_col]])
  message("Final counts for expanded subset categories:")
  print(expanded_counts)
  
  # Calculate total for CopyKAT
  expanded_total <- sum(expanded_counts[names(expanded_counts) != "Other"])
  message(paste("Total cells for CopyKAT analysis:", expanded_total))
  
  return(seurat_obj)
}

#' Create Count Matrix for CopyKAT Analysis
#'
#' This function creates a count matrix suitable for CopyKAT analysis from a Seurat object.
#'
#' @param seurat_obj A Seurat object with cells to analyze
#' @param assay The assay to use for the count matrix. Default is "Spatial"
#' @param subset_column Column in metadata to use for subsetting. Default is "ExpandedHighSubset"
#' @param include_categories Categories to include in the analysis. Default is 
#'        c("High_Hypoxia_Dominant", "High_Hypoxia", "High_Proliferation_Dominant", 
#'          "High_Proliferation", "Medium_High")
#' @param output_file Optional path to save the count matrix as a tab-delimited file
#'
#' @return A matrix of raw counts for CopyKAT analysis
#' @export
#'
#' @examples
#' \dontrun{
#' copykat_matrix <- prepareCopyKATMatrix(seurat_obj, output_file = "copykat_input.txt")
#' }
prepareCopyKATMatrix <- function(seurat_obj,
                                assay = "Spatial",
                                subset_column = "ExpandedHighSubset",
                                include_categories = c("High_Hypoxia_Dominant", "High_Hypoxia", 
                                                    "High_Proliferation_Dominant", "High_Proliferation", 
                                                    "Medium_High"),
                                output_file = NULL) {
  
  # Check if subset column exists
  if (!subset_column %in% colnames(seurat_obj@meta.data)) {
    stop("Subset column not found in Seurat metadata")
  }
  
  # Create subset if requested
  if (!is.null(subset_column) && !is.null(include_categories)) {
    subset_obj <- subset(seurat_obj, subset = get(subset_column) %in% include_categories)
    message(paste("Created subset with", ncol(subset_obj), "cells"))
  } else {
    subset_obj <- seurat_obj
    message("Using all cells for CopyKAT input")
  }
  
  # Set the assay
  Seurat::DefaultAssay(subset_obj) <- assay
  
  # Get count matrix
  copykat_input <- Seurat::GetAssayData(subset_obj, slot = "counts")
  
  # Write to file if requested
  if (!is.null(output_file)) {
    utils::write.table(copykat_input, file = output_file, sep = "\t", quote = FALSE)
    message(paste("Wrote count matrix to", output_file))
  }
  
  return(copykat_input)
}

#' Run CopyKAT Analysis
#'
#' This function runs the CopyKAT algorithm to detect CNAs from single-cell or spatial data.
#'
#' @param count_matrix A matrix of raw counts (genes in rows, cells in columns)
#' @param output_dir Directory to save CopyKAT outputs. Default is "copykat_results"
#' @param sample_name Sample name for output files. Default is "sample"
#' @param genome Genome build. Default is "hg20"
#' @param cell_line Whether data is from cell lines. Default is "no"
#' @param min_genes_per_chr Minimum number of genes per chromosome to use. Default is 5
#' @param window_size Window size of genes for CNA calculations. Default is 100
#' @param ks_cutoff KS test cutoff. Default is 0.1
#' @param distance Distance metric for clustering. Default is "euclidean"
#' @param n_cores Number of cores to use. Default is 1
#' @param reference_cells Optional vector of cell names to use as normal reference
#' @param max_cells Maximum number of cells to analyze (will randomly sample if exceeded). Default is 20000
#'
#' @return CopyKAT results object
#' @export
#'
#' @examples
#' \dontrun{
#' copykat_results <- runCopyKAT(copykat_matrix, output_dir = "results/copykat")
#' }
runCopyKAT <- function(count_matrix,
                      output_dir = "copykat_results",
                      sample_name = "sample",
                      genome = "hg20",
                      cell_line = "no",
                      min_genes_per_chr = 5,
                      window_size = 100,
                      ks_cutoff = 0.1,
                      distance = "euclidean",
                      n_cores = 1,
                      reference_cells = NULL,
                      max_cells = 20000) {
  
  # Check if copykat is installed
  if (!requireNamespace("copykat", quietly = TRUE)) {
    stop("Package copykat must be installed")
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message(paste("Created output directory:", output_dir))
  }
  
  # Sample cells if there are too many
  if (ncol(count_matrix) > max_cells) {
    message(paste("Sampling", max_cells, "cells from", ncol(count_matrix), "total cells"))
    set.seed(42)  # For reproducibility
    sampled_cells <- sample(colnames(count_matrix), max_cells)
    count_matrix_subset <- count_matrix[, sampled_cells]
  } else {
    count_matrix_subset <- count_matrix
  }
  
  # Save current directory and change to output directory (copykat writes files to current dir)
  current_dir <- getwd()
  setwd(output_dir)
  
  # Run copykat
  message("Running copykat analysis (this may take a while)...")
  copykat_results <- copykat::copykat(
    rawmat = count_matrix_subset,
    id.type = "S",  # Use gene symbol as id
    cell.line = cell_line,
    ngene.chr = min_genes_per_chr,
    win.size = window_size,
    KS.cut = ks_cutoff,
    sam.name = sample_name,
    distance = distance,
    norm.cell.names = reference_cells,
    output.seg = FALSE,
    plot.genes = TRUE,
    genome = genome,
    n.cores = n_cores
  )
  
  # Return to original directory
  setwd(current_dir)
  
  # Save copykat results object
  saveRDS(copykat_results, file = file.path(output_dir, paste0(sample_name, "_results.rds")))
  message(paste("Saved copykat results to", file.path(output_dir, paste0(sample_name, "_results.rds"))))
  
  return(copykat_results)
}

#' Add CopyKAT Results to a Seurat Object
#'
#' This function adds CopyKAT predictions and CNA values to a Seurat object.
#'
#' @param seurat_obj A Seurat object containing cells that were analyzed with CopyKAT
#' @param copykat_results A CopyKAT results object from runCopyKAT()
#' @param prediction_col Name for the column to store CopyKAT predictions. Default is "copykat_prediction"
#' @param has_cna_col Name for the column indicating if a cell has CNA data. Default is "has_cna_data"
#' @param min_bound Minimum bound for CNA values (for scaling). Default is -2.6
#' @param max_bound Maximum bound for CNA values (for scaling). Default is 2.6
#'
#' @return A Seurat object with CopyKAT results added
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_obj <- addCopyKATResults(seurat_obj, copykat_results)
#' }
addCopyKATResults <- function(seurat_obj,
                             copykat_results,
                             prediction_col = "copykat_prediction",
                             has_cna_col = "has_cna_data",
                             min_bound = -2.6,
                             max_bound = 2.6) {
  
  # Extract the CNA matrix without the metadata columns
  cna_mat_cells <- as.matrix(copykat_results$CNAmat[, -(1:3)])
  
  # Convert column names from dots to dashes to match Seurat object
  colnames(cna_mat_cells) <- gsub("\\.", "-", colnames(cna_mat_cells))
  
  # Check for common cells between CNA matrix and Seurat object
  common_cells <- intersect(colnames(cna_mat_cells), colnames(seurat_obj))
  message(paste("Number of cells with CNA data:", length(common_cells)))
  
  # Apply bounds to CNA values
  message("Applying post-processing to CNA values...")
  bounded_cna <- pmin(pmax(cna_mat_cells, min_bound), max_bound)
  
  # Subset to cells present in the Seurat object
  bounded_cna_subset <- bounded_cna[, common_cells, drop = FALSE]
  
  # Create CNA assay and add to Seurat object
  cna_assay <- Seurat::CreateAssayObject(data = bounded_cna_subset)
  seurat_obj[["cna"]] <- cna_assay
  message("Added CNA data to Seurat object")
  
  # Add CopyKAT predictions to Seurat metadata
  seurat_obj@meta.data[[prediction_col]] <- "unknown"
  cells_with_pred <- intersect(copykat_results$prediction$cell.names, colnames(seurat_obj))
  
  seurat_obj@meta.data[[prediction_col]][match(cells_with_pred, colnames(seurat_obj))] <- 
    copykat_results$prediction$copykat.pred[match(cells_with_pred, copykat_results$prediction$cell.names)]
  
  # Create a flag for cells with CNA data
  seurat_obj@meta.data[[has_cna_col]] <- colnames(seurat_obj) %in% common_cells
  
  # Report summary statistics
  pred_table <- table(seurat_obj@meta.data[[prediction_col]])
  message("CopyKAT prediction counts:")
  print(pred_table)
  
  return(seurat_obj)
}

#' Create Simplified Functional Classes
#'
#' Simplifies the detailed functional categories into three main groups for easier analysis.
#'
#' @param seurat_obj A Seurat object with detailed functional categories
#' @param detailed_category_col Name of the column with detailed categories. Default is "ExpandedHighSubset"
#' @param output_col Name for the simplified functional class column. Default is "functional_class"
#' @param hypoxic_categories Categories to classify as "Hypoxic". Default is 
#'        c("High_Hypoxia_Dominant", "High_Hypoxia")
#' @param proliferative_categories Categories to classify as "Proliferative". Default is 
#'        c("High_Proliferation_Dominant", "High_Proliferation")
#'
#' @return A Seurat object with simplified functional class column added to metadata
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_obj <- createSimplifiedClasses(seurat_obj)
#' }
createSimplifiedClasses <- function(seurat_obj,
                                  detailed_category_col = "ExpandedHighSubset",
                                  output_col = "functional_class",
                                  hypoxic_categories = c("High_Hypoxia_Dominant", "High_Hypoxia"),
                                  proliferative_categories = c("High_Proliferation_Dominant", "High_Proliferation")) {
  
  if (!detailed_category_col %in% colnames(seurat_obj@meta.data)) {
    stop("Detailed category column not found in Seurat metadata")
  }
  
  # Display the existing categories and their counts
  message("Existing functional categories:")
  print(table(seurat_obj@meta.data[[detailed_category_col]]))
  
  # Create a new classification with three categories
  seurat_obj@meta.data[[output_col]] <- "Undefined"
  
  # Map hypoxic categories
  seurat_obj@meta.data[[output_col]][
    seurat_obj@meta.data[[detailed_category_col]] %in% hypoxic_categories
  ] <- "Hypoxic"
  
  # Map proliferative categories
  seurat_obj@meta.data[[output_col]][
    seurat_obj@meta.data[[detailed_category_col]] %in% proliferative_categories
  ] <- "Proliferative"
  
  # Convert to factor for better plotting
  seurat_obj@meta.data[[output_col]] <- factor(
    seurat_obj@meta.data[[output_col]],
    levels = c("Hypoxic", "Proliferative", "Undefined")
  )
  
  # Display the new classification
  message("Mapped to three categories:")
  print(table(seurat_obj@meta.data[[output_col]]))
  
  return(seurat_obj)
}

#' Extract CNA Data for Analysis
#'
#' Extracts and processes CNA data from a Seurat object for further analysis.
#'
#' @param seurat_obj A Seurat object with CNA data in a dedicated assay
#' @param assay_name Name of the assay containing CNA data. Default is "cna"
#' @param functional_class_col Name of the column with functional classes. Default is "functional_class"
#'
#' @return A list containing the sorted CNA matrix and chromosome annotations
#' @export
#'
#' @examples
#' \dontrun{
#' cna_data <- extractCNAData(seurat_obj)
#' }
extractCNAData <- function(seurat_obj,
                          assay_name = "cna",
                          functional_class_col = "functional_class") {
  
  # Check if CNA assay exists
  if (!assay_name %in% names(seurat_obj@assays)) {
    stop("CNA assay not found in Seurat object")
  }
  
  # Check if functional class column exists
  if (!functional_class_col %in% colnames(seurat_obj@meta.data)) {
    stop("Functional class column not found in Seurat metadata")
  }
  
  # Get CNA matrix
  cna_matrix <- Seurat::GetAssayData(seurat_obj, assay = assay_name, slot = "data")
  
  # Get chromosome information
  gene_locations <- rownames(cna_matrix)
  
  # Parse chromosome data
  # Assuming format is like "chr1_123456" or similar
  chrom_data <- data.frame(
    gene = gene_locations,
    chrom = gsub("_.*$", "", gene_locations),
    chrompos = as.numeric(gsub("^.*_", "", gene_locations)),
    stringsAsFactors = FALSE
  )
  
  # Clean chromosome names to keep just the number/letter
  chrom_data$chrom <- gsub("chr", "", chrom_data$chrom)
  
  # Calculate absolute position (for visualization)
  chrom_data$abspos <- 1:nrow(chrom_data)
  
  # Create chromosome arm annotations
  chrom_anno <- data.frame(
    gene = chrom_data$gene,
    chrom = as.character(chrom_data$chrom),
    arm = rep("p", nrow(chrom_data)),  # Default to p
    stringsAsFactors = FALSE
  )
  
  # Define chromosome lengths and centromere positions (hg38)
  chrom_lengths <- data.frame(
    chrom = c(1:22, "X", "Y"),
    size = c(248956422, 242193529, 198295559, 190214555, 181538259, 
             170805979, 159345973, 145138636, 138394717, 133797422, 
             135086622, 133275309, 114364328, 107043718, 101991189, 
             90338345, 83257441, 80373285, 58617616, 64444167, 
             46709983, 50818468, 156040895, 57227415),
    centromere = c(125000000, 93300000, 91000000, 50400000, 48400000,
                   61000000, 59900000, 45600000, 49000000, 40200000,
                   53700000, 35800000, 17900000, 17600000, 19000000,
                   36600000, 24000000, 17200000, 26500000, 27500000,
                   13200000, 14700000, 60600000, 10400000),
    stringsAsFactors = FALSE
  )
  
  # Assign chromosome arms based on position
  for (i in 1:nrow(chrom_anno)) {
    chr <- chrom_anno$chrom[i]
    pos <- chrom_data$chrompos[i]
    
    # Find the matching chromosome in reference data
    idx <- which(chrom_lengths$chrom == chr)
    if (length(idx) > 0) {
      centro_pos <- chrom_lengths$centromere[idx]
      chrom_anno$arm[i] <- ifelse(pos < centro_pos, "p", "q")
    }
  }
  
  # Add combined chrom.arm annotation
  chrom_anno$chrom_arm <- paste0(chrom_anno$chrom, chrom_anno$arm)
  
  # Sort the matrix by chromosomes and positions
  sorted_idx <- order(
    as.numeric(factor(chrom_data$chrom, levels = c(1:22, "X", "Y"))), 
    chrom_data$chrompos
  )
  sorted_cna <- cna_matrix[sorted_idx, ]
  sorted_anno <- chrom_anno[sorted_idx, ]
  
  # Create metadata dataframe with functional classification
  cell_meta <- data.frame(
    cell_id = colnames(cna_matrix),
    functional_class = seurat_obj@meta.data[[functional_class_col]][match(colnames(cna_matrix), colnames(seurat_obj))],
    stringsAsFactors = FALSE
  )
  
  # Replace NAs with "Unknown"
  cell_meta$functional_class[is.na(cell_meta$functional_class)] <- "Unknown"
  
  return(list(
    sorted_cna = sorted_cna,
    sorted_anno = sorted_anno,
    cell_meta = cell_meta,
    chrom_lengths = chrom_lengths
  ))
}

#' Calculate CNA Profiles by Functional Class
#'
#' Calculates average CNA profiles for each functional class.
#'
#' @param cna_data Output from extractCNAData()
#' @param classes Vector of class names to calculate profiles for. Default is 
#'        c("Hypoxic", "Proliferative", "Undefined")
#'
#' @return Data frame with CNA profiles by functional class and chromosome arm
#' @export
#'
#' @examples
#' \dontrun{
#' cna_data <- extractCNAData(seurat_obj)
#' profiles <- calculateCNAProfiles(cna_data)
#' }
calculateCNAProfiles <- function(cna_data,
                               classes = c("Hypoxic", "Proliferative", "Undefined")) {
  
  # Extract data from the list
  sorted_cna <- cna_data$sorted_cna
  sorted_anno <- cna_data$sorted_anno
  cell_meta <- cna_data$cell_meta
  
  # Filter out cells with unknown class
  cell_meta_known <- cell_meta[cell_meta$functional_class != "Unknown", ]
  sorted_cna_known <- sorted_cna[, colnames(sorted_cna) %in% cell_meta_known$cell_id, drop = FALSE]
  
  # Initialize the result data frame
  result <- sorted_anno
  
  # Calculate mean profiles for each functional class
  for (class in classes) {
    # Identify cells in this class
    class_cells <- cell_meta_known$cell_id[cell_meta_known$functional_class == class]
    
    if (length(class_cells) > 0) {
      # Calculate mean profile
      class_profile <- rowMeans(sorted_cna_known[, colnames(sorted_cna_known) %in% class_cells, drop = FALSE], na.rm = TRUE)
      
      # Add to result
      result[[tolower(class)]] <- class_profile
    }
  }
  
  # Calculate differences between classes
  if (all(c("hypoxic", "proliferative") %in% colnames(result))) {
    result$diff_hypoxic_proliferative <- result$hypoxic - result$proliferative
  }
  
  if (all(c("hypoxic", "undefined") %in% colnames(result))) {
    result$diff_hypoxic_undefined <- result$hypoxic - result$undefined
  }
  
  if (all(c("proliferative", "undefined") %in% colnames(result))) {
    result$diff_proliferative_undefined <- result$proliferative - result$undefined
  }
  
  return(result)
}

#' Aggregate CNA Profiles by Chromosome Arm
#'
#' Aggregates gene-level CNA profiles to chromosome arm level.
#'
#' @param cna_profiles Output from calculateCNAProfiles()
#' @param value_columns Columns containing CNA values to aggregate. If NULL (default),
#'        aggregates all columns except gene, chrom, arm, and chrom_arm.
#'
#' @return Data frame with CNA profiles by chromosome arm
#' @export
#'
#' @examples
#' \dontrun{
#' cna_data <- extractCNAData(seurat_obj)
#' profiles <- calculateCNAProfiles(cna_data)
#' arm_summary <- aggregateByChromosomeArm(profiles)
#' }
aggregateByChromosomeArm <- function(cna_profiles, value_columns = NULL) {
  
  # Identify columns to aggregate if not specified
  if (is.null(value_columns)) {
    value_columns <- setdiff(
      colnames(cna_profiles),
      c("gene", "chrom", "arm", "chrom_arm")
    )
  }
  
  # Check that value columns exist
  missing_cols <- setdiff(value_columns, colnames(cna_profiles))
  if (length(missing_cols) > 0) {
    stop("The following columns are missing from cna_profiles: ", paste(missing_cols, collapse = ", "))
  }
  
  # Create formula for aggregation
  formula_str <- paste("chrom_arm ~", paste(value_columns, collapse = " + "))
  
  # Aggregate by chromosome arm
  arm_summary <- aggregate(
    formula = as.formula(formula_str),
    data = cna_profiles,
    FUN = mean,
    na.rm = TRUE
  )
  
  # Extract chromosome and arm for sorting
  arm_summary$chrom <- sub("([0-9XY]+)[pq]", "\\1", arm_summary$chrom_arm)
  arm_summary$arm <- sub("[0-9XY]+([pq])", "\\1", arm_summary$chrom_arm)
  
  # Sort by chromosome and arm
  arm_summary <- arm_summary[order(
    as.numeric(factor(arm_summary$chrom, levels = c(1:22, "X", "Y"))),
    arm_summary$arm
  ), ]
  
  return(arm_summary)
}

#' Plot CNA Heatmap by Functional Class
#'
#' Creates a heatmap visualization of CNA values by chromosome arm and functional class.
#'
#' @param arm_summary Output from aggregateByChromosomeArm()
#' @param class_columns Names of columns for each functional class
#' @param min_value Minimum value for color scale. If NULL, will use data minimum.
#' @param max_value Maximum value for color scale. If NULL, will use data maximum.
#' @param color_palette Vector of colors for the heatmap. Default is colorRampPalette(c("blue", "white", "red"))(100)
#' @param cluster_rows Whether to cluster rows. Default is FALSE
#' @param title Plot title. Default is "CNA Profiles by Chromosome Arm"
#' @param annotation_row Optional data frame for row annotation
#'
#' @return A pheatmap object
#' @export
#'
#' @examples
#' \dontrun{
#' cna_data <- extractCNAData(seurat_obj)
#' profiles <- calculateCNAProfiles(cna_data)
#' arm_summary <- aggregateByChromosomeArm(profiles)
#' plotCNAHeatmap(arm_summary, class_columns = c("hypoxic", "proliferative", "undefined"))
#' }
plotCNAHeatmap <- function(arm_summary,
                          class_columns = c("hypoxic", "proliferative", "undefined"),
                          min_value = NULL,
                          max_value = NULL,
                          color_palette = NULL,
                          cluster_rows = FALSE,
                          title = "CNA Profiles by Chromosome Arm",
                          annotation_row = NULL) {
  
  # Check if required packages are available
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package pheatmap must be installed")
  }
  
  # Extract columns for heatmap
  heatmap_data <- as.matrix(arm_summary[, class_columns, drop = FALSE])
  rownames(heatmap_data) <- arm_summary$chrom_arm
  
  # Convert column names to title case for display
  colnames(heatmap_data) <- sapply(class_columns, function(x) {
    paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  })
  
  # Set color scale range
  if (is.null(min_value)) {
    min_value <- min(heatmap_data, na.rm = TRUE)
  }
  if (is.null(max_value)) {
    max_value <- max(heatmap_data, na.rm = TRUE)
  }
  
  # Use symmetric range if appropriate
  if (min_value < 0 && max_value > 0) {
    abs_max <- max(abs(min_value), abs(max_value))
    min_value <- -abs_max
    max_value <- abs_max
  }
  
  # Create color palette if not provided
  if (is.null(color_palette)) {
    color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  }
  
  # Set breaks
  breaks <- seq(min_value, max_value, length.out = length(color_palette) + 1)
  
  # Create the heatmap
  heatmap_plot <- pheatmap::pheatmap(
    heatmap_data,
    cluster_rows = cluster_rows,
    cluster_cols = FALSE,
    annotation_row = annotation_row,
    color = color_palette,
    breaks = breaks,
    cellwidth = 50,
    cellheight = 12,
    main = title,
    silent = TRUE  # Return the object instead of plotting directly
  )
  
  return(heatmap_plot)
}

#' Plot CNA Chromosome Ideogram
#'
#' Creates a chromosome ideogram visualization showing CNA values by functional class.
#'
#' @param cna_profiles Output from calculateCNAProfiles()
#' @param classes Vector of class names to visualize. Default is 
#'        c("Hypoxic", "Proliferative", "Undefined")
#' @param min_value Minimum value for color scale. Default is -0.07
#' @param max_value Maximum value for color scale. Default is 0.09
#' @param colors Vector of colors for each class. Default is 
#'        c("#E41A1C", "#377EB8", "#9E9E9E") (red, blue, gray)
#' @param output_file Path to save the PDF. If NULL, will display the plot.
#' @param width PDF width in inches. Default is 15
#' @param height PDF height in inches. Default is 12
#'
#' @return NULL (creates a plot)
#' @export
#'
#' @examples
#' \dontrun{
#' cna_data <- extractCNAData(seurat_obj)
#' profiles <- calculateCNAProfiles(cna_data)
#' plotChromosomeIdeogram(profiles, output_file = "chromosome_ideogram.pdf")
#' }
plotChromosomeIdeogram <- function(cna_profiles,
                                  classes = c("Hypoxic", "Proliferative", "Undefined"),
                                  min_value = -0.07,
                                  max_value = 0.09,
                                  colors = c("#E41A1C", "#377EB8", "#9E9E9E"),
                                  output_file = NULL,
                                  width = 15,
                                  height = 12) {
  
  # Set up color palettes
  heatmap_colors <- colorRampPalette(c("darkblue", "blue", "lightblue", 
                                      "white", "pink", "red", "darkred"))(100)
  
  # Open PDF if output file is specified
  if (!is.null(output_file)) {
    pdf(output_file, width = width, height = height)
    on.exit(dev.off())
  }
  
  # Create a layout with 24 chromosomes (5x5 grid) plus legend
  layout_matrix <- matrix(c(1:24, rep(25, 1)), nrow = 5, byrow = TRUE)
  layout(layout_matrix)
  
  # Convert class names to lowercase for accessing data
  class_cols <- tolower(classes)
  
  # Split chromosome annotation by chromosome
  chrom_list <- split(cna_profiles, cna_profiles$chrom)
  
  # Function to plot each chromosome with multiple categories
  plot_chromosome_multi_categories <- function(chr_data, title) {
    # Split by arm
    p_arm <- chr_data[chr_data$arm == "p", ]
    q_arm <- chr_data[chr_data$arm == "q", ]
    
    # Set up the plot
    par(mar = c(4, 4, 3, 1))
    
    # Create empty plot
    plot(0, 0, type = "n", xlim = c(-1.5, 1.5), ylim = c(0, 100),
         xlab = "", ylab = "", main = title, axes = FALSE)
    
    # Plot rectangles for the categories
    offset_vals <- seq(-1, 1, length.out = length(classes))
    cat_cols <- colors
    
    # Add chromosome outlines first
    for(cat_idx in 1:length(classes)) {
      rect(offset_vals[cat_idx] - 0.3, 0, offset_vals[cat_idx] + 0.3, 100, 
           border = cat_cols[cat_idx], lwd = 2)
    }
    
    # Add centromeres
    for(cat_idx in 1:length(classes)) {
      rect(offset_vals[cat_idx] - 0.3, 49, offset_vals[cat_idx] + 0.3, 51, 
           col = "gray50", border = "black")
    }
    
    # Plot the p arm (top half)
    if(nrow(p_arm) > 0) {
      # Calculate band heights
      p_heights <- rep(49/nrow(p_arm), nrow(p_arm))
      
      # For each category, plot rectangles
      for(cat_idx in 1:length(classes)) {
        val_col <- class_cols[cat_idx]
        
        for(i in 1:nrow(p_arm)) {
          if(!is.na(p_arm[[val_col]][i])) {
            # Scale the value to 1-100 for color index
            val <- p_arm[[val_col]][i]
            color_idx <- round(((val - min_value) / (max_value - min_value)) * 100)
            color_idx <- max(1, min(100, color_idx))
            
            # Simplified band calculation
            band_bottom <- 51 + ((i-1) * (49/nrow(p_arm)))
            band_top <- 51 + (i * (49/nrow(p_arm)))
            rect(offset_vals[cat_idx] - 0.3, band_bottom, 
                 offset_vals[cat_idx] + 0.3, band_top,
                 col = heatmap_colors[color_idx],
                 border = NA)
          }
        }
      }
    }
    
    # Plot the q arm (bottom half)
    if(nrow(q_arm) > 0) {
      # Calculate band heights
      q_heights <- rep(49/nrow(q_arm), nrow(q_arm))
      
      # For each category, plot rectangles
      for(cat_idx in 1:length(classes)) {
        val_col <- class_cols[cat_idx]
        
        for(i in 1:nrow(q_arm)) {
          if(!is.na(q_arm[[val_col]][i])) {
            # Scale the value to 1-100 for color index
            val <- q_arm[[val_col]][i]
            color_idx <- round(((val - min_value) / (max_value - min_value)) * 100)
            color_idx <- max(1, min(100, color_idx))
            
            # Simplified band calculation
            band_top <- 49 - ((i-1) * (49/nrow(q_arm)))
            band_bottom <- 49 - (i * (49/nrow(q_arm)))
            rect(offset_vals[cat_idx] - 0.3, band_bottom, 
                 offset_vals[cat_idx] + 0.3, band_top,
                 col = heatmap_colors[color_idx],
                 border = NA)
          }
        }
      }
    }
    
    # Add labels
    for (i in 1:length(classes)) {
      text(offset_vals[i], -10, classes[i], col = cat_cols[i], cex = 0.7)
    }
  }
  
  # Plot each chromosome
  for(chr in names(chrom_list)) {
    plot_chromosome_multi_categories(chrom_list[[chr]], paste("Chr", chr))
  }
  
  # Add a comprehensive legend
  par(mar = c(4, 4, 4, 4))
  plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(0, 10), axes = FALSE, xlab = "", ylab = "")
  title("CNA Value Legend", font = 2)
  
  # Add color gradient with clearer boundaries
  gradient_width <- 6
  gradient_x <- seq(2, 2 + gradient_width, length.out = 101)
  for(i in 1:100) {
    rect(gradient_x[i], 7, gradient_x[i+1], 8, 
         col = heatmap_colors[i], 
         border = NA)
  }
  
  # Add border to the gradient
  rect(2, 7, 2 + gradient_width, 8, border = "black", lwd = 1)
  
  # Add labels with better positioning
  text(2, 8.2, "Loss", pos = 4, cex = 0.9)
  text(8, 8.2, "Gain", pos = 4, cex = 0.9)
  
  # Add tick marks for reference values
  tick_positions <- c(0, 0.25, 0.5, 0.75, 1)
  tick_x <- 2 + gradient_width * tick_positions
  segments(tick_x, 6.8, tick_x, 7, lwd = 1)
  tick_labels <- c(round(min_value, 2), round(min_value/2, 2), 0, 
                   round(max_value/2, 2), round(max_value, 2))
  text(tick_x, 6.6, labels = tick_labels, cex = 0.8)
  
  # Add category legend
  legend_y <- 5
  for (i in 1:length(classes)) {
    rect(3, legend_y - 0.3, 3.3, legend_y + 0.3, col = "white", 
         border = colors[i], lwd = 2)
    text(3.5, legend_y, paste0(classes[i], " (", colors[i], " outline)"), 
         pos = 4, cex = 0.8)
    legend_y <- legend_y - 0.7
  }
  
  # Add title
  text(5, 3.2, "Functional Categories", cex = 0.9, font = 2)
}

#' Create Summary Excel File of CNA Results
#'
#' Creates a comprehensive Excel file with multiple sheets summarizing CNA findings.
#'
#' @param arm_summary Output from aggregateByChromosomeArm()
#' @param class_columns Names of columns for each functional class. Default is 
#'        c("hypoxic", "proliferative", "undefined")
#' @param output_file Path to the output Excel file. Default is "CNA_Analysis_Findings.xlsx"
#' @param significance_threshold Threshold for considering a CNA significant. Default is 0.02
#'
#' @return NULL (creates a file)
#' @export
#'
#' @examples
#' \dontrun{
#' cna_data <- extractCNAData(seurat_obj)
#' profiles <- calculateCNAProfiles(cna_data)
#' arm_summary <- aggregateByChromosomeArm(profiles)
#' createCNASummaryExcel(arm_summary, output_file = "results/CNA_summary.xlsx")
#' }
createCNASummaryExcel <- function(arm_summary,
                                 class_columns = c("hypoxic", "proliferative", "undefined"),
                                 output_file = "CNA_Analysis_Findings.xlsx",
                                 significance_threshold = 0.02) {
  
  # Check if required packages are available
  if (!requireNamespace("writexl", quietly = TRUE)) {
    stop("Package writexl must be installed")
  }
  
  # Create title case function
  title_case <- function(x) {
    paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  }
  
  # Ensure output directory exists
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir) && output_dir != ".") {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # ===== SHEET 1-3: TOP CNA REGIONS BY CATEGORY =====
  top_tables <- list()
  
  for (class in class_columns) {
    # Get top 10 altered regions for this category
    top_altered <- arm_summary[order(abs(arm_summary[[class]]), decreasing = TRUE), ][1:10, ]
    
    # Prepare top regions for Excel
    top_table <- data.frame(
      Rank = 1:10,
      Chromosome_Arm = top_altered$chrom_arm,
      CNA_Value = round(top_altered[[class]], 4),
      Status = ifelse(top_altered[[class]] > 0, "Gain", "Loss"),
      Chromosome = top_altered$chrom,
      Arm = top_altered$arm
    )
    
    top_tables[[paste0("Top_", title_case(class), "_CNAs")]] <- top_table
  }
  
  # ===== SHEET 4: ALL CHROMOSOME ARMS WITH CNA VALUES =====
  # Prepare complete list of chromosome arms with their CNA values
  all_arms_table <- data.frame(
    Chromosome_Arm = arm_summary$chrom_arm,
    Chromosome = arm_summary$chrom,
    Arm = arm_summary$arm
  )
  
  # Add columns for each class
  for (class in class_columns) {
    all_arms_table[[paste0(title_case(class), "_CNA")]] <- round(arm_summary[[class]], 4)
    all_arms_table[[paste0(title_case(class), "_Status")]] <- ifelse(
      arm_summary[[class]] > significance_threshold, "Gain", 
      ifelse(arm_summary[[class]] < -significance_threshold, "Loss", "Neutral")
    )
  }
  
  # Sort by chromosome and arm
  all_arms_table <- all_arms_table[order(all_arms_table$Chromosome, all_arms_table$Arm), ]
  
  # ===== SHEET 5-6: SIGNIFICANT GAINS & LOSSES =====
  all_gains <- data.frame()
  all_losses <- data.frame()
  
  for (class in class_columns) {
    # Filter for significant gains (CNA > threshold)
    gains <- all_arms_table[all_arms_table[[paste0(title_case(class), "_CNA")]] > significance_threshold, 
                           c("Chromosome_Arm", paste0(title_case(class), "_CNA"))]
    colnames(gains)[2] <- "CNA_Value"
    gains$Category <- title_case(class)
    
    # Filter for significant losses (CNA < -threshold)
    losses <- all_arms_table[all_arms_table[[paste0(title_case(class), "_CNA")]] < -significance_threshold, 
                            c("Chromosome_Arm", paste0(title_case(class), "_CNA"))]
    colnames(losses)[2] <- "CNA_Value"
    losses$Category <- title_case(class)
    
    # Add to combined tables
    all_gains <- rbind(all_gains, gains)
    all_losses <- rbind(all_losses, losses)
  }
  
  # Sort the tables
  all_gains <- all_gains[order(all_gains$CNA_Value, decreasing = TRUE), ]
  all_losses <- all_losses[order(all_losses$CNA_Value), ]  # Ascending for losses (most negative first)
  
  # ===== SHEET 7: KEY CANCER GENE REGIONS =====
  # Define key cancer genes by chromosome arm
  key_cancer_regions <- data.frame(
    Chromosome_Arm = c("1p", "1q", "2p", "3p", "3q", "4p", "6p", "7p", "8q", "9p", "10q", "11q", "12p", "16q", "17p", "17q", "19q", "20q", "21q"),
    Key_Genes = c(
      "SDHB, ARID1A, CDKN2C",               # 1p
      "MCL1, MDM4, ABL2",                   # 1q
      "MYCN, REL, EPAS1 (HIF2A)",           # 2p
      "VHL, SETD2, BAP1",                   # 3p
      "PIK3CA, SOX2, TP63",                 # 3q
      "FGFR3, TACC3",                       # 4p
      "CDKN1A, PIM1, VEGFA",                # 6p
      "EGFR, TWIST1, HOXA9",                # 7p
      "MYC, RECQL4, RAD21",                 # 8q
      "CDKN2A, CDKN2B, JAK2",               # 9p
      "PTEN, FAS, FGFR2",                   # 10q
      "ATM, BIRC3, CCND1",                  # 11q
      "KRAS, CCND2, SOX5",                  # 12p
      "CDH1, CTCF",                         # 16q
      "TP53, MAP2K4",                       # 17p
      "BRCA1, ERBB2, STAT3",                # 17q
      "AKT2, PPP2R1A, CCNE1",               # 19q
      "AURKA, SRC, TOP1",                   # 20q
      "RUNX1, ERG"                          # 21q
    ),
    Potential_Role = c(
      "Tumor suppressors",
      "Oncogenes",
      "Oncogenes",
      "Tumor suppressors",
      "Oncogenes",
      "Growth factors",
      "Cell cycle regulation",
      "Growth factors",
      "Oncogenes",
      "Tumor suppressors",
      "Tumor suppressors",
      "DNA damage & cell cycle",
      "Oncogenes",
      "Cell adhesion",
      "Tumor suppressors",
      "DNA repair & growth",
      "Cell cycle & growth",
      "Mitosis & growth",
      "Transcription factors"
    )
  )
  
  # Merge with CNA data
  key_regions_cna <- merge(key_cancer_regions, all_arms_table, by = "Chromosome_Arm")
  
  # Sort by average absolute CNA value
  avg_abs_cna <- rowMeans(abs(as.matrix(key_regions_cna[, grep("_CNA$", colnames(key_regions_cna)), drop = FALSE])))
  key_regions_cna$Average_Abs_CNA <- avg_abs_cna
  key_regions_cna <- key_regions_cna[order(key_regions_cna$Average_Abs_CNA, decreasing = TRUE), ]
  
  # ===== SHEET 8: DIFFERENCES BETWEEN CATEGORIES =====
  # Calculate absolute differences between categories
  category_differences <- data.frame(
    Chromosome_Arm = arm_summary$chrom_arm
  )
  
  # Add all pairwise differences
  if (length(class_columns) > 1) {
    for (i in 1:(length(class_columns)-1)) {
      for (j in (i+1):length(class_columns)) {
        class1 <- class_columns[i]
        class2 <- class_columns[j]
        diff_col <- paste0(class1, "_vs_", class2)
        
        if (diff_col %in% colnames(arm_summary)) {
          category_differences[[paste0(title_case(class1), "_vs_", title_case(class2))]] <- 
            round(arm_summary[[diff_col]], 4)
        } else {
          category_differences[[paste0(title_case(class1), "_vs_", title_case(class2))]] <- 
            round(arm_summary[[class1]] - arm_summary[[class2]], 4)
        }
      }
    }
  }
  
  # Add significance indicators
  for (col in setdiff(colnames(category_differences), "Chromosome_Arm")) {
    sig_col <- paste0(col, "_Significant")
    category_differences[[sig_col]] <- abs(category_differences[[col]]) > significance_threshold
  }
  
  # Calculate maximum difference for sorting
  diff_cols <- grep("_vs_", colnames(category_differences), value = TRUE)
  diff_cols <- diff_cols[!grepl("_Significant$", diff_cols)]
  
   if (length(diff_cols) > 0) {
    category_differences$Max_Difference <- apply(
      abs(category_differences[, diff_cols, drop = FALSE]), 
      1, max
    )
    category_differences <- category_differences[order(category_differences$Max_Difference, decreasing = TRUE), ]
    
    # Clean up by removing temporary column
    category_differences$Max_Difference <- NULL
  }
  
  # ===== COMBINE ALL SHEETS INTO EXCEL FILE =====
  sheets_list <- c(
    top_tables,
    list(
      "All_Chromosome_Arms" = all_arms_table,
      "Significant_Gains" = all_gains,
      "Significant_Losses" = all_losses,
      "Key_Cancer_Gene_Regions" = key_regions_cna
    )
  )
  
  # Add differences sheet if it has content
  if (ncol(category_differences) > 1) {
    sheets_list$Category_Differences <- category_differences
  }
  
  # Write to Excel
  writexl::write_xlsx(sheets_list, path = output_file)
  
  message(paste("Created comprehensive CNA analysis summary Excel file at:", output_file))
}


