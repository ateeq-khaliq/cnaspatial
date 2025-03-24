#' Plot Functional Score Scatter Plot
#'
#' Creates a scatter plot of hypoxia vs proliferation scores with cells colored by category.
#'
#' @param seurat_obj A Seurat object with hypoxia and proliferation scores
#' @param hypoxia_score_col Name of column containing hypoxia scores. Default is "HypoxiaScore1"
#' @param prolif_score_col Name of column containing proliferation scores. Default is "ProliferationScore1"
#' @param category_col Column to use for coloring points. Default is "HypoxProlifCategory"
#' @param point_size Size of points in the plot. Default is 1
#' @param title Plot title. Default is "Hypoxia vs Proliferation Scores"
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' score_plot <- plotFunctionalScores(seurat_obj)
#' print(score_plot)
#' }
plotFunctionalScores <- function(seurat_obj,
                               hypoxia_score_col = "HypoxiaScore1",
                               prolif_score_col = "ProliferationScore1",
                               category_col = "HypoxProlifCategory",
                               point_size = 1,
                               title = "Hypoxia vs Proliferation Scores") {
  
  # Check if required packages are installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 must be installed")
  }
  
  # Check if columns exist
  if (!all(c(hypoxia_score_col, prolif_score_col, category_col) %in% colnames(seurat_obj@meta.data))) {
    stop("One or more required columns not found in Seurat metadata")
  }
  
  # Extract data for plotting
  plot_data <- data.frame(
    Hypoxia = seurat_obj@meta.data[[hypoxia_score_col]],
    Proliferation = seurat_obj@meta.data[[prolif_score_col]],
    Category = seurat_obj@meta.data[[category_col]]
  )
  
  # Define category colors based on the categories present
  categories <- unique(plot_data$Category)
  
  if (all(c("High_Both", "High_Hypoxia", "High_Proliferation", "Medium", "Undefined") %in% categories)) {
    # Standard 5-category color scheme
    colors <- c(
      "High_Both" = "#8C2981", # Purple
      "High_Hypoxia" = "#D6604D", # Red-orange
      "High_Proliferation" = "#2166AC", # Blue
      "Medium" = "#4D9221", # Green
      "Undefined" = "#999999" # Gray
    )
  } else if (all(c("High_Hypoxia_Dominant", "High_Proliferation_Dominant", "High_Hypoxia", 
                 "High_Proliferation", "Medium", "Undefined") %in% categories)) {
    # Refined 6-category color scheme
    colors <- c(
      "High_Hypoxia_Dominant" = "#D73027", # Bright red
      "High_Proliferation_Dominant" = "#4575B4", # Bright blue
      "High_Hypoxia" = "#FDAE61", # Orange
      "High_Proliferation" = "#ABD9E9", # Light blue
      "Medium" = "#4D9221", # Green
      "Undefined" = "#999999" # Gray
    )
  } else if (all(c("Hypoxic", "Proliferative", "Undefined") %in% categories)) {
    # Simplified 3-category color scheme
    colors <- c(
      "Hypoxic" = "#E41A1C", # Red
      "Proliferative" = "#377EB8", # Blue
      "Undefined" = "#999999" # Gray
    )
  } else {
    # Default to rainbow colors if categories don't match expected sets
    colors <- setNames(
      rainbow(length(categories)),
      categories
    )
  }
  
  # Create scatter plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Hypoxia, y = Proliferation, color = Category)) +
    ggplot2::geom_point(size = point_size, alpha = 0.7) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(
      title = title,
      x = "Hypoxia Score",
      y = "Proliferation Score",
      color = "Category"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  return(p)
}

#' Plot CNA Barplot by Chromosome Arm
#'
#' Creates a barplot of CNA values by chromosome arm for different functional categories.
#'
#' @param arm_summary Data frame from aggregateByChromosomeArm()
#' @param class_columns Names of columns for each functional class. Default is 
#'        c("hypoxic", "proliferative", "undefined")
#' @param n_top Number of top arms to display. Default is 15
#' @param sort_by Class to sort arms by. Default is the first class in class_columns
#' @param sort_abs Whether to sort by absolute CNA values. Default is TRUE
#' @param colors Vector of colors for classes. If NULL, will use default colors.
#' @param horizontal Whether to create a horizontal barplot. Default is TRUE
#' @param title Plot title. Default is "Top Chromosome Arms by CNA Value"
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' cna_data <- extractCNAData(seurat_obj)
#' profiles <- calculateCNAProfiles(cna_data)
#' arm_summary <- aggregateByChromosomeArm(profiles)
#' barplot <- plotCNABarplot(arm_summary)
#' print(barplot)
#' }
plotCNABarplot <- function(arm_summary,
                          class_columns = c("hypoxic", "proliferative", "undefined"),
                          n_top = 15,
                          sort_by = NULL,
                          sort_abs = TRUE,
                          colors = NULL,
                          horizontal = TRUE,
                          title = "Top Chromosome Arms by CNA Value") {
  
  # Check if required packages are installed
  if (!requireNamespace("ggplot2", quietly = TRUE) || 
      !requireNamespace("tidyr", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE)) {
    stop("Packages ggplot2, tidyr, and dplyr must be installed")
  }
  
  # Set sort_by to first class column if not specified
  if (is.null(sort_by)) {
    sort_by <- class_columns[1]
  }
  
  # Check if sort_by is valid
  if (!sort_by %in% class_columns) {
    stop("sort_by must be one of the class_columns")
  }
  
  # Define colors if not provided
  if (is.null(colors)) {
    colors <- c(
      "hypoxic" = "#E41A1C",      # Red
      "proliferative" = "#377EB8", # Blue
      "undefined" = "#999999"      # Gray
    )
  }
  
  # Make sure colors are properly named
  if (is.null(names(colors)) || !all(class_columns %in% names(colors))) {
    colors <- setNames(colors[1:length(class_columns)], class_columns)
  }
  
  # Convert class column names to title case for display
  class_display_names <- sapply(class_columns, function(x) {
    paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  })
  names(class_display_names) <- class_columns
  
  # Select chromosome arms
  if (sort_abs) {
    sorted_arms <- arm_summary[order(abs(arm_summary[[sort_by]]), decreasing = TRUE), "chrom_arm"][1:min(n_top, nrow(arm_summary))]
  } else {
    sorted_arms <- arm_summary[order(arm_summary[[sort_by]], decreasing = TRUE), "chrom_arm"][1:min(n_top, nrow(arm_summary))]
  }
  
  # Filter arm_summary to top arms and selected columns
  plot_data <- arm_summary[arm_summary$chrom_arm %in% sorted_arms, c("chrom_arm", class_columns)]
  
  # Reshape data for plotting
  plot_data_long <- tidyr::pivot_longer(
    plot_data,
    cols = class_columns,
    names_to = "Class",
    values_to = "CNA_Value"
  )
  
  # Add display names
  plot_data_long$Class_Display <- class_display_names[plot_data_long$Class]
  
  # Order arms by absolute value of sort_by class
  if (sort_abs) {
    arm_order <- plot_data[order(abs(plot_data[[sort_by]]), decreasing = TRUE), "chrom_arm"]
  } else {
    arm_order <- plot_data[order(plot_data[[sort_by]], decreasing = TRUE), "chrom_arm"]
  }
  plot_data_long$chrom_arm <- factor(plot_data_long$chrom_arm, levels = arm_order)
  
  # Create barplot
  if (horizontal) {
    p <- ggplot2::ggplot(plot_data_long, ggplot2::aes(x = CNA_Value, y = chrom_arm, fill = Class_Display)) +
      ggplot2::geom_col(position = "dodge", width = 0.7, color = "black", linewidth = 0.2) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.5) +
      ggplot2::scale_fill_manual(values = colors, name = "Functional Class") +
      ggplot2::labs(
        title = title,
        x = "Copy Number Alteration (CNA) Value",
        y = "Chromosome Arm"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.title = ggplot2::element_text(face = "bold"),
        axis.text.y = ggplot2::element_text(face = "bold"),
        legend.position = "bottom"
      )
  } else {
    p <- ggplot2::ggplot(plot_data_long, ggplot2::aes(x = chrom_arm, y = CNA_Value, fill = Class_Display)) +
      ggplot2::geom_col(position = "dodge", width = 0.7, color = "black", linewidth = 0.2) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.5) +
      ggplot2::scale_fill_manual(values = colors, name = "Functional Class") +
      ggplot2::labs(
        title = title,
        x = "Chromosome Arm",
        y = "Copy Number Alteration (CNA) Value"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.title = ggplot2::element_text(face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold"),
        legend.position = "bottom"
      )
  }
  
  return(p)
}

#' Plot CNA Circos Plot
#'
#' Creates a circos plot visualizing CNA values by chromosome and functional class.
#'
#' @param arm_summary Data frame from aggregateByChromosomeArm()
#' @param class_columns Names of columns for each functional class. Default is 
#'        c("hypoxic", "proliferative", "undefined")
#' @param min_val Minimum value for color scale. If NULL, uses data minimum.
#' @param max_val Maximum value for color scale. If NULL, uses data maximum.
#' @param colors Vector of colors for classes. If NULL, will use default colors.
#' @param title Plot title. Default is "Chromosome CNA Profiles by Functional Class"
#'
#' @return NULL (creates a plot)
#' @export
#'
#' @examples
#' \dontrun{
#' cna_data <- extractCNAData(seurat_obj)
#' profiles <- calculateCNAProfiles(cna_data)
#' arm_summary <- aggregateByChromosomeArm(profiles)
#' plotCNACircos(arm_summary)
#' }
plotCNACircos <- function(arm_summary,
                         class_columns = c("hypoxic", "proliferative", "undefined"),
                         min_val = NULL,
                         max_val = NULL,
                         colors = NULL,
                         title = "Chromosome CNA Profiles by Functional Class") {
  
  # Check if required packages are installed
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package circlize must be installed")
  }
  
  # Calculate data-driven min/max range for scaling if not provided
  if (is.null(min_val) || is.null(max_val)) {
    actual_min <- min(as.matrix(arm_summary[, class_columns]), na.rm = TRUE)
    actual_max <- max(as.matrix(arm_summary[, class_columns]), na.rm = TRUE)
    
    # Use slightly padded values for better visualization
    if (is.null(min_val)) {
      min_val <- floor(actual_min * 100) / 100  # Round down to nearest 0.01
    }
    if (is.null(max_val)) {
      max_val <- ceiling(actual_max * 100) / 100  # Round up to nearest 0.01
    }
    
    message(paste("Using data-driven color scale range:", min_val, "to", max_val))
  }
  
  # Define colors if not provided
  if (is.null(colors)) {
    colors <- c(
      "hypoxic" = "#E41A1C",      # Red
      "proliferative" = "#377EB8", # Blue
      "undefined" = "#9E9E9E"      # Gray
    )
  }
  
  # Make sure colors are properly named
  if (is.null(names(colors)) || !all(class_columns %in% names(colors))) {
    colors <- setNames(colors[1:length(class_columns)], class_columns)
  }
  
  # Define a color palette for CNA values
  blue_red_palette <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", 
                        "#F7F7F7", "#F4A582