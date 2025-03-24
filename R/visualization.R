```r
#' Plot Functional Score Scatter Plot
#'
#' Creates a scatter plot of any two functional scores with cells colored by category.
#'
#' @param seurat_obj A Seurat object with functional scores
#' @param x_score_col Name of column containing the first score to plot on x-axis
#' @param y_score_col Name of column containing the second score to plot on y-axis
#' @param category_col Column to use for coloring points. Default is "FunctionalCategory"
#' @param point_size Size of points in the plot. Default is 1
#' @param title Plot title. Default is "Functional Score Comparison"
#' @param custom_colors Optional named vector of colors for categories
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Default hypoxia vs proliferation plot
#' score_plot <- plotFunctionalScores(seurat_obj, 
#'                                   x_score_col = "HypoxiaScore1", 
#'                                   y_score_col = "ProliferationScore1")
#' 
#' # Custom functional signatures
#' score_plot <- plotFunctionalScores(seurat_obj, 
#'                                  x_score_col = "EMTScore1", 
#'                                  y_score_col = "StemnessScore1",
#'                                  category_col = "EMTStemCategory")
#' print(score_plot)
#' }
plotFunctionalScores <- function(seurat_obj,
                               x_score_col,
                               y_score_col,
                               category_col = "FunctionalCategory",
                               point_size = 1,
                               title = "Functional Score Comparison",
                               custom_colors = NULL) {
  
  # Check if required packages are installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 must be installed")
  }
  
  # Check if columns exist
  if (!all(c(x_score_col, y_score_col, category_col) %in% colnames(seurat_obj@meta.data))) {
    stop("One or more required columns not found in Seurat metadata")
  }
  
  # Extract column labels for axis titles (remove number suffix if present)
  x_label <- gsub("([A-Za-z]+)\\d*$", "\\1", x_score_col)
  y_label <- gsub("([A-Za-z]+)\\d*$", "\\1", y_score_col)
  
  # Extract data for plotting
  plot_data <- data.frame(
    x_score = seurat_obj@meta.data[[x_score_col]],
    y_score = seurat_obj@meta.data[[y_score_col]],
    Category = seurat_obj@meta.data[[category_col]]
  )
  
  # Use provided custom colors or generate based on categories
  categories <- unique(plot_data$Category)
  
  if (!is.null(custom_colors)) {
    # Check if all categories have colors
    missing_cats <- setdiff(categories, names(custom_colors))
    if (length(missing_cats) > 0) {
      warning("Missing colors for categories: ", paste(missing_cats, collapse = ", "), 
              ". Using default colors for these.")
      # Add default colors for missing categories
      custom_colors[missing_cats] <- rainbow(length(missing_cats))
    }
    colors <- custom_colors
  } else {
    # Auto-generate colors based on number of categories
    n_cats <- length(categories)
    if (n_cats <= 8) {
      # Use predefined color schemes for typical cases
      if (any(grepl("High_", categories)) && any(grepl("Medium", categories))) {
        # Likely a standard classification scheme
        colors <- c(
          "High_Both" = "#8C2981",            # Purple
          "High_First" = "#D6604D",           # Red-orange
          "High_Second" = "#2166AC",          # Blue
          "High_First_Dominant" = "#D73027",  # Bright red
          "High_Second_Dominant" = "#4575B4", # Bright blue
          "Medium" = "#4D9221",               # Green
          "Undefined" = "#999999"             # Gray
        )
        # Rename generically if needed (in case of custom signatures)
        if (!"High_First" %in% categories && any(grepl("High_.*_Dominant", categories))) {
          names(colors)[3] <- "High_First_Dominant"
          names(colors)[4] <- "High_Second_Dominant"
        }
      } else {
        # Use a simple color palette for unknown schemes
        colors <- setNames(
          c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")[1:n_cats],
          categories
        )
      }
    } else {
      # For many categories, use a rainbow palette
      colors <- setNames(
        rainbow(n_cats),
        categories
      )
    }
  }
  
  # Create scatter plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x_score, y = y_score, color = Category)) +
    ggplot2::geom_point(size = point_size, alpha = 0.7) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(
      title = title,
      x = paste(x_label, "Score"),
      y = paste(y_label, "Score"),
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
#' @param class_columns Names of columns for each functional class
#' @param n_top Number of top arms to display. Default is 15
#' @param sort_by Class to sort arms by. Default is the first class in class_columns
#' @param sort_abs Whether to sort by absolute CNA values. Default is TRUE
#' @param custom_colors Optional named vector of colors for classes
#' @param horizontal Whether to create a horizontal barplot. Default is TRUE
#' @param title Plot title. Default is "Top Chromosome Arms by CNA Value"
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Standard hypoxic/proliferative barplot
#' barplot <- plotCNABarplot(arm_summary, 
#'                          class_columns = c("hypoxic", "proliferative", "undefined"))
#'                          
#' # Custom functional signatures
#' barplot <- plotCNABarplot(arm_summary,
#'                         class_columns = c("emt_high", "stemness_high", "basal"),
#'                         custom_colors = c("emt_high" = "#D73027", 
#'                                          "stemness_high" = "#4575B4",
#'                                          "basal" = "#4D9221"))
#' print(barplot)
#' }
plotCNABarplot <- function(arm_summary,
                          class_columns,
                          n_top = 15,
                          sort_by = NULL,
                          sort_abs = TRUE,
                          custom_colors = NULL,
                          horizontal = TRUE,
                          title = "Top Chromosome Arms by CNA Value") {
  
  # Check if required packages are installed
  if (!requireNamespace("ggplot2", quietly = TRUE) || 
      !requireNamespace("tidyr", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE)) {
    stop("Packages ggplot2, tidyr, and dplyr must be installed")
  }
  
  # Check if we have class columns
  if (length(class_columns) == 0) {
    stop("No class columns specified")
  }
  
  # Verify all class columns exist in the data
  missing_cols <- setdiff(class_columns, colnames(arm_summary))
  if (length(missing_cols) > 0) {
    stop("The following class columns are missing from arm_summary: ", 
         paste(missing_cols, collapse = ", "))
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
  if (is.null(custom_colors)) {
    # Standard color palette
    std_colors <- c(
      "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
      "#FF7F00", "#FFFF33", "#A65628", "#F781BF"
    )
    
    # Apply colors to class columns
    colors <- setNames(
      std_colors[1:length(class_columns)],
      class_columns
    )
  } else {
    # Use provided custom colors
    colors <- custom_colors
    
    # Check if all class columns have colors
    missing_cols <- setdiff(class_columns, names(colors))
    if (length(missing_cols) > 0) {
      warning("Missing colors for columns: ", paste(missing_cols, collapse = ", "), 
              ". Using default colors for these.")
      
      # Add default colors for missing classes
      std_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
                     "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
      colors[missing_cols] <- std_colors[1:length(missing_cols)]
    }
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
#' @param class_columns Names of columns for each functional class
#' @param min_val Minimum value for color scale. If NULL, uses data minimum.
#' @param max_val Maximum value for color scale. If NULL, uses data maximum.
#' @param custom_colors Optional named vector of colors for classes
#' @param title Plot title. Default is "Chromosome CNA Profiles by Functional Class"
#'
#' @return NULL (creates a plot)
#' @export
#'
#' @examples
#' \dontrun{
#' # Standard hypoxic/proliferative circos plot
#' plotCNACircos(arm_summary, class_columns = c("hypoxic", "proliferative", "undefined"))
#' 
#' # Custom functional signatures
#' plotCNACircos(arm_summary, 
#'              class_columns = c("emt_high", "stemness_high", "basal"),
#'              custom_colors = c("emt_high" = "#D73027", 
#'                               "stemness_high" = "#4575B4",
#'                               "basal" = "#4D9221"))
#' }
plotCNACircos <- function(arm_summary,
                         class_columns,
                         min_val = NULL,
                         max_val = NULL,
                         custom_colors = NULL,
                         title = "Chromosome CNA Profiles by Functional Class") {
  
  # Check if required packages are installed
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package circlize must be installed")
  }
  
  # Check if we have class columns
  if (length(class_columns) == 0) {
    stop("No class columns specified")
  }
  
  # Verify all class columns exist in the data
  missing_cols <- setdiff(class_columns, colnames(arm_summary))
  if (length(missing_cols) > 0) {
    stop("The following class columns are missing from arm_summary: ", 
         paste(missing_cols, collapse = ", "))
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
  if (is.null(custom_colors)) {
    # Standard color palette
    std_colors <- c(
      "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
      "#FF7F00", "#FFFF33", "#A65628", "#F781BF"
    )
    
    # Apply colors to class columns
    colors <- setNames(
      std_colors[1:length(class_columns)],
      class_columns
    )
  } else {
    # Use provided custom colors
    colors <- custom_colors
    
    # Check if all class columns have colors
    missing_cols <- setdiff(class_columns, names(colors))
    if (length(missing_cols) > 0) {
      warning("Missing colors for columns: ", paste(missing_cols, collapse = ", "), 
              ". Using default colors for these.")
      
      # Add default colors for missing classes
      std_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
                     "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
      colors[missing_cols] <- std_colors[1:length(missing_cols)]
    }
  }
  
  # Define a color palette for CNA values
  blue_red_palette <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", 
                        "#F7F7F7", "#F4A582", "#D6604D", "#B2182B", "#67001F")
  
  # Prepare data for circos plot
  # Create a data frame with chromosome info
  circos_data <- data.frame(
    chr = paste0("chr", arm_summary$chrom),
    start = ifelse(arm_summary$arm == "p", 0, 50),
    end = ifelse(arm_summary$arm == "p", 50, 100)
  )
  
  # Add CNA values for all classes
  for (class in class_columns) {
    circos_data[[class]] <- arm_summary[[class]]
  }
  
  # Initialize circos layout
  circlize::circos.clear()
  circlize::circos.par(
    gap.degree = 5, 
    start.degree = 90,
    cell.padding = c(0, 0, 0, 0)
  )
  
  # Prepare chromosome factors and sizes
  chr_names <- unique(circos_data$chr)
  chr_lengths <- rep(100, length(chr_names))
  names(chr_lengths) <- chr_names
  
  # Initialize circos with chromosome information
  circlize::circos.initialize(
    factors = chr_names,
    xlim = matrix(rep(c(0, 100), length(chr_names)), ncol = 2, byrow = TRUE)
  )
  
  # Add chromosome ideogram
  circlize::circos.track(
    ylim = c(0, 1),
    bg.col = "grey95",
    bg.border = NA,
    track.height = 0.05
  )
  
  # Add chromosome labels
  circlize::circos.track(
    ylim = c(0, 1),
    bg.col = NA,
    bg.border = NA,
    track.height = 0.05,
    panel.fun = function(x, y) {
      chr_name <- gsub("chr", "", circlize::CELL_META$sector.index)
      circlize::circos.text(
        x = circlize::CELL_META$xcenter, 
        y = 0.5, 
        labels = chr_name, 
        cex = 0.8,
        facing = "bending.inside"
      )
    }
  )
  
  # Add a track for each functional class
  color_fun <- circlize::colorRamp2(
    seq(min_val, max_val, length = length(blue_red_palette)),
    blue_red_palette
  )
  
  # Convert class column names to title case for display
  class_display_names <- sapply(class_columns, function(x) {
    paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  })
  names(class_display_names) <- class_columns
  
  for (i in seq_along(class_columns)) {
    class <- class_columns[i]
    display_name <- class_display_names[class]
    
    # Create track for this class
    circlize::circos.track(
      factors = chr_names,
      ylim = c(0, 1),
      bg.border = NA,
      track.height = 0.1,
      panel.fun = function(x, y) {
        # Get data for this chromosome
        chr <- circlize::CELL_META$sector.index
        chr_data <- circos_data[circos_data$chr == chr, ]
        
        # Iterate through each arm (if present) in this chromosome 
        if(nrow(chr_data) > 0) {
          for(j in 1:nrow(chr_data)) {
            # Get CNA value and calculate color
            cna_value <- chr_data[[class]][j]
            fill_col <- color_fun(cna_value)
            
            # Draw a rectangle for this segment
            circlize::circos.rect(
              chr_data$start[j], 0,
              chr_data$end[j], 1,
              col = fill_col,
              border = colors[class],
              lwd = 1
            )
          }
        }
      }
    )
    
    # Add class label
    circlize::circos.text(
      sector.index = chr_names[1],
      x = -20,  # Negative value puts it outside the circle
      y = 0.5,  # Center vertically
      labels = display_name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0.5, 0),
      cex = 0.8,
      col = colors[class],
      track.index = i + 2  # +2 because we already have two tracks (ideogram + labels)
    )
  }
  
  # Add title
  title(title, cex.main = 1.2)
  
  # Add color legend for CNA values
  lgd <- circlize::color.legend(
    x = "bottom",
    at = seq(min_val, max_val, length.out = 9),
    col = blue_red_palette,
    title = "CNA Value",
    grid.line.lwd = 0,
    title.adj = c(0.5, 0)
  )
  
  # Return to default par settings
  circlize::circos.clear()
  
  # Return NULL as indicated in function definition
  return(NULL)
}

#' Plot CNA Chromosome Ideogram
#'
#' Creates a chromosome ideogram visualization showing CNA values by functional class.
#'
#' @param cna_profiles Output from calculateCNAProfiles()
#' @param classes Vector of class names to visualize
#' @param min_value Minimum value for color scale. Default is -0.07
#' @param max_value Maximum value for color scale. Default is 0.09
#' @param custom_colors Optional named vector of colors for classes
#' @param output_file Path to save the PDF. If NULL, will display the plot.
#' @param width PDF width in inches. Default is 15
#' @param height PDF height in inches. Default is 12
#'
#' @return NULL (creates a plot)
#' @export
#'
#' @examples
#' \dontrun{
#' # Standard hypoxic/proliferative/undefined ideogram
#' cna_data <- extractCNAData(seurat_obj)
#' profiles <- calculateCNAProfiles(cna_data)
#' plotChromosomeIdeogram(profiles, classes = c("Hypoxic", "Proliferative", "Undefined"))
#' 
#' # Custom functional signatures
#' plotChromosomeIdeogram(profiles, 
#'                       classes = c("EMT", "Stemness", "Basal"),
#'                       custom_colors = c("EMT" = "#D73027", 
#'                                         "Stemness" = "#4575B4",
#'                                         "Basal" = "#4D9221"))
#' }
plotChromosomeIdeogram <- function(cna_profiles,
                                  classes,
                                  min_value = -0.07,
                                  max_value = 0.09,
                                  custom_colors = NULL,
                                  output_file = NULL,
                                  width = 15,
                                  height = 12) {
  
  # Check if we have classes
  if (length(classes) == 0) {
    stop("No classes specified")
  }
  
  # Convert class names to lowercase for accessing data
  class_cols <- tolower(classes)
  
  # Verify all class columns exist in the data
  missing_cols <- setdiff(class_cols, colnames(cna_profiles))
  if (length(missing_cols) > 0) {
    stop("The following class columns are missing from cna_profiles: ", 
         paste(missing_cols, collapse = ", "))
  }
  
  # Set up color palettes
  heatmap_colors <- colorRampPalette(c("darkblue", "blue", "lightblue", 
                                      "white", "pink", "red", "darkred"))(100)
  
  # Define colors if not provided
  if (is.null(custom_colors)) {
    # Standard color palette
    if (length(classes) <= 8) {
      colors <- c(
        "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
        "#FF7F00", "#FFFF33", "#A65628", "#F781BF"
      )[1:length(classes)]
    } else {
      colors <- rainbow(length(classes))
    }
  } else {
    # Check if custom colors match classes
    if (length(custom_colors) < length(classes)) {
      warning("Not enough custom colors provided. Using default colors for remaining classes.")
      std_colors <- c(
        "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
        "#FF7F00", "#FFFF33", "#A65628", "#F781BF"
      )
      # Add missing colors
      missing_classes <- setdiff(classes, names(custom_colors))
      custom_colors[missing_classes] <- std_colors[1:length(missing_classes)]
    }
    colors <- custom_colors[classes]
  }
  
  # Open PDF if output file is specified
  if (!is.null(output_file)) {
    pdf(output_file, width = width, height = height)
    on.exit(dev.off())
  }
  
  # Create a layout with 24 chromosomes (5x5 grid) plus legend
  layout_matrix <- matrix(c(1:24, rep(25, 1)), nrow = 5, byrow = TRUE)
  layout(layout_matrix)
  
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
  
  return(NULL)
}

#' Plot Chromosome Arm Difference Heatmap
#'
#' Creates a heatmap showing differences in CNA values between functional classes at the chromosome arm level.
#'
#' @param arm_summary Output from aggregateByChromosomeArm()
#' @param class_columns Names of columns for each functional class
#' @param reference_class Optional reference class to compare others against. If NULL, all pairwise comparisons are shown.
#' @param min_value Minimum value for color scale. If NULL, will use data minimum.
#' @param max_value Maximum value for color scale. If NULL, will use data maximum.
#' @param cluster_rows Whether to cluster rows. Default is FALSE
#' @param title Plot title. Default is "CNA Differences by Chromosome Arm"
#'
#' @return A pheatmap object
#' @export
#'
#' @examples
#' \dontrun{
#' cna_data <- extractCNAData(seurat_obj)
#' profiles <- calculateCNAProfiles(cna_data)
#' arm_summary <- aggregateByChromosomeArm(profiles)
#' 
#' # Compare all classes pairwise
#' plotChromosomeDifferenceHeatmap(arm_summary, 
#'                               class_columns = c("hypoxic", "proliferative", "undefined"))
#'                              
#' # Compare against a reference class
#' plotChromosomeDifferenceHeatmap(arm_summary, 
#'                               class_columns = c("hypoxic", "proliferative", "undefined"),
#'                               reference_class = "undefined")
#' }
plotChromosomeDifferenceHeatmap <- function(arm_summary,
                                          class_columns,
                                          reference_class = NULL,
                                          min_value = NULL,
                                          max_value = NULL,
                                          cluster_rows = FALSE,
                                          title = "CNA Differences by Chromosome Arm") {
  
  # Check if required packages are installed
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package pheatmap must be installed")
  }
  
  # Check if we have class columns
  if (length(class_columns) < 2) {
    stop("At least two class columns are needed for difference calculation")
  }
  
  # Verify all class columns exist in the data
  missing_cols <- setdiff(class_columns, colnames(arm_summary))
  if (length(missing_cols) > 0) {
    stop("The following class columns are missing from arm_summary: ", 
         paste(missing_cols, collapse = ", "))
  }
  
  # If reference class is specified, check if it's valid
  if (!is.null(reference_class)) {
    if (!reference_class %in% class_columns) {
      stop("reference_class must be one of the class_columns")
    }
  }
  
  # Create title case function for column names
  title_case <- function(x) {
    paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  }
  
  # Prepare the difference matrix
  difference_data <- data.frame(
    chrom_arm = arm_summary$chrom_arm,
    stringsAsFactors = FALSE
  )
  
  # Calculate differences based on reference class or all pairwise
  if (!is.null(reference_class)) {
    # Compare all classes to the reference class
    classes_to_compare <- setdiff(class_columns, reference_class)
    for (class in classes_to_compare) {
      col_name <- paste0(title_case(class), "_vs_", title_case(reference_class))
      difference_data[[col_name]] <- arm_summary[[class]] - arm_summary[[reference_class]]
    }
  } else {
    # All pairwise comparisons
    for (i in 1:(length(class_columns)-1)) {
      for (j in (i+1):length(class_columns)) {
        class1 <- class_columns[i]
        class2 <- class_columns[j]
        col_name <- paste0(title_case(class1), "_vs_", title_case(class2))
        difference_data[[col_name]] <- arm_summary[[class1]] - arm_summary[[class2]]
      }
    }
  }
  
  # Extract the matrix for heatmap
  diff_cols <- setdiff(colnames(difference_data), "chrom_arm")
  heatmap_data <- as.matrix(difference_data[, diff_cols, drop = FALSE])
  rownames(heatmap_data) <- difference_data$chrom_arm
  
  # Set color scale range
  if (is.null(min_value) || is.null(max_value)) {
    max_abs <- max(abs(heatmap_data), na.rm = TRUE)
    
    if (is.null(min_value)) min_value <- -max_abs
    if (is.null(max_value)) max_value <- max_abs
  }
  
  # Create annotation for chromosome arms
  annotation_row <- data.frame(
    Chromosome = sub("([0-9XY]+)[pq]", "\\1", rownames(heatmap_data)),
    Arm = sub("[0-9XY]+([pq])", "\\1", rownames(heatmap_data))
  )
  rownames(annotation_row) <- rownames(heatmap_data)
  
  # Define annotation colors
  annotation_colors <- list(
    Arm = c(p = "#f08080", q = "#87cefa")
  )
  
  # Create the heatmap
  heatmap_plot <- pheatmap::pheatmap(
    heatmap_data,
    cluster_rows = cluster_rows,
    cluster_cols = FALSE,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = seq(min_value, max_value, length.out = 101),
    cellwidth = 50,
    cellheight = 12,
    main = title,
    silent = TRUE  # Return the object instead of plotting directly
  )
  
  return(heatmap_plot)
}

#' Create Summary Excel File of CNA Results
#'
#' Creates a comprehensive Excel file with multiple sheets summarizing CNA findings.
#'
#' @param arm_summary Output from aggregateByChromosomeArm()
#' @param class_columns Names of columns for each functional class
#' @param output_file Path to the output Excel file. Default is "CNA_Analysis_Findings.xlsx"
#' @param significance_threshold Threshold for considering a CNA significant. Default is 0.02
#' @param reference_class Optional reference class to compare others against
#' @param key_genes Optional data frame with key genes by chromosome arm
#'
#' @return NULL (creates a file)
#' @export
#'
#' @examples
#' \dontrun{
#' cna_data <- extractCNAData(seurat_obj)
#' profiles <- calculateCNAProfiles(cna_data, 
#'                                 classes = c("Hypoxic", "Proliferative", "Undefined"))
#' arm_summary <- aggregateByChromosomeArm(profiles)
#' 
#' # Standard analysis
#' createCNASummaryExcel(arm_summary, 
#'                      class_columns = c("hypoxic", "proliferative", "undefined"),
#'                      output_file = "results/CNA_summary.xlsx")
#' 
#' # With custom reference class
#' createCNASummaryExcel(arm_summary, 
#'                      class_columns = c("emt_high", "stemness_high", "basal"),
#'                      reference_class = "basal",
#'                      output_file = "results/CNA_custom_summary.xlsx")
#' }
createCNASummaryExcel <- function(arm_summary,
                                 class_columns,
                                 output_file = "CNA_Analysis_Findings.xlsx",
                                 significance_threshold = 0.02,
                                 reference_class = NULL,
                                 key_genes = NULL) {
  
  # Check if required packages are available
  if (!requireNamespace("writexl", quietly = TRUE)) {
    stop("Package writexl must be installed")
  }
  
  # Check if we have class columns
  if (length(class_columns) == 0) {
    stop("No class columns specified")
  }
  
  # Verify all class columns exist in the data
  missing_cols <- setdiff(class_columns, colnames(arm_summary))
  if (length(missing_cols) > 0) {
    stop("The following class columns are missing from arm_summary: ", 
         paste(missing_cols, collapse = ", "))
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
  
  # ===== SHEET 1-N: TOP CNA REGIONS BY CATEGORY =====
  top_tables <- list()
  
  for (class in class_columns) {
    # Get top 10 altered regions for this category (or fewer if not enough data)
    n_regions <- min(10, nrow(arm_summary))
    top_altered <- arm_summary[order(abs(arm_summary[[class]]), decreasing = TRUE), ][1:n_regions, ]
    
    # Prepare top regions for Excel
    top_table <- data.frame(
      Rank = 1:n_regions,
      Chromosome_Arm = top_altered$chrom_arm,
      CNA_Value = round(top_altered[[class]], 4),
      Status = ifelse(top_altered[[class]] > 0, "Gain", "Loss"),
      Chromosome = top_altered$chrom,
      Arm = top_altered$arm
    )
    
    top_tables[[paste0("Top_", title_case(class), "_CNAs")]] <- top_table
  }
  
  # ===== SHEET N+1: ALL CHROMOSOME ARMS WITH CNA VALUES =====
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
  all_arms_table <- all_arms_table[order(
    as.numeric(factor(all_arms_table$Chromosome, levels = c(1:22, "X", "Y"))),
    all_arms_table$Arm), ]
  
  # ===== SHEET N+2,N+3: SIGNIFICANT GAINS & LOSSES =====
  all_gains <- data.frame()
  all_losses <- data.frame()
  
  for (class in class_columns) {
    # Filter for significant gains (CNA > threshold)
    gains <- all_arms_table[all_arms_table[[paste0(title_case(class), "_CNA")]] > significance_threshold, 
                           c("Chromosome_Arm", paste0(title_case(class), "_CNA"))]
    if (nrow(gains) > 0) {
      colnames(gains)[2] <- "CNA_Value"
      gains$Category <- title_case(class)
      
      # Add to combined table
      all_gains <- rbind(all_gains, gains)
    }
    
    # Filter for significant losses (CNA < -threshold)
    losses <- all_arms_table[all_arms_table[[paste0(title_case(class), "_CNA")]] < -significance_threshold, 
                            c("Chromosome_Arm", paste0(title_case(class), "_CNA"))]
    if (nrow(losses) > 0) {
      colnames(losses)[2] <- "CNA_Value"
      losses$Category <- title_case(class)
      
      # Add to combined table
      all_losses <- rbind(all_losses, losses)
    }
  }
  
  # Sort the tables
  if (nrow(all_gains) > 0) {
    all_gains <- all_gains[order(all_gains$CNA_Value, decreasing = TRUE), ]
  }
  if (nrow(all_losses) > 0) {
    all_losses <- all_losses[order(all_losses$CNA_Value), ]  # Ascending for losses (most negative first)
  }
  
  # ===== SHEET N+4: KEY CANCER GENE REGIONS =====
  # Use provided key genes table or create default
  if (is.null(key_genes)) {
    # Define default key cancer genes by chromosome arm
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
      ),
      stringsAsFactors = FALSE
    )
  } else {
    key_cancer_regions <- key_genes
  }
  
  # Merge with CNA data
  key_regions_cna <- merge(key_cancer_regions, all_arms_table, by = "Chromosome_Arm", all.x = TRUE)
  
  # Sort by average absolute CNA value
  cna_cols <- grep("_CNA$", colnames(key_regions_cna), value = TRUE)
  if (length(cna_cols) > 0) {
    avg_abs_cna <- rowMeans(abs(as.matrix(key_regions_cna[, cna_cols, drop = FALSE])), na.rm = TRUE)
    key_regions_cna$Average_Abs_CNA <- avg_abs_cna
    key_regions_cna <- key_regions_cna[order(key_regions_cna$Average_Abs_CNA, decreasing = TRUE), ]
  }
  
  # ===== SHEET N+5: DIFFERENCES BETWEEN CATEGORIES =====
  # Calculate differences between categories
  if (length(class_columns) > 1) {
    category_differences <- data.frame(
      Chromosome_Arm = arm_summary$chrom_arm,
      stringsAsFactors = FALSE
    )
    
    # Calculate differences based on reference class or all pairwise
    if (!is.null(reference_class) && reference_class %in% class_columns) {
      # Compare all classes to the reference class
      classes_to_compare <- setdiff(class_columns, reference_class)
      for (class in classes_to_compare) {
        col_name <- paste0(title_case(class), "_vs_", title_case(reference_class))
        category_differences[[col_name]] <- round(arm_summary[[class]] - arm_summary[[reference_class]], 4)
      }
    } else {
      # All pairwise comparisons
      for (i in 1:(length(class_columns)-1)) {
        for (j in (i+1):length(class_columns)) {
          class1 <- class_columns[i]
          class2 <- class_columns[j]
          col_name <- paste0(title_case(class1), "_vs_", title_case(class2))
          
          if (paste0(class1, "_vs_", class2) %in% colnames(arm_summary)) {
            category_differences[[col_name]] <- round(arm_summary[[paste0(class1, "_vs_", class2)]], 4)
          } else {
            category_differences[[col_name]] <- round(arm_summary[[class1]] - arm_summary[[class2]], 4)
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
  }
  
  # ===== COMBINE ALL SHEETS INTO EXCEL FILE =====
  sheets_list <- c(
    top_tables,
    list(
      "All_Chromosome_Arms" = all_arms_table,
      "Key_Cancer_Gene_Regions" = key_regions_cna
    )
  )
  
  # Add gains and losses sheets if they have content
  if (nrow(all_gains) > 0) {
    sheets_list$Significant_Gains <- all_gains
  }
  if (nrow(all_losses) > 0) {
    sheets_list$Significant_Losses <- all_losses
  }
  
  # Add differences sheet if it exists and has content
  if (exists("category_differences") && ncol(category_differences) > 1) {
    sheets_list$Category_Differences <- category_differences
  }
  
  # Write to Excel
  writexl::write_xlsx(sheets_list, path = output_file)
  
  message(paste("Created comprehensive CNA analysis summary Excel file at:", output_file))
  
  return(NULL)
}

#' Plot Interactive CNA Visualization Dashboard
#'
#' Creates an interactive dashboard for exploring CNA data by functional class.
#'
#' @param seurat_obj A Seurat object with CNA data
#' @param class_column Name of the column containing functional classifications
#' @param cna_assay Name of the assay containing CNA data. Default is "cna"
#' @param title Dashboard title. Default is "CNA Analysis by Functional Class"
#'
#' @return A Shiny app object that can be run with shinyApp()
#' @export
#'
 #' @examples
#' \dontrun{
#' # Create and run the dashboard
#' cna_dashboard <- plotInteractiveCNADashboard(seurat_obj, 
#'                                            class_column = "functional_class")
#' if (interactive()) {
#'   shiny::runApp(cna_dashboard)
#' }
#' }
plotInteractiveCNADashboard <- function(seurat_obj,
                                       class_column,
                                       cna_assay = "cna",
                                       title = "CNA Analysis by Functional Class") {
  
  # Check if required packages are installed
  if (!requireNamespace("shiny", quietly = TRUE) || 
      !requireNamespace("Seurat", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("DT", quietly = TRUE)) {
    stop("Packages shiny, Seurat, ggplot2, and DT must be installed")
  }
  
  # Extract CNA data and functional class information
  if (!cna_assay %in% names(seurat_obj@assays)) {
    stop("CNA assay not found in Seurat object")
  }
  
  if (!class_column %in% colnames(seurat_obj@meta.data)) {
    stop("Class column not found in Seurat metadata")
  }
  
  # Define UI
  ui <- shiny::fluidPage(
    shiny::titlePanel(title),
    
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::selectInput("chromosome", "Select Chromosome:", 
                          c("All", as.character(1:22), "X", "Y")),
        shiny::checkboxGroupInput("classes", "Select Functional Classes:",
                                 choices = unique(seurat_obj@meta.data[[class_column]]),
                                 selected = unique(seurat_obj@meta.data[[class_column]])[1]),
        shiny::sliderInput("cna_range", "CNA Value Range:",
                          min = -3, max = 3, value = c(-0.5, 0.5)),
        shiny::hr(),
        shiny::downloadButton("download_data", "Download Data"),
        shiny::hr(),
        shiny::helpText("This dashboard allows you to explore copy number alterations (CNAs) by functional class. 
                        Use the controls to filter and analyze the data.")
      ),
      
      shiny::mainPanel(
        shiny::tabsetPanel(
          shiny::tabPanel("Summary", 
                         shiny::plotOutput("summary_plot"),
                         shiny::hr(),
                         shiny::h4("Class Statistics"),
                         DT::dataTableOutput("class_stats")),
          shiny::tabPanel("Chromosome View",
                         shiny::plotOutput("chromosome_plot"),
                         shiny::hr(),
                         shiny::h4("Top Altered Genes"),
                         DT::dataTableOutput("top_genes")),
          shiny::tabPanel("Gene-Level Analysis",
                         shiny::selectizeInput("gene_search", "Search for Gene:", 
                                              choices = NULL, options = list(create = TRUE)),
                         shiny::plotOutput("gene_plot"),
                         shiny::hr(),
                         shiny::h4("Gene Details"),
                         shiny::verbatimTextOutput("gene_details")),
          shiny::tabPanel("Data Table",
                         DT::dataTableOutput("full_data"))
        )
      )
    )
  )
  
  # Define server logic
  server <- function(input, output, session) {
    # This is a placeholder for the actual server logic
    # In a real implementation, you would process the data and create all the plots and tables
    
    # Example output
    output$summary_plot <- shiny::renderPlot({
      ggplot2::ggplot() + 
        ggplot2::geom_text(ggplot2::aes(x = 0.5, y = 0.5), 
                          label = "CNA Summary Plot\n(Implementation details would depend on specific analysis needs)") +
        ggplot2::theme_void()
    })
    
    # Placeholder for data download
    output$download_data <- shiny::downloadHandler(
      filename = function() {
        paste("cna_analysis_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        # In a real implementation, you would save the processed data here
        write.csv(data.frame(Note = "Placeholder for CNA data"), file)
      }
    )
  }
  
  # Return the Shiny app
  return(shiny::shinyApp(ui = ui, server = server))
}