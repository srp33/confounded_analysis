# plot_reduced.R
#
# This script reads all coordinate files in a given dataset directory,
# and uses the metadata already present in those files to generate plots.
# It saves all plots for the dataset to a single PDF.

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(argparse)
  library(tools)
  library(tibble)
})

# This tibble centralizes metadata for different datasets.
# It is used to dynamically set column names for plotting.
all_info_df <- tribble(
  ~dataset,    ~title,                            ~batch_label,      ~true_label,
  "gse20194", "GSE 20194 ER",                   "meta_batch",      "meta_er_status",
  "gse20194", "GSE 20194 HER2",                 "meta_batch",      "meta_her2_status",
  "gse20194", "GSE 20194 PR",                   "meta_batch",      "meta_pr_status",
  "gse24080",  "GSE 24080 Cytogenetic Abnormality",    "meta_batch",      "meta_cytogenetic_abnormality",
  "gse49711",  "GSE 49711 Stage 3 4",           "meta_Sex",      "meta_INSS_Stage_Split_3_4"
)

# --- Main Execution Block ---
main <- function() {
  parser <- ArgumentParser(description = "Generate plots from pre-computed coordinate files.")
  parser$add_argument("-i", "--input-dir", required = TRUE, help = "Directory containing coordinate CSV files for a single dataset.")
  parser$add_argument("-o", "--output-pdf", required = TRUE, help = "Path to save the output PDF.")
  args <- parser$parse_args()

  # The dataset ID is the name of the input directory
  dataset_id <- basename(args$input_dir)

  # Find all coordinate files in the input directory
  coord_files <- list.files(path = args$input_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
  if (length(coord_files) == 0) {
    stop(sprintf("No coordinate CSV files found in: %s.", args$input_dir))
  }
  
  # --- Find Relevant Metadata Rows ---
  dataset_info_rows <- all_info_df[all_info_df$dataset == dataset_id, ]
  if (nrow(dataset_info_rows) == 0) {
      stop(sprintf("FATAL: No entry for dataset '%s' in the R script's all_info_df.", dataset_id))
  }

  plot_list <- list()

  for (file_path in coord_files) {
    plot_data <- readr::read_csv(file_path, show_col_types = FALSE)

    # A single file can be used to generate multiple plots (e.g., for ER, PR, HER2)
    for (i in seq_len(nrow(dataset_info_rows))) {
      current_info <- dataset_info_rows[i, ]
      
      message(paste("Generating plot for:", basename(file_path), "| Title:", current_info$title))

      batch_col_name <- current_info$batch_label
      signal_col_name <- current_info$true_label
      plot_title <- current_info$title
      clean_batch_label <- tools::toTitleCase(gsub("_", " ", batch_col_name))
      clean_signal_label <- tools::toTitleCase(gsub("_", " ", signal_col_name))
      
      filename_parts <- strsplit(file_path_sans_ext(basename(file_path)), "-")[[1]]
      method_name <- filename_parts[1]
      plot_type <- toupper(filename_parts[2])

      # Check if required columns exist before plotting
      required_cols <- c(batch_col_name, signal_col_name)
      if (!all(required_cols %in% colnames(plot_data))) {
        missing_cols <- required_cols[!required_cols %in% colnames(plot_data)]
        message(paste("Warning: Skipping plot '", plot_title, "' for file", basename(file_path), "because required column(s) are missing:", paste(missing_cols, collapse=", ")))
        next
      }

      # Since the true labels are binary, we can pick two highly distinct shapes.
      # A solid circle (16) and a bold cross (4) are an excellent, clear pair.
      binary_shapes <- c(16, 4)

      p <- ggplot(plot_data, aes(x = Dim1, y = Dim2, color = factor(.data[[batch_col_name]]), shape = factor(.data[[signal_col_name]]))) +
        # Use a larger size and a thick stroke to make the cross shape bold.
        geom_point(size = 5, alpha = 0.8, stroke = 1.5) +
        # Apply the custom shape palette for two classes.
        scale_shape_manual(values = binary_shapes) +
        labs(
          title = paste(plot_title, "-", tools::toTitleCase(method_name), paste0("(", plot_type, ")")),
          subtitle = paste("Colored by", clean_batch_label, ", Shaped by", clean_signal_label),
          x = paste0(plot_type, "-1"),
          y = paste0(plot_type, "-2"),
          color = clean_batch_label,
          shape = clean_signal_label
        ) +
        theme_bw(base_size = 14) +
        theme(
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)
        ) +
        # Ensure the legend glyphs are large and clear.
        guides(color = guide_legend(override.aes = list(size=5)),
               shape = guide_legend(override.aes = list(size=5, stroke = 1.5)))
      
      plot_list[[paste(basename(file_path), plot_title)]] <- p
    }
  }

  if (length(plot_list) == 0) {
    stop("No plots were generated. Check for warnings about missing columns.")
  }

  message(paste("Saving", length(plot_list), "plots to:", args$output_pdf))
  
  pdf(args$output_pdf, width = 11, height = 8.5)
  for (plot in plot_list) {
    print(plot)
  }
  dev.off()

  message("---")
  message(paste("Successfully generated:", args$output_pdf))
}

if (sys.nframe() == 0) {
  main()
}