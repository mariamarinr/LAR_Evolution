# Load required libraries
library(tidyverse)

# File paths — update these for your data
gene_list_file   <- "data/gene_list.csv"        # One gene ID per line
count_data_file  <- "data/count_matrix.tsv"     # Tab-delimited count matrix: Gene x Samples
metadata_file    <- "data/sample_metadata.csv"  # First column = sample IDs

# Output files
expression_output_file <- "results/filtered_expression_matrix.csv"
merged_output_file     <- "results/merged_expression_metadata.csv"

# Step 1: Filter expression data by gene list
filter_expression_data <- function(gene_list_path, count_data_path, output_path) {
  gene_list <- read_lines(gene_list_path) %>% str_trim()
  
  count_data <- read_tsv(count_data_path, show_col_types = FALSE) %>%
    rename(GeneID = 1) %>%  # Rename first column to 'GeneID'
    mutate(GeneID = str_trim(as.character(GeneID))) %>%
    filter(GeneID %in% gene_list)
  
  transposed <- count_data %>%
    column_to_rownames("GeneID") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Sample")
  
  write_csv(transposed, output_path)
  return(transposed)
}

# Step 2: Merge expression matrix with metadata
merge_expression_with_metadata <- function(expression_path, metadata_path, output_path) {
  expression_data <- read_csv(expression_path, show_col_types = FALSE)
  
  metadata <- read_csv(metadata_path, show_col_types = FALSE) %>%
    rename(Sample = 1)  # Rename first column to 'Sample'
  
  merged <- inner_join(expression_data, metadata, by = "Sample")
  
  write_csv(merged, output_path)
  return(merged)
}

# Create output folders if needed
dir.create("results", showWarnings = FALSE, recursive = TRUE)

filtered_expression <- filter_expression_data(
  gene_list_path   = gene_list_file,
  count_data_path  = count_data_file,
  output_path      = expression_output_file
)

merged_data <- merge_expression_with_metadata(
  expression_path = expression_output_file,
  metadata_path   = metadata_file,
  output_path     = merged_output_file
)
