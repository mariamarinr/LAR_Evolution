# Load required library
library(tidyverse)

# File paths — update these for your data
genelist_path <- "/path/to/gene_list.csv"         # # One gene ID per line
countdata_path <- "/path/to/count_table.txt"       # Tab-delimited count matrix with 'gene' column and sample IDs
metadata_path <- "/path/to/metadata.csv"          # e.g., CSV with 'Run' and 'tissue' columns

# Output files
expression_output <- "filtered_expression_data.csv" # Output expression matrix
merged_output <- "merged_expression_metadata.csv"   # Output merged data

# Step 1: Extract gene expression data for selected genes
extract_expression_data <- function(genelist_path, countdata_path, output_path) {
  # Read gene list and count data
  gene_list <- read_csv(genelist_path, col_names = c("GeneID", "Alias")) %>%
    mutate(GeneID = str_trim(GeneID))

  count_data <- read_tsv(countdata_path, show_col_types = FALSE) %>%
    rename(GeneID = gene) %>%
    mutate(GeneID = str_trim(GeneID))

  # Filter for selected genes
  filtered <- count_data %>%
    filter(GeneID %in% gene_list$GeneID) %>%
    column_to_rownames("GeneID") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Sample")

  write_csv(filtered, output_path)
  return(filtered)
}

# Step 2: Merge expression data with metadata by Sample ID ("Run")
merge_with_metadata <- function(expression_path, metadata_path, output_path) {
  expression_data <- read_csv(expression_path, show_col_types = FALSE)

  metadata <- read_csv(metadata_path, show_col_types = FALSE) %>%
    rename(Sample = Run) %>%  # Assumes 'Run' column holds sample IDs matching expression data
    select(Sample, tissue)    # Keep only relevant columns

  merged <- inner_join(expression_data, metadata, by = "Sample")
  write_csv(merged, output_path)
  return(merged)
}

# Run the workflow
expression_data <- extract_expression_data(genelist_path, countdata_path, expression_output)
merged_data <- merge_with_metadata(expression_output, metadata_path, merged_output)
