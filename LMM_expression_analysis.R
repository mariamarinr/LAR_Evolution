# Gene Expression Analysis with Linear Mixed-Effects Models
# Author: Maria F. Marin-Recinos
# Description: This script performs log-transformed expression analysis using linear mixed-effects models, 
# post-hoc comparisons, and visualizes the results using boxplots and violin plots.
# Adaptable for any gene set across tissues and samples.

# Load necessary libraries
library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(tidyr)
library(emmeans)

# Path to input CSV (must contain: 'Sample', 'tissue', and gene expression columns)
input_data <- "path/to/expression_data.csv"  # Replace with your actual path

# Gene prefix for filtering (e.g., "LAR", "MYB", etc.)
gene_prefix <- "LAR"

# Output directory
output_dir <- "output"

# Plot title (e.g., species name)
plot_title <- "Species name"

# -----------------------------

# Output file paths
output_plot_png <- file.path(output_dir, "GE-plot.png")
output_plot_svg <- file.path(output_dir, "GE-plot.svg")
output_model_summary <- file.path(output_dir, "model_summary.txt")
output_emmeans_results <- file.path(output_dir, "emmeans_results.txt")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load data
data <- read.csv(input_data)
head(data)

# Reshape to long format based on selected gene prefix
data_long <- data %>%
  pivot_longer(cols = starts_with(gene_prefix), names_to = "Gene", values_to = "TPM") %>%
  mutate(GeneGroup = case_when(
    grepl("1", Gene) ~ paste0(gene_prefix, "1"),
    grepl("2", Gene) ~ paste0(gene_prefix, "2"),
    TRUE ~ Gene  # fallback
  ),
  tissue_Gene = paste(tissue, GeneGroup, sep = "-"))

# Log transformation
data_long <- data_long %>%
  mutate(TPM_log = log(TPM + 1))

# Fit Linear Mixed-Effects Model
model <- lmer(TPM_log ~ GeneGroup * tissue + (1 | Sample), data = data_long)
model_summary <- summary(model)

# Post-hoc pairwise comparisons
emmeans_results <- emmeans(model, pairwise ~ GeneGroup | tissue)

# Save model outputs
capture.output(model_summary, file = output_model_summary)
capture.output(summary(emmeans_results), file = output_emmeans_results)

# Plot boxplot and violin plot
plot_species <- ggplot(data_long, aes(x = tissue, y = TPM_log, fill = GeneGroup)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.25, alpha = 0.7) +
  geom_violin(position = position_dodge(width = 0.75), alpha = 0.3) +
  labs(title = plot_title, x = "Tissue", y = "Log(TPM + 1)") +
  theme_classic() +
  theme(
    plot.title = element_text(face = "italic", hjust = 0.5, size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 11),
    legend.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

# Save plots
ggsave(output_plot_png, plot = plot_species, width = 12, height = 5, dpi = 300, units = "in")
ggsave(output_plot_svg, plot = plot_species, width = 12, height = 5, units = "in")

# Model diagnostics (optional visualization)
par(mfrow = c(2, 2))
plot(model)

# Predicted vs. Actual
data_long <- data_long %>% mutate(predicted = predict(model))

ggplot(data_long, aes(x = predicted, y = TPM_log)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "red") +
  labs(title = "Predicted vs. Actual Log-Transformed TPM Values", 
       x = "Predicted Values", y = "Log(TPM + 1)")
