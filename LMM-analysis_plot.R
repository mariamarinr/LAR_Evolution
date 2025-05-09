# Load necessary libraries
library(lme4)       # Mixed-effects models
library(lmerTest)   # P-values in mixed models
library(ggplot2)    # Visualization
library(dplyr)      # Data manipulation
library(tidyr)      # Data reshaping
library(emmeans)    # Post-hoc analysis

# Define file paths
input_data <- "path/to/LAR-GE_data.csv" #table containing sample, tissue-type, genes with count data
output_plot_png <- "output/GE-plot.png"
output_plot_svg <- "output/GE-plot.svg"
output_model_summary <- "output/model_summary.txt"
output_emmeans_results <- "output/emmeans_results.txt"

# Load data
data <- read.csv(input_data)
head(data)  # Inspect first few rows

# Step 1: Reshape Data (Long Format)
data_long <- data %>%
  pivot_longer(cols = starts_with("LAR"), names_to = "Gene", values_to = "TPM") %>%
  mutate(Genes = ifelse(grepl("LAR1", Gene), "LAR1", "LAR2"),
         tissue_Gene = paste(tissue, Genes, sep = "-"))

# Step 2: Log Transformation
data_long <- data_long %>% mutate(TPM_log = log(TPM + 1))

# Step 3: Linear Mixed-Effects Model (LMM)
model <- lmer(TPM_log ~ Genes * tissue + (1 | Sample), data = data_long)
summary(model)

# Step 4: Post-Hoc Analysis
emmeans_results <- emmeans(model, pairwise ~ Genes | tissue)
summary(emmeans_results)

# Save Model Results
capture.output(summary(model), file = output_model_summary)
capture.output(summary(emmeans_results), file = output_emmeans_results)

# Step 5: Create Boxplot & Violin Plot
plot_species <- ggplot(data_long, aes(x = tissue, y = TPM_log, fill = Genes)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.25, alpha = 0.7) +
  geom_violin(position = position_dodge(width = 0.75), alpha = 0.3) +
  scale_fill_manual(values = c("LAR1" = "#cb4154", "LAR2" = "#4682b4"),
                    labels = c(expression(italic("LAR1")), expression(italic("LAR2")))) +
  labs(title = "Species name", x = "tissue", y = "Log(TPM + 1)") +
  theme_classic() +
  theme(plot.title = element_text(face = "italic", hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

ggsave(output_plot_png, plot = plot_species, width = 12, height = 5, dpi = 300, units = "in")
ggsave(output_plot_svg, plot = plot_species, width = 12, height = 5, units = "in")

# Step 6: Model Diagnostics
par(mfrow = c(2, 2))
plot(model)

# Step 7: Visualize Predicted vs. Actual
predicted_values <- predict(model)
data_long <- data_long %>% mutate(predicted = predicted_values)

ggplot(data_long, aes(x = predicted, y = TPM_log)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "red") +
  labs(title = "Predicted vs. Actual Log-Transformed TPM Values", 
       x = "Predicted Values", y = "Log(TPM + 1)")

