head(meta_result$Gene_ID)

gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = meta_result$Gene_ID,
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL")
)

head(gene_map)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db", ask = FALSE, update = FALSE)
}

library(org.Hs.eg.db)
library(AnnotationDbi)
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = meta_result$Gene_ID,
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL")
)
head(gene_map)

gene_map_clean <- gene_map |>
  dplyr::filter(!is.na(SYMBOL), SYMBOL != "") |>
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)
meta_result_symbol <- meta_result |>
  dplyr::left_join(
    gene_map_clean |>
      dplyr::select(ENSEMBL, SYMBOL),
    by = c("Gene_ID" = "ENSEMBL")
  ) |>
  dplyr::rename(
    Gene_Symbol = SYMBOL
  )

head(
  meta_result_symbol[, c("Gene_ID", "Gene_Symbol")]
)

# ============================================================
# Volcano plot visualizing meta-analysis results
# Colorectal cancer
# Random-effects meta-analysis
# ============================================================

library(tidyverse)
library(ggrepel)
library(rio)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ------------------------------------------------------------
# Create output directories
# ------------------------------------------------------------

dir.create(
  "results/figures/meta-analysis",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "results/tables/meta-analysis",
  showWarnings = FALSE,
  recursive = TRUE
)

# ------------------------------------------------------------
# Load random-effects meta-analysis results
# ------------------------------------------------------------

meta_result <- import(
  "results/tables/meta-analysis/random_effect_model.csv"
)

# ------------------------------------------------------------
# Keep required meta-analysis columns
# ------------------------------------------------------------

meta_result <- meta_result |>
  dplyr::select(
    Gene_ID,
    randomSummary,
    randomP
  ) |>
  dplyr::rename(
    log2FC = randomSummary,
    P.Value = randomP
  )

# ------------------------------------------------------------
# Remove invalid values
# ------------------------------------------------------------

meta_result <- meta_result |>
  dplyr::filter(
    !is.na(log2FC),
    !is.na(P.Value),
    P.Value > 0,
    P.Value <= 1
  )

# ------------------------------------------------------------
# Adjust P-values using Benjamini-Hochberg method
# ------------------------------------------------------------

meta_result <- meta_result |>
  dplyr::mutate(
    P.adj = p.adjust(
      P.Value,
      method = "BH"
    )
  )

# ------------------------------------------------------------
# Gene ID -> Gene Symbol
# ------------------------------------------------------------

gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = meta_result$Gene_ID,
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL")
)

# Remove duplicate mappings
gene_map_clean <- gene_map |>
  dplyr::filter(
    !is.na(SYMBOL),
    SYMBOL != ""
  ) |>
  dplyr::distinct(
    ENSEMBL,
    .keep_all = TRUE
  )

# Add Gene Symbol
meta_result <- meta_result |>
  dplyr::left_join(
    gene_map_clean |>
      dplyr::select(
        ENSEMBL,
        SYMBOL
      ),
    by = c("Gene_ID" = "ENSEMBL")
  ) |>
  dplyr::rename(
    Gene_Symbol = SYMBOL
  )

# ------------------------------------------------------------
# Define significance
# ------------------------------------------------------------

meta_result <- meta_result |>
  dplyr::mutate(
    Significance = dplyr::case_when(
      P.adj < 0.05 & log2FC > 1 ~ "Up",
      P.adj < 0.05 & log2FC < -1 ~ "Down",
      TRUE ~ "Not Significant"
    )
  )

# ------------------------------------------------------------
# Top significant genes for labeling
# ------------------------------------------------------------

top_genes <- meta_result |>
  dplyr::filter(
    !is.na(Gene_Symbol),
    Gene_Symbol != "",
    P.adj < 0.05,
    abs(log2FC) > 1
  ) |>
  dplyr::slice_max(
    order_by = abs(log2FC),
    n = 30
  )

# ------------------------------------------------------------
# Symmetric x-axis
# ------------------------------------------------------------

x_lim <- ceiling(
  max(
    abs(meta_result$log2FC),
    na.rm = TRUE
  )
)

x_lim <- max(x_lim, 2)

x_break <- if (x_lim > 10) 2 else 1

# ------------------------------------------------------------
# Volcano plot
# ------------------------------------------------------------

volcano <- ggplot(
  meta_result,
  aes(
    x = log2FC,
    y = -log10(P.adj),
    color = Significance
  )
) +
  
  geom_point(
    alpha = 0.8,
    size = 2.5
  ) +
  
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "black"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "black"
  ) +
  
  scale_color_manual(
    values = c(
      "Down" = "#2c7fb8",
      "Not Significant" = "#636363",
      "Up" = "#e34a33"
    )
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Meta-analysis Volcano Plot",
    x = "log2 (Fold Change)",
    y = "-log10 (Adjusted P-value)",
    color = "Significance"
  ) +
  
  theme(
    legend.title = element_text(
      size = 15,
      face = "bold"
    ),
    legend.text = element_text(
      size = 14
    ),
    legend.position = "right",
    axis.title.x = element_text(
      size = 14
    ),
    axis.title.y = element_text(
      size = 14
    ),
    axis.text.x = element_text(
      colour = "black",
      size = 12
    ),
    axis.text.y = element_text(
      colour = "black",
      size = 12
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    )
  ) +
  
  geom_text_repel(
    data = top_genes,
    aes(
      label = Gene_Symbol
    ),
    colour = "black",
    size = 3,
    max.overlaps = 30,
    direction = "both",
    max.time = 5,
    force = 5,
    force_pull = 5,
    point.padding = 0.5,
    seed = 40
  ) +
  
  scale_x_continuous(
    limits = c(-x_lim, x_lim),
    breaks = seq(
      -x_lim,
      x_lim,
      by = x_break
    )
  )

# ------------------------------------------------------------
# Display plot
# ------------------------------------------------------------

print(volcano)

# ------------------------------------------------------------
# Save volcano plot
# ------------------------------------------------------------

ggsave(
  "results/figures/meta-analysis/volcano_plot.png",
  plot = volcano,
  width = 12,
  height = 12,
  dpi = 600,
  bg = "white"
)

# ------------------------------------------------------------
# Save result table
# ------------------------------------------------------------

write.csv(
  meta_result,
  "results/tables/meta-analysis/volcano_plot_results.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

cat("\nVolcano plot completed successfully.\n")

cat(
  "Figure saved to:\n",
  "results/figures/meta-analysis/volcano_plot.png\n"
)

cat(
  "Result table saved to:\n",
  "results/tables/meta-analysis/volcano_plot_results.csv\n"
)

cat(
  "\nNumber of Up-regulated genes:",
  sum(meta_result$Significance == "Up"),
  "\n"
)

cat(
  "Number of Down-regulated genes:",
  sum(meta_result$Significance == "Down"),
  "\n"
)

cat(
  "Number of Not Significant genes:",
  sum(meta_result$Significance == "Not Significant"),
  "\n"
)

source("scripts/07_volcano_plot.R")


browseURL(
  normalizePath(
    "results/figures/meta-analysis/volcano_plot.png"
  )
)


