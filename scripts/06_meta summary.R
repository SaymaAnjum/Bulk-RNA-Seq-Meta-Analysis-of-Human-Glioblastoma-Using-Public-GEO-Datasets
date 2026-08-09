# ============================================================
# 06_meta summary.R
# ============================================================
# Summary of Random-Effects REML Meta-analysis
# Glioblastoma cancer: Tumor vs Normal

# Input:
#   results/tables/meta-analysis/random_effect_model.csv
#
# Output:
#   - FDR-adjusted meta-analysis results
#   - Significant meta-DEG tables
#   - Up/Down summary statistics
#   - Meta-analysis volcano plot
#   - Venn intersection cross-reference (if available)
# ============================================================


# ============================================================
# 01. Load packages
# ============================================================

library(tidyverse)
library(rio)


# ============================================================
# 02. Pin dplyr verbs
# ============================================================
# Prevent conflicts if plyr or other packages mask dplyr verbs

mutate <- dplyr::mutate
summarise <- dplyr::summarise
summarize <- dplyr::summarize
arrange <- dplyr::arrange
rename <- dplyr::rename
count <- dplyr::count
desc <- dplyr::desc
select <- dplyr::select
filter <- dplyr::filter


# ============================================================
# 03. Define directories
# ============================================================

meta_dir <- "results/tables/meta-analysis"
figure_dir <- "results/figures/meta-analysis"

dir.create(
  meta_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  figure_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 04. Load Random-Effects REML results
# ============================================================
# IMPORTANT:
# random_effect_model.csv contains:
#
# Gene_ID
# randomSummary
# randomP
# heterogeneity statistics
# error
#
# Gene_Symbol is NOT directly present in this file.
# Therefore annotation is added separately below.


meta_result <- import(
  file.path(
    meta_dir,
    "random_effect_model.csv"
  )
) |>
  select(
    Gene_ID,
    randomSummary,
    randomP
  ) |>
  rename(
    log2FC = randomSummary,
    P.Value = randomP
  ) |>
  mutate(
    Gene_ID = as.character(Gene_ID),
    log2FC = as.numeric(log2FC),
    P.Value = as.numeric(P.Value)
  )


# ============================================================
# 05. Multiple-testing correction
# ============================================================
# BH correction across all 10,128 genes

meta_result <- meta_result |>
  mutate(
    FDR = p.adjust(
      P.Value,
      method = "BH"
    )
  )


# Check number of genes

cat(
  "\nTotal genes analyzed: ",
  nrow(meta_result),
  "\n"
)


# ============================================================
# 06. Load gene annotation
# ============================================================

anno_dir <- "results/tables/annotated"

if (dir.exists(anno_dir)) {
  
  anno_files <- list.files(
    anno_dir,
    pattern = "[.]csv$",
    full.names = TRUE
  )
  
  if (length(anno_files) > 0) {
    
    gene_annotation <- anno_files |>
      lapply(function(f) {
        
        dat <- import(f)
        
        dat |>
          select(
            any_of(
              c(
                "Gene_ID",
                "Gene_Symbol",
                "Gene_Description"
              )
            )
          )
        
      }) |>
      bind_rows() |>
      mutate(
        Gene_ID = as.character(Gene_ID)
      ) |>
      distinct(
        Gene_ID,
        .keep_all = TRUE
      )
    
  } else {
    
    gene_annotation <- tibble(
      Gene_ID = character(),
      Gene_Symbol = character(),
      Gene_Description = character()
    )
    
  }
  
} else {
  
  gene_annotation <- tibble(
    Gene_ID = character(),
    Gene_Symbol = character(),
    Gene_Description = character()
  )
  
}


# ============================================================
# 07. Add annotation to meta-analysis results
# ============================================================

meta_result <- meta_result |>
  left_join(
    gene_annotation,
    by = "Gene_ID"
  )


# ============================================================
# 08. Significance classification
# ============================================================
# Final criteria:
#
#   FDR < 0.05
#   AND
#   |log2FC| >= 1
#
# Up:
#   log2FC >= 1
#
# Down:
#   log2FC <= -1


meta_result <- meta_result |>
  mutate(
    Significance = case_when(
      
      FDR < 0.05 &
        log2FC >= 1 ~ "Up",
      
      FDR < 0.05 &
        log2FC <= -1 ~ "Down",
      
      TRUE ~ "NS"
      
    )
  )


# ============================================================
# 09. Up / Down regulation summary
# ============================================================

gene_stats <- meta_result |>
  filter(
    Significance != "NS"
  ) |>
  count(
    Significance
  ) |>
  mutate(
    Percentage =
      n / nrow(meta_result) * 100
  )


print(gene_stats)


# Save summary

export(
  gene_stats,
  file.path(
    meta_dir,
    "meta_degs_summary_stats.csv"
  )
)


# ============================================================
# 10. Annotated significant DEGs only
# ============================================================

annotated_genes <- meta_result |>
  mutate(
    Gene_Symbol = na_if(
      Gene_Symbol,
      ""
    )
  ) |>
  filter(
    !is.na(Gene_Symbol),
    Significance != "NS"
  )


export(
  annotated_genes,
  file.path(
    meta_dir,
    "filtered_meta_degs_annotated_only.csv"
  )
)


# ============================================================
# 11. Final significant meta-DEG regulation table
# ============================================================

meta_key_results <- meta_result |>
  filter(
    Significance != "NS"
  ) |>
  select(
    Gene_ID,
    Gene_Symbol,
    Gene_Description,
    log2FoldChange = log2FC,
    P.Value,
    FDR,
    Regulation = Significance
  )


export(
  meta_key_results,
  file.path(
    meta_dir,
    "meta_degs_regulation.csv"
  )
)


# ============================================================
# 12. Save complete meta-analysis results with FDR
# ============================================================
# This file contains all 10,128 genes.

export(
  meta_result,
  file.path(
    meta_dir,
    "random_effect_model_with_FDR.csv"
  )
)


# ============================================================
# 13. Print final summary
# ============================================================

total_genes <- nrow(meta_result)

up_genes <- sum(
  meta_result$Significance == "Up",
  na.rm = TRUE
)

down_genes <- sum(
  meta_result$Significance == "Down",
  na.rm = TRUE
)

significant_genes <- up_genes + down_genes


cat("\n")
cat("====================================================\n")
cat(" RANDOM-EFFECTS REML META-ANALYSIS SUMMARY\n")
cat("====================================================\n")

cat(
  "Total genes tested       : ",
  total_genes,
  "\n"
)

cat(
  "Significant meta-DEGs    : ",
  significant_genes,
  "\n"
)

cat(
  "Up-regulated genes       : ",
  up_genes,
  "\n"
)

cat(
  "Down-regulated genes     : ",
  down_genes,
  "\n"
)

cat(
  "Criteria                 : FDR < 0.05 & |log2FC| >= 1\n"
)

cat("====================================================\n")


# ============================================================
# 14. Meta-analysis volcano plot data
# ============================================================

volcano_data <- meta_result |>
  filter(
    !is.na(FDR),
    FDR > 0,
    is.finite(FDR),
    is.finite(log2FC)
  ) |>
  mutate(
    neg_log10_FDR = -log10(FDR)
  )


# ============================================================
# 15. Meta-analysis volcano plot
# ============================================================

volcano_plot <- ggplot(
  volcano_data,
  aes(
    x = log2FC,
    y = neg_log10_FDR
  )
) +
  
  geom_point(
    aes(
      color = Significance
    ),
    alpha = 0.6,
    size = 1.5
  ) +
  
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  
  labs(
    title = "Meta-Analysis Volcano Plot",
    
    subtitle =
      "Random-Effects REML Meta-Analysis",
    
    x =
      "Pooled Random-Effects log2 Fold Change",
    
    y =
      "-log10(FDR)",
    
    color =
      "Significance"
  ) +
  
  theme_bw()


# Display plot

print(volcano_plot)


# ============================================================
# 16. Save volcano plot
# ============================================================

ggsave(
  filename = file.path(
    figure_dir,
    "meta_analysis_volcano.pdf"
  ),
  plot = volcano_plot,
  width = 8,
  height = 6
)


ggsave(
  filename = file.path(
    figure_dir,
    "meta_analysis_volcano.png"
  ),
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================
# 17. Cross-reference Venn intersection genes
# ============================================================

venn_common <-
  "results/tables/venn/common_genes_all_datasets.csv"

comb_file <-
  file.path(
    meta_dir,
    "meta_combining_mean.csv"
  )


if (
  file.exists(venn_common) &&
  file.exists(comb_file)
) {
  
  intersect_ids <- import(
    venn_common
  ) |>
    pull(Gene_ID) |>
    as.character()
  
  
  combined_meta <- import(
    comb_file
  ) |>
    mutate(
      Gene_ID =
        as.character(Gene_ID)
    )
  
  
  intersects <- combined_meta |>
    filter(
      Gene_ID %in% intersect_ids
    )
  
  
  export(
    intersects,
    "results/tables/venn/intersect_genes_meta_analysis.csv"
  )
  
  
  message(
    "Venn-intersect genes with meta stats: ",
    nrow(intersects)
  )
  
  
} else {
  
  message(
    "Missing Venn or combining-approach file; ",
    "skipping cross-reference."
  )
  
}


# ============================================================
# 18. Final message
# ============================================================

cat("\n")
cat("Meta-analysis summary completed successfully.\n")
cat(
  "Volcano plot saved in: ",
  figure_dir,
  "\n"
)
cat(
  "Meta-DEG tables saved in: ",
  meta_dir,
  "\n"
)
source("scripts/06_meta summary.R")

volcano_plot
file.exists("results/figures/meta-analysis/meta_analysis_volcano.png")
list.files(
  "results/figures/meta-analysis",
  full.names = TRUE
)