# ============================================================
# 04_common_genes_venn.R
# Common DEGs across Glioblastoma datasets
# ============================================================


# ------------------------------------------------------------
# Step 1: Load required packages
# ------------------------------------------------------------

library(tidyverse)
library(VennDiagram)
library(ggVennDiagram)
library(rio)


# ------------------------------------------------------------
# Step 2: Pin dplyr verbs
# ------------------------------------------------------------

mutate <- dplyr::mutate
summarise <- dplyr::summarise
summarize <- dplyr::summarize
arrange <- dplyr::arrange
rename <- dplyr::rename
count <- dplyr::count
desc <- dplyr::desc
select <- dplyr::select
filter <- dplyr::filter


# ------------------------------------------------------------
# Step 3: Create output directories
# ------------------------------------------------------------

dir.create(
  "results/figures/venn",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "results/tables/venn",
  showWarnings = FALSE,
  recursive = TRUE
)


# ------------------------------------------------------------
# Step 4: Define GEO datasets
# ------------------------------------------------------------

geo_ids <- c(
  "GSE147352",
  "GSE151352"
)


# ------------------------------------------------------------
# Step 5: Locate annotated DEG tables
# ------------------------------------------------------------

input_dir <- "results/tables/annotated"

deg_files <- file.path(
  input_dir,
  paste0(
    geo_ids,
    ".csv"
  )
)

names(deg_files) <- geo_ids


# ------------------------------------------------------------
# Step 6: Check that files exist
# ------------------------------------------------------------

missing_files <- deg_files[
  !file.exists(deg_files)
]

if (length(missing_files) > 0) {
  
  stop(
    "The following annotation files are missing:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}


# ------------------------------------------------------------
# Step 7: Define DEG thresholds
# ------------------------------------------------------------

sig_cutoff <- 0.05
lfc_cutoff <- 1


# ------------------------------------------------------------
# Step 8: Load and filter significant DEGs
# ------------------------------------------------------------

gene_sets <- map(
  deg_files,
  function(f) {
    
    import(f) |>
      filter(
        !is.na(padj),
        padj < sig_cutoff,
        abs(log2FoldChange) >= lfc_cutoff
      ) |>
      pull(Gene_ID) |>
      unique() |>
      as.character()
  }
)


# ------------------------------------------------------------
# Step 9: Name the gene sets
# ------------------------------------------------------------

names(gene_sets) <- geo_ids


# ------------------------------------------------------------
# Step 10: Print number of significant DEGs
# ------------------------------------------------------------

message("Number of significant DEGs:")

for (geo_id in names(gene_sets)) {
  
  message(
    geo_id,
    ": ",
    length(gene_sets[[geo_id]])
  )
}


# ------------------------------------------------------------
# Step 11: Remove empty gene sets
# ------------------------------------------------------------

gene_sets <- gene_sets[
  map_lgl(
    gene_sets,
    ~ length(.x) > 0
  )
]

n_sets <- length(gene_sets)


# ------------------------------------------------------------
# Step 12: Check number of datasets
# ------------------------------------------------------------

if (n_sets < 2) {
  
  stop(
    "At least two datasets with significant DEGs are required."
  )
}


# ------------------------------------------------------------
# Step 13: Create Venn diagram
# ------------------------------------------------------------

if (n_sets >= 2 && n_sets <= 5) {
  
  fill_cols <- c(
    "#99d8c9",
    "#addd8e",
    "#bcbddc",
    "#fec44f",
    "#fc9272"
  )[seq_len(n_sets)]
  
  cat_cols <- c(
    "#2ca25f",
    "#31a354",
    "#756bb1",
    "#e6550d",
    "#de2d26"
  )[seq_len(n_sets)]
  
  
  png(
    "results/figures/venn/venn_diagram.png",
    units = "in",
    width = 7,
    height = 7,
    res = 600
  )
  
  
  venn.plot <- venn.diagram(
    x = gene_sets,
    filename = NULL,
    fill = fill_cols,
    alpha = 0.6,
    cex = 0.8,
    cat.cex = 1,
    cat.col = cat_cols,
    margin = 0.12,
    lwd = 2,
    disable.logging = TRUE
  )
  
  
  grid::grid.draw(
    venn.plot
  )
  
  dev.off()
  
  
} else if (n_sets > 5) {
  
  upset <- ggVennDiagram::ggVennDiagram(
    gene_sets,
    force_upset = TRUE
  )
  
  ggsave(
    "results/figures/venn/upset_plot.png",
    upset,
    width = 12,
    height = 7,
    units = "in",
    bg = "white",
    dpi = 600
  )
}


# ------------------------------------------------------------
# Step 14: Calculate all intersection regions
# ------------------------------------------------------------

regions <- ggVennDiagram::process_region_data(
  ggVennDiagram::Venn(
    gene_sets
  )
)


# ------------------------------------------------------------
# Step 15: Save intersection summary
# ------------------------------------------------------------

regions_out <- regions |>
  select(
    id,
    name,
    count
  ) |>
  arrange(
    desc(count)
  )


export(
  as.data.frame(regions_out),
  "results/tables/venn/intersection_summary.csv"
)


# ------------------------------------------------------------
# Step 16: Find genes shared by ALL datasets
# ------------------------------------------------------------

common_genes <- Reduce(
  intersect,
  gene_sets
)


# ------------------------------------------------------------
# Step 17: Save common genes
# ------------------------------------------------------------

common_df <- data.frame(
  Gene_ID = as.character(
    common_genes
  )
)


export(
  common_df,
  "results/tables/venn/common_genes_all_datasets.csv"
)


# ------------------------------------------------------------
# Step 18: Create DEG frequency table
# ------------------------------------------------------------

deg_frequency <- gene_sets |>
  unlist() |>
  table() |>
  as.data.frame() |>
  setNames(
    c(
      "Gene_ID",
      "n_datasets"
    )
  ) |>
  arrange(
    desc(n_datasets)
  )


# ------------------------------------------------------------
# Step 19: Save DEG frequency table
# ------------------------------------------------------------

export(
  deg_frequency,
  "results/tables/venn/deg_frequency.csv"
)


# ------------------------------------------------------------
# Step 20: Print final summary
# ------------------------------------------------------------

message(
  "Total common DEGs across all datasets: ",
  length(common_genes)
)

message(
  "Finished common DEG analysis."
)