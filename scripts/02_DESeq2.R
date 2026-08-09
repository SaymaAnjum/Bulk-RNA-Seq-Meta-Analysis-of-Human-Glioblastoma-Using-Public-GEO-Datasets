
# ============================================================
# 02_DESeq2.R
# DESeq2 + SVA + PCA
# Datasets: GSE147352 and GSE151352
# ============================================================


# ------------------------------------------------------------
# Step 1: Load packages
# ------------------------------------------------------------

library(tidyverse)
library(DESeq2)
library(sva)
library(rio)
library(tximport)
library(EnsDb.Hsapiens.v86)


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
  "results/tables/DESeq2",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "results/figures/PCA",
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
# Step 5: Define PCA colors
# ------------------------------------------------------------

custom_colors <- c(
  "Cancer" = "#de2d26",
  "Normal" = "#2171b5",
  "cancer" = "#de2d26",
  "normal" = "#2171b5",
  "Tumor" = "#de2d26",
  "tumor" = "#de2d26"
)


# ============================================================
# Process each GEO dataset
# ============================================================

for (geo_id in geo_ids) {
  
  
  # ----------------------------------------------------------
  # Step 6: Load metadata
  # ----------------------------------------------------------
  
  metadata <- import(
    paste0(
      "data/metadata/",
      geo_id,
      "_metadata.csv"
    )
  )
  
  colnames(metadata) <- tolower(colnames(metadata))
  
  
  # ----------------------------------------------------------
  # Step 7: Build Salmon quant.sf file paths
  # ----------------------------------------------------------
  
  salmon_dir <- paste0(
    "data/",
    geo_id,
    "/salmon"
  )
  
  files <- file.path(
    salmon_dir,
    metadata$sample,
    "quant.sf"
  )
  
  
  # Check that all Salmon files exist
  
  if (!all(file.exists(files))) {
    
    missing_files <- files[!file.exists(files)]
    
    stop(
      paste(
        "Missing Salmon quant.sf files:",
        paste(missing_files, collapse = "\n")
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Step 8: Create transcript-to-gene mapping
  # ----------------------------------------------------------
  
  tx <- transcripts(
    EnsDb.Hsapiens.v86,
    return.type = "DataFrame"
  )
  
  tx2gene <- as.data.frame(
    tx[, c("tx_id", "gene_id")]
  )
  
  
  # ----------------------------------------------------------
  # Step 9: Import Salmon data with tximport
  # ----------------------------------------------------------
  
  txi <- tximport(
    files,
    type = "salmon",
    tx2gene = tx2gene,
    ignoreTxVersion = TRUE,
    ignoreAfterBar = TRUE
  )
  
  
  # ----------------------------------------------------------
  # Step 10: Prepare sample information
  # ----------------------------------------------------------
  
  colData <- data.frame(
    condition = relevel(
      factor(
        ifelse(
          grepl(
            "normal|control|adjacent|healthy|non.?tumou?r",
            tolower(as.character(metadata$condition))
          ),
          "Normal",
          "Tumor"
        )
      ),
      ref = "Normal"
    ),
    row.names = metadata$sample
  )
  
  
  # ----------------------------------------------------------
  # Step 11: Check sample order
  # ----------------------------------------------------------
  
  if (!all(
    rownames(colData) == colnames(txi$counts)
  )) {
    
    stop(
      paste(
        "Sample order mismatch in",
        geo_id
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Step 12: Create DESeq2 object
  # ----------------------------------------------------------
  
  dds <- DESeqDataSetFromTximport(
    txi = txi,
    colData = colData,
    design = ~ condition
  )
  
  
  # ----------------------------------------------------------
  # Step 13: Filter low-count genes
  # ----------------------------------------------------------
  
  keep <- rowSums(
    counts(dds) >= 10
  ) >= round(
    nrow(metadata) / 2
  )
  
  dds <- dds[keep, ]
  
  
  # ----------------------------------------------------------
  # Step 14: Estimate surrogate variables
  # ----------------------------------------------------------
  
  mod <- model.matrix(
    ~ condition,
    data = colData(dds)
  )
  
  mod0 <- model.matrix(
    ~ 1,
    data = colData(dds)
  )
  
  svobj <- svaseq(
    counts(dds),
    mod,
    mod0
  )
  
  n.sv <- svobj$n.sv
  
  
  # ----------------------------------------------------------
  # Step 15: Add surrogate variables
  # ----------------------------------------------------------
  
  if (n.sv > 0) {
    
    for (i in seq_len(n.sv)) {
      
      colData(dds)[[paste0("SV", i)]] <-
        svobj$sv[, i]
    }
    
    
    # --------------------------------------------------------
    # Step 16: Update design formula
    # --------------------------------------------------------
    
    sv_terms <- paste0(
      "SV",
      seq_len(n.sv),
      collapse = " + "
    )
    
    design(dds) <- as.formula(
      paste(
        "~",
        sv_terms,
        "+ condition"
      )
    )
    
  } else {
    
    design(dds) <- ~ condition
    
  }
  
  
  # ----------------------------------------------------------
  # Step 17: Run DESeq2
  # ----------------------------------------------------------
  
  dds_sva <- DESeq(dds)
  
  
  # ----------------------------------------------------------
  # Step 18: Extract differential expression results
  # ----------------------------------------------------------
  
  res <- results(
    dds_sva,
    contrast = c(
      "condition",
      "Tumor",
      "Normal"
    )
  )
  
  
  # ----------------------------------------------------------
  # Step 19: Convert results to dataframe
  # ----------------------------------------------------------
  
  res_df <- res |>
    as.data.frame() |>
    rownames_to_column("Gene_ID")
  
  
  # ----------------------------------------------------------
  # Step 20: Export DESeq2 results
  # ----------------------------------------------------------
  
  output_file <- paste0(
    "results/tables/DESeq2/",
    geo_id,
    ".csv"
  )
  
  export(
    res_df,
    output_file
  )
  
  
  # ----------------------------------------------------------
  # Step 21: VST before and after SVA
  # ----------------------------------------------------------
  
  dds_original <- dds
  
  design(dds_original) <- ~ condition
  
  vsd_pre <- vst(
    dds_original,
    blind = TRUE
  )
  
  vsd_post <- vst(
    dds_sva,
    blind = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Step 22: PCA data
  # ----------------------------------------------------------
  
  pca_data_pre <- plotPCA(
    vsd_pre,
    intgroup = "condition",
    returnData = TRUE
  )
  
  pca_data_post <- plotPCA(
    vsd_post,
    intgroup = "condition",
    returnData = TRUE
  )
  
  
  # ----------------------------------------------------------
  # Step 23: Variance explained
  # ----------------------------------------------------------
  
  percent_pre <- round(
    100 * attr(
      pca_data_pre,
      "percentVar"
    )
  )
  
  percent_post <- round(
    100 * attr(
      pca_data_post,
      "percentVar"
    )
  )
  
  
  # ----------------------------------------------------------
  # Step 24: Pre-correction PCA
  # ----------------------------------------------------------
  
  pca_pre <- ggplot(
    pca_data_pre,
    aes(
      PC1,
      PC2,
      color = condition
    )
  ) +
    
    geom_point(size = 3) +
    
    xlab(
      paste0(
        "PC1: ",
        percent_pre[1],
        "% variance"
      )
    ) +
    
    ylab(
      paste0(
        "PC2: ",
        percent_pre[2],
        "% variance"
      )
    ) +
    
    ggtitle(
      paste0(
        geo_id,
        " (Pre-Correction)"
      )
    ) +
    
    scale_color_manual(
      values = custom_colors
    ) +
    
    theme_minimal() +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 10
      ),
      
      legend.title = element_text(
        face = "bold",
        size = 10
      ),
      
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 1
      ),
      
      legend.position = "right",
      
      aspect.ratio = 1
    )
  
  
  # ----------------------------------------------------------
  # Step 25: Post-correction PCA
  # ----------------------------------------------------------
  
  pca_post <- ggplot(
    pca_data_post,
    aes(
      PC1,
      PC2,
      color = condition
    )
  ) +
    
    geom_point(size = 3) +
    
    xlab(
      paste0(
        "PC1: ",
        percent_post[1],
        "% variance"
      )
    ) +
    
    ylab(
      paste0(
        "PC2: ",
        percent_post[2],
        "% variance"
      )
    ) +
    
    ggtitle(
      paste0(
        geo_id,
        " (Post-Correction)"
      )
    ) +
    
    scale_color_manual(
      values = custom_colors
    ) +
    
    theme_minimal() +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 10
      ),
      
      legend.title = element_text(
        face = "bold",
        size = 10
      ),
      
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 1
      ),
      
      legend.position = "right",
      
      aspect.ratio = 1
    )
  
  
  # ----------------------------------------------------------
  # Step 26: Save PCA plots
  # ----------------------------------------------------------
  
  pca_dir <- paste0(
    "results/figures/PCA/",
    geo_id
  )
  
  dir.create(
    pca_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  
  ggsave(
    filename = paste0(
      pca_dir,
      "/",
      geo_id,
      "_pre.png"
    ),
    
    plot = pca_pre,
    
    width = 5,
    height = 5,
    
    units = "in",
    
    bg = "white",
    
    dpi = 600
  )
  
  
  ggsave(
    filename = paste0(
      pca_dir,
      "/",
      geo_id,
      "_post.png"
    ),
    
    plot = pca_post,
    
    width = 5,
    height = 5,
    
    units = "in",
    
    bg = "white",
    
    dpi = 600
  )
  
  
  # ----------------------------------------------------------
  # Step 27: Print progress
  # ----------------------------------------------------------
  
  message(
    "Finished: ",
    geo_id
  )
  
}