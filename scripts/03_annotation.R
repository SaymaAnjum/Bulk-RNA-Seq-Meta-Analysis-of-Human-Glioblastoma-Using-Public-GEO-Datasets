# ============================================================
# 03_annotation.R
# Gene annotation of DESeq2 results
# ============================================================


# ------------------------------------------------------------
# Step 1: Load packages
# ------------------------------------------------------------

library(genekitr)
library(tidyverse)
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
# Step 3: Create output directory
# ------------------------------------------------------------

dir.create(
  "results/tables/annotated",
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
# Step 5: Process each dataset
# ------------------------------------------------------------

for (geo_id in geo_ids) {
  
  message("Processing ", geo_id, " ...")
  
  
  # ----------------------------------------------------------
  # Step 6: Load DESeq2 results
  # ----------------------------------------------------------
  
  input_file <- paste0(
    "results/tables/DESeq2/",
    geo_id,
    ".csv"
  )
  
  
  if (!file.exists(input_file)) {
    
    message(
      "Skipping ",
      geo_id,
      ": DESeq2 output not found"
    )
    
    next
  }
  
  
  deseq_results <- import(
    input_file
  )
  
  
  # ----------------------------------------------------------
  # Step 7: Make Gene_ID character
  # ----------------------------------------------------------
  
  deseq_results <- deseq_results |>
    mutate(
      Gene_ID = as.character(Gene_ID)
    )
  
  
  # ----------------------------------------------------------
  # Step 8: Remove Ensembl version suffix
  # ----------------------------------------------------------
  
  deseq_results <- deseq_results |>
    mutate(
      Gene_ID = sub(
        "\\..*$",
        "",
        Gene_ID
      )
    )
  
  
  # ----------------------------------------------------------
  # Step 9: Get gene annotation
  # ----------------------------------------------------------
  
  gene_info <- tryCatch(
    
    genInfo(
      id = deseq_results$Gene_ID,
      org = "hs",
      unique = TRUE,
      keepNA = FALSE
    ),
    
    error = function(e) {
      
      message(
        "genInfo failed for ",
        geo_id,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  
  # ----------------------------------------------------------
  # Step 10: Check annotation
  # ----------------------------------------------------------
  
  if (is.null(gene_info)) {
    
    next
  }
  
  
  gene_info <- as.data.frame(
    gene_info
  )
  
  
  # ----------------------------------------------------------
  # Step 11: Join annotation with DESeq2 results
  # ----------------------------------------------------------
  
  annotated <- deseq_results |>
    left_join(
      gene_info,
      by = c(
        "Gene_ID" = "input_id"
      )
    )
  
  
  # ----------------------------------------------------------
  # Step 12: Rename annotation columns
  # ----------------------------------------------------------
  
  colnames(annotated)[
    colnames(annotated) == "symbol"
  ] <- "Gene_Symbol"
  
  
  colnames(annotated)[
    colnames(annotated) == "gene_name"
  ] <- "Gene_Description"
  
  
  colnames(annotated)[
    colnames(annotated) == "gene_biotype"
  ] <- "Gene_Biotype"
  
  
  # ----------------------------------------------------------
  # Step 13: Keep relevant columns
  # ----------------------------------------------------------
  
  annotated <- annotated |>
    select(
      any_of(
        c(
          "Gene_ID",
          "Gene_Symbol",
          "Gene_Description",
          "Gene_Biotype",
          "baseMean",
          "log2FoldChange",
          "lfcSE",
          "stat",
          "pvalue",
          "padj"
        )
      )
    )
  
  
  # ----------------------------------------------------------
  # Step 14: Keep protein-coding genes
  # ----------------------------------------------------------
  
  if (
    "Gene_Biotype" %in%
    colnames(annotated)
  ) {
    
    annotated <- annotated |>
      filter(
        Gene_Biotype ==
          "protein_coding"
      )
  }
  
  
  # ----------------------------------------------------------
  # Step 15: Export annotation results
  # ----------------------------------------------------------
  
  output_file <- paste0(
    "results/tables/annotated/",
    geo_id,
    ".csv"
  )
  
  
  export(
    annotated,
    output_file
  )
  
  
  # ----------------------------------------------------------
  # Step 16: Completion message
  # ----------------------------------------------------------
  
  message(
    "Finished annotation: ",
    geo_id
  )
}

list.files(
  "results/tables/annotated",
  full.names = TRUE
)

ann147 <- rio::import(
  "results/tables/annotated/GSE147352.csv"
)

dim(ann147)
head(ann147)


ann151 <- rio::import(
  "results/tables/annotated/GSE151352.csv"
)

dim(ann151)
head(ann151)

colnames(ann147)