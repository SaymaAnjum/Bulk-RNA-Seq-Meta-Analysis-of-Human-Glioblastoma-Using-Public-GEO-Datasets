# ============================================================
# Functional Enrichment Analysis
# ORA + GSEA
# Glioblastoma Meta-analysis
# ============================================================

# ------------------------------------------------------------
# Load packages
# ------------------------------------------------------------

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ReactomePA)
library(msigdbr)
library(tidyverse)
library(rio)
library(AnnotationDbi)

# ------------------------------------------------------------
# Pin dplyr functions
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
# Thresholds
# ------------------------------------------------------------

PADJ_CUTOFF <- 0.05
LFC_CUTOFF <- 1

# ------------------------------------------------------------
# Input / Output
# ------------------------------------------------------------

META_FILE <- "results/tables/meta-analysis/random_effect_model.csv"

OUT_CSV <- "results/tables/enrichment"
OUT_FIG <- "results/figures/enrichment"

PREFIX <- "glioblastoma"

dir.create(
  OUT_CSV,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  OUT_FIG,
  showWarnings = FALSE,
  recursive = TRUE
)

# ============================================================
# Ensembl -> Entrez conversion
# ============================================================

convert_ensembl_to_entrez <- function(gene_ids) {
  
  gene_ids <- as.character(gene_ids)
  
  gene_ids <- gene_ids[
    !is.na(gene_ids) &
      gene_ids != ""
  ]
  
  gene_ids <- unique(gene_ids)
  
  mapping <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = gene_ids,
    keytype = "ENSEMBL",
    columns = c(
      "ENSEMBL",
      "ENTREZID"
    )
  )
  
  mapping <- mapping |>
    dplyr::filter(
      !is.na(ENTREZID),
      ENTREZID != ""
    ) |>
    dplyr::distinct(
      ENSEMBL,
      .keep_all = TRUE
    )
  
  return(mapping)
}

# ============================================================
# Build ranked list for GSEA
# ============================================================

build_ranked_list <- function(meta_df) {
  
  df <- meta_df |>
    dplyr::filter(
      !is.na(randomP),
      !is.na(randomSummary),
      randomP > 0,
      !is.na(ENTREZID)
    ) |>
    dplyr::mutate(
      rank_score =
        sign(randomSummary) *
        -log10(randomP),
      ENTREZID = as.character(ENTREZID)
    ) |>
    dplyr::group_by(ENTREZID) |>
    dplyr::slice_max(
      order_by = abs(rank_score),
      n = 1,
      with_ties = FALSE
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(
      dplyr::desc(rank_score)
    )
  
  ranked <- df$rank_score
  
  names(ranked) <- df$ENTREZID
  
  return(ranked)
}

# ============================================================
# ORA
# ============================================================

run_ora <- function(
    entrez_ids,
    universe_entrez
) {
  
  if (length(entrez_ids) < 5) {
    
    message(
      "    Too few genes for ORA (n=",
      length(entrez_ids),
      "), skipping."
    )
    
    return(NULL)
  }
  
  results <- list()
  
  # ----------------------------------------------------------
  # GO
  # ----------------------------------------------------------
  
  for (ont in c("BP", "MF", "CC")) {
    
    results[[paste0("GO_", ont)]] <-
      tryCatch(
        
        enrichGO(
          gene = entrez_ids,
          universe = universe_entrez,
          OrgDb = org.Hs.eg.db,
          keyType = "ENTREZID",
          ont = ont,
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.20,
          readable = TRUE
        ),
        
        error = function(e) {
          
          message(
            "    GO-",
            ont,
            " error: ",
            e$message
          )
          
          NULL
        }
      )
  }
  
  # ----------------------------------------------------------
  # KEGG
  # ----------------------------------------------------------
  
  results[["KEGG"]] <-
    tryCatch(
      
      enrichKEGG(
        gene = entrez_ids,
        universe = universe_entrez,
        organism = "hsa",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05
      ),
      
      error = function(e) {
        
        message(
          "    KEGG error: ",
          e$message
        )
        
        NULL
      }
    )
  
  # ----------------------------------------------------------
  # Reactome
  # ----------------------------------------------------------
  
  results[["Reactome"]] <-
    tryCatch(
      
      enrichPathway(
        gene = entrez_ids,
        universe = universe_entrez,
        organism = "human",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        readable = TRUE
      ),
      
      error = function(e) {
        
        message(
          "    Reactome error: ",
          e$message
        )
        
        NULL
      }
    )
  
  return(results)
}

# ============================================================
# GSEA
# ============================================================

run_gsea <- function(ranked_list) {
  
  if (length(ranked_list) < 100) {
    
    message(
      "    Too few ranked genes (n=",
      length(ranked_list),
      "), skipping GSEA."
    )
    
    return(NULL)
  }
  
  results <- list()
  
  # ----------------------------------------------------------
  # GO Biological Process
  # ----------------------------------------------------------
  
  results[["GSEA_GO_BP"]] <-
    tryCatch(
      
      gseGO(
        geneList = ranked_list,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        verbose = FALSE
      ),
      
      error = function(e) {
        
        message(
          "    GSEA GO error: ",
          e$message
        )
        
        NULL
      }
    )
  
  # ----------------------------------------------------------
  # KEGG
  # ----------------------------------------------------------
  
  results[["GSEA_KEGG"]] <-
    tryCatch(
      
      gseKEGG(
        geneList = ranked_list,
        organism = "hsa",
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        verbose = FALSE
      ),
      
      error = function(e) {
        
        message(
          "    GSEA KEGG error: ",
          e$message
        )
        
        NULL
      }
    )
  
  # ----------------------------------------------------------
  # Hallmark
  # ----------------------------------------------------------
  
  hallmark <- msigdbr(
    species = "Homo sapiens",
    collection = "H"
  ) |>
    dplyr::select(
      gs_name,
      ncbi_gene
    ) |>
    dplyr::mutate(
      ncbi_gene = as.character(ncbi_gene)
    )
  
  results[["GSEA_Hallmark"]] <-
    tryCatch(
      
      GSEA(
        geneList = ranked_list,
        TERM2GENE = hallmark,
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        verbose = FALSE
      ),
      
      error = function(e) {
        
        message(
          "    GSEA Hallmark error: ",
          e$message
        )
        
        NULL
      }
    )
  
  return(results)
}

# ============================================================
# Export results
# ============================================================

export_enrichment <- function(
    result_list,
    out_dir,
    prefix
) {
  
  if (is.null(result_list)) {
    return(invisible(NULL))
  }
  
  dir.create(
    out_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  for (name in names(result_list)) {
    
    res <- result_list[[name]]
    
    if (is.null(res)) {
      next
    }
    
    df <- tryCatch(
      as.data.frame(res),
      error = function(e) NULL
    )
    
    if (is.null(df) || nrow(df) == 0) {
      next
    }
    
    export(
      df,
      file.path(
        out_dir,
        paste0(
          prefix,
          "_",
          name,
          ".csv"
        )
      )
    )
  }
}

# ============================================================
# Dotplot
# ============================================================

save_dotplot <- function(
    res,
    title,
    path,
    width = 10,
    height = 8,
    show = 20
) {
  
  if (is.null(res)) {
    return(invisible(NULL))
  }
  
  df <- tryCatch(
    as.data.frame(res),
    error = function(e) NULL
  )
  
  if (is.null(df) || nrow(df) == 0) {
    return(invisible(NULL))
  }
  
  is_gsea <- all(
    c("NES", "setSize") %in%
      colnames(df)
  )
  
  if (is_gsea) {
    
    df <- df[
      order(df$p.adjust),
      ,
      drop = FALSE
    ]
    
    df <- head(
      df,
      show
    )
    
    df$Label <-
      stringr::str_wrap(
        df$Description,
        34
      )
    
    p <- ggplot(
      df,
      aes(
        x = NES,
        y = reorder(
          Label,
          NES
        )
      )
    ) +
      
      geom_vline(
        xintercept = 0,
        color = "grey70",
        linewidth = 0.4
      ) +
      
      geom_segment(
        aes(
          x = 0,
          xend = NES,
          yend = reorder(
            Label,
            NES
          )
        ),
        color = "grey80",
        linewidth = 0.5
      ) +
      
      geom_point(
        aes(
          size = setSize,
          fill = p.adjust
        ),
        shape = 21,
        color = "grey30"
      ) +
      
      scale_fill_gradient(
        low = "#b2182b",
        high = "#2166ac",
        trans = "log10",
        name = "p.adjust"
      ) +
      
      scale_size_continuous(
        range = c(3, 10),
        name = "Set size"
      ) +
      
      labs(
        title = title,
        x = "Normalized Enrichment Score (NES)",
        y = NULL
      ) +
      
      theme_bw(
        base_size = 13
      ) +
      
      theme(
        plot.title = element_text(
          face = "bold",
          hjust = 0.5
        ),
        axis.text.y = element_text(
          color = "black"
        ),
        panel.grid.minor =
          element_blank()
      )
    
  } else {
    
    p <-
      dotplot(
        res,
        showCategory = show
      ) +
      ggtitle(title)
  }
  
  ggsave(
    path,
    plot = p,
    width = width,
    height = height,
    dpi = 600,
    bg = "white"
  )
  
  invisible(p)
}

# ============================================================
# GSEA enrichment plot
# ============================================================

save_gsea_plot <- function(
    res,
    title,
    path,
    n = 3
) {
  
  if (is.null(res)) {
    return(invisible(NULL))
  }
  
  df <- tryCatch(
    as.data.frame(res),
    error = function(e) NULL
  )
  
  if (is.null(df) || nrow(df) == 0) {
    return(invisible(NULL))
  }
  
  ids <- head(
    df$ID,
    n
  )
  
  p <- gseaplot2(
    res,
    geneSetID = ids,
    title = title
  )
  
  ggsave(
    path,
    plot = p,
    width = 12,
    height = 9,
    dpi = 300,
    bg = "white"
  )
  
  invisible(p)
}

# ============================================================
# DRIVER
# ============================================================

message(
  "\n=== ",
  toupper(PREFIX),
  " ENRICHMENT ==="
)

# ------------------------------------------------------------
# Load meta-analysis results
# ------------------------------------------------------------

meta <- import(
  META_FILE
) |>
  dplyr::mutate(
    Gene_ID = as.character(
      Gene_ID
    )
  )

message(
  "  Genes in meta: ",
  nrow(meta)
)

# ------------------------------------------------------------
# Required columns
# ------------------------------------------------------------

required_columns <- c(
  "Gene_ID",
  "randomSummary",
  "randomP"
)

missing_columns <-
  setdiff(
    required_columns,
    colnames(meta)
  )

if (length(missing_columns) > 0) {
  
  stop(
    "Missing required columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------
# BH adjustment
# ------------------------------------------------------------

meta <- meta |>
  dplyr::mutate(
    rem_padj =
      p.adjust(
        randomP,
        method = "BH"
      )
  )

# ------------------------------------------------------------
# Ensembl -> Entrez
# ------------------------------------------------------------

message(
  "  Converting Ensembl IDs to Entrez IDs..."
)

gene_mapping <-
  convert_ensembl_to_entrez(
    meta$Gene_ID
  )

message(
  "  Successfully mapped genes: ",
  nrow(gene_mapping)
)

# ------------------------------------------------------------
# Add Entrez ID
# ------------------------------------------------------------

meta <- meta |>
  dplyr::left_join(
    gene_mapping |>
      dplyr::select(
        ENSEMBL,
        ENTREZID
      ),
    by = c(
      "Gene_ID" = "ENSEMBL"
    )
  )

message(
  "  Genes with Entrez ID: ",
  sum(
    !is.na(
      meta$ENTREZID
    )
  )
)

# ------------------------------------------------------------
# Universe
# ------------------------------------------------------------

universe_entrez <-
  meta |>
  dplyr::filter(
    !is.na(ENTREZID)
  ) |>
  dplyr::pull(
    ENTREZID
  ) |>
  unique()

message(
  "  Universe genes for ORA: ",
  length(universe_entrez)
)

# ============================================================
# Significant DEGs
# ============================================================

sig <- meta |>
  dplyr::filter(
    rem_padj < PADJ_CUTOFF,
    abs(randomSummary) >= LFC_CUTOFF,
    !is.na(ENTREZID)
  )

message(
  "  Significant DEGs (padj < ",
  PADJ_CUTOFF,
  ", |log2FC| >= ",
  LFC_CUTOFF,
  "): ",
  nrow(sig)
)

# ------------------------------------------------------------
# Up / Down
# ------------------------------------------------------------

entrez_up <-
  sig |>
  dplyr::filter(
    randomSummary > 0
  ) |>
  dplyr::pull(
    ENTREZID
  ) |>
  unique()

entrez_down <-
  sig |>
  dplyr::filter(
    randomSummary < 0
  ) |>
  dplyr::pull(
    ENTREZID
  ) |>
  unique()

entrez_all <-
  sig |>
  dplyr::pull(
    ENTREZID
  ) |>
  unique()

message(
  "  Up: ",
  length(entrez_up),
  "  Down: ",
  length(entrez_down)
)

# ============================================================
# ORA
# ============================================================

message(
  "  Running ORA..."
)

ora_up <-
  run_ora(
    entrez_up,
    universe_entrez
  )

ora_down <-
  run_ora(
    entrez_down,
    universe_entrez
  )

ora_all <-
  run_ora(
    entrez_all,
    universe_entrez
  )

# ------------------------------------------------------------
# Export ORA
# ------------------------------------------------------------

export_enrichment(
  ora_up,
  OUT_CSV,
  paste0(
    PREFIX,
    "_ora_up"
  )
)

export_enrichment(
  ora_down,
  OUT_CSV,
  paste0(
    PREFIX,
    "_ora_down"
  )
)

export_enrichment(
  ora_all,
  OUT_CSV,
  paste0(
    PREFIX,
    "_ora_all"
  )
)

# ------------------------------------------------------------
# ORA plots
# ------------------------------------------------------------

if (!is.null(ora_all)) {
  
  for (name in names(ora_all)) {
    
    save_dotplot(
      ora_all[[name]],
      paste(
        PREFIX,
        name
      ),
      file.path(
        OUT_FIG,
        paste0(
          PREFIX,
          "_ora_all_",
          name,
          "_dotplot.png"
        )
      )
    )
  }
}

# ============================================================
# GSEA
# ============================================================

message(
  "  Running GSEA..."
)

ranked <-
  build_ranked_list(
    meta
  )

message(
  "  Ranked genes for GSEA: ",
  length(ranked)
)

gsea_res <-
  run_gsea(
    ranked
  )

# ------------------------------------------------------------
# Export GSEA
# ------------------------------------------------------------

export_enrichment(
  gsea_res,
  OUT_CSV,
  paste0(
    PREFIX,
    "_gsea"
  )
)

# ------------------------------------------------------------
# GSEA plots
# ------------------------------------------------------------

if (!is.null(gsea_res)) {
  
  for (name in names(gsea_res)) {
    
    save_dotplot(
      gsea_res[[name]],
      paste(
        PREFIX,
        name
      ),
      file.path(
        OUT_FIG,
        paste0(
          PREFIX,
          "_",
          name,
          "_dotplot.png"
        )
      )
    )
    
    save_gsea_plot(
      gsea_res[[name]],
      paste(
        PREFIX,
        name,
        "- Top Pathways"
      ),
      file.path(
        OUT_FIG,
        paste0(
          PREFIX,
          "_",
          name,
          "_enrichplot.png"
        )
      )
    )
  }
}

# ============================================================
# Final summary
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "Functional enrichment analysis completed successfully.\n"
)

cat(
  "============================================================\n"
)

cat(
  "Meta-analysis genes: ",
  nrow(meta),
  "\n"
)

cat(
  "Genes mapped to Entrez: ",
  sum(
    !is.na(
      meta$ENTREZID
    )
  ),
  "\n"
)

cat(
  "Significant DEGs: ",
  nrow(sig),
  "\n"
)

cat(
  "Up-regulated: ",
  length(entrez_up),
  "\n"
)

cat(
  "Down-regulated: ",
  length(entrez_down),
  "\n"
)

cat(
  "\nEnrichment tables saved to:\n",
  OUT_CSV,
  "\n"
)

cat(
  "\nEnrichment figures saved to:\n",
  OUT_FIG,
  "\n"
)

cat(
  "============================================================\n"
)

parse("scripts/08_enrichment_analysis.R")
file.edit("scripts/08_enrichment_analysis.R")
file.exists("scripts/08_enrichment_analysis.R")
file.exists("scripts/08_enrichment_analysis.R")
source("scripts/08_enrichment_analysis.R")
list.files(
  "results/tables/enrichment",
  full.names = TRUE
)

list.files(
  "results/figures/enrichment",
  full.names = TRUE
)