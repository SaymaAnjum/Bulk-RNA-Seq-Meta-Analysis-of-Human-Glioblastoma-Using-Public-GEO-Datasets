# ============================================================
# GLIOBLASTOMA META-ANALYSIS PPI NETWORK ANALYSIS
# Ensembl Gene ID -> Gene Symbol -> STRING PPI
# ============================================================

library(httr)
library(igraph)
library(ggraph)
library(ggplot2)
library(tidyverse)
library(rio)
library(scales)
library(AnnotationDbi)
library(org.Hs.eg.db)

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
distinct <- dplyr::distinct
transmute <- dplyr::transmute

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------

PADJ_CUTOFF <- 0.05
LFC_CUTOFF <- 1

STRING_SCORE <- 700
STRING_SPECIES <- 9606
STRING_MAX_NODES <- 2000

HUB_QUANTILE <- 0.90
MIN_MODULE_SIZE <- 5

# ------------------------------------------------------------
# Input / Output
# ------------------------------------------------------------

META_FILE <- "results/tables/meta-analysis/random_effect_model.csv"

OUT_CSV <- "results/tables/network"
OUT_FIG <- "results/figures/network"

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
# 1. Convert Ensembl IDs to Gene Symbols
# ============================================================

convert_ensembl_to_symbol <- function(gene_ids) {
  
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
      "SYMBOL"
    )
  )
  
  mapping <- mapping |>
    dplyr::filter(
      !is.na(SYMBOL),
      SYMBOL != ""
    ) |>
    dplyr::distinct(
      ENSEMBL,
      .keep_all = TRUE
    )
  
  mapping
}

# ============================================================
# 2. STRING REST API
# ============================================================

string_api_network <- function(
    symbols,
    required_score = STRING_SCORE,
    species = STRING_SPECIES
) {
  
  symbols <- unique(
    symbols[
      !is.na(symbols) &
        symbols != ""
    ]
  )
  
  if (length(symbols) < 2) {
    return(NULL)
  }
  
  message(
    "  Sending ",
    length(symbols),
    " genes to STRING..."
  )
  
  response <- tryCatch(
    
    httr::POST(
      "https://string-db.org/api/tsv/network",
      
      body = list(
        identifiers = paste(
          symbols,
          collapse = "\n"
        ),
        species = species,
        required_score = required_score,
        caller_identity = "glioblastoma_meta_pipeline"
      ),
      
      encode = "form",
      
      httr::timeout(120)
    ),
    
    error = function(e) {
      
      message(
        "  STRING API error: ",
        e$message
      )
      
      NULL
    }
  )
  
  if (is.null(response)) {
    return(NULL)
  }
  
  if (httr::http_error(response)) {
    
    message(
      "  STRING API returned HTTP error."
    )
    
    return(NULL)
  }
  
  text_result <- httr::content(
    response,
    as = "text",
    encoding = "UTF-8"
  )
  
  network <- tryCatch(
    
    utils::read.delim(
      text = text_result,
      stringsAsFactors = FALSE
    ),
    
    error = function(e) NULL
  )
  
  if (is.null(network)) {
    return(NULL)
  }
  
  required_columns <- c(
    "preferredName_A",
    "preferredName_B",
    "score"
  )
  
  if (!all(required_columns %in% colnames(network))) {
    
    message(
      "  STRING returned no usable interaction table."
    )
    
    return(NULL)
  }
  
  network
}

# ============================================================
# 3. Build igraph network
# ============================================================

build_ppi_graph <- function(
    network_df,
    lfc_map
) {
  
  if (
    is.null(network_df) ||
    nrow(network_df) == 0
  ) {
    return(NULL)
  }
  
  edges <- network_df |>
    
    dplyr::transmute(
      gene_from = preferredName_A,
      gene_to = preferredName_B,
      combined_score = score
    ) |>
    
    dplyr::filter(
      gene_from %in% names(lfc_map),
      gene_to %in% names(lfc_map)
    )
  
  if (nrow(edges) == 0) {
    return(NULL)
  }
  
  nodes <- unique(
    c(
      edges$gene_from,
      edges$gene_to
    )
  )
  
  vertex_df <- tibble(
    name = nodes,
    pooled_log2FC = unname(
      lfc_map[nodes]
    )
  )
  
  graph <- igraph::graph_from_data_frame(
    d = edges,
    vertices = vertex_df,
    directed = FALSE
  )
  
  graph <- igraph::simplify(
    graph,
    remove.multiple = TRUE,
    remove.loops = TRUE
  )
  
  graph
}

# ============================================================
# 4. Network metrics
# ============================================================

compute_metrics <- function(graph) {
  
  tibble(
    Symbol = igraph::V(graph)$name,
    
    Degree = igraph::degree(graph),
    
    Betweenness =
      igraph::betweenness(
        graph,
        normalized = TRUE
      ),
    
    Closeness =
      igraph::closeness(
        graph,
        normalized = TRUE
      ),
    
    Eigenvector =
      igraph::eigen_centrality(
        graph
      )$vector,
    
    Log2FC =
      igraph::V(graph)$pooled_log2FC
    
  ) |>
    
    dplyr::arrange(
      dplyr::desc(Degree)
    )
}

# ============================================================
# 5. Hub genes
# ============================================================

get_hub_genes <- function(
    metrics_df,
    quantile_cut = HUB_QUANTILE
) {
  
  threshold <- quantile(
    metrics_df$Degree,
    quantile_cut,
    na.rm = TRUE
  )
  
  metrics_df |>
    dplyr::filter(
      Degree >= threshold
    )
}

# ============================================================
# 6. Community detection
# ============================================================

detect_communities <- function(graph) {
  
  communities <- igraph::cluster_louvain(
    graph
  )
  
  igraph::V(graph)$community <-
    igraph::membership(
      communities
    )
  
  graph
}

# ============================================================
# 7. Network plot
# ============================================================

save_network_plot <- function(
    graph,
    title,
    path
) {
  
  if (
    is.null(graph) ||
    igraph::vcount(graph) < 2
  ) {
    return(invisible(NULL))
  }
  
  png(
    path,
    width = 12,
    height = 12,
    units = "in",
    res = 300
  )
  
  set.seed(42)
  
  degree_values <-
    igraph::degree(graph)
  
  vertex_sizes <-
    scales::rescale(
      degree_values,
      to = c(4, 16)
    )
  
  direction <-
    ifelse(
      igraph::V(graph)$pooled_log2FC >= 0,
      "Upregulated",
      "Downregulated"
    )
  
  vertex_colors <-
    ifelse(
      direction == "Upregulated",
      "#de2d26",
      "#2171b5"
    )
  
  top_nodes <-
    names(
      sort(
        degree_values,
        decreasing = TRUE
      )
    )
  
  top_nodes <-
    head(top_nodes, 15)
  
  labels <-
    ifelse(
      igraph::V(graph)$name %in% top_nodes,
      igraph::V(graph)$name,
      NA
    )
  
  plot(
    graph,
    
    layout = igraph::layout_with_fr(graph),
    
    vertex.size = vertex_sizes,
    
    vertex.color = vertex_colors,
    
    vertex.label = labels,
    
    vertex.label.cex = 0.7,
    
    vertex.label.color = "black",
    
    edge.width = 0.5,
    
    edge.color = "gray75",
    
    main = title
  )
  
  legend(
    "topright",
    
    legend = c(
      "Upregulated",
      "Downregulated"
    ),
    
    col = c(
      "#de2d26",
      "#2171b5"
    ),
    
    pch = 19,
    
    bty = "n"
  )
  
  dev.off()
}

# ============================================================
# 8. Hub gene barplot
# ============================================================

save_hub_barplot <- function(
    hub_df,
    path,
    n = 25
) {
  
  if (
    is.null(hub_df) ||
    nrow(hub_df) == 0
  ) {
    return(invisible(NULL))
  }
  
  df <- head(
    hub_df,
    n
  )
  
  p <- ggplot(
    df,
    aes(
      x = reorder(Symbol, Degree),
      y = Degree,
      fill = Log2FC
    )
  ) +
    
    geom_col() +
    
    scale_fill_gradient2(
      low = "#2171b5",
      mid = "white",
      high = "#de2d26",
      midpoint = 0,
      name = "Log2FC"
    ) +
    
    coord_flip() +
    
    labs(
      title = "Glioblastoma Top Hub Genes",
      x = NULL,
      y = "Degree Centrality"
    ) +
    
    theme_bw(
      base_size = 12
    ) +
    
    theme(
      plot.title =
        element_text(
          face = "bold",
          hjust = 0.5
        )
    )
  
  ggsave(
    path,
    plot = p,
    width = 9,
    height = 8,
    dpi = 300,
    bg = "white"
  )
}

# ============================================================
# 9. MAIN ANALYSIS
# ============================================================

message(
  "\n============================================================"
)

message(
  "\nGLIOBLASTOMA PPI NETWORK ANALYSIS"
)

message(
  "\n============================================================"
)

# ------------------------------------------------------------
# Read meta-analysis result
# ------------------------------------------------------------

meta <- rio::import(
  META_FILE
)

message(
  "\nGenes in meta-analysis: ",
  nrow(meta)
)

# ------------------------------------------------------------
# Check required columns
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

if (
  length(missing_columns) > 0
) {
  
  stop(
    "Missing required columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------
# BH adjusted p-value
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
# Ensembl -> Gene Symbol
# ------------------------------------------------------------

message(
  "\nConverting Ensembl IDs to Gene Symbols..."
)

gene_mapping <-
  convert_ensembl_to_symbol(
    meta$Gene_ID
  )

message(
  "Genes mapped to symbols: ",
  nrow(gene_mapping)
)

# ------------------------------------------------------------
# Add SYMBOL to meta-analysis result
# ------------------------------------------------------------

meta <- meta |>
  
  dplyr::left_join(
    
    gene_mapping |>
      
      dplyr::select(
        ENSEMBL,
        SYMBOL
      ),
    
    by = c(
      "Gene_ID" = "ENSEMBL"
    )
  )

message(
  "Genes with Gene Symbol: ",
  sum(
    !is.na(meta$SYMBOL)
  )
)

# ------------------------------------------------------------
# Significant DEGs
# ------------------------------------------------------------

sig <- meta |>
  
  dplyr::filter(
    
    rem_padj < PADJ_CUTOFF,
    
    abs(randomSummary) >= LFC_CUTOFF,
    
    !is.na(SYMBOL),
    
    SYMBOL != ""
  ) |>
  
  dplyr::arrange(
    randomP,
    dplyr::desc(
      abs(randomSummary)
    )
  ) |>
  
  dplyr::distinct(
    SYMBOL,
    .keep_all = TRUE
  ) |>
  
  dplyr::transmute(
    
    Symbol = SYMBOL,
    
    pooled_log2FC =
      randomSummary
  )

message(
  "\nSignificant DEGs: ",
  nrow(sig)
)

if (
  nrow(sig) < 5
) {
  
  stop(
    "Too few significant DEGs for PPI analysis."
  )
}

# ------------------------------------------------------------
# STRING limit
# ------------------------------------------------------------

if (
  nrow(sig) > STRING_MAX_NODES
) {
  
  message(
    "More than 2000 genes detected."
  )
  
  message(
    "Keeping top ",
    STRING_MAX_NODES,
    " significant genes."
  )
  
  sig <-
    head(
      sig,
      STRING_MAX_NODES
    )
}

# ------------------------------------------------------------
# Direction summary
# ------------------------------------------------------------

message(
  "Upregulated: ",
  sum(
    sig$pooled_log2FC > 0
  )
)

message(
  "Downregulated: ",
  sum(
    sig$pooled_log2FC < 0
  )
)

# ------------------------------------------------------------
# LFC map
# ------------------------------------------------------------

lfc_map <-
  setNames(
    sig$pooled_log2FC,
    sig$Symbol
  )

# ------------------------------------------------------------
# STRING
# ------------------------------------------------------------

message(
  "\nQuerying STRING database..."
)

network_df <-
  string_api_network(
    sig$Symbol
  )

if (
  is.null(network_df)
) {
  
  stop(
    "STRING API returned no network."
  )
}

message(
  "STRING interactions returned: ",
  nrow(network_df)
)

# ------------------------------------------------------------
# Build graph
# ------------------------------------------------------------

message(
  "\nBuilding PPI network..."
)

ppi_graph <-
  build_ppi_graph(
    network_df,
    lfc_map
  )

if (
  is.null(ppi_graph)
) {
  
  stop(
    "No STRING interactions found."
  )
}

message(
  "Network nodes: ",
  igraph::vcount(ppi_graph)
)

message(
  "Network edges: ",
  igraph::ecount(ppi_graph)
)

# ------------------------------------------------------------
# Metrics
# ------------------------------------------------------------

metrics <-
  compute_metrics(
    ppi_graph
  )

hub_genes <-
  get_hub_genes(
    metrics
  )

message(
  "Hub genes: ",
  nrow(hub_genes)
)

# ------------------------------------------------------------
# Communities
# ------------------------------------------------------------

ppi_graph <-
  detect_communities(
    ppi_graph
  )

modules_df <- tibble(
  
  Symbol =
    igraph::V(ppi_graph)$name,
  
  Module =
    igraph::V(ppi_graph)$community,
  
  Degree =
    igraph::degree(ppi_graph),
  
  Log2FC =
    igraph::V(ppi_graph)$pooled_log2FC
  
) |>
  
  dplyr::arrange(
    Module,
    dplyr::desc(Degree)
  )

# ============================================================
# 10. Export tables
# ============================================================

rio::export(
  metrics,
  file.path(
    OUT_CSV,
    paste0(
      PREFIX,
      "_network_metrics.csv"
    )
  )
)

rio::export(
  hub_genes,
  file.path(
    OUT_CSV,
    paste0(
      PREFIX,
      "_hub_genes.csv"
    )
  )
)

rio::export(
  modules_df,
  file.path(
    OUT_CSV,
    paste0(
      PREFIX,
      "_network_modules.csv"
    )
  )
)

rio::export(
  sig,
  file.path(
    OUT_CSV,
    paste0(
      PREFIX,
      "_network_input_DEGs.csv"
    )
  )
)

rio::export(
  network_df,
  file.path(
    OUT_CSV,
    paste0(
      PREFIX,
      "_STRING_interactions.csv"
    )
  )
)

# ============================================================
# 11. Save figures
# ============================================================

message(
  "\nSaving network figures..."
)

save_network_plot(
  
  ppi_graph,
  
  "Glioblastoma PPI Network",
  
  file.path(
    OUT_FIG,
    paste0(
      PREFIX,
      "_ppi_network.png"
    )
  )
)

save_hub_barplot(
  
  hub_genes,
  
  file.path(
    OUT_FIG,
    paste0(
      PREFIX,
      "_hub_genes_barplot.png"
    )
  )
)

# ============================================================
# 12. Final summary
# ============================================================

message(
  "\n============================================================"
)

message(
  "\nPPI NETWORK ANALYSIS COMPLETED SUCCESSFULLY"
)

message(
  "\n============================================================"
)

message(
  "\nMeta-analysis genes: ",
  nrow(meta)
)

message(
  "\nSignificant DEGs: ",
  nrow(sig)
)

message(
  "\nNetwork nodes: ",
  igraph::vcount(ppi_graph)
)

message(
  "\nNetwork edges: ",
  igraph::ecount(ppi_graph)
)

message(
  "\nHub genes: ",
  nrow(hub_genes)
)

message(
  "\nTables saved to:"
)

message(
  OUT_CSV
)

message(
  "\nFigures saved to:"
)

message(
  OUT_FIG
)

message(
  "\n============================================================\n"
)


parse("scripts/09_network_analysis.R")

source("scripts/09_network_analysis.R")

list.files("results/tables/network")
list.files("results/figures/network")

hub_genes <- rio::import(
  "results/tables/network/glioblastoma_hub_genes.csv"
)

head(hub_genes, 20)

hub_genes


network_metrics <- rio::import(
  "results/tables/network/glioblastoma_network_metrics.csv"
)

head(network_metrics, 20)
