
pak::pkg_install(c("tximport", "rtracklayer", "rio", "tidyverse"))

# ==============================
# Load packages
# ==============================

library(tximport)
library(rtracklayer)
library(tidyverse)
library(rio)

# ==============================
# Create output folder
# ==============================

dir.create("data/raw_counts",
           recursive = TRUE,
           showWarnings = FALSE)

# ==============================
# GEO datasets
# ==============================

geo_ids <- c(
  "GSE147352",
  "GSE151352"
)

# ==============================
# Build tx2gene from GENCODE v49
# ==============================
tx2gene <- gtf |>
  as.data.frame() |>
  dplyr::filter(type == "transcript") |>
  dplyr::select(transcript_id, gene_id) |>
  dplyr::mutate(
    transcript_id = sub("\\..*$", "", transcript_id)
  ) |>
  distinct()


# ==============================
# Process each dataset
# ==============================

for (geo in geo_ids) {
  
  message("Processing: ", geo)
  
  salmon_dir <- file.path("data", geo, "salmon")
  
  samples <- list.dirs(
    salmon_dir,
    recursive = FALSE,
    full.names = FALSE
  )
  
  files <- file.path(
    salmon_dir,
    samples,
    "quant.sf"
  )
  
  names(files) <- samples
  
  txi <- tximport(
    files,
    type = "salmon",
    tx2gene = tx2gene,
    ignoreTxVersion = TRUE
  )
  
  counts <- as.data.frame(txi$counts)
  counts <- tibble::rownames_to_column(counts, "GeneID")
  
  export(
    counts,
    file.path(
      "data/raw_counts",
      paste0(geo, "_raw_counts.csv")
    )
  )
  
  message(geo, " completed")
}

message("All datasets finished.")




counts <- as.data.frame(txi$counts)
counts <- tibble::rownames_to_column(counts, "GeneID")

rio::export(
  counts,
  file.path(
    "data/raw_counts",
    paste0(geo, "_raw_counts.csv")
  )
)
list.files("data/raw_counts")

#RAW COUNTS FOR GSE151352

geo <- "GSE151352"

salmon_dir <- file.path("data", geo, "salmon")

samples <- list.dirs(
  salmon_dir,
  recursive = FALSE,
  full.names = FALSE
)

files <- file.path(
  salmon_dir,
  samples,
  "quant.sf"
)

names(files) <- samples

files

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE
)

dim(txi$counts)
head(txi$counts)

counts <- as.data.frame(txi$counts)
counts <- tibble::rownames_to_column(counts, "GeneID")

rio::export(
  counts,
  "data/raw_counts/GSE151352_raw_counts.csv"
)

list.files("data/raw_counts")








head(tx2gene)


