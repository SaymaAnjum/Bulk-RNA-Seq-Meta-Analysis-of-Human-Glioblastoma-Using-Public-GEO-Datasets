# ============================================================
# 05_Meta-Analysis of RNA-seq data using MetaVolcanoR
# Glioblastoma: Tumor vs Normal
# GEO datasets: GSE147352 and GSE151352
# ============================================================
#Random-effects REML meta-analysis

list.files("scripts")

source("scripts/05_meta_analysis_metavolcano.R")

names(studies)
colnames(studies$GSE147352)
colnames(studies$GSE151352)


summary(studies$GSE147352$lfcSE)

summary(studies$GSE151352$lfcSE)


sum(is.na(studies$GSE147352$lfcSE))
sum(is.na(studies$GSE151352$lfcSE))


head(
  studies$GSE147352[, c(
    "Gene_ID",
    "log2FoldChange",
    "lfcSE",
    "padj"
  )]
)

head(
  studies$GSE151352[, c(
    "Gene_ID",
    "log2FoldChange",
    "lfcSE",
    "padj"
  )]
)


packageVersion("MetaVolcanoR")

args(rem_mv)


summary(test_studies$GSE147352$vi)
summary(test_studies$GSE151352$vi)

test_studies <- lapply(
  studies,
  function(df) {
    df |>
      mutate(
        Gene_ID = as.character(Gene_ID),
        vi = lfcSE^2
      )
  }
)


summary(test_studies$GSE147352$vi)

summary(test_studies$GSE151352$vi)


test_meta <- rem_mv(
  diffexp       = test_studies,
  pcriteria     = "padj",
  foldchangecol = "log2FoldChange",
  genenamecol   = "Gene_ID",
  geneidcol     = NULL,
  collaps       = TRUE,
  vcol          = "vi",
  cvar          = FALSE,
  metathr       = 0.01,
  jobname       = "MetaVolcano_REM",
  outputfolder  = "results/figures/meta-analysis/",
  draw          = "PDF",
  ncores        = 1
)



test_one <- test_studies$GSE147352

test_one <- test_one |>
  dplyr::select(Gene_ID, log2FoldChange, vi, padj)

colnames(test_one)


test_two <- test_studies$GSE151352

test_two <- test_two |>
  dplyr::select(Gene_ID, log2FoldChange, vi, padj)

colnames(test_two)




MetaVolcanoR::remodel(
  test_one[1, ],
  "log2FoldChange",
  "vi"
)



MetaVolcanoR::remodel(
  test_two[1, ],
  "log2FoldChange",
  "vi"
)




common_ids <- intersect(
  test_studies$GSE147352$Gene_ID,
  test_studies$GSE151352$Gene_ID
)

length(common_ids)



meta_studies <- lapply(
  test_studies,
  function(df) {
    df |>
      filter(Gene_ID %in% common_ids) |>
      distinct(Gene_ID, .keep_all = TRUE)
  }
)


sapply(meta_studies, nrow)

sapply(meta_studies, function(x) sum(is.na(x$vi)))

sapply(meta_studies, function(x) length(unique(x$Gene_ID)))


lapply(
  meta_studies,
  function(x) {
    x |>
      select(Gene_ID, log2FoldChange, vi, padj) |>
      head()
  }
)


test_meta <- MetaVolcanoR::rem_mv(
  diffexp       = meta_studies,
  pcriteria     = "padj",
  foldchangecol = "log2FoldChange",
  genenamecol   = "Gene_ID",
  geneidcol     = NULL,
  collaps       = FALSE,
  vcol          = "vi",
  cvar          = FALSE,
  metathr       = 0.01,
  jobname       = "MetaVolcano_REM",
  outputfolder  = "results/figures/meta-analysis/",
  draw          = "PDF",
  ncores        = 1
)

one_gene <- merge(
  meta_studies$GSE147352[, c("Gene_ID", "log2FoldChange", "vi")],
  meta_studies$GSE151352[, c("Gene_ID", "log2FoldChange", "vi")],
  by = "Gene_ID",
  suffixes = c("_GSE147352", "_GSE151352")
)

colnames(one_gene)


one_gene[1, ]

MetaVolcanoR::remodel(
  one_gene[1, ],
  "log2FoldChange",
  "vi"
)

meta_results_list <- lapply(
  seq_len(nrow(one_gene)),
  function(i) {
    res <- tryCatch(
      MetaVolcanoR::remodel(
        one_gene[i, , drop = FALSE],
        "log2FoldChange",
        "vi"
      ),
      error = function(e) NULL
    )
    
    if (is.null(res)) {
      return(NULL)
    }
    
    as.data.frame(res)
  }
)
meta_results <- dplyr::bind_rows(meta_results_list)
dim(meta_results)
head(meta_results)
meta_results$Gene_ID <- one_gene$Gene_ID
head(meta_results[, c(
  "Gene_ID",
  "randomSummary",
  "randomCi.lb",
  "randomCi.ub",
  "randomP",
  "het_QE",
  "het_QEp",
  "error"
)])

meta_key_results <- meta_results |>
  filter(
    randomP < 0.05,
    abs(randomSummary) >= 1,
    error == FALSE
  )

dim(meta_key_results)


head(meta_key_results)
dir.create(
  "results/tables/meta-analysis",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "results/figures/meta-analysis",
  showWarnings = FALSE,
  recursive = TRUE
)
export(
  meta_results,
  "results/tables/meta-analysis/random_effect_model.csv"
)

export(
  meta_key_results,
  "results/tables/meta-analysis/filtered_meta_degs.csv"
)
anno_files <- list.files(
  "results/tables/annotated",
  pattern = "[.]csv$",
  full.names = TRUE
)

gene_annotation <- anno_files |>
  lapply(function(f) {
    import(f) |>
      select(any_of(c(
        "Gene_ID",
        "Gene_Symbol",
        "Gene_Description"
      )))
  }) |>
  bind_rows() |>
  mutate(Gene_ID = as.character(Gene_ID)) |>
  distinct(Gene_ID, .keep_all = TRUE)



meta_annotated <- meta_results |>
  mutate(Gene_ID = as.character(Gene_ID)) |>
  left_join(
    gene_annotation,
    by = "Gene_ID"
  )

colnames(meta_annotated)


head(
  meta_annotated[, c(
    "Gene_ID",
    "Gene_Symbol",
    "randomSummary",
    "randomP"
  )]
)



meta_key_annotated <- meta_key_annotated |>
  mutate(
    Regulation = case_when(
      randomSummary > 0 ~ "UP",
      randomSummary < 0 ~ "DOWN",
      TRUE ~ "NS"
    )
  )

export(
  meta_key_annotated,
  "results/tables/meta-analysis/filtered_meta_degs_annotated.csv"
)

meta_key_annotated |>
  count(Regulation)

meta_key_summary <- meta_key_annotated |>
  count(Regulation) |>
  mutate(
    Percentage = n / sum(n) * 100
  )

print(meta_key_summary)

export(
  meta_key_summary,
  "results/tables/meta-analysis/meta_degs_summary_stats.csv"
)


volcano_data <- meta_annotated |>
  mutate(
    neg_log10_p = -log10(randomP),
    Significance = case_when(
      randomP < 0.05 & randomSummary >= 1 ~ "UP",
      randomP < 0.05 & randomSummary <= -1 ~ "DOWN",
      TRUE ~ "NS"
    )
  )


volcano_plot <- ggplot(
  volcano_data,
  aes(
    x = randomSummary,
    y = neg_log10_p
  )
) +
  geom_point(
    aes(color = Significance),
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
    subtitle = "Glioblastoma: Tumor vs Normal",
    x = "Pooled Random-Effects log2 Fold Change",
    y = "-log10(FDR)",
    color = "Significance"
  )
   +
  theme_bw()

print(volcano_plot)
ggsave(
  "results/figures/meta-analysis/Glioblastoma_MetaVolcano_REML.png",
  plot = volcano_plot,
  width = 9,
  height = 7,
  dpi = 600
)
file.exists(
  "results/figures/meta-analysis/Glioblastoma_MetaVolcano_REML.png"
)
#MetaVolcano plot
meta_annotated <- meta_annotated |>
  mutate(
    FDR = p.adjust(randomP, method = "BH")
  )
head(meta_annotated[, c(
  "Gene_ID",
  "Gene_Symbol",
  "randomSummary",
  "randomP",
  "FDR"
)])
volcano_data <- meta_annotated |>
  mutate(
    neg_log10_FDR = -log10(FDR),
    
    Significance = case_when(
      FDR < 0.05 & randomSummary >= 1 ~ "UP",
      FDR < 0.05 & randomSummary <= -1 ~ "DOWN",
      TRUE ~ "NS"
    )
  )

table(volcano_data$Significance)

meta_sig <- volcano_data |>
  filter(
    FDR < 0.05,
    abs(randomSummary) >= 1,
    error == FALSE
  )
nrow(meta_sig)


table(meta_sig$Significance)


library(ggplot2)

volcano_plot <- ggplot(
  volcano_data,
  aes(
    x = randomSummary,
    y = neg_log10_FDR
  )
) +
  geom_point(
    aes(color = Significance),
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
  volcano_plot <- volcano_plot +
  labs(
    title = "Meta-Analysis Volcano Plot",
    subtitle = "Glioblastoma: Tumor vs Normal — Random-Effects REML",
    x = "Pooled Random-Effects log2 Fold Change",
    y = "-log10(FDR)"
  )

print(volcano_plot)

volcano_plot


ggsave(
  "results/figures/meta-analysis/Glioblastoma_Meta_Analysis_Volcano.png",
  plot = volcano_plot,
  width = 9,
  height = 7,
  dpi = 600
)

file.exists(
  "results/figures/meta-analysis/Glioblastoma_Meta_Analysis_Volcano.png"
)
meta_sig <- meta_annotated |>
  dplyr::filter(
    FDR < 0.05,
    abs(randomSummary) >= 1,
    error == FALSE
  ) |>
  dplyr::mutate(
    Regulation = dplyr::case_when(
      randomSummary >= 1 ~ "UP",
      randomSummary <= -1 ~ "DOWN",
      TRUE ~ "NS"
    )
  ) |>
  dplyr::select(
    Gene_ID,
    Gene_Symbol,
    Gene_Description,
    randomSummary,
    randomCi.lb,
    randomCi.ub,
    randomP,
    FDR,
    het_QE,
    het_QEp,
    Regulation
  )

rio::export(
  meta_sig,
  "results/tables/meta-analysis/Glioblastoma_Final_Meta_DEGs_FDR.csv"
)

dim(meta_sig)
table(meta_sig$Regulation)

















































































































































