# Scoring single-cell data with PhenoMapR

## Overview

This vignette demonstrates using PhenoMapR on a single-cell RNA-seq
dataset. Here, we chose to highlight a pancreatic adenocarcinoma (PAAD)
dataset from Peng et al. (*Cell Research* 2019; CRA001160) due to it’s
broad cell type representation, sample number, and inclusion of normal
controls. Pre-processed expression and metadata were obtained from
[TISCH2](https://tisch.compbio.cn/home/). We score cells with both
**TCGA PAAD** and **PRECOG Pancreatic** meta-z scores, compare results
between the two references, then define prognostic groups, identify
marker genes, and visualize proportions by sample and cell type using
the original cell type annotations from the metadata file. Finally, we
compare our results with those of [Jolasun et
al.](https://www.nature.com/articles/s41467-025-66162-4/figures/2) using
their method
**\[SIDISH\]**(<https://www.nature.com/articles/s41467-025-66162-4>) on
the same dataset.

## CRA001160: Load data from Google Drive

We download the expression matrix (10X H5) and cell metadata (TSV) from
Google Drive, then build a Seurat object for scoring and visualization.

### Download and load

``` r
suppressPackageStartupMessages(library(PhenoMapR))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
knitr::opts_chunk$set(fig.width = 12, out.width = "100%")

report_timing <- function(step_name, t0, obj = NULL) {
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  mem_mb <- if (!is.null(obj)) format(object.size(obj), units = "MB") else "-"
  message(sprintf("[%s] Runtime: %.2f s | Memory: %s", step_name, elapsed, mem_mb))
}

if (!requireNamespace("googledrive", quietly = TRUE)) {
  stop("The 'googledrive' package is required. Install with: install.packages('googledrive')")
}
vignette_dir <- if (dir.exists("vignettes")) "vignettes" else if (dir.exists("Vignettes")) "Vignettes" else "."
if (!dir.exists(vignette_dir)) dir.create(vignette_dir, recursive = TRUE, showWarnings = FALSE)

options(googledrive_quiet = TRUE)
googledrive::drive_deauth()
# Expression matrix (10X H5)
path_h5 <- file.path(vignette_dir, "PAAD_CRA001160_expression.h5")
googledrive::drive_download(googledrive::as_id("1PolTXggREz8XmhutCLTQJGCfKxFAzqMl"), path_h5, overwrite = TRUE)
# Cell metadata
path_meta <- file.path(vignette_dir, "PAAD_CRA001160_CellMetainfo_table.tsv")
googledrive::drive_download(googledrive::as_id("17mqxnKOZJn0jW2iD9RV0wZeWsilAIwdu"), path_meta, overwrite = TRUE)

has_data <- file.exists(path_h5) && file.exists(path_meta)
knitr::opts_chunk$set(eval = has_data)
if (!has_data) {
  message("Could not download CRA001160 expression or metadata from Google Drive.")
} else {
  t0 <- Sys.time()
  expr_mat <- Seurat::Read10X_h5(path_h5)
  meta <- read.delim(path_meta, stringsAsFactors = FALSE, check.names = FALSE)
  report_timing("Read H5 + metadata", t0)

  # Align metadata to matrix columns (cell barcodes)
  cell_ids <- colnames(expr_mat)
  id_col <- NULL
  for (cand in c("Barcode", "barcode", "cell_id", "Cell", "cell_barcode", names(meta)[1])) {
    if (cand %in% names(meta) && length(intersect(meta[[cand]], cell_ids)) > 0) {
      id_col <- cand
      break
    }
  }
  if (is.null(id_col)) id_col <- names(meta)[1]
  meta <- meta[meta[[id_col]] %in% cell_ids, , drop = FALSE]
  meta <- meta[match(cell_ids, meta[[id_col]]), , drop = FALSE]
  rownames(meta) <- cell_ids

  # Detect cell type column (character/factor with 2--500 unique values)
  celltype_col <- NULL
  for (col in names(meta)) {
    if (!(is.character(meta[[col]]) || is.factor(meta[[col]]))) next
    nlev <- length(unique(meta[[col]][!is.na(meta[[col]])]))
    if (nlev >= 2L && nlev <= 500L) {
      celltype_col <- col
      break
    }
  }
  if (is.null(celltype_col)) celltype_col <- names(meta)[2]
  # Original cell type annotation (TISCH2 metadata) when present
  celltype_original_col <- if ("Celltype (original)" %in% names(meta)) "Celltype (original)" else NULL

  # Subsample cells for vignette runtime (avoid very large Seurat objects in CI/pkgdown)
  max_cells <- 25000L
  if (ncol(expr_mat) > max_cells) {
    set.seed(1)
    keep_cells <- sample(colnames(expr_mat), max_cells)
    expr_mat <- expr_mat[, keep_cells, drop = FALSE]
    meta <- meta[keep_cells, , drop = FALSE]
    cell_ids <- colnames(expr_mat)
    message(sprintf("Subsampled cells from %d to %d for vignette.", length(rownames(meta)), length(keep_cells)))
  }

  # Build Seurat object; avoid double-normalization when data is already normalized (e.g. TISCH2)
  t0 <- Sys.time()
  already_norm <- PhenoMapR:::is_likely_normalized(expr_mat)
  seurat <- Seurat::CreateSeuratObject(counts = expr_mat, meta.data = meta, assay = "RNA")
  if (already_norm) {
    Seurat::SetAssayData(seurat, layer = "data", new.data = expr_mat, assay = "RNA")
    message("Data detected as already normalized (e.g. TISCH2); skipped NormalizeData to avoid double normalization.")
  } else {
    seurat <- Seurat::NormalizeData(seurat, verbose = FALSE)
  }
  report_timing(paste0("Create Seurat", if (already_norm) " (data as-is)" else " + normalize"), t0, seurat)
  message(sprintf("Cells: %d | Genes: %d | Cell type column: %s", ncol(seurat), nrow(seurat), celltype_col))
}
```

    ## [Read H5 + metadata] Runtime: 7.42 s | Memory: -

    ## Subsampled cells from 25000 to 25000 for vignette.

    ## Warning in .M2v(x): sparse->dense coercion: allocating vector of size 3.9 GiB

    ## Data detected as already normalized (e.g. TISCH2); skipped NormalizeData to avoid double normalization.

    ## [Create Seurat (data as-is)] Runtime: 3.93 s | Memory: 717.9 Mb

    ## Cells: 25000 | Genes: 21066 | Cell type column: Celltype (malignancy)

## Score cells with TCGA and PRECOG Pancreatic

We score each cell using both the **TCGA** (PAAD) and **PRECOG**
(Pancreatic) meta-z references, then add both score columns to the
Seurat object.

``` r
t0 <- Sys.time()

# Ensure the RNA data layer has gene names (required by PhenoMapR). If the
# data layer is empty or lacks rownames, copy from the counts layer.
expr_data <- tryCatch(
  Seurat::GetAssayData(seurat, layer = "data", assay = "RNA"),
  error = function(e) NULL
)
```

    ## Warning: Layer 'data' is empty

``` r
if (!is.null(expr_data) && (nrow(expr_data) == 0 || is.null(rownames(expr_data)))) {
  counts_data <- tryCatch(
    Seurat::GetAssayData(seurat, layer = "counts", assay = "RNA"),
    error = function(e) NULL
  )
  if (!is.null(counts_data) && !is.null(rownames(counts_data))) {
    expr_fix <- counts_data
    seurat <- Seurat::SetAssayData(seurat, layer = "data", new.data = expr_fix, assay = "RNA")
  }
}

scores_tcga <- PhenoMap(
  expression = seurat,
  reference = "tcga",
  cancer_type = "PAAD",
  assay = "RNA",
  slot = "data",
  verbose = TRUE
)
```

    ## Detected input type: seurat

    ## 4844 genes used for scoring against PAADCalculating scores...
    ## Completed scoring for PAAD

``` r
scores_precog <- PhenoMap(
  expression = seurat,
  reference = "precog",
  cancer_type = "Pancreatic",
  assay = "RNA",
  slot = "data",
  verbose = TRUE
)
```

    ## Detected input type: seurat

    ## 6556 genes used for scoring against PancreaticCalculating scores...
    ## Completed scoring for Pancreatic

``` r
seurat <- add_scores_to_seurat(seurat, scores_tcga)
```

    ## Added 1 score column(s) to Seurat metadata

``` r
seurat <- add_scores_to_seurat(seurat, scores_precog)
```

    ## Added 1 score column(s) to Seurat metadata

``` r
report_timing("Score TCGA + PRECOG", t0, seurat)
```

    ## [Score TCGA + PRECOG] Runtime: 4.65 s | Memory: 1423.6 Mb

## Scatterplot: TCGA vs PRECOG scores

Comparison of per-cell scores from the two references. Cells are colored
by the original cell type from the CellMetainfo file.

``` r
score_cols <- names(seurat@meta.data)[grep("weighted_sum_score", names(seurat@meta.data))]
score_tcga_col <- score_cols[grep("PAAD", score_cols, ignore.case = TRUE)][1]
score_precog_col <- score_cols[grep("Pancreatic", score_cols, ignore.case = TRUE)][1]
if (is.na(score_tcga_col)) score_tcga_col <- score_cols[1]
if (is.na(score_precog_col)) score_precog_col <- score_cols[2]

df_scatter <- seurat@meta.data
df_scatter$cell_type_anno <- df_scatter[[celltype_col]]
pal_ct <- PhenoMapR::get_celltype_palette(unique(as.character(df_scatter$cell_type_anno)))
ggplot(df_scatter, aes(x = .data[[score_precog_col]], y = .data[[score_tcga_col]], color = cell_type_anno)) +
  geom_point(alpha = 0.3, size = 0.5) +
  scale_color_manual(values = pal_ct, name = "Cell type") +
  theme_minimal() +
  labs(
    x = "PRECOG Pancreatic score",
    y = "TCGA PAAD score",
    title = "CRA001160: TCGA vs PRECOG pancreatic scores"
  ) +
  theme(legend.position = "right", legend.key.size = unit(0.4, "cm"))
```

![](single-cell_files/figure-html/scatter-tcga-precog-1.png)

## Score distribution: TCGA and PRECOG by cell type (ordered by PRECOG)

Boxplots of score distributions for both references, with cell types
ordered by median PRECOG score so that the most adverse (high PRECOG)
cell types appear toward one side.

``` r
# Order cell types by median PRECOG score (ascending: favorable left, adverse right)
med_precog <- setNames(
  tapply(seurat@meta.data[[score_precog_col]], seurat@meta.data[[celltype_col]], median, na.rm = TRUE),
  levels(factor(seurat@meta.data[[celltype_col]]))
)
med_precog <- med_precog[!is.na(med_precog)]
ct_order <- names(sort(med_precog))
seurat@meta.data[[celltype_col]] <- factor(seurat@meta.data[[celltype_col]], levels = ct_order)

# Long format: Reference (TCGA / PRECOG), Score, Cell type
dl <- rbind(
  data.frame(
    Reference = "TCGA PAAD",
    Score = seurat@meta.data[[score_tcga_col]],
    Cell_type = seurat@meta.data[[celltype_col]],
    stringsAsFactors = FALSE
  ),
  data.frame(
    Reference = "PRECOG Pancreatic",
    Score = seurat@meta.data[[score_precog_col]],
    Cell_type = seurat@meta.data[[celltype_col]],
    stringsAsFactors = FALSE
  )
)
dl$Reference <- factor(dl$Reference, levels = c("TCGA PAAD", "PRECOG Pancreatic"))

pal_ct <- PhenoMapR::get_celltype_palette(ct_order)
ggplot(dl, aes(x = Cell_type, y = Score, fill = Cell_type)) +
  geom_boxplot(outlier.alpha = 0.2) +
  facet_wrap(~ Reference, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = pal_ct, name = "Cell type") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    strip.text = element_text(size = 12)
  ) +
  labs(x = "Cell type (ordered by PRECOG median)", y = "PhenoMapR score", title = "CRA001160: Score distribution by reference and cell type")
```

![](single-cell_files/figure-html/score-distribution-boxplots-1.png)

## Score distribution by Celltype (original)

When the metadata includes **Celltype (original)** (TISCH2 original
annotations), we show the same score distributions grouped by that
column.

``` r
if (!is.null(celltype_original_col) && celltype_original_col %in% names(seurat@meta.data)) {
  med_precog_orig <- setNames(
    tapply(seurat@meta.data[[score_precog_col]], seurat@meta.data[[celltype_original_col]], median, na.rm = TRUE),
    levels(factor(seurat@meta.data[[celltype_original_col]]))
  )
  med_precog_orig <- med_precog_orig[!is.na(med_precog_orig)]
  ct_order_orig <- names(sort(med_precog_orig))
  seurat@meta.data$celltype_original <- factor(seurat@meta.data[[celltype_original_col]], levels = ct_order_orig)

  dl_orig <- rbind(
    data.frame(
      Reference = "TCGA PAAD",
      Score = seurat@meta.data[[score_tcga_col]],
      Cell_type = seurat@meta.data$celltype_original,
      stringsAsFactors = FALSE
    ),
    data.frame(
      Reference = "PRECOG Pancreatic",
      Score = seurat@meta.data[[score_precog_col]],
      Cell_type = seurat@meta.data$celltype_original,
      stringsAsFactors = FALSE
    )
  )
  dl_orig$Reference <- factor(dl_orig$Reference, levels = c("TCGA PAAD", "PRECOG Pancreatic"))
  pal_ct_orig <- PhenoMapR::get_celltype_palette(ct_order_orig)
  print(ggplot(dl_orig, aes(x = Cell_type, y = Score, fill = Cell_type)) +
    geom_boxplot(outlier.alpha = 0.2) +
    facet_wrap(~ Reference, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = pal_ct_orig, name = "Celltype (original)") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      strip.text = element_text(size = 12)
    ) +
    labs(x = "Celltype (original)", y = "PhenoMapR score", title = "CRA001160: Score by reference and Celltype (original)"))
} else {
  message("Celltype (original) not in metadata; skipping boxplots by original annotations.")
}
```

![](single-cell_files/figure-html/score-distribution-celltype-original-1.png)

## Prognostic groups and marker genes

We define Most Adverse and Most Favorable groups from the **PRECOG**
score (5th/95th percentiles), then run marker identification (adverse vs
rest, favorable vs rest) using the original cell type annotation for
context.

``` r
scores_df <- seurat@meta.data[, c(score_tcga_col, score_precog_col), drop = FALSE]
groups <- define_prognostic_groups(scores_df, percentile = 0.05, score_columns = score_precog_col)
seurat <- AddMetaData(seurat, groups)
group_col <- grep("prognostic_group", names(seurat@meta.data), value = TRUE)[1]

markers <- NULL
if (!is.na(group_col)) {
  markers <- find_prognostic_markers(
    seurat,
    group_labels = seurat@meta.data[[group_col]],
    group_column = NULL,
    assay = "RNA",
    slot = "data",
    max_cells_per_ident = 5000L
  )
  if (!is.null(markers)) {
    message("Adverse markers (top 5):")
    print(head(markers$adverse_markers, 5))
    message("Favorable markers (top 5):")
    print(head(markers$favorable_markers, 5))
  }
}
```

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 3.9 GiB

    ## Using Seurat FindMarkers: Most Adverse n=1250, Most Favorable n=1250

    ## Subsampled Other from 22500 to 5000 cells (memory limit)

    ## Adverse markers (top 5):

    ##   p_val avg_log2FC pct_in_group pct_rest     gene p_adj
    ## 1     0   2.916575        0.822    0.158    LAMB3     0
    ## 2     0   2.729147        0.770    0.139      SFN     0
    ## 3     0   2.450000        0.800    0.172   GPRC5A     0
    ## 4     0   2.302042        0.899    0.272   PHLDA2     0
    ## 5     0   1.652242        0.898    0.276 C19orf33     0

    ## Favorable markers (top 5):

    ##   p_val avg_log2FC pct_in_group pct_rest   gene p_adj
    ## 1     0   6.135377        0.491    0.053 CELA2B     0
    ## 2     0   6.359580        0.454    0.039   CTRL     0
    ## 3     0  -2.492736        0.520    0.917  ANXA2     0
    ## 4     0  -1.911960        0.540    0.932   RAC1     0
    ## 5     0   4.599349        0.411    0.046  PDIA2     0

## Marker heatmap

Top adverse and favorable marker genes; columns (cells) ordered by
PRECOG score. Cell type from CellMetainfo is shown in the annotation.

``` r
if (is.null(markers)) {
  message("Markers not available; skipping heatmap.")
} else {
  n_top <- 15
  expr <- tryCatch(
    Seurat::GetAssayData(seurat, layer = "data", assay = "RNA"),
    error = function(e) Seurat::GetAssayData(seurat, slot = "data", assay = "RNA")
  )
  adverse_pos <- markers$adverse_markers[markers$adverse_markers$avg_log2FC > 0, ]
  favorable_pos <- markers$favorable_markers[markers$favorable_markers$avg_log2FC > 0, ]
  top_genes <- c(
    head(adverse_pos$gene[order(adverse_pos$p_adj)], n_top),
    head(favorable_pos$gene[order(favorable_pos$p_adj)], n_top)
  )
  top_genes <- unique(top_genes[top_genes %in% rownames(expr)])
  if (length(top_genes) == 0) top_genes <- head(rownames(expr), 20)

  mat <- as.matrix(expr[top_genes, , drop = FALSE])
  mat_scaled <- t(scale(t(mat)))
  mat_scaled[mat_scaled < -3] <- -3
  mat_scaled[mat_scaled > 3] <- 3

  ord <- order(seurat@meta.data[[score_precog_col]])
  mat_scaled <- mat_scaled[, ord, drop = FALSE]
  meta_ord <- seurat@meta.data[ord, ]

  score_vals <- meta_ord[[score_precog_col]]
  # Use Celltype (original) for heatmap annotation when available
  ct_anno_col <- if (!is.null(celltype_original_col) && celltype_original_col %in% names(meta_ord)) celltype_original_col else celltype_col
  ann_col <- data.frame(
    `PRECOG Score` = score_vals,
    `Cell type` = factor(meta_ord[[ct_anno_col]]),
    `Prognostic group` = factor(meta_ord[[group_col]], levels = c("Most Adverse", "Other", "Most Favorable")),
    check.names = FALSE
  )
  rownames(ann_col) <- colnames(mat_scaled)

  pal_group <- c(`Most Adverse` = "#B2182B", Other = "#F7F7F7", `Most Favorable` = "#2166AC")
  pal_celltype <- PhenoMapR::get_celltype_palette(as.character(ann_col$`Cell type`))
  ann_colors <- list(
    `PRECOG Score` = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
    `Cell type` = pal_celltype,
    `Prognostic group` = pal_group
  )
  heatmap_colors <- if (requireNamespace("paletteer", quietly = TRUE)) {
    colorRampPalette(paletteer::paletteer_d("MexBrewer::Vendedora"))(100)
  } else {
    colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
  }
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    pheatmap::pheatmap(
      mat_scaled,
      scale = "none",
      cluster_cols = FALSE,
      cluster_rows = TRUE,
      show_colnames = FALSE,
      annotation_col = ann_col,
      annotation_colors = ann_colors,
      color = heatmap_colors,
      breaks = seq(-3, 3, length.out = 101),
      main = "Top adverse & favorable marker genes (CRA001160)",
      fontsize_row = 8
    )
  } else {
    heatmap(mat_scaled, scale = "none", Colv = NA, col = heatmap_colors, labCol = FALSE,
            main = "Top marker genes (CRA001160)")
  }
}
```

![](single-cell_files/figure-html/heatmap-markers-1.png)

## Proportion of prognostic cells by sample and cell type

Proportion of Most Adverse and Most Favorable cells (defined by PRECOG)
within each sample and cell type, using the original cell type
annotation from CellMetainfo.

``` r
sample_col <- NULL
for (cand in c("Sample", "sample", "Patient", "patient", "Sample_ID", "orig.ident")) {
  if (cand %in% names(seurat@meta.data)) {
    sample_col <- cand
    break
  }
}
if (is.null(sample_col)) sample_col <- names(seurat@meta.data)[1]

meta_plot <- seurat@meta.data
meta_plot$sample <- meta_plot[[sample_col]]
meta_plot$cell_type <- meta_plot[[celltype_col]]
meta_plot$prognostic_grp <- meta_plot[[group_col]]
meta_plot <- meta_plot[meta_plot$prognostic_grp %in% c("Most Adverse", "Most Favorable"), ]
meta_plot$sample <- factor(meta_plot$sample)
meta_plot$cell_type <- factor(meta_plot$cell_type)

if (nrow(meta_plot) > 0) {
  counts <- as.data.frame(table(meta_plot$sample, meta_plot$cell_type, meta_plot$prognostic_grp, dnn = c("sample", "cell_type", "pg")))
  totals <- aggregate(counts$Freq, by = list(sample = counts$sample, pg = counts$pg), FUN = sum)
  names(totals)[3] <- "total"
  counts <- merge(counts, totals, by = c("sample", "pg"))
  counts$proportion <- counts$Freq / counts$total
  counts$proportion[counts$total == 0] <- 0
  sample_lev <- levels(counts$sample)
  counts$x_num <- as.numeric(counts$sample) + ifelse(counts$pg == "Most Adverse", -0.18, 0.18)
  print(ggplot(counts, aes(x = x_num, y = cell_type, size = proportion, color = pg)) +
    geom_point(alpha = 0.85) +
    scale_color_manual(values = c(`Most Adverse` = "#B2182B", `Most Favorable` = "#2166AC"), name = "Prognostic group") +
    scale_size_continuous(range = c(0, 8), name = "Proportion") +
    scale_x_continuous(breaks = seq_along(sample_lev), labels = sample_lev) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_line(color = "grey90"),
      axis.title = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    ) +
    labs(x = "Sample", y = "Cell type", title = "Proportion of adverse vs favorable cells by cell type and sample (CRA001160)"))
} else {
  message("No Most Adverse or Most Favorable cells found; proportion plot skipped.")
}
```

![](single-cell_files/figure-html/proportion-by-sample-celltype-1.png)

## Conclusions

Here, we demonstrate that using PhenoMapR on single-cell RNAseq data in
a pancreatic cancer dataset successfully identified cells known to be
associated with disease (malignant and ductal type 1). These results
agree with those of Jolasun et al., however our PhenoMapR results
provide an increased level of granularity compared to other methods,
since we retain absolute and rank-orderd information regarding all
cell’s phenotype association. We also nominate those most associated
with favorable outcomes in PAAD, highlighting potential areas for
additional therapeutic focus.

## References

- **CRA001160**: Peng J, Sun B-F, Chen C-Y, et al. Single-cell RNA-seq
  highlights intra-tumoral heterogeneity and malignant progression in
  pancreatic ductal adenocarcinoma. *Cell Research*. 2019.
  <https://doi.org/10.1038/s41422-019-0195-y>. [GSA:
  CRA001160](https://ngdc.cncb.ac.cn/gsa/browse/CRA001160).

- **PRECOG 2.0**: Benard B, Lalgudi S, et al. PRECOG 2.0: an updated
  resource of pan-cancer gene-level prognostic meta-z scores. *Nucleic
  Acids Research*. 2026.

## Session Info

``` r
sessionInfo()
```

    ## R version 4.5.3 (2026-03-11)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] ggplot2_4.0.2      Seurat_5.4.0       SeuratObject_5.3.0 sp_2.2-1          
    ## [5] PhenoMapR_0.1.0   
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3     jsonlite_2.0.0         magrittr_2.0.4        
    ##   [4] spatstat.utils_3.2-2   farver_2.1.2           rmarkdown_2.30        
    ##   [7] fs_1.6.7               ragg_1.5.1             vctrs_0.7.1           
    ##  [10] ROCR_1.0-12            spatstat.explore_3.7-0 paletteer_1.7.0       
    ##  [13] htmltools_0.5.9        curl_7.0.0             sass_0.4.10           
    ##  [16] sctransform_0.4.3      parallelly_1.46.1      KernSmooth_2.23-26    
    ##  [19] bslib_0.10.0           htmlwidgets_1.6.4      desc_1.4.3            
    ##  [22] ica_1.0-3              plyr_1.8.9             plotly_4.12.0         
    ##  [25] zoo_1.8-15             cachem_1.1.0           igraph_2.2.2          
    ##  [28] mime_0.13              lifecycle_1.0.5        pkgconfig_2.0.3       
    ##  [31] Matrix_1.7-4           R6_2.6.1               fastmap_1.2.0         
    ##  [34] fitdistrplus_1.2-6     future_1.69.0          shiny_1.13.0          
    ##  [37] digest_0.6.39          rematch2_2.1.2         patchwork_1.3.2       
    ##  [40] tensor_1.5.1           prismatic_1.1.2        RSpectra_0.16-2       
    ##  [43] irlba_2.3.7            textshaping_1.0.5      labeling_0.4.3        
    ##  [46] progressr_0.18.0       spatstat.sparse_3.1-0  httr_1.4.8            
    ##  [49] polyclip_1.10-7        abind_1.4-8            compiler_4.5.3        
    ##  [52] gargle_1.6.1           bit64_4.6.0-1          withr_3.0.2           
    ##  [55] S7_0.2.1               fastDummies_1.7.5      MASS_7.3-65           
    ##  [58] tools_4.5.3            lmtest_0.9-40          otel_0.2.0            
    ##  [61] googledrive_2.1.2      httpuv_1.6.16          future.apply_1.20.2   
    ##  [64] goftest_1.2-3          glue_1.8.0             nlme_3.1-168          
    ##  [67] promises_1.5.0         grid_4.5.3             Rtsne_0.17            
    ##  [70] cluster_2.1.8.2        reshape2_1.4.5         generics_0.1.4        
    ##  [73] hdf5r_1.3.12           gtable_0.3.6           spatstat.data_3.1-9   
    ##  [76] tidyr_1.3.2            data.table_1.18.2.1    spatstat.geom_3.7-0   
    ##  [79] RcppAnnoy_0.0.23       ggrepel_0.9.7          RANN_2.6.2            
    ##  [82] pillar_1.11.1          stringr_1.6.0          limma_3.66.0          
    ##  [85] spam_2.11-3            RcppHNSW_0.6.0         later_1.4.8           
    ##  [88] splines_4.5.3          dplyr_1.2.0            lattice_0.22-9        
    ##  [91] survival_3.8-6         bit_4.6.0              deldir_2.0-4          
    ##  [94] tidyselect_1.2.1       miniUI_0.1.2           pbapply_1.7-4         
    ##  [97] knitr_1.51             gridExtra_2.3          scattermore_1.2       
    ## [100] xfun_0.56              statmod_1.5.1          matrixStats_1.5.0     
    ## [103] pheatmap_1.0.13        stringi_1.8.7          lazyeval_0.2.2        
    ## [106] yaml_2.3.12            evaluate_1.0.5         codetools_0.2-20      
    ## [109] tibble_3.3.1           cli_3.6.5              uwot_0.2.4            
    ## [112] xtable_1.8-8           reticulate_1.45.0      systemfonts_1.3.2     
    ## [115] jquerylib_0.1.4        Rcpp_1.1.1             globals_0.19.0        
    ## [118] spatstat.random_3.4-4  png_0.1-8              spatstat.univar_3.1-6 
    ## [121] parallel_4.5.3         pkgdown_2.2.0          presto_1.0.0          
    ## [124] dotCall64_1.2          listenv_0.10.1         viridisLite_0.4.3     
    ## [127] scales_1.4.0           ggridges_0.5.7         purrr_1.2.1           
    ## [130] rlang_1.1.7            cowplot_1.2.0
