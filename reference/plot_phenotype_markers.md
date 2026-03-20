# Plot phenotype marker gene heatmaps (ComplexHeatmap)

Draws a **cell-type-agnostic** (global) or **cell-type-specific**
heatmap from the output of
[`find_phenotype_markers()`](https://brooksbenard.github.io/PhenoMapR/reference/find_phenotype_markers.md).
Expression is subset to the selected marker genes and ordered cells
**before** row-wise scaling, to keep transpose/scale costs small on
large matrices.

## Usage

``` r
plot_phenotype_markers(
  markers,
  expr_mat,
  meta,
  cell_id_col = "Cell",
  group_col,
  score_col,
  celltype_col = "celltype_original",
  celltype_palette = NULL,
  heatmap_type = c("global", "cell_type_specific"),
  top_n_markers = 20L,
  n_mark_labels = 5L,
  p_adj_threshold = 0.05,
  scale_clip = NULL,
  column_title = NULL,
  draw = TRUE,
  use_raster = FALSE
)
```

## Arguments

- markers:

  List returned by
  [`find_phenotype_markers()`](https://brooksbenard.github.io/PhenoMapR/reference/find_phenotype_markers.md)
  with elements `adverse_markers` and `favorable_markers`.

- expr_mat:

  Numeric matrix, genes \\\times\\ cells (column names = cell IDs).

- meta:

  Data frame of cell metadata with at least `cell_id_col`, `group_col`,
  `score_col`, and `celltype_col`.

- cell_id_col:

  Column in `meta` with cell IDs matching `colnames(expr_mat)`.

- group_col:

  Column with phenotype groups (`Most Favorable`, `Other`,
  `Most Adverse`).

- score_col:

  Column with continuous phenotype scores for the top color bar.

- celltype_col:

  Column with cell type labels.

- celltype_palette:

  Named vector of colors for cell types. If `NULL`,
  [`get_celltype_palette()`](https://brooksbenard.github.io/PhenoMapR/reference/get_celltype_palette.md)
  is used.

- heatmap_type:

  `"global"` (cell-type agnostic markers) or `"cell_type_specific"`
  (markers per cell type from `marker_scope = "cell_type_specific"`).

- top_n_markers:

  Maximum number of genes to keep per contrast block (per tail for
  global; per phenotype bin \\\times\\ cell type for
  cell-type-specific).

- n_mark_labels:

  Number of row labels to draw per block via
  [`ComplexHeatmap::anno_mark`](https://rdrr.io/pkg/ComplexHeatmap/man/anno_mark.html)
  (top genes by `avg_log2FC` within each block).

- p_adj_threshold:

  Only genes with `p_adj` below this value (and positive `avg_log2FC`)
  are candidates before ordering by log fold change.

- scale_clip:

  Length-2 numeric vector `c(lo, hi)` applied after row scaling (values
  outside are clipped). If `NULL`, uses `c(-3, 3)` for
  `heatmap_type = "global"` and `c(-5, 5)` for cell-type-specific.

- column_title:

  Optional title above the heatmap.

- draw:

  If `TRUE` (default), calls
  [`ComplexHeatmap::draw()`](https://rdrr.io/pkg/ComplexHeatmap/man/draw-dispatch.html).
  If `FALSE`, returns the `Heatmap` object invisibly.

- use_raster:

  Passed to `Heatmap()` (default `FALSE`).

## Value

Invisibly, the
[`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
object (or `NULL` if nothing was plotted).

## Details

Requires suggested packages **ComplexHeatmap** and **circlize**.

For `heatmap_type = "global"`, rows are favorable-marker genes then
adverse-only genes (no duplicate genes in the adverse block). For
`cell_type_specific`, rows follow phenotype bin then cell type, matching
column order.

## See also

[`find_phenotype_markers()`](https://brooksbenard.github.io/PhenoMapR/reference/find_phenotype_markers.md)
