#!/usr/bin/env Rscript

# WGCNA_118553.R
# Temporal cortex WGCNA prioritization for GSE118553.
#
# This script is designed to run after:
#   1. PCGS/run_pcgs_temporal.R
#   2. PCA_118553.R
#
# It loads the PCA-expanded temporal gene set, constructs a signed WGCNA
# network, computes intramodular connectivity to module eigengenes (kME),
# and selects a final gene panel of 200 to 300 genes for downstream analyses.
# The script also exports publication-ready visual summaries, including
# a gene dendrogram with module color annotation and an expression heatmap.

suppressPackageStartupMessages({
  library(WGCNA)
})

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

params <- list(
  input_pcgs_rds = file.path("outputs", "pcgs_temporal_pipeline.rds"),
  input_pca_csv = file.path("outputs", "PCA_118553_results.csv"),
  out_dir = "outputs",
  target_min = 200,
  target_max = 300,
  target_n = 250,
  kme_grid = seq(0.40, 0.95, by = 0.005),
  network_type = "signed",
  tom_type = "signed",
  min_module_size = 30,
  merge_cut_height = 0.25,
  max_block_size = 20000,
  scale_free_r2 = 0.80,
  heatmap_top_var_genes = 75
)

dir.create(params$out_dir, showWarnings = FALSE, recursive = TRUE)

log_step <- function(text) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), text))
}

require_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required input not found: ", path)
  }
}

diagnose_and_prepare_wgcna_inputs <- function(expr_matrix, sample_groups, min_genes_required) {
  issues <- character()

  if (!is.matrix(expr_matrix)) {
    expr_matrix <- as.matrix(expr_matrix)
    issues <- c(issues, "Expression object was coerced to a numeric matrix.")
  }

  mode(expr_matrix) <- "numeric"

  if (is.null(rownames(expr_matrix)) || any(rownames(expr_matrix) == "")) {
    stop("The expression matrix must contain non-empty gene symbols as row names.")
  }
  if (is.null(colnames(expr_matrix)) || any(colnames(expr_matrix) == "")) {
    stop("The expression matrix must contain non-empty sample identifiers as column names.")
  }

  duplicated_genes <- duplicated(rownames(expr_matrix))
  if (any(duplicated_genes)) {
    expr_matrix <- expr_matrix[!duplicated_genes, , drop = FALSE]
    issues <- c(
      issues,
      sprintf("Removed %d duplicated gene symbols before WGCNA.", sum(duplicated_genes))
    )
  }

  zero_variance <- apply(expr_matrix, 1, function(x) sd(x, na.rm = TRUE) == 0)
  zero_variance[is.na(zero_variance)] <- TRUE
  if (any(zero_variance)) {
    expr_matrix <- expr_matrix[!zero_variance, , drop = FALSE]
    issues <- c(
      issues,
      sprintf("Removed %d zero-variance genes before WGCNA.", sum(zero_variance))
    )
  }

  if (nrow(expr_matrix) < min_genes_required) {
    stop(
      "Too few genes remain after preprocessing for WGCNA. Remaining genes: ",
      nrow(expr_matrix)
    )
  }

  if (length(sample_groups) != ncol(expr_matrix)) {
    if (!is.null(names(sample_groups))) {
      matched_groups <- sample_groups[colnames(expr_matrix)]
      if (sum(!is.na(matched_groups)) == ncol(expr_matrix)) {
        sample_groups <- unname(matched_groups)
        issues <- c(issues, "Realigned sample group labels using sample identifiers.")
      } else {
        stop(
          "Sample group labels do not match the expression matrix columns after alignment. ",
          "Labels: ", length(sample_groups), ", samples: ", ncol(expr_matrix)
        )
      }
    } else {
      stop(
        "Sample group labels do not match the expression matrix columns. ",
        "Labels: ", length(sample_groups), ", samples: ", ncol(expr_matrix)
      )
    }
  }

  datExpr <- as.data.frame(t(expr_matrix))

  if (anyNA(datExpr) || any(!is.finite(as.matrix(datExpr)))) {
    datExpr <- as.data.frame(lapply(datExpr, function(v) {
      v[!is.finite(v)] <- NA_real_
      if (all(is.na(v))) {
        v[] <- 0
      } else {
        v[is.na(v)] <- median(v, na.rm = TRUE)
      }
      v
    }))
    rownames(datExpr) <- colnames(expr_matrix)
    issues <- c(issues, "Imputed NA/Inf values using gene-wise medians.")
  }

  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  if (!gsg$allOK) {
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
    sample_groups <- sample_groups[gsg$goodSamples]
    issues <- c(
      issues,
      sprintf(
        "Removed %d samples and %d genes flagged by goodSamplesGenes().",
        sum(!gsg$goodSamples),
        sum(!gsg$goodGenes)
      )
    )
  }

  if (ncol(datExpr) < min_genes_required) {
    stop(
      "Too few genes remain after WGCNA quality control. Remaining genes: ",
      ncol(datExpr)
    )
  }

  list(
    expr_matrix = expr_matrix,
    datExpr = datExpr,
    sample_groups = sample_groups,
    issues = unique(issues)
  )  
}

run_blockwise_modules_safe <- function(datExpr, params) {
  tryCatch(
    blockwiseModules(
      datExpr,
      power = params$soft_power,
      networkType = params$network_type,
      TOMType = params$tom_type,
      minModuleSize = params$min_module_size,
      mergeCutHeight = params$merge_cut_height,
      numericLabels = TRUE,
      pamRespectsDendro = FALSE,
      saveTOMs = FALSE,
      maxBlockSize = params$max_block_size,
      verbose = 3
    ),
    error = function(e) {
      stop(
        paste(
          "blockwiseModules() failed after input diagnostics.",
          "Likely causes include insufficient variance structure,",
          "an overly small effective sample size, or an unstable soft-threshold.",
          "Original error:",
          conditionMessage(e)
        )
      )
    }
  )
}

select_gene_panel <- function(kme_table, target_min, target_max, target_n, kme_grid) {
  scan_df <- data.frame(
    kME_threshold = kme_grid,
    n_genes = vapply(
      kme_grid,
      function(threshold) sum(kme_table$abs_kME >= threshold, na.rm = TRUE),
      integer(1)
    )
  )

  in_range <- which(scan_df$n_genes >= target_min & scan_df$n_genes <= target_max)

  if (length(in_range) > 0) {
    best_idx <- in_range[which.min(abs(scan_df$n_genes[in_range] - target_n))]
    threshold_selected <- scan_df$kME_threshold[best_idx]
    gene_panel <- kme_table[kme_table$abs_kME >= threshold_selected, , drop = FALSE]
    selection_mode <- "threshold_within_target_range"
  } else {
    under_cap <- which(scan_df$n_genes <= target_max)

    if (length(under_cap) > 0) {
      best_idx <- under_cap[which.min(abs(scan_df$n_genes[under_cap] - target_n))]
      threshold_selected <- scan_df$kME_threshold[best_idx]
      gene_panel <- kme_table[kme_table$abs_kME >= threshold_selected, , drop = FALSE]
      selection_mode <- "threshold_under_target_max_without_hard_cap"
    } else {
      threshold_distance <- abs(scan_df$n_genes - target_n)
      best_idx <- which.min(threshold_distance)
      threshold_selected <- scan_df$kME_threshold[best_idx]
      gene_panel <- kme_table[kme_table$abs_kME >= threshold_selected, , drop = FALSE]
      selection_mode <- "rank_adjusted_to_target_n"

      if (nrow(kme_table) < target_min) {
        stop(
          "Fewer than ", target_min, " genes are available after WGCNA scoring. ",
          "Current total: ", nrow(kme_table)
        )
      }

      gene_panel <- kme_table[seq_len(min(target_n, nrow(kme_table))), , drop = FALSE]
    }
  }

  gene_panel <- gene_panel[order(-gene_panel$abs_kME, gene_panel$gene_symbol), , drop = FALSE]
  rownames(gene_panel) <- NULL

  list(
    scan_df = scan_df,
    threshold_selected = threshold_selected,
    selection_mode = selection_mode,
    gene_panel = gene_panel
  )
}

plot_gene_dendrogram <- function(dat_expr_selected, module_colors_selected, out_file) {
  gene_tree <- hclust(dist(t(dat_expr_selected)), method = "average")

  grDevices::pdf(out_file, width = 14, height = 8)
  plotDendroAndColors(
    gene_tree,
    colors = module_colors_selected,
    groupLabels = "Module",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = "Temporal Cortex WGCNA Gene Dendrogram"
  )
  grDevices::dev.off()
}

plot_expression_heatmap <- function(dat_expr_selected, sample_groups, out_file, top_n = 75) {
  gene_variance <- apply(dat_expr_selected, 2, var, na.rm = TRUE)
  top_n <- min(top_n, ncol(dat_expr_selected))
  selected_columns <- names(sort(gene_variance, decreasing = TRUE))[seq_len(top_n)]
  heatmap_matrix <- t(scale(dat_expr_selected[, selected_columns, drop = FALSE]))
  heatmap_matrix[!is.finite(heatmap_matrix)] <- 0

  group_colors <- labels2colors(as.numeric(factor(sample_groups)))
  heatmap_limit <- max(abs(heatmap_matrix), na.rm = TRUE)
  if (!is.finite(heatmap_limit) || heatmap_limit == 0) {
    heatmap_limit <- 1
  }

  grDevices::pdf(out_file, width = 14, height = 10)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  }, add = FALSE)
  layout(matrix(c(1, 2), nrow = 2), heights = c(4.5, 1))
  par(mar = c(8, 8, 4, 2))
  labeledHeatmap(
    Matrix = heatmap_matrix,
    xLabels = rownames(dat_expr_selected),
    yLabels = rownames(heatmap_matrix),
    ySymbols = rownames(heatmap_matrix),
    colorLabels = FALSE,
    colors = blueWhiteRed(100),
    textMatrix = NULL,
    setStdMargins = FALSE,
    cex.lab.x = 0.6,
    cex.lab.y = 0.7,
    zlim = c(-heatmap_limit, heatmap_limit),
    main = "Temporal Cortex WGCNA Expression Heatmap"
  )
  par(mar = c(6, 8, 2, 2))
  plot.new()
  plot.window(xlim = c(0, nrow(dat_expr_selected)), ylim = c(0, 1))
  rect(
    xleft = seq(0, nrow(dat_expr_selected) - 1),
    ybottom = 0,
    xright = seq_len(nrow(dat_expr_selected)),
    ytop = 1,
    col = group_colors,
    border = NA
  )
  axis(1, at = seq_len(nrow(dat_expr_selected)) - 0.5, labels = rownames(dat_expr_selected), las = 2, cex.axis = 0.6)
  title(main = "Sample Group Annotation", line = 0.5)
}

log_step("Loading temporal PCGS and PCA outputs")
require_file(params$input_pcgs_rds)
require_file(params$input_pca_csv)

pcgs_temporal <- readRDS(params$input_pcgs_rds)
pca_results <- read.csv(params$input_pca_csv, stringsAsFactors = FALSE, check.names = FALSE)

if (!all(c("expr_gene_all", "group_temporal") %in% names(pcgs_temporal))) {
  stop("The PCGS session object must contain 'expr_gene_all' and 'group_temporal'.")
}
if (!("gene_symbol" %in% colnames(pca_results))) {
  stop("The PCA results table must contain a 'gene_symbol' column.")
}

expr_gene_all <- as.matrix(pcgs_temporal$expr_gene_all)
mode(expr_gene_all) <- "numeric"

expanded_genes <- unique(na.omit(as.character(pca_results$gene_symbol)))
expanded_genes <- intersect(expanded_genes, rownames(expr_gene_all))

if (length(expanded_genes) < params$target_min) {
  stop(
    "Too few PCA-expanded genes overlap with the temporal expression matrix. ",
    "Matched genes: ", length(expanded_genes)
  )
}

expr_temporal_expanded <- expr_gene_all[expanded_genes, , drop = FALSE]
sample_group_temporal <- as.character(pcgs_temporal$group_temporal)
if (length(sample_group_temporal) != ncol(expr_temporal_expanded)) {
  stop("Sample group labels do not match the temporal expression matrix columns.")
}

log_step(sprintf(
  "Matched %d PCA-expanded genes across %d temporal samples",
  nrow(expr_temporal_expanded), ncol(expr_temporal_expanded)
))

log_step("Running deep input diagnostics and automatic repair for WGCNA")
input_diagnostics <- diagnose_and_prepare_wgcna_inputs(
  expr_matrix = expr_temporal_expanded,
  sample_groups = sample_group_temporal,
  min_genes_required = params$target_min
)
expr_temporal_expanded <- input_diagnostics$expr_matrix
datExpr <- input_diagnostics$datExpr
sample_group_temporal <- input_diagnostics$sample_groups

if (length(input_diagnostics$issues) > 0) {
  for (issue in input_diagnostics$issues) {
    log_step(paste("Input repair:", issue))
  }
} else {
  log_step("No structural WGCNA input issues were detected")
}

log_step(sprintf("WGCNA input matrix: %d samples x %d genes", nrow(datExpr), ncol(datExpr)))

log_step("Selecting the soft-thresholding power")
power_candidates <- c(1:10, seq(12, 30, by = 2))
soft_threshold_fit <- pickSoftThreshold(
  datExpr,
  powerVector = power_candidates,
  networkType = params$network_type,
  verbose = 5
)

fit_table <- soft_threshold_fit$fitIndices
eligible_powers <- fit_table$Power[fit_table$SFT.R.sq >= params$scale_free_r2]
if (length(eligible_powers) > 0) {
  soft_power <- min(eligible_powers)
} else {
  soft_power <- fit_table$Power[which.max(fit_table$SFT.R.sq)]
  warning(
    sprintf(
      "No power achieved SFT.R^2 >= %.2f; using the best available power (%d).",
      params$scale_free_r2, soft_power
    )
  )
}

log_step(sprintf("Selected soft-thresholding power: %d", soft_power))

log_step("Constructing the signed co-expression network")
network_params <- params
network_params$soft_power <- soft_power
network <- run_blockwise_modules_safe(datExpr, network_params)

module_colors_all <- labels2colors(network$colors)
names(module_colors_all) <- colnames(datExpr)
module_eigengenes <- orderMEs(network$MEs)

log_step("Computing module membership (kME)")
kme_matrix <- as.data.frame(cor(datExpr, module_eigengenes, use = "p"))
colnames(kme_matrix) <- paste0("kME_", colnames(module_eigengenes))

gene_names <- colnames(datExpr)
gene_modules <- module_colors_all[gene_names]
self_kme <- numeric(length(gene_names))
names(self_kme) <- gene_names

for (i in seq_along(gene_names)) {
  module_name <- gene_modules[i]
  module_column <- paste0("kME_ME", module_name)

  if (module_column %in% colnames(kme_matrix)) {
    self_kme[i] <- kme_matrix[i, module_column]
  } else {
    self_kme[i] <- kme_matrix[i, which.max(abs(as.numeric(kme_matrix[i, ])))]
  }
}

wgcna_ranked_table <- data.frame(
  gene_symbol = gene_names,
  module_color = unname(gene_modules),
  kME = unname(self_kme),
  abs_kME = abs(unname(self_kme)),
  stringsAsFactors = FALSE
)
wgcna_ranked_table <- wgcna_ranked_table[order(-wgcna_ranked_table$abs_kME, wgcna_ranked_table$gene_symbol), , drop = FALSE]
rownames(wgcna_ranked_table) <- NULL

log_step("Selecting a final WGCNA gene panel of 200 to 300 genes")
selection <- select_gene_panel(
  kme_table = wgcna_ranked_table,
  target_min = params$target_min,
  target_max = params$target_max,
  target_n = params$target_n,
  kme_grid = params$kme_grid
)

WGCNA_118553_results <- selection$gene_panel
selected_gene_symbols <- WGCNA_118553_results$gene_symbol
selected_gene_modules <- WGCNA_118553_results$module_color
selected_datExpr <- datExpr[, selected_gene_symbols, drop = FALSE]

log_step(sprintf(
  "Selected %d genes using %s (kME threshold %.3f)",
  nrow(WGCNA_118553_results),
  selection$selection_mode,
  selection$threshold_selected
))

write.csv(
  selection$scan_df,
  file.path(params$out_dir, "temporal_wgcna_kme_threshold_scan.csv"),
  row.names = FALSE
)
write.csv(
  wgcna_ranked_table,
  file.path(params$out_dir, "temporal_wgcna_ranked_genes.csv"),
  row.names = FALSE
)
write.csv(
  WGCNA_118553_results,
  file.path(params$out_dir, "WGCNA_118553_results.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(gene_symbol = selected_gene_symbols, stringsAsFactors = FALSE),
  file.path(params$out_dir, "temporal_wgcna_selected_genes.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    issue = if (length(input_diagnostics$issues) > 0) input_diagnostics$issues else "No issues detected",
    stringsAsFactors = FALSE
  ),
  file.path(params$out_dir, "temporal_wgcna_input_diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    parameter = c(
      "input_pcgs_rds",
      "input_pca_csv",
      "out_dir",
      "target_min",
      "target_max",
      "target_n",
      "kme_grid_start",
      "kme_grid_end",
      "kme_grid_step",
      "network_type",
      "tom_type",
      "min_module_size",
      "merge_cut_height",
      "max_block_size",
      "scale_free_r2",
      "heatmap_top_var_genes",
      "soft_power",
      "kME_threshold_selected",
      "selection_mode",
      "n_selected_genes"
    ),
    value = c(
      params$input_pcgs_rds,
      params$input_pca_csv,
      params$out_dir,
      params$target_min,
      params$target_max,
      params$target_n,
      min(params$kme_grid),
      max(params$kme_grid),
      unique(diff(params$kme_grid))[1],
      params$network_type,
      params$tom_type,
      params$min_module_size,
      params$merge_cut_height,
      params$max_block_size,
      params$scale_free_r2,
      params$heatmap_top_var_genes,
      soft_power,
      selection$threshold_selected,
      selection$selection_mode,
      nrow(WGCNA_118553_results)
    ),
    stringsAsFactors = FALSE
  ),
  file.path(params$out_dir, "temporal_wgcna_parameters.csv"),
  row.names = FALSE
)

log_step("Saving publication-ready WGCNA figures")
plot_gene_dendrogram(
  dat_expr_selected = selected_datExpr,
  module_colors_selected = selected_gene_modules,
  out_file = file.path(params$out_dir, "temporal_wgcna_gene_dendrogram.pdf")
)
plot_expression_heatmap(
  dat_expr_selected = selected_datExpr,
  sample_groups = sample_group_temporal,
  out_file = file.path(params$out_dir, "temporal_wgcna_expression_heatmap.pdf"),
  top_n = params$heatmap_top_var_genes
)

saveRDS(
  list(
    params = params,
    soft_power = soft_power,
    module_colors_all = module_colors_all,
    module_eigengenes = module_eigengenes,
    kme_matrix = kme_matrix,
    wgcna_ranked_table = wgcna_ranked_table,
    selection = selection,
    WGCNA_118553_results = WGCNA_118553_results
  ),
  file = file.path(params$out_dir, "temporal_wgcna_diagnostics.rds")
)

log_step("Temporal WGCNA prioritization completed")
