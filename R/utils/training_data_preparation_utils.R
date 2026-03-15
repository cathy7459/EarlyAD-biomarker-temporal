#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)
})

source(file.path("R/utils", "PCGS.R"))

log_step <- function(text) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), text))
}

resolve_input_path <- function(path_candidates) {
  existing <- path_candidates[file.exists(path_candidates)]
  if (length(existing) == 0) {
    stop(
      "None of the expected input files were found: ",
      paste(path_candidates, collapse = ", ")
    )
  }
  existing[1]
}

normalize_string <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\s+", " ", x)
  x
}

normalize_column_name <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", "", x)
  x
}

normalize_metadata_value <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "", x)
  x
}

resolve_metadata_column <- function(pdat, column_candidates, column_role) {
  if (length(column_candidates) == 0) {
    stop("No candidate columns were provided for ", column_role, ".")
  }

  normalized_colnames <- normalize_column_name(colnames(pdat))
  normalized_candidates <- unique(normalize_column_name(column_candidates))

  matched_idx <- match(normalized_candidates, normalized_colnames)
  matched_idx <- matched_idx[!is.na(matched_idx)]

  if (length(matched_idx) == 0) {
    stop(
      "Missing required ", column_role, " column. Candidates checked: ",
      paste(column_candidates, collapse = ", "),
      ". Available columns: ",
      paste(colnames(pdat), collapse = ", ")
    )
  }

  colnames(pdat)[matched_idx[1]]
}

filter_temporal_ad_control_samples <- function(pdat, tissue_column_candidates, tissue_value_candidates, disease_column_candidates,
                                               control_label_candidates, case_label_candidates) {
  tissue_column <- resolve_metadata_column(pdat, tissue_column_candidates, "tissue")
  disease_column <- resolve_metadata_column(pdat, disease_column_candidates, "disease")

  tissue_values <- normalize_string(pdat[[tissue_column]])
  disease_values <- normalize_string(pdat[[disease_column]])

  normalized_tissue_values <- normalize_metadata_value(tissue_values)
  normalized_disease_values <- normalize_metadata_value(disease_values)
  normalized_tissue_candidates <- unique(normalize_metadata_value(tissue_value_candidates))
  normalized_control_candidates <- unique(normalize_metadata_value(control_label_candidates))
  normalized_case_candidates <- unique(normalize_metadata_value(case_label_candidates))

  keep_tissue <- normalized_tissue_values %in% normalized_tissue_candidates
  keep_group <- normalized_disease_values %in% c(normalized_control_candidates, normalized_case_candidates)
  keep_samples <- keep_tissue & keep_group

  group_temporal_training <- ifelse(
    normalized_disease_values[keep_samples] %in% normalized_control_candidates,
    "Control",
    ifelse(
      normalized_disease_values[keep_samples] %in% normalized_case_candidates,
      "AD",
      NA_character_
    )
  )

  if (sum(keep_samples) == 0) {
    stop(
      "No samples matched the requested temporal and disease filters. ",
      "Resolved tissue column: ", tissue_column,
      "; resolved disease column: ", disease_column,
      ". Unique tissue values include: ", paste(utils::head(unique(tissue_values), 10), collapse = ", "),
      ". Unique disease values include: ", paste(utils::head(unique(disease_values), 10), collapse = ", ")
    )
  }
  if (anyNA(group_temporal_training)) {
    stop("Some filtered samples could not be mapped to Control or AD.")
  }

  list(
    keep_samples = keep_samples,
    tissue_column = tissue_column,
    disease_column = disease_column,
    tissue_values = tissue_values,
    disease_values = disease_values,
    tissue_value_counts = sort(table(tissue_values), decreasing = TRUE),
    disease_value_counts = sort(table(disease_values), decreasing = TRUE),
    group_temporal_training = factor(group_temporal_training, levels = c("Control", "AD"))
  )
}

prune_correlated_features <- function(x, ranked_genes, correlation_threshold = 0.90, max_genes = NULL) {
  x <- as.matrix(x)
  mode(x) <- "numeric"

  ranked_genes <- ranked_genes[ranked_genes %in% colnames(x)]
  if (length(ranked_genes) == 0) {
    stop("No ranked genes overlapped with the projected training matrix.")
  }

  selected_genes <- character(0)
  exclusion_reason <- setNames(rep("not_reached", length(ranked_genes)), ranked_genes)

  for (gene_symbol in ranked_genes) {
    if (!is.null(max_genes) && length(selected_genes) >= max_genes) {
      exclusion_reason[gene_symbol] <- "max_genes_reached"
      next
    }

    if (length(selected_genes) == 0) {
      selected_genes <- c(selected_genes, gene_symbol)
      exclusion_reason[gene_symbol] <- "selected"
      next
    }

    correlation_with_selected <- suppressWarnings(cor(
      x[, gene_symbol],
      x[, selected_genes, drop = FALSE],
      use = "pairwise.complete.obs"
    ))
    correlation_with_selected <- abs(as.numeric(correlation_with_selected))
    correlation_with_selected <- correlation_with_selected[is.finite(correlation_with_selected)]

    if (length(correlation_with_selected) == 0 || max(correlation_with_selected) < correlation_threshold) {
      selected_genes <- c(selected_genes, gene_symbol)
      exclusion_reason[gene_symbol] <- "selected"
    } else {
      exclusion_reason[gene_symbol] <- "pruned_by_correlation"
    }
  }

  selection_table <- data.frame(
    gene_symbol = ranked_genes,
    selection_status = unname(exclusion_reason[ranked_genes]),
    rank_in_feature_space = seq_along(ranked_genes),
    stringsAsFactors = FALSE
  )

  list(
    selected_genes = selected_genes,
    selection_table = selection_table
  )
}

build_metadata_diagnostics <- function(sample_filter, sample_ids_kept) {
  data.frame(
    resolved_tissue_column = sample_filter$tissue_column,
    resolved_disease_column = sample_filter$disease_column,
    n_samples_kept = length(sample_ids_kept),
    kept_samples_preview = paste(utils::head(sample_ids_kept, 10), collapse = " | "),
    tissue_values_preview = paste(utils::head(names(sample_filter$tissue_value_counts), 10), collapse = " | "),
    disease_values_preview = paste(utils::head(names(sample_filter$disease_value_counts), 10), collapse = " | "),
    stringsAsFactors = FALSE
  )
}

run_univariate_filter <- function(expr_gene_by_sample, group_factor,
                                  univariate_fdr_threshold = 0.10,
                                  univariate_raw_p_threshold = 0.05,
                                  univariate_min_genes = 5) {
  expr_gene_by_sample <- as.matrix(expr_gene_by_sample)
  mode(expr_gene_by_sample) <- "numeric"

  design <- model.matrix(~ 0 + group_factor)
  colnames(design) <- make.names(colnames(design))
  fit <- limma::lmFit(expr_gene_by_sample, design)
  contrast_matrix <- limma::makeContrasts(
    AD_vs_Control = group_factorAD - group_factorControl,
    levels = design
  )
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, contrast_matrix))

  univariate_stats <- limma::topTable(fit2, coef = "AD_vs_Control", number = Inf, sort.by = "P")
  univariate_stats$gene_symbol <- rownames(univariate_stats)
  univariate_stats <- univariate_stats[, c("gene_symbol", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

  keep_fdr <- univariate_stats$adj.P.Val < univariate_fdr_threshold
  keep_raw_p <- univariate_stats$P.Value < univariate_raw_p_threshold

  if (sum(keep_fdr, na.rm = TRUE) >= univariate_min_genes) {
    selected_table <- univariate_stats[keep_fdr, , drop = FALSE]
    selection_mode <- "fdr_filtered_training_only"
  } else if (sum(keep_raw_p, na.rm = TRUE) >= univariate_min_genes) {
    selected_table <- univariate_stats[keep_raw_p, , drop = FALSE]
    selection_mode <- "raw_p_filtered_training_only"
  } else {
    selected_table <- univariate_stats
    selection_mode <- "ranked_by_univariate_p_training_only"
  }

  selected_table <- selected_table[order(selected_table$adj.P.Val, selected_table$P.Value, -abs(selected_table$logFC)), , drop = FALSE]
  rownames(selected_table) <- NULL

  list(
    stats_table = univariate_stats,
    selected_table = selected_table,
    selected_genes = selected_table$gene_symbol,
    selection_mode = selection_mode
  )
}

collapse_probes_to_genes <- function(expr_probe, gpl_tbl) {
  symbol_col <- pick_symbol_col_from_gpl(gpl_tbl)
  if (is.na(symbol_col)) {
    symbol_candidates <- grep("gene|symbol|ilmn_gene", colnames(gpl_tbl), ignore.case = TRUE, value = TRUE)
    if (length(symbol_candidates) == 0) {
      stop("Unable to locate a gene symbol column in the GPL table.")
    }
    symbol_col <- symbol_candidates[1]
  }

  probe_col <- colnames(gpl_tbl)[1]
  map_df <- gpl_tbl[, c(probe_col, symbol_col)]
  colnames(map_df) <- c("probe", "symbol_raw")
  map_df$symbol <- clean_symbol(map_df$symbol_raw)
  map_df <- map_df[!is.na(map_df$symbol), c("probe", "symbol"), drop = FALSE]
  map_df <- map_df[map_df$probe %in% rownames(expr_probe), , drop = FALSE]

  if (nrow(map_df) == 0) {
    stop("No overlap was found between the GPL annotation and the expression probes.")
  }

  expr_mapped <- expr_probe[map_df$probe, , drop = FALSE]
  expr_gene <- limma::avereps(expr_mapped, ID = map_df$symbol)
  expr_gene <- as.matrix(expr_gene)
  mode(expr_gene) <- "numeric"

  list(
    expr_gene = expr_gene,
    probe_to_gene = map_df,
    symbol_col = symbol_col
  )
}

impute_and_scale_matrix <- function(x) {
  x <- as.matrix(x)
  mode(x) <- "numeric"
  gene_medians <- rep(NA_real_, ncol(x))
  names(gene_medians) <- colnames(x)

  for (j in seq_len(ncol(x))) {
    x[!is.finite(x[, j]), j] <- NA_real_
    if (all(is.na(x[, j]))) {
      x[, j] <- 0
      gene_medians[j] <- 0
    } else {
      gene_medians[j] <- median(x[, j], na.rm = TRUE)
      x[is.na(x[, j]), j] <- median(x[, j], na.rm = TRUE)
    }
  }

  gene_means <- colMeans(x, na.rm = TRUE)
  gene_sds <- apply(x, 2, sd, na.rm = TRUE)
  gene_sds[!is.finite(gene_sds) | gene_sds == 0] <- 1

  x_scaled <- sweep(x, 2, gene_means, FUN = "-")
  x_scaled <- sweep(x_scaled, 2, gene_sds, FUN = "/")
  x_scaled[!is.finite(x_scaled)] <- 0

  list(
    x_scaled = x_scaled,
    medians = gene_medians,
    means = gene_means,
    sds = gene_sds
  )
}

apply_saved_scaling <- function(x, scaling_stats) {
  x <- as.matrix(x)
  mode(x) <- "numeric"

  required_columns <- c("gene_symbol", "mean_training", "sd_training")
  if (!all(required_columns %in% colnames(scaling_stats))) {
    stop(
      "Scaling statistics must contain columns: ",
      paste(required_columns, collapse = ", ")
    )
  }

  target_genes <- as.character(scaling_stats$gene_symbol)
  x_aligned <- matrix(
    NA_real_,
    nrow = nrow(x),
    ncol = length(target_genes),
    dimnames = list(rownames(x), target_genes)
  )

  overlapping_genes <- intersect(colnames(x), target_genes)
  x_aligned[, overlapping_genes] <- x[, overlapping_genes, drop = FALSE]

  gene_means <- setNames(as.numeric(scaling_stats$mean_training), target_genes)
  gene_sds <- setNames(as.numeric(scaling_stats$sd_training), target_genes)
  gene_sds[!is.finite(gene_sds) | gene_sds == 0] <- 1
  if ("median_training" %in% colnames(scaling_stats)) {
    gene_medians <- setNames(as.numeric(scaling_stats$median_training), target_genes)
  } else {
    gene_medians <- gene_means
  }

  for (gene_symbol in target_genes) {
    column_values <- x_aligned[, gene_symbol]
    column_values[!is.finite(column_values)] <- NA_real_
    column_values[is.na(column_values)] <- gene_medians[gene_symbol]
    x_aligned[, gene_symbol] <- column_values
  }

  x_scaled <- sweep(x_aligned, 2, gene_means[target_genes], FUN = "-")
  x_scaled <- sweep(x_scaled, 2, gene_sds[target_genes], FUN = "/")
  x_scaled[!is.finite(x_scaled)] <- 0

  list(
    x_scaled = x_scaled,
    overlapping_genes = overlapping_genes,
    missing_genes = setdiff(target_genes, overlapping_genes)
  )
}
