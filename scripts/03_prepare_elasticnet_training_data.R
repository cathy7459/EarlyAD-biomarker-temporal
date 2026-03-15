#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GEOquery)
})

#source(file.path("config", "training_data_preparation_config.R"))
#source(file.path("utils", "training_data_preparation_utils.R"))

run_training_data_preparation <- function(params = training_data_preparation_config) {
  dir.create(params$out_dir, showWarnings = FALSE, recursive = TRUE)

  log_step("Resolving the temporal WGCNA feature space")
  feature_space_path <- resolve_input_path(c(
    params$feature_space_csv,
    file.path("Code_test", "outputs", "WGCNA_118553_results.csv"),
    file.path("Code_test", "outputs", "temporal_wgcna_selected_genes.csv")
  ))

  feature_space_table <- read.csv(feature_space_path, stringsAsFactors = FALSE, check.names = FALSE)
  feature_column <- if ("gene_symbol" %in% colnames(feature_space_table)) {
    "gene_symbol"
  } else {
    grep("gene|symbol", colnames(feature_space_table), ignore.case = TRUE, value = TRUE)[1]
  }
  if (is.na(feature_column) || is.null(feature_column)) {
    stop("The feature-space file does not contain a gene-symbol column.")
  }

  temporal_wgcna_genes <- unique(trimws(as.character(feature_space_table[[feature_column]])))
  temporal_wgcna_genes <- temporal_wgcna_genes[!is.na(temporal_wgcna_genes) & temporal_wgcna_genes != ""]

  log_step(sprintf("Loaded %d temporal WGCNA genes", length(temporal_wgcna_genes)))

  log_step(sprintf("Downloading temporal training cohort %s", params$gse_id))
  gset <- getGEO(params$gse_id, GSEMatrix = TRUE, getGPL = FALSE)[[1]]
  expr_probe_all <- exprs(gset)
  pdat <- pData(gset)

  sample_filter <- filter_temporal_ad_control_samples(
    pdat = pdat,
    tissue_column_candidates = params$tissue_column_candidates,
    tissue_value_candidates = params$tissue_value_candidates,
    disease_column_candidates = params$disease_column_candidates,
    control_label_candidates = params$control_label_candidates,
    case_label_candidates = params$case_label_candidates
  )

  expr_probe_temporal_training <- expr_probe_all[, sample_filter$keep_samples, drop = FALSE]
  sample_ids_temporal_training <- colnames(expr_probe_temporal_training)
  group_temporal_training <- sample_filter$group_temporal_training

  log_step(sprintf(
    "Retained %d temporal-cortex samples after filtering disease state to AD and control",
    ncol(expr_probe_temporal_training)
  ))
  log_step(sprintf(
    "Temporal training group counts: %s",
    paste(names(table(group_temporal_training)), table(group_temporal_training), collapse = ", ")
  ))

  gpl_id <- annotation(gset)
  if (is.na(gpl_id) || gpl_id == "") {
    stop("The GEO series does not provide a valid GPL identifier.")
  }

  log_step(sprintf("Collapsing probes to gene symbols using %s", gpl_id))
  gpl_tbl <- Table(getGEO(gpl_id))
  collapse_result <- collapse_probes_to_genes(expr_probe_temporal_training, gpl_tbl)
  expr_gene_temporal_training <- collapse_result$expr_gene

  log_step(sprintf(
    "Collapsed the temporal training matrix to %d genes across %d samples",
    nrow(expr_gene_temporal_training), ncol(expr_gene_temporal_training)
  ))

  mapped_temporal_wgcna_genes <- intersect(temporal_wgcna_genes, rownames(expr_gene_temporal_training))
  if (length(mapped_temporal_wgcna_genes) == 0) {
    stop(
      "Too few temporal WGCNA genes were mapped into GSE132903 Temporal_Cortex. ",
      "Mapped genes: ", length(mapped_temporal_wgcna_genes)
    )
  }

  expr_temporal_projected <- expr_gene_temporal_training[mapped_temporal_wgcna_genes, , drop = FALSE]
  expr_temporal_projected <- expr_temporal_projected[
    temporal_wgcna_genes[temporal_wgcna_genes %in% rownames(expr_temporal_projected)],
    ,
    drop = FALSE
  ]

  log_step(sprintf(
    "Projected the training cohort into the temporal WGCNA feature space with %d genes",
    nrow(expr_temporal_projected)
  ))

  temporal_projected_expression_raw <- t(expr_temporal_projected)
  pruning_result <- prune_correlated_features(
    x = temporal_projected_expression_raw,
    ranked_genes = temporal_wgcna_genes,
    correlation_threshold = params$correlation_pruning_threshold
  )

  temporal_training_genes <- pruning_result$selected_genes
  expr_temporal_training_pruned <- expr_temporal_projected[temporal_training_genes, , drop = FALSE]

  univariate_result <- run_univariate_filter(
    expr_gene_by_sample = expr_temporal_training_pruned,
    group_factor = group_temporal_training,
    univariate_fdr_threshold = params$univariate_fdr_threshold,
    univariate_raw_p_threshold = params$univariate_raw_p_threshold,
    univariate_min_genes = params$univariate_min_genes
  )
  temporal_training_genes <- univariate_result$selected_genes
  expr_temporal_training_selected <- expr_temporal_training_pruned[temporal_training_genes, , drop = FALSE]
  temporal_training_feature_panel <- univariate_result$selected_table

  log_step(sprintf(
    paste(
      "Retained %d genes after correlation pruning and univariate filtering",
      "(correlation threshold %.2f, final max %d genes, mode = %s)"
    ),
    nrow(temporal_training_feature_panel),
    params$correlation_pruning_threshold,
    params$max_genes_final,
    univariate_result$selection_mode
  ))

  temporal_training_expression_raw <- t(expr_temporal_training_selected)
  scale_result <- impute_and_scale_matrix(temporal_training_expression_raw)
  temporal_elasticnet_training_expression <- scale_result$x_scaled

  temporal_training_labels <- data.frame(
    sample_id = rownames(temporal_elasticnet_training_expression),
    group_temporal_training = as.character(group_temporal_training[
      match(rownames(temporal_elasticnet_training_expression), colnames(expr_temporal_training_selected))
    ]),
    y_binary = ifelse(
      as.character(group_temporal_training[
        match(rownames(temporal_elasticnet_training_expression), colnames(expr_temporal_training_selected))
      ]) == "AD",
      1,
      0
    ),
    stringsAsFactors = FALSE
  )

  temporal_training_scaling_stats <- data.frame(
    gene_symbol = colnames(temporal_elasticnet_training_expression),
    median_training = unname(scale_result$medians[colnames(temporal_elasticnet_training_expression)]),
    mean_training = unname(scale_result$means[colnames(temporal_elasticnet_training_expression)]),
    sd_training = unname(scale_result$sds[colnames(temporal_elasticnet_training_expression)]),
    stringsAsFactors = FALSE
  )

  temporal_training_parameters <- data.frame(
    parameter = c(
      "gse_id",
      "tissue_column",
      "tissue_value_candidates",
      "disease_column",
      "control_label_candidates",
      "case_label_candidates",
      "feature_space_csv",
      "correlation_pruning_threshold",
      "max_genes_final_for_model_interpretation",
      "univariate_fdr_threshold",
      "univariate_raw_p_threshold",
      "univariate_min_genes",
      "selection_mode",
      "n_projected_genes",
      "n_selected_genes",
      "n_temporal_training_samples"
    ),
    value = c(
      params$gse_id,
      sample_filter$tissue_column,
      paste(params$tissue_value_candidates, collapse = " | "),
      sample_filter$disease_column,
      paste(params$control_label_candidates, collapse = " | "),
      paste(params$case_label_candidates, collapse = " | "),
      feature_space_path,
      params$correlation_pruning_threshold,
      params$max_genes_final,
      params$univariate_fdr_threshold,
      params$univariate_raw_p_threshold,
      params$univariate_min_genes,
      univariate_result$selection_mode,
      nrow(expr_temporal_projected),
      ncol(temporal_elasticnet_training_expression),
      nrow(temporal_elasticnet_training_expression)
    ),
    stringsAsFactors = FALSE
  )
  temporal_training_metadata_diagnostics <- build_metadata_diagnostics(
    sample_filter = sample_filter,
    sample_ids_kept = sample_ids_temporal_training
  )

  write.csv(
    collapse_result$probe_to_gene,
    file.path(params$out_dir, "GSE132903_temporal_probe_to_gene_map.csv"),
    row.names = FALSE
  )
  write.csv(
    temporal_training_feature_panel,
    file.path(params$out_dir, "temporal_training_feature_panel.csv"),
    row.names = FALSE
  )
  write.csv(
    univariate_result$stats_table,
    file.path(params$out_dir, "temporal_training_univariate_statistics.csv"),
    row.names = FALSE
  )
  write.csv(
    pruning_result$selection_table,
    file.path(params$out_dir, "temporal_training_correlation_pruning.csv"),
    row.names = FALSE
  )
  write.csv(
    data.frame(gene_symbol = temporal_training_genes, stringsAsFactors = FALSE),
    file.path(params$out_dir, "temporal_training_gene_list.csv"),
    row.names = FALSE
  )
  write.csv(
    temporal_elasticnet_training_expression,
    file.path(params$out_dir, "temporal_elasticnet_training_expression_sample_x_gene.csv"),
    row.names = TRUE
  )
  write.csv(
    temporal_training_labels,
    file.path(params$out_dir, "temporal_elasticnet_training_labels.csv"),
    row.names = FALSE
  )
  write.csv(
    temporal_training_scaling_stats,
    file.path(params$out_dir, "temporal_elasticnet_training_scaling_stats.csv"),
    row.names = FALSE
  )
  write.csv(
    temporal_training_parameters,
    file.path(params$out_dir, "temporal_elasticnet_training_parameters.csv"),
    row.names = FALSE
  )
  write.csv(
    temporal_training_metadata_diagnostics,
    file.path(params$out_dir, "temporal_elasticnet_training_metadata_diagnostics.csv"),
    row.names = FALSE
  )

  saveRDS(
    list(
      params = params,
      feature_space_path = feature_space_path,
      temporal_wgcna_genes = temporal_wgcna_genes,
      mapped_temporal_wgcna_genes = mapped_temporal_wgcna_genes,
      group_temporal_training = group_temporal_training,
      temporal_training_metadata_diagnostics = temporal_training_metadata_diagnostics,
      pruning_result = pruning_result,
      univariate_result = univariate_result,
      temporal_training_feature_panel = temporal_training_feature_panel,
      temporal_elasticnet_training_expression = temporal_elasticnet_training_expression,
      temporal_training_labels = temporal_training_labels,
      temporal_training_scaling_stats = temporal_training_scaling_stats
    ),
    file = file.path(params$out_dir, "temporal_training_preparation_diagnostics.rds")
  )

  log_step(sprintf(
    "Temporal training preparation completed: %d samples x %d genes ready for Elastic Net",
    nrow(temporal_elasticnet_training_expression),
    ncol(temporal_elasticnet_training_expression)
  ))

  invisible(
    list(
      temporal_elasticnet_training_expression = temporal_elasticnet_training_expression,
      temporal_training_labels = temporal_training_labels,
      temporal_training_scaling_stats = temporal_training_scaling_stats,
      temporal_training_feature_panel = temporal_training_feature_panel
    )
  )
}

run_training_data_preparation()
