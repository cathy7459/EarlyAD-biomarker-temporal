#!/usr/bin/env Rscript

# test set preparation.R
# External test-set preparation for Elastic Net evaluation.
#
# This script harmonizes GSE122063 with the temporal training pipeline by:
#   - retaining only Alzheimer's disease and control samples
#   - collapsing probes to gene symbols using the same strategy as the training set
#   - projecting the test cohort into the final training gene space
#   - applying the training-set scaling statistics to standardize the test matrix
#   - exporting a sample x gene matrix that can be used directly for Elastic Net inference

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(stringr)
})

# source(file.path("utils", "training_data_preparation_utils.R"))

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

resolve_training_reference_pair <- function(gene_list_candidates, scaling_stats_candidates) {
  gene_list_existing <- gene_list_candidates[file.exists(gene_list_candidates)]
  scaling_stats_existing <- scaling_stats_candidates[file.exists(scaling_stats_candidates)]

  if (length(gene_list_existing) == 0) {
    stop(
      "None of the expected training gene-list files were found: ",
      paste(gene_list_candidates, collapse = ", ")
    )
  }
  if (length(scaling_stats_existing) == 0) {
    stop(
      "None of the expected training scaling-stat files were found: ",
      paste(scaling_stats_candidates, collapse = ", ")
    )
  }

  candidate_pairs <- expand.grid(
    gene_list_path = gene_list_existing,
    scaling_stats_path = scaling_stats_existing,
    stringsAsFactors = FALSE
  )

  pair_summaries <- vector("list", nrow(candidate_pairs))

  for (i in seq_len(nrow(candidate_pairs))) {
    gene_table <- read.csv(candidate_pairs$gene_list_path[i], stringsAsFactors = FALSE, check.names = FALSE)
    scaling_table <- read.csv(candidate_pairs$scaling_stats_path[i], stringsAsFactors = FALSE, check.names = FALSE)

    if (!("gene_symbol" %in% colnames(gene_table)) || !("gene_symbol" %in% colnames(scaling_table))) {
      next
    }

    gene_symbols <- unique(as.character(gene_table$gene_symbol))
    gene_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]
    scaling_symbols <- unique(as.character(scaling_table$gene_symbol))
    scaling_symbols <- scaling_symbols[!is.na(scaling_symbols) & scaling_symbols != ""]

    if (length(gene_symbols) == 0 || length(scaling_symbols) == 0) {
      next
    }

    pair_summaries[[i]] <- data.frame(
      gene_list_path = candidate_pairs$gene_list_path[i],
      scaling_stats_path = candidate_pairs$scaling_stats_path[i],
      n_gene_list = length(gene_symbols),
      n_scaling_stats = length(scaling_symbols),
      same_gene_set = setequal(gene_symbols, scaling_symbols),
      stringsAsFactors = FALSE
    )
  }

  pair_summary_table <- do.call(rbind, pair_summaries)
  if (is.null(pair_summary_table) || nrow(pair_summary_table) == 0) {
    stop("No valid training reference pair could be constructed from the available files.")
  }

  valid_pairs <- pair_summary_table[
    pair_summary_table$same_gene_set &
      pair_summary_table$n_gene_list == pair_summary_table$n_scaling_stats,
    ,
    drop = FALSE
  ]

  if (nrow(valid_pairs) == 0) {
    write.csv(
      pair_summary_table,
      file.path("outputs", "temporal_test_training_reference_pair_diagnostics.csv"),
      row.names = FALSE
    )
    stop(
      "No consistent training gene-list/scaling-stat pair was found. ",
      "Diagnostics were written to outputs/temporal_test_training_reference_pair_diagnostics.csv"
    )
  }

  valid_pairs <- valid_pairs[order(-valid_pairs$n_gene_list, valid_pairs$gene_list_path), , drop = FALSE]
  chosen_pair <- valid_pairs[1, , drop = FALSE]

  list(
    gene_list_path = chosen_pair$gene_list_path[[1]],
    scaling_stats_path = chosen_pair$scaling_stats_path[[1]],
    diagnostics = pair_summary_table
  )
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

clean_symbol <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == "" | is.na(x)] <- NA_character_
  x <- str_split_fixed(x, "\\s*///\\s*|\\s*;\\s*|\\s*,\\s*", 2)[, 1]
  x <- trimws(x)
  x[x == "" | is.na(x)] <- NA_character_
  x
}

pick_symbol_col_from_gpl <- function(gpl_tbl) {
  cn <- colnames(gpl_tbl)
  c0 <- tolower(cn)

  preferred_exact <- c("GENE_SYMBOL", "Gene Symbol", "gene_symbol", "SYMBOL", "Symbol")
  exact_hit <- preferred_exact[preferred_exact %in% cn][1]
  if (!is.na(exact_hit) && !is.null(exact_hit)) {
    return(exact_hit)
  }

  idx1 <- which(
    c0 %in% c("gene symbol", "genesymbol", "symbol") |
      str_detect(c0, "gene\\s*symbol") |
      str_detect(c0, "\\bsymbol\\b")
  )
  if (length(idx1) > 0) return(cn[idx1[1]])

  idx2 <- which(str_detect(c0, "gene") & str_detect(c0, "symbol"))
  if (length(idx2) > 0) return(cn[idx2[1]])

  NA_character_
}

fetch_gpl_annotation <- function(gpl_id) {
  gpl <- tryCatch(
    getGEO(gpl_id, AnnotGPL = TRUE),
    error = function(e) NULL
  )

  if (is.null(gpl)) {
    log_step(sprintf(
      "AnnotGPL retrieval failed for %s; falling back to the submitter GPL table",
      gpl_id
    ))
    gpl <- getGEO(gpl_id)
  } else {
    log_step(sprintf(
      "Loaded GPL annotation for %s using AnnotGPL = TRUE",
      gpl_id
    ))
  }

  Table(gpl)
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

  # ------------------------------------------------------------
  # Resolve probe column explicitly
  # ------------------------------------------------------------
  probe_candidates <- c("ID", "ProbeID", "PROBE_ID", "probe_id", "ILMN_ID", "illumina_probe_id")
  probe_col <- probe_candidates[probe_candidates %in% colnames(gpl_tbl)][1]

  if (is.na(probe_col) || !nzchar(probe_col)) {
    stop(
      "Unable to locate a probe ID column in the GPL table.\nAvailable columns:\n",
      paste(colnames(gpl_tbl), collapse = ", ")
    )
  }

  # ------------------------------------------------------------
  # Keep only non-control probes when CONTROL_TYPE is available
  # ------------------------------------------------------------
  gpl_use <- gpl_tbl
  if ("CONTROL_TYPE" %in% colnames(gpl_use)) {
    gpl_use <- gpl_use[
      is.na(gpl_use$CONTROL_TYPE) |
        gpl_use$CONTROL_TYPE == "FALSE" |
        gpl_use$CONTROL_TYPE == FALSE,
      ,
      drop = FALSE
    ]
  }

  # ------------------------------------------------------------
  # Build mapping table
  # ------------------------------------------------------------
  map_df <- gpl_use[, c(probe_col, symbol_col), drop = FALSE]
  colnames(map_df) <- c("probe_raw", "symbol_raw")

  # Clean gene symbols
  map_df$symbol <- clean_symbol(map_df$symbol_raw)

  # Normalize probe IDs in GPL annotation
  map_df$probe_norm <- trimws(as.character(map_df$probe_raw))
  map_df$probe_norm <- gsub('^"|"$', "", map_df$probe_norm)
  map_df$probe_norm <- sub("\\.0+$", "", map_df$probe_norm)

  # Normalize probe IDs in expression matrix rownames
  expr_probe_rownames_raw <- rownames(expr_probe)
  expr_probe_rownames_norm <- trimws(as.character(expr_probe_rownames_raw))
  expr_probe_rownames_norm <- gsub('^"|"$', "", expr_probe_rownames_norm)
  expr_probe_rownames_norm <- sub("\\.0+$", "", expr_probe_rownames_norm)

  # Remove invalid rows
  map_df <- map_df[
    !is.na(map_df$probe_norm) & nzchar(map_df$probe_norm) &
      !is.na(map_df$symbol) & nzchar(map_df$symbol),
    c("probe_norm", "symbol"),
    drop = FALSE
  ]

  # Keep only probes present in the expression matrix
  map_df <- map_df[map_df$probe_norm %in% expr_probe_rownames_norm, , drop = FALSE]

  if (nrow(map_df) == 0) {
    stop("No overlap was found between the GPL annotation and the expression probes.")
  }

  # Re-map normalized probe IDs back to the exact rownames of expr_probe
  probe_match_idx <- match(map_df$probe_norm, expr_probe_rownames_norm)
  map_df$probe <- expr_probe_rownames_raw[probe_match_idx]

  # Keep unique probe-symbol pairs only
  map_df <- unique(map_df[, c("probe", "symbol"), drop = FALSE])

  if (nrow(map_df) == 0) {
    stop("No valid probe-to-gene mappings remained after filtering.")
  }

  log_step(sprintf(
    "Probe collapse used probe column '%s' and symbol column '%s' (overlap n = %d)",
    probe_col, symbol_col, nrow(map_df)
  ))

  # ------------------------------------------------------------
  # Subset expression matrix and collapse probes to genes
  # ------------------------------------------------------------
  expr_mapped <- expr_probe[map_df$probe, , drop = FALSE]

  # Ensure symbol order exactly matches the expression row order
  symbol_vec <- map_df$symbol[match(rownames(expr_mapped), map_df$probe)]

  expr_gene <- limma::avereps(expr_mapped, ID = symbol_vec)
  expr_gene <- as.matrix(expr_gene)
  mode(expr_gene) <- "numeric"

  list(
    expr_gene = expr_gene,
    probe_to_gene = map_df,
    symbol_col = symbol_col,
    probe_col = probe_col
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

params <- list(
  gse_id = "GSE122063",
  out_dir = "outputs",
  disease_column_candidates = c(
    "patient diagnosis:ch1",
    "diagnosis:ch1",
    "diagnosis",
    "disease state:ch1",
    "disease state",
    "disease:ch1",
    "disease",
    "group:ch1",
    "group",
    "status:ch1",
    "status"
  ),
  control_label_candidates = c(
    "control",
    "Control",
    "healthy control",
    "Healthy Control",
    "ND",
    "nd"
  ),
  case_label_candidates = c(
    "AD",
    "ad",
    "Alzheimer's disease",
    "alzheimers disease"
  ),
  training_gene_list_csv = file.path("outputs", "temporal_training_gene_list.csv"),
  training_scaling_stats_csv = file.path("outputs", "temporal_elasticnet_training_scaling_stats.csv")
)

normalize_test_group <- function(disease_values, control_label_candidates, case_label_candidates) {
  normalized_disease_values <- normalize_metadata_value(disease_values)
  normalized_control_candidates <- unique(normalize_metadata_value(control_label_candidates))
  normalized_case_candidates <- unique(normalize_metadata_value(case_label_candidates))

  ifelse(
    normalized_disease_values %in% normalized_control_candidates,
    "Control",
    ifelse(
      normalized_disease_values %in% normalized_case_candidates,
      "AD",
      NA_character_
    )
  )
}

log_step(sprintf("Downloading external test cohort %s", params$gse_id))
dir.create(params$out_dir, showWarnings = FALSE, recursive = TRUE)

training_gene_list_path <- resolve_input_path(c(
  params$training_gene_list_csv,
  file.path("Code_test", "outputs", "temporal_training_gene_list.csv")
))
training_scaling_stats_path <- resolve_input_path(c(
  params$training_scaling_stats_csv,
  file.path("Code_test", "outputs", "temporal_elasticnet_training_scaling_stats.csv")
))

training_reference_pair <- resolve_training_reference_pair(
  gene_list_candidates = c(
    params$training_gene_list_csv,
    file.path("Code_test", "outputs", "temporal_training_gene_list.csv")
  ),
  scaling_stats_candidates = c(
    params$training_scaling_stats_csv,
    file.path("Code_test", "outputs", "temporal_elasticnet_training_scaling_stats.csv")
  )
)
training_gene_list_path <- training_reference_pair$gene_list_path
training_scaling_stats_path <- training_reference_pair$scaling_stats_path

training_gene_table <- read.csv(training_gene_list_path, stringsAsFactors = FALSE, check.names = FALSE)
training_scaling_stats <- read.csv(training_scaling_stats_path, stringsAsFactors = FALSE, check.names = FALSE)

if (!("gene_symbol" %in% colnames(training_gene_table))) {
  stop("The training gene-list file must contain a 'gene_symbol' column.")
}

temporal_training_genes <- unique(as.character(training_gene_table$gene_symbol))
temporal_training_genes <- temporal_training_genes[!is.na(temporal_training_genes) & temporal_training_genes != ""]

gset <- getGEO(params$gse_id, GSEMatrix = TRUE, getGPL = FALSE)[[1]]
expr_probe_all <- exprs(gset)
pdat <- pData(gset)

disease_column <- resolve_metadata_column(
  pdat = pdat,
  column_candidates = params$disease_column_candidates,
  column_role = "disease"
)
disease_values <- normalize_string(pdat[[disease_column]])
group_test <- normalize_test_group(
  disease_values = disease_values,
  control_label_candidates = params$control_label_candidates,
  case_label_candidates = params$case_label_candidates
)
keep_test <- !is.na(group_test)

if (sum(keep_test) == 0) {
  stop(
    "No AD/control samples were identified in the test cohort. ",
    "Resolved disease column: ", disease_column,
    ". Unique disease values include: ", paste(utils::head(unique(disease_values), 10), collapse = ", ")
  )
}

expr_probe_test <- expr_probe_all[, keep_test, drop = FALSE]
group_test <- factor(group_test[keep_test], levels = c("Control", "AD"))

log_step(sprintf(
  "Retained %d external test samples after filtering to AD and control",
  ncol(expr_probe_test)
))
log_step(sprintf(
  "External test group counts: %s",
  paste(names(table(group_test)), table(group_test), collapse = ", ")
))

gpl_id <- annotation(gset)
if (is.na(gpl_id) || gpl_id == "") {
  stop("The GEO series does not provide a valid GPL identifier.")
}

log_step(sprintf("Collapsing probes to gene symbols using %s", gpl_id))
gpl_tbl <- fetch_gpl_annotation(gpl_id)
collapse_result <- collapse_probes_to_genes(expr_probe_test, gpl_tbl)
expr_gene_test <- collapse_result$expr_gene

mapped_test_genes <- intersect(temporal_training_genes, rownames(expr_gene_test))
expr_temporal_test_projected <- expr_gene_test[mapped_test_genes, , drop = FALSE]
expr_temporal_test_projected <- expr_temporal_test_projected[
  temporal_training_genes[temporal_training_genes %in% rownames(expr_temporal_test_projected)],
  ,
  drop = FALSE
]

log_step(sprintf(
  "Projected the external cohort into the final training gene space with %d overlapping genes",
  nrow(expr_temporal_test_projected)
))

X_test_raw <- t(expr_temporal_test_projected)
scale_result <- apply_saved_scaling(
  x = X_test_raw,
  scaling_stats = training_scaling_stats
)
temporal_elasticnet_test_expression <- scale_result$x_scaled

temporal_test_labels <- data.frame(
  sample_id = rownames(temporal_elasticnet_test_expression),
  group_temporal_test = as.character(group_test[match(
    rownames(temporal_elasticnet_test_expression),
    colnames(expr_probe_test)
  )]),
  y_binary = ifelse(
    as.character(group_test[match(
      rownames(temporal_elasticnet_test_expression),
      colnames(expr_probe_test)
    )]) == "AD",
    1,
    0
  ),
  stringsAsFactors = FALSE
)

temporal_test_parameters <- data.frame(
  parameter = c(
    "gse_id",
    "disease_column",
    "training_gene_list_csv",
    "training_scaling_stats_csv",
    "n_test_samples",
    "n_training_genes",
    "n_overlapping_genes",
    "n_missing_genes_imputed_to_training_mean"
  ),
  value = c(
    params$gse_id,
    disease_column,
    training_gene_list_path,
    training_scaling_stats_path,
    nrow(temporal_elasticnet_test_expression),
    length(temporal_training_genes),
    length(scale_result$overlapping_genes),
    length(scale_result$missing_genes)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  collapse_result$probe_to_gene,
  file.path(params$out_dir, "GSE122063_probe_to_gene_map.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(gene_symbol = temporal_training_genes, stringsAsFactors = FALSE),
  file.path(params$out_dir, "temporal_test_expected_gene_list.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(gene_symbol = scale_result$overlapping_genes, stringsAsFactors = FALSE),
  file.path(params$out_dir, "temporal_test_overlapping_gene_list.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(gene_symbol = scale_result$missing_genes, stringsAsFactors = FALSE),
  file.path(params$out_dir, "temporal_test_missing_gene_list.csv"),
  row.names = FALSE
)
write.csv(
  temporal_elasticnet_test_expression,
  file.path(params$out_dir, "temporal_elasticnet_test_expression_sample_x_gene.csv"),
  row.names = TRUE
)
write.csv(
  temporal_test_labels,
  file.path(params$out_dir, "temporal_elasticnet_test_labels.csv"),
  row.names = FALSE
)
write.csv(
  temporal_test_parameters,
  file.path(params$out_dir, "temporal_elasticnet_test_parameters.csv"),
  row.names = FALSE
)
write.csv(
  training_reference_pair$diagnostics,
  file.path(params$out_dir, "temporal_test_training_reference_pair_diagnostics.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    params = params,
    temporal_training_genes = temporal_training_genes,
    disease_column = disease_column,
    scale_result = scale_result,
    temporal_elasticnet_test_expression = temporal_elasticnet_test_expression,
    temporal_test_labels = temporal_test_labels
  ),
  file = file.path(params$out_dir, "temporal_test_preparation_diagnostics.rds")
)

log_step(sprintf(
  "External test preparation completed: %d samples x %d genes ready for Elastic Net inference",
  nrow(temporal_elasticnet_test_expression),
  ncol(temporal_elasticnet_test_expression)
))
