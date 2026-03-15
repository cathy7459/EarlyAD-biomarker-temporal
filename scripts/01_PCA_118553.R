#!/usr/bin/env Rscript

# PCA_118553.R
# Temporal cortex PCA-based feature expansion for GSE118553.
#
# This script is designed to run after `PCGS/run_pcgs_temporal.R`.
# It loads the temporal PCGS session object, performs PCA on the
# temporal core-gene expression space, quantifies gene-PC associations
# across the full temporal transcriptome, and expands the feature set
# using statistically supported simple and partial correlations.
#
# Selection strategy:
#   1. Search statistically meaningful parameter combinations.
#   2. Prefer combinations that already yield 1,000 to 3,000 genes.
#   3. Apply a hard cap at 3,000 genes only when the selected
#      combination exceeds that limit.
#
# Main outputs:
#   - outputs/PCA_118553_results.csv
#   - outputs/temporal_pca_expanded_genes.txt
#   - outputs/temporal_pca_expanded_expression_gene_x_sample.csv
#   - outputs/temporal_pca_expansion_diagnostics.rds

options(stringsAsFactors = FALSE)

params <- list(
  input_rds = file.path("outputs", "pcgs_temporal_pipeline.rds"),
  out_dir = "outputs",
  pc_cumvar = 0.70,
  min_pcs = 4,
  max_pcs = 12,
  drop_pc1 = FALSE,
  target_min = 1000,
  target_max = 3000,
  m_and_grid = c(4, 3, 2),
  alpha_grid = c(0.005, 0.01, 0.02),
  alpha_soft_grid = c(0.01, 0.02, 0.05),
  min_abs_simple_cor_grid = c(0.30, 0.25, 0.20),
  min_abs_partial_cor_grid = c(0.20, 0.15, 0.10)
)

dir.create(params$out_dir, showWarnings = FALSE, recursive = TRUE)

log_step <- function(text) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), text))
}

require_named_elements <- function(x, required_names, object_name) {
  missing_names <- setdiff(required_names, names(x))
  if (length(missing_names) > 0) {
    stop(
      object_name, " is missing required entries: ",
      paste(missing_names, collapse = ", ")
    )
  }
}

safe_scale <- function(x) {
  x <- as.matrix(x)
  mode(x) <- "numeric"

  if (nrow(x) == 0 || ncol(x) == 0) {
    stop("safe_scale() received an empty matrix.")
  }

  col_means <- colMeans(x, na.rm = TRUE)
  col_sds <- apply(x, 2, sd, na.rm = TRUE)

  col_means[!is.finite(col_means)] <- 0
  col_sds[!is.finite(col_sds) | col_sds == 0] <- 1

  x_scaled <- sweep(x, 2, col_means, FUN = "-")
  x_scaled <- sweep(x_scaled, 2, col_sds, FUN = "/")

  x_scaled[!is.finite(x_scaled)] <- 0
  x_scaled <- as.matrix(x_scaled)

  x_scaled
}


cor_to_p <- function(r, n) {
  r <- pmin(pmax(r, -0.999999), 0.999999)
  t_stat <- r * sqrt((n - 2) / (1 - r^2))
  2 * pt(-abs(t_stat), df = n - 2)
}

cor_with_pc <- function(x, y) {
  suppressWarnings(cor(x, y, use = "pairwise.complete.obs"))
}

partial_cor_onegene <- function(g, score_matrix) {
  n_pc <- ncol(score_matrix)
  out <- rep(NA_real_, n_pc)

  for (k in seq_len(n_pc)) {
    if (n_pc == 1) {
      out[k] <- cor_with_pc(g, score_matrix[, k])
    } else {
      others <- score_matrix[, -k, drop = FALSE]
      g_resid <- residuals(lm(g ~ others))
      pc_resid <- residuals(lm(score_matrix[, k] ~ others))
      out[k] <- cor_with_pc(g_resid, pc_resid)
    }
  }

  out
}

build_candidate_table <- function(Cmat, PCmat, FDR_simple, FDR_partial, core_genes, combo) {
  pass_partial <- rowSums(
    FDR_partial < combo$alpha & abs(PCmat) >= combo$min_abs_partial_cor,
    na.rm = TRUE
  ) >= combo$m_and
  pass_simple <- rowSums(
    FDR_simple < combo$alpha_soft & abs(Cmat) >= combo$min_abs_simple_cor,
    na.rm = TRUE
  ) >= 1
  keep <- pass_partial & pass_simple

  candidate_genes <- setdiff(rownames(Cmat)[keep], core_genes)
  expanded_genes <- unique(c(core_genes, candidate_genes))

  candidate_table <- data.frame(
    gene_symbol = expanded_genes,
    is_core = expanded_genes %in% core_genes,
    max_abs_simple_cor = apply(abs(Cmat[expanded_genes, , drop = FALSE]), 1, max, na.rm = TRUE),
    max_abs_partial_cor = apply(abs(PCmat[expanded_genes, , drop = FALSE]), 1, max, na.rm = TRUE),
    min_fdr_simple = apply(FDR_simple[expanded_genes, , drop = FALSE], 1, min, na.rm = TRUE),
    min_fdr_partial = apply(FDR_partial[expanded_genes, , drop = FALSE], 1, min, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  candidate_table$evidence_score <- with(
    candidate_table,
    max_abs_simple_cor * max_abs_partial_cor * (-log10(pmax(min_fdr_partial, 1e-300)))
  )
  candidate_table <- candidate_table[order(!candidate_table$is_core, -candidate_table$evidence_score), , drop = FALSE]

  list(
    genes = expanded_genes,
    table = candidate_table,
    n_genes = length(expanded_genes),
    n_non_core = sum(!expanded_genes %in% core_genes)
  )
}

select_parameter_combo <- function(selection_grid, target_min, target_max) {
  in_range <- selection_grid$n_genes >= target_min & selection_grid$n_genes <= target_max
  if (any(in_range)) {
    return(selection_grid[which(in_range)[1], , drop = FALSE])
  }

  distance_to_range <- ifelse(
    selection_grid$n_genes < target_min,
    target_min - selection_grid$n_genes,
    selection_grid$n_genes - target_max
  )
  selection_grid[which.min(distance_to_range), , drop = FALSE]
}

cap_expanded_table <- function(candidate_table, core_genes, target_max) {
  core_table <- candidate_table[candidate_table$is_core, , drop = FALSE]
  non_core_table <- candidate_table[!candidate_table$is_core, , drop = FALSE]

  n_keep_non_core <- max(target_max - nrow(core_table), 0)
  if (nrow(non_core_table) > n_keep_non_core) {
    non_core_table <- non_core_table[seq_len(n_keep_non_core), , drop = FALSE]
  }

  final_table <- rbind(core_table, non_core_table)
  final_table <- final_table[match(c(core_genes, setdiff(final_table$gene_symbol, core_genes)), final_table$gene_symbol), , drop = FALSE]
  rownames(final_table) <- NULL
  final_table
}

log_step("Loading temporal PCGS session object")
if (!file.exists(params$input_rds)) {
  stop(
    "Required input not found: ", params$input_rds, "\n",
    "Run PCGS/run_pcgs_temporal.R first, or update params$input_rds to the correct RDS path."
  )
}

pcgs_temporal <- readRDS(params$input_rds)
require_named_elements(
  pcgs_temporal,
  c("params", "expr_gene_all", "expr_core", "core_genes", "group_temporal"),
  "pcgs_temporal"
)

expr_gene_all <- as.matrix(pcgs_temporal$expr_gene_all)
expr_core <- as.matrix(pcgs_temporal$expr_core)
core_genes <- unique(as.character(pcgs_temporal$core_genes))
core_genes <- core_genes[!is.na(core_genes) & core_genes != ""]
group_temporal <- pcgs_temporal$group_temporal
pcgs_params <- pcgs_temporal$params

if (is.null(rownames(expr_gene_all)) || is.null(colnames(expr_gene_all))) {
  stop("expr_gene_all must be a genes x samples matrix with row and column names.")
}
if (is.null(rownames(expr_core)) || is.null(colnames(expr_core))) {
  stop("expr_core must be a genes x samples matrix with row and column names.")
}

common_samples <- intersect(colnames(expr_gene_all), colnames(expr_core))
if (length(common_samples) < 5) {
  stop("Too few overlapping temporal samples between expr_gene_all and expr_core.")
}

expr_gene_all <- expr_gene_all[, common_samples, drop = FALSE]
expr_core <- expr_core[, common_samples, drop = FALSE]
core_genes <- intersect(core_genes, rownames(expr_gene_all))
core_genes <- intersect(core_genes, rownames(expr_core))

if (length(core_genes) < 10) {
  stop("Fewer than 10 mapped core genes are available for temporal PCA.")
}

log_step(sprintf("Temporal matrix loaded: %d genes x %d samples", nrow(expr_gene_all), ncol(expr_gene_all)))
log_step(sprintf("Mapped temporal core genes: %d", length(core_genes)))

log_step("Running PCA on the temporal core-gene matrix")
X_core <- t(expr_core[core_genes, , drop = FALSE])
X_core <- safe_scale(X_core)

pca <- prcomp(X_core, center = TRUE, scale. = TRUE)
variance_explained <- summary(pca)$importance[2, ]
K_fixed <- which(cumsum(variance_explained) >= params$pc_cumvar)[1]
K_fixed <- max(params$min_pcs, min(K_fixed, params$max_pcs, ncol(pca$x)))

pc_idx <- seq_len(K_fixed)
if (params$drop_pc1) {
  pc_idx <- setdiff(pc_idx, 1)
}
if (length(pc_idx) == 0) {
  stop("No principal components were retained after applying the PCA settings.")
}

scores <- pca$x[, pc_idx, drop = FALSE]
log_step(sprintf("Retained %d principal components for expansion", ncol(scores)))

log_step("Computing simple and partial correlations against temporal PCA scores")
expr_gene_all <- expr_gene_all[, rownames(scores), drop = FALSE]
X_all <- safe_scale(t(expr_gene_all))

Cmat <- sapply(seq_len(ncol(scores)), function(k) {
  apply(X_all, 2, function(g) cor_with_pc(g, scores[, k]))
})
Cmat <- as.matrix(Cmat)
colnames(Cmat) <- colnames(scores)
rownames(Cmat) <- colnames(X_all)

n_samples <- nrow(scores)
P_simple <- apply(Cmat, 2, function(r) cor_to_p(r, n_samples))
P_simple <- as.matrix(P_simple)
rownames(P_simple) <- rownames(Cmat)
colnames(P_simple) <- colnames(Cmat)

FDR_simple <- apply(P_simple, 2, p.adjust, method = "BH")
FDR_simple <- as.matrix(FDR_simple)
rownames(FDR_simple) <- rownames(Cmat)
colnames(FDR_simple) <- colnames(Cmat)

PCmat <- t(apply(X_all, 2, partial_cor_onegene, score_matrix = scores))
colnames(PCmat) <- colnames(scores)
rownames(PCmat) <- colnames(X_all)

df_partial <- max(n_samples - (ncol(scores) - 1) - 2, 5)
P_partial <- apply(PCmat, 2, function(r) cor_to_p(r, df_partial + 2))
P_partial <- as.matrix(P_partial)
rownames(P_partial) <- rownames(PCmat)
colnames(P_partial) <- colnames(PCmat)

FDR_partial <- apply(P_partial, 2, p.adjust, method = "BH")
FDR_partial <- as.matrix(FDR_partial)
rownames(FDR_partial) <- rownames(PCmat)
colnames(FDR_partial) <- colnames(PCmat)

log_step("Searching statistically supported parameter combinations")
selection_grid <- expand.grid(
  m_and = params$m_and_grid,
  alpha = params$alpha_grid,
  alpha_soft = params$alpha_soft_grid,
  min_abs_simple_cor = params$min_abs_simple_cor_grid,
  min_abs_partial_cor = params$min_abs_partial_cor_grid,
  stringsAsFactors = FALSE
)
selection_grid <- selection_grid[order(
  -selection_grid$m_and,
  selection_grid$alpha,
  selection_grid$alpha_soft,
  -selection_grid$min_abs_partial_cor,
  -selection_grid$min_abs_simple_cor
), , drop = FALSE]

candidate_results <- vector("list", nrow(selection_grid))
for (i in seq_len(nrow(selection_grid))) {
  combo <- selection_grid[i, , drop = FALSE]
  candidate_results[[i]] <- build_candidate_table(
    Cmat = Cmat,
    PCmat = PCmat,
    FDR_simple = FDR_simple,
    FDR_partial = FDR_partial,
    core_genes = core_genes,
    combo = combo
  )
  selection_grid$n_genes[i] <- candidate_results[[i]]$n_genes
  selection_grid$n_non_core[i] <- candidate_results[[i]]$n_non_core
}

best_combo <- select_parameter_combo(selection_grid, params$target_min, params$target_max)
best_idx <- which(
  selection_grid$m_and == best_combo$m_and &
    selection_grid$alpha == best_combo$alpha &
    selection_grid$alpha_soft == best_combo$alpha_soft &
    selection_grid$min_abs_simple_cor == best_combo$min_abs_simple_cor &
    selection_grid$min_abs_partial_cor == best_combo$min_abs_partial_cor
)[1]

best_result <- candidate_results[[best_idx]]
selected_within_range <- best_result$n_genes >= params$target_min && best_result$n_genes <= params$target_max
applied_cap <- FALSE

if (best_result$n_genes > params$target_max) {
  final_results <- cap_expanded_table(best_result$table, core_genes, params$target_max)
  applied_cap <- TRUE
} else {
  final_results <- best_result$table
}

expanded_genes <- final_results$gene_symbol
expr_temporal_expanded <- expr_gene_all[expanded_genes, , drop = FALSE]

sample_metadata_temporal <- data.frame(
  sample = colnames(expr_temporal_expanded),
  group_temporal = as.character(group_temporal),
  stringsAsFactors = FALSE
)
if (nrow(sample_metadata_temporal) != ncol(expr_temporal_expanded)) {
  sample_metadata_temporal$group_temporal <- as.character(group_temporal)[seq_len(ncol(expr_temporal_expanded))]
}

PCA_118553_results <- final_results
PCA_118553_results$rank <- seq_len(nrow(PCA_118553_results))
PCA_118553_results <- PCA_118553_results[, c(
  "rank", "gene_symbol", "is_core", "max_abs_simple_cor",
  "max_abs_partial_cor", "min_fdr_simple", "min_fdr_partial", "evidence_score"
)]

log_step(sprintf(
  paste(
    "Selected parameter set:",
    "m_and=%d, alpha=%.3f, alpha_soft=%.3f,",
    "min_abs_simple_cor=%.2f, min_abs_partial_cor=%.2f"
  ),
  best_combo$m_and,
  best_combo$alpha,
  best_combo$alpha_soft,
  best_combo$min_abs_simple_cor,
  best_combo$min_abs_partial_cor
))
log_step(sprintf("Expanded temporal feature space: %d genes", nrow(PCA_118553_results)))
if (applied_cap) {
  log_step(sprintf("Applied evidence-based cap at %d genes", params$target_max))
} else if (selected_within_range) {
  log_step("No hard cap was applied because the selected parameter set already fell within the target range")
} else {
  log_step("No hard cap was required, but the selected statistically supported result remained outside the preferred range")
}

write.csv(PCA_118553_results, file.path(params$out_dir, "PCA_118553_results.csv"), row.names = FALSE)
writeLines(expanded_genes, file.path(params$out_dir, "temporal_pca_expanded_genes.txt"))
write.csv(
  data.frame(gene_symbol = rownames(expr_temporal_expanded), expr_temporal_expanded, check.names = FALSE),
  file.path(params$out_dir, "temporal_pca_expanded_expression_gene_x_sample.csv"),
  row.names = FALSE
)
write.csv(
  sample_metadata_temporal,
  file.path(params$out_dir, "sample_metadata_temporal_for_pca.csv"),
  row.names = FALSE
)
write.csv(
  selection_grid,
  file.path(params$out_dir, "temporal_pca_parameter_scan.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    params = params,
    pcgs_params = pcgs_params,
    selected_combo = best_combo,
    selected_within_range = selected_within_range,
    applied_cap = applied_cap,
    selection_grid = selection_grid,
    core_genes = core_genes,
    expanded_genes = expanded_genes,
    sample_metadata_temporal = sample_metadata_temporal,
    pca = pca,
    scores = scores,
    Cmat = Cmat,
    PCmat = PCmat,
    FDR_simple = FDR_simple,
    FDR_partial = FDR_partial,
    PCA_118553_results = PCA_118553_results
  ),
  file = file.path(params$out_dir, "temporal_pca_expansion_diagnostics.rds")
)

log_step("Temporal PCA expansion completed")
