#!/usr/bin/env Rscript

# enet_train_test.R
# Elastic Net model training and external validation.
#
# This script assumes that the training and test cohorts have already been
# preprocessed independently and exported as sample x gene matrices.
# The training cohort is used exclusively for repeated cross-validation,
# model fitting, and alpha selection. The external test cohort is used
# exclusively for out-of-sample validation.

suppressPackageStartupMessages({
  library(glmnet)
  library(pROC)
})

params <- list(
  out_dir = "outputs",
  train_expression_csv_candidates = c(
    file.path("outputs", "temporal_elasticnet_training_expression_sample_x_gene.csv"),
    file.path("Code_test", "outputs", "temporal_elasticnet_training_expression_sample_x_gene.csv")
  ),
  train_labels_csv_candidates = c(
    file.path("outputs", "temporal_elasticnet_training_labels.csv"),
    file.path("Code_test", "outputs", "temporal_elasticnet_training_labels.csv")
  ),
  test_expression_csv_candidates = c(
    file.path("outputs", "temporal_elasticnet_test_expression_sample_x_gene.csv"),
    file.path("Code_test", "outputs", "temporal_elasticnet_test_expression_sample_x_gene.csv")
  ),
  test_labels_csv_candidates = c(
    file.path("outputs", "temporal_elasticnet_test_labels.csv"),
    file.path("Code_test", "outputs", "temporal_elasticnet_test_labels.csv")
  ),
  alpha_grid = c(0.4, 0.55, 0.7, 0.85),
  nfolds = 5,
  n_repeats = 30,
  seed = 20260313
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

resolve_input_path <- function(path_candidates, label) {
  existing_paths <- path_candidates[file.exists(path_candidates)]
  if (length(existing_paths) == 0) {
    stop(
      "Required ", label, " file was not found. Checked paths: ",
      paste(path_candidates, collapse = ", ")
    )
  }
  existing_paths[1]
}

read_expression_matrix <- function(path) {
  x <- read.csv(path, row.names = 1, check.names = FALSE)
  x <- as.matrix(x)
  mode(x) <- "numeric"
  x
}

read_binary_labels <- function(path, group_column) {
  label_table <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required_columns <- c("sample_id", group_column, "y_binary")
  if (!all(required_columns %in% colnames(label_table))) {
    stop(
      "Label file is missing required columns: ",
      paste(setdiff(required_columns, colnames(label_table)), collapse = ", ")
    )
  }
  label_table
}

align_expression_and_labels <- function(x, label_table) {
  common_samples <- intersect(rownames(x), label_table$sample_id)
  if (length(common_samples) == 0) {
    stop("No overlapping samples were found between the expression matrix and the label table.")
  }

  x <- x[common_samples, , drop = FALSE]
  label_table <- label_table[match(common_samples, label_table$sample_id), , drop = FALSE]

  if (!identical(rownames(x), label_table$sample_id)) {
    stop("Sample ordering failed during matrix/label alignment.")
  }

  list(
    x = x,
    y = as.numeric(label_table$y_binary),
    labels = label_table
  )
}

save_selected_genes <- function(fit, out_prefix) {
  beta <- as.matrix(coef(fit))
  beta <- beta[setdiff(rownames(beta), "(Intercept)"), , drop = FALSE]
  beta <- beta[beta[, 1] != 0, , drop = FALSE]

  selected_gene_table <- data.frame(
    gene_symbol = rownames(beta),
    coefficient = as.numeric(beta[, 1]),
    direction = ifelse(beta[, 1] > 0, "up_in_AD", "down_in_AD"),
    stringsAsFactors = FALSE
  )
  selected_gene_table <- selected_gene_table[order(-abs(selected_gene_table$coefficient), selected_gene_table$gene_symbol), , drop = FALSE]

  write.csv(
    selected_gene_table,
    file.path(params$out_dir, paste0(out_prefix, "_selected_genes.csv")),
    row.names = FALSE
  )
  write.csv(
    data.frame(gene_symbol = selected_gene_table$gene_symbol, stringsAsFactors = FALSE),
    file.path(params$out_dir, paste0(out_prefix, "_selected_gene_list.csv")),
    row.names = FALSE
  )
  write.csv(
    selected_gene_table[seq_len(min(30, nrow(selected_gene_table))), , drop = FALSE],
    file.path(params$out_dir, paste0(out_prefix, "_top30_for_interpretation.csv")),
    row.names = FALSE
  )
  write.csv(
    data.frame(
      gene_symbol = selected_gene_table$gene_symbol[seq_len(min(30, nrow(selected_gene_table)))],
      stringsAsFactors = FALSE
    ),
    file.path(params$out_dir, paste0(out_prefix, "_top30_gene_list_for_interpretation.csv")),
    row.names = FALSE
  )

  selected_gene_table
}

fit_alpha_grid <- function(x, y, alpha_grid, nfolds, n_repeats, seed) {
  cv_results <- vector("list", length(alpha_grid))
  repeat_level_results <- vector("list", length(alpha_grid))
  summary_table <- data.frame(
    alpha = alpha_grid,
    lambda_min_mean = NA_real_,
    lambda_1se_mean = NA_real_,
    cv_auc_lambda_min_mean = NA_real_,
    cv_auc_lambda_1se_mean = NA_real_,
    cv_auc_lambda_1se_sd = NA_real_,
    n_selected_genes_lambda_1se_mean = NA_real_,
    n_selected_genes_lambda_1se_median = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(alpha_grid)) {
    alpha_value <- alpha_grid[i]
    alpha_repeat_results <- vector("list", n_repeats)

    for (repeat_idx in seq_len(n_repeats)) {
      set.seed(seed + i * 1000 + repeat_idx)
      foldid <- sample(rep(seq_len(nfolds), length.out = nrow(x)))

      cv_fit <- cv.glmnet(
        x = x,
        y = y,
        family = "binomial",
        alpha = alpha_value,
        foldid = foldid,
        nfolds = nfolds,
        type.measure = "auc",
        standardize = FALSE,
        keep = TRUE
      )

      idx_lambda_min <- which.min(abs(log(cv_fit$lambda) - log(cv_fit$lambda.min)))
      idx_lambda_1se <- which.min(abs(log(cv_fit$lambda) - log(cv_fit$lambda.1se)))
      pred_lambda_min <- as.numeric(cv_fit$fit.preval[, idx_lambda_min])
      pred_lambda_1se <- as.numeric(cv_fit$fit.preval[, idx_lambda_1se])

      auc_lambda_min <- as.numeric(auc(roc(y, pred_lambda_min, quiet = TRUE, direction = "<")))
      auc_lambda_1se <- as.numeric(auc(roc(y, pred_lambda_1se, quiet = TRUE, direction = "<")))

      fit_lambda_1se <- glmnet(
        x = x,
        y = y,
        family = "binomial",
        alpha = alpha_value,
        lambda = cv_fit$lambda.1se,
        standardize = FALSE
      )
      beta_lambda_1se <- as.matrix(coef(fit_lambda_1se))
      n_selected_genes_lambda_1se <- sum(
        rownames(beta_lambda_1se) != "(Intercept)" & beta_lambda_1se[, 1] != 0
      )

      alpha_repeat_results[[repeat_idx]] <- data.frame(
        repeat_id = repeat_idx,
        alpha = alpha_value,
        lambda_min = cv_fit$lambda.min,
        lambda_1se = cv_fit$lambda.1se,
        cv_auc_lambda_min = auc_lambda_min,
        cv_auc_lambda_1se = auc_lambda_1se,
        n_selected_genes_lambda_1se = n_selected_genes_lambda_1se,
        stringsAsFactors = FALSE
      )
    }

    alpha_repeat_table <- do.call(rbind, alpha_repeat_results)
    repeat_level_results[[i]] <- alpha_repeat_table

    cv_fit_full <- cv.glmnet(
      x = x,
      y = y,
      family = "binomial",
      alpha = alpha_value,
      nfolds = nfolds,
      type.measure = "auc",
      standardize = FALSE,
      keep = TRUE
    )

    cv_results[[i]] <- cv_fit_full
    summary_table$lambda_min_mean[i] <- mean(alpha_repeat_table$lambda_min)
    summary_table$lambda_1se_mean[i] <- mean(alpha_repeat_table$lambda_1se)
    summary_table$cv_auc_lambda_min_mean[i] <- mean(alpha_repeat_table$cv_auc_lambda_min)
    summary_table$cv_auc_lambda_1se_mean[i] <- mean(alpha_repeat_table$cv_auc_lambda_1se)
    summary_table$cv_auc_lambda_1se_sd[i] <- stats::sd(alpha_repeat_table$cv_auc_lambda_1se)
    summary_table$n_selected_genes_lambda_1se_mean[i] <- mean(alpha_repeat_table$n_selected_genes_lambda_1se)
    summary_table$n_selected_genes_lambda_1se_median[i] <- stats::median(alpha_repeat_table$n_selected_genes_lambda_1se)
  }

  best_auc <- max(summary_table$cv_auc_lambda_1se_mean)
  auc_tolerance <- summary_table$cv_auc_lambda_1se_sd[which.max(summary_table$cv_auc_lambda_1se_mean)]
  candidate_idx <- which(summary_table$cv_auc_lambda_1se_mean >= best_auc - auc_tolerance)
  candidate_table <- summary_table[candidate_idx, , drop = FALSE]
  best_idx <- candidate_idx[order(
    candidate_table$n_selected_genes_lambda_1se_median,
    -candidate_table$alpha
  )[1]]

  list(
    best_idx = best_idx,
    best_alpha = alpha_grid[best_idx],
    best_cvfit = cv_results[[best_idx]],
    summary_table = summary_table,
    repeat_level_table = do.call(rbind, repeat_level_results)
  )
}

log_step("Loading preprocessed training and external test matrices")
train_expression_csv <- resolve_input_path(params$train_expression_csv_candidates, "training expression")
train_labels_csv <- resolve_input_path(params$train_labels_csv_candidates, "training labels")
test_expression_csv <- resolve_input_path(params$test_expression_csv_candidates, "test expression")
test_labels_csv <- resolve_input_path(params$test_labels_csv_candidates, "test labels")

log_step(sprintf("Training expression file: %s", train_expression_csv))
log_step(sprintf("Training labels file: %s", train_labels_csv))
log_step(sprintf("Test expression file: %s", test_expression_csv))
log_step(sprintf("Test labels file: %s", test_labels_csv))

X_train <- read_expression_matrix(train_expression_csv)
X_test <- read_expression_matrix(test_expression_csv)
train_labels <- read_binary_labels(train_labels_csv, "group_temporal_training")
test_labels <- read_binary_labels(test_labels_csv, "group_temporal_test")

aligned_train <- align_expression_and_labels(X_train, train_labels)
aligned_test <- align_expression_and_labels(X_test, test_labels)

if (!identical(colnames(aligned_train$x), colnames(aligned_test$x))) {
  stop("Training and test gene spaces are not aligned. Re-run preprocessing to enforce the same gene order.")
}

X_train <- aligned_train$x
X_test <- aligned_test$x
y_train <- aligned_train$y
y_test <- aligned_test$y

if (!all(y_train %in% c(0, 1)) || !all(y_test %in% c(0, 1))) {
  stop("Binary outcome labels must be coded as 0 or 1.")
}

log_step(sprintf(
  "Elastic Net input prepared: training %d x %d, test %d x %d",
  nrow(X_train), ncol(X_train), nrow(X_test), ncol(X_test)
))

set.seed(params$seed)

log_step("Selecting alpha by repeated 5-fold cross-validation on the training cohort only")
alpha_search <- fit_alpha_grid(
  x = X_train,
  y = y_train,
  alpha_grid = params$alpha_grid,
  nfolds = min(params$nfolds, nrow(X_train)),
  n_repeats = params$n_repeats,
  seed = params$seed
)
cvfit <- alpha_search$best_cvfit

final_fit <- glmnet(
  x = X_train,
  y = y_train,
  family = "binomial",
  alpha = alpha_search$best_alpha,
  lambda = cvfit$lambda.1se,
  standardize = FALSE
)

idx_lambda_1se <- which.min(abs(log(cvfit$lambda) - log(cvfit$lambda.1se)))
train_cv_predictions <- as.numeric(cvfit$fit.preval[, idx_lambda_1se])
train_cv_roc <- roc(y_train, train_cv_predictions, quiet = TRUE, direction = "<")
train_cv_auc <- as.numeric(auc(train_cv_roc))

test_predictions <- as.numeric(predict(final_fit, newx = X_test, type = "response"))
test_roc <- roc(y_test, test_predictions, quiet = TRUE, direction = "<")
test_auc <- as.numeric(auc(test_roc))

log_step(sprintf(
  "Selected alpha = %.2f | Cross-validated training AUC = %.4f | External test AUC = %.4f",
  alpha_search$best_alpha, train_cv_auc, test_auc
))

grDevices::pdf(file.path(params$out_dir, "temporal_elasticnet_coefficient_path.pdf"), width = 9, height = 7)
plot(
  cvfit$glmnet.fit,
  xvar = "lambda",
  main = "Elastic Net Coefficient Path\nTraining Cohort Only\n"
)
grDevices::dev.off()

grDevices::pdf(file.path(params$out_dir, "temporal_elasticnet_external_roc.pdf"), width = 8, height = 7)
plot(
  test_roc,
  main = sprintf("External Validation ROC Curve\nAUC = %.4f\n", test_auc),
  col = "#1f4e79",
  lwd = 2
)
abline(a = 0, b = 1, lty = 2, col = "grey60")
grDevices::dev.off()

grDevices::pdf(file.path(params$out_dir, "temporal_elasticnet_training_cv_roc.pdf"), width = 8, height = 7)
plot(
  train_cv_roc,
  main = sprintf("Cross-Validated Training ROC Curve\nAUC = %.4f\n", train_cv_auc),
  col = "#7a0019",
  lwd = 2
)
abline(a = 0, b = 1, lty = 2, col = "grey60")
grDevices::dev.off()

selected_gene_table <- save_selected_genes(final_fit, "temporal_elasticnet_lambda_1se")

write.csv(
  data.frame(
    sample_id = rownames(X_test),
    y_true = y_test,
    y_pred_probability = test_predictions,
    stringsAsFactors = FALSE
  ),
  file.path(params$out_dir, "temporal_elasticnet_external_predictions.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    sample_id = rownames(X_train),
    y_true = y_train,
    y_pred_probability = train_cv_predictions,
    stringsAsFactors = FALSE
  ),
  file.path(params$out_dir, "temporal_elasticnet_training_cv_predictions.csv"),
  row.names = FALSE
)

write.csv(
  alpha_search$summary_table,
  file.path(params$out_dir, "temporal_elasticnet_alpha_grid_search.csv"),
  row.names = FALSE
)
write.csv(
  alpha_search$repeat_level_table,
  file.path(params$out_dir, "temporal_elasticnet_alpha_grid_search_repeats.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    parameter = c(
      "train_expression_csv",
      "train_labels_csv",
      "test_expression_csv",
      "test_labels_csv",
      "alpha_grid",
      "selected_alpha",
      "nfolds_requested",
      "nfolds_used",
      "n_repeats",
      "seed",
      "lambda_min",
      "lambda_1se",
      "n_genes_train_and_test",
      "n_selected_genes_lambda_1se",
      "training_cv_auc",
      "external_test_auc"
    ),
    value = c(
      train_expression_csv,
      train_labels_csv,
      test_expression_csv,
      test_labels_csv,
      paste(params$alpha_grid, collapse = " | "),
      alpha_search$best_alpha,
      params$nfolds,
      min(params$nfolds, nrow(X_train)),
      params$n_repeats,
      params$seed,
      cvfit$lambda.min,
      cvfit$lambda.1se,
      ncol(X_train),
      nrow(selected_gene_table),
      train_cv_auc,
      test_auc
    ),
    stringsAsFactors = FALSE
  ),
  file.path(params$out_dir, "temporal_elasticnet_model_parameters.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    params = params,
    alpha_search = alpha_search,
    common_genes = colnames(X_train),
    cvfit = cvfit,
    final_fit = final_fit,
    train_cv_auc = train_cv_auc,
    test_auc = test_auc,
    selected_gene_table = selected_gene_table
  ),
  file = file.path(params$out_dir, "temporal_elasticnet_model_diagnostics.rds")
)

# ------------------------------------------------------------------------------
# Additional plot: Horizontal barplot for selected gene regulation direction
#   - logFC-based direction on the training cohort only
#   - Up in AD   : mean(training AD) - mean(training Control) > 0
#   - Down in AD : mean(training AD) - mean(training Control) < 0
#   - External test cohort remains validation-only
# ------------------------------------------------------------------------------

log_step("Saving selected gene regulation direction barplot (training cohort logFC-based)")

selected_gene_direction_barplot_file <- file.path(
  params$out_dir,
  "temporal_elasticnet_selected_gene_direction_barplot.pdf"
)

if (!exists("selected_gene_table") || nrow(selected_gene_table) == 0) {
  warning("selected_gene_table is empty; skipping selected gene direction barplot.")
} else {
  selected_genes <- selected_gene_table$gene_symbol
  selected_genes <- intersect(selected_genes, colnames(X_train))
  
  if (length(selected_genes) == 0) {
    warning("No selected genes were found in X_train; skipping selected gene direction barplot.")
  } else {
    ad_idx_train <- which(y_train == 1)
    control_idx_train <- which(y_train == 0)
    
    if (length(ad_idx_train) == 0 || length(control_idx_train) == 0) {
      warning("Both AD and Control samples are required in X_train to compute training logFC; skipping selected gene direction barplot.")
    } else {
      plot_table <- data.frame(
        gene_symbol = selected_genes,
        mean_AD_training = apply(X_train[ad_idx_train, selected_genes, drop = FALSE], 2, mean, na.rm = TRUE),
        mean_Control_training = apply(X_train[control_idx_train, selected_genes, drop = FALSE], 2, mean, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      
      plot_table$logFC_training <- plot_table$mean_AD_training - plot_table$mean_Control_training
      plot_table <- plot_table[is.finite(plot_table$logFC_training) & plot_table$logFC_training != 0, , drop = FALSE]
      
      if (nrow(plot_table) == 0) {
        warning("All computed training logFC values were zero or non-finite; skipping selected gene direction barplot.")
      } else {
        plot_table$direction <- ifelse(plot_table$logFC_training > 0, "up_in_AD", "down_in_AD")
        plot_table$plot_value <- plot_table$logFC_training
        
        plot_table <- merge(
          plot_table,
          selected_gene_table[, c("gene_symbol", "coefficient"), drop = FALSE],
          by = "gene_symbol",
          all.x = TRUE,
          sort = FALSE
        )
        
        plot_table <- plot_table[order(plot_table$plot_value), , drop = FALSE]
        
        bar_colors <- ifelse(
          plot_table$direction == "up_in_AD",
          "firebrick3",
          "royalblue3"
        )
        
        grDevices::pdf(
          selected_gene_direction_barplot_file,
          width = 12,
          height = max(7, nrow(plot_table) * 0.32 + 3),
          useDingbats = FALSE
        )
        
        par(mar = c(5, 14, 5, 3))
        bar_mid <- barplot(
          height = plot_table$plot_value,
          names.arg = plot_table$gene_symbol,
          horiz = TRUE,
          las = 1,
          col = bar_colors,
          border = NA,
          xlab = "Training cohort logFC (mean AD - mean Control)",
          main = "Regulation direction of Elastic Net-selected genes\nTraining cohort logFC-based\n",
          cex.names = if (nrow(plot_table) <= 40) 0.9 else if (nrow(plot_table) <= 80) 0.75 else 0.6
        )
        
        abline(v = 0, lty = 2, lwd = 1.2, col = "grey40")
        
        text(
          x = plot_table$plot_value,
          y = bar_mid,
          labels = sprintf("%.3f", plot_table$plot_value),
          pos = ifelse(plot_table$plot_value > 0, 4, 2),
          xpd = TRUE,
          cex = 0.8
        )
        
        legend(
          "topright",
          legend = c("Up in AD", "Down in AD"),
          fill = c("firebrick3", "royalblue3"),
          border = NA,
          bty = "n",
          cex = 0.95
        )
        
        box()
        grDevices::dev.off()
        
        write.csv(
          plot_table,
          file.path(params$out_dir, "temporal_elasticnet_selected_gene_direction_logFC_table.csv"),
          row.names = FALSE
        )
        
        log_step(sprintf(
          "Saved selected gene direction barplot to: %s",
          selected_gene_direction_barplot_file
        ))
      }
    }
  }
}

log_step("Elastic Net training and external validation completed")
