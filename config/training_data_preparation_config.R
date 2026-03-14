#!/usr/bin/env Rscript

training_data_preparation_config <- list(
  gse_id = "GSE132903",
  tissue_column_candidates = c(
    "tissue:ch1",
    "tissue",
    "brain region:ch1",
    "brain region",
    "region:ch1",
    "region"
  ),
  tissue_value_candidates = c(
    "Temporal_Cortex",
    "Temporal Cortex",
    "temporal cortex",
    "middle temporal gyrus",
    "Middle temporal gyrus",
    "MTG",
    "tissue: Temporal_Cortex",
    "tissue: Temporal Cortex"
  ),
  disease_column_candidates = c(
    "disease state:ch1",
    "disease state",
    "disease:ch1",
    "disease",
    "diagnosis:ch1",
    "diagnosis",
    "group:ch1",
    "group",
    "status:ch1",
    "status"
  ),
  control_label = "control",
  case_label = "AD",
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
  feature_space_csv = file.path("outputs", "WGCNA_118553_results.csv"),
  out_dir = "outputs",
  max_genes_final = 30,
  correlation_pruning_threshold = 0.90,
  univariate_fdr_threshold = 0.10,
  univariate_raw_p_threshold = 0.05,
  univariate_min_genes = 5
)
