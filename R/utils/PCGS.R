# Utility helpers for the PCGS (Preclinical Core Gene Set) pipeline

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
})

pick_col <- function(df, patterns, required = TRUE) {
  cn <- colnames(df)
  hit <- NA_character_
  for (p in patterns) {
    idx <- which(str_detect(tolower(cn), p))
    if (length(idx) > 0) { hit <- cn[idx[1]]; break }
  }
  if (is.na(hit) && required) {
    stop(
      "Could not find a column matching patterns: ", paste(patterns, collapse = " | "), "\n",
      "Available columns: ", paste(colnames(df), collapse = ", ")
    )
  }
  if (is.na(hit)) return(NULL)
  hit
}

normalize_group <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  x0 <- str_replace_all(x0, "\\s+", "")
  dplyr::case_when(
    str_detect(x0, "control|ctl|normal|nonad|healthy") ~ "Control",
    str_detect(x0, "asym") ~ "AsymAD",
    str_detect(x0, "^ad$|alzheimer|symad|dementia") ~ "AD",
    TRUE ~ NA_character_
  )
}

clean_symbol <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == "" | is.na(x)] <- NA_character_
  # "A /// B", "A;B", "A,B" -> keep the first token
  x <- str_split_fixed(x, "\\s*///\\s*|\\s*;\\s*|\\s*,\\s*", 2)[, 1]
  x <- trimws(x)
  x[x == "" | is.na(x)] <- NA_character_
  x
}

pick_symbol_col_from_gpl <- function(gpl_tbl) {
  cn <- colnames(gpl_tbl)
  c0 <- tolower(cn)

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

choose_cor_threshold <- function(max_abs_cor, core_genes, cor_grid,
                                 target_min, target_max) {
  genes <- names(max_abs_cor)

  best_thr <- NA_real_
  best_n   <- NA_integer_

  for (thr in cor_grid) {
    expanded <- genes[which(max_abs_cor >= thr)]
    expanded <- union(expanded, core_genes)
    nn <- length(expanded)
    if (nn >= target_min && nn <= target_max) {
      best_thr <- thr
      best_n <- nn
      break
    }
  }

  if (is.na(best_thr)) {
    cand_n <- sapply(cor_grid, function(thr) {
      expanded <- genes[which(max_abs_cor >= thr)]
      length(union(expanded, core_genes))
    })
    dist <- ifelse(cand_n < target_min, target_min - cand_n,
                   ifelse(cand_n > target_max, cand_n - target_max, 0))
    idx <- which.min(dist)
    best_thr <- cor_grid[idx]
    best_n <- cand_n[idx]
  }

  list(threshold = best_thr, expanded_n = best_n)
}
