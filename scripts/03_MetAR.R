# ================================================================
# 03_run_meta_analysis.R
# RNA-seq variable-expression meta-analysis pipeline
#
# Purpose:
#   1. Read log2 TMM-CPM matrices produced by preprocessing script 02.
#   2. Extract target gene expression values for each dataset.
#   3. Merge target gene expression with manually curated phenotype data.
#   4. Fit per-dataset variable-expression models.
#   5. Pool per-dataset effects using random-effects meta-analysis.
#   6. Generate diagnostic plots, publication tables, and summary outputs.
#
# Input requirements:
#   - inputs/tmm_data/raw_tmm_data/*_log2_tmm_cpm.csv
#   - inputs/pdata/proc_pdata/*_proc_seq_pdata.csv
#
# Notes:
#   - TMM normalisation and expression filtering are performed in script 02.
#   - This script treats input expression matrices as already log2 TMM-CPM.
# ================================================================


# ================================================================
# 0. Package setup
# ================================================================

required_packages <- c(
  "dplyr",
  "readr",
  "tibble",
  "tidyr",
  "stringr",
  "purrr",
  "forcats",
  "ggplot2",
  "ggtext",
  "patchwork",
  "metafor",
  "mgcv",
  "future",
  "future.apply",
  "here",
  "scales",
  "grid"
)

install_missing_packages <- function(packages) {
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) > 0) {
    install.packages(missing_packages, dependencies = TRUE)
  }

  invisible(lapply(packages, library, character.only = TRUE))
}

install_missing_packages(required_packages)


# ================================================================
# 1. User-adjustable analysis settings
# ================================================================

# Primary phenotype/metadata variable to analyse.
# This column must exist in each processed pData file:
#   inputs/pdata/proc_pdata/*_proc_seq_pdata.csv

variable_of_interest <- "age"

# User-facing label used in plots and tables.
variable_label <- "Age"

# Unit label for the variable of interest.
# Examples:
#   "years", "kg/m²", "nmol/L", "pmol/L"
variable_units <- "years"

# Scaling used when calculating interpretable expression change.
# For age, 10 gives percent change per decade.
# For other variables, set to 1 for percent change per 1 unit,
# or another biologically meaningful increment.
effect_scale <- 10

effect_scale_label <- ifelse(
  effect_scale == 1,
  paste0("per 1 ", variable_units),
  paste0("per ", effect_scale, " ", variable_units)
)

min_samples   <- 6       # Minimum samples required per dataset model
min_k         <- 5       # Minimum datasets required to retain a gene for meta-analysis/forest plot
edf_threshold <- 1.7     # GAM non-linearity threshold
aic_delta     <- 5       # GAM must improve AIC by this amount to replace LM
n_boot        <- 1000    # Bootstrap iterations for GAM-derived average slope SE
alpha         <- 0.05    # Significance threshold

parallel_workers <- max(1, future::availableCores() - 1)

future::plan(future::multisession, workers = parallel_workers)

# Jungle plot category settings
# These should match values in inputs/target_gene_list.csv -> category column.
# Edit these to suit your analysis.

jungle_groupings <- list(
  "All target genes" = NULL
)

jungle_facet_groupings <- list(
  "All target genes" = NULL
)


# ================================================================
# 2. Project directories
# ================================================================

project_root <- here::here()

dirs <- list(
  inputs              = file.path(project_root, "inputs"),
  outputs             = file.path(project_root, "outputs"),
  logs                = file.path(project_root, "outputs", "logs"),

  raw_tmm_data        = file.path(project_root, "inputs", "tmm_data", "raw_tmm_data"),
  proc_tmm_data       = file.path(project_root, "inputs", "tmm_data", "proc_tmm_data"),
  pre_proc_data       = file.path(project_root, "inputs", "pre_proc_data"),

  pdata               = file.path(project_root, "inputs", "pdata"),
  proc_pdata          = file.path(project_root, "inputs", "pdata", "proc_pdata"),

  figures             = file.path(project_root, "outputs", "figures"),
  baujat_plots        = file.path(project_root, "outputs", "figures", "baujat_plots"),
  funnel_plots        = file.path(project_root, "outputs", "figures", "funnel_plots"),
  forest_plots        = file.path(project_root, "outputs", "figures", "forest_plots"),
  scatterplot_grids   = file.path(project_root, "outputs", "figures", "scatterplot_grids"),
  jungle_plots        = file.path(project_root, "outputs", "figures", "jungle_plots"),

  tables              = file.path(project_root, "outputs", "tables"),
  gene_effects        = file.path(project_root, "outputs", "tables", "gene_effects"),
  meta_effects        = file.path(project_root, "outputs", "tables", "meta_effects"),
  publication_tables  = file.path(project_root, "outputs", "tables", "publication_tables"),
  loo_tables          = file.path(project_root, "outputs", "tables", "loo")
)

invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

message("Project root: ", project_root)
message("Directories created or already exist.")


# ================================================================
# 3. Helper functions
# ================================================================

get_dataset_id <- function(file_path) {
  basename(file_path) |>
    stringr::str_remove("_log2_tmm_cpm\\.csv$") |>
    stringr::str_remove("_tmm\\.csv$") |>
    stringr::str_remove("_raw_counts(_\\d+)?\\.csv$") |>
    stringr::str_remove("_raw_counts(_\\d+)?$") |>
    stringr::str_remove("_proc_tmm\\.csv$") |>
    stringr::str_remove("_pre_proc_data\\.csv$") |>
    stringr::str_remove("_proc_seq_pdata\\.csv$")
}

strip_ensembl_version <- function(x) {
  stringr::str_remove(as.character(x), "\\..*$")
}

format_beta <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    NA_character_,
    formatC(round(as.numeric(x), digits), format = "f", digits = digits)
  )
}

format_p <- function(p, digits = 4, threshold = 0.0001) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(
      p < threshold,
      paste0("<", threshold),
      formatC(as.numeric(p), format = "f", digits = digits)
    )
  )
}

format_p_inline <- function(p) {
  ifelse(
    is.na(p),
    "p = NA",
    ifelse(p < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p)))
  )
}

percent_change_scaled <- function(beta, scale = effect_scale) {
  # beta is the per-unit change on a log2 expression scale.
  # scale defines the reporting increment, for example 10 years for age.
  round((2^as.numeric(beta * scale) - 1) * 100, 3)
}

safe_mean_sd_range <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return("Not reported")
  }

  sprintf(
    "%.1f ± %.1f (%.1f–%.1f)",
    mean(x),
    stats::sd(x),
    min(x),
    max(x)
  )
}

first_matching_column <- function(df, patterns) {
  nms <- names(df)
  idx <- which(stringr::str_detect(tolower(nms), paste(patterns, collapse = "|")))

  if (length(idx) == 0) {
    NA_character_
  } else {
    nms[idx[1]]
  }
}

detect_most_common_value <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]

  if (length(x) == 0) {
    return("Not reported")
  }

  tab <- sort(table(x), decreasing = TRUE)
  names(tab)[1]
}

# ================================================================
# 4. Target gene list
# ================================================================

target_gene_file <- file.path(project_root, "inputs", "target_gene_list.csv")

if (!file.exists(target_gene_file)) {
  template_gene_table <- tibble::tibble(
    symbol = character(),
    ensembl_id = character(),
    category = character()
  )
  
  readr::write_csv(template_gene_table, target_gene_file)
  
  stop(
    "Created template target gene file: inputs/target_gene_list.csv\n\n",
    "MANUAL STEP REQUIRED:\n",
    "  1. Open inputs/target_gene_list.csv.\n",
    "  2. Add target genes using the columns:\n",
    "     symbol, ensembl_id, category\n",
    "  3. Save the file.\n",
    "  4. Re-run this script."
  )
}

gene_table <- readr::read_csv(
  target_gene_file,
  show_col_types = FALSE
) |>
  dplyr::mutate(
    symbol = stringr::str_trim(symbol),
    ensembl_id = stringr::str_trim(ensembl_id),
    category = stringr::str_trim(category)
  ) |>
  dplyr::filter(
    !is.na(symbol),
    symbol != "",
    !is.na(ensembl_id),
    ensembl_id != "",
    !is.na(category),
    category != ""
  ) |>
  dplyr::distinct(symbol, ensembl_id, .keep_all = TRUE)

required_gene_cols <- c("symbol", "ensembl_id", "category")

missing_gene_cols <- setdiff(required_gene_cols, colnames(gene_table))

if (length(missing_gene_cols) > 0) {
  stop(
    "inputs/target_gene_list.csv is missing required columns: ",
    paste(missing_gene_cols, collapse = ", ")
  )
}

if (nrow(gene_table) == 0) {
  stop(
    "No target genes found in inputs/target_gene_list.csv.\n\n",
    "Add genes using columns: symbol, ensembl_id, category"
  )
}

message("Loaded ", nrow(gene_table), " target genes from inputs/target_gene_list.csv")

# ================================================================
# 5. Validate required input files
# ================================================================

tmm_files <- list.files(
  dirs$raw_tmm_data,
  pattern = "(_log2_tmm_cpm|_tmm)\\.csv$",
  full.names = TRUE
)

pdata_files <- list.files(
  dirs$proc_pdata,
  pattern = "_proc_seq_pdata\\.csv$",
  full.names = TRUE
)

if (length(tmm_files) == 0) {
  stop(
    "No log2 TMM-CPM files found in inputs/tmm_data/raw_tmm_data/.\n",
    "Run 02_download_counts_and_tmm_normalise.R before this script."
  )
}

if (length(pdata_files) == 0) {
  stop(
    "No processed pData files found in inputs/pdata/proc_pdata/.\n",
    "Complete the manual phenotype curation step before this script."
  )
}


# Check that the variable of interest is present in at least one processed pData file.
variable_presence <- purrr::map_lgl(
  pdata_files,
  function(f) {
    pdata_names <- names(readr::read_csv(f, n_max = 0, show_col_types = FALSE))
    variable_of_interest %in% pdata_names
  }
)

if (!any(variable_presence)) {
  stop(
    "Variable of interest not found in any processed pData file: ",
    variable_of_interest, "\n",
    "Check that inputs/pdata/proc_pdata/*_proc_seq_pdata.csv contains a column named exactly: ",
    variable_of_interest
  )
}

message(
  "Variable of interest: ", variable_of_interest,
  " (", variable_label, ", ", variable_units, ")"
)


# ================================================================
# 6. Extract target gene expression from log2 TMM-CPM matrices
# ================================================================

extract_target_gene <- function(tmm_file, gene_symbol, ensembl_id, output_dir) {
  dataset_id <- get_dataset_id(tmm_file)

  tryCatch({
    expr <- read.csv(tmm_file, row.names = 1, check.names = FALSE)

    if (nrow(expr) == 0 || ncol(expr) == 0) {
      stop("Empty expression matrix.")
    }

    row_ids <- rownames(expr)
    row_ids_clean <- strip_ensembl_version(row_ids)
    target_clean <- strip_ensembl_version(ensembl_id)

    idx_ensembl <- which(row_ids == ensembl_id | row_ids_clean == target_clean)
    idx_symbol  <- which(row_ids == gene_symbol)

    idx <- if (length(idx_ensembl) > 0) {
      idx_ensembl[1]
    } else if (length(idx_symbol) > 0) {
      idx_symbol[1]
    } else {
      integer(0)
    }

    if (length(idx) == 0) {
      stop("Target gene not found.")
    }

    matched <- expr[idx, , drop = FALSE]

    out <- tibble::tibble(
      sample = colnames(matched),
      target_gene_TMM = as.numeric(matched[1, ])
    )

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    output_file <- file.path(output_dir, paste0(dataset_id, "_proc_tmm.csv"))
    readr::write_csv(out, output_file)

    tibble::tibble(
      dataset = dataset_id,
      gene = gene_symbol,
      status = "success",
      matched_id = rownames(matched)[1],
      output_file = output_file,
      reason = NA_character_
    )
  }, error = function(e) {
    tibble::tibble(
      dataset = dataset_id,
      gene = gene_symbol,
      status = "failed",
      matched_id = NA_character_,
      output_file = NA_character_,
      reason = conditionMessage(e)
    )
  })
}

message("\nExtracting target gene expression values...")

target_extraction_log <- purrr::map_dfr(seq_len(nrow(gene_table)), function(i) {
  gene_symbol <- gene_table$symbol[i]
  ensembl_id  <- gene_table$ensembl_id[i]

  message("  Gene: ", gene_symbol)

  gene_output_dir <- file.path(dirs$proc_tmm_data, gene_symbol)

  gene_log <- purrr::map_dfr(
    tmm_files,
    extract_target_gene,
    gene_symbol = gene_symbol,
    ensembl_id = ensembl_id,
    output_dir = gene_output_dir
  )

  gene_log_dir <- file.path(dirs$logs, gene_symbol)
  dir.create(gene_log_dir, recursive = TRUE, showWarnings = FALSE)

  readr::write_csv(
    gene_log,
    file.path(gene_log_dir, "target_gene_extraction_log.csv")
  )

  gene_log
})

readr::write_csv(
  target_extraction_log,
  file.path(dirs$logs, "target_gene_extraction_log_all_genes.csv")
)


# ================================================================
# 7. Merge target gene expression with processed phenotype data
# ================================================================

merge_expression_with_pdata <- function(gene_symbol, pdata_file) {
  dataset_id <- get_dataset_id(pdata_file)
  expr_file <- file.path(dirs$proc_tmm_data, gene_symbol, paste0(dataset_id, "_proc_tmm.csv"))

  if (!file.exists(expr_file)) {
    return(tibble::tibble(
      gene = gene_symbol,
      dataset = dataset_id,
      status = "no_matching_expression_file",
      merged_rows = NA_integer_,
      excluded_rows = NA_integer_
    ))
  }

  tryCatch({
    pdata <- readr::read_csv(pdata_file, show_col_types = FALSE)
    expr_data <- readr::read_csv(expr_file, show_col_types = FALSE)

    if (!"sample" %in% colnames(expr_data)) {
      stop("Expression file is missing sample column.")
    }

    tmm_samples <- as.character(expr_data$sample)

    match_counts <- vapply(
      names(pdata),
      function(col) sum(as.character(pdata[[col]]) %in% tmm_samples, na.rm = TRUE),
      integer(1)
    )

    best_col <- names(match_counts)[which.max(match_counts)]
    max_matches <- max(match_counts)

    if (is.na(best_col) || max_matches == 0) {
      return(tibble::tibble(
        gene = gene_symbol,
        dataset = dataset_id,
        status = "no_matching_samples",
        merged_rows = 0L,
        excluded_rows = nrow(expr_data)
      ))
    }

    pdata <- pdata |>
      dplyr::mutate(merge_id = as.character(.data[[best_col]]))

    expr_data <- expr_data |>
      dplyr::mutate(merge_id = as.character(sample))

    merged_data <- dplyr::inner_join(pdata, expr_data, by = "merge_id")

    if (variable_of_interest %in% colnames(merged_data)) {
      merged_data <- merged_data |>
        dplyr::mutate(
          "{variable_of_interest}" := suppressWarnings(as.numeric(.data[[variable_of_interest]]))
        )
    }

    out_dir <- file.path(dirs$pre_proc_data, gene_symbol)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    out_file <- file.path(out_dir, paste0(dataset_id, "_pre_proc_data.csv"))
    readr::write_csv(merged_data, out_file)

    tibble::tibble(
      gene = gene_symbol,
      dataset = dataset_id,
      status = "merged",
      merged_rows = nrow(merged_data),
      excluded_rows = nrow(expr_data) - nrow(merged_data)
    )
  }, error = function(e) {
    tibble::tibble(
      gene = gene_symbol,
      dataset = dataset_id,
      status = paste0("failed: ", conditionMessage(e)),
      merged_rows = NA_integer_,
      excluded_rows = NA_integer_
    )
  })
}

message("\nMerging expression values with processed phenotype data...")

merge_log_all <- purrr::map_dfr(gene_table$symbol, function(gene_symbol) {
  message("  Gene: ", gene_symbol)

  gene_log <- purrr::map_dfr(
    pdata_files,
    ~ merge_expression_with_pdata(gene_symbol, .x)
  )

  log_dir <- file.path(dirs$logs, gene_symbol)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  readr::write_csv(
    gene_log,
    file.path(log_dir, "data_merge_log.csv")
  )

  gene_log
})

readr::write_csv(
  merge_log_all,
  file.path(dirs$logs, "data_merge_log_all_genes.csv")
)


# ================================================================
# 8. Per-dataset model fitting
# ================================================================

fit_dataset_model <- function(file_path, gene_symbol, category) {
  dataset_id <- get_dataset_id(file_path)

  df <- tryCatch(
    readr::read_csv(file_path, show_col_types = FALSE),
    error = function(e) NULL
  )

  if (is.null(df)) {
    return(list(result = NULL, plot = NULL, message = "read_failed"))
  }

  if (!all(c(variable_of_interest, "target_gene_TMM") %in% colnames(df))) {
    return(list(result = NULL, plot = NULL, message = "missing_required_columns"))
  }

  df <- df |>
    dplyr::transmute(
      variable_value = suppressWarnings(as.numeric(.data[[variable_of_interest]])),
      target_gene_TMM = suppressWarnings(as.numeric(target_gene_TMM))
    ) |>
    dplyr::filter(!is.na(variable_value), !is.na(target_gene_TMM))

  if (nrow(df) < min_samples) {
    return(list(result = NULL, plot = NULL, message = "too_few_samples"))
  }

  if (stats::var(df$target_gene_TMM, na.rm = TRUE) == 0) {
    return(list(result = NULL, plot = NULL, message = "zero_expression_variance"))
  }

  lm_fit <- tryCatch(
    stats::lm(target_gene_TMM ~ variable_value, data = df),
    error = function(e) NULL
  )

  if (is.null(lm_fit)) {
    return(list(result = NULL, plot = NULL, message = "lm_failed"))
  }

  lm_summary <- summary(lm_fit)

  if (!"variable_value" %in% rownames(coef(lm_summary))) {
    return(list(result = NULL, plot = NULL, message = "lm_variable_coef_missing"))
  }

  lm_beta <- as.numeric(coef(lm_summary)["variable_value", "Estimate"])
  lm_se   <- as.numeric(coef(lm_summary)["variable_value", "Std. Error"])
  lm_p    <- as.numeric(coef(lm_summary)["variable_value", "Pr(>|t|)"])

  chosen_model <- "LM"
  effect_beta  <- lm_beta
  effect_se    <- lm_se
  effect_p     <- lm_p
  edf_val      <- NA_real_
  gam_fit      <- NULL

  n_unique_variable <- length(unique(df$variable_value))

  if (n_unique_variable >= 3) {
    k_proposed <- ceiling(nrow(df) / 4)
    k_use <- min(k_proposed, n_unique_variable - 1)
    k_use <- max(k_use, 3)
    k_use <- min(k_use, 15)

    gam_fit <- tryCatch(
      mgcv::gam(
        target_gene_TMM ~ s(variable_value, k = k_use, bs = "tp"),
        data = df,
        method = "REML"
      ),
      error = function(e) NULL
    )

    if (!is.null(gam_fit)) {
      gam_summary <- summary(gam_fit)

      edf_val <- if (!is.null(gam_summary$s.table) && nrow(gam_summary$s.table) >= 1) {
        as.numeric(gam_summary$s.table[1, "edf"])
      } else {
        NA_real_
      }

      aic_lm <- tryCatch(stats::AIC(lm_fit), error = function(e) NA_real_)
      aic_gam <- tryCatch(stats::AIC(gam_fit), error = function(e) NA_real_)

      use_gam <- !is.na(edf_val) &&
        edf_val > edf_threshold &&
        !is.na(aic_lm) &&
        !is.na(aic_gam) &&
        aic_gam <= (aic_lm - aic_delta)

      if (use_gam) {
        variable_grid <- seq(min(df$variable_value), max(df$variable_value), length.out = 200)

        preds <- tryCatch(
          predict(gam_fit, newdata = data.frame(variable_value = variable_grid), se.fit = TRUE),
          error = function(e) NULL
        )

        if (!is.null(preds)) {
          derivs <- diff(preds$fit) / diff(variable_grid)

          variable_density <- stats::density(
            df$variable_value,
            from = min(variable_grid),
            to = max(variable_grid),
            n = length(variable_grid) - 1
          )$y

          weights <- variable_density / sum(variable_density, na.rm = TRUE)

          effect_beta <- sum(weights * derivs)

          boot_effects <- future.apply::future_sapply(
            seq_len(n_boot),
            function(b) {
              boot_df <- df[sample(seq_len(nrow(df)), replace = TRUE), ]

              boot_fit <- tryCatch(
                mgcv::gam(
                  target_gene_TMM ~ s(variable_value, k = k_use, bs = "tp"),
                  data = boot_df,
                  method = "REML"
                ),
                error = function(e) NULL
              )

              if (is.null(boot_fit)) {
                return(NA_real_)
              }

              boot_preds <- tryCatch(
                predict(boot_fit, newdata = data.frame(variable_value = variable_grid)),
                error = function(e) NULL
              )

              if (is.null(boot_preds)) {
                return(NA_real_)
              }

              boot_derivs <- diff(boot_preds) / diff(variable_grid)

              boot_variable_density <- stats::density(
                boot_df$variable_value,
                from = min(variable_grid),
                to = max(variable_grid),
                n = length(variable_grid) - 1
              )$y

              boot_weights <- boot_variable_density / sum(boot_variable_density, na.rm = TRUE)

              sum(boot_weights * boot_derivs)
            },
            future.seed = TRUE
          )

          boot_effects <- boot_effects[!is.na(boot_effects)]

          effect_se <- if (length(boot_effects) > 1) {
            stats::sd(boot_effects)
          } else {
            NA_real_
          }

          effect_p <- if (!is.na(effect_se) && effect_se > 0) {
            2 * stats::pt(-abs(effect_beta / effect_se), df = nrow(df) - 2)
          } else {
            NA_real_
          }

          chosen_model <- "GAM"
        }
      }
    }
  }

  # -----------------------------
  # Scatterplot
  # -----------------------------
  df_plot <- df |>
    dplyr::arrange(variable_value)

  variable_min <- min(df_plot$variable_value, na.rm = TRUE)
  variable_max <- max(df_plot$variable_value, na.rm = TRUE)
  variable_grid <- seq(variable_min, variable_max, length.out = 200)

  pred_df <- NULL

  if (identical(chosen_model, "GAM") && !is.null(gam_fit)) {
    p <- tryCatch(
      predict(gam_fit, newdata = data.frame(variable_value = variable_grid), se.fit = TRUE),
      error = function(e) NULL
    )

    if (!is.null(p)) {
      pred_df <- tibble::tibble(
        variable_value = variable_grid,
        fit = as.numeric(p$fit),
        fit_lo = as.numeric(p$fit - 1.96 * p$se.fit),
        fit_hi = as.numeric(p$fit + 1.96 * p$se.fit)
      )
    }
  } else {
    p <- tryCatch(
      predict(
        lm_fit,
        newdata = data.frame(variable_value = variable_grid),
        interval = "confidence",
        level = 0.95
      ),
      error = function(e) NULL
    )

    if (!is.null(p)) {
      pred_df <- tibble::tibble(
        variable_value = variable_grid,
        fit = as.numeric(p[, "fit"]),
        fit_lo = as.numeric(p[, "lwr"]),
        fit_hi = as.numeric(p[, "upr"])
      )
    }
  }

  scatter_plot <- ggplot2::ggplot(df_plot, ggplot2::aes(variable_value, target_gene_TMM)) +
    ggplot2::geom_point(size = 2.2, alpha = 0.6, colour = "grey30") +
    ggplot2::labs(
      title = dataset_id,
      subtitle = paste0(
        "n = ", nrow(df_plot),
        " | ", tolower(variable_label), " = ", variable_min, "–", variable_max,
        " | β = ", sprintf("%.3f", effect_beta),
        " | ", format_p_inline(effect_p)
      ),
      x = paste0(variable_label, " (", variable_units, ")"),
      y = paste0(gene_symbol, " expression (log2 TMM-CPM)")
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(pred_df)) {
    scatter_plot <- scatter_plot +
      ggplot2::geom_ribbon(
        data = pred_df,
        ggplot2::aes(x = variable_value, ymin = fit_lo, ymax = fit_hi),
        fill = "grey70",
        alpha = 0.25,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_line(
        data = pred_df,
        ggplot2::aes(x = variable_value, y = fit),
        colour = "#D1495B",
        linewidth = 1,
        inherit.aes = FALSE
      )
  }

  result <- tibble::tibble(
    dataset = dataset_id,
    gene = gene_symbol,
    category = category,
    chosen_model = chosen_model,
    edf = edf_val,
    beta = as.numeric(effect_beta),
    se = as.numeric(effect_se),
    pval = as.numeric(effect_p),
    n = nrow(df),
    var_min = min(df$variable_value, na.rm = TRUE),
    var_max = max(df$variable_value, na.rm = TRUE)
  )

  list(result = result, plot = scatter_plot, message = "success")
}


# ================================================================
# 9. Meta-analysis loop
# ================================================================

all_meta_diagnostics <- tibble::tibble(
  gene = character(),
  k = integer(),
  Q = numeric(),
  df = numeric(),
  Q_pval = numeric(),
  I_squared = numeric(),
  tau2 = numeric(),
  H_squared = numeric(),
  egger_pval = numeric(),
  begg_pval = numeric()
)

meta_gene_summary <- tibble::tibble(
  Gene = character(),
  Category = character(),
  k = integer(),
  Variable_range = character(),
  `Beta (p-value)` = character(),
  `95% CI` = character(),
  `95% PI` = character(),
  beta_num = double(),
  p_value_num = double(),
  percent_change_scaled = double(),
  Direction = character(),
  Heterogeneity = character(),
  Publication_bias = character(),
  Influential_datasets = character()
)

loo_all_genes <- tibble::tibble(
  Gene = character(),
  Category = character(),
  `Excluded dataset` = character(),
  Full_beta_num = double(),
  LOO_beta_num = double(),
  Full_p_num = double(),
  LOO_p_num = double(),
  `Full model Beta, p-value` = character(),
  `LOO model Beta, p-value` = character(),
  `Significance change?` = logical(),
  `Directional change?` = logical(),
  Interpretation = character(),
  influential = logical()
)

gene_table <- gene_table |>
  dplyr::mutate(
    available_datasets = vapply(
      symbol,
      function(gene) {
        input_dir <- file.path(dirs$pre_proc_data, gene)
        if (!dir.exists(input_dir)) {
          return(0L)
        }
        length(list.files(input_dir, pattern = "_pre_proc_data\\.csv$", full.names = TRUE))
      },
      integer(1)
    )
  ) |>
  dplyr::filter(available_datasets >= min_k)

message("\nGenes retained after dataset-availability filter: ", nrow(gene_table))

if (nrow(gene_table) == 0) {
  stop("No genes had enough available datasets for analysis. Check merge logs.")
}

for (i in seq_len(nrow(gene_table))) {
  gene <- gene_table$symbol[i]
  category <- gene_table$category[i]

  message("\n============================================================")
  message("Processing gene: ", gene, " (", category, ")")
  message("============================================================")

  input_dir <- file.path(dirs$pre_proc_data, gene)
  files <- list.files(input_dir, pattern = "_pre_proc_data\\.csv$", full.names = TRUE)

  model_outputs <- purrr::map(files, fit_dataset_model, gene_symbol = gene, category = category)

  model_log <- tibble::tibble(
    dataset = vapply(files, get_dataset_id, character(1)),
    status = vapply(model_outputs, function(x) x$message, character(1))
  )

  gene_log_dir <- file.path(dirs$logs, gene)
  dir.create(gene_log_dir, recursive = TRUE, showWarnings = FALSE)

  readr::write_csv(model_log, file.path(gene_log_dir, "dataset_model_log.csv"))

  results <- purrr::map(model_outputs, "result") |>
    purrr::compact()

  plot_list <- purrr::map(model_outputs, "plot") |>
    purrr::compact()

  if (length(results) == 0) {
    message("No valid models for gene ", gene)
    next
  }

  effects_df <- dplyr::bind_rows(results)

  gene_n_total <- sum(effects_df$n, na.rm = TRUE)
  gene_n_min <- min(effects_df$n, na.rm = TRUE)
  gene_n_max <- max(effects_df$n, na.rm = TRUE)

  effects_df <- effects_df |>
    dplyr::mutate(
      gene_n_total = gene_n_total,
      gene_n_min = gene_n_min,
      gene_n_max = gene_n_max,
      n_range = paste0(gene_n_min, "–", gene_n_max)
    )

  effects_file <- file.path(dirs$gene_effects, paste0(gene, "_model_effects.csv"))
  readr::write_csv(effects_df, effects_file)

  # -----------------------------
  # Scatterplot grid
  # -----------------------------
  if (length(plot_list) > 0) {
    names(plot_list) <- effects_df$dataset[seq_along(plot_list)]

    n_plots <- length(plot_list)
    ncol_grid <- 3
    nrow_grid <- ceiling(n_plots / ncol_grid)

    grid_plot <- patchwork::wrap_plots(plot_list, ncol = ncol_grid) +
      patchwork::plot_annotation(
        title = paste0(gene, " expression vs ", tolower(variable_label), " by dataset"),
        theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold"))
      )

    grid_file <- file.path(dirs$scatterplot_grids, paste0(gene, "_scatterplots.png"))

    tryCatch({
      ggplot2::ggsave(
        grid_file,
        grid_plot,
        width = 12,
        height = max(4, 3.5 * nrow_grid),
        dpi = 300
      )
      message("Scatterplot grid saved: ", grid_file)
    }, error = function(e) {
      warning("Failed to save scatterplot grid for ", gene, ": ", conditionMessage(e))
    })
  }

  # -----------------------------
  # Meta-analysis
  # -----------------------------
  meta_df <- effects_df |>
    dplyr::filter(!is.na(beta), !is.na(se), se > 0)

  meta_res <- NULL
  meta_variable_range <- NA_character_

  if (nrow(meta_df) >= 2) {
    meta_res <- tryCatch(
      metafor::rma(yi = beta, sei = se, data = meta_df, method = "REML"),
      error = function(e) {
        warning("Meta-analysis failed for ", gene, ": ", conditionMessage(e))
        NULL
      }
    )
  }

  if (!is.null(meta_res)) {
    pi_res <- tryCatch(
      predict(meta_res),
      error = function(e) list(pi.lb = NA_real_, pi.ub = NA_real_)
    )

    PI_lb <- as.numeric(pi_res$pi.lb)
    PI_ub <- as.numeric(pi_res$pi.ub)

    meta_variable_min <- min(meta_df$var_min, na.rm = TRUE)
    meta_variable_max <- max(meta_df$var_max, na.rm = TRUE)
    meta_variable_range <- paste0(meta_variable_min, "–", meta_variable_max, " ", variable_units)

    egger_pval <- NA_real_
    begg_pval <- NA_real_

    if (meta_res$k >= 3 && stats::var(meta_df$se, na.rm = TRUE) != 0) {
      egger_res <- tryCatch(
        metafor::regtest(meta_res, model = "rma", predictor = "sei"),
        error = function(e) NULL
      )

      begg_res <- tryCatch(
        metafor::ranktest(meta_res),
        error = function(e) NULL
      )

      if (!is.null(egger_res)) egger_pval <- as.numeric(egger_res$pval)
      if (!is.null(begg_res)) begg_pval <- as.numeric(begg_res$pval)
    }

    heterogeneity <- paste0(
      "I² = ", round(as.numeric(meta_res$I2), 1), "%, ",
      "τ² = ", signif(as.numeric(meta_res$tau2), 3), ", ",
      "H² = ", signif(as.numeric(meta_res$H2), 3)
    )

    pub_bias <- paste0(
      "Egger p = ", signif(egger_pval, 3),
      "; Begg p = ", signif(begg_pval, 3)
    )

    all_meta_diagnostics <- dplyr::bind_rows(
      all_meta_diagnostics,
      tibble::tibble(
        gene = gene,
        k = meta_res$k,
        Q = as.numeric(meta_res$QE),
        df = meta_res$k - 1,
        Q_pval = as.numeric(meta_res$QEp),
        I_squared = as.numeric(meta_res$I2),
        tau2 = as.numeric(meta_res$tau2),
        H_squared = as.numeric(meta_res$H2),
        egger_pval = egger_pval,
        begg_pval = begg_pval
      )
    )

    meta_raw_tbl <- tibble::tibble(
      Gene = gene,
      Category = category,
      k = meta_res$k,
      beta = as.numeric(meta_res$b[[1]]),
      se = as.numeric(meta_res$se[[1]]),
      zval = as.numeric(meta_res$zval[[1]]),
      pval = as.numeric(meta_res$pval[[1]]),
      ci_lb = as.numeric(meta_res$ci.lb[[1]]),
      ci_ub = as.numeric(meta_res$ci.ub[[1]]),
      QE = as.numeric(meta_res$QE),
      QEp = as.numeric(meta_res$QEp),
      I2 = as.numeric(meta_res$I2),
      tau2 = as.numeric(meta_res$tau2),
      H2 = as.numeric(meta_res$H2),
      PI_lb = PI_lb,
      PI_ub = PI_ub,
      Egger_p = egger_pval,
      Begg_p = begg_pval
    )

    meta_file <- file.path(dirs$meta_effects, paste0(gene, "_meta_effects.csv"))
    readr::write_csv(meta_raw_tbl, meta_file)

    # -----------------------------
    # Funnel plot
    # -----------------------------
    funnel_file <- file.path(dirs$funnel_plots, paste0(gene, "_funnel.png"))

    tryCatch({
      grDevices::png(funnel_file, width = 1800, height = 1200, res = 300)
      metafor::funnel(
        meta_res,
        xlab = paste0("β coefficient (log2 TMM-CPM change per 1 ", variable_units, ")"),
        ylab = "Standard error",
        main = paste0("Funnel plot for ", tolower(variable_label), " effects on ", gene, " expression")
      )
      grDevices::dev.off()
    }, error = function(e) {
      if (grDevices::dev.cur() > 1) grDevices::dev.off()
      warning("Failed to generate funnel plot for ", gene, ": ", conditionMessage(e))
    })

    # -----------------------------
    # Baujat plot
    # -----------------------------
    baujat_file <- file.path(dirs$baujat_plots, paste0(gene, "_baujat.png"))

    tryCatch({
      baujat_data <- metafor::baujat(meta_res, plot = FALSE)

      grDevices::png(baujat_file, width = 1800, height = 1200, res = 300)
      graphics::plot(
        baujat_data$x,
        baujat_data$y,
        pch = 16,
        col = "black",
        xlab = "Contribution to heterogeneity (Q)",
        ylab = "Influence on pooled effect",
        main = paste0("Baujat plot for ", tolower(variable_label), " effects on ", gene, " expression")
      )
      grDevices::dev.off()
    }, error = function(e) {
      if (grDevices::dev.cur() > 1) grDevices::dev.off()
      warning("Failed to generate Baujat plot for ", gene, ": ", conditionMessage(e))
    })

    # -----------------------------
    # Leave-one-out analysis
    # -----------------------------
    influential_datasets <- "None"

    if (meta_res$k >= 3) {
      loo_res <- tryCatch(metafor::leave1out(meta_res), error = function(e) NULL)

      if (!is.null(loo_res)) {
        loo_df <- as.data.frame(loo_res)

        full_est <- as.numeric(meta_res$b)
        full_p <- as.numeric(meta_res$pval)
        full_sig <- full_p < alpha
        full_dir <- sign(full_est)

        loo_df$dataset_ID <- as.character(meta_df$dataset[seq_len(nrow(loo_df))])

        loo_long_tbl <- loo_df |>
          dplyr::mutate(
            loo_sig = pval < alpha,
            loo_dir = sign(estimate),
            significance_change = loo_sig != full_sig,
            direction_change = loo_dir != full_dir & (loo_sig | full_sig),
            influential = significance_change | direction_change,
            Interpretation = dplyr::case_when(
              significance_change & direction_change ~ "Significance and direction changed",
              significance_change ~ "Significance changed",
              direction_change ~ "Direction changed",
              TRUE ~ "No meaningful change"
            )
          ) |>
          dplyr::transmute(
            Gene = as.character(gene),
            Category = as.character(category),
            `Excluded dataset` = as.character(dataset_ID),
            Full_beta_num = as.numeric(full_est),
            LOO_beta_num = as.numeric(estimate),
            Full_p_num = as.numeric(full_p),
            LOO_p_num = as.numeric(pval),
            `Full model Beta, p-value` = paste0(
              format_beta(full_est),
              " (",
              format_p(full_p),
              ")"
            ),
            `LOO model Beta, p-value` = paste0(
              format_beta(estimate),
              " (",
              format_p(pval),
              ")"
            ),
            `Significance change?` = as.logical(significance_change),
            `Directional change?` = as.logical(direction_change),
            Interpretation = as.character(Interpretation),
            influential = as.logical(influential)
          )

        loo_all_genes <- dplyr::bind_rows(loo_all_genes, loo_long_tbl)

        influential_ids <- loo_long_tbl |>
          dplyr::filter(influential) |>
          dplyr::pull(`Excluded dataset`)

        influential_datasets <- if (length(influential_ids) == 0) {
          "None"
        } else {
          paste(influential_ids, collapse = ", ")
        }

        loo_file <- file.path(dirs$loo_tables, paste0(gene, "_loo.csv"))
        readr::write_csv(loo_long_tbl, loo_file)
      }
    }

    gene_summary_row <- tibble::tibble(
      Gene = as.character(gene),
      Category = as.character(category),
      k = as.integer(meta_res$k),
      Variable_range = meta_variable_range,
      `Beta (p-value)` = paste0(
        format_beta(meta_res$b[[1]]),
        " (",
        format_p(meta_res$pval[[1]]),
        ")"
      ),
      `95% CI` = paste0(
        round(as.numeric(meta_res$ci.lb[[1]]), 4),
        ", ",
        round(as.numeric(meta_res$ci.ub[[1]]), 4)
      ),
      `95% PI` = paste0(
        round(as.numeric(PI_lb), 4),
        ", ",
        round(as.numeric(PI_ub), 4)
      ),
      beta_num = as.numeric(meta_res$b[[1]]),
      p_value_num = as.numeric(meta_res$pval[[1]]),
      percent_change_scaled = percent_change_scaled(meta_res$b[[1]]),
      Direction = ifelse(as.numeric(meta_res$b[[1]]) > 0, "↑", "↓"),
      Heterogeneity = as.character(heterogeneity),
      Publication_bias = as.character(pub_bias),
      Influential_datasets = as.character(influential_datasets)
    )

    meta_gene_summary <- dplyr::bind_rows(meta_gene_summary, gene_summary_row)
  }

  # -----------------------------
  # Forest plot
  # -----------------------------
  if (nrow(effects_df) >= min_k) {
    beta_min <- min(effects_df$beta - 1.96 * effects_df$se, na.rm = TRUE)
    beta_max <- max(effects_df$beta + 1.96 * effects_df$se, na.rm = TRUE)
    beta_range <- beta_max - beta_min
    padding <- ifelse(beta_range == 0, 0.01, beta_range * 0.05)

    effects_plot_df <- effects_df |>
      dplyr::mutate(
        Percent_Change_scaled = percent_change_scaled(beta),
        CI_Lower = beta - 1.96 * se,
        CI_Upper = beta + 1.96 * se,
        pval_text = dplyr::if_else(
          pval < 0.001,
          "< 0.001",
          paste0("= ", sprintf("%.3f", pval))
        ),
        Significance = dplyr::if_else(pval < alpha, "Significant", "Not Significant"),
        Delta_expr_label = paste0(
          ifelse(Percent_Change_scaled > 0, "+", ""),
          Percent_Change_scaled,
          "% ", effect_scale_label
        ),
        Label = paste0(
          dataset,
          "\n n = ", n,
          ", ", tolower(variable_label), " range = ", var_min, "–", var_max,
          "\nβ = ", sprintf("%.3f", beta),
          ", p ", pval_text,
          "\nΔ expression = ", Delta_expr_label
        ),
        Type = "Dataset"
      )

    plot_data <- effects_plot_df

    if (!is.null(meta_res)) {
      meta_row <- tibble::tibble(
        dataset = "Pooled Meta Effect",
        beta = as.numeric(meta_res$b),
        se = as.numeric(meta_res$se),
        pval = as.numeric(meta_res$pval),
        n = meta_res$k,
        Percent_Change_scaled = percent_change_scaled(meta_res$b),
        CI_Lower = as.numeric(meta_res$ci.lb),
        CI_Upper = as.numeric(meta_res$ci.ub),
        pval_text = ifelse(
          meta_res$pval < 0.001,
          "< 0.001",
          paste0("= ", sprintf("%.3f", meta_res$pval))
        ),
        Significance = ifelse(meta_res$pval < alpha, "Significant", "Not Significant"),
        Delta_expr_label = paste0(
          ifelse(Percent_Change_scaled > 0, "+", ""),
          Percent_Change_scaled,
          "% ", effect_scale_label
        ),
        Label = paste0(
          "**Pooled Meta Effect**",
          "\n", variable_label, " range = ", meta_variable_range,
          "\nβ = ", sprintf("%.3f", beta),
          ", p ", pval_text,
          "\nΔ expression = ", Delta_expr_label
        ),
        Type = "Meta"
      )

      plot_data <- dplyr::bind_rows(effects_plot_df, meta_row)
    }

    dataset_levels <- effects_plot_df |>
      dplyr::arrange(beta) |>
      dplyr::pull(Label)

    all_levels <- if ("Meta" %in% plot_data$Type) {
      c(dataset_levels, plot_data$Label[plot_data$Type == "Meta"])
    } else {
      dataset_levels
    }

    plot_data$Label <- factor(plot_data$Label, levels = rev(all_levels))

    forest_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = beta, y = Label)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = CI_Lower, xmax = CI_Upper),
        height = 0.2,
        linewidth = 0.8
      ) +
      ggplot2::geom_point(ggplot2::aes(color = Significance), size = 3) +
      ggplot2::scale_color_manual(
        values = c("Significant" = "#D1495B", "Not Significant" = "#2E86AB")
      ) +
      ggplot2::scale_x_continuous(
        name = paste0("Effect size β (", gene, " log2 TMM-CPM change per 1 ", variable_units, ")"),
        limits = c(beta_min - padding, beta_max + padding)
      ) +
      ggplot2::ylab("") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        axis.text.y = ggtext::element_markdown(),
        axis.title.y = ggplot2::element_text(size = 12),
        legend.position = "none",
        panel.grid.minor = ggplot2::element_blank()
      ) +
      ggplot2::coord_cartesian(clip = "off")

    forest_file <- file.path(dirs$forest_plots, paste0(gene, "_forestplot.png"))

    ggplot2::ggsave(
      forest_file,
      forest_plot,
      width = 8,
      height = max(4, 0.8 * nrow(plot_data)),
      dpi = 300
    )
  }
}


# ================================================================
# 10. Save publication tables
# ================================================================

readr::write_csv(
  all_meta_diagnostics,
  file.path(dirs$tables, "meta_diagnostics_all_genes.csv")
)

readr::write_csv(
  meta_gene_summary,
  file.path(dirs$meta_effects, "meta_gene_summary_with_numeric_columns.csv")
)

if (nrow(loo_all_genes) > 0) {
  loo_full_table <- loo_all_genes

  loo_influential_table <- loo_all_genes |>
    dplyr::filter(influential)

  pub_ready_loo_full_tbl <- loo_full_table |>
    dplyr::select(-Full_beta_num, -Full_p_num, -LOO_beta_num, -LOO_p_num)

  pub_ready_loo_influential_tbl <- loo_influential_table |>
    dplyr::select(-Full_beta_num, -Full_p_num, -LOO_beta_num, -LOO_p_num)

  readr::write_csv(
    pub_ready_loo_full_tbl,
    file.path(dirs$publication_tables, "loo_full_all_genes.csv")
  )

  readr::write_csv(
    pub_ready_loo_influential_tbl,
    file.path(dirs$publication_tables, "loo_influential_only.csv")
  )
}

if (nrow(meta_gene_summary) > 0) {
  meta_category_summary <- meta_gene_summary |>
    dplyr::group_by(Category) |>
    dplyr::summarise(
      `Genes analysed` = dplyr::n(),
      `Genes associated with variable` = sum(p_value_num < alpha, na.rm = TRUE),
      `Direction of association` = {
        sig_betas <- beta_num[p_value_num < alpha]

        n_pos <- sum(sig_betas > 0, na.rm = TRUE)
        n_neg <- sum(sig_betas < 0, na.rm = TRUE)

        if (length(sig_betas) == 0) {
          "–"
        } else if (n_pos > 0 && n_neg > 0) {
          paste0("↑ (", n_pos, "), ↓ (", n_neg, ")")
        } else if (n_pos > 0) {
          paste0("↑ (", n_pos, ")")
        } else {
          paste0("↓ (", n_neg, ")")
        }
      },
      `Beta range` = paste0(
        round(min(beta_num, na.rm = TRUE), 3),
        " – ",
        round(max(beta_num, na.rm = TRUE), 3)
      ),
      `Sensitivity issues?` = ifelse(
        any(Influential_datasets != "None"),
        "Yes",
        "No"
      ),
      .groups = "drop"
    )

  readr::write_csv(
    meta_category_summary,
    file.path(dirs$publication_tables, "meta_category_summary.csv")
  )

  pub_ready_gene_summary <- meta_gene_summary |>
    dplyr::select(-beta_num, -p_value_num)

  readr::write_csv(
    pub_ready_gene_summary,
    file.path(dirs$publication_tables, "meta_gene_summary.csv")
  )
}


# ================================================================
# 11. Jungle plots
# ================================================================

if (nrow(meta_gene_summary) > 0) {
  
  all_categories <- unique(meta_gene_summary$Category)
  
  # If user leaves groupings as NULL, include all categories
  jungle_groupings <- purrr::map(
    jungle_groupings,
    ~ if (is.null(.x)) all_categories else .x
  )
  
  jungle_facet_groupings <- purrr::map(
    jungle_facet_groupings,
    ~ if (is.null(.x)) all_categories else .x
  )
  
  meta_plot_df <- meta_gene_summary |>
    tidyr::separate(`95% CI`, into = c("ci_lb", "ci_ub"), sep = ",\\s*", convert = TRUE) |>
    dplyr::mutate(
      is_sig = p_value_num < alpha,
      forest_label = paste0(
        Gene,
        "\n",
        "k = ", k,
        ", β = ", sprintf("%.3f", beta_num),
        ", ",
        format_p_inline(p_value_num),
        "\nΔ expression = ",
        sprintf("%+.1f", percent_change_scaled),
        "% ", effect_scale_label
      ),
      category = factor(Category),
      group_facet = purrr::map_chr(
        as.character(category),
        function(cat) {
          found <- names(jungle_facet_groupings)[
            vapply(jungle_facet_groupings, function(x) cat %in% x, logical(1))
          ]
          
          if (length(found) == 0) "Other" else found[1]
        }
      )
    ) |>
    dplyr::mutate(
      group_facet = factor(group_facet, levels = c(names(jungle_facet_groupings), "Other")),
      group_facet_wrapped = stringr::str_wrap(group_facet, width = 20)
    )
  
  purrr::walk2(names(jungle_groupings), jungle_groupings, function(plot_name, cats) {
    
    df_group <- meta_plot_df |>
      dplyr::filter(as.character(category) %in% cats)
    
    if (nrow(df_group) == 0) {
      return(NULL)
    }
    
    df_group <- df_group |>
      dplyr::arrange(beta_num) |>
      dplyr::mutate(forest_label = forcats::fct_inorder(forest_label))
    
    p <- ggplot2::ggplot(df_group, ggplot2::aes(x = beta_num, y = forest_label, color = is_sig)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = ci_lb, xmax = ci_ub),
        height = 0.3,
        colour = "grey40"
      ) +
      ggplot2::geom_point(size = 3) +
      ggplot2::facet_grid(
        group_facet_wrapped ~ .,
        scales = "free_y",
        space = "free_y"
      ) +
      ggplot2::scale_color_manual(
        values = c("TRUE" = "#D1495B", "FALSE" = "#2E86AB")
      ) +
      ggplot2::labs(
        title = NULL,
        x = paste0("Pooled effect size β (log2 TMM-CPM change per 1 ", variable_units, ")"),
        y = ""
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 7),
        strip.text.y = ggplot2::element_text(face = "bold", size = 7),
        strip.background = ggplot2::element_rect(fill = "grey90", colour = NA),
        legend.position = "none",
        panel.spacing = grid::unit(1.2, "lines"),
        panel.grid.minor = ggplot2::element_blank()
      )
    
    save_path <- file.path(
      dirs$jungle_plots,
      paste0(
        "jungle_plot_",
        stringr::str_replace_all(tolower(plot_name), "[^a-z0-9]+", "_"),
        ".png"
      )
    )
    
    ggplot2::ggsave(
      save_path,
      p,
      width = 10,
      height = max(6, 0.5 * nrow(df_group)),
      dpi = 600
    )
  })
}

# ================================================================
# 12. Dataset characteristics table
# ================================================================

dataset_characteristics <- purrr::map_dfr(pdata_files, function(pdata_file) {
  dataset_id <- get_dataset_id(pdata_file)

  pdata <- tryCatch(
    readr::read_csv(pdata_file, show_col_types = FALSE),
    error = function(e) NULL
  )

  if (is.null(pdata)) {
    return(tibble::tibble(
      Dataset_ID = dataset_id,
      Sample_size = NA_integer_,
      Variable = "Failed to read pData",
      BMI = "Not reported",
      Tissue = "Not reported",
      Instrument_model = "Not reported",
      Platform_id = "Not reported",
      Reference_genome = "Not reported",
      Library_prep = "",
      Read_type_length = "",
      Sequencing_depth = ""
    ))
  }

  variable_col <- first_matching_column(pdata, c(paste0("^", tolower(variable_of_interest), "$")))
  bmi_col <- first_matching_column(pdata, c("^bmi$"))
  tissue_col <- first_matching_column(pdata, c("tissue", "source", "organ"))
  instrument_col <- first_matching_column(pdata, c("instrument", "sequencer", "model"))
  platform_col <- first_matching_column(pdata, c("platform", "gpl"))
  genome_col <- first_matching_column(pdata, c("genome", "reference", "assembly"))

  tibble::tibble(
    Dataset_ID = dataset_id,
    Sample_size = nrow(pdata),
    Variable = if (!is.na(variable_col)) safe_mean_sd_range(pdata[[variable_col]]) else "Not reported",
    BMI = if (!is.na(bmi_col)) safe_mean_sd_range(pdata[[bmi_col]]) else "Not reported",
    Tissue = if (!is.na(tissue_col)) detect_most_common_value(pdata[[tissue_col]]) else "Not reported",
    Instrument_model = if (!is.na(instrument_col)) detect_most_common_value(pdata[[instrument_col]]) else "Not reported",
    Platform_id = if (!is.na(platform_col)) detect_most_common_value(pdata[[platform_col]]) else "Not reported",
    Reference_genome = if (!is.na(genome_col)) detect_most_common_value(pdata[[genome_col]]) else "Not reported",
    Library_prep = "",
    Read_type_length = "",
    Sequencing_depth = ""
  )
})

readr::write_csv(
  dataset_characteristics,
  file.path(dirs$tables, "dataset_characteristics.csv")
)


# ================================================================
# 13. Finish
# ================================================================

future::plan(future::sequential)

message(
  "\n============================================================\n",
  "META-ANALYSIS PIPELINE COMPLETE\n",
  "============================================================\n\n",
  "Key outputs:\n",
  "  outputs/tables/gene_effects/\n",
  "  outputs/tables/meta_effects/\n",
  "  outputs/tables/publication_tables/\n",
  "  outputs/figures/scatterplot_grids/\n",
  "  outputs/figures/forest_plots/\n",
  "  outputs/figures/funnel_plots/\n",
  "  outputs/figures/baujat_plots/\n",
  "  outputs/figures/jungle_plots/\n",
  "  outputs/logs/\n",
  "============================================================\n"
)
