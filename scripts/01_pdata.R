# ================================================================
# 01_download_and_prepare_pdata.R
# RNA-seq GEO phenotype metadata download script
#
# Purpose:
#   1. Create the required project directory structure.
#   2. Create a template inputs/potential_seq_datalist.csv if absent.
#   3. Import screened GEO Series IDs.
#   4. Download raw GEO phenotype metadata (pData).
#   5. Save a pData download log.
#
# Manual step after this script:
#   Review and manually format phenotype metadata files before running:
#   02_download_counts_and_tmm_normalise.R
# ================================================================


# ================================================================
# 0. Package setup
# ================================================================

cran_packages <- c(
  "dplyr",
  "readr",
  "tibble",
  "stringr",
  "here"
)

bioconductor_packages <- c(
  "GEOquery"
)

install_cran_packages <- function(packages) {
  missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing_packages) > 0) {
    install.packages(missing_packages, dependencies = TRUE)
  }
}

install_bioconductor_packages <- function(packages) {
  missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing_packages) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }

    BiocManager::install(missing_packages, ask = FALSE, update = FALSE)
  }
}

install_cran_packages(cran_packages)
install_bioconductor_packages(bioconductor_packages)

invisible(lapply(c(cran_packages, bioconductor_packages), library, character.only = TRUE))


# ================================================================
# 1. Project directories
# ================================================================

project_root <- here::here()

dirs <- list(
  inputs     = file.path(project_root, "inputs"),
  outputs    = file.path(project_root, "outputs"),
  logs       = file.path(project_root, "outputs", "logs"),
  pdata_raw  = file.path(project_root, "inputs", "pdata", "raw_pdata"),
  pdata_proc = file.path(project_root, "inputs", "pdata", "proc_pdata")
)

invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

message("Project directories created or already exist.")


# ================================================================
# 2. Helper functions
# ================================================================

clean_gse_id <- function(x) {
  x |>
    stringr::str_trim() |>
    stringr::str_replace_all("E0", "E")
}

download_pdata <- function(gse_id, output_dir) {
  tryCatch({
    geo_obj <- GEOquery::getGEO(gse_id)

    if (is.null(geo_obj) || length(geo_obj) == 0 || is.null(geo_obj[[1]])) {
      stop("Invalid or empty GEO object.")
    }

    # pData() is made available when GEOquery and its dependencies are loaded.
    # This avoids explicit Biobase::pData() namespacing for systems where
    # Biobase is not directly available under the expected package name.
    if (!exists("pData", mode = "function")) {
      stop("pData() is not available. Restart R and load GEOquery before retrying.")
    }

    sample_info <- pData(geo_obj[[1]]) |>
      as.data.frame()

    output_file <- file.path(output_dir, paste0(gse_id, "_raw_seq_pdata.csv"))

    utils::write.csv(
      sample_info,
      output_file,
      row.names = TRUE,
      quote = TRUE,
      fileEncoding = "UTF-8"
    )

    tibble::tibble(
      GSE_ID = gse_id,
      success = TRUE,
      output_file = output_file,
      error_message = NA_character_
    )
  }, error = function(e) {
    tibble::tibble(
      GSE_ID = gse_id,
      success = FALSE,
      output_file = NA_character_,
      error_message = conditionMessage(e)
    )
  })
}

download_pdata_with_retries <- function(gse_ids, output_dir, max_retries = 2) {
  results <- list()
  pending <- gse_ids

  for (attempt in seq_len(max_retries)) {
    if (length(pending) == 0) break

    message("pData download pass ", attempt, " of ", max_retries)

    attempt_results <- lapply(pending, download_pdata, output_dir = output_dir) |>
      dplyr::bind_rows() |>
      dplyr::mutate(attempt = attempt)

    results[[attempt]] <- attempt_results

    pending <- attempt_results |>
      dplyr::filter(!success) |>
      dplyr::pull(GSE_ID)

    if (length(pending) > 0 && attempt < max_retries) {
      Sys.sleep(2^attempt)
    }
  }

  dplyr::bind_rows(results) |>
    dplyr::group_by(GSE_ID) |>
    dplyr::slice_tail(n = 1) |>
    dplyr::ungroup()
}


# ================================================================
# 3. Create or import screened GEO dataset list
# ================================================================

dataset_list_path <- file.path(dirs$inputs, "potential_seq_datalist.csv")

if (!file.exists(dataset_list_path)) {
  readr::write_csv(
    tibble::tibble(GSE_ID = character()),
    dataset_list_path
  )

  stop(
    "Created template file: inputs/potential_seq_datalist.csv\n\n",
    "MANUAL STEP REQUIRED:\n",
    "  1. Open inputs/potential_seq_datalist.csv.\n",
    "  2. Add screened NCBI GEO Series IDs under the GSE_ID column.\n",
    "  3. Save the file.\n",
    "  4. Re-run this script."
  )
}

potential_seq_datalist <- readr::read_csv(
  dataset_list_path,
  show_col_types = FALSE
) |>
  dplyr::mutate(GSE_ID = clean_gse_id(GSE_ID)) |>
  dplyr::filter(!is.na(GSE_ID), GSE_ID != "") |>
  dplyr::distinct(GSE_ID)

if (nrow(potential_seq_datalist) == 0) {
  stop(
    "No GSE IDs found in inputs/potential_seq_datalist.csv.\n\n",
    "MANUAL STEP REQUIRED:\n",
    "  1. Add screened NCBI GEO Series IDs under the GSE_ID column.\n",
    "  2. Save the file.\n",
    "  3. Re-run this script."
  )
}

message("Loaded ", nrow(potential_seq_datalist), " GEO dataset IDs.")


# ================================================================
# 4. Download raw GEO phenotype metadata
# ================================================================

pdata_log <- download_pdata_with_retries(
  gse_ids = potential_seq_datalist$GSE_ID,
  output_dir = dirs$pdata_raw,
  max_retries = 2
)

pdata_log_path <- file.path(dirs$logs, "pdata_download_log.csv")

readr::write_csv(
  pdata_log,
  pdata_log_path
)

message("pData download log saved to: ", pdata_log_path)


# ================================================================
# 5. Manual checkpoint before script 02
# ================================================================

message(
  "\n============================================================\n",
  "MANUAL CHECKPOINT REQUIRED BEFORE RUNNING SCRIPT 02\n",
  "============================================================\n\n",
  "1. Review the pData download log:\n",
  "   outputs/logs/pdata_download_log.csv\n\n",
  "2. Review raw phenotype metadata files in:\n",
  "   inputs/pdata/raw_pdata/\n\n",
  "3. For each included dataset, manually format the phenotype file so that\n",
  "   required variables use consistent column names across datasets.\n",
  "   Example:\n",
  "     age\n",
  "     sex\n",
  "     tissue\n",
  "     dataset-specific inclusion/exclusion variables\n\n",
  "4. Save each processed phenotype file to:\n",
  "   inputs/pdata/proc_pdata/\n\n",
  "5. Use the required filename format:\n",
  "   GSEID_proc_seq_pdata.csv\n\n",
  "Example:\n",
  "   GSE111010_proc_seq_pdata.csv\n\n",
  "6. Once all included pData files have been formatted and saved,\n",
  "   run:\n",
  "   02_download_counts_and_tmm_normalise.R\n",
  "============================================================\n"
)
