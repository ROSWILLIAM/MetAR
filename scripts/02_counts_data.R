# ================================================================
# 02_download_counts_and_tmm_normalise.R
# RNA-seq GEO count download and TMM normalisation script
#
# Purpose:
#   1. Detect manually processed pData files.
#   2. Build the final dataset list from processed pData filenames.
#   3. Download GEO supplementary count-like files.
#   4. Convert supported count tables to CSV.
#   5. Perform TMM normalisation and export log2 TMM CPM matrices.
#   6. Save reproducibility logs.
#
# Prerequisite:
#   Run 01_download_and_prepare_pdata.R first and complete the manual
#   phenotype metadata formatting checkpoint.
# ================================================================


# ================================================================
# 0. Package setup
# ================================================================

cran_packages <- c(
  "dplyr",
  "readr",
  "readxl",
  "tibble",
  "stringr",
  "tools",
  "xml2",
  "httr",
  "R.utils",
  "here"
)

bioconductor_packages <- c(
  "edgeR"
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

  pdata_proc = file.path(project_root, "inputs", "pdata", "proc_pdata"),

  counts_raw = file.path(project_root, "inputs", "counts_data", "raw_counts_data"),

  tmm_raw    = file.path(project_root, "inputs", "tmm_data", "raw_tmm_data"),
  tmm_proc   = file.path(project_root, "inputs", "tmm_data", "proc_tmm_data")
)

invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

message("Project directories created or already exist.")


# ================================================================
# 2. Helper functions
# ================================================================

is_raw_counts_file <- function(filename) {
  filename_lower <- tolower(filename)

  grepl("count", filename_lower) &&
    grepl("\\.(txt|csv|tsv|xlsx|xls)$", filename_lower)
}

is_archive_file <- function(filename) {
  grepl("\\.(gz|tar|tar\\.gz|tgz|zip)$", tolower(filename))
}

detect_normalisation <- function(filenames) {
  text <- tolower(paste(filenames, collapse = " "))

  tmm_flag <- grepl("\\btmm\\b", text)
  norm_flag <- grepl("norm|normalised|normalized|deseq|vst|rlog|rpkm|fpkm|tpm|cpm", text)

  list(
    is_normalised = norm_flag || tmm_flag,
    is_tmm = tmm_flag
  )
}

extract_archive <- function(archive_path, extract_dir) {
  dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

  if (grepl("\\.(tar\\.gz|tgz)$", archive_path, ignore.case = TRUE)) {
    utils::untar(archive_path, exdir = extract_dir)
  } else if (grepl("\\.tar$", archive_path, ignore.case = TRUE)) {
    utils::untar(archive_path, exdir = extract_dir)
  } else if (grepl("\\.zip$", archive_path, ignore.case = TRUE)) {
    utils::unzip(archive_path, exdir = extract_dir)
  } else if (grepl("\\.gz$", archive_path, ignore.case = TRUE)) {
    dest_path <- file.path(
      extract_dir,
      sub("\\.gz$", "", basename(archive_path), ignore.case = TRUE)
    )

    R.utils::gunzip(
      archive_path,
      destname = dest_path,
      overwrite = TRUE,
      remove = FALSE
    )
  }
}

read_tabular_file <- function(file_path) {
  ext <- tolower(tools::file_ext(file_path))

  if (ext %in% c("txt", "tsv")) {
    readr::read_delim(
      file_path,
      delim = "\t",
      show_col_types = FALSE,
      name_repair = "minimal"
    )
  } else if (ext == "csv") {
    readr::read_csv(
      file_path,
      show_col_types = FALSE,
      name_repair = "minimal"
    )
  } else if (ext %in% c("xls", "xlsx")) {
    readxl::read_excel(file_path)
  } else {
    stop("Unsupported file type: ", ext)
  }
}

convert_to_csv <- function(file_path, output_dir, gse_id, index = NULL) {
  output_name <- if (is.null(index)) {
    paste0(gse_id, "_raw_counts.csv")
  } else {
    paste0(gse_id, "_raw_counts_", index, ".csv")
  }

  output_path <- file.path(output_dir, output_name)

  tryCatch({
    df <- read_tabular_file(file_path)
    readr::write_csv(df, output_path)
    output_path
  }, error = function(e) {
    message("Failed to convert: ", basename(file_path), " | ", conditionMessage(e))
    NULL
  })
}


# ================================================================
# 3. Build dataset list from manually processed pData files
# ================================================================

proc_pdata_files <- list.files(
  dirs$pdata_proc,
  pattern = "_proc_seq_pdata\\.csv$",
  full.names = FALSE
)

if (length(proc_pdata_files) == 0) {
  stop(
    "No processed pData files found in inputs/pdata/proc_pdata.\n\n",
    "MANUAL STEP REQUIRED:\n",
    "  1. Run 01_download_and_prepare_pdata.R.\n",
    "  2. Review raw pData files in inputs/pdata/raw_pdata/.\n",
    "  3. Format included phenotype files with consistent variable names.\n",
    "  4. Save them to inputs/pdata/proc_pdata/ using:\n",
    "     GSEID_proc_seq_pdata.csv\n",
    "  5. Re-run this script."
  )
}

seq_datalist <- tibble::tibble(
  GSE_ID = stringr::str_remove(proc_pdata_files, "_proc_seq_pdata\\.csv$")
)

message("Found ", nrow(seq_datalist), " processed phenotype files.")


# ================================================================
# 4. Download GEO supplementary count files
# ================================================================

download_gse_counts <- function(gse_id, output_dir, max_size_mb = 200) {
  geo_url <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gse_id)

  page <- tryCatch({
    response <- httr::GET(
      geo_url,
      httr::user_agent("Mozilla/5.0"),
      httr::timeout(60)
    )

    if (httr::status_code(response) != 200) {
      stop("HTTP status ", httr::status_code(response))
    }

    xml2::read_html(httr::content(response, as = "text", encoding = "UTF-8"))
  }, error = function(e) {
    return(NULL)
  })

  if (is.null(page)) {
    return(list(
      success = FALSE,
      error = "Failed to read GEO page.",
      original_filenames = character()
    ))
  }

  supp_links <- page |>
    xml2::xml_find_all("//a[contains(@href, 'geo/download/')]") |>
    xml2::xml_attr("href")

  if (length(supp_links) == 0) {
    return(list(
      success = FALSE,
      error = "No supplementary files found.",
      original_filenames = character()
    ))
  }

  supp_links <- supp_links[
    !grepl("\\.(sra|fastq|fq|bam|cram)$", tolower(supp_links))
  ]

  if (length(supp_links) == 0) {
    return(list(
      success = FALSE,
      error = "Supplementary files contained only raw sequencing files.",
      original_filenames = character()
    ))
  }

  supp_links_full <- paste0("https://www.ncbi.nlm.nih.gov", supp_links)

  kept_files <- character()
  original_filenames <- character()

  for (link in supp_links_full) {
    tryCatch({
      head_response <- httr::HEAD(
        link,
        httr::user_agent("Mozilla/5.0"),
        httr::timeout(60)
      )

      file_size <- as.numeric(httr::headers(head_response)[["content-length"]])

      if (!is.na(file_size) && file_size / 1024^2 > max_size_mb) {
        message(
          "Skipping large file for ", gse_id, ": ",
          round(file_size / 1024^2, 1), " MB"
        )
        next
      }

      response <- httr::GET(
        link,
        httr::user_agent("Mozilla/5.0"),
        httr::timeout(120)
      )

      if (httr::status_code(response) != 200) next

      content_disposition <- httr::headers(response)[["content-disposition"]]

      file_name <- if (is.null(content_disposition)) {
        basename(link)
      } else {
        sub('.*filename="?([^";]+)"?.*', "\\1", content_disposition)
      }

      dest_path <- file.path(output_dir, paste0(gse_id, "_", file_name))
      writeBin(httr::content(response, "raw"), dest_path)

      if (is_archive_file(file_name)) {
        extract_dir <- tempfile(pattern = paste0(gse_id, "_extracted_"))
        dir.create(extract_dir, recursive = TRUE)

        extract_archive(dest_path, extract_dir)

        extracted_files <- list.files(
          extract_dir,
          full.names = TRUE,
          recursive = TRUE
        )

        raw_files <- extracted_files[vapply(
          extracted_files,
          is_raw_counts_file,
          logical(1)
        )]

        if (length(raw_files) > 0) {
          for (i in seq_along(raw_files)) {
            converted <- convert_to_csv(
              raw_files[i],
              output_dir = output_dir,
              gse_id = gse_id,
              index = i
            )

            if (!is.null(converted)) {
              kept_files <- c(kept_files, converted)
              original_filenames <- c(original_filenames, basename(raw_files[i]))
            }
          }
        }

        unlink(extract_dir, recursive = TRUE, force = TRUE)
        file.remove(dest_path)

      } else if (is_raw_counts_file(file_name)) {
        converted <- convert_to_csv(
          dest_path,
          output_dir = output_dir,
          gse_id = gse_id
        )

        if (!is.null(converted)) {
          kept_files <- c(kept_files, converted)
          original_filenames <- c(original_filenames, file_name)
        }
      } else {
        file.remove(dest_path)
      }
    }, error = function(e) {
      message("Failed to process supplementary file for ", gse_id, ": ", conditionMessage(e))
    })
  }

  if (length(kept_files) > 0) {
    list(
      success = TRUE,
      error = NA_character_,
      original_filenames = original_filenames
    )
  } else {
    list(
      success = FALSE,
      error = "No raw count-like files found.",
      original_filenames = original_filenames
    )
  }
}

download_counts_with_retries <- function(gse_id, output_dir, max_attempts = 3) {
  success <- FALSE
  attempts <- 0
  result <- NULL

  while (!success && attempts < max_attempts) {
    attempts <- attempts + 1

    message("Downloading counts for ", gse_id, " | attempt ", attempts)

    result <- download_gse_counts(
      gse_id = gse_id,
      output_dir = output_dir
    )

    success <- isTRUE(result$success)

    if (!success && attempts < max_attempts) {
      Sys.sleep(5)
    }
  }

  flags <- detect_normalisation(result$original_filenames)

  tibble::tibble(
    GSE_ID = gse_id,
    success = success,
    attempts = attempts,
    original_filenames = paste(result$original_filenames, collapse = "; "),
    is_normalised = flags$is_normalised,
    is_tmm = flags$is_tmm,
    error_message = ifelse(is.na(result$error), "", result$error)
  )
}

counts_log <- lapply(
  seq_datalist$GSE_ID,
  download_counts_with_retries,
  output_dir = dirs$counts_raw,
  max_attempts = 3
) |>
  dplyr::bind_rows()

counts_log_path <- file.path(dirs$logs, "counts_download_log.csv")

readr::write_csv(
  counts_log,
  counts_log_path
)

message("Counts download log saved to: ", counts_log_path)


# ================================================================
# 5. TMM normalisation
# ================================================================

normalise_tmm <- function(file_path, output_dir) {
  dataset_id <- tools::file_path_sans_ext(basename(file_path))

  message("TMM normalising: ", dataset_id)

  tryCatch({
    counts <- read.csv(
      file_path,
      row.names = 1,
      check.names = FALSE
    )

    counts <- counts[, vapply(counts, is.numeric, logical(1)), drop = FALSE]

    if (nrow(counts) == 0 || ncol(counts) == 0) {
      stop("No numeric count matrix detected.")
    }

    if (anyNA(counts)) {
      stop("Count matrix contains missing values.")
    }

    if (any(counts < 0)) {
      stop("Count matrix contains negative values.")
    }

    if (any(counts %% 1 != 0)) {
      message("Warning: count matrix contains non-integer values.")
    }

    dge <- edgeR::DGEList(counts = counts)

    keep <- edgeR::filterByExpr(dge)

    if (sum(keep) < 5) {
      stop("Too few genes retained after filtering.")
    }

    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- edgeR::calcNormFactors(dge, method = "TMM")

    log2_tmm_cpm <- edgeR::cpm(
      dge,
      normalized.lib.sizes = TRUE,
      log = TRUE
    )

    output_file <- file.path(
      output_dir,
      paste0(dataset_id, "_log2_tmm_cpm.csv")
    )

    write.csv(
      as.data.frame(log2_tmm_cpm),
      output_file,
      row.names = TRUE
    )

    tibble::tibble(
      dataset_id = dataset_id,
      status = "success",
      n_genes_input = nrow(counts),
      n_genes_retained = sum(keep),
      n_samples = ncol(counts),
      output_file = output_file,
      message = "TMM normalisation successful."
    )
  }, error = function(e) {
    tibble::tibble(
      dataset_id = dataset_id,
      status = "failed",
      n_genes_input = NA_integer_,
      n_genes_retained = NA_integer_,
      n_samples = NA_integer_,
      output_file = NA_character_,
      message = conditionMessage(e)
    )
  })
}

count_files <- list.files(
  dirs$counts_raw,
  pattern = "\\.csv$",
  full.names = TRUE
)

if (length(count_files) == 0) {
  stop(
    "No raw count CSV files found in inputs/counts_data/raw_counts_data.\n",
    "Check outputs/logs/counts_download_log.csv for download failures."
  )
}

tmm_log <- lapply(
  count_files,
  normalise_tmm,
  output_dir = dirs$tmm_raw
) |>
  dplyr::bind_rows()

tmm_log_path <- file.path(dirs$logs, "tmm_normalisation_log.csv")

readr::write_csv(
  tmm_log,
  tmm_log_path
)

message("TMM normalisation log saved to: ", tmm_log_path)


# ================================================================
# 6. Pipeline complete
# ================================================================

message(
  "\n============================================================\n",
  "PRE-PROCESSING PIPELINE COMPLETE\n",
  "============================================================\n\n",
  "Review the following logs:\n",
  "  outputs/logs/counts_download_log.csv\n",
  "  outputs/logs/tmm_normalisation_log.csv\n\n",
  "Raw count CSV files are saved in:\n",
  "  inputs/counts_data/raw_counts_data/\n\n",
  "Log2 TMM CPM files are saved in:\n",
  "  inputs/tmm_data/raw_tmm_data/\n",
  "============================================================\n"
)
