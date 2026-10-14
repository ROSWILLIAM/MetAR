# MetAR: Meta-Analysis of RNA-seq based associations

MetAR is a reproducible R-based pipeline for performing meta-analysis of RNA-seq datasets from the NCBI Gene Expression Omnibus (GEO). It enables analysis of relationships between a user-defined variable (e.g. age, hormone levels) and gene expression across multiple studies.

This guide provides a clear, step-by-step walkthrough for first-time users.

---------------------------------------------------------------------

## OVERVIEW

MetAR consists of three scripts:

1. 01_pdata.R        	→ Downloads phenotype (pData)
2. 02_counts_data.R     → Processes RNA-seq count data
3. 03_MetAR.R     	→ Runs modelling and meta-analysis

IMPORTANT:
MetAR does NOT select datasets for you.
You must first identify potentially suitable GEO datasets manually through systematic searching.

---------------------------------------------------------------------

## STEP 0 — IDENTIFY DATASETS (MANUAL)

Before using MetAR, you must perform a manual or systematic search of GEO.

Go to:
https://www.ncbi.nlm.nih.gov/geo/

Search using relevant terms, for example:
- "muscle AND age"

Select datasets based on:
- RNA-seq (NOT microarray)
- Availability of raw counts
- Availability of phenotype data (e.g. age)
- Sufficient sample size

OUTPUT:
Create a list of GSE IDs.

Example:

GSE_ID
GSE111010
GSE113165
GSE125892

This will be imported into MetAR in script 01_pdata.R
---------------------------------------------------------------------

## STEP 1 — DOWNLOAD PHENOTYPE DATA

Run in R:

source("scripts/01_pdata.R")

This downloads metadata for each dataset and saves it to:

inputs/pdata/raw_pdata/

---------------------------------------------------------------------

## STEP 2 — CLEAN PHENOTYPE DATA (MANUAL)

You must manually clean and standardise the phenotype data.

For each dataset:
- Open files in inputs/pdata/raw_pdata/
- Identify relevant variables (e.g. age)
- Standardise column names across datasets
- Convert values to numeric where needed
- Remove missing values

Save cleaned files to:

inputs/pdata/proc_pdata/

Each file MUST contain:
- Sample ID
- Variable of interest (e.g. age)

---------------------------------------------------------------------

## STEP 3 — PROCESS RNA-SEQ COUNT DATA

Run:

source("scripts/02_counts_data.R")

This performs:
- Filtering of low-expression genes
- TMM normalisation
- Conversion to log2 TMM-CPM

---------------------------------------------------------------------

## STEP 4 — VERIFY DATA ALIGNMENT (MANUAL)

Check:
- Sample IDs match between counts and pData
- No duplicates
- Variable of interest is numeric

If this step is incorrect, the analysis will fail.

---------------------------------------------------------------------

## STEP 5 — DEFINE TARGET GENES

Create:

inputs/target_gene_list.csv

Format:

symbol,ensembl_id,category
AR,ENSG00000169083,AR signalling
STAR,ENSG00000177706,Steroidogenic enzyme

Notes:
- symbol = gene name
- ensembl_id = preferred for matching
- category = used for plotting/grouping

---------------------------------------------------------------------

## STEP 6 — SET VARIABLE OF INTEREST

Open:

scripts/03_MetAR.R

Edit the following:

variable_of_interest <- "age"
variable_label <- "Age"
variable_units <- "years"
effect_scale <- 10

Examples:
- Age → effect_scale = 10 (per decade)
- Hormone → effect_scale = 1

---------------------------------------------------------------------

## STEP 7 — RUN META-ANALYSIS

Run:

source("scripts/03_MetAR.R")

This will:
- Fit models within each dataset
- Compare LM vs GAM
- Perform random-effects meta-analysis
- Generate all publication ready plots and tables

---------------------------------------------------------------------

## OUTPUTS

Outputs are saved in:

outputs/

Key folders:

dataset_level/
- Model outputs per dataset

meta_analysis/
- Pooled results

plots/
- Forest plots
- Funnel plots
- Jungle plots
- Baujat plots
- Scatterplots

tables/
- model summary tables
- publication ready results tables

logs/
- Run summaries

---------------------------------------------------------------------

## COMMON ERRORS

Problem: Variable not found
→ Check column name in pData

Problem: Empty results
→ Too few samples

Problem: Sample mismatch
→ IDs not aligned

---------------------------------------------------------------------

## EXAMPLE APPLICATION

Williams et al. (2026)
MetAR: A semi-automated meta-analysis of skeletal muscle androgen receptors association with age
https://www.biorxiv.org/content/10.64898/2026.04.21.719741v1

---------------------------------------------------------------------
