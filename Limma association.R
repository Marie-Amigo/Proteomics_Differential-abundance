# =============================================================
# Unified Disease-Associated Proteomics Pipeline
#
# PURPOSE
# -------
# This script performs a differential protein abundance analysis
# per condition in comparison with control.
# -------
# Statistical framework:
#   - limma pairwise contrasts
#   - optional covariates
#   - presence filtering
#   - LOSO stability
#   - label-randomization permutation sanity check
#   - heatmap-only imputation for visualisation
#   - correctness
#   - consistency
#   - robustness
#   - reproducibility
#   - auditability
#   - code cleanliness
#
# IMPORTANT:
# Change the toggle accordingly to your dataset
# =============================================================

# ========================= 0) SETUP ==========================

options(stringsAsFactors = FALSE)
set.seed(123)

# ---------- Package loading ----------
required_pkgs <- c(
  "conflicted",
  "readxl", "dplyr", "tidyr", "stringr", "tibble", "purrr",
  "limma", "openxlsx", "ggplot2", "ggpubr", "RColorBrewer",
  "jsonlite", "digest", "readr", "pheatmap", "tools", "ggrepel"
)

to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)

invisible(lapply(required_pkgs, require, character.only = TRUE))

if ("package:conflicted" %in% search()) {
  conflicted::conflicts_prefer(dplyr::desc)
  conflicted::conflicts_prefer(base::unname)
  conflicted::conflicts_prefer(dplyr::intersect)
  conflicted::conflicts_prefer(dplyr::setdiff)
  conflicted::conflicts_prefer(dplyr::filter)
  conflicted::conflicts_prefer(dplyr::lag)
  conflicted::conflicts_prefer(dplyr::select)
  conflicted::conflicts_prefer(base::union)
}

# ====================== 1) USER SETTINGS =====================

# ---------------------- INPUT FILES --------------------------
# CHANGE THESE to your current input files / sheets
input_path       <- "C:/Users/90957842/Documents/20250730_R scripts/00_Inputs/20260330_DAP_MS Biomarkers.xlsx"
expression_sheet <- "DAP_matrix"
metadata_sheet   <- "Metadata"

# FASTA used only for accession -> gene symbol mapping
# CHANGE to your local FASTA if needed
gene_map_fasta   <- "C:/Users/90957842/Documents/20250730_R scripts/03_Library/Uniprot_Human_revised-FINAL.fasta"

# ---------------------- OUTPUT PATHS -------------------------
# CHANGE output_dir if you want outputs elsewhere
output_dir <- "C:/Users/90957842/Documents/20250730_R scripts/02_Processed/Disease correlations"

# -------------------- DATA STATUS TOGGLES --------------------
# Set TRUE if your matrix is already log2-transformed
already_log2 <- TRUE

# Set TRUE if your matrix is already normalized
already_normalised <- TRUE

# IMPORTANT:
# Set TRUE only if zeros in your input matrix mean "missing / undetected"
# Set FALSE if zero can be a real numeric value in your processed matrix
zero_is_missing <- TRUE

# ---------------- PRESENCE FILTER TOGGLES --------------------
# Primary presence threshold used for the main analysis
min_presence_frac <- 0.5

# Thresholds used for sensitivity analysis (Jaccard comparison)
presence_thresholds <- sort(unique(c(min_presence_frac, 0.30, 0.50, 0.60, 0.70)))

# Keep TRUE for group-aware presence filtering
require_presence_by_group <- TRUE

# Presence mode options:
#   "all_groups" = protein must pass threshold in every group
#   "any_group"  = protein must pass threshold in at least one group
#   "overall"    = threshold applied across all samples pooled
presence_mode <- "all_group"

# ----------------- GROUP / CONTRAST TOGGLES ------------------
# Reference / baseline group label as written in metadata$Condition
control_label <- "Control"

# TRUE = also test disease-vs-disease contrasts
# FALSE = only disease-vs-control contrasts
include_disease_vs_disease <- TRUE

# Not implemented in this script version
paired_mode <- FALSE

# ------------------ COVARIATE TOGGLES ------------------------
# TRUE = include covariates in limma model
use_covariates <- TRUE

# TRUE = drop samples with missing covariates before modelling
# Recommended for clean modelling
drop_missing_covariates <- TRUE

# CHANGE covariates here if needed
# Example: c("Age", "Sex")
covariates <- c("Age", "Sex")

# ---------------- SIGNIFICANCE TOGGLES -----------------------
alpha_fdr_primary   <- 0.05
alpha_p_exploratory <- 0.05

# Optional secondary global BH over all Protein x Contrast tests
add_global_FDR <- FALSE

# ---------------------- QC TOGGLES ---------------------------
run_QC     <- FALSE
run_LOSO   <- FALSE
make_figures <- TRUE

# Heatmap z-score cap
HEATMAP_ZLIM <- 3

# ------------------ PERMUTATION TOGGLES ----------------------
# This is a label-randomization benchmark / sanity check
# It is NOT a covariate-preserving permutation test
run_permutation <- FALSE
n_perm          <- 1000
perm_alpha      <- 0.05

# ------------------- VOLCANO LABEL TOGGLES -------------------
volcano_label_mode <- "top"   # "top" or "all_fdr"
volcano_label_top  <- 25
volcano_label_fdr  <- 0.05
volcano_label_cap  <- Inf

# -------------------- FOREST TOGGLES -------------------------
use_only_significant_forest <- TRUE
max_forest_labels <- 60

# -------------------- BOXPLOT TOGGLES ------------------------
# Display order for groups on boxplots
# Edit this if your project uses a different preferred order
box_order <- c("Control", "RRMS", "Before", "After", "Stable", "SPMS", "R", "NR", "Fibromyalgia", "RA", "CPC")

# --------------------- SAFETY CHECKS -------------------------
if (paired_mode) {
  stop("paired_mode=TRUE is not implemented in this script version.")
}

# ======================= 2) PATH SETUP =======================

date_str  <- format(Sys.Date(), "%Y-%m-%d")
infile_nm <- tools::file_path_sans_ext(basename(input_path))
base_dir  <- file.path(output_dir, date_str, infile_nm)

vis_dir        <- file.path(base_dir, "Visualisations")
heatmap_dir    <- file.path(vis_dir, "Heatmap")
boxplot_dir    <- file.path(vis_dir, "Boxplot")
forest_dir     <- file.path(vis_dir, "Forest")
qc_dir         <- file.path(base_dir, "QC")
perm_dir       <- file.path(base_dir, "Permutation")
volcano_dir    <- file.path(vis_dir, "Volcano")

box_p_sig_dir    <- file.path(vis_dir, "Boxplot_P_lt005")
box_fdr_sig_dir  <- file.path(vis_dir, "Boxplot_FDR_lt005")
box_fdr_near_dir <- file.path(vis_dir, "Boxplot_FDR_lt01")

dir_list <- c(
  base_dir, vis_dir, heatmap_dir, boxplot_dir, forest_dir,
  qc_dir, perm_dir, volcano_dir,
  box_p_sig_dir, box_fdr_sig_dir, box_fdr_near_dir
)

for (p in dir_list) dir.create(p, recursive = TRUE, showWarnings = FALSE)

# ======================= 3) HELPERS ==========================

as_one_string <- function(x, fallback = "") {
  if (is.null(x) || length(x) != 1) return(fallback)
  x <- as.character(x)
  if (is.na(x) || !nzchar(x)) return(fallback)
  x
}

safe_stub <- function(x, fallback = "item") {
  x <- as_one_string(x, fallback = fallback)
  if (!nzchar(x)) x <- fallback
  y <- gsub("[^A-Za-z0-9]+", "_", x)
  y <- gsub("^_+|_+$", "", y)
  if (!nzchar(y)) y <- fallback
  substr(y, 1, 120)
}

safe_path <- function(...) {
  p <- file.path(...)
  p <- as_one_string(p, fallback = "")
  if (!nzchar(p)) stop("Invalid output path produced (empty/NA).")
  p
}

safe_ggsave <- function(filename, plot, ...) {
  filename <- as_one_string(filename, fallback = "")
  if (!nzchar(filename)) stop("Invalid filename for ggsave (empty/NA).")
  try(ggplot2::ggsave(filename = filename, plot = plot, ..., limitsize = FALSE), silent = TRUE)
}

robust_save_png <- function(plot_obj, file, width = 3.5, height = 4.5, dpi = 300, bg = "white") {
  file <- as_one_string(file, fallback = "")
  if (!nzchar(file)) return(invisible(NULL))
  try(
    ggplot2::ggsave(file, plot_obj, width = width, height = height, dpi = dpi, bg = bg, limitsize = FALSE),
    silent = TRUE
  )
}

excel_safe_sheet_name <- function(x, existing = character(0), max_len = 31) {
  nm <- gsub("[:\\\\/\\?\\*\\[\\]]", "_", as_one_string(x, "Sheet"))
  nm <- gsub("\\s+", " ", nm)
  nm <- trimws(nm)
  if (!nzchar(nm)) nm <- "Sheet"
  nm <- substr(nm, 1, max_len)
  
  base_nm <- nm
  k <- 1L
  while (nm %in% existing) {
    suffix <- paste0("_", k)
    nm <- substr(base_nm, 1, max_len - nchar(suffix))
    nm <- paste0(nm, suffix)
    k <- k + 1L
  }
  nm
}

wb_add_autosized <- function(wb, sheet, df) {
  existing <- unique(c(
    tryCatch(names(wb$worksheets), error = function(e) character(0)),
    tryCatch(wb$sheet_names,       error = function(e) character(0))
  ))
  existing <- existing[nzchar(existing)]
  
  safe_nm <- excel_safe_sheet_name(sheet, existing = existing)
  openxlsx::addWorksheet(wb, safe_nm)
  openxlsx::writeData(wb, safe_nm, df)
  
  if (ncol(df) > 0) {
    openxlsx::setColWidths(wb, safe_nm, cols = 1:ncol(df), widths = "auto")
  }
}

keep_positive_sd <- function(mat) {
  sdv <- apply(mat, 1, function(x) stats::sd(x, na.rm = TRUE))
  is.finite(sdv) & sdv > 0
}

recode_sex <- function(x) {
  x <- trimws(as.character(x))
  x <- toupper(x)
  x <- dplyr::recode(
    x,
    "FEMALE" = "F", "WOMAN" = "F", "WOMEN" = "F", "F" = "F",
    "MALE"   = "M", "MAN"   = "M", "MEN"   = "M", "M" = "M",
    .default = x
  )
  x <- ifelse(x %in% c("F", "M"), x, NA_character_)
  x
}

assign_condition_colors <- function(condition,
                                    control_label = "Control",
                                    overrides = NULL,
                                    control_color = "#9CA3AF",
                                    fallback_color = "#BDBDBD") {
  levs <- levels(condition)
  disease_levels <- setdiff(levs, control_label)
  disease_map <- if (!is.null(overrides)) overrides else character(0)
  disease_map <- disease_map[intersect(names(disease_map), disease_levels)]
  
  missing <- setdiff(disease_levels, names(disease_map))
  if (length(missing) > 0) {
    disease_map <- c(disease_map, setNames(rep(fallback_color, length(missing)), missing))
  }
  
  c(setNames(control_color, control_label), disease_map)
}

# ================== 4) COLOUR DEFINITIONS ====================

overrides <- c(
  "RRMS" = "#FF69B4",
  "Fibromyalgia" = "#BCEE68",
  "RA" = "#FF7F00",
  "SPMS" = "#ADD8E6",
  "Resp" = "#8B0000",
  "Before" = "#8B3A62",
  "After" = "#FF69B4",
  "Active" = "#FF69B4",
  "Stable" = "#D15FEE",
  "CPC" = "#825740",
  "R" = "#D15FEE",
  "NR" = "#ADD8E6",
  
  "Phenotype_1" = "#FF69B4",
  "Phenotype_2" = "#BCEE68",
  "Phenotype_3" = "#FF7F00",
  "Phenotype_4" = "#ADD8E6"
  )

# ==================== 5) LOAD INPUT DATA =====================

if (!file.exists(input_path)) stop("Input file not found: ", input_path)
if (!file.exists(gene_map_fasta)) stop("FASTA file not found: ", gene_map_fasta)

ncols_expr <- ncol(readxl::read_excel(input_path, sheet = expression_sheet, n_max = 0))
col_types  <- c("text", rep("numeric", max(0, ncols_expr - 1)))

expr_raw <- suppressMessages(
  readxl::read_excel(
    input_path,
    sheet = expression_sheet,
    col_types = col_types,
    col_names = TRUE,
    na = c("", "NA", "NaN", "nan", "#N/A")
  )
)

meta <- suppressMessages(
  readxl::read_excel(input_path, sheet = metadata_sheet, col_names = TRUE)
) %>% as.data.frame()

if (!"Protein" %in% names(expr_raw)) stop("Expression sheet must contain a 'Protein' column.")
if (!all(c("Sample_ID", "Condition") %in% names(meta))) {
  stop("Metadata sheet must contain columns: Sample_ID and Condition")
}

# ==================== 6) BASIC CLEANING ======================

expr_raw$Protein <- trimws(as.character(expr_raw$Protein))
expr_raw$Protein[expr_raw$Protein == ""] <- NA

meta$Sample_ID  <- trimws(as.character(meta$Sample_ID))
meta$Condition  <- trimws(as.character(meta$Condition))

if (anyDuplicated(meta$Sample_ID)) {
  stop("Duplicate Sample_ID values found in metadata.")
}

# Drop completely empty rows with missing Protein and all missing numeric values
num_cols <- setdiff(names(expr_raw), "Protein")
all_num_na <- if (length(num_cols)) {
  apply(expr_raw[, num_cols, drop = FALSE], 1, function(z) all(is.na(z)))
} else {
  rep(TRUE, nrow(expr_raw))
}

drop_idx <- which(is.na(expr_raw$Protein) & all_num_na)
if (length(drop_idx)) expr_raw <- expr_raw[-drop_idx, , drop = FALSE]

# Keep original accession, make rownames unique for modelling
expr_wide <- expr_raw %>%
  mutate(
    Protein_Accession = Protein,
    Protein = make.unique(Protein_Accession, sep = "__dup")
  ) %>%
  relocate(Protein_Accession, .after = Protein)

if (anyDuplicated(expr_wide$Protein)) {
  stop("Duplicate Protein keys after make.unique() -- unexpected.")
}

expr0 <- expr_wide[, setdiff(names(expr_wide), c("Protein", "Protein_Accession")), drop = FALSE] %>%
  as.data.frame()

rownames(expr0) <- expr_wide$Protein
expr0[] <- lapply(expr0, function(x) suppressWarnings(as.numeric(as.character(x))))

# Standardise non-finite values to NA
nonfinite_before <- sum(!is.finite(as.matrix(expr0)), na.rm = TRUE)
expr0[] <- lapply(expr0, function(v) {
  v[!is.finite(v)] <- NA_real_
  v
})
nonfinite_after <- sum(!is.finite(as.matrix(expr0)), na.rm = TRUE)
nonfinite_converted <- nonfinite_before - nonfinite_after

cat("Rows:", nrow(expr_wide), " | Cols:", ncol(expr_wide), "\n")

# Optional cohort audit
if ("Cohort" %in% names(meta)) {
  meta$Cohort <- trimws(as.character(meta$Cohort))
  cat("\n=== Condition x Cohort table ===\n")
  print(table(meta$Condition, meta$Cohort, useNA = "ifany"))
  
  cat("\n=== Row-wise proportions (within Condition) ===\n")
  print(round(prop.table(table(meta$Condition, meta$Cohort), margin = 1), 3))
} else {
  message("No 'Cohort' column found in metadata.")
}

# Align samples
expr_cols <- colnames(expr0)

if (length(setdiff(meta$Sample_ID, expr_cols))) {
  warning("Dropping metadata samples not in expression: ",
          paste(setdiff(meta$Sample_ID, expr_cols), collapse = ", "))
}

if (length(setdiff(expr_cols, meta$Sample_ID))) {
  warning("Expression columns with no metadata (dropped): ",
          paste(setdiff(expr_cols, meta$Sample_ID), collapse = ", "))
}

keep_ids <- intersect(meta$Sample_ID, expr_cols)
if (!length(keep_ids)) stop("No overlapping samples between expression and metadata.")

meta  <- meta[match(keep_ids, meta$Sample_ID), , drop = FALSE]
expr0 <- expr0[, keep_ids, drop = FALSE]
expr0 <- expr0[, meta$Sample_ID, drop = FALSE]

# ================= 7) LOG / ZERO / COVARIATES ================

expr <- expr0
log2_applied <- FALSE

if (!already_log2) {
  expr[] <- lapply(expr, function(x) ifelse(is.na(x), NA_real_, log2(x + 1)))
  log2_applied <- TRUE
} else {
  message("Data already log2-transformed -- skipped.")
}

nan_to_NA_count <- sum(unlist(lapply(expr, function(v) sum(is.nan(v), na.rm = TRUE))), na.rm = TRUE)
expr[] <- lapply(expr, function(v) {
  v[is.nan(v)] <- NA_real_
  v
})

if (zero_is_missing) {
  zeros_to_NA_count <- sum(unlist(lapply(expr, function(v) sum(v == 0, na.rm = TRUE))), na.rm = TRUE)
  expr[] <- lapply(expr, function(v) {
    v[v == 0] <- NA_real_
    v
  })
  message("Zero values were treated as missing and converted to NA.")
} else {
  zeros_to_NA_count <- 0
  message("Zero values were NOT converted to NA.")
}

cat("Converted", zeros_to_NA_count, "zeros and", nan_to_NA_count, "NaNs to NA.\n")

# Covariate-ready metadata
meta_model <- meta %>%
  dplyr::mutate(
    Age = if ("Age" %in% names(meta)) as.numeric(Age) else NA_real_,
    Sex = if ("Sex" %in% names(meta)) factor(recode_sex(Sex), levels = c("F", "M")) else factor(NA, levels = c("F", "M"))
  )

if (use_covariates) {
  missing_covariate_cols <- setdiff(covariates, names(meta_model))
  if (length(missing_covariate_cols)) {
    stop("Missing required covariate columns in metadata: ",
         paste(missing_covariate_cols, collapse = ", "))
  }
  
  if (drop_missing_covariates) {
    cc <- complete.cases(meta_model[, covariates, drop = FALSE])
    dropped_samples <- meta_model$Sample_ID[!cc]
    if (any(!cc)) {
      message("Covariates complete-case: dropping ", sum(!cc), " samples with missing covariates.")
    }
    keep_samps <- meta_model$Sample_ID[cc]
  } else {
    dropped_samples <- character(0)
    keep_samps <- meta_model$Sample_ID
    message("Covariates ON but drop_missing_covariates=FALSE.")
  }
} else {
  dropped_samples <- character(0)
  keep_samps <- meta_model$Sample_ID
  message("Covariates OFF: model uses Condition only.")
}

if (length(keep_samps) < 3) stop("Fewer than 3 samples remain after covariate filtering.")

expr_sub <- expr[, keep_samps, drop = FALSE]
meta_sub <- meta_model[match(keep_samps, meta_model$Sample_ID), , drop = FALSE]

condition0 <- factor(meta_sub$Condition)
if (control_label %in% levels(condition0)) {
  condition0 <- stats::relevel(condition0, ref = control_label)
}
names(condition0) <- keep_samps

box_levels <- box_order[box_order %in% levels(condition0)]
if (!length(box_levels)) box_levels <- levels(condition0)

condition <- factor(as.character(condition0), levels = box_levels)
names(condition) <- keep_samps

if (anyNA(condition)) stop("Condition contains NA after factor handling.")

# Warn if covariate degenerates after filtering
if (use_covariates) {
  for (cv in covariates) {
    vv <- meta_sub[[cv]]
    if (is.factor(vv) || is.character(vv)) {
      nlev <- length(unique(na.omit(vv)))
      if (nlev < 2) {
        warning("Covariate '", cv, "' has <2 non-missing levels after filtering.")
      }
    }
  }
}

condition_colors <- assign_condition_colors(
  condition,
  control_label = control_label,
  overrides = overrides
)

# ================ 8) PRESENCE FILTER FUNCTIONS ===============

apply_presence_filter <- function(expr_in, condition, frac, by_group = TRUE,
                                  mode = c("all_groups", "any_group", "overall")) {
  mode <- match.arg(mode)
  
  condition <- condition[colnames(expr_in)]
  condition <- droplevels(condition)
  
  stopifnot(length(condition) == ncol(expr_in))
  stopifnot(identical(names(condition), colnames(expr_in)))
  stopifnot(!anyNA(condition))
  
  if (mode == "overall" || !by_group) {
    req <- ceiling(ncol(expr_in) * frac)
    keep <- rowSums(!is.na(expr_in)) >= req
  } else {
    grps <- levels(condition)
    sizes <- table(condition)[grps]
    req <- ceiling(as.numeric(sizes) * frac)
    
    keep_by_grp <- sapply(seq_along(grps), function(i) {
      idx <- which(condition == grps[i])
      rowSums(!is.na(expr_in[, idx, drop = FALSE])) >= req[i]
    })
    
    if (is.vector(keep_by_grp)) {
      keep_by_grp <- matrix(keep_by_grp, ncol = length(grps))
    }
    
    if (mode == "all_groups") {
      keep <- apply(keep_by_grp, 1, all)
    } else {
      keep <- apply(keep_by_grp, 1, any)
    }
  }
  
  expr_in[keep, , drop = FALSE]
}

# =================== 9) MAIN FILTERING =======================

n_start <- nrow(expr_sub)

expr_pf <- apply_presence_filter(
  expr_sub,
  condition,
  frac = min_presence_frac,
  by_group = require_presence_by_group,
  mode = presence_mode
)

n_after_presence <- nrow(expr_pf)

keep_var <- apply(expr_pf, 1, function(x) {
  sdv <- stats::sd(x, na.rm = TRUE)
  is.finite(sdv) && sdv > 0
})

expr_pf <- expr_pf[keep_var, , drop = FALSE]
n_after_variance <- nrow(expr_pf)

expr_norm <- if (already_normalised) {
  expr_pf
} else {
  meds <- apply(expr_pf, 2, stats::median, na.rm = TRUE)
  sweep(expr_pf, 2, meds, "-")
}

expr_mat <- as.matrix(expr_norm)

if (nrow(expr_mat) < 2) stop("Too few proteins remain after filtering.")
if (ncol(expr_mat) < 3) stop("Too few samples remain for modelling.")

# Heatmap-only imputation
expr_imp_viz <- as.data.frame(expr_mat)
for (g in levels(condition)) {
  idx <- which(condition == g)
  grp <- expr_imp_viz[, idx, drop = FALSE]
  grp_meds <- apply(grp, 1, stats::median, na.rm = TRUE)
  miss <- is.na(grp)
  if (any(miss)) grp[miss] <- grp_meds[row(miss)[miss]]
  expr_imp_viz[, idx] <- grp
}

saveRDS(
  expr_imp_viz,
  safe_path(base_dir, paste0(date_str, "_expr_for_heatmap_visualisation_only.rds"))
)

# =================== 10) FASTA GENE MAPPING ==================

parse_fasta_gene_map <- function(fasta_path) {
  fasta_path <- as_one_string(fasta_path, "")
  if (!nzchar(fasta_path) || !file.exists(fasta_path)) {
    stop("FASTA file not found: ", fasta_path)
  }
  
  hdrs <- readLines(fasta_path, warn = FALSE)
  hdrs <- hdrs[startsWith(hdrs, ">")]
  if (!length(hdrs)) stop("No FASTA headers found in ", fasta_path)
  
  get_acc <- function(h) {
    m <- regexec("^>[^|]*\\|([^|]+)\\|", h)
    r <- regmatches(h, m)[[1]]
    if (length(r) >= 2) return(r[2])
    sub("^>(\\S+).*", "\\1", h)
  }
  
  get_gene <- function(h) {
    g1 <- sub(".*\\bGN=([A-Za-z0-9._-]+)\\b.*", "\\1", h)
    if (!identical(g1, h)) return(g1)
    g2 <- sub(".*\\bGene=([A-Za-z0-9._-]+)\\b.*", "\\1", h)
    if (!identical(g2, h)) return(g2)
    NA_character_
  }
  
  acc  <- vapply(hdrs, get_acc,  character(1))
  gene <- vapply(hdrs, get_gene, character(1))
  
  acc <- sub("-\\d+$", "", acc)
  keep <- !is.na(gene) & nzchar(gene)
  acc  <- acc[keep]
  gene <- gene[keep]
  
  dupe <- duplicated(acc)
  setNames(gene[!dupe], acc[!dupe])
}

acc_for_map <- sub("-\\d+$", "", expr_wide$Protein_Accession)
names(acc_for_map) <- expr_wide$Protein

fasta_map <- parse_fasta_gene_map(gene_map_fasta)
gene_vec  <- base::unname(fasta_map[acc_for_map])

id_map <- tibble::tibble(
  Protein   = expr_wide$Protein,
  Accession = expr_wide$Protein_Accession,
  Gene      = gene_vec
)

message(sprintf("FASTA mapping: %d/%d accessions mapped to Gene.",
                sum(!is.na(id_map$Gene)), nrow(id_map)))

label_from_map_base <- function(keys, id_map) {
  df <- id_map[match(keys, id_map$Protein), c("Gene", "Accession")]
  gene <- as.character(df$Gene)
  acc  <- as.character(df$Accession)
  
  ifelse(!is.na(gene) & nzchar(gene), gene,
         ifelse(!is.na(acc) & nzchar(acc), acc, keys))
}

label_from_map_unique <- function(keys, id_map) {
  base_lab <- label_from_map_base(keys, id_map)
  acc <- id_map$Accession[match(keys, id_map$Protein)]
  dup <- duplicated(base_lab) | duplicated(base_lab, fromLast = TRUE)
  base_lab[dup] <- paste0(base_lab[dup], " [", acc[dup], "]")
  make.unique(base_lab, sep = " ")
}

# ==================== 11) LIMMA FUNCTIONS ====================

build_contrasts <- function(condition, control_label, include_dd = TRUE) {
  levs <- levels(condition)
  if (length(levs) < 2) return(character(0))
  
  if (!(control_label %in% levs)) {
    warning("control_label not found in Condition levels; using first level as baseline.")
    control_label <- levs[1]
  }
  
  dis <- setdiff(levs, control_label)
  out <- character(0)
  
  for (g in dis) out <- c(out, paste0(g, "_vs_", control_label))
  
  if (include_dd && length(dis) >= 2) {
    cmb <- combn(dis, 2, simplify = FALSE)
    for (p in cmb) out <- c(out, paste0(p[1], "_vs_", p[2]))
  }
  
  unique(out)
}

build_design0 <- function(condition, meta, covariates = NULL, use_covariates = TRUE) {
  condition <- droplevels(condition)
  mm <- meta[match(names(condition), meta$Sample_ID), , drop = FALSE]
  stopifnot(all(mm$Sample_ID == names(condition)))
  
  if (isTRUE(use_covariates) && !is.null(covariates) && length(covariates)) {
    f <- stats::as.formula(paste("~ 0 + condition +", paste(covariates, collapse = " + ")))
  } else {
    f <- ~ 0 + condition
  }
  
  d <- stats::model.matrix(f, data = cbind(mm, condition = condition))
  colnames(d) <- make.names(colnames(d))
  rownames(d) <- names(condition)
  d
}

run_limma <- function(X, condition, meta, covariates, control_label,
                      include_dd = TRUE, use_covariates = TRUE) {
  
  condition <- droplevels(condition)
  design <- build_design0(condition, meta, covariates, use_covariates = use_covariates)
  fit0 <- limma::lmFit(X, design)
  
  res_df <- fit0$df.residual
  ok_df <- is.finite(res_df) & res_df > 1
  if (any(!ok_df)) {
    message("Dropping ", sum(!ok_df), " proteins with insufficient residual df before eBayes.")
    fit0 <- fit0[ok_df, ]
  }
  
  fit0e <- limma::eBayes(fit0, trend = TRUE, robust = TRUE)
  
  cond_cols <- grep("^condition", colnames(design))
  
  if (length(cond_cols) < 2) {
    omni0 <- limma::topTable(fit0e, coef = cond_cols, number = Inf,
                             adjust.method = "BH", sort.by = "P")
  } else {
    omni0 <- limma::topTable(fit0e, coef = cond_cols, number = Inf,
                             adjust.method = "BH", sort.by = "F")
  }
  
  omnibus <- omni0 %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Protein")
  
  if ("F" %in% names(omnibus))         { omnibus$Omnibus_F   <- omnibus$F;         omnibus$F <- NULL }
  if ("t" %in% names(omnibus))         { omnibus$Omnibus_t   <- omnibus$t;         omnibus$t <- NULL }
  if ("P.Value" %in% names(omnibus))   { omnibus$Omnibus_P   <- omnibus$P.Value;   omnibus$P.Value <- NULL }
  if ("adj.P.Val" %in% names(omnibus)) { omnibus$Omnibus_FDR <- omnibus$adj.P.Val; omnibus$adj.P.Val <- NULL }
  
  contrast_names <- build_contrasts(condition, control_label, include_dd)
  if (!length(contrast_names)) {
    return(list(
      omnibus = omnibus,
      results = list(),
      fit2 = NULL,
      contrast_mat = NULL,
      design = design
    ))
  }
  
  coef_for_level <- function(level) {
    nm <- make.names(paste0("condition", level))
    if (!nm %in% colnames(design)) {
      stop("Could not find coefficient for level '", level, "' -> '", nm, "'.")
    }
    nm
  }
  
  contrast_terms <- vapply(contrast_names, function(nm) {
    parts <- strsplit(nm, "_vs_")[[1]]
    g1 <- parts[1]
    g2 <- parts[2]
    paste0(coef_for_level(g1), " - ", coef_for_level(g2))
  }, character(1))
  
  contrast_mat <- limma::makeContrasts(contrasts = contrast_terms, levels = design)
  fit2 <- limma::contrasts.fit(fit0, contrast_mat)
  fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)
  
  results <- list()
  for (i in seq_len(ncol(contrast_mat))) {
    nm <- contrast_names[i]
    
    tt <- limma::topTable(fit2, coef = i, number = Inf,
                          adjust.method = "BH", sort.by = "none") %>%
      tibble::rownames_to_column("Protein") %>%
      dplyr::rename(
        limma_logFC = logFC,
        limma_P     = P.Value,
        limma_FDR   = adj.P.Val,
        Limma_t     = t,
        B           = B
      )
    
    idx <- match(tt$Protein, rownames(fit2$coefficients))
    tt$limma_stdev_unscaled <- fit2$stdev.unscaled[, i][idx]
    tt$limma_s2_post        <- fit2$s2.post[idx]
    tt$limma_df_total       <- fit2$df.total[idx]
    
    tt$Contrast <- nm
    tt$Group1 <- sub("_vs_.*", "", nm)
    tt$Group2 <- sub(".*_vs_", "", nm)
    
    results[[nm]] <- tt
  }
  
  list(
    omnibus = omnibus,
    results = results,
    fit2 = fit2,
    contrast_mat = contrast_mat,
    design = design
  )
}

add_limma_ci <- function(df) {
  SE_post <- df$limma_stdev_unscaled * sqrt(df$limma_s2_post)
  tcrit <- ifelse(
    is.finite(df$limma_df_total) & df$limma_df_total > 0,
    stats::qt(0.975, df = df$limma_df_total),
    1.96
  )
  
  dplyr::mutate(
    df,
    limma_logFC_SE_post = SE_post,
    logFC_CI_low  = limma_logFC - tcrit * SE_post,
    logFC_CI_high = limma_logFC + tcrit * SE_post
  )
}

attach_labels <- function(df) {
  df %>%
    dplyr::left_join(id_map, by = "Protein") %>%
    dplyr::relocate(Accession, Gene, .after = Protein)
}

# ==================== 12) MAIN LIMMA RUN =====================

limma_pack <- run_limma(
  expr_mat,
  condition,
  meta_sub,
  covariates,
  control_label,
  include_dd = include_disease_vs_disease,
  use_covariates = use_covariates
)

omnibus_tbl    <- limma_pack$omnibus
limma_results  <- limma_pack$results
contrast_names <- names(limma_results)

saveRDS(
  list(
    design = limma_pack$design,
    contrast_mat = limma_pack$contrast_mat,
    use_covariates = use_covariates,
    covariates = if (use_covariates) covariates else character(0),
    drop_missing_covariates = drop_missing_covariates,
    dropped_samples_missing_covariates = dropped_samples
  ),
  safe_path(base_dir, paste0(date_str, "_limma_design_contrasts_audit.rds"))
)

merged_all_flagged <- if (length(limma_results)) {
  dplyr::bind_rows(lapply(limma_results, function(df) {
    out <- df %>%
      attach_labels() %>%
      dplyr::left_join(omnibus_tbl, by = "Protein") %>%
      dplyr::mutate(
        Is_Significant_FDR = is.finite(limma_FDR) & (limma_FDR < alpha_fdr_primary),
        Is_Significant_P   = is.finite(limma_P)   & (limma_P   < alpha_p_exploratory),
        Abs_log2FC         = abs(limma_logFC)
      )
    add_limma_ci(out)
  }))
} else {
  tibble::tibble()
}

if (add_global_FDR && nrow(merged_all_flagged)) {
  merged_all_flagged <- merged_all_flagged %>%
    dplyr::mutate(
      Global_BH_FDR = stats::p.adjust(limma_P, method = "BH"),
      Is_Significant_GlobalFDR = is.finite(Global_BH_FDR) & (Global_BH_FDR < alpha_fdr_primary)
    )
}

sig_hits_fdr <- merged_all_flagged %>%
  dplyr::filter(Is_Significant_FDR) %>%
  dplyr::arrange(limma_FDR, dplyr::desc(abs(limma_logFC)))

# ============ 13) PRESENCE SENSITIVITY (JACCARD) =============

sens_list <- list()

for (pf in presence_thresholds) {
  expr_pf2 <- apply_presence_filter(
    expr_sub,
    condition,
    frac = pf,
    by_group = require_presence_by_group,
    mode = presence_mode
  )
  
  expr_pf2 <- expr_pf2[keep_positive_sd(expr_pf2), , drop = FALSE]
  
  expr_pf2 <- if (already_normalised) {
    expr_pf2
  } else {
    meds <- apply(expr_pf2, 2, median, na.rm = TRUE)
    sweep(expr_pf2, 2, meds, "-")
  }
  
  lp_pf <- run_limma(
    as.matrix(expr_pf2),
    condition,
    meta_sub,
    covariates,
    control_label,
    include_dd = include_disease_vs_disease,
    use_covariates = use_covariates
  )
  
  sens_list[[as.character(pf)]] <- lp_pf$results
}

jaccard <- purrr::map_dfr(contrast_names, function(nm) {
  sets <- lapply(names(sens_list), function(pf) {
    df <- sens_list[[pf]][[nm]]
    if (is.null(df)) return(character(0))
    df %>%
      dplyr::filter(is.finite(limma_FDR) & limma_FDR < alpha_fdr_primary) %>%
      dplyr::pull(Protein) %>%
      unique()
  })
  names(sets) <- names(sens_list)
  
  base_set <- sets[[as.character(min_presence_frac)]]
  other_ps <- setdiff(names(sets), as.character(min_presence_frac))
  
  purrr::map_dfr(other_ps, function(pf) {
    s <- sets[[pf]]
    inter <- length(intersect(base_set, s))
    uni   <- length(union(base_set, s))
    
    tibble::tibble(
      Contrast = nm,
      Baseline = min_presence_frac,
      Compare  = as.numeric(pf),
      Jaccard  = ifelse(uni == 0, NA_real_, inter / uni)
    )
  })
})

# ================= 14) LOSO STABILITY CHECK ==================

if (run_LOSO && length(contrast_names)) {
  loso_flags <- list()
  
  for (nm in contrast_names) {
    parts <- strsplit(nm, "_vs_")[[1]]
    g1 <- parts[1]
    g2 <- parts[2]
    
    sam <- names(condition)[condition %in% c(g1, g2)]
    if (length(sam) < 4) next
    
    X_base <- expr_mat[, sam, drop = FALSE]
    m_base <- meta_sub[match(sam, meta_sub$Sample_ID), , drop = FALSE]
    c_base <- droplevels(condition[sam])
    
    lp_base <- run_limma(
      X_base, c_base, m_base, covariates,
      control_label = g2,
      include_dd = FALSE,
      use_covariates = use_covariates
    )
    
    target_nm <- paste0(g1, "_vs_", g2)
    if (!target_nm %in% names(lp_base$results)) next
    
    base_tab <- lp_base$results[[target_nm]] %>%
      dplyr::select(Protein, limma_logFC, limma_FDR) %>%
      dplyr::mutate(
        Is_Significant_FDR = is.finite(limma_FDR) & (limma_FDR < alpha_fdr_primary)
      )
    
    stability <- rep(TRUE, nrow(base_tab))
    names(stability) <- base_tab$Protein
    
    for (u in sam) {
      keep_u <- setdiff(sam, u)
      if (length(keep_u) < 4) next
      
      X_loo <- X_base[, keep_u, drop = FALSE]
      m_loo <- m_base[match(keep_u, m_base$Sample_ID), , drop = FALSE]
      c_loo <- droplevels(c_base[keep_u])
      
      lp_loo <- run_limma(
        X_loo, c_loo, m_loo, covariates,
        control_label = g2,
        include_dd = FALSE,
        use_covariates = use_covariates
      )
      
      if (!target_nm %in% names(lp_loo$results)) next
      
      tt_loo <- lp_loo$results[[target_nm]] %>%
        dplyr::select(Protein, limma_logFC, limma_FDR) %>%
        dplyr::rename(
          limma_logFC.y = limma_logFC,
          limma_FDR.y   = limma_FDR
        )
      
      chk <- base_tab %>% dplyr::left_join(tt_loo, by = "Protein")
      
      flip <- is.finite(chk$limma_logFC) &
        is.finite(chk$limma_logFC.y) &
        (sign(chk$limma_logFC) != sign(chk$limma_logFC.y))
      
      lost <- is.na(chk$limma_FDR.y) |
        (!is.finite(chk$limma_FDR.y)) |
        (chk$limma_FDR.y >= alpha_fdr_primary)
      
      instability <- chk$Is_Significant_FDR & (flip | lost)
      stability[chk$Protein] <- stability[chk$Protein] & !instability
    }
    
    loso_flags[[nm]] <- tibble::tibble(
      Protein = names(stability),
      Stable_LOSO = as.logical(stability),
      Contrast = nm
    )
  }
  
  loso_tbl <- dplyr::bind_rows(loso_flags)
  
  merged_all_flagged <- merged_all_flagged %>%
    dplyr::left_join(loso_tbl, by = c("Protein", "Contrast"))
}

# ===================== 15) EXPORT TABLES =====================

wb <- openxlsx::createWorkbook()

wb_add_autosized(wb, "AllMerged", merged_all_flagged)
wb_add_autosized(wb, "Omnibus_Ftest", omnibus_tbl)
if (nrow(jaccard)) wb_add_autosized(wb, "PresenceSensitivity_Jaccard", jaccard)

if (length(contrast_names)) {
  for (nm in contrast_names) {
    df <- merged_all_flagged %>%
      dplyr::filter(Contrast == nm) %>%
      dplyr::arrange(limma_FDR, dplyr::desc(abs(limma_logFC)))
    wb_add_autosized(wb, paste0(nm, "_limma"), df)
  }
}

wb_add_autosized(wb, "SignificantHits_FDR", sig_hits_fdr)

out_xlsx <- safe_path(base_dir, paste0(date_str, "_DiseaseAssoc_LIMMA_Results.xlsx"))
openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

readr::write_csv(sig_hits_fdr, safe_path(base_dir, paste0(date_str, "_SignificantHits_FDR.csv")))
readr::write_csv(omnibus_tbl,  safe_path(base_dir, paste0(date_str, "_Omnibus_Ftest.csv")))

# ======================= 16) FIGURES =========================

if (make_figures) {
  
  # -------------------- 16.1 HEATMAPS -----------------------
  make_heatmap <- function(protein_ids, title_text, file_stub) {
    protein_ids <- unique(protein_ids)
    
    if (length(protein_ids) < 2) {
      message("Heatmap skipped (", file_stub, "): need >=2 proteins; found ", length(protein_ids), ".")
      return(invisible(NULL))
    }
    
    Xh <- expr_imp_viz[rownames(expr_imp_viz) %in% protein_ids, , drop = FALSE]
    if (nrow(Xh) < 2) {
      message("Heatmap skipped (", file_stub, "): matched <2 proteins in matrix.")
      return(invisible(NULL))
    }
    
    rownames(Xh) <- label_from_map_unique(rownames(Xh), id_map)
    heat_scaled <- t(scale(t(as.matrix(Xh))))
    
    ann_col <- data.frame(Condition = factor(as.character(condition), levels = levels(condition)))
    rownames(ann_col) <- colnames(heat_scaled)
    
    pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(7, "RdBu")))(100)
    
    pdf_file <- safe_path(heatmap_dir, paste0(date_str, "_Heatmap_", file_stub, ".pdf"))
    png_file <- safe_path(heatmap_dir, paste0(date_str, "_Heatmap_", file_stub, ".png"))
    
    n_rows <- nrow(heat_scaled)
    n_cols <- ncol(heat_scaled)
    
    # # fixed figure size
    # PDF_W <- 14
    # PDF_H <- 10
    # 
    # # adaptive cell sizes
    # cell_w <- max(3, min(8, 160 / n_cols))
    # cell_h <- max(3, min(8, 120 / n_rows))
    # 
    # # layout tuning
    # tree_row <- 20
    # tree_col <- 20
    # fs_row <- if (n_rows <= 50) 9 else if (n_rows <= 100) 8 else 6
    # fs_col <- if (n_cols <= 30) 8 else if (n_cols <= 60) 7 else 5
    
    pdf_file <- safe_path(heatmap_dir, paste0(date_str, "_Heatmap_", file_stub, ".pdf"))
    png_file <- safe_path(heatmap_dir, paste0(date_str, "_Heatmap_", file_stub, ".png"))
    
    PDF_W <- 10
    PDF_H <- 10
    cell_h <- 4
    cell_w <- 3.5
    tree_row <- 30
    tree_col <- 30
    fs_row <- 8
    fs_col <- 8
    
    grDevices::pdf(pdf_file, width = PDF_W, height = PDF_H, useDingbats = FALSE)
    pheatmap::pheatmap(
      heat_scaled,
      show_colnames = TRUE,
      show_rownames = TRUE,
      annotation_col = ann_col,
      annotation_colors = list(Condition = condition_colors),
      annotation_legend = TRUE,
      cluster_cols = TRUE,
      cluster_rows = TRUE,
      breaks = seq(-HEATMAP_ZLIM, HEATMAP_ZLIM, length.out = 101),
      main = title_text,
      color = pal,
      border_color = NA,
      cellheight = cell_h,
      cellwidth = cell_w,
      treeheight_row = tree_row,
      treeheight_col = tree_col,
      fontsize_row = fs_row,
      fontsize_col = fs_col
    )
    grDevices::dev.off()
    
    grDevices::png(filename = png_file, width = PDF_W, height = PDF_H, units = "in", res = 600)
    pheatmap::pheatmap(
      heat_scaled,
      show_colnames = FALSE,
      show_rownames = TRUE,
      annotation_col = ann_col,
      annotation_colors = list(Condition = condition_colors),
      annotation_legend = TRUE,
      cluster_cols = TRUE,
      cluster_rows = TRUE,
      breaks = seq(-HEATMAP_ZLIM, HEATMAP_ZLIM, length.out = 101),
      main = title_text,
      color = pal,
      border_color = NA,
      cellheight = cell_h,
      cellwidth = cell_w,
      treeheight_row = tree_row,
      treeheight_col = tree_col,
      fontsize_row = fs_row,
      fontsize_col = fs_col
    )
    grDevices::dev.off()
    
    invisible(NULL)
  }
  
  heat_ids_all <- sig_hits_fdr %>% dplyr::pull(Protein) %>% unique()
  make_heatmap(
    protein_ids = heat_ids_all,
    title_text = paste0(
      "Significant proteins across any contrast (limma BH-FDR<", alpha_fdr_primary, ")\n",
      "Row z-score; viz-only imputation used for heatmap display only"
    ),
    file_stub = paste0("AllContrasts_FDRlt", gsub("\\.", "", as.character(alpha_fdr_primary)))
  )
  
  major_contrasts <- contrast_names[grepl(paste0("_vs_", control_label, "$"), contrast_names)]
  
  if (length(major_contrasts)) {
    for (cn in major_contrasts) {
      ids_cn <- merged_all_flagged %>%
        dplyr::filter(Contrast == cn, is.finite(limma_FDR), limma_FDR < alpha_fdr_primary) %>%
        dplyr::pull(Protein) %>%
        unique()
      
      make_heatmap(
        protein_ids = ids_cn,
        title_text = paste0(
          "Significant proteins: ", gsub("_vs_", " vs ", cn),
          " (limma BH-FDR<", alpha_fdr_primary, ")\n",
          "Row z-score; viz-only imputation used for heatmap display only"
        ),
        file_stub = paste0("Major_", safe_stub(cn, "contrast"),
                           "_FDRlt", gsub("\\.", "", as.character(alpha_fdr_primary)))
      )
    }
  } else {
    message("No major (Disease vs Control) contrasts found for per-contrast heatmaps.")
  }
  
  # ------------------- 16.2 VOLCANO PLOTS -------------------
  NS_col <- "#9CA3AF"
  
  contrast_color <- function(contrast_name, condition_colors,
                             control_label = "Control",
                             dd_color = "black",
                             fallback = "#6B7280") {
    parts <- strsplit(contrast_name, "_vs_")[[1]]
    if (length(parts) != 2) return(fallback)
    g1 <- parts[1]
    g2 <- parts[2]
    
    if (g2 == control_label && g1 %in% names(condition_colors)) {
      return(unname(condition_colors[g1]))
    }
    if (g1 != control_label && g2 != control_label) {
      return(dd_color)
    }
    fallback
  }
  
  make_volcano_df <- function(df, fdr_thr = 0.05, p_thr = 0.05) {
    df %>%
      dplyr::mutate(
        neglog10P = -log10(pmax(limma_P, 1e-300)),
        SigClass = dplyr::case_when(
          is.finite(limma_FDR) & limma_FDR < fdr_thr ~ paste0("FDR<", fdr_thr),
          is.finite(limma_P)   & limma_P   < p_thr  ~ paste0("Nominal p<", p_thr),
          TRUE ~ "NS"
        ),
        Label = dplyr::case_when(
          !is.na(Gene) & nzchar(Gene) ~ Gene,
          !is.na(Accession) & nzchar(Accession) ~ Accession,
          TRUE ~ Protein
        )
      )
  }
  
  plot_volcano_one <- function(df_contrast, contrast_name,
                               condition_colors,
                               control_label = "Control",
                               fdr_thr = 0.05, p_thr = 0.05,
                               label_mode = "top",
                               label_top = 10,
                               label_fdr = 0.05,
                               label_cap = Inf) {
    
    if (!label_mode %in% c("top", "all_fdr")) {
      stop("label_mode must be 'top' or 'all_fdr'")
    }
    
    parts <- strsplit(contrast_name, "_vs_")[[1]]
    g1 <- if (length(parts) == 2) parts[1] else "Group1"
    g2 <- if (length(parts) == 2) parts[2] else "Group2"
    
    d <- make_volcano_df(df_contrast, fdr_thr = fdr_thr, p_thr = p_thr)
    
    sig_col <- contrast_color(contrast_name, condition_colors, control_label = control_label, dd_color = "black")
    d <- d %>% dplyr::mutate(PlotColor = ifelse(SigClass == "NS", NS_col, sig_col))
    
    lab <- if (label_mode == "all_fdr") {
      d %>%
        dplyr::filter(is.finite(limma_FDR) & limma_FDR < label_fdr) %>%
        dplyr::arrange(limma_FDR, dplyr::desc(abs(limma_logFC)))
    } else {
      d %>%
        dplyr::filter(SigClass != "NS") %>%
        dplyr::arrange(
          dplyr::case_when(grepl("^FDR<", SigClass) ~ 0, TRUE ~ 1),
          limma_FDR, limma_P, dplyr::desc(abs(limma_logFC))
        ) %>%
        dplyr::slice(1:label_top)
    }
    
    if (is.finite(label_cap) && nrow(lab) > label_cap) {
      lab <- lab %>% dplyr::slice(1:label_cap)
    }
    
    p <- ggplot(d, aes(x = limma_logFC, y = neglog10P)) +
      geom_hline(yintercept = -log10(p_thr), linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_point(
        data = d %>% dplyr::filter(SigClass == "NS"),
        aes(color = PlotColor), shape = 16, size = 1.2, alpha = 0.35
      ) +
      geom_point(
        data = d %>% dplyr::filter(grepl("^FDR<", SigClass)),
        aes(color = PlotColor), shape = 16, size = 2, alpha = 0.95
      ) +
      geom_point(
        data = d %>% dplyr::filter(SigClass == paste0("Nominal p<", p_thr)),
        aes(color = PlotColor), shape = 1, size = 2.3, alpha = 0.95, stroke = 1.1
      ) +
      scale_color_identity() +
      ggrepel::geom_text_repel(data = lab, aes(label = Label), size = 3, max.overlaps = Inf) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      labs(
        title = paste0("Volcano: ", g1, " vs ", g2),
        subtitle = paste0(
          "Filled = BH-FDR<", fdr_thr,
          "; Open = nominal p<", p_thr,
          " (FDR≥", fdr_thr, "); dashed line = p=", p_thr
        ),
        x = paste0("limma adjusted log2FC (", g1, " - ", g2, ")"),
        y = "-log10(p-value)"
      )
    
    stub <- safe_stub(contrast_name, "contrast")
    safe_ggsave(
      safe_path(volcano_dir, paste0(date_str, "_Volcano_", stub, ".png")),
      plot = p, width = 6.2, height = 4.5, dpi = 300
    )
    safe_ggsave(
      safe_path(volcano_dir, paste0(date_str, "_Volcano_", stub, ".pdf")),
      plot = p, width = 6.2, height = 4.5
    )
    
    invisible(p)
  }
  
  if (nrow(merged_all_flagged)) {
    for (cn in unique(merged_all_flagged$Contrast)) {
      dfc <- merged_all_flagged %>% dplyr::filter(Contrast == cn)
      plot_volcano_one(
        dfc, cn, condition_colors, control_label,
        fdr_thr = alpha_fdr_primary,
        p_thr = alpha_p_exploratory,
        label_mode = volcano_label_mode,
        label_top = volcano_label_top,
        label_fdr = volcano_label_fdr,
        label_cap = volcano_label_cap
      )
    }
    message("Volcano plots saved in: ", volcano_dir)
  }
  
  # -------------------- 16.3 BOXPLOTS -----------------------
  meta_simple <- meta_sub %>%
    dplyr::select(Sample_ID, Condition) %>%
    dplyr::distinct()
  
  expr_long_all <- expr_mat %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Protein") %>%
    tidyr::pivot_longer(-Protein, names_to = "Sample_ID", values_to = "Expression") %>%
    dplyr::left_join(meta_simple, by = "Sample_ID") %>%
    dplyr::mutate(Condition = factor(as.character(Condition), levels = box_levels))
  
  readr::write_csv(
    expr_long_all,
    safe_path(boxplot_dir, paste0(date_str, "_Boxplot_LongData_NO_IMPUTATION.csv"))
  )
  
  ids_fdr_005 <- merged_all_flagged %>%
    dplyr::filter(is.finite(limma_FDR), limma_FDR < 0.05) %>%
    dplyr::pull(Protein) %>%
    unique()
  
  ids_fdr_01 <- merged_all_flagged %>%
    dplyr::filter(is.finite(limma_FDR), limma_FDR < 0.10) %>%
    dplyr::pull(Protein) %>%
    unique()
  
  ids_p_005 <- merged_all_flagged %>%
    dplyr::filter(is.finite(limma_P), limma_P < 0.05) %>%
    dplyr::pull(Protein) %>%
    unique()
  
  p_to_star <- function(p) {
    if (!is.finite(p)) return("")
    if (p < 1e-4) "****"
    else if (p < 1e-3) "***"
    else if (p < 1e-2) "**"
    else if (p < 5e-2) "*"
    else ""
  }
  
  build_star_pairs_for_protein <- function(prot, levels_present, stat_df, use = c("P", "FDR"),
                                           y_max_data, y_range) {
    use <- match.arg(use)
    
    d <- stat_df %>%
      dplyr::filter(
        Protein == prot,
        Group1 %in% levels_present,
        Group2 %in% levels_present
      )
    
    if (!nrow(d)) return(NULL)
    
    d <- d %>%
      dplyr::mutate(
        pair_a = pmin(Group1, Group2),
        pair_b = pmax(Group1, Group2),
        pair_key = paste(pair_a, pair_b, sep = "__")
      ) %>%
      dplyr::group_by(pair_key) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
    
    pcol <- if (use == "P") "limma_P" else "limma_FDR"
    
    d <- d %>%
      dplyr::mutate(label = vapply(.data[[pcol]], p_to_star, character(1))) %>%
      dplyr::filter(nzchar(label))
    
    if (!nrow(d)) return(NULL)
    
    base <- y_max_data + 0.06 * y_range
    step <- 0.06 * y_range
    
    d %>%
      dplyr::arrange(.data[[pcol]]) %>%
      dplyr::mutate(y.position = base + (dplyr::row_number() - 1) * step) %>%
      dplyr::transmute(
        group1 = Group1,
        group2 = Group2,
        y.position = y.position,
        p.adj.signif = label
      )
  }
  
  plot_box_set <- function(id_set, out_dir, pdf_stub, star_mode = c("P", "FDR")) {
    star_mode <- match.arg(star_mode)
    if (!length(id_set)) return(invisible(NULL))
    
    df_long <- expr_long_all %>% dplyr::filter(Protein %in% id_set)
    pdf_file <- safe_path(out_dir, paste0(date_str, "_", pdf_stub, ".pdf"))
    grDevices::pdf(pdf_file, width = 6.5, height = 4.8)
    
    produced <- 0L
    skipped  <- 0L
    
    for (prot in unique(df_long$Protein)) {
      dfp <- df_long %>%
        dplyr::filter(Protein == prot, is.finite(Expression), !is.na(Condition))
      
      if (all(is.na(dfp$Expression))) {
        skipped <- skipped + 1L
        next
      }
      
      yr_p <- range(dfp$Expression, na.rm = TRUE)
      if (!all(is.finite(yr_p))) {
        skipped <- skipped + 1L
        next
      }
      
      y_range <- diff(yr_p)
      if (!is.finite(y_range) || y_range <= 0) y_range <- 1
      y_max_data <- max(dfp$Expression, na.rm = TRUE)
      
      spairs <- build_star_pairs_for_protein(
        prot = prot,
        levels_present = levels(condition),
        stat_df = merged_all_flagged,
        use = star_mode,
        y_max_data = y_max_data,
        y_range = y_range
      )
      
      y_top <- yr_p[2] + 0.12 * y_range
      if (!is.null(spairs) && nrow(spairs)) {
        y_top <- max(y_top, max(spairs$y.position, na.rm = TRUE) + 0.06 * y_range)
      }
      y_lim_p <- c(yr_p[1] - 0.05 * y_range, y_top)
      
      lab_idx <- match(prot, id_map$Protein)
      title_txt <- if (is.na(lab_idx)) {
        prot
      } else {
        gene <- id_map$Gene[lab_idx]
        acc  <- id_map$Accession[lab_idx]
        if (!is.na(gene) && nzchar(gene) && !is.na(acc) && nzchar(acc)) paste0(gene, " [", acc, "]")
        else if (!is.na(gene) && nzchar(gene)) gene
        else if (!is.na(acc) && nzchar(acc)) acc
        else prot
      }
      
      p <- ggplot(dfp, aes(x = Condition, y = Expression, fill = Condition)) +
        geom_boxplot(width = 0.3, outlier.shape = NA, linewidth = 0.6, alpha = 1) +
        geom_jitter(width = 0.12, height = 0, alpha = 0.60, size = 1.2) +
        ggtitle(paste0(title_txt, " expression")) +
        ylab("Protein expression (Log2)") +
        xlab(NULL) +
        scale_fill_manual(values = condition_colors, drop = FALSE) +
        scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
        coord_cartesian(ylim = y_lim_p) +
        theme_minimal(base_size = 14) +
        theme(
          legend.position = "none",
          axis.text.x = element_text(size = 13, angle = 0, vjust = 0.8),
          axis.text.y = element_text(size = 13),
          axis.title.y = element_text(size = 14),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_line(linewidth = 0.5),
          plot.title = element_text(size = 14, face = "bold"),
          plot.margin = margin(t = 10, r = 10, b = 5, l = 5)
        )
      
      if (!is.null(spairs) && nrow(spairs)) {
        p <- p + ggpubr::stat_pvalue_manual(
          spairs,
          label = "p.adj.signif",
          tip.length = 0.01,
          hide.ns = TRUE
        )
      }
      
      print(p)
      produced <- produced + 1L
      
      fname_stub <- paste0(safe_stub(title_txt, "protein"), "_", safe_stub(prot, "id"))
      robust_save_png(
        p,
        safe_path(out_dir, paste0(date_str, "_", fname_stub, "_boxplot.png"))
      )
    }
    
    grDevices::dev.off()
    
    message(sprintf(
      "Boxplots (%s): requested=%d, produced=%d, skipped=%d. Output: %s",
      star_mode, length(unique(df_long$Protein)), produced, skipped, out_dir
    ))
    
    invisible(NULL)
  }
  
  plot_box_set(ids_p_005,   box_p_sig_dir,    "Boxplots_P_lt005",   star_mode = "P")
  plot_box_set(ids_fdr_005, box_fdr_sig_dir,  "Boxplots_FDR_lt005", star_mode = "FDR")
  plot_box_set(ids_fdr_01,  box_fdr_near_dir, "Boxplots_FDR_lt01",  star_mode = "FDR")
  
  # -------------------- 16.4 FOREST PLOTS -------------------
  for (cn in unique(merged_all_flagged$Contrast)) {
    d <- merged_all_flagged %>% dplyr::filter(Contrast == cn)
    if (use_only_significant_forest) d <- d %>% dplyr::filter(Is_Significant_FDR)
    if (!nrow(d)) next
    
    d <- d %>%
      dplyr::mutate(Label = label_from_map_unique(Protein, id_map)) %>%
      dplyr::arrange(limma_FDR, dplyr::desc(abs(limma_logFC))) %>%
      dplyr::slice(1:min(max_forest_labels, dplyr::n()))
    
    p <- ggplot(d, aes(x = limma_logFC, y = reorder(Label, limma_logFC))) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_errorbarh(aes(xmin = logFC_CI_low, xmax = logFC_CI_high), height = 0.15, alpha = 0.85) +
      geom_point(size = 2) +
      labs(
        title = paste0("Forest: ", cn),
        subtitle = "Points = limma adjusted log2FC; whiskers = 95% CI (empirical-Bayes posterior SE)",
        x = "Adjusted log2 Fold-change (Group1 vs Group2)",
        y = NULL
      ) +
      theme_minimal(base_size = 11)
    
    stub <- safe_stub(cn, "contrast")
    safe_ggsave(
      safe_path(forest_dir, paste0(date_str, "_Forest_", stub, ".pdf")),
      plot = p,
      width = 7.8,
      height = max(3.5, 0.22 * nrow(d) + 1.5)
    )
    safe_ggsave(
      safe_path(forest_dir, paste0(date_str, "_Forest_", stub, ".png")),
      plot = p,
      width = 7.8,
      height = max(3.5, 0.22 * nrow(d) + 1.5),
      dpi = 300
    )
  }
}

# ============ 17) PERMUTATION LABEL RANDOMIZATION ============

if (run_permutation && length(contrast_names)) {
  message("Running permutation label-randomization benchmark...")
  
  perm_summary <- list()
  
  for (nm in contrast_names) {
    parts <- strsplit(nm, "_vs_")[[1]]
    g1 <- parts[1]
    g2 <- parts[2]
    
    sam <- names(condition)[condition %in% c(g1, g2)]
    if (length(sam) < 4) next
    
    Xsub <- expr_mat[, sam, drop = FALSE]
    msub <- meta_sub[match(sam, meta_sub$Sample_ID), , drop = FALSE]
    csub <- droplevels(condition[sam])
    
    obs <- merged_all_flagged %>% dplyr::filter(Contrast == nm)
    obs_nsig <- sum(obs$limma_FDR < perm_alpha, na.rm = TRUE)
    
    nsig_perm <- rep(NA_real_, n_perm)
    
    for (b in seq_len(n_perm)) {
      cperm <- sample(csub, replace = FALSE)
      names(cperm) <- names(csub)
      cperm <- droplevels(cperm)
      
      lp <- run_limma(
        Xsub, cperm, msub, covariates,
        control_label = g2,
        include_dd = FALSE,
        use_covariates = use_covariates
      )
      
      target_nm <- paste0(g1, "_vs_", g2)
      if (target_nm %in% names(lp$results)) {
        nsig_perm[b] <- sum(lp$results[[target_nm]]$limma_FDR < perm_alpha, na.rm = TRUE)
      }
    }
    
    n_finite <- sum(is.finite(nsig_perm))
    p_perm <- (sum(nsig_perm >= obs_nsig, na.rm = TRUE) + 1) / (n_finite + 1)
    
    perm_summary[[nm]] <- tibble::tibble(
      Contrast = nm,
      Observed_Nsig = obs_nsig,
      Perm_Median_Nsig = median(nsig_perm, na.rm = TRUE),
      Perm_05pct_Nsig = as.numeric(stats::quantile(nsig_perm, 0.05, na.rm = TRUE)),
      Perm_95pct_Nsig = as.numeric(stats::quantile(nsig_perm, 0.95, na.rm = TRUE)),
      Perm_Mean_Nsig = mean(nsig_perm, na.rm = TRUE),
      Perm_NA = sum(is.na(nsig_perm)),
      Perm_N_finite = n_finite,
      Perm_P_value_Nsig = p_perm
    )
  }
  
  perm_df <- dplyr::bind_rows(perm_summary)
  
  readr::write_csv(
    perm_df,
    safe_path(perm_dir, paste0(date_str, "_Permutation_label_randomization_summary.csv"))
  )
  
  wbperm <- openxlsx::createWorkbook()
  wb_add_autosized(wbperm, "Permutation_label_randomization", perm_df)
  openxlsx::saveWorkbook(
    wbperm,
    safe_path(perm_dir, paste0(date_str, "_Permutation_label_randomization.xlsx")),
    overwrite = TRUE
  )
}

# ===================== 18) QC & AUDIT ========================

if (run_QC) {
  wbq <- openxlsx::createWorkbook()
  
  samp_df <- tibble::tibble(
    Sample    = colnames(expr_mat),
    Condition = as.character(condition),
    Median    = apply(expr_mat, 2, median, na.rm = TRUE),
    MAD       = apply(expr_mat, 2, mad, na.rm = TRUE),
    PctNA     = apply(expr_mat, 2, function(x) mean(is.na(x)) * 100)
  )
  wb_add_autosized(wbq, "Samples", samp_df)
  
  miss_df <- tibble::tibble(Protein = rownames(expr_mat))
  for (g in levels(condition)) {
    miss_df[[paste0("PctNA_", g)]] <- apply(
      expr_mat[, condition == g, drop = FALSE],
      1,
      function(x) mean(is.na(x)) * 100
    )
  }
  wb_add_autosized(wbq, "Missingness", miss_df)
  
  gm_summ <- id_map %>%
    dplyr::mutate(Mapped = !is.na(Gene)) %>%
    dplyr::count(Mapped, name = "N")
  wb_add_autosized(wbq, "GeneMapSummary", gm_summ)
  
  cov_audit <- tibble::tibble(
    Use_Covariates = use_covariates,
    Covariates = if (use_covariates) paste(covariates, collapse = ", ") else "None",
    Drop_Missing_Covariates = drop_missing_covariates,
    N_total_samples_after_alignment = nrow(meta),
    N_samples_used_for_modeling = ncol(expr_mat),
    N_dropped_missing_covariates = length(dropped_samples),
    Dropped_Sample_IDs = paste(dropped_samples, collapse = ", ")
  )
  wb_add_autosized(wbq, "CovariateCompleteCase", cov_audit)
  
  if (nrow(jaccard)) wb_add_autosized(wbq, "Presence_Jaccard", jaccard)
  
  wb_add_autosized(
    wbq,
    "ZerosNaNs",
    data.frame(
      Zeros_Converted = zeros_to_NA_count,
      NaNs_Converted = nan_to_NA_count,
      Zero_is_missing = zero_is_missing
    )
  )
  
  openxlsx::saveWorkbook(
    wbq,
    safe_path(qc_dir, paste0(date_str, "_QC.xlsx")),
    overwrite = TRUE
  )
}

# ============== 19) PARAMETERS + SESSION INFO ================

params <- list(
  input_path = input_path,
  expression_sheet = expression_sheet,
  metadata_sheet = metadata_sheet,
  gene_map_fasta = gene_map_fasta,
  output_dir = output_dir,
  already_log2 = already_log2,
  already_normalised = already_normalised,
  zero_is_missing = zero_is_missing,
  min_presence_frac = min_presence_frac,
  presence_thresholds = presence_thresholds,
  require_presence_by_group = require_presence_by_group,
  presence_mode = presence_mode,
  control_label = control_label,
  include_disease_vs_disease = include_disease_vs_disease,
  paired_mode = paired_mode,
  use_covariates = use_covariates,
  drop_missing_covariates = drop_missing_covariates,
  covariates = covariates,
  alpha_fdr_primary = alpha_fdr_primary,
  alpha_p_exploratory = alpha_p_exploratory,
  add_global_FDR = add_global_FDR,
  run_QC = run_QC,
  run_LOSO = run_LOSO,
  make_figures = make_figures,
  HEATMAP_ZLIM = HEATMAP_ZLIM,
  run_permutation = run_permutation,
  n_perm = n_perm,
  perm_alpha = perm_alpha,
  volcano_label_mode = volcano_label_mode,
  volcano_label_top = volcano_label_top,
  volcano_label_fdr = volcano_label_fdr,
  volcano_label_cap = volcano_label_cap,
  use_only_significant_forest = use_only_significant_forest,
  max_forest_labels = max_forest_labels,
  box_order = box_order,
  seed = 123
)

jsonlite::write_json(
  params,
  safe_path(base_dir, paste0(date_str, "_parameters.json")),
  pretty = TRUE,
  auto_unbox = TRUE
)

writeLines(
  capture.output(sessionInfo()),
  con = safe_path(base_dir, paste0(date_str, "_sessionInfo.txt"))
)

# ====================== 20) SUMMARY ==========================

by_group_n <- table(condition)

sig_n_rows   <- nrow(sig_hits_fdr)
sig_n_unique <- dplyr::n_distinct(sig_hits_fdr$Protein)

box_fdr_005_n <- if (exists("ids_fdr_005")) length(ids_fdr_005) else 0L
box_fdr_01_n  <- if (exists("ids_fdr_01"))  length(ids_fdr_01)  else 0L
box_p_005_n   <- if (exists("ids_p_005"))   length(ids_p_005)   else 0L

summary_lines <- c(
  paste0("Initial proteins (pre-filter, modeled samples): ", n_start),
  paste0("After presence filter (>= ", min_presence_frac * 100, "%): ", n_after_presence),
  paste0("After variance filter: ", n_after_variance),
  paste0("Final modeled matrix: ", nrow(expr_mat), " proteins x ", ncol(expr_mat), " samples"),
  paste0("Groups (n): ", paste(names(by_group_n), as.integer(by_group_n), sep = "=", collapse = ", ")),
  paste0("Covariates used: ", if (use_covariates) paste(covariates, collapse = ", ") else "None"),
  paste0("Dropped samples (missing covariates): ", length(dropped_samples)),
  paste0("Contrasts: ", paste(contrast_names, collapse = ", ")),
  "Primary inference: limma pairwise contrasts with BH-FDR per contrast",
  "Backup omnibus: limma moderated F-test of Condition terms, BH-FDR",
  if (add_global_FDR) {
    "Secondary: Global BH-FDR across all Protein x Contrast tests included as Global_BH_FDR"
  } else {
    "Global BH-FDR: OFF"
  },
  paste0("Significant hits (Protein x Contrast, BH-FDR<", alpha_fdr_primary, "): ", sig_n_rows),
  paste0("Unique significant proteins (any contrast): ", sig_n_unique),
  paste0("Boxplot proteins (FDR<0.05): ", box_fdr_005_n),
  paste0("Boxplot proteins (FDR<0.1): ", box_fdr_01_n),
  paste0("Boxplot proteins (P<0.05): ", box_p_005_n),
  paste0("Normalization: ", if (already_normalised) "none (already normalized)" else "sample median-centering"),
  "Imputation: heatmap-only (group median) -- NOT used for modeling or boxplots",
  paste0("Zero handling: ", if (zero_is_missing) "zeros converted to NA" else "zeros retained"),
  paste0("Zeros converted to NA: ", zeros_to_NA_count, " | NaNs converted: ", nan_to_NA_count),
  paste0("Non-finite values standardised to NA (Inf/NaN): ", nonfinite_converted),
  paste0("Log2 transformation applied: ", log2_applied),
  paste0("Permutation sanity check: ", if (run_permutation) paste0("ON (n=", n_perm, ")") else "OFF"),
  "Permutation mode: label-randomization benchmark (not covariate-preserving permutation inference)",
  paste0("FASTA mapping: ", sum(!is.na(id_map$Gene)), " mapped / ", nrow(id_map), " total.")
)

writeLines(
  summary_lines,
  con = safe_path(base_dir, paste0(date_str, "_summary.txt"))
)

cat(
  "\n=== Summary ===\n",
  "Significant hits (Protein x Contrast, BH-FDR<", alpha_fdr_primary, "): ", sig_n_rows, "\n",
  "Unique significant proteins (any contrast): ", sig_n_unique, "\n",
  "Boxplot proteins (FDR<0.05): ", box_fdr_005_n, "\n",
  "Boxplot proteins (FDR<0.1): ", box_fdr_01_n, "\n",
  "Boxplot proteins (P<0.05): ", box_p_005_n, "\n",
  "Contrasts: ", paste(contrast_names, collapse = ", "), "\n",
  "Plots in: ", vis_dir, "\n",
  "Excel workbook: ", basename(out_xlsx), "\n"
)

message("Done. Outputs -> ", normalizePath(base_dir, winslash = "/", mustWork = FALSE))