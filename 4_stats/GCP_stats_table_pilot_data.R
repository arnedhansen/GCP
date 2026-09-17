#!/usr/bin/env Rscript
# Pilot descriptive statistics table (Gaze + EEG).
# Loads GCP_merged_data.csv (Include == 1), computes Mdn (MAD) by contrast,
# and writes CSV + Word table with AOC-matched flextable aesthetics.
# Hierarchy: section -> measure (full stimulus window 0-2000 ms only).
# MAD uses stats::mad() default constant (1.4826).

suppressPackageStartupMessages({
  library(officer)
  library(flextable)
})

## Paths -------------------------------------------------------------------

features_csv <- "/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/GCP_merged_data.csv"
stats_root <- "/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/stats"
tables_root <- file.path(stats_root, "tables")
stats_pilot_dir <- file.path(stats_root, "pilot")
tables_pilot_dir <- file.path(tables_root, "pilot")

for (d in c(stats_root, tables_root, stats_pilot_dir, tables_pilot_dir)) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}

## Aesthetics (matched to AOC summary tables) ------------------------------

TABLE_WIDTH_IN <- 6.5
FONT_SIZE <- 6
RULE_WIDTH <- 0.75
LINE_SPACING <- 0.75
ROW_HEIGHT_IN <- 0.14
PAD_MEASURE_IN <- 12

## Measure specs -----------------------------------------------------------

cond_levels <- c(1L, 2L, 3L, 4L)
cond_labels <- c("25%", "50%", "75%", "100%")
cond_map <- setNames(cond_labels, cond_levels)

measure_specs <- data.frame(
  section = c(
    "Gaze Measures", "Gaze Measures", "Gaze Measures", "Gaze Measures",
    "EEG Measures", "EEG Measures"
  ),
  column = c(
    "PupilSize_bl", "MSRate_bl", "BCEA_bl", "Vel2D_bl",
    "Power", "Frequency"
  ),
  label = c(
    "Pupil Size [%]",
    "Microsaccade Rate [%]",
    "BCEA [%]",
    "Eye Velocity [%]",
    "Gamma Peak Power [dB]",
    "Gamma Peak Frequency [Hz]"
  ),
  digits = c(2L, 2L, 2L, 2L, 2L, 1L),
  stringsAsFactors = FALSE
)

## Helpers -----------------------------------------------------------------

empty_contrast_row <- function(measure_label) {
  data.frame(
    Measure = measure_label,
    `25%` = "",
    `50%` = "",
    `75%` = "",
    `100%` = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

fmt_mdn_mad <- function(mdn_val, mad_val, digits = 2L) {
  if (!is.finite(mdn_val) || !is.finite(mad_val)) {
    return("")
  }
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mdn_val, mad_val)
}

doc_add_title <- function(doc, text, size = FONT_SIZE) {
  body_add_fpar(
    doc,
    fpar(ftext(text, fp_text(bold = TRUE, font.size = size)))
  )
}

table_section_spacer_fpar <- function(font_size = FONT_SIZE) {
  fpar(
    ftext("\u00a0", fp_text(font.size = font_size)),
    fp_p = fp_par(line_spacing = 1, padding.top = 0, padding.bottom = 0)
  )
}

style_summary_ft <- function(ft, col_widths) {
  ft <- theme_booktabs(ft)
  rule_border <- fp_border(color = "#666666", width = RULE_WIDTH)
  ft <- hline_top(ft, border = rule_border, part = "header")
  ft <- hline_bottom(ft, border = rule_border, part = "header")
  ft <- hline_bottom(ft, border = rule_border, part = "body")
  ft <- fontsize(ft, size = FONT_SIZE, part = "body")
  ft <- bold(ft, part = "header")
  ft <- padding(ft, padding = 1, part = "all")
  ft <- align(ft, align = "left", j = 1, part = "all")
  if (ncol(ft$body$dataset) > 1) {
    ft <- align(ft, align = "center", j = 2:ncol(ft$body$dataset), part = "all")
  }
  ft <- line_spacing(ft, space = LINE_SPACING, part = "all")
  ft <- hrule(ft, rule = "exact", part = "all")
  ft <- height_all(ft, height = ROW_HEIGHT_IN, part = "all")
  col_widths <- col_widths * (TABLE_WIDTH_IN / sum(col_widths))
  ft <- width(ft, j = seq_along(col_widths), width = col_widths)
  ft <- set_table_properties(ft, layout = "fixed", align = "left")
  ft <- align(ft, align = "center", part = "header")
  align(ft, j = 1, align = "left", part = "header")
}

add_contrast_header_rule <- function(ft, title_row = 1L, j_contrast) {
  bdr <- fp_border(color = "#666666", width = RULE_WIDTH * 0.85)
  ft <- valign(ft, i = title_row, valign = "bottom", part = "header")
  ft <- hline(ft, i = title_row, j = j_contrast, border = bdr, part = "header")
  fontsize(ft, size = FONT_SIZE, part = "header")
}

## Load data ---------------------------------------------------------------

if (!file.exists(features_csv)) {
  stop("Missing merged features CSV: ", features_csv)
}

raw <- utils::read.csv(features_csv, stringsAsFactors = FALSE, check.names = FALSE)
if (!("Include" %in% names(raw))) {
  stop("Column 'Include' not found in ", features_csv)
}

dat <- raw[raw$Include == 1, , drop = FALSE]
if (nrow(dat) == 0) {
  stop("No rows with Include == 1 in ", features_csv)
}

n_subj <- length(unique(dat$ID))
message(sprintf("[GCP PILOT TABLE] Include == 1: N = %d subjects, %d rows", n_subj, nrow(dat)))

needed_cols <- unique(measure_specs$column)
missing_cols <- setdiff(needed_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing measure columns: ", paste(missing_cols, collapse = ", "))
}

## Compute descriptives ----------------------------------------------------

tidy_rows <- list()
wide_rows <- list()
row_i <- 0L

for (sec in unique(measure_specs$section)) {
  row_i <- row_i + 1L
  wide_rows[[row_i]] <- empty_contrast_row(sec)
  attr(wide_rows[[row_i]], "row_type") <- "section"

  specs_sec <- measure_specs[measure_specs$section == sec, , drop = FALSE]
  for (r in seq_len(nrow(specs_sec))) {
    col_name <- specs_sec$column[r]
    lab <- specs_sec$label[r]
    digits <- specs_sec$digits[r]

    formatted <- character(length(cond_levels))
    names(formatted) <- cond_labels

    for (ci in seq_along(cond_levels)) {
      cond <- cond_levels[ci]
      x <- as.numeric(dat[[col_name]][dat$Condition == cond])
      x <- x[is.finite(x)]
      mdn_val <- if (length(x) > 0) stats::median(x) else NA_real_
      mad_val <- if (length(x) > 0) stats::mad(x, constant = 1.4826) else NA_real_
      formatted[ci] <- fmt_mdn_mad(mdn_val, mad_val, digits = digits)
      tidy_rows[[length(tidy_rows) + 1L]] <- data.frame(
        Section = sec,
        Measure = lab,
        Window = "Full Window (0-2000 ms)",
        Column = col_name,
        Condition = cond_map[[as.character(cond)]],
        ConditionCode = cond,
        N = length(x),
        Median = mdn_val,
        MAD = mad_val,
        Formatted = formatted[ci],
        stringsAsFactors = FALSE
      )
    }

    row_i <- row_i + 1L
    wide_rows[[row_i]] <- data.frame(
      Measure = lab,
      `25%` = formatted[["25%"]],
      `50%` = formatted[["50%"]],
      `75%` = formatted[["75%"]],
      `100%` = formatted[["100%"]],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    attr(wide_rows[[row_i]], "row_type") <- "measure"
  }
}

tidy_tbl <- do.call(rbind, tidy_rows)
wide_tbl <- do.call(rbind, wide_rows)
row_types <- vapply(wide_rows, function(x) attr(x, "row_type"), character(1))
section_idx <- which(row_types == "section")
measure_idx <- which(row_types == "measure")

## Write CSV ---------------------------------------------------------------

tidy_csv <- file.path(stats_pilot_dir, "GCP_pilot_descriptives.csv")
wide_csv <- file.path(tables_pilot_dir, "GCP_pilot_descriptives_table.csv")
utils::write.csv(tidy_tbl, tidy_csv, row.names = FALSE)
utils::write.csv(wide_tbl, wide_csv, row.names = FALSE)
message("[GCP PILOT TABLE] Wrote stats CSV -> ", tidy_csv)
message("[GCP PILOT TABLE] Wrote table CSV -> ", wide_csv)

## Build flextable ---------------------------------------------------------

ft <- flextable(wide_tbl)
ft <- delete_part(ft, part = "header")
ft <- add_header_row(
  ft,
  values = c("", "Mdn (MAD)", "Mdn (MAD)", "Mdn (MAD)", "Mdn (MAD)")
)
ft <- add_header_row(
  ft,
  values = c("", "25%", "50%", "75%", "100%"),
  top = TRUE
)
ft <- style_summary_ft(ft, col_widths = c(2.8, 0.925, 0.925, 0.925, 0.925))
ft <- add_contrast_header_rule(ft, title_row = 1L, j_contrast = 2:5)

if (length(section_idx) > 0) {
  ft <- bold(ft, i = section_idx, j = 1, part = "body")
}
if (length(measure_idx) > 0) {
  ft <- bold(ft, i = measure_idx, j = 1, part = "body")
  ft <- padding(ft, i = measure_idx, j = 1, padding.left = PAD_MEASURE_IN, part = "body")
}

## Write Word --------------------------------------------------------------

docx_path <- file.path(tables_pilot_dir, "GCP_pilot_descriptives_table.docx")
doc <- read_docx()
doc <- doc_add_title(doc, "Table 1: Pilot Data Descriptive Statistics", size = FONT_SIZE)
doc <- body_add_flextable(doc, value = ft, align = "left")
doc <- body_add_fpar(doc, table_section_spacer_fpar(FONT_SIZE))
doc <- body_add_fpar(
  doc,
  fpar(
    ftext("Note. ", fp_text(bold = TRUE, font.size = FONT_SIZE)),
    ftext(
      paste(
        "Values are median (Mdn) and median absolute deviation (MAD) across participants.",
        "All measures are from the full stimulus window (0-2000 ms)."
      ),
      fp_text(font.size = FONT_SIZE)
    )
  )
)

print(doc, target = docx_path)
message("[GCP PILOT TABLE] Wrote Word table -> ", docx_path)
message("[GCP PILOT TABLE] Done.")
