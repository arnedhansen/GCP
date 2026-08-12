# GCP Effect Sizes Table
# Literature effect sizes used for GCP power analysis (contrast on gamma power and frequency).
#
# Key output (data/power_analysis/):
#   GCP_effect_sizes_table.docx
#
# Run: Rscript GCP_effect_sizes_table.R

ensure_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
}

ensure_packages(c("flextable", "officer", "readxl"))
suppressPackageStartupMessages({
  library(flextable)
  library(officer)
  library(readxl)
})

xlsx_path <- "/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/power_analysis/GCP_reported_effect_sizes_table.xlsx"
out_dir <- "/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/power_analysis"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

TABLE_FONT <- "Arial"
TABLE_FONT_SIZE <- 10
EM_DASH <- "\u2014"
DAGGER <- "\u2020"

# %% Excel parsing -------------------------------------------------------------

read_literature_sheet <- function(path) {
  raw <- readxl::read_excel(path, sheet = 1, col_names = FALSE, .name_repair = "minimal")
  raw <- as.data.frame(raw, stringsAsFactors = FALSE)
  names(raw) <- paste0("V", seq_len(ncol(raw)))

  section_starts <- which(grepl("^CONTRAST on Gamma", raw$V1, ignore.case = TRUE))
  if (length(section_starts) < 2L) {
    stop("Expected at least two gamma sections in ", path)
  }

  parse_section <- function(start_row, end_row) {
    rows <- list()
    for (i in seq.int(start_row + 2L, end_row)) {
      ref <- raw$V1[i]
      if (is.na(ref) || !nzchar(trimws(ref))) next
      rows[[length(rows) + 1L]] <- data.frame(
        reference_raw = ref,
        n_raw = as.character(raw$V2[i]),
        conditions_raw = as.character(raw$V3[i]),
        stimuli_raw = as.character(raw$V4[i]),
        sf_raw = as.character(raw$V5[i]),
        effect_size_raw = raw$V7[i],
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  }

  list(
    gamma_power = parse_section(section_starts[1L], section_starts[2L] - 1L),
    gamma_frequency = parse_section(section_starts[2L], nrow(raw))
  )
}

match_reference <- function(df, pattern) {
  hit <- grepl(pattern, df$reference_raw, ignore.case = TRUE)
  if (sum(hit) != 1L) {
    stop("Expected exactly one match for '", pattern, "', found ", sum(hit))
  }
  df[hit, , drop = FALSE]
}

select_references <- function(df, patterns) {
  do.call(rbind, lapply(patterns, function(p) match_reference(df, p)))
}

# %% Formatting helpers --------------------------------------------------------

fmt_reference <- function(x) {
  x <- trimws(x)
  x <- sub("^(Lowet) E\\. et al", "\\1 et al.", x, ignore.case = TRUE)
  x <- sub("^(.+ & .+?) (\\d{4})$", "\\1 (\\2)", x)
  x <- sub(" et al\\. (\\d{4})$", " et al. (\\1)", x)
  x <- sub(" et al (\\d{4})$", " et al. (\\1)", x)
  x <- sub("^(Van Pelt|van Pelt)", "Van Pelt", x)
  x
}

fmt_n <- function(n_raw) {
  n_raw <- ifelse(is.na(n_raw), "", n_raw)
  nhp_only <- grepl("monkey|primate|non-human|rhesus", n_raw, ignore.case = TRUE) &
    !grepl("participants", n_raw, ignore.case = TRUE)
  num <- suppressWarnings(as.numeric(sub("^([0-9]+).*", "\\1", n_raw)))
  out <- ifelse(is.finite(num), as.character(num), EM_DASH)
  ifelse(nhp_only, paste0(out, DAGGER), out)
}

fmt_conditions <- function(conditions_raw, reference_raw) {
  ref <- tolower(reference_raw)
  if (grepl("^hall", ref)) {
    return("Six levels, random order")
  }
  if (grepl("henrie", ref)) {
    before_paren <- sub("\\s*\\(.*", "", conditions_raw)
    nums <- regmatches(before_paren, gregexpr("[0-9]+(?:\\.[0-9]+)?", before_paren))[[1]]
    return(paste(nums, collapse = ", "))
  }
  if (grepl("^orekhova", ref)) {
    return("50, 100")
  }
  nums <- regmatches(conditions_raw, gregexpr("[0-9]+(?:\\.[0-9]+)?", conditions_raw))[[1]]
  if (length(nums) == 0L) {
    return(EM_DASH)
  }
  paste(nums, collapse = ", ")
}

fmt_stimuli <- function(stimuli_raw) {
  x <- tolower(trimws(stimuli_raw))
  mapping <- c(
    "static square-wave grating" = "Static square-wave",
    "static sine-wave grating" = "Static sine-wave",
    "drifting sinusoidal gratings" = "Drifting sinusoidals",
    "drifting sine wave gratings" = "Drifting sine-waves",
    "dynamic concentric square-wave grating" = "Concentric dynamic square-wave",
    "inward-moving, concentric sine-wave grating" = "Concentric dynamic sine-wave",
    "dynamic concentric sine-wave gratings" = "Concentric dynamic sine-wave",
    "sinusoidal luminance gratings" = "Sinusoidal luminances"
  )
  hit <- names(mapping)[vapply(names(mapping), function(k) grepl(k, x, fixed = TRUE), logical(1))]
  if (length(hit) > 0L) {
    return(unname(mapping[[hit[1L]]]))
  }
  if (is.na(stimuli_raw) || !nzchar(trimws(stimuli_raw))) {
    return(EM_DASH)
  }
  stimuli_raw
}

fmt_spatial_frequency <- function(sf_raw, reference_raw) {
  ref <- tolower(reference_raw)
  if (grepl("^orekhova", ref)) {
    return("1.66")
  }
  if (is.na(sf_raw) || grepl("not reported|^-$", sf_raw, ignore.case = TRUE)) {
    return(EM_DASH)
  }
  num <- regmatches(sf_raw, regexpr("[0-9]+(?:\\.[0-9]+)?", sf_raw))
  if (length(num) == 0L) {
    return(EM_DASH)
  }
  num[[1L]]
}

fmt_effect_size <- function(x, reference_raw = "") {
  ref <- tolower(reference_raw)
  if (grepl("murty", ref) && is.character(x) && grepl("/", x)) {
    x <- sub("/.*$", "", x)
  }
  if (is.na(x) || (is.character(x) && (!nzchar(trimws(x)) || trimws(x) %in% c("-", EM_DASH)))) {
    return(EM_DASH)
  }
  if (is.numeric(x)) {
    return(format(x, trim = TRUE, scientific = FALSE, digits = 3))
  }
  x <- trimws(as.character(x))
  if (x %in% c("-", EM_DASH)) {
    return(EM_DASH)
  }
  suppressWarnings(num <- as.numeric(x))
  if (is.finite(num)) {
    return(format(num, trim = TRUE, scientific = FALSE, digits = 3))
  }
  x
}

prepare_contrast_table <- function(df) {
  data.frame(
    Reference = fmt_reference(df$reference_raw),
    N = fmt_n(df$n_raw),
    `Conditions (Michelson Contrast [%])` = mapply(
      fmt_conditions, df$conditions_raw, df$reference_raw, USE.NAMES = FALSE
    ),
    `Grating Stimuli` = vapply(df$stimuli_raw, fmt_stimuli, character(1)),
    `Spatial Frequency (cpd)` = mapply(fmt_spatial_frequency, df$sf_raw, df$reference_raw, USE.NAMES = FALSE),
    `Effect Size` = mapply(fmt_effect_size, df$effect_size_raw, df$reference_raw, USE.NAMES = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# %% Table formatting -----------------------------------------------------------

style_effect_size_ft <- function(ft, col_widths = NULL) {
  ft <- theme_booktabs(ft)
  rule_border <- fp_border(color = "#666666", width = 0.75)
  ft <- hline_top(ft, border = rule_border, part = "header")
  ft <- hline_bottom(ft, border = rule_border, part = "header")
  ft <- hline_bottom(ft, border = rule_border, part = "body")
  ft <- font(ft, fontname = TABLE_FONT, part = "all")
  ft <- fontsize(ft, size = TABLE_FONT_SIZE, part = "all")
  ft <- bold(ft, part = "header")
  ft <- padding(ft, padding = 3, part = "all")
  ft <- align(ft, align = "left", j = 1, part = "all")
  if (ncol(ft$body$dataset) > 1) {
    ft <- align(ft, align = "center", j = 2:ncol(ft$body$dataset), part = "all")
  }
  ft <- valign(ft, valign = "center", part = "all")
  ft <- set_table_properties(ft, layout = "fixed")
  if (!is.null(col_widths) && length(col_widths) == ncol(ft$body$dataset)) {
    ft <- width(ft, j = seq_along(col_widths), width = col_widths)
  } else {
    ft <- autofit(ft)
  }
  ft
}

make_contrast_flextable <- function(df) {
  body <- prepare_contrast_table(df)
  ft <- flextable(body)
  style_effect_size_ft(
    ft,
    col_widths = c(1.8, 0.35, 1.35, 1.35, 0.75, 0.55)
  )
}

add_section <- function(doc, title, ft) {
  doc <- body_add_fpar(doc, fpar(
    ftext(title, prop = fp_text(bold = TRUE, font.family = TABLE_FONT, font.size = TABLE_FONT_SIZE))
  ))
  doc <- body_add_flextable(doc, ft)
  body_add_par(doc, "")
}

# %% Export --------------------------------------------------------------------

literature <- read_literature_sheet(xlsx_path)

gamma_power <- select_references(
  literature$gamma_power,
  c(
    "Bartoli", "Gebodh", "Hall", "Henrie", "Jia", "Koch",
    "Murty", "Orekhova", "Schadow", "Van Pelt et al 2018"
  )
)

gamma_frequency <- select_references(
  literature$gamma_frequency,
  c(
    "Bartoli", "Hadjipapas", "Jia", "Lowet", "Orekhova",
    "Ray", "Roberts", "van Pelt et al 2018"
  )
)

doc <- read_docx()
doc <- add_section(
  doc,
  "Studies reporting on the effect of contrast on gamma power",
  make_contrast_flextable(gamma_power)
)
doc <- add_section(
  doc,
  "Studies reporting on the effect of contrast on gamma frequency",
  make_contrast_flextable(gamma_frequency)
)

docx_path <- file.path(out_dir, "GCP_effect_sizes_table.docx")
print(doc, target = docx_path)
message("Saved: ", docx_path)
