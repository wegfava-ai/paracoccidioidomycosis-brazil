# ============================================================
# Script: 24d_create_combined_oni_pcm_burden_figure.R
# Project: Ecoepidemiology of acute/subacute PCM in Brazil
# Purpose: Create a single publication-ready figure combining
#          annual mean ONI and annual PCM burden.
# Language: English
# ============================================================

rm(list = ls())
gc()

# ------------------------------------------------------------
# 0. Performance configuration
# ------------------------------------------------------------

n_cores_requested <- 10

if (requireNamespace("parallel", quietly = TRUE) &&
    requireNamespace("data.table", quietly = TRUE)) {
  available_cores <- parallel::detectCores(logical = TRUE)
  data.table::setDTthreads(threads = min(n_cores_requested, available_cores))
}

# ------------------------------------------------------------
# 1. Package management
# ------------------------------------------------------------

required_packages <- c(
  "data.table",
  "ggplot2",
  "patchwork",
  "scales"
)

install_if_missing <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing_pkgs) > 0) {
    install.packages(missing_pkgs, dependencies = TRUE)
  }
}

install_if_missing(required_packages)

library(data.table)
library(ggplot2)
library(patchwork)
library(scales)

# ------------------------------------------------------------
# 2. User-defined paths
# ------------------------------------------------------------

datasus_root <- "D:/DATASUS"
project_root <- "D:/DATASUS/pcm_acute_subacute_brazil"

data_processed_dir <- file.path(project_root, "data_processed")
models_dir <- file.path(data_processed_dir, "models")
enso_analysis_dir <- file.path(models_dir, "enso_pcm_incidence")
enso_analysis_tables_dir <- file.path(enso_analysis_dir, "tables")
enso_analysis_figures_dir <- file.path(enso_analysis_dir, "figures")
enso_analysis_data_dir <- file.path(enso_analysis_dir, "data")

outputs_dir <- file.path(project_root, "outputs")
outputs_figures_dir <- file.path(outputs_dir, "figures")
publication_ready_outputs_dir <- file.path(project_root, "publication_ready_outputs")

dir.create(enso_analysis_figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(outputs_figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(publication_ready_outputs_dir, recursive = TRUE, showWarnings = FALSE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
report_file <- file.path(
  project_root,
  "logs",
  paste0("24d_create_combined_oni_pcm_burden_figure_report_", timestamp, ".txt")
)

dir.create(dirname(report_file), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

normalize_path <- function(path) {
  ifelse(file.exists(path), normalizePath(path, winslash = "/", mustWork = FALSE), path)
}

find_optional_file <- function(candidate_paths) {
  candidate_paths <- unique(candidate_paths[!is.na(candidate_paths)])
  existing <- candidate_paths[file.exists(candidate_paths)]
  if (length(existing) == 0) return(NA_character_)
  existing[1]
}

find_required_file <- function(candidate_paths, description) {
  file_path <- find_optional_file(candidate_paths)
  if (is.na(file_path)) {
    stop(
      "Could not find required file: ", description,
      "\nCandidate paths:\n",
      paste(candidate_paths, collapse = "\n")
    )
  }
  file_path
}

find_column <- function(dt, candidates, required = TRUE, description = "column") {
  matched <- candidates[candidates %in% names(dt)]
  if (length(matched) > 0) return(matched[1])

  if (required) {
    stop(
      "Could not find required ", description,
      ". Available columns are: ",
      paste(names(dt), collapse = ", ")
    )
  }

  NA_character_
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

classify_enso_phase_plot <- function(oni_values) {
  out <- ifelse(
    is.na(oni_values),
    "NEUTRAL",
    ifelse(
      oni_values >= 0.5, "EL_NINO",
      ifelse(oni_values <= -0.5, "LA_NINA", "NEUTRAL")
    )
  )
  factor(out, levels = c("LA_NINA", "NEUTRAL", "EL_NINO"))
}

write_text_report <- function(path, lines) {
  writeLines(lines, con = path, useBytes = TRUE)
}

save_plot_pair <- function(plot_object, stem, width, height, dpi = 400) {
  png_file <- paste0(stem, ".png")
  pdf_file <- paste0(stem, ".pdf")

  ggsave(
    filename = png_file,
    plot = plot_object,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )

  ggsave(
    filename = pdf_file,
    plot = plot_object,
    width = width,
    height = height,
    bg = "white"
  )

  list(png = png_file, pdf = pdf_file)
}

copy_if_exists <- function(from, to_dir) {
  if (file.exists(from)) {
    file.copy(from, file.path(to_dir, basename(from)), overwrite = TRUE)
  }
}

# ------------------------------------------------------------
# 4. Locate and read annual dataset from Script 24 v5
# ------------------------------------------------------------

annual_dataset_file <- find_required_file(
  candidate_paths = c(
    file.path(enso_analysis_data_dir, "enso_pcm_annual_national_dataset_latest.rds"),
    file.path(enso_analysis_tables_dir, "enso_pcm_annual_national_dataset_latest.csv")
  ),
  description = "annual national PCM + ENSO dataset"
)

if (grepl("\\.rds$", annual_dataset_file, ignore.case = TRUE)) {
  annual_dt <- as.data.table(readRDS(annual_dataset_file))
} else {
  annual_dt <- as.data.table(fread(annual_dataset_file, encoding = "UTF-8"))
}

# ------------------------------------------------------------
# 5. Standardize annual dataset columns
# ------------------------------------------------------------

year_col <- find_column(annual_dt, c("year", "ano", "YEAR"), description = "year column")
cases_col <- find_column(annual_dt, c("cases", "n_cases", "pcm_cases"), description = "cases column")
incidence_col <- find_column(
  annual_dt,
  c("incidence_100k", "annual_incidence_100k", "incidence_per_100k"),
  description = "incidence column"
)
oni_col <- find_column(
  annual_dt,
  c("oni_current", "enso_oni_mean", "oni_mean", "annual_oni_mean"),
  description = "annual mean ONI column"
)
phase_col <- find_column(
  annual_dt,
  c("enso_phase_plot", "enso_phase_current", "enso_phase"),
  required = FALSE,
  description = "ENSO phase column"
)

setnames(annual_dt, year_col, "year", skip_absent = TRUE)
setnames(annual_dt, cases_col, "cases", skip_absent = TRUE)
setnames(annual_dt, incidence_col, "incidence_100k", skip_absent = TRUE)
setnames(annual_dt, oni_col, "oni_current", skip_absent = TRUE)
if (!is.na(phase_col)) setnames(annual_dt, phase_col, "enso_phase_raw", skip_absent = TRUE)

annual_dt[
  ,
  `:=`(
    year = as.integer(year),
    cases = as.integer(round(safe_numeric(cases))),
    incidence_100k = safe_numeric(incidence_100k),
    oni_current = safe_numeric(oni_current)
  )
]

annual_dt <- annual_dt[order(year)]

if ("enso_phase_raw" %in% names(annual_dt)) {
  annual_dt[, enso_phase_plot := factor(as.character(enso_phase_raw), levels = c("LA_NINA", "NEUTRAL", "EL_NINO"))]
} else {
  annual_dt[, enso_phase_plot := classify_enso_phase_plot(oni_current)]
}

annual_dt[, enso_phase_label := fifelse(
  enso_phase_plot == "EL_NINO", "El Niño",
  fifelse(enso_phase_plot == "LA_NINA", "La Niña", "Neutral")
)]

annual_dt[, year_label := as.character(year)]

# ------------------------------------------------------------
# 6. Prepare plotting data
# ------------------------------------------------------------

phase_colors <- c(
  "LA_NINA" = "#2F6DB3",
  "NEUTRAL" = "#B9B9B9",
  "EL_NINO" = "#C8102E"
)

max_cases <- max(annual_dt$cases, na.rm = TRUE)
max_incidence <- max(annual_dt$incidence_100k, na.rm = TRUE)

if (!is.finite(max_incidence) || max_incidence <= 0) {
  stop("Could not calculate a positive maximum annual incidence value.")
}

incidence_scale_factor <- max_cases / max_incidence * 0.85

annual_dt[, incidence_scaled := incidence_100k * incidence_scale_factor]

# For cleaner labels, only label all case bars if the number of years is manageable.
label_all_cases <- nrow(annual_dt) <= 30

# ------------------------------------------------------------
# 7. Panel A: Annual mean ONI
# ------------------------------------------------------------

p_oni <- ggplot(annual_dt, aes(x = factor(year), y = oni_current, fill = enso_phase_plot)) +
  geom_col(width = 0.78, color = NA) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey55") +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.45, color = "black") +
  scale_fill_manual(
    values = phase_colors,
    breaks = c("LA_NINA", "NEUTRAL", "EL_NINO"),
    labels = c("La Niña", "Neutral", "El Niño"),
    name = "Annual ENSO phase"
  ) +
  labs(
    title = "A. Annual mean Oceanic Niño Index",
    subtitle = "ENSO phase defined by annual mean ONI: El Niño ≥ 0.5; La Niña ≤ -0.5",
    x = NULL,
    y = "Annual mean ONI"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# ------------------------------------------------------------
# 8. Panel B: Annual PCM cases + incidence
# ------------------------------------------------------------

p_pcm <- ggplot(annual_dt, aes(x = factor(year))) +
  geom_col(
    aes(y = cases, fill = enso_phase_plot),
    width = 0.78,
    alpha = 0.95,
    color = NA
  ) +
  geom_line(
    aes(y = incidence_scaled, group = 1),
    color = "black",
    linewidth = 0.95
  ) +
  geom_point(
    aes(y = incidence_scaled),
    color = "black",
    size = 1.9
  ) +
  {
    if (label_all_cases) {
      geom_text(
        aes(y = cases, label = cases),
        vjust = -0.25,
        size = 3.0,
        color = "grey10"
      )
    }
  } +
  scale_fill_manual(
    values = phase_colors,
    breaks = c("LA_NINA", "NEUTRAL", "EL_NINO"),
    labels = c("La Niña", "Neutral", "El Niño"),
    name = "Annual ENSO phase"
  ) +
  scale_y_continuous(
    name = "Annual PCM cases (n)",
    sec.axis = sec_axis(
      trans = ~ . / incidence_scale_factor,
      name = "Annual incidence per 100,000 inhabitants"
    ),
    expand = expansion(mult = c(0, 0.10))
  ) +
  labs(
    title = "B. Annual national probable acute/subacute PCM burden",
    subtitle = "Bars show yearly case counts; black line shows annual incidence",
    x = "Year",
    y = "Annual PCM cases (n)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# ------------------------------------------------------------
# 9. Combined figure
# ------------------------------------------------------------

combined_title <- "Annual ENSO pattern and annual PCM burden in Brazil"
combined_subtitle <- "Probable acute/subacute PCM, 1998-2024"

combined_plot <- p_oni / p_pcm +
  plot_layout(heights = c(1, 1.15), guides = "collect") +
  plot_annotation(
    title = combined_title,
    subtitle = combined_subtitle,
    caption = paste(
      "Top panel: annual mean ONI categorized as El Niño, Neutral, or La Niña.",
      "Bottom panel: annual probable acute/subacute PCM cases (bars) and annual incidence per 100,000 inhabitants (black line).",
      "Both panels use the same ENSO phase color coding."
    )
  ) &
  theme(
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 13),
    plot.caption = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

# ------------------------------------------------------------
# 10. Save outputs
# ------------------------------------------------------------

figure_stem_main <- file.path(
  enso_analysis_figures_dir,
  paste0("24d_combined_oni_pcm_burden_figure_", timestamp)
)

figure_stem_latest <- file.path(
  enso_analysis_figures_dir,
  "combined_oni_pcm_burden_figure_latest"
)

saved_main <- save_plot_pair(combined_plot, figure_stem_main, width = 16, height = 12, dpi = 420)
saved_latest <- save_plot_pair(combined_plot, figure_stem_latest, width = 16, height = 12, dpi = 420)

# Save the plotting dataset used in the figure
figure_data_file <- file.path(
  enso_analysis_tables_dir,
  paste0("24d_combined_oni_pcm_burden_figure_data_", timestamp, ".csv")
)
figure_data_latest <- file.path(
  enso_analysis_tables_dir,
  "combined_oni_pcm_burden_figure_data_latest.csv"
)

fwrite(annual_dt, figure_data_file)
fwrite(annual_dt, figure_data_latest)

# Copy to user-friendly output folders
for (file_i in c(saved_main$png, saved_main$pdf, saved_latest$png, saved_latest$pdf)) {
  copy_if_exists(file_i, outputs_figures_dir)
  copy_if_exists(file_i, publication_ready_outputs_dir)
}

copy_if_exists(figure_data_file, outputs_figures_dir)

# ------------------------------------------------------------
# 11. Build textual report
# ------------------------------------------------------------

report_lines <- c(
  "============================================================",
  "24d_create_combined_oni_pcm_burden_figure.R",
  "============================================================",
  paste("Run timestamp:", timestamp),
  "",
  "Input file",
  "------------------------------------------------------------",
  paste("Annual dataset:", normalize_path(annual_dataset_file)),
  "",
  "Figure concept",
  "------------------------------------------------------------",
  "One single figure combining annual mean ONI with annual PCM burden.",
  "Top panel: annual mean ONI as bars, colored by ENSO phase.",
  "Bottom panel: annual PCM cases as bars, colored by the same ENSO phase.",
  "Bottom panel also shows annual incidence as a black line with a secondary y-axis.",
  "",
  "Key settings",
  "------------------------------------------------------------",
  paste("Requested CPU threads:", n_cores_requested),
  paste("data.table threads actually used:", data.table::getDTthreads()),
  paste("Years in the figure:", min(annual_dt$year, na.rm = TRUE), "-", max(annual_dt$year, na.rm = TRUE)),
  paste("Total years plotted:", nrow(annual_dt)),
  paste("Total PCM cases across all years:", format(sum(annual_dt$cases, na.rm = TRUE), big.mark = ",")),
  paste("Mean annual cases:", round(mean(annual_dt$cases, na.rm = TRUE), 2)),
  paste("Median annual cases:", round(median(annual_dt$cases, na.rm = TRUE), 2)),
  paste("Mean annual incidence per 100,000:", round(mean(annual_dt$incidence_100k, na.rm = TRUE), 4)),
  paste("Median annual incidence per 100,000:", round(median(annual_dt$incidence_100k, na.rm = TRUE), 4)),
  paste("Incidence scale factor for the secondary axis:", round(incidence_scale_factor, 4)),
  "",
  "Saved outputs",
  "------------------------------------------------------------",
  paste("Main PNG:", normalize_path(saved_main$png)),
  paste("Main PDF:", normalize_path(saved_main$pdf)),
  paste("Latest PNG:", normalize_path(saved_latest$png)),
  paste("Latest PDF:", normalize_path(saved_latest$pdf)),
  paste("Figure data CSV:", normalize_path(figure_data_file)),
  "",
  "Additional copies",
  "------------------------------------------------------------",
  paste("Outputs figures folder:", normalize_path(outputs_figures_dir)),
  paste("Publication-ready outputs folder:", normalize_path(publication_ready_outputs_dir)),
  "",
  "Interpretation note",
  "------------------------------------------------------------",
  "This figure is descriptive and is designed to unify ENSO phase classification,",
  "annual ONI magnitude, annual PCM case counts, and annual incidence in a single visual.",
  "Formal inference about ENSO effects should still rely on the regression models from Script 24.",
  "============================================================"
)

write_text_report(report_file, report_lines)

# ------------------------------------------------------------
# 12. Console summary
# ------------------------------------------------------------

message("Combined figure created successfully.")
message("Input dataset: ", normalize_path(annual_dataset_file))
message("Total PCM cases across all years: ", format(sum(annual_dt$cases, na.rm = TRUE), big.mark = ","))
message("Latest PNG: ", normalize_path(saved_latest$png))
message("Latest PDF: ", normalize_path(saved_latest$pdf))
message("Report: ", normalize_path(report_file))
