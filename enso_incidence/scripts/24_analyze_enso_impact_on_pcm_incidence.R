# ============================================================
# Project: Ecoepidemiology of acute/subacute paracoccidioidomycosis in Brazil
# Script: 24_analyze_enso_impact_on_pcm_incidence.R
# Purpose: Analyze the association between El Niño/ENSO and annual incidence of probable acute/subacute PCM
# Author: Wellington Fava
# Date: 2026-06-26
# Version: v1_enso_pcm_incidence
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

datasus_root <- "D:/DATASUS"
project_root <- file.path(datasus_root, "pcm_acute_subacute_brazil")

logs_dir <- file.path(project_root, "logs")
outputs_dir <- file.path(project_root, "outputs")
tables_dir <- file.path(outputs_dir, "tables")
figures_dir <- file.path(outputs_dir, "figures")

data_processed_dir <- file.path(project_root, "data_processed")
modeling_dir <- file.path(data_processed_dir, "modeling")
final_modeling_dir <- file.path(data_processed_dir, "final_modeling")
enso_processed_dir <- file.path(data_processed_dir, "enso")
models_dir <- file.path(data_processed_dir, "models")

enso_analysis_dir <- file.path(models_dir, "enso_pcm_incidence")
enso_analysis_tables_dir <- file.path(enso_analysis_dir, "tables")
enso_analysis_figures_dir <- file.path(enso_analysis_dir, "figures")
enso_analysis_data_dir <- file.path(enso_analysis_dir, "data")

for (dir_path in c(
  logs_dir,
  outputs_dir,
  tables_dir,
  figures_dir,
  data_processed_dir,
  modeling_dir,
  final_modeling_dir,
  enso_processed_dir,
  models_dir,
  enso_analysis_dir,
  enso_analysis_tables_dir,
  enso_analysis_figures_dir,
  enso_analysis_data_dir
)) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("Created folder: ", dir_path)
  }
}

# ------------------------------------------------------------
# 2. Settings
# ------------------------------------------------------------

primary_period_start <- 1998
primary_period_end <- 2024
age_cutoff_years <- 20

outcome_definition <- "probable acute/subacute paracoccidioidomycosis"
outcome_short <- "probable acute/subacute PCM"

el_nino_threshold <- 0.5
la_nina_threshold <- -0.5
include_linear_year_trend <- TRUE

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# ------------------------------------------------------------
# 3. Packages
# ------------------------------------------------------------

required_packages <- c("data.table", "ggplot2", "MASS")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "The following required packages are not installed: ",
    paste(missing_packages, collapse = ", "),
    "\nPlease install them before running this script:\n",
    "install.packages(c(",
    paste(sprintf("'%s'", missing_packages), collapse = ", "),
    "))"
  )
}

library(data.table)
library(ggplot2)
library(MASS)

# ------------------------------------------------------------
# 4. Helper functions
# ------------------------------------------------------------

safe_character <- function(x) {
  x <- as.character(x)
  x <- iconv(x, from = "", to = "UTF-8", sub = "")
  x <- trimws(x)
  x[x == ""] <- NA_character_
  x
}

safe_numeric <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x <- safe_character(x)
  has_comma <- grepl(",", x)
  x_out <- x
  x_out[has_comma] <- gsub("\\.", "", x_out[has_comma])
  x_out[has_comma] <- gsub(",", ".", x_out[has_comma])
  x_out <- gsub("\\s+", "", x_out)
  x_out <- gsub("[^0-9\\.-]", "", x_out)
  suppressWarnings(as.numeric(x_out))
}

normalize_path <- function(x) {
  x <- safe_character(x)
  x <- gsub("\\\\", "/", x)
  x
}

find_latest_file <- function(directory, pattern) {
  if (!dir.exists(directory)) return(NA_character_)
  files <- list.files(
    directory,
    pattern = pattern,
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(files) == 0) return(NA_character_)
  files[which.max(file.info(files)$mtime)]
}

find_optional_file <- function(candidate_paths) {
  candidate_paths <- candidate_paths[!is.na(candidate_paths)]
  existing_paths <- candidate_paths[file.exists(candidate_paths)]
  if (length(existing_paths) == 0) return(NA_character_)
  existing_paths[1]
}

find_required_file <- function(candidate_paths, file_description) {
  file_path <- find_optional_file(candidate_paths)
  if (is.na(file_path)) {
    stop(
      "Could not find required file: ",
      file_description,
      "\nCandidate paths:\n",
      paste(candidate_paths, collapse = "\n")
    )
  }
  file_path
}

find_column <- function(dt, candidates, required = TRUE, description = "column") {
  dt <- as.data.table(dt)
  name_map <- data.table(
    original_name = names(dt),
    lower_name = tolower(names(dt))
  )

  for (candidate_i in candidates) {
    hit_i <- name_map[lower_name == tolower(candidate_i), original_name]
    if (length(hit_i) > 0) return(hit_i[1])
  }

  for (candidate_i in candidates) {
    hit_i <- name_map[grepl(tolower(candidate_i), lower_name, fixed = TRUE), original_name]
    if (length(hit_i) > 0) return(hit_i[1])
  }

  if (required) {
    stop(
      "Could not find required ",
      description,
      ". Candidate names searched: ",
      paste(candidates, collapse = ", "),
      "\nAvailable columns:\n",
      paste(names(dt), collapse = ", ")
    )
  }

  NA_character_
}

classify_enso_phase <- function(oni_value) {
  oni_value <- suppressWarnings(as.numeric(oni_value))
  out <- rep(NA_character_, length(oni_value))
  out[is.finite(oni_value) & oni_value >= el_nino_threshold] <- "El Nino"
  out[is.finite(oni_value) & oni_value <= la_nina_threshold] <- "La Nina"
  out[is.finite(oni_value) & oni_value > la_nina_threshold & oni_value < el_nino_threshold] <- "Neutral"
  factor(out, levels = c("Neutral", "El Nino", "La Nina"))
}

fmt_num <- function(x, digits = 4) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || is.na(x[1]) || !is.finite(x[1])) return("NA")
  formatC(x[1], format = "f", digits = digits)
}

fmt_int <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || is.na(x[1]) || !is.finite(x[1])) return("NA")
  formatC(round(x[1]), format = "d", big.mark = ",")
}

write_text_file <- function(text, file_path) {
  writeLines(text, file_path, useBytes = TRUE)
}

save_plot_pair <- function(plot_object, file_stem, width = 11, height = 7, dpi = 320) {
  png_file <- file.path(enso_analysis_figures_dir, paste0(file_stem, ".png"))
  pdf_file <- file.path(enso_analysis_figures_dir, paste0(file_stem, ".pdf"))

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

  file.copy(
    from = png_file,
    to = file.path(figures_dir, basename(png_file)),
    overwrite = TRUE
  )

  file.copy(
    from = pdf_file,
    to = file.path(figures_dir, basename(pdf_file)),
    overwrite = TRUE
  )

  data.table(
    figure_stem = file_stem,
    png_file = normalize_path(png_file),
    pdf_file = normalize_path(pdf_file),
    png_file_outputs = normalize_path(file.path(figures_dir, basename(png_file))),
    pdf_file_outputs = normalize_path(file.path(figures_dir, basename(pdf_file)))
  )
}

tidy_model <- function(model_object, model_name, exposure_label, model_family) {
  if (is.null(model_object)) return(data.table())

  coef_matrix <- coef(summary(model_object))
  coef_dt <- as.data.table(coef_matrix, keep.rownames = "term")

  estimate_col <- find_column(coef_dt, c("Estimate"), description = "estimate column")
  se_col <- find_column(coef_dt, c("Std. Error"), description = "standard error column")
  statistic_col <- find_column(coef_dt, c("z value", "t value"), required = FALSE, description = "statistic column")
  p_col <- find_column(coef_dt, c("Pr(>|z|)", "Pr(>|t|)", "p value"), required = FALSE, description = "p-value column")

  coef_dt[, estimate_log := safe_numeric(get(estimate_col))]
  coef_dt[, std_error := safe_numeric(get(se_col))]

  if (!is.na(statistic_col)) {
    coef_dt[, statistic := safe_numeric(get(statistic_col))]
  } else {
    coef_dt[, statistic := NA_real_]
  }

  if (!is.na(p_col)) {
    coef_dt[, p_value := safe_numeric(get(p_col))]
  } else {
    coef_dt[, p_value := NA_real_]
  }

  coef_dt[
    ,
    `:=`(
      conf_low_log = estimate_log - 1.96 * std_error,
      conf_high_log = estimate_log + 1.96 * std_error,
      irr = exp(estimate_log),
      conf_low_95 = exp(estimate_log - 1.96 * std_error),
      conf_high_95 = exp(estimate_log + 1.96 * std_error),
      model_name = model_name,
      exposure_label = exposure_label,
      model_family = model_family
    )
  ]

  coef_dt[
    ,
    .(
      model_name,
      exposure_label,
      model_family,
      term,
      estimate_log,
      std_error,
      statistic,
      p_value,
      irr,
      conf_low_95,
      conf_high_95
    )
  ]
}

fit_quasipoisson <- function(model_data, formula_object, model_name, exposure_label) {
  out <- tryCatch(
    glm(
      formula = formula_object,
      family = quasipoisson(link = "log"),
      data = model_data
    ),
    error = function(e) {
      message("Quasi-Poisson model failed: ", model_name, " | ", conditionMessage(e))
      NULL
    }
  )

  tidy_model(
    model_object = out,
    model_name = model_name,
    exposure_label = exposure_label,
    model_family = "quasi_poisson"
  )
}

fit_negative_binomial <- function(model_data, formula_object, model_name, exposure_label) {
  out <- tryCatch(
    MASS::glm.nb(
      formula = formula_object,
      data = model_data,
      control = glm.control(maxit = 100)
    ),
    error = function(e) {
      message("Negative binomial model failed: ", model_name, " | ", conditionMessage(e))
      NULL
    }
  )

  tidy_model(
    model_object = out,
    model_name = model_name,
    exposure_label = exposure_label,
    model_family = "negative_binomial"
  )
}

# ------------------------------------------------------------
# 5. Locate and read input files
# ------------------------------------------------------------

primary_dataset_file <- find_required_file(
  candidate_paths = c(
    file.path(final_modeling_dir, "pcm_final_municipality_year_modeling_dataset_1998_2024_latest.rds"),
    file.path(modeling_dir, "pcm_final_municipality_year_modeling_dataset_1998_2024_latest.rds"),
    find_latest_file(final_modeling_dir, "pcm_final_municipality_year_modeling_dataset.*\\.rds$"),
    find_latest_file(modeling_dir, "pcm_final_municipality_year_modeling_dataset.*\\.rds$")
  ),
  file_description = "final municipality-year PCM modeling dataset"
)

enso_year_file <- find_required_file(
  candidate_paths = c(
    file.path(enso_processed_dir, "oni_modeling_year_indices_latest.rds"),
    find_latest_file(enso_processed_dir, "oni_modeling_year_indices.*\\.rds$")
  ),
  file_description = "annual ENSO modeling-year indices"
)

primary_dt <- as.data.table(readRDS(primary_dataset_file))
enso_year <- as.data.table(readRDS(enso_year_file))

input_file_summary <- data.table(
  input_name = c("primary_dataset_file", "enso_year_file"),
  file_path = c(primary_dataset_file, enso_year_file)
)

# ------------------------------------------------------------
# 6. Standardize core columns
# ------------------------------------------------------------

year_col <- find_column(primary_dt, c("year", "ano", "event_year", "YEAR"), description = "year column")
municipality_col <- find_column(primary_dt, c("municipality_code_6", "codmun6", "municipality_code", "codigo_municipio"), description = "municipality code column")
outcome_col <- find_column(primary_dt, c("pcm_event_count_main_final", "pcm_event_count_main", "pcm_event_count", "event_count", "cases", "n_cases"), description = "PCM event count column")
population_col <- find_column(primary_dt, c("population_final", "population", "pop", "population_denominator", "populacao"), description = "population denominator column")

enso_year_col <- find_column(enso_year, c("year", "ano", "event_year", "YEAR"), description = "ENSO year column")

setnames(primary_dt, year_col, "year", skip_absent = TRUE)
setnames(primary_dt, municipality_col, "municipality_code_6", skip_absent = TRUE)
setnames(primary_dt, outcome_col, "pcm_cases", skip_absent = TRUE)
setnames(primary_dt, population_col, "population", skip_absent = TRUE)
setnames(enso_year, enso_year_col, "year", skip_absent = TRUE)

primary_dt[
  ,
  `:=`(
    year = as.integer(year),
    municipality_code_6 = safe_character(municipality_code_6),
    pcm_cases = safe_numeric(pcm_cases),
    population = safe_numeric(population)
  )
]

enso_year[, year := as.integer(year)]

analysis_dt <- primary_dt[
  year >= primary_period_start &
    year <= primary_period_end
]

if (nrow(analysis_dt) == 0) {
  stop(
    "No rows were available for the primary analysis period ",
    primary_period_start,
    "-",
    primary_period_end,
    "."
  )
}

# ------------------------------------------------------------
# 7. Build annual analysis dataset
# ------------------------------------------------------------

annual_dt <- analysis_dt[
  ,
  .(
    cases = sum(pcm_cases, na.rm = TRUE),
    population = sum(population, na.rm = TRUE),
    municipalities = uniqueN(municipality_code_6),
    rows = .N,
    event_municipality_years = sum(pcm_cases > 0, na.rm = TRUE)
  ),
  by = year
][order(year)]

annual_dt[
  ,
  incidence_100k := fifelse(
    is.finite(population) & population > 0,
    cases / population * 100000,
    NA_real_
  )
]

oni_current_col <- find_column(
  enso_year,
  c("enso_oni_mean", "oni_mean", "annual_oni_mean", "oni_annual_mean"),
  required = FALSE,
  description = "current annual ONI mean column"
)

oni_lag1_col <- find_column(
  enso_year,
  c("enso_oni_mean_lag1y", "oni_mean_lag1y", "oni_lag1", "oni_mean_lag1"),
  required = FALSE,
  description = "ONI mean lag 1 year column"
)

oni_lag2_col <- find_column(
  enso_year,
  c("enso_oni_mean_lag2y", "oni_mean_lag2y", "oni_lag2", "oni_mean_lag2"),
  required = FALSE,
  description = "ONI mean lag 2 years column"
)

oni_lag3_col <- find_column(
  enso_year,
  c("enso_oni_mean_lag3y", "oni_mean_lag3y", "oni_lag3", "oni_mean_lag3"),
  required = FALSE,
  description = "ONI mean lag 3 years column"
)

if (is.na(oni_current_col)) {
  stop(
    "Could not identify an annual ONI mean column in the ENSO file. Available columns:\n",
    paste(names(enso_year), collapse = ", ")
  )
}

enso_join_cols <- unique(c("year", oni_current_col, oni_lag1_col, oni_lag2_col, oni_lag3_col))
enso_join_cols <- enso_join_cols[!is.na(enso_join_cols)]

enso_join <- enso_year[, ..enso_join_cols]

setnames(enso_join, oni_current_col, "oni_current", skip_absent = TRUE)
if (!is.na(oni_lag1_col)) setnames(enso_join, oni_lag1_col, "oni_lag1y", skip_absent = TRUE)
if (!is.na(oni_lag2_col)) setnames(enso_join, oni_lag2_col, "oni_lag2y", skip_absent = TRUE)
if (!is.na(oni_lag3_col)) setnames(enso_join, oni_lag3_col, "oni_lag3y", skip_absent = TRUE)

annual_dt <- merge(annual_dt, enso_join, by = "year", all.x = TRUE)

for (col_i in c("oni_current", "oni_lag1y", "oni_lag2y", "oni_lag3y")) {
  if (col_i %in% names(annual_dt)) {
    annual_dt[, (col_i) := safe_numeric(get(col_i))]
  }
}

annual_dt[, enso_phase_current := classify_enso_phase(oni_current)]
if ("oni_lag1y" %in% names(annual_dt)) annual_dt[, enso_phase_lag1y := classify_enso_phase(oni_lag1y)]
if ("oni_lag2y" %in% names(annual_dt)) annual_dt[, enso_phase_lag2y := classify_enso_phase(oni_lag2y)]
if ("oni_lag3y" %in% names(annual_dt)) annual_dt[, enso_phase_lag3y := classify_enso_phase(oni_lag3y)]

annual_dt[, year_centered := year - mean(year, na.rm = TRUE)]
annual_dt[, log_population_offset := log(population)]

annual_model_dt <- annual_dt[
  is.finite(cases) &
    is.finite(population) &
    population > 0 &
    is.finite(log_population_offset) &
    !is.na(oni_current)
]

# ------------------------------------------------------------
# 8. Descriptive summaries
# ------------------------------------------------------------

annual_summary <- annual_dt[
  ,
  .(
    years = .N,
    first_year = min(year, na.rm = TRUE),
    last_year = max(year, na.rm = TRUE),
    total_cases = sum(cases, na.rm = TRUE),
    mean_annual_cases = mean(cases, na.rm = TRUE),
    median_annual_cases = median(cases, na.rm = TRUE),
    total_population_person_years = sum(population, na.rm = TRUE),
    mean_annual_incidence_100k = mean(incidence_100k, na.rm = TRUE),
    median_annual_incidence_100k = median(incidence_100k, na.rm = TRUE),
    min_annual_incidence_100k = min(incidence_100k, na.rm = TRUE),
    max_annual_incidence_100k = max(incidence_100k, na.rm = TRUE)
  )
]

phase_summary_current <- annual_dt[
  !is.na(enso_phase_current),
  .(
    years = .N,
    years_list = paste(year, collapse = ", "),
    total_cases = sum(cases, na.rm = TRUE),
    mean_cases = mean(cases, na.rm = TRUE),
    median_cases = median(cases, na.rm = TRUE),
    mean_incidence_100k = mean(incidence_100k, na.rm = TRUE),
    median_incidence_100k = median(incidence_100k, na.rm = TRUE),
    min_incidence_100k = min(incidence_100k, na.rm = TRUE),
    max_incidence_100k = max(incidence_100k, na.rm = TRUE),
    mean_oni = mean(oni_current, na.rm = TRUE)
  ),
  by = enso_phase_current
][order(enso_phase_current)]

phase_summary_lag1 <- data.table()

if ("enso_phase_lag1y" %in% names(annual_dt)) {
  phase_summary_lag1 <- annual_dt[
    !is.na(enso_phase_lag1y),
    .(
      years = .N,
      years_list = paste(year, collapse = ", "),
      total_cases = sum(cases, na.rm = TRUE),
      mean_cases = mean(cases, na.rm = TRUE),
      median_cases = median(cases, na.rm = TRUE),
      mean_incidence_100k = mean(incidence_100k, na.rm = TRUE),
      median_incidence_100k = median(incidence_100k, na.rm = TRUE),
      min_incidence_100k = min(incidence_100k, na.rm = TRUE),
      max_incidence_100k = max(incidence_100k, na.rm = TRUE),
      mean_oni_lag1y = mean(oni_lag1y, na.rm = TRUE)
    ),
    by = enso_phase_lag1y
  ][order(enso_phase_lag1y)]
}

# ------------------------------------------------------------
# 9. Model fitting
# ------------------------------------------------------------

model_terms <- data.table()

build_formula <- function(exposure_variable) {
  if (isTRUE(include_linear_year_trend)) {
    as.formula(paste0("cases ~ ", exposure_variable, " + year_centered + offset(log_population_offset)"))
  } else {
    as.formula(paste0("cases ~ ", exposure_variable, " + offset(log_population_offset)"))
  }
}

exposure_specs <- data.table(
  model_name = c("continuous_ONI_current", "phase_current"),
  exposure_variable = c("oni_current", "enso_phase_current"),
  exposure_label = c("Annual mean ONI, same year", "ENSO phase, same year")
)

if ("oni_lag1y" %in% names(annual_model_dt)) {
  exposure_specs <- rbindlist(
    list(exposure_specs, data.table(model_name = "continuous_ONI_lag1y", exposure_variable = "oni_lag1y", exposure_label = "Annual mean ONI, 1-year lag")),
    fill = TRUE
  )
}

if ("oni_lag2y" %in% names(annual_model_dt)) {
  exposure_specs <- rbindlist(
    list(exposure_specs, data.table(model_name = "continuous_ONI_lag2y", exposure_variable = "oni_lag2y", exposure_label = "Annual mean ONI, 2-year lag")),
    fill = TRUE
  )
}

if ("oni_lag3y" %in% names(annual_model_dt)) {
  exposure_specs <- rbindlist(
    list(exposure_specs, data.table(model_name = "continuous_ONI_lag3y", exposure_variable = "oni_lag3y", exposure_label = "Annual mean ONI, 3-year lag")),
    fill = TRUE
  )
}

if ("enso_phase_lag1y" %in% names(annual_model_dt)) {
  exposure_specs <- rbindlist(
    list(exposure_specs, data.table(model_name = "phase_lag1y", exposure_variable = "enso_phase_lag1y", exposure_label = "ENSO phase, 1-year lag")),
    fill = TRUE
  )
}

if ("enso_phase_lag2y" %in% names(annual_model_dt)) {
  exposure_specs <- rbindlist(
    list(exposure_specs, data.table(model_name = "phase_lag2y", exposure_variable = "enso_phase_lag2y", exposure_label = "ENSO phase, 2-year lag")),
    fill = TRUE
  )
}

if ("enso_phase_lag3y" %in% names(annual_model_dt)) {
  exposure_specs <- rbindlist(
    list(exposure_specs, data.table(model_name = "phase_lag3y", exposure_variable = "enso_phase_lag3y", exposure_label = "ENSO phase, 3-year lag")),
    fill = TRUE
  )
}

for (i in seq_len(nrow(exposure_specs))) {
  exposure_i <- exposure_specs$exposure_variable[i]
  model_name_i <- exposure_specs$model_name[i]
  exposure_label_i <- exposure_specs$exposure_label[i]

  model_dt_i <- annual_model_dt[!is.na(get(exposure_i))]

  if (nrow(model_dt_i) < 10) next

  formula_i <- build_formula(exposure_i)

  qp_i <- fit_quasipoisson(
    model_data = model_dt_i,
    formula_object = formula_i,
    model_name = model_name_i,
    exposure_label = exposure_label_i
  )

  nb_i <- fit_negative_binomial(
    model_data = model_dt_i,
    formula_object = formula_i,
    model_name = model_name_i,
    exposure_label = exposure_label_i
  )

  model_terms <- rbindlist(list(model_terms, qp_i, nb_i), fill = TRUE)
}

model_summary <- model_terms[
  term == "(Intercept)",
  .(
    model_name,
    exposure_label,
    model_family
  )
]

model_summary <- merge(model_summary, exposure_specs, by = c("model_name", "exposure_label"), all.x = TRUE)

model_summary[
  ,
  `:=`(
    years_used = NA_integer_,
    total_cases_used = NA_real_,
    total_population_used = NA_real_
  )
]

for (i in seq_len(nrow(model_summary))) {
  exposure_i <- model_summary$exposure_variable[i]
  model_dt_i <- annual_model_dt[!is.na(get(exposure_i))]

  model_summary[
    i,
    `:=`(
      years_used = nrow(model_dt_i),
      total_cases_used = sum(model_dt_i$cases, na.rm = TRUE),
      total_population_used = sum(model_dt_i$population, na.rm = TRUE)
    )
  ]
}

enso_effect_terms <- model_terms[
  term != "(Intercept)" &
    term != "year_centered" &
    !grepl("offset", term, ignore.case = TRUE)
]

enso_effect_terms[
  ,
  interpretation := fifelse(
    grepl("^oni_", term),
    "IRR per +1.0 unit increase in annual mean ONI",
    fifelse(
      grepl("El Nino", term, ignore.case = TRUE),
      "IRR compared with Neutral years",
      fifelse(
        grepl("La Nina", term, ignore.case = TRUE),
        "IRR compared with Neutral years",
        "IRR for model term"
      )
    )
  )
]

# ------------------------------------------------------------
# 10. Figures
# ------------------------------------------------------------

figure_log <- data.table()

plot_dt <- copy(annual_dt)
oni_range <- range(plot_dt$oni_current, na.rm = TRUE)
inc_range <- range(plot_dt$incidence_100k, na.rm = TRUE)

if (all(is.finite(oni_range)) && diff(oni_range) > 0 && all(is.finite(inc_range)) && diff(inc_range) > 0) {
  plot_dt[, oni_scaled := (oni_current - oni_range[1]) / diff(oni_range) * diff(inc_range) + inc_range[1]]
} else {
  plot_dt[, oni_scaled := NA_real_]
}

p_time <- ggplot(plot_dt, aes(x = year)) +
  geom_col(aes(y = incidence_100k), fill = "grey75", color = "grey35", linewidth = 0.2) +
  geom_line(aes(y = oni_scaled), color = "firebrick", linewidth = 1.1, na.rm = TRUE) +
  geom_point(aes(y = oni_scaled), color = "firebrick", size = 1.8, na.rm = TRUE) +
  geom_hline(yintercept = mean(plot_dt$incidence_100k, na.rm = TRUE), linetype = "dashed", color = "grey35") +
  scale_x_continuous(breaks = seq(primary_period_start, primary_period_end, by = 2)) +
  labs(
    title = "Annual probable acute/subacute PCM incidence and ENSO",
    subtitle = paste0(
      "Integrated SIM + SIH period: ",
      primary_period_start,
      "-",
      primary_period_end,
      "; red line = annual mean ONI scaled to the incidence axis"
    ),
    x = "Year",
    y = "Annual incidence per 100,000 inhabitants"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

figure_log <- rbindlist(
  list(
    figure_log,
    save_plot_pair(
      p_time,
      paste0("24_annual_pcm_incidence_and_oni_", timestamp),
      width = 12,
      height = 7
    )
  ),
  fill = TRUE
)

p_phase <- ggplot(
  annual_dt[!is.na(enso_phase_current)],
  aes(x = enso_phase_current, y = incidence_100k)
) +
  geom_boxplot(fill = "grey85", color = "grey25", outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.75, size = 2) +
  labs(
    title = "Annual probable acute/subacute PCM incidence by ENSO phase",
    subtitle = paste0(
      "ENSO phase defined by annual mean ONI: El Niño ≥ ",
      el_nino_threshold,
      ", La Niña ≤ ",
      la_nina_threshold
    ),
    x = "ENSO phase, same year",
    y = "Annual incidence per 100,000 inhabitants"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

figure_log <- rbindlist(
  list(
    figure_log,
    save_plot_pair(
      p_phase,
      paste0("24_incidence_by_enso_phase_current_", timestamp),
      width = 9,
      height = 7
    )
  ),
  fill = TRUE
)

continuous_terms <- enso_effect_terms[
  model_family == "quasi_poisson" &
    grepl("^oni_", term)
]

if (nrow(continuous_terms) > 0) {
  continuous_terms[
    ,
    exposure_label := factor(
      exposure_label,
      levels = c(
        "Annual mean ONI, same year",
        "Annual mean ONI, 1-year lag",
        "Annual mean ONI, 2-year lag",
        "Annual mean ONI, 3-year lag"
      )
    )
  ]

  p_oni_irr <- ggplot(
    continuous_terms,
    aes(x = exposure_label, y = irr)
  ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey35") +
    geom_pointrange(aes(ymin = conf_low_95, ymax = conf_high_95), linewidth = 0.6) +
    coord_flip() +
    labs(
      title = "Association between annual mean ONI and probable acute/subacute PCM incidence",
      subtitle = "Quasi-Poisson models with population offset and linear calendar-year term",
      x = NULL,
      y = "Incidence rate ratio per +1.0 ONI unit"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )

  figure_log <- rbindlist(
    list(
      figure_log,
      save_plot_pair(
        p_oni_irr,
        paste0("24_irr_continuous_oni_lag_models_", timestamp),
        width = 10,
        height = 6
      )
    ),
    fill = TRUE
  )
}

phase_terms <- enso_effect_terms[
  model_family == "quasi_poisson" &
    grepl("phase", model_name, ignore.case = TRUE)
]

if (nrow(phase_terms) > 0) {
  phase_terms[, term_clean := gsub("enso_phase_current|enso_phase_lag1y|enso_phase_lag2y|enso_phase_lag3y", "", term)]

  p_phase_irr <- ggplot(
    phase_terms,
    aes(x = exposure_label, y = irr, shape = term_clean)
  ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey35") +
    geom_pointrange(
      aes(ymin = conf_low_95, ymax = conf_high_95),
      position = position_dodge(width = 0.55),
      linewidth = 0.6
    ) +
    coord_flip() +
    labs(
      title = "Association between ENSO phase and probable acute/subacute PCM incidence",
      subtitle = "Quasi-Poisson models with Neutral years as reference",
      x = NULL,
      y = "Incidence rate ratio compared with Neutral years",
      shape = "ENSO phase"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )

  figure_log <- rbindlist(
    list(
      figure_log,
      save_plot_pair(
        p_phase_irr,
        paste0("24_irr_enso_phase_lag_models_", timestamp),
        width = 10,
        height = 7
      )
    ),
    fill = TRUE
  )
}

# ------------------------------------------------------------
# 11. Save outputs
# ------------------------------------------------------------

annual_dataset_file <- file.path(enso_analysis_data_dir, paste0("24_enso_pcm_annual_dataset_", timestamp, ".rds"))
annual_dataset_latest <- file.path(enso_analysis_data_dir, "enso_pcm_annual_dataset_latest.rds")
annual_dataset_csv <- file.path(enso_analysis_tables_dir, paste0("24_enso_pcm_annual_dataset_", timestamp, ".csv"))
annual_dataset_csv_latest <- file.path(enso_analysis_tables_dir, "enso_pcm_annual_dataset_latest.csv")

saveRDS(annual_dt, annual_dataset_file)
saveRDS(annual_dt, annual_dataset_latest)
fwrite(annual_dt, annual_dataset_csv)
fwrite(annual_dt, annual_dataset_csv_latest)

model_terms_file <- file.path(enso_analysis_tables_dir, paste0("24_enso_pcm_model_terms_", timestamp, ".csv"))
model_terms_latest <- file.path(enso_analysis_tables_dir, "enso_pcm_model_terms_latest.csv")
model_summary_file <- file.path(enso_analysis_tables_dir, paste0("24_enso_pcm_model_summary_", timestamp, ".csv"))
model_summary_latest <- file.path(enso_analysis_tables_dir, "enso_pcm_model_summary_latest.csv")
effect_terms_file <- file.path(enso_analysis_tables_dir, paste0("24_enso_pcm_effect_terms_", timestamp, ".csv"))
effect_terms_latest <- file.path(enso_analysis_tables_dir, "enso_pcm_effect_terms_latest.csv")
phase_summary_current_file <- file.path(enso_analysis_tables_dir, paste0("24_enso_phase_summary_current_", timestamp, ".csv"))
phase_summary_current_latest <- file.path(enso_analysis_tables_dir, "enso_phase_summary_current_latest.csv")
phase_summary_lag1_file <- file.path(enso_analysis_tables_dir, paste0("24_enso_phase_summary_lag1_", timestamp, ".csv"))
phase_summary_lag1_latest <- file.path(enso_analysis_tables_dir, "enso_phase_summary_lag1_latest.csv")
annual_summary_file <- file.path(enso_analysis_tables_dir, paste0("24_enso_pcm_annual_summary_", timestamp, ".csv"))
annual_summary_latest <- file.path(enso_analysis_tables_dir, "enso_pcm_annual_summary_latest.csv")
figure_log_file <- file.path(enso_analysis_tables_dir, paste0("24_enso_pcm_figure_log_", timestamp, ".csv"))
figure_log_latest <- file.path(enso_analysis_tables_dir, "enso_pcm_figure_log_latest.csv")
input_summary_file <- file.path(enso_analysis_tables_dir, paste0("24_enso_pcm_input_file_summary_", timestamp, ".csv"))
input_summary_latest <- file.path(enso_analysis_tables_dir, "enso_pcm_input_file_summary_latest.csv")

fwrite(model_terms, model_terms_file)
fwrite(model_terms, model_terms_latest)
fwrite(model_summary, model_summary_file)
fwrite(model_summary, model_summary_latest)
fwrite(enso_effect_terms, effect_terms_file)
fwrite(enso_effect_terms, effect_terms_latest)
fwrite(phase_summary_current, phase_summary_current_file)
fwrite(phase_summary_current, phase_summary_current_latest)
fwrite(phase_summary_lag1, phase_summary_lag1_file)
fwrite(phase_summary_lag1, phase_summary_lag1_latest)
fwrite(annual_summary, annual_summary_file)
fwrite(annual_summary, annual_summary_latest)
fwrite(figure_log, figure_log_file)
fwrite(figure_log, figure_log_latest)
fwrite(input_file_summary, input_summary_file)
fwrite(input_file_summary, input_summary_latest)

for (file_i in c(
  annual_dataset_csv,
  model_terms_file,
  model_summary_file,
  effect_terms_file,
  phase_summary_current_file,
  phase_summary_lag1_file,
  annual_summary_file,
  figure_log_file,
  input_summary_file
)) {
  if (!is.na(file_i) && file.exists(file_i)) {
    file.copy(from = file_i, to = file.path(tables_dir, basename(file_i)), overwrite = TRUE)
  }
}

# ------------------------------------------------------------
# 12. Build interpretation text
# ------------------------------------------------------------

get_effect_row <- function(effect_dt, model_name_i, term_pattern_i, family_i = "quasi_poisson") {
  effect_dt[
    model_name == model_name_i &
      model_family == family_i &
      grepl(term_pattern_i, term, ignore.case = TRUE)
  ][1]
}

current_oni_effect <- get_effect_row(enso_effect_terms, "continuous_ONI_current", "^oni_current$")
lag1_oni_effect <- get_effect_row(enso_effect_terms, "continuous_ONI_lag1y", "^oni_lag1y$")
phase_el_nino_effect <- get_effect_row(enso_effect_terms, "phase_current", "El Nino")
phase_lag1_el_nino_effect <- get_effect_row(enso_effect_terms, "phase_lag1y", "El Nino")

interpretation_text <- paste0(
  "ENSO and probable acute/subacute PCM incidence analysis\n",
  "=======================================================\n\n",
  "This analysis evaluated the association between El Niño/ENSO and annual incidence of ",
  outcome_definition,
  " in Brazil using the common integrated SIM + SIH period ",
  primary_period_start,
  "-",
  primary_period_end,
  ".\n\n",
  "Important methodological note\n",
  "-----------------------------\n",
  "ENSO is a year-level climatic exposure. Therefore, the primary analysis was performed at the national annual level rather than by treating municipality-year rows as independent ENSO observations. ",
  "Counts were modeled with a population offset, and the main models included a linear calendar-year term to reduce confounding by long-term temporal trends. ",
  "Because the annual time series is short, results should be interpreted as exploratory evidence of association, not as causal proof.\n\n",
  "Annual dataset summary\n",
  "----------------------\n",
  "Years analyzed: ", fmt_int(annual_summary$years), " (", fmt_int(annual_summary$first_year), "-", fmt_int(annual_summary$last_year), ")\n",
  "Total events: ", fmt_int(annual_summary$total_cases), "\n",
  "Median annual incidence: ", fmt_num(annual_summary$median_annual_incidence_100k, 4), " per 100,000\n",
  "Mean annual incidence: ", fmt_num(annual_summary$mean_annual_incidence_100k, 4), " per 100,000\n\n",
  "Main quasi-Poisson estimates\n",
  "----------------------------\n",
  "Annual mean ONI, same year: IRR ", fmt_num(current_oni_effect$irr, 3), " (95% CI ", fmt_num(current_oni_effect$conf_low_95, 3), "-", fmt_num(current_oni_effect$conf_high_95, 3), "), p = ", fmt_num(current_oni_effect$p_value, 4), "\n",
  "Annual mean ONI, 1-year lag: IRR ", fmt_num(lag1_oni_effect$irr, 3), " (95% CI ", fmt_num(lag1_oni_effect$conf_low_95, 3), "-", fmt_num(lag1_oni_effect$conf_high_95, 3), "), p = ", fmt_num(lag1_oni_effect$p_value, 4), "\n",
  "El Niño phase, same year versus Neutral: IRR ", fmt_num(phase_el_nino_effect$irr, 3), " (95% CI ", fmt_num(phase_el_nino_effect$conf_low_95, 3), "-", fmt_num(phase_el_nino_effect$conf_high_95, 3), "), p = ", fmt_num(phase_el_nino_effect$p_value, 4), "\n",
  "El Niño phase, 1-year lag versus Neutral: IRR ", fmt_num(phase_lag1_el_nino_effect$irr, 3), " (95% CI ", fmt_num(phase_lag1_el_nino_effect$conf_low_95, 3), "-", fmt_num(phase_lag1_el_nino_effect$conf_high_95, 3), "), p = ", fmt_num(phase_lag1_el_nino_effect$p_value, 4), "\n\n",
  "Interpretation guidance\n",
  "-----------------------\n",
  "If IRR values are above 1, they suggest higher annual incidence during or after warmer ENSO conditions; if below 1, they suggest lower annual incidence. ",
  "Confidence intervals crossing 1 indicate substantial statistical uncertainty. ",
  "Lagged ONI models are biologically plausible because PCM clinical records may occur after environmental exposure, infection, disease progression, diagnosis, and hospitalization or death registration.\n"
)

interpretation_file <- file.path(enso_analysis_dir, paste0("24_enso_pcm_interpretation_text_", timestamp, ".txt"))
interpretation_latest <- file.path(enso_analysis_dir, "enso_pcm_interpretation_text_latest.txt")
write_text_file(interpretation_text, interpretation_file)
write_text_file(interpretation_text, interpretation_latest)

# ------------------------------------------------------------
# 13. Save report
# ------------------------------------------------------------

report_file <- file.path(logs_dir, paste0("24_analyze_enso_impact_on_pcm_incidence_report_", timestamp, ".txt"))

sink(report_file)

cat("============================================================\n")
cat("ENSO impact on probable acute/subacute PCM incidence report\n")
cat("Project: Ecoepidemiology of acute/subacute paracoccidioidomycosis in Brazil\n")
cat("Generated at:", as.character(Sys.time()), "\n")
cat("Script: 24_analyze_enso_impact_on_pcm_incidence.R\n")
cat("Version: v1_enso_pcm_incidence\n")
cat("============================================================\n\n")

cat("Purpose:\n")
cat("Analyze the association between El Niño/ENSO and annual incidence of probable acute/subacute PCM using the available integrated record period.\n\n")

cat("Primary analytical period:\n")
cat(primary_period_start, "-", primary_period_end, "\n\n")

cat("Methodological note:\n")
cat("The integrated SIM + SIH outcome uses the common period available for both systems.\n")
cat("ENSO is a year-level exposure; therefore, the primary analysis is national annual rather than municipality-year independent modeling.\n")
cat("Counts are modeled with a population offset and a linear year term.\n")
cat("The analysis is exploratory and should not be interpreted as causal proof.\n\n")

cat("Input files:\n")
print(input_file_summary)
cat("\n\n")

cat("Detected columns:\n")
cat("Year column: year\n")
cat("Municipality column: municipality_code_6\n")
cat("Outcome column used as pcm_cases:", outcome_col, "\n")
cat("Population column used as population:", population_col, "\n")
cat("Current ONI column:", oni_current_col, "\n")
cat("Lag 1 ONI column:", oni_lag1_col, "\n")
cat("Lag 2 ONI column:", oni_lag2_col, "\n")
cat("Lag 3 ONI column:", oni_lag3_col, "\n\n")

cat("Annual summary:\n")
print(annual_summary)
cat("\n\n")

cat("Annual dataset:\n")
print(annual_dt)
cat("\n\n")

cat("ENSO phase summary, same year:\n")
print(phase_summary_current)
cat("\n\n")

cat("ENSO phase summary, 1-year lag:\n")
print(phase_summary_lag1)
cat("\n\n")

cat("Model summary:\n")
print(model_summary)
cat("\n\n")

cat("ENSO effect terms:\n")
print(enso_effect_terms)
cat("\n\n")

cat("Figure log:\n")
print(figure_log)
cat("\n\n")

cat("Interpretation text:\n")
cat(interpretation_text)
cat("\n\n")

cat("Output folders:\n")
cat("Analysis folder:", enso_analysis_dir, "\n")
cat("Tables:", enso_analysis_tables_dir, "\n")
cat("Figures:", enso_analysis_figures_dir, "\n")
cat("Data:", enso_analysis_data_dir, "\n")
cat("General output tables:", tables_dir, "\n")
cat("General output figures:", figures_dir, "\n")

sink()

# ------------------------------------------------------------
# 14. Console summary
# ------------------------------------------------------------

message("============================================================")
message("ENSO impact analysis completed.")
message("Script version: v1_enso_pcm_incidence")
message("Annual years analyzed: ", nrow(annual_dt))
message("Total events: ", sum(annual_dt$cases, na.rm = TRUE))
message("Median annual incidence per 100,000: ", fmt_num(annual_summary$median_annual_incidence_100k, 4))
message("Effect terms rows: ", nrow(enso_effect_terms))
message("Figures created: ", nrow(figure_log))
message("Report saved to: ", report_file)
message("============================================================")

print(
  data.table(
    output = c(
      "annual_dataset_csv",
      "effect_terms_latest",
      "phase_summary_current_latest",
      "interpretation_latest",
      "figure_log_latest",
      "report_file"
    ),
    path = c(
      annual_dataset_csv_latest,
      effect_terms_latest,
      phase_summary_current_latest,
      interpretation_latest,
      figure_log_latest,
      report_file
    )
  )
)
