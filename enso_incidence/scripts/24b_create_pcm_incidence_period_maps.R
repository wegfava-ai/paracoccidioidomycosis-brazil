# ============================================================
# Project: Ecoepidemiology of acute/subacute paracoccidioidomycosis in Brazil
# Script: 24b_create_pcm_incidence_period_maps.R
# Purpose: Create national maps of annualized municipal incidence of probable acute/subacute PCM for the overall period and selected subperiods
# Author: Wellington Fava
# Date: 2026-06-30
# Version: v1_incidence_period_maps
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
maps_dir <- file.path(outputs_dir, "maps")

data_processed_dir <- file.path(project_root, "data_processed")
spatial_processed_dir <- file.path(data_processed_dir, "spatial")
models_dir <- file.path(data_processed_dir, "models")

incidence_maps_dir <- file.path(models_dir, "pcm_incidence_period_maps")
incidence_maps_tables_dir <- file.path(incidence_maps_dir, "tables")
incidence_maps_maps_dir <- file.path(incidence_maps_dir, "maps")
incidence_maps_data_dir <- file.path(incidence_maps_dir, "data")

for (dir_path in c(
  logs_dir,
  outputs_dir,
  tables_dir,
  maps_dir,
  data_processed_dir,
  spatial_processed_dir,
  models_dir,
  incidence_maps_dir,
  incidence_maps_tables_dir,
  incidence_maps_maps_dir,
  incidence_maps_data_dir
)) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("Created folder: ", dir_path)
  }
}

# ------------------------------------------------------------
# 2. Settings
# ------------------------------------------------------------

period_definitions <- data.frame(
  period_label = c(
    "Overall period | 1998-2024",
    "Period 1 | 1998-2004",
    "Period 2 | 2005-2011",
    "Period 3 | 2012-2018",
    "Period 4 | 2019-2024"
  ),
  period_id = c(
    "overall_1998_2024",
    "period_1_1998_2004",
    "period_2_2005_2011",
    "period_3_2012_2018",
    "period_4_2019_2024"
  ),
  start_year = c(1998, 1998, 2005, 2012, 2019),
  end_year = c(2024, 2004, 2011, 2018, 2024),
  stringsAsFactors = FALSE
)

# Incidence is computed as annualized incidence per 100,000 person-years:
# sum(cases during period) / sum(population denominators during period) * 100,000.
incidence_multiplier <- 100000

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# ------------------------------------------------------------
# 3. Packages
# ------------------------------------------------------------

required_packages <- c(
  "data.table",
  "sf",
  "ggplot2"
)

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
library(sf)
library(ggplot2)

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

safe_file_stem <- function(x) {
  x <- safe_character(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
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
  dt_names <- names(dt)

  name_map <- data.table(
    original_name = dt_names,
    lower_name = tolower(dt_names)
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
      paste(dt_names, collapse = ", ")
    )
  }

  NA_character_
}

read_spatial_object <- function(file_path) {
  if (grepl("\\.rds$", file_path, ignore.case = TRUE)) {
    obj <- readRDS(file_path)
  } else if (grepl("\\.gpkg$|\\.shp$", file_path, ignore.case = TRUE)) {
    obj <- sf::st_read(file_path, quiet = TRUE)
  } else {
    stop("Unsupported spatial file type: ", file_path)
  }

  if (!inherits(obj, "sf")) {
    stop("The selected file does not contain an sf object: ", file_path)
  }

  obj
}

standardize_crs <- function(x, target_crs = 4674) {
  if (is.na(sf::st_crs(x))) {
    sf::st_crs(x) <- target_crs
  }
  if (sf::st_crs(x)$epsg != target_crs) {
    x <- sf::st_transform(x, target_crs)
  }
  x
}

save_map_pair <- function(plot_object, file_stem, width = 11, height = 9, dpi = 320) {
  png_file <- file.path(incidence_maps_maps_dir, paste0(file_stem, ".png"))
  pdf_file <- file.path(incidence_maps_maps_dir, paste0(file_stem, ".pdf"))

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
    to = file.path(maps_dir, basename(png_file)),
    overwrite = TRUE
  )

  file.copy(
    from = pdf_file,
    to = file.path(maps_dir, basename(pdf_file)),
    overwrite = TRUE
  )

  data.table(
    map_stem = file_stem,
    png_file = normalize_path(png_file),
    pdf_file = normalize_path(pdf_file),
    png_file_outputs = normalize_path(file.path(maps_dir, basename(png_file))),
    pdf_file_outputs = normalize_path(file.path(maps_dir, basename(pdf_file)))
  )
}

# ------------------------------------------------------------
# 5. Locate and read national municipality-year spatial dataset
# ------------------------------------------------------------

national_spatial_file <- find_required_file(
  candidate_paths = c(
    file.path(spatial_processed_dir, "pcm_municipality_spatial_main_year_1998_2024_latest.rds"),
    find_latest_file(spatial_processed_dir, "pcm_municipality_spatial_main_year_1998_2024.*\\.rds$"),
    file.path(spatial_processed_dir, "pcm_municipality_spatial_main_year_1998_2024_latest.gpkg"),
    find_latest_file(spatial_processed_dir, "pcm_municipality_spatial_main_year_1998_2024.*\\.gpkg$")
  ),
  file_description = "full national municipality-year spatial PCM dataset"
)

national_sf <- read_spatial_object(national_spatial_file)
national_sf <- standardize_crs(national_sf, target_crs = 4674)

input_file_summary <- data.table(
  input_name = "national_spatial_file",
  file_path = national_spatial_file
)

# ------------------------------------------------------------
# 6. Standardize columns
# ------------------------------------------------------------

year_col <- find_column(
  national_sf,
  c("year", "ano", "event_year", "YEAR"),
  description = "year column"
)

municipality_col <- find_column(
  national_sf,
  c("municipality_code_6", "codmun6", "municipality_code", "codigo_municipio", "code_muni"),
  description = "municipality code column"
)

municipality_name_col <- find_column(
  national_sf,
  c("municipality_name_final", "municipality_name", "name_muni", "nome_municipio", "NM_MUN"),
  required = FALSE,
  description = "municipality name column"
)

uf_col <- find_column(
  national_sf,
  c("uf_abbrev", "abbrev_state", "state_abbrev", "sigla_uf", "SIGLA_UF", "uf"),
  required = FALSE,
  description = "state abbreviation column"
)

state_code_col <- find_column(
  national_sf,
  c("state_code", "code_state", "cod_uf", "CD_UF", "uf_code"),
  required = FALSE,
  description = "state code column"
)

region_col <- find_column(
  national_sf,
  c("region", "region_name", "nome_regiao", "NM_REGIAO"),
  required = FALSE,
  description = "region column"
)

outcome_col <- find_column(
  national_sf,
  c(
    "pcm_event_count_main_final",
    "pcm_event_count_main",
    "pcm_event_count",
    "pcm_records_under20_main_period",
    "pcm_records_main_period",
    "main_records",
    "event_count",
    "cases",
    "n_cases"
  ),
  description = "PCM event count column"
)

population_col <- find_column(
  national_sf,
  c(
    "population_final",
    "population",
    "pop",
    "population_denominator",
    "populacao"
  ),
  description = "population denominator column"
)

setnames(national_sf, year_col, "year")
setnames(national_sf, municipality_col, "municipality_code_6")
setnames(national_sf, outcome_col, "pcm_cases")
setnames(national_sf, population_col, "population")

if (!is.na(municipality_name_col) && municipality_name_col != "municipality_name") {
  setnames(national_sf, municipality_name_col, "municipality_name")
}

if (!is.na(uf_col) && uf_col != "uf_abbrev") {
  setnames(national_sf, uf_col, "uf_abbrev")
}

if (!is.na(state_code_col) && state_code_col != "state_code") {
  setnames(national_sf, state_code_col, "state_code")
}

if (!is.na(region_col) && region_col != "region") {
  setnames(national_sf, region_col, "region")
}

national_sf$year <- as.integer(national_sf$year)
national_sf$municipality_code_6 <- safe_character(national_sf$municipality_code_6)
national_sf$pcm_cases <- safe_numeric(national_sf$pcm_cases)
national_sf$population <- safe_numeric(national_sf$population)
national_sf$pcm_cases[is.na(national_sf$pcm_cases)] <- 0

if (!("municipality_name" %in% names(national_sf))) {
  national_sf$municipality_name <- national_sf$municipality_code_6
}

if (!("uf_abbrev" %in% names(national_sf))) {
  national_sf$uf_abbrev <- NA_character_
}

if (!("state_code" %in% names(national_sf))) {
  national_sf$state_code <- NA_character_
}

if (!("region" %in% names(national_sf))) {
  national_sf$region <- NA_character_
}

# ------------------------------------------------------------
# 7. Build stable municipality geometry and boundaries
# ------------------------------------------------------------

# Geometry for each municipality should be repeated by year; keep one geometry per municipality.
municipality_base_sf <- national_sf[
  !duplicated(national_sf$municipality_code_6),
  c("municipality_code_6", "municipality_name", "uf_abbrev", "state_code", "region", attr(national_sf, "sf_column"))
]

municipality_base_sf <- standardize_crs(municipality_base_sf, target_crs = 4674)

# State boundaries from municipality union.
boundary_group_col <- NA_character_

if ("state_code" %in% names(municipality_base_sf) &&
    any(!is.na(municipality_base_sf$state_code))) {
  boundary_group_col <- "state_code"
} else if ("uf_abbrev" %in% names(municipality_base_sf) &&
           any(!is.na(municipality_base_sf$uf_abbrev))) {
  boundary_group_col <- "uf_abbrev"
}

if (!is.na(boundary_group_col)) {
  state_boundaries <- aggregate(
    municipality_base_sf[, boundary_group_col],
    by = list(municipality_base_sf[[boundary_group_col]]),
    FUN = function(x) x[1]
  )

  names(state_boundaries)[1] <- "state_group"
  state_boundaries <- standardize_crs(state_boundaries, target_crs = 4674)
} else {
  state_boundaries <- sf::st_sf(
    state_group = character(0),
    geometry = sf::st_sfc(crs = sf::st_crs(municipality_base_sf))
  )
}

brazil_boundary <- sf::st_as_sf(
  data.frame(country = "Brazil"),
  geometry = sf::st_union(sf::st_geometry(municipality_base_sf))
)

brazil_boundary <- standardize_crs(brazil_boundary, target_crs = 4674)

# ------------------------------------------------------------
# 8. Aggregate incidence by period
# ------------------------------------------------------------

national_dt <- as.data.table(sf::st_drop_geometry(national_sf))

period_results <- list()

for (i in seq_len(nrow(period_definitions))) {

  period_i <- period_definitions[i, ]

  period_dt_i <- national_dt[
    year >= period_i$start_year &
      year <= period_i$end_year
  ]

  if (nrow(period_dt_i) == 0) {
    warning("No rows found for period: ", period_i$period_label)
    next
  }

  agg_i <- period_dt_i[
    ,
    .(
      cases = sum(pcm_cases, na.rm = TRUE),
      person_years = sum(population, na.rm = TRUE),
      years_in_period = uniqueN(year),
      rows = .N,
      event_years = sum(pcm_cases > 0, na.rm = TRUE),
      missing_population_rows = sum(is.na(population) | !is.finite(population) | population <= 0),
      municipality_name = municipality_name[which.max(!is.na(municipality_name))][1],
      uf_abbrev = uf_abbrev[which.max(!is.na(uf_abbrev))][1],
      state_code = state_code[which.max(!is.na(state_code))][1],
      region = region[which.max(!is.na(region))][1]
    ),
    by = municipality_code_6
  ]

  agg_i[
    ,
    `:=`(
      annualized_incidence_100k = fifelse(
        is.finite(person_years) & person_years > 0,
        cases / person_years * incidence_multiplier,
        NA_real_
      ),
      cumulative_cases = cases,
      period_label = period_i$period_label,
      period_id = period_i$period_id,
      start_year = period_i$start_year,
      end_year = period_i$end_year
    )
  ]

  period_results[[period_i$period_id]] <- agg_i
}

period_incidence_dt <- rbindlist(period_results, fill = TRUE)

period_incidence_sf <- merge(
  municipality_base_sf,
  period_incidence_dt,
  by = "municipality_code_6",
  all.y = TRUE
)

period_summary <- period_incidence_dt[
  ,
  .(
    municipalities = .N,
    municipalities_with_valid_incidence = sum(!is.na(annualized_incidence_100k)),
    municipalities_with_cases = sum(cases > 0, na.rm = TRUE),
    total_cases = sum(cases, na.rm = TRUE),
    total_person_years = sum(person_years, na.rm = TRUE),
    national_annualized_incidence_100k = sum(cases, na.rm = TRUE) / sum(person_years, na.rm = TRUE) * incidence_multiplier,
    median_municipal_incidence_100k = median(annualized_incidence_100k, na.rm = TRUE),
    mean_municipal_incidence_100k = mean(annualized_incidence_100k, na.rm = TRUE),
    max_municipal_incidence_100k = max(annualized_incidence_100k, na.rm = TRUE)
  ),
  by = .(period_id, period_label, start_year, end_year)
][order(start_year, end_year)]

top_municipalities <- period_incidence_dt[
  !is.na(annualized_incidence_100k),
  .SD[order(-annualized_incidence_100k, -cases)][1:min(.N, 25)],
  by = .(period_id, period_label)
][
  ,
  .(
    period_id,
    period_label,
    municipality_code_6,
    municipality_name,
    uf_abbrev,
    region,
    cases,
    person_years,
    annualized_incidence_100k,
    years_in_period
  )
]

# ------------------------------------------------------------
# 9. Save tables before plotting
# ------------------------------------------------------------

period_incidence_csv <- file.path(
  incidence_maps_tables_dir,
  paste0("24b_pcm_municipal_annualized_incidence_by_period_", timestamp, ".csv")
)

period_incidence_latest <- file.path(
  incidence_maps_tables_dir,
  "pcm_municipal_annualized_incidence_by_period_latest.csv"
)

period_summary_csv <- file.path(
  incidence_maps_tables_dir,
  paste0("24b_pcm_incidence_period_summary_", timestamp, ".csv")
)

period_summary_latest <- file.path(
  incidence_maps_tables_dir,
  "pcm_incidence_period_summary_latest.csv"
)

top_municipalities_csv <- file.path(
  incidence_maps_tables_dir,
  paste0("24b_top_municipal_incidence_by_period_", timestamp, ".csv")
)

top_municipalities_latest <- file.path(
  incidence_maps_tables_dir,
  "top_municipal_incidence_by_period_latest.csv"
)

input_summary_csv <- file.path(
  incidence_maps_tables_dir,
  paste0("24b_incidence_maps_input_file_summary_", timestamp, ".csv")
)

input_summary_latest <- file.path(
  incidence_maps_tables_dir,
  "incidence_maps_input_file_summary_latest.csv"
)

fwrite(period_incidence_dt, period_incidence_csv)
fwrite(period_incidence_dt, period_incidence_latest)
fwrite(period_summary, period_summary_csv)
fwrite(period_summary, period_summary_latest)
fwrite(top_municipalities, top_municipalities_csv)
fwrite(top_municipalities, top_municipalities_latest)
fwrite(input_file_summary, input_summary_csv)
fwrite(input_file_summary, input_summary_latest)

for (file_i in c(
  period_incidence_csv,
  period_summary_csv,
  top_municipalities_csv,
  input_summary_csv
)) {
  file.copy(
    from = file_i,
    to = file.path(tables_dir, basename(file_i)),
    overwrite = TRUE
  )
}

period_incidence_gpkg <- file.path(
  incidence_maps_data_dir,
  paste0("24b_pcm_incidence_period_maps_", timestamp, ".gpkg")
)

period_incidence_gpkg_latest <- file.path(
  incidence_maps_data_dir,
  "pcm_incidence_period_maps_latest.gpkg"
)

sf::st_write(period_incidence_sf, period_incidence_gpkg, delete_dsn = TRUE, quiet = TRUE)
sf::st_write(period_incidence_sf, period_incidence_gpkg_latest, delete_dsn = TRUE, quiet = TRUE)

# ------------------------------------------------------------
# 10. Create incidence maps
# ------------------------------------------------------------

map_log <- data.table()

# Common upper limit improves comparability across period maps.
global_max <- max(period_incidence_sf$annualized_incidence_100k, na.rm = TRUE)

if (!is.finite(global_max) || global_max <= 0) {
  global_max <- 1
}

for (i in seq_len(nrow(period_definitions))) {

  period_i <- period_definitions[i, ]

  map_sf_i <- period_incidence_sf[
    period_incidence_sf$period_id == period_i$period_id,
  ]

  summary_i <- period_summary[
    period_id == period_i$period_id
  ]

  title_i <- paste0(
    "Annualized incidence of probable acute/subacute PCM, ",
    period_i$start_year,
    "-",
    period_i$end_year
  )

  subtitle_i <- paste0(
    "Cases = ",
    formatC(summary_i$total_cases, format = "d", big.mark = ","),
    "; national annualized incidence = ",
    formatC(summary_i$national_annualized_incidence_100k, format = "f", digits = 4),
    " per 100,000 person-years"
  )

  file_stem_i <- paste0(
    "24b_pcm_annualized_incidence_100k_",
    period_i$period_id,
    "_",
    timestamp
  )

  p_map_i <- ggplot() +
    geom_sf(
      data = map_sf_i,
      aes(fill = annualized_incidence_100k),
      color = NA
    ) +
    geom_sf(
      data = state_boundaries,
      fill = NA,
      color = "grey25",
      linewidth = 0.22
    ) +
    geom_sf(
      data = brazil_boundary,
      fill = NA,
      color = "black",
      linewidth = 0.45
    ) +
    scale_fill_viridis_c(
      option = "magma",
      direction = -1,
      trans = "sqrt",
      limits = c(0, global_max),
      na.value = "grey92",
      name = "Annualized\nincidence\nper 100,000"
    ) +
    labs(
      title = title_i,
      subtitle = subtitle_i,
      caption = "Incidence = sum of probable acute/subacute PCM records divided by summed population person-years × 100,000. State and national borders are shown in grey/black."
    ) +
    coord_sf(datum = NA) +
    theme_void(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold", size = 15, hjust = 0),
      plot.subtitle = element_text(size = 11, hjust = 0, margin = margin(b = 8)),
      plot.caption = element_text(size = 8.5, hjust = 0, margin = margin(t = 8)),
      legend.position = "right"
    )

  map_log <- rbindlist(
    list(
      map_log,
      cbind(
        data.table(
          period_id = period_i$period_id,
          period_label = period_i$period_label,
          start_year = period_i$start_year,
          end_year = period_i$end_year,
          map_type = "annualized_incidence_100k",
          border_layers = "Brazil boundary and state boundaries"
        ),
        save_map_pair(
          p_map_i,
          file_stem_i,
          width = 10.5,
          height = 9,
          dpi = 320
        )
      )
    ),
    fill = TRUE
  )
}

# ------------------------------------------------------------
# 11. Create a multi-panel period map
# ------------------------------------------------------------

period_only_sf <- period_incidence_sf[
  period_incidence_sf$period_id != "overall_1998_2024",
]

if (nrow(period_only_sf) > 0) {

  p_panel <- ggplot() +
    geom_sf(
      data = period_only_sf,
      aes(fill = annualized_incidence_100k),
      color = NA
    ) +
    geom_sf(
      data = state_boundaries,
      fill = NA,
      color = "grey25",
      linewidth = 0.18
    ) +
    geom_sf(
      data = brazil_boundary,
      fill = NA,
      color = "black",
      linewidth = 0.35
    ) +
    facet_wrap(~ period_label, ncol = 2) +
    scale_fill_viridis_c(
      option = "magma",
      direction = -1,
      trans = "sqrt",
      limits = c(0, global_max),
      na.value = "grey92",
      name = "Annualized\nincidence\nper 100,000"
    ) +
    labs(
      title = "Annualized incidence of probable acute/subacute PCM by period",
      subtitle = "Brazilian municipalities, 1998-2024",
      caption = "Incidence = sum of probable acute/subacute PCM records divided by summed population person-years × 100,000. State and national borders are shown in grey/black."
    ) +
    coord_sf(datum = NA) +
    theme_void(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold", size = 16, hjust = 0),
      plot.subtitle = element_text(size = 11, hjust = 0, margin = margin(b = 8)),
      plot.caption = element_text(size = 8.5, hjust = 0, margin = margin(t = 8)),
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "right"
    )

  map_log <- rbindlist(
    list(
      map_log,
      cbind(
        data.table(
          period_id = "period_panel_1998_2024",
          period_label = "Period panel | 1998-2024",
          start_year = 1998,
          end_year = 2024,
          map_type = "annualized_incidence_100k_period_panel",
          border_layers = "Brazil boundary and state boundaries"
        ),
        save_map_pair(
          p_panel,
          paste0("24b_pcm_annualized_incidence_100k_period_panel_1998_2024_", timestamp),
          width = 13,
          height = 10,
          dpi = 320
        )
      )
    ),
    fill = TRUE
  )
}

# ------------------------------------------------------------
# 12. Save map log and captions
# ------------------------------------------------------------

map_log_file <- file.path(
  incidence_maps_tables_dir,
  paste0("24b_pcm_incidence_map_file_log_", timestamp, ".csv")
)

map_log_latest <- file.path(
  incidence_maps_tables_dir,
  "pcm_incidence_map_file_log_latest.csv"
)

fwrite(map_log, map_log_file)
fwrite(map_log, map_log_latest)

file.copy(
  from = map_log_file,
  to = file.path(tables_dir, basename(map_log_file)),
  overwrite = TRUE
)

map_captions <- map_log[
  ,
  .(
    period_id,
    period_label,
    map_type,
    map_caption = paste0(
      "Annualized municipal incidence of probable acute/subacute paracoccidioidomycosis in Brazil for ",
      period_label,
      ". Incidence was calculated as the total number of records divided by summed population person-years and multiplied by 100,000. ",
      "Municipal polygons are colored by annualized incidence; state boundaries and the national boundary are overlaid to improve spatial orientation."
    ),
    png_file,
    pdf_file
  )
]

map_captions_file <- file.path(
  incidence_maps_tables_dir,
  paste0("24b_pcm_incidence_map_captions_", timestamp, ".csv")
)

map_captions_latest <- file.path(
  incidence_maps_tables_dir,
  "pcm_incidence_map_captions_latest.csv"
)

fwrite(map_captions, map_captions_file)
fwrite(map_captions, map_captions_latest)

# ------------------------------------------------------------
# 13. Save report
# ------------------------------------------------------------

report_file <- file.path(
  logs_dir,
  paste0("24b_create_pcm_incidence_period_maps_report_", timestamp, ".txt")
)

sink(report_file)

cat("============================================================\n")
cat("PCM incidence period maps report\n")
cat("Project: Ecoepidemiology of acute/subacute paracoccidioidomycosis in Brazil\n")
cat("Generated at:", as.character(Sys.time()), "\n")
cat("Script: 24b_create_pcm_incidence_period_maps.R\n")
cat("Version: v1_incidence_period_maps\n")
cat("============================================================\n\n")

cat("Purpose:\n")
cat("Create national maps of annualized municipal incidence of probable acute/subacute PCM for the overall period and selected subperiods.\n\n")

cat("Period definitions:\n")
print(as.data.table(period_definitions))
cat("\n\n")

cat("Methodological note:\n")
cat("Incidence is annualized within each period: sum of cases divided by summed population person-years multiplied by 100,000.\n")
cat("The maps use the full national municipality-year spatial dataset when available.\n")
cat("All maps include municipality fills, state borders, and the national Brazil border.\n\n")

cat("Input files:\n")
print(input_file_summary)
cat("\n\n")

cat("Detected columns:\n")
cat("Year column: year\n")
cat("Municipality code column: municipality_code_6\n")
cat("Municipality name column:", ifelse("municipality_name" %in% names(national_sf), "municipality_name", NA_character_), "\n")
cat("UF column:", ifelse("uf_abbrev" %in% names(national_sf), "uf_abbrev", NA_character_), "\n")
cat("State code column:", ifelse("state_code" %in% names(national_sf), "state_code", NA_character_), "\n")
cat("Region column:", ifelse("region" %in% names(national_sf), "region", NA_character_), "\n")
cat("Outcome column used as pcm_cases:", outcome_col, "\n")
cat("Population column used as population:", population_col, "\n")
cat("Boundary group column:", boundary_group_col, "\n\n")

cat("Spatial summary:\n")
cat("National spatial rows:", nrow(national_sf), "\n")
cat("Municipality base rows:", nrow(municipality_base_sf), "\n")
cat("State boundary rows:", nrow(state_boundaries), "\n")
cat("Brazil boundary rows:", nrow(brazil_boundary), "\n")
cat("CRS:", sf::st_crs(municipality_base_sf)$input, "\n\n")

cat("Period summary:\n")
print(period_summary)
cat("\n\n")

cat("Top municipalities by period:\n")
print(top_municipalities)
cat("\n\n")

cat("Map log:\n")
print(map_log)
cat("\n\n")

cat("Output files:\n")
cat("Municipal incidence table:", period_incidence_latest, "\n")
cat("Period summary:", period_summary_latest, "\n")
cat("Top municipalities:", top_municipalities_latest, "\n")
cat("Spatial GeoPackage latest:", period_incidence_gpkg_latest, "\n")
cat("Map log:", map_log_latest, "\n")
cat("Map captions:", map_captions_latest, "\n\n")

cat("Output folders:\n")
cat("Analysis folder:", incidence_maps_dir, "\n")
cat("Tables:", incidence_maps_tables_dir, "\n")
cat("Maps:", incidence_maps_maps_dir, "\n")
cat("Data:", incidence_maps_data_dir, "\n")
cat("General output tables:", tables_dir, "\n")
cat("General output maps:", maps_dir, "\n")

sink()

# ------------------------------------------------------------
# 14. Console summary
# ------------------------------------------------------------

message("============================================================")
message("PCM incidence period maps completed.")
message("Script version: v1_incidence_period_maps")
message("Period maps created: ", nrow(map_log))
message("Municipality-period rows: ", nrow(period_incidence_dt))
message("Report saved to: ", report_file)
message("============================================================")

print(
  data.table(
    output = c(
      "period_incidence_latest",
      "period_summary_latest",
      "top_municipalities_latest",
      "period_incidence_gpkg_latest",
      "map_log_latest",
      "map_captions_latest",
      "report_file"
    ),
    path = c(
      period_incidence_latest,
      period_summary_latest,
      top_municipalities_latest,
      period_incidence_gpkg_latest,
      map_log_latest,
      map_captions_latest,
      report_file
    )
  )
)
