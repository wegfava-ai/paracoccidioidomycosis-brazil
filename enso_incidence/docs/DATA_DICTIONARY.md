# Data dictionary

## `enso_pcm_annual_national_dataset_latest.csv`

- `year`: calendar year.
- `cases`: annual number of probable acute/subacute PCM records.
- `population`: annual national population denominator after correction using the broad IBGE municipal population table.
- `municipalities`: number of municipalities in the annual grid.
- `event_municipality_years`: number of municipality-year rows with at least one PCM record.
- `missing_population_rows`: number of municipality-year rows without valid population denominator.
- `incidence_100k`: annual national incidence per 100,000 inhabitants.
- `oni_current`: annual mean Oceanic Niño Index in the same year.
- `oni_lag1y`, `oni_lag2y`, `oni_lag3y`: lagged annual mean ONI variables.
- `enso_phase_current`: ENSO phase based on current annual mean ONI.

## `pcm_municipal_annualized_incidence_by_period_latest.csv`

- `municipality_code_6`: six-digit municipality code used in the analytical grid.
- `cases`: total probable acute/subacute PCM records in the period.
- `person_years`: summed population denominator across years in the period.
- `annualized_incidence_100k`: cases divided by person-years multiplied by 100,000.
- `period_label`: period label.
- `period_id`: machine-readable period ID.
- `start_year`: first year of the period.
- `end_year`: final year of the period.

## `pcm_incidence_period_summary_latest.csv`

- Summary of total cases, person-years, and national annualized incidence by period.

## `enso_pcm_national_effect_terms_latest.csv`

- Regression model effect estimates.
- `irr`: incidence rate ratio.
- `conf_low_95` and `conf_high_95`: approximate 95% confidence interval.
- `p_value`: model-based p-value.

## Important interpretation note

All outputs are based on aggregate records. Incidence estimates are proxies based on administrative health records and should be interpreted in light of underdiagnosis, underreporting, and access-to-care differences.
