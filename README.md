# Ecoepidemiology of probable acute/subacute paracoccidioidomycosis in Brazil

## Overview

This repository folder contains tables, figures, maps, captions, and manuscript text generated from a validated analytical workflow. The analysis modeled the spatial distribution and future climatic suitability of **probable acute/subacute paracoccidioidomycosis (PCM)** in Brazil using municipality-level administrative health records, environmental predictors, and future climate projections.

**Important interpretation note:** all current and future values shown here are **model-predicted annualized suitability rates per 100,000 inhabitants**. They should not be interpreted as observed incidence or deterministic forecasts of future case counts.

## Key outputs

- Retrospective modeling period: **1998-2024**.
- Main projection model: **M3 climate-soil-elevation negative binomial model**.
- Complementary full retrospective model: **M6 full retrospective negative binomial model including ENSO**.
- Main future scenarios: **SSP245** and **SSP585**, periods **2041-2060** and **2061-2080**.
- Highest median future suitability scenario: **SSP585 | 2061-2080**.
- Current median predicted suitability: **0.0320 per 100,000**.
- Median predicted suitability under SSP585 | 2061-2080: **0.0797 per 100,000**.
- Median percent change under SSP585 | 2061-2080: **119.8%**.
- Top regional increase under SSP585 | 2061-2080: **North; change=0.1382**.
- Top state-level increase under SSP585 | 2061-2080: **PI; change=0.2570**.
- Municipalities with very large increase under SSP585 | 2061-2080: **3,424**.

## Folder structure

```text
github_publication_package/
├── README.md
├── tables/              # publication-ready CSV tables
├── figures/             # publication-ready figures, PNG and PDF
├── maps/                # publication-ready maps, PNG and PDF
├── captions/            # figure, table, and map captions
├── manuscript_text/     # editable Methods, Results, limitations, and callouts
└── data/                # file manifests and compact key-result tables
```

## Main modeling results

The BIC-supported model used for current and future suitability projections was the M3 climate-soil-elevation negative binomial model. The AIC-supported full retrospective model was M6, which included ENSO and additional predictor groups. M3 was used for projection because its predictors were suitable for static spatial and future climate analyses, while M6 included temporal ENSO predictors that were not projected under future climate scenarios.

### Current and future median suitability

<p align="center">
  <img src="figures/figure_1_current_future_median_suitability_latest.png" alt="Current and future median suitability" width="850">
</p>

The figure above compares the current baseline with future projections across climate scenario-periods. The median suitability increased progressively, with the largest increase under **SSP585 | 2061-2080**.

### Future-minus-current suitability

<p align="center">
  <img src="figures/figure_2_median_future_minus_current_latest.png" alt="Median future-minus-current suitability" width="850">
</p>

Positive values indicate higher future model-predicted suitability compared with the current baseline.

### Regional trajectories

<p align="center">
  <img src="figures/figure_3_regional_change_across_scenarios_latest.png" alt="Regional change across scenarios" width="850">
</p>

Regional trajectories show that increases were not spatially uniform. The largest median increase under **SSP585 | 2061-2080** was observed in **North; change=0.1382**.

### Municipality-level change categories

<p align="center">
  <img src="figures/figure_4_change_categories_by_scenario_latest.png" alt="Municipality-level suitability change categories" width="850">
</p>

Categories summarize percent change relative to the current baseline. Under **SSP585 | 2061-2080**, **3,424** municipalities were classified as having a very large increase in modeled suitability.

### State-level change under the most severe scenario

<p align="center">
  <img src="figures/figure_5_uf_median_change_severe_scenario_latest.png" alt="State-level median change under the severe scenario" width="850">
</p>

The state-level summary highlights federative units with the largest median future-minus-current suitability under **SSP585 | 2061-2080**.

### Priority municipalities

<p align="center">
  <img src="figures/figure_6_priority_municipalities_severe_scenario_latest.png" alt="Priority municipalities under the severe scenario" width="850">
</p>

Priority municipalities were ranked using a combined score integrating high future predicted suitability and high absolute increase relative to the current baseline.

## Publication-ready maps

### Current predicted suitability

<p align="center">
  <img src="maps/map_1_current_suitability_predicted_annual_rate_100k.png" alt="Current predicted suitability map" width="850">
</p>

This map shows current model-predicted suitability for probable acute/subacute PCM in Brazil using the M3 climate-soil-elevation negative binomial model.

### Future predicted suitability maps

#### SSP245 | 2041-2060

<p align="center">
  <img src="maps/map_2_future_suitability_predicted_annual_rate_100k_SSP245_2041_2060.png" alt="Future predicted suitability SSP245 2041-2060" width="850">
</p>

#### SSP245 | 2061-2080

<p align="center">
  <img src="maps/map_3_future_suitability_predicted_annual_rate_100k_SSP245_2061_2080.png" alt="Future predicted suitability SSP245 2061-2080" width="850">
</p>

#### SSP585 | 2041-2060

<p align="center">
  <img src="maps/map_4_future_suitability_predicted_annual_rate_100k_SSP585_2041_2060.png" alt="Future predicted suitability SSP585 2041-2060" width="850">
</p>

#### SSP585 | 2061-2080

<p align="center">
  <img src="maps/map_5_future_suitability_predicted_annual_rate_100k_SSP585_2061_2080.png" alt="Future predicted suitability SSP585 2061-2080" width="850">
</p>

### Future-minus-current change under SSP585 | 2061-2080

<p align="center">
  <img src="maps/map_change_absolute_future_minus_current_rate_100k_SSP585_2061_2080.png" alt="Future-minus-current suitability SSP585 2061-2080" width="850">
</p>

Positive values indicate higher future suitability compared with the current baseline.

### Suitability change categories under SSP585 | 2061-2080

<p align="center">
  <img src="maps/map_change_categories_SSP585_2061_2080.png" alt="Suitability change categories SSP585 2061-2080" width="850">
</p>

This map shows municipality-level categories of future suitability change based on percent change from current model-predicted suitability.

## Summary tables

### Current and future suitability summary

| scenario | municipalities | complete_predictions | median_rate_100k | median_future_minus_current_100k | median_percent_change |
| --- | --- | --- | --- | --- | --- |
| Current baseline | 5570 | 5356 | 0.032 | 0 | 0 |
| SSP245 \| 2041-2060 | 5570 | 5285 | 0.0523 | 0.0178 | 52.8466 |
| SSP245 \| 2061-2080 | 5570 | 5285 | 0.0566 | 0.0214 | 64.8897 |
| SSP585 \| 2041-2060 | 5570 | 5285 | 0.0603 | 0.0255 | 71.7926 |
| SSP585 \| 2061-2080 | 5570 | 5285 | 0.0797 | 0.0438 | 119.7828 |

### Top state-level increases under SSP585 | 2061-2080

| uf | region | median_future_minus_current_100k | median_percent_change |
| --- | --- | --- | --- |
| PI | Northeast | 0.257 | 333.6907 |
| RR | North | 0.2474 | 303.6722 |
| TO | North | 0.2289 | 269.2721 |
| MT | Central-West | 0.1897 | 205.9339 |
| RO | North | 0.1698 | 405.8602 |
| DF | Central-West | 0.1154 | 183.7921 |
| PA | North | 0.1019 | 211.4294 |
| MG | Southeast | 0.1017 | 135.6854 |
| MA | Northeast | 0.0746 | 129.8662 |
| GO | Central-West | 0.0615 | 225.4306 |

### Priority municipalities under SSP585 | 2061-2080

| municipality | uf | future_rate_100k | future_minus_current_100k | percent_change | priority_score |
| --- | --- | --- | --- | --- | --- |
| Itaúba | MT | 1.6929 | 1.6009 | 1739.8917 | 2 |
| Ipiranga do Norte | MT | 1.2102 | 1.094 | 941.8439 | 5 |
| Colíder | MT | 1.2885 | 1.0029 | 351.145 | 7 |
| Tabaporã | MT | 1.159 | 1.03 | 798.6876 | 7 |
| Cláudia | MT | 1.1405 | 1.0037 | 733.56 | 9 |
| Porto Alegre do Norte | MT | 1.0381 | 0.9762 | 1577.3884 | 13 |
| Nova Santa Helena | MT | 1.0229 | 0.9582 | 1481.474 | 16 |
| Novo Horizonte do Norte | MT | 1.0257 | 0.9229 | 897.4682 | 17 |
| Novo Acordo | TO | 1.0011 | 0.9274 | 1257.5806 | 18 |
| Marcelândia | MT | 0.9978 | 0.9119 | 1062.1148 | 21 |
| Terra Nova do Norte | MT | 1.0524 | 0.6798 | 182.4303 | 23 |
| Porto dos Gaúchos | MT | 0.937 | 0.8337 | 807.3191 | 24 |
| Santa Terezinha | MT | 0.9417 | 0.8313 | 752.9449 | 24 |
| Sinop | MT | 0.908 | 0.8063 | 792.5768 | 27 |
| Lajeado | TO | 0.7796 | 0.7357 | 1675.5958 | 33 |

## Available files

### Tables

- `tables/table_1_retrospective_model_selection_latest.csv`
- `tables/table_2_retrospective_prediction_performance_latest.csv`
- `tables/table_3_predictor_group_effects_latest.csv`
- `tables/table_4_current_future_overall_summary_latest.csv`
- `tables/table_5_regional_future_change_latest.csv`
- `tables/table_6_uf_future_change_severe_scenario_latest.csv`
- `tables/table_7_future_change_categories_latest.csv`
- `tables/table_8_priority_municipalities_severe_scenario_latest.csv`
- `tables/table_9_top_future_rate_municipalities_severe_scenario_latest.csv`
- `tables/table_10_top_absolute_increase_municipalities_severe_scenario_latest.csv`

### Manuscript text

- `manuscript_text/methods_results_text_latest.txt`
- `manuscript_text/methods_section_latest.txt`
- `manuscript_text/results_section_latest.txt`
- `manuscript_text/limitations_paragraph_latest.txt`
- `manuscript_text/suggested_table_figure_map_callouts_latest.txt`

### Captions

- `captions/publication_figure_and_table_captions_latest.txt`
- `captions/publication_map_captions_latest.txt`
- `captions/publication_captions_latest.csv`

## Methodological notes

- The case definition represents **probable acute/subacute PCM**, not clinically confirmed acute/subacute PCM.
- The age restriction of **≤20 years** was used as a conservative proxy for the acute/subacute clinical form.
- Soil predictors were extracted from SoilGrids using municipality centroids, not polygon-level zonal means.
- Future projections held soil and elevation constant and replaced current BIOCLIM predictors with CMIP6 ensemble future predictors.
- Values are model-predicted suitability rates per 100,000 inhabitants and should not be interpreted as observed incidence.
- M6 included temporal ENSO predictors and was retained as a full retrospective model, but M3 was used for future projections.

## Reproducibility

This folder was generated automatically by:

```text
23_build_publication_output_package_and_readme.R
```

The package includes only publication-ready outputs. The complete analytical workflow should be documented separately with the numbered R scripts used to generate the source datasets, models, predictions, figures, and maps.

## Citation and reuse

This repository is intended to support transparent review, reproducibility, and publication preparation. Please cite the associated manuscript when available.

