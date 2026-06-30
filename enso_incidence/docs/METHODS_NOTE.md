# Methods note

## Outcome definition

Probable acute/subacute paracoccidioidomycosis was defined using DATASUS records with ICD-10 code B41 and age <= 20 years at the event. This is an epidemiological proxy and should not be interpreted as confirmed clinical form.

## ENSO exposure

Annual mean Oceanic Niño Index (ONI) was used as the main continuous ENSO variable. Annual ENSO phase was defined descriptively using thresholds: El Niño >= 0.5, La Niña <= -0.5, and Neutral between -0.5 and 0.5.

## National annual ENSO analysis

Annual national PCM counts were aggregated across municipalities. Models used a population offset and linear calendar-year term. Because ENSO is a year-level exposure and the time series has 27 years, regression results should be interpreted as exploratory.

## Municipal incidence maps

Annualized municipal incidence was calculated as total cases in the period divided by summed municipal population person-years in the same period, multiplied by 100,000.

Periods used:

- 1998-2024 overall
- 1998-2004
- 2005-2011
- 2012-2018
- 2019-2024

## Population denominator correction

Population denominators were updated using the widest valid municipal population file found in the project environment. The validated file was `D:/DATASUS/data_raw/population/ibge_official_municipal_population_sources.rds`, reducing missing population rows from 85,791 to 60 in the national municipality-year grid.

## Mapping note

Maps were generated without municipal borders to avoid visual artifacts. State and national borders were drawn from a light local state boundary layer when available.
