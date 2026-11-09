# =============================================================================
# R/data_prep.R
# Builds the analytic dataset for the Maine municipal housing-permit analysis.
#
# Source workbook: "2025 Statewide Housing Data.xlsx" (sheet "All").
# The workbook mixes three reporting streams:
#   * Local ............ 182 municipalities that answered the state survey
#   * Local/LUPC ....... 38 LUPC municipalities (permit totals only, no breakdown)
#   * Census ........... 220 municipalities covered only by the Census BP survey
#   * Not Available .... 90 municipalities with no permit data at all
# Building-permit (BP) columns are populated for local respondents, the Census
# columns for Census-covered towns, so the two streams are harmonised below.
#
# Municipal characteristics (median household income, median home value, land
# area/density, tenure, housing-stock composition) come from the 2020-2024 ACS
# 5-year estimates via R/fetch_acs.R -> data/town_covariates.csv.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(janitor)
  library(stringr)
})

prepare_housing_data <- function(path = "2025 Statewide Housing Data.xlsx") {

  raw <- read_excel(path, sheet = "All") |>
    clean_names() |>
    filter(!is.na(town))

  stopifnot(nrow(raw) == 530)   # 529 municipalities + 1 UNORGANIZED TERRITORIES row

  # ---- 1. Recode / tidy ------------------------------------------------------
  towns <- raw |>
    mutate(
      town        = str_to_title(str_trim(town)),
      county      = str_to_title(str_trim(county)),
      data_source = case_when(
        str_detect(data_source, "Local / LUPC") ~ "Local/LUPC",
        TRUE ~ data_source
      ),
      data_source = factor(data_source,
                           levels = c("Local", "Local/LUPC", "Census", "Not Available")),
      required    = required == "Yes",
      # LUPC permits are locally reported as well -> same reporting route
      route       = factor(if_else(data_source %in% c("Local", "Local/LUPC"),
                                   "Local", "Census"),
                           levels = c("Census", "Local")),
      pop         = census_population_24
    ) |>
    # the UNORGANIZED TERRITORIES (OTHER) row aggregates 71 unincorporated
    # places and has no population; it is kept for state totals but excluded
    # from the municipal analysis
    filter(town != "Unorganized Territories (Other)")

  # ---- 2. Harmonise the two permit reporting streams -------------------------
  # BP_Total = ADU + single-family + 2-unit + 3-4-unit + 5+-unit (verified)
  towns <- towns |>
    mutate(
      perm_total = coalesce(bp_total, census_yearly_permits_total),
      perm_sf    = coalesce(bp_single_family_houses, census_single_family_houses),
      perm_adu   = bp_accessory_dwelling_units,
      perm_2     = coalesce(bp_2_unit_buildings,
                            census_two_unit_multifamily_structures),
      perm_34    = coalesce(bp_3_4_unit_buildings,
                            census_three_four_unit_multifamily_structures),
      perm_5p    = coalesce(bp_5_or_more_unit_buildings,
                            census_five_or_more_unit_multifamily_structures),
      # multifamily = 2-unit + 3-4-unit + 5+-unit buildings
      perm_mf    = case_when(
        !is.na(bp_2_unit_buildings) ~ bp_2_unit_buildings +
          bp_3_4_unit_buildings + bp_5_or_more_unit_buildings,
        TRUE ~ census_all_multifamily_structures
      ),
      has_local_bp  = !is.na(bp_total),
      has_census_bp = !is.na(census_yearly_permits_total)
    )

  # ---- 3. Feature engineering ------------------------------------------------
  towns <- towns |>
    mutate(
      log_pop          = if_else(pop > 0, log(pop), NA_real_),
      permits_per_1000 = if_else(pop > 0, 1000 * perm_total / pop, NA_real_),
      mf_share         = if_else(perm_total > 0, perm_mf / perm_total, NA_real_),
      mf_any           = case_when(
        is.na(perm_mf)     ~ NA,
        perm_mf > 0        ~ "Yes",
        TRUE               ~ "No"
      ) |> factor(levels = c("No", "Yes")),
      size_band = cut(pop,
                      breaks = c(0, 500, 1000, 2000, 5000, 10000, Inf),
                      labels = c("<500", "500-1k", "1k-2k", "2k-5k",
                                 "5k-10k", "10k+"),
                      right = FALSE)
    )

  # ---- 3b. ACS municipal characteristics ------------------------------------
  cov_path <- file.path(dirname(path), "data", "town_covariates.csv")
  cov_path <- if (file.exists(cov_path)) cov_path else "data/town_covariates.csv"
  if (file.exists(cov_path)) {
    cov <- read_csv(cov_path, show_col_types = FALSE) |>
      mutate(town = str_to_title(TOWN)) |>
      select(town, inc, home_value, rent, land_sqmi, density, owner_share,
             mf_stock_share, med_year_built, hh_size, hu)
    towns <- towns |> left_join(cov, by = "town")
  } else {
    warning("data/town_covariates.csv not found: run R/fetch_acs.R first")
    towns <- towns |>
      mutate(inc = NA_real_, home_value = NA_real_, rent = NA_real_,
             land_sqmi = NA_real_, density = NA_real_, owner_share = NA_real_,
             mf_stock_share = NA_real_, med_year_built = NA_real_,
             hh_size = NA_real_, hu = NA_real_)
  }

  towns <- towns |>
    mutate(
      log_inc   = if_else(inc > 0, log(inc), NA_real_),
      log_hv    = if_else(home_value > 0, log(home_value), NA_real_),
      log_dens  = if_else(density > 0, log(density), NA_real_),
      log_rate  = if_else(permits_per_1000 > 0, log(permits_per_1000), NA_real_),
      log1p_rate = if_else(permits_per_1000 >= 0,
                           log1p(permits_per_1000), NA_real_),
      built_post_1980 = if_else(!is.na(med_year_built),
                                pmax(0, med_year_built - 1980), NA_real_)
    )

  # ---- 4. Analysis sample ----------------------------------------------------
  # Excluded: 90 municipalities with no permit data (structural missingness)
  # and 1 municipality with zero reported population.
  analytic <- towns |>
    filter(!is.na(perm_total), !is.na(log_pop))

  # Sample used by the covariate models: the analytic sample with a complete
  # set of ACS predictors (income and home value are suppressed by the Census
  # for a handful of very small places).
  modelled <- analytic |>
    filter(!is.na(log_inc), !is.na(log_hv), !is.na(log_dens),
           !is.na(owner_share), !is.na(mf_stock_share), !is.na(med_year_built))

  # Towns observed under *both* reporting streams -> measurement comparison
  overlap <- towns |>
    filter(has_local_bp, has_census_bp, !is.na(pop), pop > 0) |>
    mutate(
      ratio      = bp_total / census_yearly_permits_total,
      diff_units = bp_total - census_yearly_permits_total
    )

  list(all = towns, data = analytic, modelled = modelled, overlap = overlap)
}
