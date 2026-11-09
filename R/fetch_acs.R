# =============================================================================
# R/fetch_acs.R
# Builds data/town_covariates.csv: municipal characteristics for every Maine
# county subdivision, from the 2020-2024 ACS 5-year estimates plus the Census
# Gazetteer (land area -> population density).
#
# Why not the API? api.census.gov now requires an API key for every request, so
# the tables are read from the public ACS table-based summary file instead:
#   https://www2.census.gov/programs-surveys/acs/summary_file/2024/table-based-SF/
# Each table is a pipe-delimited file with one row per geography (GEO_ID, E, M).
#
# Variables
#   B01003 total population          B25001 housing units
#   B19013 median household income   B25077 median value, owner-occupied units
#   B25064 median gross rent         B25003 tenure (owner / renter)
#   B25024 units in structure        B25035 median year structure built
#
# Run once; the script is idempotent and skips files already on disk.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(readxl)
})

ACS_YEAR <- 2024          # 2020-2024 5-year estimates
ACS_BASE <- sprintf(
  "https://www2.census.gov/programs-surveys/acs/summary_file/%1$d/table-based-SF/data/5YRData/",
  ACS_YEAR
)
GAZ_URL <- sprintf(
  "https://www2.census.gov/geo/docs/maps-data/data/gazetteer/%1$d_Gazetteer/%1$d_gaz_cousubs_23.txt",
  ACS_YEAR
)

# FIPS -> county name, so that the join can use (name, county) and so that
# duplicated subdivision names in different counties cannot collide
fips_county <- c(
  "001" = "Androscoggin", "003" = "Aroostook", "005" = "Cumberland",
  "007" = "Franklin", "009" = "Hancock", "011" = "Kennebec",
  "013" = "Knox", "015" = "Lincoln", "017" = "Oxford", "019" = "Penobscot",
  "021" = "Piscataquis", "023" = "Sagadahoc", "025" = "Somerset",
  "027" = "Waldo", "029" = "Washington", "031" = "York"
)

# Municipal suffixes are dropped before matching, and these two names lose
# their meaning when the suffix is stripped ("Old Town", "Monhegan Island")
suffix_rx <- "\\s+(CITY|TOWN|PLANTATION|UT|GORE|RESERVATION|TOWNSHIP)$"
clean_name <- function(x) {
  x <- str_replace_all(toupper(str_trim(x)), "[.-]", " ")
  x <- str_squish(x)
  x <- str_replace_all(x, suffix_rx, "")
  str_squish(x)
}
name_fix <- c("OLD" = "OLD TOWN", "MONHEGAN" = "MONHEGAN ISLAND")

# ---- 1. download helpers ----------------------------------------------------
fetch <- function(url, dest) {
  dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
  if (!file.exists(dest) || file.size(dest) == 0) {
    message("downloading ", basename(dest))
    utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  }
  dest
}

acs_table <- function(table, cols) {
  fn <- sprintf("acsdt5y%d-%s.dat", ACS_YEAR, table)
  dest <- fetch(paste0(ACS_BASE, fn), file.path("data", "acs", fn))
  read_delim(dest, delim = "|", show_col_types = FALSE, progress = FALSE) |>
    rename_with(str_trim) |>
    filter(str_starts(GEO_ID, "0600000US23")) |>
    mutate(GEOID = str_remove(GEO_ID, "^0600000US"), .after = GEO_ID) |>
    select(GEOID, all_of(cols)) |>
    mutate(across(all_of(cols), ~ {
      v <- suppressWarnings(as.numeric(.x))
      v[v <= -100000] <- NA          # ACS annotation codes (-666666666 etc.)
      v
    }))
}

# ---- 2. geography ------------------------------------------------------------
gaz <- read_tsv(fetch(GAZ_URL, file.path("data", "acs", "gaz_cousubs_24.txt")),
                show_col_types = FALSE, progress = FALSE) |>
  filter(USPS == "ME") |>
  transmute(
    GEOID      = str_trim(GEOID),
    acs_name   = str_trim(NAME),
    land_sqmi  = as.numeric(ALAND_SQMI),
    county     = unname(fips_county[str_sub(GEOID, 3, 5)]),
    key        = clean_name(acs_name)
  ) |>
  mutate(key = coalesce(unname(name_fix[key]), key))

# ---- 3. ACS tables ----------------------------------------------------------
pop  <- acs_table("b01003", "B01003_E001")
hu   <- acs_table("b25001", "B25001_E001")
inc  <- acs_table("b19013", "B19013_E001")
val  <- acs_table("b25077", "B25077_E001")
rent <- acs_table("b25064", "B25064_E001")
ten  <- acs_table("b25003", c("B25003_E001", "B25003_E002", "B25003_E003"))
strc <- acs_table("b25024", sprintf("B25024_E%03d", 1:9))
yr   <- acs_table("b25035", "B25035_E001")

covariates <- gaz |>
  left_join(pop,  by = "GEOID") |>
  left_join(hu,   by = "GEOID") |>
  left_join(inc,  by = "GEOID") |>
  left_join(val,  by = "GEOID") |>
  left_join(rent, by = "GEOID") |>
  left_join(ten,  by = "GEOID") |>
  left_join(strc, by = "GEOID") |>
  left_join(yr,   by = "GEOID") |>
  transmute(
    key, GEOID, acs_name, county,
    land_sqmi,
    pop_acs       = B01003_E001,
    hu            = B25001_E001,
    inc           = B19013_E001,
    home_value    = B25077_E001,
    rent          = B25064_E001,
    occ           = B25003_E001,
    owner_share   = B25003_E002 / B25003_E001,
    mf_stock_share = rowSums(across(sprintf("B25024_E%03d", 3:9)), na.rm = FALSE) /
      B25024_E001,
    med_year_built = B25035_E001,
    density       = pop_acs / land_sqmi,
    hh_size       = if_else(occ > 0, pop_acs / occ, NA_real_)
  )

# ---- 4. attach to the workbook's municipality list --------------------------
towns <- read_excel("2025 Statewide Housing Data.xlsx", sheet = "All") |>
  janitor::clean_names() |>
  filter(!is.na(town)) |>
  transmute(
    TOWN       = str_trim(town),
    COUNTY     = str_to_title(str_trim(county)),
    key        = coalesce(unname(name_fix[clean_name(TOWN)]), clean_name(TOWN)),
    county_key = str_remove_all(toupper(COUNTY), "[^A-Z]")
  ) |>
  distinct(TOWN, .keep_all = TRUE)

acs_lookup <- covariates |>
  mutate(county_key = str_remove_all(toupper(county), "[^A-Z]")) |>
  select(-county) |>
  # reservation land is tabulated once per county; one row per (name, county)
  distinct(key, county_key, .keep_all = TRUE)

out <- towns |>
  left_join(acs_lookup, by = c("key", "county_key")) |>
  distinct(TOWN, .keep_all = TRUE) |>
  select(TOWN, COUNTY, GEOID, land_sqmi, pop_acs, hu, inc, home_value, rent,
         owner_share, mf_stock_share, med_year_built, density, hh_size)

dir.create("data", showWarnings = FALSE)
write_csv(out, "data/town_covariates.csv")
message(sprintf("wrote data/town_covariates.csv: %d towns, %d matched to ACS",
                nrow(out), sum(!is.na(out$GEOID))))
invisible(out)
