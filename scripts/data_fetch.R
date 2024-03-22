library(tidyverse)
library(ggplot2)
library(sf)
library(leaflet)
library(tidycensus)
library(ipumsr)

city_list <- list("Dallas, TX", "San Diego, CA", "Austin, TX", "Charlotte, NC", "Washington, DC", "Phoenix, AZ", "San Antonio, TX", "Atlanta, GA", "Houston, TX", "Miami, FL", "Philadelphia, PA")

# city_geometry <- tigris::places(year = 2022, cb = TRUE) %>%
#   left_join(get_acs(
#     "place",
#     year = 2022,
#     variables = c("pop" = "B01001_001", 
#                   "pop_black" = "B01001B_001", 
#                   "pop_native" = "B01001C_001", 
#                   "pop_asian" = "B01001D_001", 
#                   "pop_pacific" = "B01001E_001", 
#                   "pop_otherrace" = "B01001F_001", 
#                   "pop_multirace" = "B01001G_001", 
#                   "pop_white" = "B01001H_001", 
#                   "pop_latino" = "B01001I_001"),
#     output = "wide"
#   ), by = "GEOID", suffix = c("", "_acs"))
# 
# dallas_remove <- sf::st_sfc(sf::st_polygon(list(matrix(c(
#   -96.6, 32.75,
#   -96.6, 33.1,
#   -96.3, 33.1,
#   -96.3, 32.75,
#   -96.6, 32.75
# ), ncol = 2, byrow = TRUE))), crs = sf::st_crs(city_geometry))
# 
# target_cities <- city_geometry %>%
#   mutate(cityst = paste0(NAME, ", ", STUSPS)) %>%
#   filter(cityst %in% city_list) %>%
#   mutate(geometry = if_else(cityst == "Dallas, TX", 
#                             sf::st_difference(geometry, dallas_remove), 
#                             geometry)) %>%
#   mutate(geometry = sf::st_simplify(geometry, dTolerance = 1))
# st_write(target_cities, "data/target_cities.geojson")
target_cities <- st_read("data/target_cities.geojson")

# city_plots <- lapply(target_cities$cityst, function(city) {
#   city_data <- target_cities %>% filter(cityst == city)
#   ggplot(city_data) + 
#     geom_sf(fill = "grey10", color = NA) +
#     theme_void()
# })

fips_codes <- tidycensus::fips_codes %>%
  mutate(CountyFIPS = paste0(state_code, county_code))

# target_counties <- tigris::counties(year = 2021) %>%
#   sf::st_join(target_cities, sf::st_intersects) %>%
#   filter(!is.na(GEOID.y)) %>%
#   rename(CountyFIPS = GEOID.x) %>%
#   left_join(fips_codes, "CountyFIPS")
# st_write(target_counties, "data/target_counties.geojson")
target_counties <- st_read("data/target_counties.geojson")

# target_metros <- tigris::core_based_statistical_areas(year = 2021) %>%
#   sf::st_join(target_cities, sf::st_intersects, suffix = c("", "_city")) %>%
#   filter(!is.na(GEOID_city)) %>%
#   left_join(get_acs(
#     "cbsa",
#     year = 2022,
#     variables = c("pop" = "B01001_001", 
#                   "pop_black" = "B01001B_001", 
#                   "pop_native" = "B01001C_001", 
#                   "pop_asian" = "B01001D_001", 
#                   "pop_pacific" = "B01001E_001", 
#                   "pop_otherrace" = "B01001F_001", 
#                   "pop_multirace" = "B01001G_001", 
#                   "pop_white" = "B01001H_001", 
#                   "pop_latino" = "B01001I_001"),
#     output = "wide"
#   ), by = "GEOID", suffix = c("", "_acs"))
# st_write(target_metros, "data/target_metros.geojson")
target_metros <- st_read("data/target_metros.geojson")

# target_tracts <- unique(target_counties$state) %>%
#   map(~ tigris::tracts(state = .x, year = 2021) %>% 
#           sf::st_join(target_cities, sf::st_intersects) %>%
#           filter(!is.na(GEOID.y))) %>%
#   bind_rows() %>%
#   rename(GEOID = GEOID.x) %>%
#   mutate(TractFIPS = GEOID) %>%
#   distinct(TractFIPS, .keep_all=TRUE) %>%
#   left_join(get_acs(
#     "tract",
#     year = 2022,
#     state = unique(target_counties$state),
#     variables = c("pop" = "B01001_001", 
#                   "pop_black" = "B01001B_001", 
#                   "pop_native" = "B01001C_001", 
#                   "pop_asian" = "B01001D_001", 
#                   "pop_pacific" = "B01001E_001", 
#                   "pop_otherrace" = "B01001F_001", 
#                   "pop_multirace" = "B01001G_001", 
#                   "pop_white" = "B01001H_001", 
#                   "pop_latino" = "B01001I_001"),
#     output = "wide"
#   ), by = "GEOID", suffix = c("_city", ""))
# st_write(target_tracts, "data/target_tracts.geojson")
target_tracts <- st_read("data/target_tracts.geojson")
dallas_tracts <- target_tracts %>%
  filter(NAME.y == "Dallas")

# target_2019tracts <- unique(target_counties$state) %>%
#   map(~ tigris::tracts(state = .x, year = 2019) %>% 
#         sf::st_join(target_cities, sf::st_intersects) %>%
#         filter(!is.na(GEOID.y))) %>%
#   bind_rows() %>%
#   rename(GEOID = GEOID.x) %>%
#   mutate(TractFIPS = GEOID) %>%
#   distinct(TractFIPS, .keep_all=TRUE) %>%
#   left_join(get_acs(
#     "tract",
#     year = 2019,
#     state = unique(target_counties$state),
#     variables = c("pop" = "B01001_001", 
#                   "pop_black" = "B01001B_001", 
#                   "pop_native" = "B01001C_001", 
#                   "pop_asian" = "B01001D_001", 
#                   "pop_pacific" = "B01001E_001", 
#                   "pop_otherrace" = "B01001F_001", 
#                   "pop_multirace" = "B01001G_001", 
#                   "pop_white" = "B01001H_001", 
#                   "pop_latino" = "B01001I_001"),
#     output = "wide"
#   ), by = "GEOID", suffix = c("_city", ""))
# st_write(target_2019tracts, "data/target_2019tracts.geojson")
target_2019tracts <- st_read("data/target_2019tracts.geojson")

# target_PUMAs <- unique(target_counties$state) %>%
#   map(~ tigris::pumas(state = .x, year = 2022) %>% 
#           sf::st_join(target_cities, sf::st_intersects) %>%
#           filter(!is.na(GEOID))) %>%
#   bind_rows() %>%
#   rename(GEOID.y = GEOID,
#          GEOID = GEOID20) %>%
#   mutate(PUMAFIPS = GEOID) %>%
#   distinct(PUMAFIPS, .keep_all=TRUE) %>%
#   left_join(get_acs(
#     "puma",
#     year = 2022,
#     state = unique(target_counties$state),
#     variables = c("pop" = "B01001_001", 
#                   "pop_black" = "B01001B_001", 
#                   "pop_native" = "B01001C_001", 
#                   "pop_asian" = "B01001D_001", 
#                   "pop_pacific" = "B01001E_001", 
#                   "pop_otherrace" = "B01001F_001", 
#                   "pop_multirace" = "B01001G_001", 
#                   "pop_white" = "B01001H_001", 
#                   "pop_latino" = "B01001I_001"),
#     output = "wide"
#   ), by = "GEOID", suffix = c("_city", ""))
# st_write(target_PUMAs, "data/target_PUMAs.geojson")
target_PUMAs <- st_read("data/target_PUMAs.geojson")

# target_sds_TX <- tigris::school_districts(state = "TX") %>% 
#   filter(apply(sf::st_intersects(geometry, sf::st_geometry(target_cities), sparse = FALSE), 1, any)) %>%
#   mutate(
#     sd_area = sf::st_area(geometry)
#   ) %>%
#   sf::st_intersection(target_cities) %>%
#   sf::st_make_valid() %>%
#   mutate(intersection_area = sf::st_area(geometry),
#          city_area = sf::st_area(sf::st_geometry(target_cities[match(NAME.1, target_cities$NAME),])),
#          city_prop = as.numeric((intersection_area / city_area)),
#          sd_prop = as.numeric((intersection_area / sd_area))) %>%
#   filter(sd_prop > 2/3) %>%
#   mutate(NAME = stringr::str_replace(NAME, "Independent School District", "ISD"),
#          DNAME = str_to_upper(NAME))
# st_write(target_sds_TX, "data/target_sds_TX.geojson")
target_sds_TX <- st_read("data/target_sds_TX.geojson")

# Data import
races = c("all", "black", "native", "aapi", "add", "white", "latino")
acs_races = c("black", "native", "asian", "pacific", "native", "otherrace", "multirace", "white", "latino")
keeps = c("GEOID", races)
threshold = 0.5
sum_columns <- function(df, races, max_cols) {
  for (race in races) {
    for (suffix in c("_E", "_M")) {
      column_names <- paste0(seq_len(max_cols), race, suffix)
      df <- df %>%
        mutate(!!sym(paste0(race, suffix)) := rowSums(select(., all_of(column_names)), na.rm = TRUE))
    }
  }
  return(df)
}

sum_columns2 <- function(df, races, max_cols) {
  for (race in races) {
    for (suffix in c("E", "M", "_E", "_M")) {
      column_names <- paste0(seq_len(max_cols), race, suffix)
      df <- df %>%
        mutate(!!sym(paste0(race, suffix)) := rowSums(select(., all_of(column_names)), na.rm = TRUE))
    }
  }
  return(df)
}

pums <- read_ipums_micro(read_ipums_ddi("data/IPUMS/usa_00003.xml"), data_file = "data/IPUMS/usa_00003.dat") %>%
  mutate(cityst = as_factor(CITY, levels = "labels"),
         stpuma = paste0(STATEFIP, sprintf("%05d", PUMA)),
         GEOID = stpuma) %>%
  filter(STATEFIP %in% target_counties$STATEFP.x, YEAR == 2022) %>%
  mutate(
    RACEON = case_match(
      RACE,
      1 ~ "white",
      2 ~ "black",
      3 ~ "native",
      4:6 ~ "aapi",
      7:9 ~ "add",
    ),
    hisp = if_else(HISPAN %in% c(1,2,3,4), "latino", "!latino")
  )

################################################################################
### BEGIN ###

join_list <- list()

HD_HQA_overcrowded_housing_place <- get_acs(
      "place",
      variables = c(
        "all" = "B25014_001",
        "black" = "B25014B_001",
        "native" = "B25014C_001",
        "asian" = "B25014D_001",
        "pacific" = "B25014E_001",
        "otherrace" = "B25014F_001",
        "multirace" = "B25014G_001",
        "white" = "B25014H_001",
        "latino" = "B25014I_001",
        "all_" = c("B25014_005", "B25014_006", "B25014_007", "B25014_011", "B25014_012", "B25014_013"),
        "black_" = "B25014B_003",
        "native_" = "B25014C_003",
        "asian_" = "B25014D_003",
        "pacific_" = "B25014E_003",
        "otherrace_" = "B25014F_003",
        "multirace_" = "B25014G_003",
        "white_" = "B25014H_003",
        "latino_" = "B25014I_003"
      ),
      year = 2022,
      output = "wide"
    ) %>%
    filter(GEOID %in% target_cities$GEOID) %>%
    mutate(all = case_when(allM / allE > threshold ~ NA,
                           (all_1M+all_2M+all_3M+all_4M+all_5M+all_6M) / (all_1E+all_2E+all_3E+all_4E+all_5E+all_6E) > threshold ~ NA,
                           TRUE ~ (all_1E+all_2E+all_3E+all_4E+all_5E+all_6E) / allE),
           black = case_when(blackM / blackE > threshold ~ NA,
                             black_M / black_E > threshold ~ NA,
                             TRUE ~ black_E / blackE),
           native = case_when(nativeM / nativeE > threshold ~ NA,
                             native_M / native_E > threshold ~ NA,
                             TRUE ~ native_E / nativeE),
           aapi = case_when( (asianM + pacificM) / (asianE + pacificE) > threshold ~ NA,
                             (asian_M + pacific_M) / (asian_E + pacific_E) > threshold ~ NA,
                             TRUE ~ (asian_E + pacific_E) / (asianE + pacificE)),
           add = case_when(  (otherraceM+multiraceM) / (otherraceE+multiraceE) > threshold ~ NA,
                             (otherrace_M+multirace_M) / (otherrace_E+multirace_E) > threshold ~ NA,
                             TRUE ~ (otherrace_E+multirace_E) / (otherraceE+multiraceE)),
           white = case_when(whiteM / whiteE > threshold ~ NA,
                             white_M / white_E > threshold ~ NA,
                             TRUE ~ white_E / whiteE),
           latino = case_when(latinoM / latinoE > threshold ~ NA,
                             latino_M / latino_E > threshold ~ NA,
                             TRUE ~ latino_E / latinoE)) %>%
    select(all_of(keeps)) %>%
    mutate(City = target_cities$cityst[match(.data$GEOID, target_cities$GEOID)],
           Geography = "Place",
           System = "Human Development",
           Area = "Housing Quality and Affordability",
           Metric = "Overcrowded Housing (More than 1 occupant per room)",
           Measure = "%",
           Source = "ACS 5-year",
           Year = 2022,
           BiggerBetter = 0,
           Name = "Overcrowded housing",
           Id = "overcrowded-housing"
    ); join_list <- unique(c(join_list, list(HD_HQA_overcrowded_housing_place)))


HD_HQA_overcrowded_housing_tract <- get_acs(
  "tract",
  variables = c(
    "all" = "B25014_001",
    "black" = "B25014B_001",
    "native" = "B25014C_001",
    "asian" = "B25014D_001",
    "pacific" = "B25014E_001",
    "otherrace" = "B25014F_001",
    "multirace" = "B25014G_001",
    "white" = "B25014H_001",
    "latino" = "B25014I_001",
    "all_" = c("B25014_005", "B25014_006", "B25014_007", "B25014_011", "B25014_012", "B25014_013"),
    "black_" = "B25014B_003",
    "native_" = "B25014C_003",
    "asian_" = "B25014D_003",
    "pacific_" = "B25014E_003",
    "otherrace_" = "B25014F_003",
    "multirace_" = "B25014G_003",
    "white_" = "B25014H_003",
    "latino_" = "B25014I_003"
  ),
  state = unique(target_tracts$STATE_NAME),
  year = 2022,
  output = "wide"
) %>%
  filter(GEOID %in% target_tracts$TractFIPS) %>%
  mutate(all = case_when(
          # allM / allE > threshold ~ NA,
          # (all_1M+all_2M+all_3M+all_4M+all_5M+all_6M) / (all_1E+all_2E+all_3E+all_4E+all_5E+all_6E) > threshold ~ NA,
          TRUE ~ (all_1E+all_2E+all_3E+all_4E+all_5E+all_6E) / allE),
         black = case_when(
           # blackM / blackE > threshold ~ NA,
           # black_M / black_E > threshold ~ NA,
           TRUE ~ black_E / blackE),
         native = case_when(
           # nativeM / nativeE > threshold ~ NA,
           # native_M / native_E > threshold ~ NA,
           TRUE ~ native_E / nativeE),
         aapi = case_when( 
           # (asianM + pacificM) / (asianE + pacificE) > threshold ~ NA,
           # (asian_M + pacific_M) / (asian_E + pacific_E) > threshold ~ NA,
           TRUE ~ (asian_E + pacific_E) / (asianE + pacificE)),
         add = case_when(  
           # (otherraceM+multiraceM) / (otherraceE+multiraceE) > threshold ~ NA,
           # (otherrace_M+multirace_M) / (otherrace_E+multirace_E) > threshold ~ NA,
           TRUE ~ (otherrace_E+multirace_E) / (otherraceE+multiraceE)),
         white = case_when(
           # whiteM / whiteE > threshold ~ NA,
           # white_M / white_E > threshold ~ NA,
           TRUE ~ white_E / whiteE),
         latino = case_when(
           # latinoM / latinoE > threshold ~ NA,
           # latino_M / latino_E > threshold ~ NA,
           TRUE ~ latino_E / latinoE)) %>%
  select(all_of(keeps)) %>%
  mutate(City = target_tracts$cityst[match(.data$GEOID, target_tracts$GEOID)],
         Geography = "Tract",
         System = "Human Development",
         Area = "Housing Quality and Affordability",
         Metric = "Overcrowded Housing (More than 1 occupant per room)",
         Measure = "%",
         Source = "ACS 5-year",
         Year = 2022,
         BiggerBetter = 0,
         Name = "Overcrowded housing",
         Id = "overcrowded-housing"
  ); join_list <- unique(c(join_list, list(HD_HQA_overcrowded_housing_tract)))


HD_HQA_rent_burden_place_0 <- pums %>%
  left_join(target_PUMAs %>% st_drop_geometry(), by = "GEOID", suffix = c("", "_target")) %>%
  filter(
    OWNERSHP == 2,
    !is.na(PUMACE20)
    ) %>%
  mutate(
    rentburden = (RENTGRS*12) / FTOTINC,
    rentburden30 = rentburden > 0.3
  ); HD_HQA_rent_burden_place <-
  bind_rows(
    pollster::moe_crosstab(
      HD_HQA_rent_burden_place_0,
      y = rentburden30,
      x = GEOID.y,
      weight = PERWT) %>%
      mutate(raceeth = "all", GEOID.y = as.factor(GEOID.y)),
    pollster::moe_crosstab_3way(
      HD_HQA_rent_burden_place_0,
      y = rentburden30,
      x = GEOID.y,
      z = RACEON,
      weight = PERWT) %>%
      rename(raceeth = RACEON),
    pollster::moe_crosstab_3way(
      HD_HQA_rent_burden_place_0,
      y = rentburden30,
      x = GEOID.y,
      z = hisp,
      weight = PERWT) %>%
      rename(raceeth = hisp),
  ) %>%
  filter(
    rentburden30 == T,
    raceeth != "!latino"
  ) %>%
  mutate(pct = pct/100) %>%
  select(
    raceeth,
    GEOID = GEOID.y,
    pct
  ) %>%
  pivot_wider(
    names_from = raceeth,
    values_from = pct
  ) %>%
  mutate(City = target_cities$cityst[match(.data$GEOID, target_cities$GEOID)],
         Geography = "Place",
         System = "Human Development",
         Area = "Housing Quality and Affordability",
         Metric = "Rent Burden ( > 30% of income)",
         Measure = "%",
         Source = "PUMS",
         Year = 2022,
         BiggerBetter = 0,
         Name = "Rent burden",
         Id = "rent-burden"
  ); join_list <- unique(c(join_list, list(HD_HQA_rent_burden_place)))


HD_HQA_rent_burden_PUMA_0 <- pums %>%
  left_join(target_PUMAs %>% st_drop_geometry(), by = "GEOID", suffix = c("", "_target")) %>%
  filter(
    OWNERSHP == 2,
    !is.na(PUMACE20)
  ) %>%
  mutate(
    rentburden = (RENTGRS*12) / FTOTINC,
    rentburden30 = rentburden > 0.3
  ); HD_HQA_rent_burden_PUMA <-
  bind_rows(
    pollster::moe_crosstab(
      HD_HQA_rent_burden_PUMA_0,
      y = rentburden30,
      x = GEOID,
      weight = PERWT) %>%
      mutate(raceeth = "all", GEOID = as.factor(GEOID)),
    pollster::moe_crosstab_3way(
      HD_HQA_rent_burden_PUMA_0,
      y = rentburden30,
      x = GEOID,
      z = RACEON,
      weight = PERWT) %>%
      rename(raceeth = RACEON),
    pollster::moe_crosstab_3way(
      HD_HQA_rent_burden_PUMA_0,
      y = rentburden30,
      x = GEOID,
      z = hisp,
      weight = PERWT) %>%
      rename(raceeth = hisp),
  ) %>%
  filter(
    rentburden30 == T,
    raceeth != "!latino"
  ) %>%
  mutate(pct = pct/100) %>%
  select(
    raceeth,
    GEOID,
    pct
  ) %>%
  pivot_wider(
    names_from = raceeth,
    values_from = pct
  ) %>%
  mutate(City = target_PUMAs$cityst[match(.data$GEOID, target_PUMAs$GEOID)],
         Geography = "PUMA",
         System = "Human Development",
         Area = "Housing Quality and Affordability",
         Metric = "Rent Burden ( > 30% of income)",
         Measure = "%",
         Source = "PUMS",
         Year = 2022,
         BiggerBetter = 0
  ); join_list <- unique(c(join_list, list(HD_HQA_rent_burden_PUMA)))


HD_HQA_homelessness_rates <- read_csv("C:/Users/taylo/CPAL Dropbox/Data Library/HUD/Continuum of Care (Point in Time)/2023_homeless_bycity.csv") %>%
  mutate(NAME_city = City) %>%
  left_join(target_metros, by = "NAME_city") %>%
  mutate(all = Total / popE,
         black = Black / pop_blackE,
         native = Native / pop_nativeE,
         aapi = (Asian+Pacific) / (pop_asianE+pop_pacificE),
         add = Multiple / pop_multiraceE,
         white = White / pop_whiteE,
         latino = NA) %>%
  select(all_of(keeps)) %>%
  mutate(City = target_metros$cityst[match(.data$GEOID, target_metros$GEOID)],
         Geography = "CSA",
         System = "Human Development",
         Area = "Housing Quality and Affordability",
         Metric = "Unhoused Rates",
         Measure = "%",
         Source = "HUD Point in Time Count",
         Year = 2023,
         BiggerBetter = 0
  ); join_list <- unique(c(join_list, list(HD_HQA_homelessness_rates)))

## BONUS: unsheltered rates (not going to use)
  # HD_HQA_unsheltered_rates <- read_csv("C:/Users/taylo/CPAL Dropbox/Data Library/HUD/Continuum of Care (Point in Time)/2023_homeless_bycity.csv") %>%
  #   mutate(NAME_city = City) %>%
  #   left_join(target_metros, by = "NAME_city") %>%
  #   mutate(black = UnshBlack / Black,
  #          native = UnshNative / Native,
  #          asian = UnshAsian / Asian,
  #          pacific = UnshPacific / Pacific,
  #          otherrace = NA,
  #          multirace = UnshMultiple / Multiple,
  #          white = UnshWhite / White,
  #          latino = NA) %>%
  #   select(all_of(keeps)) %>%
  #   mutate(System = "Human Development",
  #          Area = "Housing Quality and Affordability",
  #          Metric = "Unsheltered Rates")

HD_HFI_uninsured_population_place <- get_acs(
  "place",
  variables = c(
    "all" = "B27001_001",
    "black" = "C27001B_001",
    "native" = "C27001C_001",
    "asian" = "C27001D_001",
    "pacific" = "C27001E_001",
    "otherrace" = "C27001F_001",
    "multirace" = "C27001G_001",
    "white" = "C27001H_001",
    "latino" = "C27001I_001",
    "1all_" = "B27001_005",
    "2all_" = "B27001_008",
    "3all_" = "B27001_011",
    "4all_" = "B27001_014",
    "5all_" = "B27001_017",
    "6all_" = "B27001_020",
    "7all_" = "B27001_023",
    "8all_" = "B27001_026",
    "9all_" = "B27001_029",
    "10all_" = "B27001_033",
    "11all_" = "B27001_036",
    "12all_" = "B27001_039",
    "13all_" = "B27001_042",
    "14all_" = "B27001_045",
    "15all_" = "B27001_048",
    "16all_" = "B27001_051",
    "17all_" = "B27001_054",
    "18all_" = "B27001_057",
    "1black_" = "C27001B_004",
    "2black_" = "C27001B_007",
    "3black_" = "C27001B_010",
    "1native_" = "C27001C_004",
    "2native_" = "C27001C_007",
    "3native_" = "C27001C_010",
    "1asian_" = "C27001D_004",
    "2asian_" = "C27001D_007",
    "3asian_" = "C27001D_010",
    "1pacific_" = "C27001E_004",
    "2pacific_" = "C27001E_007",
    "3pacific_" = "C27001E_010",
    "1otherrace_" = "C27001F_004",
    "2otherrace_" = "C27001F_007",
    "3otherrace_" = "C27001F_010",
    "1multirace_" = "C27001G_004",
    "2multirace_" = "C27001G_007",
    "3multirace_" = "C27001G_010",
    "1white_" = "C27001H_004",
    "2white_" = "C27001H_007",
    "3white_" = "C27001H_010",
    "1latino_" = "C27001I_004",
    "2latino_" = "C27001I_007",
    "3latino_" = "C27001I_010"
  ),
  year = 2022,
  output = "wide"
) %>%
  filter(GEOID %in% target_cities$GEOID) %>%
  sum_columns(races = acs_races, max_cols = 3) %>%
  sum_columns(races = "all", max_cols = 18) %>%
  mutate(all = all_E / allE,
         black = case_when(blackM / blackE > threshold ~ NA,
                           black_M / black_E > threshold ~ NA,
                           TRUE ~ black_E / blackE),
         native = case_when(nativeM / nativeE > threshold ~ NA,
                            native_M / native_E > threshold ~ NA,
                            TRUE ~ native_E / nativeE),
         aapi = case_when( (asianM + pacificM) / (asianE + pacificE) > threshold ~ NA,
                           (asian_M + pacific_M) / (asian_E + pacific_E) > threshold ~ NA,
                           TRUE ~ (asian_E + pacific_E) / (asianE + pacificE)),
         add = case_when(  (otherraceM+multiraceM) / (otherraceE+multiraceE) > threshold ~ NA,
                           (otherrace_M+multirace_M) / (otherrace_E+multirace_E) > threshold ~ NA,
                           TRUE ~ (otherrace_E+multirace_E) / (otherraceE+multiraceE)),
         white = case_when(whiteM / whiteE > threshold ~ NA,
                           white_M / white_E > threshold ~ NA,
                           TRUE ~ white_E / whiteE),
         latino = case_when(latinoM / latinoE > threshold ~ NA,
                            latino_M / latino_E > threshold ~ NA,
                            TRUE ~ latino_E / latinoE)) %>%
  select(all_of(keeps)) %>%
  mutate(City = target_cities$cityst[match(.data$GEOID, target_cities$GEOID)],
         Geography = "Place",
         System = "Human Development",
         Area = "Healthy Family and Individuals",
         Metric = "Uninsured Population",
         Measure = "%",
         Source = "ACS 5-year",
         Year = 2022,
         BiggerBetter = 0
  ); join_list <- unique(c(join_list, list(HD_HFI_uninsured_population_place)))


HD_HFI_uninsured_population_tract <- get_acs(
  "tract",
  variables = c(
    "all" = "B27001_001",
    "black" = "C27001B_001",
    "native" = "C27001C_001",
    "asian" = "C27001D_001",
    "pacific" = "C27001E_001",
    "otherrace" = "C27001F_001",
    "multirace" = "C27001G_001",
    "white" = "C27001H_001",
    "latino" = "C27001I_001",
    "1all_" = "B27001_005",
    "2all_" = "B27001_008",
    "3all_" = "B27001_011",
    "4all_" = "B27001_014",
    "5all_" = "B27001_017",
    "6all_" = "B27001_020",
    "7all_" = "B27001_023",
    "8all_" = "B27001_026",
    "9all_" = "B27001_029",
    "10all_" = "B27001_033",
    "11all_" = "B27001_036",
    "12all_" = "B27001_039",
    "13all_" = "B27001_042",
    "14all_" = "B27001_045",
    "15all_" = "B27001_048",
    "16all_" = "B27001_051",
    "17all_" = "B27001_054",
    "18all_" = "B27001_057",
    "1black_" = "C27001B_004",
    "2black_" = "C27001B_007",
    "3black_" = "C27001B_010",
    "1native_" = "C27001C_004",
    "2native_" = "C27001C_007",
    "3native_" = "C27001C_010",
    "1asian_" = "C27001D_004",
    "2asian_" = "C27001D_007",
    "3asian_" = "C27001D_010",
    "1pacific_" = "C27001E_004",
    "2pacific_" = "C27001E_007",
    "3pacific_" = "C27001E_010",
    "1otherrace_" = "C27001F_004",
    "2otherrace_" = "C27001F_007",
    "3otherrace_" = "C27001F_010",
    "1multirace_" = "C27001G_004",
    "2multirace_" = "C27001G_007",
    "3multirace_" = "C27001G_010",
    "1white_" = "C27001H_004",
    "2white_" = "C27001H_007",
    "3white_" = "C27001H_010",
    "1latino_" = "C27001I_004",
    "2latino_" = "C27001I_007",
    "3latino_" = "C27001I_010"
  ),
  state = unique(target_tracts$STATE_NAME),
  year = 2022,
  output = "wide"
) %>%
  filter(GEOID %in% target_tracts$GEOID) %>%
  sum_columns(races = acs_races, max_cols = 3) %>%
  sum_columns(races = "all", max_cols = 18) %>%
  mutate(all = all_E / allE,
         black = case_when(blackM / blackE > threshold ~ NA,
                           black_M / black_E > threshold ~ NA,
                           TRUE ~ black_E / blackE),
         native = case_when(nativeM / nativeE > threshold ~ NA,
                            native_M / native_E > threshold ~ NA,
                            TRUE ~ native_E / nativeE),
         aapi = case_when( (asianM + pacificM) / (asianE + pacificE) > threshold ~ NA,
                           (asian_M + pacific_M) / (asian_E + pacific_E) > threshold ~ NA,
                           TRUE ~ (asian_E + pacific_E) / (asianE + pacificE)),
         add = case_when(  (otherraceM+multiraceM) / (otherraceE+multiraceE) > threshold ~ NA,
                           (otherrace_M+multirace_M) / (otherrace_E+multirace_E) > threshold ~ NA,
                           TRUE ~ (otherrace_E+multirace_E) / (otherraceE+multiraceE)),
         white = case_when(whiteM / whiteE > threshold ~ NA,
                           white_M / white_E > threshold ~ NA,
                           TRUE ~ white_E / whiteE),
         latino = case_when(latinoM / latinoE > threshold ~ NA,
                            latino_M / latino_E > threshold ~ NA,
                            TRUE ~ latino_E / latinoE)) %>%
  select(all_of(keeps)) %>%
  mutate(City = target_tracts$cityst[match(.data$GEOID, target_tracts$GEOID)],
         Geography = "Tract",
         System = "Human Development",
         Area = "Healthy Family and Individuals",
         Metric = "Uninsured Population",
         Measure = "%",
         Source = "ACS 5-year",
         Year = 2022,
         BiggerBetter = 0
  ); join_list <- unique(c(join_list, list(HD_HFI_uninsured_population_tract)))


HD_HFI_life_expectancy_place <- # missing data on race / not recent
  tibble(City = NA,
         Geography = "Place",
         System = "Human Development",
         Area = "Healthy Family and Individuals",
         Metric = "Life Expectancy",
         Measure = "%",
         Source = NA,
         Year = NA,
         BiggerBetter = 1
  ); join_list <- unique(c(join_list, list(HD_HFI_life_expectancy_place)))


HD_HFI_life_expectancy_tract <- # missing data on race / not recent
  tibble(City = NA,
         Geography = "Tract",
         System = "Human Development",
         Area = "Healthy Family and Individuals",
         Metric = "Life Expectancy",
         Measure = "%",
         Source = NA,
         Year = NA,
         BiggerBetter = 1
  ); join_list <- unique(c(join_list, list(HD_HFI_life_expectancy_tract)))

# need to add infant mortality later
# HD_HFI_health_index_place <- read_csv("C:/Users/taylo/CPAL Dropbox/Data Library/CDC/PLACES and 500 Cities/2023/place-level/PLACES__Local_Data_for_Better_Health__Place_Data_2023_release.csv") %>%
  # filter(LocationID %in% target_cities$GEOID) %>%
  # filter(DataValueTypeID == "AgeAdjPrv") %>%
  # filter(MeasureId %in% c("BINGE", "MHLTH", "OBESITY", "LPA", "CSMOKING", "DIABETES", "MAMMOUSE")) %>%
  # mutate(Data_Value = Data_Value / 100) %>%
  # select(Year, StateAbbr, LocationName, LocationID, Data_Value, MeasureId) %>%
  # pivot_wider(values_from = Data_Value, names_from = MeasureId) %>%
  # mutate(MAMNONUSE = 1 - MAMMOUSE) %>%
  # select(-MAMMOUSE) %>%
  # mutate()


HD_HFI_health_index_tract <- read_csv("C:/Users/taylo/CPAL Dropbox/Data Library/CDC/PLACES and 500 Cities/2023/tract-level/PLACES__Local_Data_for_Better_Health__Census_Tract_Data_2023_release.csv") %>%
  filter(MeasureId %in% c("BINGE", "MHLTH", "OBESITY", "LPA", "CSMOKING", "DIABETES", "MAMMOUSE")) %>%
  mutate(Data_Value = Data_Value / 100) %>%
  select(StateAbbr, CountyFIPS, LocationID, Data_Value, MeasureId) %>%
  pivot_wider(values_from = Data_Value, names_from = MeasureId) %>%
  mutate(MAMNOUSE = 1 - MAMMOUSE) %>%
  select(-MAMMOUSE) %>%
  mutate(across(all_of(c("BINGE", "MHLTH", "OBESITY", "LPA", "CSMOKING", "DIABETES", "MAMNOUSE")), 
                ~ as.vector(scale(.)))) %>%
  mutate(health_index = rowMeans(select(., c("BINGE", "MHLTH", "OBESITY", "LPA", "CSMOKING", "DIABETES", "MAMNOUSE")))) %>%
  select(GEOID = LocationID, health_index) %>%
  right_join(target_tracts, by = "GEOID") %>%
  select(-ends_with("city")) %>%
  mutate(totalpop = popE,
    
         pop_all = case_when(
           # popM / popE > 0 ~ NA,
                         TRUE ~ popE),
         pop_black = case_when(
           # pop_blackM / pop_blackE > 1 ~ NA,
                           TRUE ~ pop_blackE),
         pop_native = case_when(
           # pop_nativeM / pop_nativeE > 1 ~ NA,
                            TRUE ~ pop_nativeE),
         pop_aapi = case_when(
           # pop_asianM / pop_asianE > 1 ~ NA,
                           TRUE ~ pop_asianE + pop_pacificE),
         pop_add = case_when(
           # pop_otherraceM / pop_otherraceE > 1 ~ NA,
                               TRUE ~ pop_otherraceE + pop_multiraceE),
         pop_white = case_when(
           # pop_whiteM / pop_whiteE > 1 ~ NA,
                           TRUE ~ pop_whiteE),
         pop_latino = case_when(
           # pop_latinoM / pop_latinoE > 1 ~ NA,
                            TRUE ~ pop_latinoE)) %>%
#   select(-ends_with("E")) %>%
#   select(-ends_with("M")) %>%
  mutate(across(starts_with("pop"), ~ ./totalpop,
                .names = "{str_replace(.col, 'pop_', 'pct_')}")) %>%
  select(-ends_with("E", ignore.case = F)) %>%
  select(-ends_with("M", ignore.case = F)) %>%
  mutate(across(starts_with("pct"), 
                ~ .*health_index, 
                .names = "{str_remove(.col, 'pct_')}")) %>%
  select(all_of(c(city = "GEOID.y.y", keeps, "totalpop"))) %>%
  mutate(City = target_cities$cityst[match(.data$city, target_cities$GEOID)],
         Geography = "Tract",
         System = "Human Development",
         Area = "Healthy Family and Individuals",
         Metric = "Health Index",
         Measure = "Z",
         Source = "CDC PLACES",
         Year = 2023,
         BiggerBetter = 0
  ); HD_HFI_health_index_place <- HD_HFI_health_index_tract %>%
  pivot_longer(
    cols = races,
    names_to = "race",
    values_to = "weighted_health_index"
  ) %>%
  group_by(city, race, City) %>%
  summarize(health_index = weighted.mean(weighted_health_index, totalpop, na.rm = TRUE)) %>%
  pivot_wider(
    names_from = race,
    values_from = health_index
  ) %>%
  rename(GEOID = city)%>%
  mutate(Geography = "Place",
         System = "Human Development",
         Area = "Healthy Family and Individuals",
         Metric = "Health Index",
         Measure = "Z",
         Source = "CDC PLACES",
         Year = 2023,
         BiggerBetter = 0
  ); join_list <- unique(c(join_list, list(HD_HFI_health_index_place))
  ); join_list <- unique(c(join_list, list(HD_HFI_health_index_tract %>% select(-totalpop, -city))))

public_transit <- c(
  31,32,33,36,37,38,39
)

private_transit <- c(
  10,11,12,13,14,15,20,38,50,60
)

HD_Tr_public_transit_commute_0 <- pums %>% # (reported commute larger than 45 min)
  left_join(target_PUMAs %>% st_drop_geometry(), by = "GEOID", suffix = c("", "_target")) %>%
  filter(
    TRANWORK %in% public_transit,
    !is.na(PUMACE20)
  ) %>%
  mutate(
    longcommute = TRANTIME >= 45
  ); HD_Tr_public_transit_commute_place <-
  bind_rows(
    pollster::moe_crosstab(
      HD_Tr_public_transit_commute_0,
      y = longcommute,
      x = GEOID.y,
      weight = PERWT) %>%
      mutate(GEOID.y = as.factor(GEOID.y),
             raceeth = "all"),
    pollster::moe_crosstab_3way(
      HD_Tr_public_transit_commute_0,
      y = longcommute,
      x = GEOID.y,
      z = RACEON,
      weight = PERWT) %>%
      rename(
        raceeth = RACEON
      ),
    pollster::moe_crosstab_3way(
      HD_Tr_public_transit_commute_0,
      y = longcommute,
      x = GEOID.y,
      z = hisp,
      weight = PERWT) %>%
      rename(
        raceeth = hisp
      ),
  ) %>%
  filter(
    longcommute == T,
    raceeth != "!latino"
  ) %>%
  mutate(pct = pct/100) %>%
  select(
    raceeth,
    GEOID = GEOID.y,
    pct
  ) %>%
  pivot_wider(
    names_from = raceeth,
    values_from = pct
  ) %>%
  select(all_of(keeps)) %>%
  mutate(City = target_cities$cityst[match(.data$GEOID, target_cities$GEOID)],
         Geography = "Place",
         System = "Human Development",
         Area = "Transportation",
         Metric = "Percent of public transit commutes longer than 45 min",
         Measure = "%",
         Source = "PUMS",
         Year = 2022,
         BiggerBetter = 0
  ); HD_Tr_public_transit_commute_PUMA <-
  bind_rows(
    pollster::moe_crosstab(
      HD_Tr_public_transit_commute_0,
      y = longcommute,
      x = GEOID,
      weight = PERWT) %>%
      mutate(GEOID = as.factor(GEOID),
             raceeth = "all"),
    pollster::moe_crosstab_3way(
      HD_Tr_public_transit_commute_0,
      y = longcommute,
      x = GEOID,
      z = RACEON,
      weight = PERWT) %>%
      rename(
        raceeth = RACEON
      ),
    pollster::moe_crosstab_3way(
      HD_Tr_public_transit_commute_0,
      y = longcommute,
      x = GEOID,
      z = hisp,
      weight = PERWT) %>%
      rename(
        raceeth = hisp
      ),
  ) %>%
  filter(
    longcommute == T,
    raceeth != "!latino"
  ) %>%
  mutate(pct = pct/100) %>%
  select(
    raceeth,
    GEOID = GEOID,
    pct
  ) %>%
  pivot_wider(
    names_from = raceeth,
    values_from = pct
  ) %>%
  select(all_of(keeps)) %>%
  mutate(City = target_PUMAs$cityst[match(.data$GEOID, target_PUMAs$GEOID)],
         Geography = "PUMA",
         System = "Human Development",
         Area = "Transportation",
         Metric = "Percent of public transit commutes longer than 45 min",
         Measure = "%",
         Source = "PUMS",
         Year = 2022,
         BiggerBetter = 0
  ); join_list <- unique(c(join_list, list(HD_Tr_public_transit_commute_place), list(HD_Tr_public_transit_commute_PUMA)))

HD_Tr_private_transit_commute_0 <- pums %>%
  left_join(target_PUMAs %>% st_drop_geometry(), by = "GEOID", suffix = c("", "_target")) %>%
  filter(
    TRANWORK %in% private_transit,
    !is.na(PUMACE20)
  ) %>%
  mutate(
    longcommute = TRANTIME >= 45
  ); HD_Tr_private_transit_commute_place <-
  bind_rows(
    pollster::moe_crosstab(
      HD_Tr_private_transit_commute_0,
      y = longcommute,
      x = GEOID.y,
      weight = PERWT) %>%
      mutate(GEOID.y = as.factor(GEOID.y),
             raceeth = "all"),
    pollster::moe_crosstab_3way(
      HD_Tr_private_transit_commute_0,
      y = longcommute,
      x = GEOID.y,
      z = RACEON,
      weight = PERWT) %>%
      rename(
        raceeth = RACEON
      ),
    pollster::moe_crosstab_3way(
      HD_Tr_private_transit_commute_0,
      y = longcommute,
      x = GEOID.y,
      z = hisp,
      weight = PERWT) %>%
      rename(
        raceeth = hisp
      ),
  ) %>%
  filter(
    longcommute == T,
    raceeth != "!latino"
  ) %>%
  mutate(pct = pct/100) %>%
  select(
    raceeth,
    GEOID = GEOID.y,
    pct
  ) %>%
  pivot_wider(
    names_from = raceeth,
    values_from = pct
  ) %>%
  mutate(City = target_cities$cityst[match(.data$GEOID, target_cities$GEOID)],
         Geography = "Place",
         System = "Human Development",
         Area = "Transportation",
         Metric = "Percent of private transit commutes longer than 45 min",
         Measure = "%",
         Source = "PUMS",
         Year = 2022,
         BiggerBetter = 0
  ); HD_Tr_private_transit_commute_PUMA <-
  bind_rows(
    pollster::moe_crosstab(
      HD_Tr_private_transit_commute_0,
      y = longcommute,
      x = GEOID,
      weight = PERWT) %>%
      mutate(GEOID = as.factor(GEOID),
             raceeth = "all"),
    pollster::moe_crosstab_3way(
      HD_Tr_private_transit_commute_0,
      y = longcommute,
      x = GEOID,
      z = RACEON,
      weight = PERWT) %>%
      rename(
        raceeth = RACEON
      ),
    pollster::moe_crosstab_3way(
      HD_Tr_private_transit_commute_0,
      y = longcommute,
      x = GEOID,
      z = hisp,
      weight = PERWT) %>%
      rename(
        raceeth = hisp
      ),
  ) %>%
  filter(
    longcommute == T,
    raceeth != "!latino"
  ) %>%
  mutate(pct = pct/100) %>%
  select(
    raceeth,
    GEOID = GEOID,
    pct
  ) %>%
  pivot_wider(
    names_from = raceeth,
    values_from = pct
  ) %>%
  select(all_of(keeps)) %>%
  mutate(City = target_PUMAs$cityst[match(.data$GEOID, target_PUMAs$GEOID)],
         Geography = "PUMA",
         System = "Human Development",
         Area = "Transportation",
         Metric = "Percent of private transit commutes longer than 45 min",
         Measure = "%",
         Source = "PUMS",
         Year = 2022,
         BiggerBetter = 0
  ); join_list <- unique(c(join_list, list(HD_Tr_private_transit_commute_place), list(HD_Tr_private_transit_commute_PUMA)))


HD_Tr_pct_public_commute_0 <- pums %>%
  left_join(target_PUMAs %>% st_drop_geometry(), by = "GEOID", suffix = c("", "_target")) %>%
  filter(
    TRANWORK != 0,
    !is.na(PUMACE20)
  ) %>%
  mutate(
    pubtrans = TRANWORK %in% public_transit
  ); HD_Tr_pct_public_commute_place <-
  bind_rows(
    pollster::moe_crosstab(
      HD_Tr_pct_public_commute_0,
      y = pubtrans,
      x = GEOID.y,
      weight = PERWT) %>%
      mutate(GEOID.y = as.factor(GEOID.y),
             raceeth = "all"),
    pollster::moe_crosstab_3way(
      HD_Tr_pct_public_commute_0,
      y = pubtrans,
      x = GEOID.y,
      z = RACEON,
      weight = PERWT) %>%
      rename(
        raceeth = RACEON
      ),
    pollster::moe_crosstab_3way(
      HD_Tr_pct_public_commute_0,
      y = pubtrans,
      x = GEOID.y,
      z = hisp,
      weight = PERWT) %>%
      rename(
        raceeth = hisp
      ),
  ) %>%
  filter(
    pubtrans == T,
    raceeth != "!latino"
  ) %>%
  mutate(pct = pct/100) %>%
  select(
    raceeth,
    GEOID = GEOID.y,
    pct
  ) %>%
  pivot_wider(
    names_from = raceeth,
    values_from = pct
  ) %>%
  mutate(City = target_cities$cityst[match(.data$GEOID, target_cities$GEOID)],
         Geography = "Place",
         System = "Human Development",
         Area = "Transportation",
         Metric = "Percent commuting by public transit",
         Measure = "%",
         Source = "PUMS",
         Year = 2022,
         BiggerBetter = 1
  ); HD_Tr_pct_public_commute_PUMA <-
  bind_rows(
    pollster::moe_crosstab(
      HD_Tr_pct_public_commute_0,
      y = pubtrans,
      x = GEOID,
      weight = PERWT) %>%
      mutate(GEOID = as.factor(GEOID),
             raceeth = "all"),
    pollster::moe_crosstab_3way(
      HD_Tr_pct_public_commute_0,
      y = pubtrans,
      x = GEOID,
      z = RACEON,
      weight = PERWT) %>%
      rename(
        raceeth = RACEON
      ),
    pollster::moe_crosstab_3way(
      HD_Tr_pct_public_commute_0,
      y = pubtrans,
      x = GEOID,
      z = hisp,
      weight = PERWT) %>%
      rename(
        raceeth = hisp
      ),
  ) %>%
  filter(
    pubtrans == T,
    raceeth != "!latino"
  ) %>%
  mutate(pct = pct/100) %>%
  select(
    raceeth,
    GEOID = GEOID,
    pct
  ) %>%
  pivot_wider(
    names_from = raceeth,
    values_from = pct
  ) %>%
  mutate(City = target_PUMAs$cityst[match(.data$GEOID, target_PUMAs$GEOID)],
         Geography = "PUMA",
         System = "Human Development",
         Area = "Transportation",
         Metric = "Percent commuting by public transit",
         Measure = "%",
         Source = "PUMS",
         Year = 2022,
         BiggerBetter = 1
  ); join_list <- unique(c(join_list, list(HD_Tr_pct_public_commute_place), list(HD_Tr_pct_public_commute_PUMA)))

HD_Tr_pct_private_commute_0 <- pums %>%
  left_join(target_PUMAs %>% st_drop_geometry(), by = "GEOID", suffix = c("", "_target")) %>%
  filter(
    TRANWORK != 0,
    !is.na(PUMACE20)
  ) %>%
  mutate(
    privtrans = TRANWORK %in% private_transit
  ); HD_Tr_pct_private_commute_place <-
  bind_rows(
    pollster::moe_crosstab(
      HD_Tr_pct_private_commute_0,
      y = privtrans,
      x = GEOID.y,
      weight = PERWT) %>%
      mutate(GEOID.y = as.factor(GEOID.y),
             raceeth = "all"),
    pollster::moe_crosstab_3way(
      HD_Tr_pct_private_commute_0,
      y = privtrans,
      x = GEOID.y,
      z = RACEON,
      weight = PERWT) %>%
      rename(
        raceeth = RACEON
      ),
    pollster::moe_crosstab_3way(
      HD_Tr_pct_private_commute_0,
      y = privtrans,
      x = GEOID.y,
      z = hisp,
      weight = PERWT) %>%
      rename(
        raceeth = hisp
      ),
  ) %>%
  filter(
    privtrans == T,
    raceeth != "!latino"
  ) %>%
  mutate(pct = pct/100) %>%
  select(
    raceeth,
    GEOID = GEOID.y,
    pct
  ) %>%
  pivot_wider(
    names_from = raceeth,
    values_from = pct
  ) %>%
  mutate(City = target_cities$cityst[match(.data$GEOID, target_cities$GEOID)],
         Geography = "Place",
         System = "Human Development",
         Area = "Transportation",
         Metric = "Percent commuting by private transit",
         Measure = "%",
         Source = "PUMS",
         Year = 2022,
         BiggerBetter = 1
  ); HD_Tr_pct_private_commute_PUMA <-
  bind_rows(
    pollster::moe_crosstab(
      HD_Tr_pct_private_commute_0,
      y = privtrans,
      x = GEOID,
      weight = PERWT) %>%
      mutate(GEOID = as.factor(GEOID),
             raceeth = "all"),
    pollster::moe_crosstab_3way(
      HD_Tr_pct_private_commute_0,
      y = privtrans,
      x = GEOID,
      z = RACEON,
      weight = PERWT) %>%
      rename(
        raceeth = RACEON
      ),
    pollster::moe_crosstab_3way(
      HD_Tr_pct_private_commute_0,
      y = privtrans,
      x = GEOID,
      z = hisp,
      weight = PERWT) %>%
      rename(
        raceeth = hisp
      ),
  ) %>%
  filter(
    privtrans == T,
    raceeth != "!latino"
  ) %>%
  mutate(pct = pct/100) %>%
  select(
    raceeth,
    GEOID = GEOID,
    pct
  ) %>%
  pivot_wider(
    names_from = raceeth,
    values_from = pct
  ) %>%
  mutate(City = target_PUMAs$cityst[match(.data$GEOID, target_PUMAs$GEOID)],
         Geography = "PUMA",
         System = "Human Development",
         Area = "Transportation",
         Metric = "Percent commuting by private transit",
         Measure = "%",
         Source = "PUMS",
         Year = 2022,
         BiggerBetter = 1
  ); join_list <- unique(c(join_list, list(HD_Tr_pct_private_commute_place), list(HD_Tr_pct_private_commute_PUMA)))

HD_Tr_pct_work_from_home_0 <- pums %>%
  left_join(target_PUMAs %>% st_drop_geometry(), by = "GEOID", suffix = c("", "_target")) %>%
  filter(
    TRANWORK != 0,
    !is.na(PUMACE20)
  ) %>%
  mutate(
    work_from_home = TRANWORK == 80
  ); HD_Tr_pct_work_from_home <-
  bind_rows(
    pollster::moe_crosstab_3way(
      HD_Tr_pct_work_from_home_0,
      y = work_from_home,
      x = GEOID.y,
      z = RACEON,
      weight = PERWT) %>%
      rename(
        raceeth = RACEON
      ),
    pollster::moe_crosstab_3way(
      HD_Tr_pct_work_from_home_0,
      y = work_from_home,
      x = GEOID.y,
      z = hisp,
      weight = PERWT) %>%
      rename(
        raceeth = hisp
      ),
  ) %>%
  filter(
    work_from_home == T,
    raceeth != "!latino"
  ) %>%
  mutate(pct = pct/100) %>%
  select(
    raceeth,
    GEOID = GEOID.y,
    pct
  ) %>%
  pivot_wider(
    names_from = raceeth,
    values_from = pct
  ) %>%
  mutate(System = "Human Development",
         Area = "Transportation",
         Metric = "Percent working from home",
         Measure = "%",
         Source = "PUMS",
         Year = 2022)


HD_Tr_no_vehicle_0 <- pums %>%
  left_join(target_PUMAs %>% st_drop_geometry(), by = "GEOID", suffix = c("", "_target")) %>%
  filter(
    VEHICLES != 0,
    !is.na(PUMACE20)
  ) %>%
  mutate(
    no_vehicles = VEHICLES == 9
  )

HD_Tr_no_vehicle <-
  bind_rows(
    pollster::moe_crosstab_3way(
      HD_Tr_no_vehicle_0,
      y = no_vehicles,
      x = GEOID.y,
      z = RACEON,
      weight = PERWT) %>%
      rename(
        raceeth = RACEON
      ),
    pollster::moe_crosstab_3way(
      HD_Tr_no_vehicle_0,
      y = no_vehicles,
      x = GEOID.y,
      z = hisp,
      weight = PERWT) %>%
      rename(
        raceeth = hisp
      ),
  ) %>%
  filter(
    no_vehicles == T,
    raceeth != "!latino"
  ) %>%
  mutate(pct = pct/100) %>%
  select(
    raceeth,
    GEOID = GEOID.y,
    pct
  ) %>%
  pivot_wider(
    names_from = raceeth,
    values_from = pct
  ) %>%
  mutate(System = "Human Development",
         Area = "Transportation",
         Metric = "Percent with no vehicle",
         Measure = "%",
         Source = "PUMS",
         Year = "2022")

HD_FS_pop_without_supermarket_tract <- read_csv("C:/Users/taylo/CPAL Dropbox/Data Library/USDA/Food Access Research Atlas/2019_updated_2021/Food Access Research Atlas.csv") %>%
  mutate(TractFIPS = sprintf("%011.0f", CensusTract)) %>%
  left_join(target_2019tracts %>% st_drop_geometry()) %>%
  filter(!is.na(GEOID.y)) %>%
  group_by(GEOID.y) %>%
  summarize(
    lapop1 = sum(as.double(lapop1), na.rm = T),
    lawhite1 = sum(as.double(lawhite1), na.rm = T),
    lablack1 = sum(as.double(lablack1), na.rm = T),
    laasian1 = sum(as.double(laasian1), na.rm = T),
    lanhopi1 = sum(as.double(lanhopi1), na.rm = T),
    laaian1 = sum(as.double(laaian1), na.rm = T),
    laomultir1 = sum(as.double(laomultir1), na.rm = T),
    lahisp1 = sum(as.double(lahisp1), na.rm = T),
    TractPop = sum(as.double(Pop2010), na.rm = T),
    TractWhite = sum(as.double(TractWhite), na.rm = T),
    TractBlack = sum(as.double(TractBlack), na.rm = T),
    TractAsian = sum(as.double(TractAsian), na.rm = T),
    TractNHOPI = sum(as.double(TractNHOPI), na.rm = T),
    TractAIAN = sum(as.double(TractAIAN), na.rm = T),
    TractOMultir = sum(as.double(TractOMultir), na.rm = T),
    TractHispanic = sum(as.double(TractHispanic), na.rm = T)
  ) %>%
  mutate(
    all = lapop1 / TractPop,
    white = lawhite1 / TractWhite,
    black = lablack1 / TractBlack,
    aapi =  (laasian1+lanhopi1) / (TractAsian+TractNHOPI),
    native = laaian1 / TractAIAN,
    add = laomultir1 / TractOMultir,
    latino = lahisp1 / TractHispanic,
  ) %>%
  rename(GEOID = GEOID.y) %>%
  select(all_of(keeps)) %>%
  mutate(System = "Human Development",
         Area = "Food Security",
         Metric = "Population share beyond 1 mile of supermarket",
         Measure = "%",
         Source = "USDA Food Research Atlas",
         Year = "2019")


HD_FS_pct_SNAP <- get_acs(
  "place",
  variables = c(
    "all" = "B22001_001",
    "black" = "B22005B_001",
    "native" = "B22005C_001",
    "asian" = "B22005D_001",
    "pacific" = "B22005E_001",
    "otherrace" = "B22005F_001",
    "multirace" = "B22005G_001",
    "white" = "B22005H_001",
    "latino" = "B22005I_001",
    "all_" = "B22001_002",
    "black_" = "B22005B_002",
    "native_" = "B22005C_002",
    "asian_" = "B22005D_002",
    "pacific_" = "B22005E_002",
    "otherrace_" = "B22005F_002",
    "multirace_" = "B22005G_002",
    "white_" = "B22005H_002",
    "latino_" = "B22005I_002"
  ),
  year = 2022,
  output = "wide"
) %>%
  filter(GEOID %in% target_cities$GEOID) %>%
  # sum_columns(races = acs_races, max_cols = 3) %>%
  # sum_columns(races = "all", max_cols = 18) %>%
  mutate(all = all_E / allE,
         black = case_when(blackM / blackE > threshold ~ NA,
                           black_M / black_E > threshold ~ NA,
                           TRUE ~ black_E / blackE),
         native = case_when(nativeM / nativeE > threshold ~ NA,
                            native_M / native_E > threshold ~ NA,
                            TRUE ~ native_E / nativeE),
         aapi = case_when( (asianM + pacificM) / (asianE + pacificE) > threshold ~ NA,
                           (asian_M + pacific_M) / (asian_E + pacific_E) > threshold ~ NA,
                           TRUE ~ (asian_E + pacific_E) / (asianE + pacificE)),
         add = case_when(  (otherraceM+multiraceM) / (otherraceE+multiraceE) > threshold ~ NA,
                           (otherrace_M+multirace_M) / (otherrace_E+multirace_E) > threshold ~ NA,
                           TRUE ~ (otherrace_E+multirace_E) / (otherraceE+multiraceE)),
         white = case_when(whiteM / whiteE > threshold ~ NA,
                           white_M / white_E > threshold ~ NA,
                           TRUE ~ white_E / whiteE),
         latino = case_when(latinoM / latinoE > threshold ~ NA,
                            latino_M / latino_E > threshold ~ NA,
                            TRUE ~ latino_E / latinoE)) %>%
  select(all_of(keeps)) %>%
  mutate(System = "Human Development",
         Area = "Healthy Family and Individuals",
         Metric = "Percent of Households Receiving SNAP",
         Measure = "%",
         Source = "ACS 5-year",
         Year = 2022,
         Name = "",
         Id = "")

HD_AOGS <- NULL # no data

HD_ES_pre_k_enrollment <- get_acs(
  "place",
  variables = c(
    "1all" = "B01001_003",
    "2all" = "B01001_027",
    "1black" = "B01001B_003",
    "2black" = "B01001B_018",
    "1native" = "B01001C_003",
    "2native" = "B01001C_018",
    "1asian" = "B01001D_003",
    "2asian" = "B01001D_018",
    "1pacific" = "B01001E_003",
    "2pacific" = "B01001E_018",
    "1otherrace" = "B01001F_003",
    "2otherrace" = "B01001F_018",
    "1multirace" = "B01001G_003",
    "2multirace" = "B01001G_018",
    "1white" = "B01001H_003",
    "2white" = "B01001H_018",
    "1latino" = "B01001I_003",
    "2latino" = "B01001I_018",
    "1all_" = "B14007_003",
    "2all_" = "B14007_004",
    "1black_" = "B14007B_003",
    "2black_" = "B14007B_004",
    "1native_" = "B14007C_003",
    "2native_" = "B14007C_004",
    "1asian_" = "B14007D_003",
    "2asian_" = "B14007D_004",
    "1pacific_" = "B14007E_003",
    "2pacific_" = "B14007E_004",
    "1otherrace_" = "B14007F_003",
    "2otherrace_" = "B14007F_004",
    "1multirace_" = "B14007G_003",
    "2multirace_" = "B14007G_004",
    "1white_" = "B14007H_003",
    "2white_" = "B14007H_004",
    "1latino_" = "B14007I_003",
    "2latino_" = "B14007I_004"
  ),
  year = 2022,
  output = "wide"
) %>%
  filter(GEOID %in% target_cities$GEOID) %>%
  sum_columns2(races = acs_races, max_cols = 2) %>%
  sum_columns2(races = "all", max_cols = 2) %>%
  mutate(all = all_E / allE,
         black = case_when(blackM / blackE > threshold ~ NA,
                           black_M / black_E > threshold ~ NA,
                           TRUE ~ black_E / blackE),
         native = case_when(nativeM / nativeE > threshold ~ NA,
                            native_M / native_E > threshold ~ NA,
                            TRUE ~ native_E / nativeE),
         aapi = case_when( (asianM + pacificM) / (asianE + pacificE) > threshold ~ NA,
                           (asian_M + pacific_M) / (asian_E + pacific_E) > threshold ~ NA,
                           TRUE ~ (asian_E + pacific_E) / (asianE + pacificE)),
         add = case_when(  (otherraceM+multiraceM) / (otherraceE+multiraceE) > threshold ~ NA,
                           (otherrace_M+multirace_M) / (otherrace_E+multirace_E) > threshold ~ NA,
                           TRUE ~ (otherrace_E+multirace_E) / (otherraceE+multiraceE)),
         white = case_when(whiteM / whiteE > threshold ~ NA,
                           white_M / white_E > threshold ~ NA,
                           TRUE ~ white_E / whiteE),
         latino = case_when(latinoM / latinoE > threshold ~ NA,
                            latino_M / latino_E > threshold ~ NA,
                            TRUE ~ latino_E / latinoE)) %>%
  select(all_of(keeps)) %>%
  mutate(System = "Human Development",
         Area = "Healthy Family and Individuals",
         Metric = "Pre-K Enrollment Rate of Under 5 years",
         Measure = "%",
         Source = "ACS 5-year",
         Year = 2022)

HD_ES_bachelors_above <- get_acs(
  "place",
  variables = c(
    "all" = "B15002_001",
    "black" = "C15002B_001",
    "native" = "C15002C_001",
    "asian" = "C15002D_001",
    "pacific" = "C15002E_001",
    "otherrace" = "C15002F_001",
    "multirace" = "C15002G_001",
    "white" = "C15002H_001",
    "latino" = "C15002I_001",
    "1all_" = "B15002_015",
    "2all_" = "B15002_016",
    "3all_" = "B15002_017",
    "4all_" = "B15002_018",
    "5all_" = "B15002_032",
    "6all_" = "B15002_033",
    "7all_" = "B15002_034",
    "8all_" = "B15002_035",
    "1black_" = "C15002B_006",
    "2black_" = "C15002B_011",
    "1native_" = "C15002C_006",
    "2native_" = "C15002C_011",
    "1asian_" = "C15002D_006",
    "2asian_" = "C15002D_011",
    "1pacific_" = "C15002E_006",
    "2pacific_" = "C15002E_011",
    "1otherrace_" = "C15002F_006",
    "2otherrace_" = "C15002F_011",
    "1multirace_" = "C15002G_006",
    "2multirace_" = "C15002G_011",
    "1white_" = "C15002H_006",
    "2white_" = "C15002H_011",
    "1latino_" = "C15002I_006",
    "2latino_" = "C15002I_011"
  ),
  year = 2022,
  output = "wide"
) %>%
  filter(GEOID %in% target_cities$GEOID) %>%
  sum_columns(races = acs_races, max_cols = 2) %>%
  sum_columns(races = "all", max_cols = 8) %>%
  mutate(all = all_E / allE,
         black = case_when(blackM / blackE > threshold ~ NA,
                           black_M / black_E > threshold ~ NA,
                           TRUE ~ black_E / blackE),
         native = case_when(nativeM / nativeE > threshold ~ NA,
                            native_M / native_E > threshold ~ NA,
                            TRUE ~ native_E / nativeE),
         aapi = case_when( (asianM + pacificM) / (asianE + pacificE) > threshold ~ NA,
                           (asian_M + pacific_M) / (asian_E + pacific_E) > threshold ~ NA,
                           TRUE ~ (asian_E + pacific_E) / (asianE + pacificE)),
         add = case_when(  (otherraceM+multiraceM) / (otherraceE+multiraceE) > threshold ~ NA,
                           (otherrace_M+multirace_M) / (otherrace_E+multirace_E) > threshold ~ NA,
                           TRUE ~ (otherrace_E+multirace_E) / (otherraceE+multiraceE)),
         white = case_when(whiteM / whiteE > threshold ~ NA,
                           white_M / white_E > threshold ~ NA,
                           TRUE ~ white_E / whiteE),
         latino = case_when(latinoM / latinoE > threshold ~ NA,
                            latino_M / latino_E > threshold ~ NA,
                            TRUE ~ latino_E / latinoE)) %>%
  select(all_of(keeps)) %>%
  mutate(System = "Human Development",
         Area = "Healthy Family and Individuals",
         Metric = "Bachelors or Above Attained Over 25 years",
         Measure = "%",
         Source = "ACS 5-year",
         Year = 2022)


### SHELVING EDUCATION FOR NOW

# # for edu in TX
# 
# staar_4e <- haven::read_sas("C:/Users/taylo/CPAL Dropbox/Data Library/Texas Education Agency/STAAR 2022 District Aggregate Results/Grade 4/dfy22e4.sas7bdat") %>%
#   filter(DNAME %in% target_sds_TX$DNAME,
#          DISTRICT != 234903) %>%
#   mutate(Language = "English")
# 
# staar_4s <- haven::read_sas("C:/Users/taylo/CPAL Dropbox/Data Library/Texas Education Agency/STAAR 2022 District Aggregate Results/Grade 4/dfy22s4.sas7bdat") %>%
#   filter(DNAME %in% target_sds_TX$DNAME,
#          DISTRICT != 234903) %>%
#   mutate(Language = "Spanish")
# 
# staar_4 <- bind_rows(staar_4e, staar_4s) %>%
#     select(1:5,
#            Language,
#            r_ethh_d,
#            r_ethi_d,
#            r_etha_d,
#            r_ethb_d,
#            r_ethp_d,
#            r_ethw_d,
#            r_eth2_d,
#            m_ethh_d,
#            m_ethi_d,
#            m_etha_d,
#            m_ethb_d,
#            m_ethp_d,
#            m_ethw_d,
#            m_eth2_d,
#            r_ethh_unsatgl_nm,
#            r_ethi_unsatgl_nm,
#            r_etha_unsatgl_nm,
#            r_ethb_unsatgl_nm,
#            r_ethp_unsatgl_nm,
#            r_ethw_unsatgl_nm,
#            r_eth2_unsatgl_nm,
#            m_ethh_unsatgl_nm,
#            m_ethi_unsatgl_nm,
#            m_etha_unsatgl_nm,
#            m_ethb_unsatgl_nm,
#            m_ethp_unsatgl_nm,
#            m_ethw_unsatgl_nm,
#            m_eth2_unsatgl_nm
#            ) %>%
#   pivot_longer(
#     cols = -c(1:6),
#     names_to = c("Subject", "Ethnicity", "Metric"),
#     names_pattern = "(r|m)_(eth.)_(.*)",
#     values_to = "Value"
#   ) %>%
#   mutate(
#     Metric = case_when(
#       Metric == "d" ~ "total",
#       Metric == "unsatgl_nm" ~ "unsatisfactory"
#     )
#   ) %>%
#   pivot_wider(
#     names_from = Metric,
#     values_from = Value
#   ) %>%
#   mutate(
#     unsatisfactory = if_else(
#       is.na(unsatisfactory),
#       0,
#       unsatisfactory
#     ),
#     Subject = case_when(
#       Subject == "r" ~ "reading",
#       Subject == "m" ~ "math"
#     ),
#     Ethnicity = if_else(
#       Ethnicity == "ethp",
#       "etha",
#       Ethnicity
#     )
#   ) %>%
#   left_join(select(sf::st_drop_geometry(target_sds_TX), DNAME, NAME.1), by = "DNAME") %>%
#   group_by(NAME.1, Subject, Ethnicity) %>%
#   summarise(
#     total = sum(total, na.rm = TRUE),
#     unsatisfactory = sum(unsatisfactory, na.rm = TRUE)
#   ) %>%
#     mutate(
#       BelowBasic = unsatisfactory / total,
#       Ethnicity = case_when(
#         Ethnicity == "eth2" ~ "Two or more races",
#         Ethnicity == "ethi" ~ "American Indian/Alaska Native",
#         Ethnicity == "etha" ~ "Asian/Pacific Islander",
#         Ethnicity == "ethh" ~ "Hispanic",
#         Ethnicity == "ethb" ~ "Black",
#         Ethnicity == "ethw" ~ "White"
#       )
#     ) %>%
#   ungroup() %>%
#   select(-total, -unsatisfactory) %>%
#   pivot_wider(
#     names_from = Subject,
#     values_from = BelowBasic,
#     names_glue = "{Subject}{.value}",
#   ) %>%
#   mutate(Year = 2022) %>%
#   rename(Jurisdiction = NAME.1)
# 
# 
# staar_8 <- haven::read_sas("C:/Users/taylo/CPAL Dropbox/Data Library/Texas Education Agency/STAAR 2022 District Aggregate Results/Grade 8/dfy22e8.sas7bdat") %>%
#   filter(DNAME %in% target_sds_TX$DNAME,
#          DISTRICT != 234903) %>%    
#   select(1:5, 
#          r_ethh_d,
#          r_ethi_d,
#          r_etha_d,
#          r_ethb_d,
#          r_ethp_d,
#          r_ethw_d,
#          r_eth2_d,
#          m_ethh_d,
#          m_ethi_d,
#          m_etha_d,
#          m_ethb_d,
#          m_ethp_d,
#          m_ethw_d,
#          m_eth2_d,
#          r_ethh_unsatgl_nm,
#          r_ethi_unsatgl_nm,
#          r_etha_unsatgl_nm,
#          r_ethb_unsatgl_nm,
#          r_ethp_unsatgl_nm,
#          r_ethw_unsatgl_nm,
#          r_eth2_unsatgl_nm,
#          m_ethh_unsatgl_nm,
#          m_ethi_unsatgl_nm,
#          m_etha_unsatgl_nm,
#          m_ethb_unsatgl_nm,
#          m_ethp_unsatgl_nm,
#          m_ethw_unsatgl_nm,
#          m_eth2_unsatgl_nm
#          ) %>%
#   pivot_longer(
#     cols = -c(1:5), # Assuming the first 5 columns are not to be pivoted
#     names_to = c("Subject", "Ethnicity", "Metric"),
#     names_pattern = "(r|m)_(eth.)_(.*)",
#     values_to = "Value"
#   ) %>%
#   mutate(
#     Metric = case_when(
#       Metric == "d" ~ "total",
#       Metric == "unsatgl_nm" ~ "unsatisfactory"
#     )
#   ) %>%
#   pivot_wider(
#     names_from = Metric,
#     values_from = Value
#   ) %>%
#   mutate(
#     unsatisfactory = if_else(
#       is.na(unsatisfactory),
#       0,
#       unsatisfactory
#     ),
#     Subject = case_when(
#       Subject == "r" ~ "reading",
#       Subject == "m" ~ "math"
#     ),
#     Ethnicity = if_else(
#       Ethnicity == "ethp",
#       "etha",
#       Ethnicity
#     )
#   ) %>%
#   left_join(select(sf::st_drop_geometry(target_sds_TX), DNAME, NAME.1), by = "DNAME") %>%
#   group_by(NAME.1, Subject, Ethnicity) %>%
#   summarise(
#     total = sum(total, na.rm = TRUE),
#     unsatisfactory = sum(unsatisfactory, na.rm = TRUE)
#   ) %>%
#     mutate(
#       BelowBasic = unsatisfactory / total,
#       Ethnicity = case_when(
#         Ethnicity == "eth2" ~ "Two or more races",
#         Ethnicity == "ethi" ~ "American Indian/Alaska Native",
#         Ethnicity == "etha" ~ "Asian/Pacific Islander",
#         Ethnicity == "ethh" ~ "Hispanic",
#         Ethnicity == "ethb" ~ "Black",
#         Ethnicity == "ethw" ~ "White"
#       )
#     ) %>%
#   ungroup() %>%
#   select(-total, -unsatisfactory) %>%
#   pivot_wider(
#     names_from = Subject,
#     values_from = BelowBasic,
#     names_glue = "{Subject}{.value}",
#   ) %>%
#   mutate(Year = 2022) %>%
#   rename(Jurisdiction = NAME.1)
# 
# 
# NAEP_4 <- read_csv("C:/Users/taylo/CPAL Dropbox/Data Library/NAEP/2022 District Profiles/grade4/grade4_pctBelowBasic_byRace_math.csv") %>% 
#   filter(Year == 2022) %>%
#   rename(mathBelowBasic = `below Basic`) %>%
#   left_join(read_csv("C:/Users/taylo/CPAL Dropbox/Data Library/NAEP/2022 District Profiles/grade4/grade4_pctBelowBasic_byRace_reading.csv")) %>%
#   rename(
#     Ethnicity = `Race/ethnicity used to report trends, school-reported`,
#     readingBelowBasic = `below Basic`) %>%
#   mutate(
#     readingBelowBasic = readingBelowBasic/100,
#     mathBelowBasic = mathBelowBasic/100
#   ) 
# 
# # NORMALIZE
# grade_4_norm <- left_join(staar_4, NAEP_4, by = c("Year", "Jurisdiction", "Ethnicity")) %>%
#   mutate(mathDiff = mathBelowBasic.x - mathBelowBasic.y,
#          readingDiff = readingBelowBasic.x - readingBelowBasic.y) %>%
#   group_by(Ethnicity) %>%
#   summarize(
#     mathDiff = mean(mathDiff, na.rm = T),
#     readingDiff = mean(readingDiff, na.rm = T)
#   ) %>%
#   mutate_at(vars(ends_with("Diff")), ~ifelse(is.nan(.), 0, .))
# 
# 
# HD_ES_grade_4_below_basic <- NAEP_4 %>%
#   filter(! Jurisdiction %in% c("Austin", "Dallas", "Houston", "San Antonio")) %>%
#   bind_rows(staar_4 %>%
#               left_join(grade_4_norm, by = "Ethnicity") %>%
#               mutate(
#                 mathBelowBasic = mathBelowBasic - (mathDiff),
#                 readingBelowBasic = readingBelowBasic - (readingDiff)
#               )
#             ) %>%
#   select(-mathDiff, -readingDiff) %>%
#   mutate_at(vars(ends_with("BelowBasic")), round, digits=2)
# 
# 
#   
# NAEP_8 <- read_csv("C:/Users/taylo/CPAL Dropbox/Data Library/NAEP/2022 District Profiles/grade8/grade8_pctBelowBasic_byRace_math.csv") %>% 
#   rename(mathBelowBasic = `below Basic`) %>%
#   left_join(read_csv("C:/Users/taylo/CPAL Dropbox/Data Library/NAEP/2022 District Profiles/grade8/grade8_pctBelowBasic_byRace_reading.csv")) %>%
#   rename(
#     Ethnicity = `Race/ethnicity used to report trends, school-reported`,
#     readingBelowBasic = `below Basic`) %>%
#   mutate(
#     readingBelowBasic = readingBelowBasic/100,
#     mathBelowBasic = mathBelowBasic/100
#   ) 
# 
# 
# # NORMALIZE
# grade_8_norm <- left_join(staar_8, NAEP_8, by = c("Year", "Jurisdiction", "Ethnicity")) %>%
#   mutate(mathDiff = mathBelowBasic.x - mathBelowBasic.y,
#          readingDiff = readingBelowBasic.x - readingBelowBasic.y) %>%
#   group_by(Ethnicity) %>%
#   summarize(
#     mathDiff = mean(mathDiff, na.rm = T),
#     readingDiff = mean(readingDiff, na.rm = T)
#   ) %>%
#   mutate_at(vars(ends_with("Diff")), ~ifelse(is.nan(.), 0, .))
# 
# 
# HD_ES_grade_8_below_basic <- NAEP_8 %>%
#   filter(! Jurisdiction %in% c("Austin", "Dallas", "Houston", "San Antonio")) %>%
#   bind_rows(staar_8 %>%
#               left_join(grade_8_norm, by = "Ethnicity") %>%
#               mutate(
#                 mathBelowBasic = mathBelowBasic - (mathDiff),
#                 readingBelowBasic = readingBelowBasic - (readingDiff)
#               )
#             ) %>%
#   select(-mathDiff, -readingDiff) %>%
#   mutate_at(vars(ends_with("BelowBasic")), round, digits=2)


HD_TA_broadband_access <- get_acs(
  "place",
  variables = c(
    "all" = "B28008_001",
    "black" = "B28009B_001",
    "native" = "B28009C_001",
    "asian" = "B28009D_001",
    "pacific" = "B28009E_001",
    "otherrace" = "B28009F_001",
    "multirace" = "B28009G_001",
    "white" = "B28009H_001",
    "latino" = "B28009I_001",
    "all_" = "B28008_004",
    "black_" = "B28009B_004",
    "native_" = "B28009C_004",
    "asian_" = "B28009D_004",
    "pacific_" = "B28009E_004",
    "otherrace_" = "B28009F_004",
    "multirace_" = "B28009G_004",
    "white_" = "B28009H_004",
    "latino_" = "B28009I_004"
  ),
  year = 2022,
  output = "wide"
) %>%
  filter(GEOID %in% target_cities$GEOID) %>%
  # sum_columns(races = acs_races, max_cols = 2) %>%
  # sum_columns(races = "all", max_cols = 8) %>%
  mutate(all = all_E / allE,
         black = case_when(blackM / blackE > threshold ~ NA,
                           black_M / black_E > threshold ~ NA,
                           TRUE ~ black_E / blackE),
         native = case_when(nativeM / nativeE > threshold ~ NA,
                            native_M / native_E > threshold ~ NA,
                            TRUE ~ native_E / nativeE),
         aapi = case_when( (asianM + pacificM) / (asianE + pacificE) > threshold ~ NA,
                           (asian_M + pacific_M) / (asian_E + pacific_E) > threshold ~ NA,
                           TRUE ~ (asian_E + pacific_E) / (asianE + pacificE)),
         add = case_when(  (otherraceM+multiraceM) / (otherraceE+multiraceE) > threshold ~ NA,
                           (otherrace_M+multirace_M) / (otherrace_E+multirace_E) > threshold ~ NA,
                           TRUE ~ (otherrace_E+multirace_E) / (otherraceE+multiraceE)),
         white = case_when(whiteM / whiteE > threshold ~ NA,
                           white_M / white_E > threshold ~ NA,
                           TRUE ~ white_E / whiteE),
         latino = case_when(latinoM / latinoE > threshold ~ NA,
                            latino_M / latino_E > threshold ~ NA,
                            TRUE ~ latino_E / latinoE)) %>%
  select(all_of(keeps)) %>%
  mutate(System = "Human Development",
         Area = "Healthy Family and Individuals",
         Metric = "Broadband Internet Access in Household",
         Measure = "%",
         Source = "ACS 5-year",
         Year = 2022)

NG_SSC_incarceration_rates_county <- read_csv("C:/Users/taylo/CPAL Dropbox/Data Library/Vera Institute/Incarceration rates/incarceration_rates_2023.csv") %>%
  arrange(Year) %>%
  group_by(CityRep) %>%
  slice_max(Year, n = 1) %>%
  ungroup() %>%
  rename(NAME = CityRep) %>%
  left_join(target_cities) %>%
  mutate(
    all = Overall/100000,
    black = Black/100000,
    native = Native/100000,
    aapi = AAPI/100000,
    add = NA,
    white = White/100000,
    latino = Latinx/100000
  ) %>%
  select(all_of(keeps), Year = Year) %>%
  mutate(System = "Neighborhood and Governence",
         Area = "Stability and Social Cohesion",
         Metric = "Incarceration Rate",
         Measure = "%",
         Source = "Vera Institute")  
  
  

NG_SSC_pct_vacant_housing <- NULL #ACS

# NEED TO BRAINSTORM W/ MICHAEL AND ANTHONY -- DALLAS COLLEGE CAMPUSES
NG_SSC_near_higheredu <- read_csv("C:/Users/taylo/CPAL Dropbox/Data Library/IPEDS/2022 campuses/hd2022.csv") %>%
  left_join(read_csv("C:/Users/taylo/CPAL Dropbox/Data Library/IPEDS/2022 campuses/effy2022.csv") %>% filter(EFFYALEV == 21), by = "UNITID") %>%
  left_join(read_csv("C:/Users/taylo/CPAL Dropbox/Data Library/IPEDS/2022 campuses/effy2022.csv") %>% filter(EFFYALEV == 42), by = "UNITID", suffix = c("", "_pt")) %>%
  filter(C21BASIC %in% c(1:9, 14:23)) %>% # diverse degree programs
  filter(SECTOR %in% c(1,2,4,5)) %>% # not-for-profit, 2 yrs or more
  filter(CYACTIVE == 1) %>% # currently active
  sf::st_as_sf(coords = c("LONGITUD", "LATITUDE"), crs = sf::st_crs(target_cities)) %>%
  sf::st_join(sf::st_buffer(target_cities, dist = 5 * 1609.344)) %>%
  filter(!is.na(STATEFP)) %>%
  mutate(studentpop = EFYTOTLT + (EFYTOTLT_pt / 2),
         buffer_distance_mi = ((log(studentpop, base=10))^2)/5, # log^2(x)/5 of student pop
         buffer_distance_m = buffer_distance_mi * 1609.344) %>%
  filter(!is.na(studentpop)) %>%
  mutate(geometry = st_buffer(geometry, dist = buffer_distance_m)) %>%
  st_join(target_tracts)
  # filter(cityst == "Dallas, TX") %>%
  # arrange(buffer_distance_mi)

NG_SSC_traffic_stops <- NULL # missing data

# NEED TO BRAINSTORM W/ MICHAEL AND ANTHONY
NG_NEIS_near_library <- readxl::read_excel("C:/Users/taylo/CPAL Dropbox/Data Library/PROTECCT-GLAM/2023 galleries libraries archives museums/Master_Data_List_11032023.xlsx") %>%
  filter(TYPE_MAIN == "LIB") %>%
  sf::st_as_sf(coords = c("LONG", "LAT"), crs = sf::st_crs(target_cities)) %>%
  sf::st_join(st_buffer(target_cities, dist = 2 * 1609.344)) %>%
  filter(!is.na(STATEFP))

NG_NEIS_near_park <- list.files("C:/Users/taylo/CPAL Dropbox/Data Library/CDC/EPH Tracking/Park proximity/pct", 
                                pattern = ".csv$", full.names = TRUE) %>%
  lapply(function(file) {
    read_csv(file, col_types = cols(StateFIPS = col_character(),
                                    CensusTract = col_character()))
  }) %>%
  bind_rows() %>%
  rename(TractFIPS = CensusTract) %>%
  filter(TractFIPS %in% target_tracts$TractFIPS)

NG_PS_serious_crime <- list.files(path = "C:/Users/taylo/CPAL Dropbox/Data Library/FBI/Uniform Crime Reporting/2022", pattern = "*.csv", full.names = TRUE) %>%
  lapply(function(file) {
    parse_file_name <- function(filename) {
      parts <- str_match(filename, "^(\\d{4})_(Offender|Victim)(Race|Eth)_(\\w+).csv$")
      data.frame(Year = parts[1,2], Type = parts[1,3], Category = parts[1,4], City = parts[1,5])
    }
    file_info <- parse_file_name(basename(file))
    df <- read_csv(file)
    mutate(df, Year = file_info$Year, Type = file_info$Type, Category = file_info$Category, City = file_info$City)
  }) %>% 
  bind_rows()

FD_FIA_pct_bank_acct <- list.files(path = "C:/Users/taylo/CPAL Dropbox/Data Library/FDIC/Household Survey/2017-2021 Five-Year Estimate", pattern = "*.csv", full.names = TRUE) %>%
  lapply(function(file) read_csv(file, skip = 2, n_max = 27, 
                                 col_types = cols(Characteristics = col_character(), .default = col_double()))[21:27,] %>%
           mutate(msa = stringr::str_extract(basename(file), "(?<=msa-).*(?=\\.csv)"))) %>%
  bind_rows()

FD_FIA_mortgage_apprv <- read_csv("C:/Users/taylo/CPAL Dropbox/Data Library/FFIEC/2022_public_lar_csv/2022_msamd_filtered.csv")

FD_SW_median_home_value <- NULL #ACS

FD_SW_home_ownership_rate <- NULL #ACS

ED_SHPS_business_rate <- NULL # need to brainstorm
ED_SHPS_no_employees <- NULL # need to brainstorm
ED_SHPS_revenue_per_business <- NULL # need to brainstorm

ED_EJ_unemployment_rate <- NULL #ACS

ED_EJ_labor_force_rate <- NULL #ACS

ED_EJ_median_household_income <- NULL #ACS

ED_EJ_managers <- NULL #ACS

ED_EJ_service_occupation <- NULL #ACS

AO_Inc_top_5pct_income_share <- NULL #ACS

AO_Inc_child_poverty_rate <- NULL #ACS

AO_Inc_poverty_rate <- NULL #ACS

AO_Inc_bottom_quartile_outcome <- read_csv("C:/Users/taylo/CPAL Dropbox/Data Library/Opportunity Atlas/2024 pull/tract_kfr_rP_gP_p25.csv",
                                           col_types = c("c", "c", "n")) %>%
  mutate(Quartile = cut(Household_Income_at_Age_35_rP_gP_p25,
                        breaks = quantile(Household_Income_at_Age_35_rP_gP_p25, na.rm = TRUE, probs = 0:4 / 4),
                        include.lowest = TRUE,
                        labels = 1:4)) %>%
  filter(tract %in% target_tracts$TractFIPS)
  

################################################################################

################################################################################


df <- bind_rows(join_list)

dallas_only <- df %>%
  filter(Geography %in% c("Place", "CSA")) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  filter(!is.na(all)) %>%
  group_by(System, Area, Metric, Measure) %>%
  mutate(
    DallasRank = case_when(
      BiggerBetter == 1 ~ rank(-all),
      TRUE ~ rank(all)
    ),
    Ranked = max(DallasRank, na.rm = TRUE),
  ) %>%
  ungroup() %>%
  filter(City == "Dallas, TX") %>%
  select(System, Area, Metric, Measure, Dallas = all, DallasRank, Ranked)

peer_medians <- df %>%
  filter(Geography %in% c("Place", "CSA")) %>%
  group_by(System, Area, Metric, Measure) %>%
  summarize(
    PeerMedian = median(all, na.rm = TRUE),
    BestPeer = case_when(
      sum(!is.na(all)) == 0 ~ NA,
      BiggerBetter == 1 ~ max(all, na.rm = TRUE),
      BiggerBetter == 0 ~ min(all, na.rm = TRUE),
      TRUE ~ NA
    )
  ) %>%
  ungroup() %>%
  distinct() %>%
  left_join(dallas_only)
  
write_csv(peer_medians, "data/peer_medians.csv")



geographies <- c(
  "Block group" = 1,
  "Tract" = 2,
  "PUMA" = 3,
  "Place" = 4,
  "CSA" = 5,
  "County" = 6
)

smallgeo <- df %>%
  mutate(
    geolevel = geographies[Geography]
  ) %>%
  group_by(System, Area, Metric, Measure) %>%
  filter(geolevel == min(geolevel)) 

write_csv(smallgeo %>% distinct(Geography), "data/smallest_geo_available.csv")

dallas_smallgeo <- smallgeo %>%
  filter(City == "Dallas, TX")

write_csv(dallas_smallgeo, "data/dallas_smallgeo.csv")










  
city_vars <-  c("pop" = "B01001_001", 
                "pop_black" = "B01001B_001", 
                "pop_native" = "B01001C_001", 
                "pop_asian" = "B01001D_001", 
                "pop_pacific" = "B01001E_001", 
                "pop_otherrace" = "B01001F_001", 
                "pop_multirace" = "B01001G_001", 
                "pop_white" = "B01001H_001", 
                "pop_latino" = "B01001I_001",
                "pov_under" = "S1701_C03_001",
                "pov_under_black" = "S1701_C03_014",
                "pov_under_native" = "S1701_C03_015",
                "pov_under_asian" = "S1701_C03_016",
                "pov_under_pacific" = "S1701_C03_017",
                "pov_under_otherrace" = "S1701_C03_018",
                "pov_under_multirace" = "S1701_C03_019",
                "pov_under_white" = "S1701_C03_021",
                "pov_under_latino" = "S1701_C03_020",
                "no_insurance" = 
)



city_df <- tidycensus::get_acs("place",
                               variables = city_vars,
                               year = 2022) 


  
acs_tract <- tidycensus::get_acs("tract",
                                 variables = c("pop_total" = "S1701_C01_001", 
                                               "pov_total" = "S1701_C02_001", 
                                               "pop_black" = "B01001B_001", 
                                               "pop_native" = "B01001C_001", 
                                               "pop_asian" = "B01001D_001", 
                                               "pop_pacific" = "B01001E_001", 
                                               "pop_otherrace" = "B01001F_001", 
                                               "pop_multirace" = "B01001G_001", 
                                               "pop_white" = "B01001H_001", 
                                               "pop_latino" = "B01001I_001"),
                                 year = 2022,
                                 state = "TX",
                                 county = c("Dallas", "Denton", "Collin County", "Kaufman", "Rockwall"),
                                 geometry = TRUE,
                                 output = "wide") %>%
  select(-ends_with("M")) %>%
  rename_with(~gsub("E$", "", .), ends_with("E")) %>%
  mutate(
    poppct_black = pop_black / pop_total,
    poppct_asian = pop_asian / pop_total,
    poppct_white = pop_white / pop_total,
    poppct_latino = pop_latino / pop_total,
    tract_sqmi = as.numeric(units::set_units(sf::st_area(geometry), mi^2)),
    pop_density = pop_total/tract_sqmi,
    pop_density_cut = ifelse(pop_density > 25000, 25000, pop_density)
    # poppct_other = (pop_native + pop_pacific + pop_otherrace + pop_multirace) / pop_total,
  ) %>%
  mutate(
    max_race_info = pmap(
      list(poppct_black, poppct_asian, poppct_white, poppct_latino), 
      ~ {
        percentages <- c(...); 
        names(percentages) <- c("Black", "Asian", "White", "Hispanic/Latino");
        sorted <- sort(percentages, decreasing = TRUE); 
        max_race <- names(sorted)[1]; 
        second_max_race <- names(sorted)[2]; 
        max_diff <- sorted[1] - sorted[2]; 
        list(max_race, max_diff)
      }
    )
  ) %>%
  unnest_wider(max_race_info, names_sep = "_") %>%
  rename(
    pop_largestrace = max_race_info_1,
    pop_racediff = max_race_info_2
  ) %>%
  sf::st_as_sf() %>%
  sf::st_set_crs(4269)

acs_place <- tidycensus::get_acs("place",
                                 variables = c("pop_total" = "S1701_C01_001", 
                                               "pov_total" = "S1701_C02_001", 
                                               "pop_black" = "B01001B_001", 
                                               "pop_native" = "B01001C_001", 
                                               "pop_asian" = "B01001D_001", 
                                               "pop_pacific" = "B01001E_001", 
                                               "pop_otherrace" = "B01001F_001", 
                                               "pop_multirace" = "B01001G_001", 
                                               "pop_white" = "B01001H_001", 
                                               "pop_latino" = "B01001I_001"),
                                 year = 2022,
                                 state = "TX",
                                 geometry = TRUE,
                                 output = "wide") %>%
  filter(NAME == "Dallas city, Texas") %>%
  select(-ends_with("M")) %>%
  rename_with(~gsub("E$", "", .), ends_with("E")) %>%
  mutate(
    poppct_black = pop_black / pop_total,
    poppct_asian = pop_asian / pop_total,
    poppct_white = pop_white / pop_total,
    poppct_latino = pop_latino / pop_total,
    # poppct_other = (pop_native + pop_pacific + pop_otherrace + pop_multirace) / pop_total
  ) %>%
  mutate(geometry = sf::st_difference(geometry, dallas_remove))

dallas_tracts <- sf::st_join(acs_tract, acs_place, sf::st_intersects, suffix = c("", "_DAL")) %>%
  filter(!is.na(NAM_DAL)) %>%
  select(geometry, everything()) %>%
  sf::st_simplify(dTolerance = 1) %>%
  mutate(
    tract_area = sf::st_area(geometry)
  ) %>%
  rowwise() %>%
  mutate(
    intersection = tryCatch(
      {
        inter = sf::st_intersection(geometry, acs_place$geometry[1])
        if (length(inter) == 0) { NA } else { inter }
      },
      error = function(e) NA
    ),
    intersection_area = ifelse(is.na(intersection), 0, sf::st_area(intersection)),
    area_proportion = as.numeric(intersection_area / tract_area),
    tract_name = str_remove(NAM, ";.*")
  ) %>%
  ungroup() %>%
  filter(area_proportion >= 0.20)
# %>%  pivot_longer(4:(which(names(.) == "GEOID_DAL") - 1), names_to = "variable")

bounds <- as.vector(sf::st_bbox(dallas_tracts))