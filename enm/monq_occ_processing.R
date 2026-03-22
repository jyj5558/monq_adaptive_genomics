#gbif_download
# setup -------------------------------------------------------------------

rm(list = ls())

library(rgbif)
library(tidyverse)

setwd("C:/Users/jyj55/OneDrive - purdue.edu/DeWoody_Lab/Dissertation/Montezuma quail/enm")

# data --------------------------------------------------------------------

my_species <- 'Cyrtonyx_montezumae'

key <- 
  name_backbone(my_species) |> 
  pull(usageKey)

key

gbif_download <- 
  occ_download(
    pred('taxonKey', key), 
    format = 'SIMPLE_CSV',
    user = 'jyj5558',
    pwd = '',
    email = 'jyj5558@nate.com')

gbif_download

# save citation -----------------------------------------------------------

gbif_download |> 
  write_rds(
    paste0(
      'occurrence/raw/',
      my_species,
      '_key.rds'))

read_rds(
  paste0(
  'occurrence/raw/',
  my_species,
  '_key.rds'))

# check download processing -----------------------------------------------

occ_download_wait(gbif_download)

data <- 
  occ_download_get(
    gbif_download, 
    path = 'occurrence/raw', 
    overwrite = TRUE) |> 
  occ_download_import()

# data <- 
#   occ_download_get(
#     '0037423-240906103802322',
#     path = 'data/raw',
#     overwrite = TRUE)

# save data ---------------------------------------------------------------

data |>
  write_csv(
    paste0(
      'occurrence/raw/',
      my_species,
      '_gbif_raw.csv'))



#gbif_clean
# setup -------------------------------------------------------------------

rm(list = ls())

library(CoordinateCleaner)
library(scrubr)
library(rgbif)
library(tidyverse)

# data --------------------------------------------------------------------

my_species <- 'Cyrtonyx_montezumae'

gbif_raw <- 
  read_csv(
    paste0(
      'occurrence/raw/',
      my_species,
      '_gbif_raw.csv'))

gbif_pre_clean <- 
  gbif_raw |> 
  select(
    id = gbifID,
    datasetKey,
    species, 
    subspecies = infraspecificEpithet,
    locality = locality,
    status = occurrenceStatus,
    x = decimalLongitude,
    y = decimalLatitude,
    accuracy = coordinateUncertaintyInMeters,
    record = basisOfRecord,
    year = year, 
    country = countryCode,
    state = stateProvince,
    institution = institutionCode) # consider using "locality" and not using "stateProvince"
   
gbif_pre_clean2 <- 
  gbif_pre_clean |>
  filter(
    !if_any(
      c(
        x, 
        y, 
        #accuracy, 
        year,
        state), 
      ~ is.na(.x))) |>  
  filter(x < 0, y > 0) |> 
  filter(
    !(subspecies %in% c('montezumae', 'merriami', 'rowleyi', 'sallei'))) |>
  filter(
    !(record %in% c('MATERIAL_CITATION', 'MATERIAL_SAMPLE', 'PRESERVED_SPECIMEN'))) |> 
  filter(  
    status == 'PRESENT',
    institution != 'iNaturalist',
    year >= 1970,
    accuracy <= 1700 | is.na(accuracy)) |> # based on root(2.9km2) (maximum home range 2.9km2 from Chavarria et al., 2017; Greene et al., 2020)
  select(!status) |> 
  mutate(source = 'gbif') |> 
  distinct(x, y, year, .keep_all = TRUE)

captive <- gbif_pre_clean2 %>% filter(grepl("captive", 'locality', ignore.case = TRUE))
captive

gbif_clean <- 
  gbif_pre_clean2 |>
  clean_coordinates(
    lon = 'x',
    lat = 'y',
    tests = c(
      'capitals', 
      'centroids',
      'equal', 
      'gbif', 
      'institutions', 
      'outliers', 
      'seas', 
      'zeros'),
    value = 'clean') |> 
  coord_incomplete() |> 
  coord_imprecise() |> 
  coord_impossible() |> 
  coord_unlikely()

# save data ---------------------------------------------------------------

gbif_clean |> 
  write_csv(
    paste0(
      'occurrence/processed/',
      my_species,
      '_gbif_clean.csv'))

# create derived dataset --------------------------------------------------

# https://www.gbif.org/derived-dataset/about)

derived_data <-
  gbif_clean |>
  summarize(
    n = n(),
    .by = datasetKey)

# test derived dataset

derived_dataset_prep(
  citation_data = derived_data,
  title = 'Derived Dataset Cyrtonyx montezumae',
  description = 
    'This data was filtered using CoordinateCleaner and scrubr',
  #source_url = 
    'https://github.com/jyj5558/monq/data/processed/gbif_clean.csv',
  gbif_download_doi = '10.15468/dl.cvcrnz')

# If output looks ok, run derived_dataset 

derived_dataset(
  citation_data = derived_data,
  title = 'Derived Dataset Cyrtonyx montezumae',
  description = 
    'This data was filtered using CoordinateCleaner and scrubr',
  #source_url = 
    'https://github.com/jyj5558/monq/data/processed/gbif_clean.csv',
  gbif_download_doi = '10.15468/dl.cvcrnz',
  user = 'jyj5558',
  pwd = '')



#inat
# setup -------------------------------------------------------------------

rm(list = ls())

library(rinat)
library(tidyverse)

# data --------------------------------------------------------------------

#total of observations

inat_metadata <- 
  get_inat_obs(
    query = 'Cyrtonyx montezumae', 
    meta = TRUE) |> 
  pluck('meta')

# download data

inat_data <-
  get_inat_obs(
    query = 'Cyrtonyx montezumae',
    quality = 'research',
    geo = TRUE,
    maxresults = 10000,
    meta = FALSE) |> 
  as_tibble()

# save data ---------------------------------------------------------------

my_species <- 'Cyrtonyx_montezumae'

inat_data |>
  write_csv(
    paste0(
      'occurrence/raw/',
      my_species,
      '_inat_raw.csv'))



#inat_clean
# setup -------------------------------------------------------------------

rm(list = ls())

library(CoordinateCleaner)
library(scrubr)
library(tidyverse)

# data --------------------------------------------------------------------

my_species <- 'Cyrtonyx_montezumae'

inat_pre_clean <- 
  read_csv(
    paste0(
      'occurrence/raw/',
      my_species,
      '_inat_raw.csv')) |>
  select(
    id,
    species = scientific_name,
    date = datetime,
    locality = place_guess,
    x = longitude,
    y = latitude,
    accuracy = positional_accuracy, 
    coordinates_obscured) |> 
  mutate(
    date = as_date(date),
    year = year(date), 
    .before = date) |> 
  select(!date)

inat_pre_clean2 <- 
  inat_pre_clean |>
  filter(
    !if_any(
      c(
        x, 
        y, 
        #accuracy, 
        year), 
      is.na)) |>  
  filter(x < 0, y > 0) |> 
  filter(  
    year >= 1970,
    accuracy <= 1700 | is.na(accuracy),
    coordinates_obscured == 'FALSE') |>
  mutate(source = 'INaturalist') |> 
  distinct(x, y, year, .keep_all = TRUE)

inat_clean <- 
  inat_pre_clean2 |>
  clean_coordinates(
    lon = 'x',
    lat = 'y',
    tests = c(
      'capitals', 
      'centroids',
      'equal', 
      'gbif', 
      'institutions', 
      'outliers', 
      'seas', 
      'zeros'),
    value = 'clean') |> 
  coord_incomplete() |> 
  coord_imprecise() |> 
  coord_impossible() |> 
  coord_unlikely()

# save data ---------------------------------------------------------------

inat_clean |>
  write_csv(
    paste0(
      'occurrence/processed/',
      my_species,
      '_inat_clean.csv'))



#full_dataset
# setup -------------------------------------------------------------------

rm(list = ls())

library(sf)
library(tidyverse)

# data --------------------------------------------------------------------

list.files(
  'occurrence/processed',
  pattern = 'clean.csv$',
  full.names = TRUE) |> 
  map(~ .x |> 
        read_csv()) |> 
  set_names('gbif','inat') |> 
  list2env(.GlobalEnv)

# gbif
# inat

gbif_metadata <- 
  gbif |> 
  select(id, datasetKey, institution)

full_dataset <- 
  gbif |> 
  select(!c(datasetKey, institution)) |> 
  bind_rows(
    inat |> 
      select(!coordinates_obscured)) |> 
  mutate(
    across(
      c(id, species, year, source),
      ~ as_factor(.x)))

dir.create('occurrence/processed/final')

my_species <- 'Cyrtonyx_montezumae'

full_dataset |> 
  write_csv(
    paste0(
      'occurrence/processed/final/',
      my_species,
      '_occs_clean.csv'))
