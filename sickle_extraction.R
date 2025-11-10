library(tidyverse)
library(sf)
library(geodata) 
library(exactextractr)
library(terra)
library(raster)


DATA_DIR <- "Data"    
OUTPUT_DIR <- "Output"    
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, showWarnings = FALSE)
data_path <- function(...) file.path(DATA_DIR, ...)

##Set working directory to the project root
# setwd("path/to/project/root")

country_codes_file_path <- data_path("countries_with_regional_codes.csv")
country_codes <- read.csv(country_codes_file_path, header = TRUE)

all_country_code <- country_codes[ , 3]

# save copies of gadm SpatVector as .rds files in data folder
import_gadm <- function(ISO, level){
  geodata::gadm(country = ISO, level = level, path = "./Data/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed/gadm", version = "4.0")
}

# loop over function to import admin0 (country) shapefiles
map2(all_country_code, 0, import_gadm)

# create a list of all imported countries
countries_all <- lapply(c(all_country_code), 
                           function(x){
                             list.files(path = "./Data/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed/gadm/gadm",
                                        pattern = paste0("*", x , "_0_pk.rds"), 
                                        full.names = TRUE)
                           }) |> 
  unlist()

# import MAP raster file
# downloaded from: https://data.malariaatlas.org/maps?layers=Blood_Disorders:201201_Global_Sickle_Haemoglobin_HbS_Allele_Frequency,Malaria:202508_Global_Pf_Parasite_Rate
sickle_file_path <- data_path("2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed","2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.geotiff")
sickle <- raster(sickle_file_path)
crs(sickle)

par(mar = c(1,1,1,1)) # set up plot margins
plot(sickle) # quick view
#The legend shows a prevalence of allele frequencies


# import worldpop raster
# downloaded from: https://hub.worldpop.org/geodata/summary?id=24777
# import raster file
world_pop_path <- data_path("2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed","ppp_2010_1km_Aggregated.tif")
pop <- raster(world_pop_path)
crs(pop) # matching crs to sickle
# this is showing population counts

results_df <- data.frame(country = character(),
                         results_mean = numeric(),
                         results_median = numeric(),
                         results_weighted = numeric())

for (i in countries_all) {
  object <- readRDS(i) # read in object
  object <- terra::vect(sf::st_as_sf(object)) # unpack SpatVector
  country_sf <- st_as_sf(object)  # transform to sf object
  sub_sickle <- crop(sickle, extent(country_sf))
  sub_sickle <- mask(sub_sickle, country_sf)
  # mean of allele frequency
  results_mean <- exactextractr::exact_extract(sub_sickle, country_sf, fun = "mean", progress = FALSE)
  # median of allele frequency
  results_median <- exactextractr::exact_extract(sub_sickle, country_sf, fun = "median", progress = FALSE)
  results_df <- rbind(results_df, data.frame(Country.Code = country_sf$ID_0,
                                             results_mean = results_mean,
                                             results_median = results_median))
}


results_weighted_df <- data.frame(country = character(),
                         results_weighted = numeric())

for (i in countries_all) {
  object <- readRDS(i) # read in object
  object <- terra::vect(sf::st_as_sf(object)) # unpack SpatVector
  country_sf <- st_as_sf(object)  # transform to sf object
  # crop to shapefile
  sub_sickle <- crop(sickle, extent(country_sf))
  sub_sickle <- mask(sub_sickle, country_sf)
  sub_pop <- crop(pop, extent(country_sf))
  sub_pop <- mask(sub_pop, country_sf)
  # resample raster so that it matches the sickle one
  pop_resample <- raster::resample(sub_pop, sub_sickle, method = "bilinear")
  # population weighted mean of allele frequency
  results_weighted <- exactextractr::exact_extract(sub_sickle, country_sf, fun = "weighted_mean", weights = raster::area(pop_resample), progress = TRUE)
  results_weighted_df <- rbind(results_weighted_df, data.frame(Country.Code = country_sf$ID_0,
                                             results_weighted = results_weighted))
}

merged_results <- results_df %>%
  full_join(results_weighted_df, by = "Country.Code")
merged_results
summary(merged_results)

#write.csv(merged_results, file = "Data/Sickle_cell/sample_list.csv", row.names = FALSE)
cat("\nScripts completed successfully.\n")
