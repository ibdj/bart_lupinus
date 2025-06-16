# R scripts for short course "Species distribution models with Bayesian statistics in R"
# A. Marcia Barbosa (https://modtools.wordpress.com)


# SETUP ####

# load required packages:
library(terra)  # mapping
library(geodata)  # downloading data and maps
library(fuzzySim)  # cleaning and gridding data
library(sf) # to read geopackages
library(spocc) # to get occurrence data from many sources

# # the following command sets the working directory to the folder that contains this script:
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # if you get an error, try running this again, or do it via the RStudio menu "Session - Set Working Directory - To Source File Location"
# getwd()

# # the next command creates a folder for the output files of this course (if it doesn't already exist):
# if (!file.exists("../outputs")) dir.create("../outputs")  # '../' goes up one level from the current working directory, so this creates the 'outputs' folder just outside the 'scripts' folder

# download/import a world countries map for geographical context:
countries <- geodata::world(path = "../outputs")


# IMPORT SPECIES OCCURRENCE DATA ####

# the code below imports a sample dataset; you can place your own data in the course materials 'practicals/inputs' folder and change the path below accordingly
my_presences <- st_read("/Users/ibdj/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/gis/mapping_plants.gpkg", layer = "lupines_manual_observations") |> filter(!st_is_empty(geom)) |> select(geom)
#my_presences <- read.csv("../inputs/elephant_tracking/African_elephant_HwangeNP.csv")
nrow(my_presences)

# convert the table to a spatial object (like when you import a delimited text file into a GIS, you need to specify the names of the columns that contain the spatial coordinates (point geometry), and preferably also the geographic projection / coordinate reference system):
#head(my_presences) # see which columns contain the spatial coordinates to provide as 'geom' below
my_presences <- terra::vect(my_presences)
#my_presences <- terra::vect(my_presences, geom = c("location.long", "location.lat"), keepgeom = TRUE, crs = "EPSG:4326")
plot(my_presences, cex = 0.2)

# plot a part of the countries map and these presences on top:
plot(countries, ext = terra::ext(my_presences) + 2, col = "tan", background = "lightblue")
plot(my_presences, cex = 0.2, col = "blue", add = TRUE)  # if your points don't overlap where they should, you may need to use terra::project() - ask me for help if needed!


# DOWNLOAD ADDITIONAL OCCURRENCE DATA FROM GBIF ####

# example terrestrial species:
my_species <- "Lupinus nootkatensis"  # if you choose another species, change it here, not all over the script!


# mind that data download takes a long time when there are many occurrences! first, check how many occurrences are available for this species, without downloading them:
geodata::sp_occurrence(genus = "", species = my_species, download = FALSE)

# if it's not (much) more than 10,000 records (just to keep computing times reasonable DURING THE COURSE), you can download them all:
gbif_raw <- geodata::sp_occurrence(genus = "", species = my_species, fixnames = FALSE)

# if your species spreads over a very large area or has many occurrence records, you can download points only within a specified window of longitude and latitude coordinates (otherwise, scripts may take too long to run DURING THE COURSE!):
plot(countries, col = "tan", background = "lightblue")  # see coordinates in the plot axes

# set the coordinates of your desired extent, e.g.:
download_window <- terra::ext(-25, -10, 63, 67)  # xmin, xmax, ymin, ymax


plot(download_window, border = "green", lwd = 2, add = TRUE)  # confirm it's where you want it on the map
# plot country borders within 'download_window' only:
plot(countries, col = "tan", background = "lightblue", ext = download_window)
#plot(my_presences, cex = 0.2, col = "blue", add = TRUE)

# if global data are too much (too many occurrences or too large area for timely computation DURING THE COURSE), download GBIF data from the chosen window only:
# for the TERRESTRIAL EXAMPLE SPECIES, add args=c("year=2018,2025") to the below command, otherwise the download will take much longer!
gbif_raw <- sp_occurrence(genus = "", species = my_species, ext = download_window, fixnames = FALSE)

# check how the data are organized:
head(gbif_raw)
names(gbif_raw)

# NOTE: as per the ?sp_occurrence help file, "Before using this function, please first check the GBIF data use agreement and see the note below about how to cite these data". If you plan to use GBIF data in any report or publication, you can download the data directly from www.gbif.org (then import the .csv to R) and note down the DOI and citation for the entire download. See also the new "derived datasets" tool (https://www.gbif.org/citation-guidelines#derivedDatasets) for obtaining a DOI based on the 'datasetKey' column above. It is very important to properly cite the data sources! GBIF is not a source, just a repository for many people who put in very hard work to collect these data and make them available for you.

# convert to spatial object:
names(gbif_raw)  # see which columns contain the spatial coordinates, provided as 'geom' below:
gbif_raw <- terra::vect(gbif_raw, geom = c("decimalLongitude", "decimalLatitude"), keepgeom = TRUE, crs = "EPSG:4326")

plot(countries, ext = terra::ext(gbif_raw) + 1, col = "tan", background = "lightblue", main = paste(my_species, "in Iceland"))
plot(gbif_raw, col = "darkgreen", add = TRUE)  # compare e.g. with the range map of this species at https://www.iucnredlist.org to assess if the species' distribution is well represented in this region
## There seem to be missing some of the northern most area of its extent compaired to IUCN: https://www.iucnredlist.org/species/181008073/223031019

# create a folder and save the data there:
output_data_folder <- paste0("../outputs/", my_species, "/occ_data")
dir.create(output_data_folder, recursive = TRUE, showWarnings = FALSE)

write.csv(gbif_raw, paste0(output_data_folder, "/gbif_raw.csv"), row.names = FALSE)
write.csv(as.vector(download_window), paste0(output_data_folder, "/download_window.csv"), row.names = FALSE)


# from now on, you don't need to download these data again - you can just import them from the .csv:
gbif_raw <- read.csv(paste0(output_data_folder, "/gbif_raw.csv"))
head(gbif_raw)
nrow(gbif_raw)


# CLEAN SPECIES OCCURRENCE DATA ################################################################################################

# biodiversity data typically contain many errors; careful mapping, inspection and cleaning are necessary!

# automatically remove records of absence, or with duplicate, impossible, unlikely (equal or zero), imprecise, or overly uncertain spatial coordinates (potentially over 50 km away from the actual observation site):
gbif_clean <- fuzzySim::cleanCoords(gbif_raw, coord.cols = c("decimalLongitude", "decimalLatitude"), abs.col = "occurrenceStatus", uncert.col = "coordinateUncertaintyInMeters", uncert.limit = 50000)

plot(countries, col = "tan", alpha = 0.3, add = TRUE)

# but note that this will only discard records where coordinate uncertainty is adequately reported in the dataset, which may not always be the case! Careful mapping and visual inspection are necessary

# save the cleaned data to disk as a .csv file:
write.csv(data.frame(gbif_clean), paste0(output_data_folder, "/gbif_clean.csv"), row.names = FALSE)

# see the data you have on disk so far:
list.files(output_data_folder)

# from now on, you don't need to download and clean the data again - you can just import them from the .csv:
gbif_clean <- read.csv(paste0(output_data_folder, "/gbif_clean.csv"))
head(gbif_clean)
nrow(gbif_clean)


# COMBINE THE OCCURRENCE DATA FROM DIFFERENT SOURCES ####

#names(my_presences)
names(gbif_clean)
unique(gbif_clean$species)  # before combining (below), make sure it's the same species from your 'my_presences' data! NOT THE CASE FOR THE MARINE EXAMPLE SPECIES

# make relevant column names match (for 'rbind' within 'appendData' below):
#names(my_presences)[grep("lon", names(my_presences))]
#names(my_presences)[grep("lat", names(my_presences))]
#names(my_presences) <- gsub("location.long", "decimalLongitude", names(my_presences))
#names(my_presences) <- gsub("location.lat", "decimalLatitude", names(my_presences))

presences <- gbif_clean
# if 'my_presences' are from the SAME SPECIES (not for the marine example):
presences <- fuzzySim::appendData(gbif_clean, data.frame(my_presences), fill = FALSE)

head(presences)
tail(presences)

# convert to spatial object:
presences <- terra::vect(presences, geom = c("decimalLongitude", "decimalLatitude"), keepgeom = TRUE, crs = "EPSG:4326")

# map coloured according to the last column:
plot(presences, ncol(presences), col = c("darkgreen"), background = "lightblue", main = paste(my_species, "in Iceland"))
plot(countries, col = "sienna4", alpha = 0.3, add = TRUE)


#### ASSIGNMENT 1: if you've tried different data than the example, post the map you obtained with the above commands


# save the combined occurrence data:
write.csv(presences, paste0(output_data_folder, "/presences_all.csv"), row.names = FALSE)


# (DOWN)LOAD ENVIRONMENTAL VARIABLES ####

# NOTE: if you plan to project your model(s) to other regions or time periods, you'll need to use only variables that are available there too!

# get some TERRESTRIAL layers:
layers <- geodata::worldclim_global(var = "bio", res = 5, path = "../outputs")
layers
# for the meanings of these variables, see https://www.worldclim.org/data/bioclim.html
# see also geodata::worldclim_country() if you want to download higher-resolution variables for smaller parts of the world

# or instead get some MARINE layers:
# marine_layer_names <- c("Chlorophyll", "Current.Velocity", "Ice.thickness", "Salinity", "Temperature")  # using only those with future projections available at https://bio-oracle.org/data/2.0, because we intend to project models to the future later on
# layers <- terra::rast(lapply(marine_layer_names, function(v) geodata::bio_oracle(path = "../outputs", var = v, stat = "Mean", benthic = FALSE)))  # check help file for other values of "stat" and "benthic"

names(layers)
# simplify variable names by replacing the prefix with nothing:
names(layers) <- gsub("wc2.1_5m_", "", names(layers))  # for the terrestrial example data
#names(layers) <- gsub("Present.", "", names(layers))  # for the marine example data
names(layers)

# plot a couple of layers to see how they look:
plot(layers)
plot(layers[[1]], main = names(layers)[1])
plot(layers[[3]], main = names(layers)[3])



# BIO1 = Annual Mean Temperature
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO3 = Isothermality (BIO2/B IO7) (×100)
# BIO4 = Temperature Seasonality (standard deviation ×100)
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO8 = Mean Temperature of Wettest Quarter
# BIO9 = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation
# BIO13 = Precipitation of Wettest Month
# BIO14 = Precipitation of Driest Month
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter

