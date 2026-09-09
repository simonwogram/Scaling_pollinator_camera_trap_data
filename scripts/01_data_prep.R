library(tidyverse)
library(lubridate)
library(hms)
library(geodata)
library(terra)
library(data.table)

path_cams <- "P:/R_P1/raw_data/camera/DataFigure8/Figure8CSV"
path_temp <- "P:/R_P1/raw_data/temperature"
path_ana  <- "P:/R_P1/r_data"
path_clim <- "P:/R_P1/raw_data/climate"
path_gbif <- "P:/R_P1/raw_data/gbif"

recording_file <- "P:/R_P1/raw_data/camera/recording_dates.csv"
gbif_file <- "P:/R_P1/raw_data/gbif/0001878-260126135527185/occurrence.txt"

# Camera detections
cam_orig <- list.files(path_cams, full.names = TRUE) |>
  map_dfr(~ read.csv(.x, header = FALSE)) |>
  select(V1:V6) |>
  setNames(c("SystemID", "camera", "date_orig", "time_orig", "conf", "class")) |>
  mutate(
    SystemID = recode(
      as.character(SystemID),
      "1" = "S1_146",
      "2" = "S2_123",
      "3" = "S3_194",
      "4" = "S4_199",
      "5" = "S5_187"
    ),
    date_time = as.POSIXct(
      paste(date_orig, str_pad(time_orig, 6, pad = "0")),
      format = "%Y%m%d %H%M%OS"
    ),
    date = ymd(as.character(date_orig)),
    time = as_hms(as.POSIXct(str_pad(time_orig, 6, pad = "0"), format = "%H%M%OS")),
    time_hr = as_hms(as.POSIXct(str_pad(time_orig, 6, pad = "0"), format = "%H"))
  ) |>
  filter(conf >= 70, class != 1)

cam_hr <- cam_orig |>
  count(SystemID, date_time, date, time, time_hr, class, name = "abun") |>
  group_by(SystemID, date, time_hr, class) |>
  summarise(abun = sum(abun), .groups = "drop")

recording_dates <- read.csv(recording_file) |>
  rename(SystemID = Level1Folder) |>
  separate_rows(RecordingDates, sep = ",\\s*") |>
  transmute(SystemID, date = ymd(RecordingDates))

dates_full <- recording_dates |>
  crossing(
    time_hr = hms(hours = 5:22),
    class = sort(unique(cam_hr$class))
  ) |>
  mutate(
    dates_hr = as.POSIXct(
      paste(date, time_hr),
      format = "%Y-%m-%d %H:%M:%S"
    )
  )

cam_hr <- dates_full |>
  left_join(cam_hr, by = c("SystemID", "date", "time_hr", "class")) |>
  filter(dates_hr >= as.POSIXct("2019-06-24 14:00:00")) |>
  mutate(abun = replace_na(abun, 0))

# Temperature data
file_pattern <- list.files(path_temp, pattern = "temp", full.names = TRUE)

temp_orig <- map_dfr(file_pattern, read.csv2)
names(temp_orig)[1:4] <- c("DateTime", "temp", "temp_max", "temp_min")

temp_orig <- temp_orig |>
  mutate(date_time = as.POSIXct(DateTime, tz = "Europe/Copenhagen")) |>
  select(-DateTime) |>
  distinct()

cam_hr <- cam_hr |>
  left_join(select(temp_orig, date_time, temp), by = c("dates_hr" = "date_time"))

cam_day <- cam_hr |>
  group_by(SystemID, date, class) |>
  summarise(
    temp = max(temp, na.rm = TRUE),
    abun = sum(abun, na.rm = TRUE),
    .groups = "drop"
  )

save(cam_hr, file = file.path(path_ana, "cam_hr.RData"))
save(temp_orig, file = file.path(path_ana, "temp_orig.RData"))
save(cam_day, file = file.path(path_ana, "cam_day_tmax.RData"))

# WorldClim BIO5
bio5_crop <- worldclim_global(var = "bio", res = 0.5, path = path_clim)[[5]] |>
  crop(ext(-10.41, 54.58, 30.91667, 72.75))

names(bio5_crop) <- "bio5"
writeRaster(bio5_crop, file.path(path_clim, "bio5_crop.tif"), overwrite = TRUE)

# GBIF records
species_list <- c(
  "Apis mellifera",
  "Bombus lapidarius",
  "Bombus terrestris",
  "Eupeodes corollae",
  "Episyrphus balteatus",
  "Aglais urticae",
  "Vespula vulgaris",
  "Eristalis tenax"
)

occ <- fread(gbif_file, sep = "\t", quote = "", na.strings = c("", "NA"))
occ[, `:=`(
  year = as.integer(year),
  decimalLongitude = as.numeric(decimalLongitude),
  decimalLatitude = as.numeric(decimalLatitude),
  coordinateUncertaintyInMeters = as.numeric(coordinateUncertaintyInMeters)
)]

for (species_name in species_list) {
  gbif_clean <- occ[species == species_name] |>
    mutate(lon = decimalLongitude, lat = decimalLatitude) |>
    distinct(lon, lat, .keep_all = TRUE) |>
    distinct(occurrenceID, .keep_all = TRUE) |>
    filter(
      !is.na(lon),
      !is.na(lat),
      !is.na(year),
      year >= 1970,
      !(lon == 0 & lat == 0),
      is.na(identificationVerificationStatus) | identificationVerificationStatus != "casual",
      taxonomicStatus == "ACCEPTED",
      is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters <= 5000
    )

  fwrite(
    as.data.table(gbif_clean),
    file.path(path_gbif, paste0(str_replace_all(species_name, " ", "_"), "_clean.csv"))
  )
}
