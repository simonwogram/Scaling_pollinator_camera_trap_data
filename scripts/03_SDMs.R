library(tidyverse)
library(glmmTMB)
library(raster)
library(terra)
library(modEvA)
library(ggnewscale)

path_ana  <- "D:/R_P1/r_data"
path_clim <- "D:/R_P1/raw_data/climate"
path_gbif <- "D:/R_P1/raw_data/gbif"
path_res  <- "D:/R_P1/results"
path_pred_rasters <- file.path(path_res, "prediction_rasters")

dir.create(path_pred_rasters, recursive = TRUE, showWarnings = FALSE)

# Species metadata
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

species_classes <- setNames(2:9, species_list)

# Input data
tmax_crop <- raster(file.path(path_clim, "bio5_crop.tif"))
names(tmax_crop) <- "temp"
load(file.path(path_ana, "cam_day_tmax.RData"))

# Camera model data
daily_effort <- cam_day |>
  distinct(date, SystemID) |>
  count(date, name = "n_systems")

cam_day_model <- cam_day |>
  group_by(date, class) |>
  summarise(
    temp = max(temp, na.rm = TRUE),
    abun = sum(abun, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(daily_effort, by = "date")

# Raster scaling
scale_raster_01 <- function(raster_object) {
  raster_min <- global(raster_object, "min", na.rm = TRUE)[1, 1]
  raster_max <- global(raster_object, "max", na.rm = TRUE)[1, 1]
  (raster_object - raster_min) / (raster_max - raster_min)
}

# Species distribution models
boyce_results <- tibble(
  species = species_list,
  class = as.integer(species_classes),
  boyce_index = NA_real_
)

pred_rasters <- setNames(vector("list", length(species_list)), species_list)

for (i in seq_along(species_list)) {
  species_name <- species_list[i]
  class_number <- species_classes[i]
  file_name <- str_replace_all(species_name, " ", "_")

  model_data <- cam_day_model |>
    filter(
      class == class_number,
      is.finite(temp),
      !is.na(abun),
      n_systems > 0
    )

  model <- glmmTMB(
    abun ~ scale(temp) + scale(I(temp^2)) + offset(log(n_systems)),
    data = model_data,
    family = nbinom2
  )

  pred_1km <- raster::predict(
    object = tmax_crop,
    model = model,
    type = "response",
    const = data.frame(n_systems = 1)
  ) |>
    rast()

  pred_10km <- aggregate(pred_1km, fact = 10)
  pred_10km_scaled <- scale_raster_01(pred_10km)
  pred_rasters[[species_name]] <- pred_10km

  species_dir <- file.path(path_res, file_name)
  raster_dir <- file.path(path_pred_rasters, file_name)
  dir.create(species_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)

  writeRaster(
    pred_1km,
    file.path(raster_dir, paste0(file_name, "_camera_raw_1km.tif")),
    overwrite = TRUE
  )
  writeRaster(
    pred_10km,
    file.path(raster_dir, paste0(file_name, "_camera_raw_10km.tif")),
    overwrite = TRUE
  )
  writeRaster(
    pred_10km_scaled,
    file.path(raster_dir, paste0(file_name, "_camera_scaled_10km.tif")),
    overwrite = TRUE
  )
  writeRaster(
    pred_10km,
    file.path(species_dir, paste0("pred_", file_name, ".tif")),
    overwrite = TRUE
  )
  writeRaster(
    pred_10km_scaled,
    file.path(species_dir, paste0("pred_scaled_", file_name, ".tif")),
    overwrite = TRUE
  )

  occurrence_data <- read.csv(
    file.path(path_gbif, paste0(file_name, "_clean.csv"))
  )
  occurrence_points <- occurrence_data[, c("lon", "lat")]

  presence_raster <- rast(pred_10km)
  values(presence_raster) <- NA
  presence_cells <- unique(
    cellFromXY(presence_raster, as.matrix(occurrence_points))
  )
  presence_cells <- presence_cells[!is.na(presence_cells)]
  presence_raster[presence_cells] <- 1

  png(
    file.path(species_dir, paste0("pred_", file_name, ".png")),
    width = 2400,
    height = 2000,
    res = 300
  )
  plot(pred_10km, main = species_name)
  dev.off()

  png(
    file.path(species_dir, paste0("pred_presences_", file_name, ".png")),
    width = 2400,
    height = 2000,
    res = 300
  )
  plot(pred_10km, main = species_name)
  plot(presence_raster, add = TRUE, col = "black", legend = FALSE)
  dev.off()

  png(
    file.path(species_dir, paste0("boyce_", file_name, ".png")),
    width = 2000,
    height = 2000,
    res = 300
  )
  boyce <- Boyce(
    obs = occurrence_points,
    pred = pred_10km,
    main = species_name,
    rm.dup.points = TRUE
  )
  dev.off()

  boyce_results$boyce_index[i] <- boyce$Boyce
}

# Boyce results
boyce_table <- boyce_results |>
  transmute(
    Species = paste0("*", species, "*"),
    Class = class,
    Boyce = round(boyce_index, 3)
  )

write.csv(
  boyce_results,
  file.path(path_res, "camera_sdm_boyce_results.csv"),
  row.names = FALSE
)
write.csv(
  boyce_table,
  file.path(path_res, "camera_sdm_boyce_results_formatted.csv"),
  row.names = FALSE
)

# Multi-species map data
main_species <- c(
  "Aglais urticae",
  "Apis mellifera",
  "Bombus lapidarius",
  "Bombus terrestris",
  "Episyrphus balteatus",
  "Eristalis tenax",
  "Eupeodes corollae",
  "Vespula vulgaris"
)

map_df <- imap_dfr(pred_rasters, ~ {
  species_raster <- .x
  species_name <- .y
  as.data.frame(species_raster, xy = TRUE, na.rm = FALSE) |>
    setNames(c("x", "y", "value")) |>
    mutate(species = species_name)
}) |>
  filter(species %in% main_species) |>
  group_by(species) |>
  mutate(
    value_rel = (value - min(value, na.rm = TRUE)) /
      (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))
  ) |>
  ungroup()

species_labels <- setNames(
  paste0("italic(", str_replace(main_species, " ", "~"), ")"),
  main_species
)

# GBIF occurrence cells
gbif_cells_df <- map_dfr(main_species, ~ {
  species_name <- .x
  file_name <- str_replace_all(species_name, " ", "_")
  occurrence_points <- read.csv(
    file.path(path_gbif, paste0(file_name, "_clean.csv"))
  ) |>
    dplyr::select(lon, lat) |>
    filter(is.finite(lon), is.finite(lat))

  species_raster <- pred_rasters[[species_name]]
  cells <- unique(cellFromXY(species_raster, as.matrix(occurrence_points)))
  cells <- cells[!is.na(cells)]
  cell_centres <- xyFromCell(species_raster, cells)

  tibble(
    x = cell_centres[, 1],
    y = cell_centres[, 2],
    species = species_name,
    layer = "GBIF occurrence"
  )
})

reference_raster <- pred_rasters[[main_species[1]]]
cell_width <- res(reference_raster)[1]
cell_height <- res(reference_raster)[2]

boyce_label_df <- map_df |>
  group_by(species) |>
  summarise(
    x = min(x, na.rm = TRUE) + 2,
    y = max(y, na.rm = TRUE) - 2,
    .groups = "drop"
  ) |>
  left_join(dplyr::select(boyce_results, species, boyce_index), by = "species") |>
  mutate(
    label = if_else(
      is.na(boyce_index),
      "Boyce = NA",
      sprintf("Boyce = %.3f", boyce_index)
    )
  )

# Multi-species map
p_map_main <- ggplot() +
  geom_tile(data = map_df, aes(x, y, fill = value_rel)) +
  scale_fill_viridis_c(
    name = "Predicted climatic\nsuitability (scaled)",
    limits = c(0, 1),
    na.value = "white",
    guide = guide_colorbar(
      order = 1,
      title.position = "top",
      title.hjust = 0.5,
      label.position = "bottom",
      barwidth = unit(45, "mm"),
      barheight = unit(4, "mm")
    )
  ) +
  new_scale_fill() +
  geom_tile(
    data = gbif_cells_df,
    aes(x, y, fill = layer),
    width = cell_width,
    height = cell_height,
    alpha = 0.60,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("GBIF occurrence" = "red"),
    guide = guide_legend(
      order = 2,
      override.aes = list(alpha = 0.80)
    )
  ) +
  geom_label(
    data = boyce_label_df,
    aes(x, y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3,
    label.size = 0,
    fill = scales::alpha("white", 0.65)
  ) +
  coord_fixed(expand = FALSE) +
  facet_wrap(
    ~ species,
    ncol = 2,
    labeller = as_labeller(species_labels, default = label_parsed)
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.height = unit(5, "mm"),
    legend.key.width = unit(8, "mm"),
    panel.spacing = unit(12, "pt"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    plot.margin = margin(6, 6, 6, 6)
  )

ggsave(
  file.path(path_res, "Fig_SDM_multispecies_8sp_MAIN.png"),
  p_map_main,
  width = 150,
  height = 300,
  units = "mm",
  dpi = 600,
  bg = "white"
)
