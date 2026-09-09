# README

This repository contains the R scripts used to process camera-trap detections, fit temperature-response models, generate camera-based species distribution models, and validate temporal predictions against independent pan-trap data.

## Running the scripts

Run the scripts in numerical order. Before starting each script, restart R and clear the workspace, for example with "rm(list = ls())". The scripts assume the local folder structure defined in the "path_*" objects at the top of each file.

## Scripts

- "01_data_prep.R": Imports and cleans camera detections, joins temperature data, creates hourly and daily camera datasets, crops WorldClim BIO5 data, and filters GBIF records.
- "02_activity_models.R": Fits species-specific negative-binomial temperature-activity models, generates prediction curves, model diagnostics, figures, and model-comparison tables.
- "03_SDMs.R": Predicts pollinator activity onto the BIO5 rasters using the camera-based temperature-activity models, saves species-specific prediction rasters and maps, overlays GBIF occurrences, and calculates Boyce index values.
- "04_temporal.R": Validates predictions by camera-based temperature-activity models against TERENO pan-trap data for *Bombus terrestris* and *Bombus lapidarius*, and produces validation figures, temporal skill tables, and habitat-effect summaries.

## Main outputs

Processed data are saved to "r_data/". Figures, tables, prediction rasters, model summaries, and validation outputs are saved to "results/".

