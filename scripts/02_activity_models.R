library(tidyverse)
library(glmmTMB)
library(gt)
library(DHARMa)
library(performance)

path_ana <- "D:/R_P1/r_data"
path_res <- "D:/R_P1/results"

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

# Camera model data
load(file.path(path_ana, "cam_day_tmax.RData"))

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

# Temperature-response models
model_results <- map2(species_list, species_classes, ~ {
  species_name <- .x
  class_number <- .y
  model_data <- cam_day_model |>
    filter(
      class == class_number,
      is.finite(temp),
      !is.na(abun),
      n_systems > 0
    )

  full_model <- glmmTMB(
    abun ~ scale(temp) + scale(I(temp^2)) + offset(log(n_systems)),
    data = model_data,
    family = nbinom2
  )
  
  conv <- performance::check_convergence(full_model)
  print(conv)
  
  simres <- simulateResiduals(full_model, n = 1000)
  plot(simres)

  null_model <- glmmTMB(
    abun ~ 1 + offset(log(n_systems)),
    data = model_data,
    family = nbinom2
  )

  list(
    data = model_data,
    model = full_model,
    null_model = null_model,
    LRT = anova(null_model, full_model),
    delta_AIC = AIC(null_model) - AIC(full_model),
    coefficients = summary(full_model)$coefficients$cond
  )
}) |>
  setNames(species_list)

# Temperature-response predictions
pred_df <- imap_dfr(model_results, ~ {
  result <- .x
  species_name <- .y
  new_data <- tibble(
    temp = seq(
      min(result$data$temp, na.rm = TRUE),
      max(result$data$temp, na.rm = TRUE),
      length.out = 200
    ),
    n_systems = 1
  )

  prediction <- predict(
    result$model,
    newdata = new_data,
    type = "response",
    se.fit = TRUE
  )

  new_data |>
    mutate(
      fit = as.numeric(prediction$fit),
      se = as.numeric(prediction$se.fit),
      lwr = pmax(0, fit - 1.96 * se),
      upr = fit + 1.96 * se,
      species = species_name
    )
})

raw_df <- imap_dfr(model_results, ~ {
  result <- .x
  species_name <- .y
  result$data |>
    transmute(
      temp,
      abun = abun / n_systems,
      species = species_name
    )
})

# Temperature-response figure
species_labels <- setNames(
  paste0("italic(", str_replace(species_list, " ", "~"), ")"),
  species_list
)

p_temp_eff <- ggplot() +
  geom_ribbon(
    data = pred_df,
    aes(temp, ymin = lwr, ymax = upr),
    alpha = 0.25
  ) +
  geom_line(
    data = pred_df,
    aes(temp, fit),
    linewidth = 0.9,
    lineend = "round"
  ) +
  geom_point(
    data = raw_df,
    aes(temp, abun),
    size = 0.9,
    alpha = 0.35
  ) +
  facet_wrap(
    ~ species,
    scales = "free_y",
    ncol = 4,
    labeller = as_labeller(species_labels, default = label_parsed)
  ) +
  labs(
    x = "Daily maximum temperature (°C)",
    y = "Daily camera-trap detections per active system (count)",
    title = "Temperature–response curves for daily camera-trap detections"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11, margin = margin(b = 6)),
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    panel.spacing = unit(8, "pt")
  )

ggsave(
  file.path(path_res, "Fig_temp_response_facets.pdf"),
  p_temp_eff,
  width = 202,
  height = 140,
  units = "mm",
  device = cairo_pdf
)

ggsave(
  file.path(path_res, "Fig_temp_response_facets.tiff"),
  p_temp_eff,
  width = 202,
  height = 140,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

# Model comparison table
model_table <- imap_dfr(model_results, ~ {
  result <- .x
  species_name <- .y
  lrt <- result$LRT
  table_null_model <- glmmTMB(abun ~ 1 + offset(log(n_systems)), data = result$data, family = nbinom2)

  tibble(
    Species = species_name,
    AIC_null = AIC(table_null_model),
    AIC_temp = AIC(result$model),
    `ΔAIC` = AIC_null - AIC_temp,
    `LRT χ²` = as.numeric(lrt$Chisq[2]),
    df = as.integer(lrt$Df[2]),
    p = as.numeric(lrt$`Pr(>Chisq)`[2])
  )
}) |>
  mutate(
    across(c(AIC_null, AIC_temp, `ΔAIC`, `LRT χ²`), ~ round(.x, 2)),
    p = if_else(p < 0.001, "< 0.001", sprintf("%.3f", p))
  )

model_gt <- model_table |>
  mutate(Species = paste0("*", Species, "*")) |>
  gt() |>
  fmt_markdown(columns = Species)

gtsave(
  model_gt,
  filename = file.path(path_res, "Table_AIC_LRT_temperature_models.html")
)
