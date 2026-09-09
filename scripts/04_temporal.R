library(tidyverse)
library(glmmTMB)
library(patchwork)
library(DHARMa)

path_ana  <- "D:/R_P1/r_data"
path_data <- "D:/R_P1/raw_data/tereno"
path_res  <- "D:/R_P1/results"

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

# TERENO and temperature data
tmax_orig <- read.csv(file.path(path_data, "dwd_tmax.csv"), row.names = NULL)
names(tmax_orig) <- c(
  "Produkt_Code",
  "station",
  "day",
  "temp",
  "Qualitaet_Byte",
  "Qualitaet_Niveau"
)

tereno <- read.csv(file.path(path_data, "tereno_filtered.csv")) |>
  mutate(
    StartDate = as.Date(StartDate),
    EndDate = as.Date(EndDate)
  )

# Camera models
fit_camera_model <- function(class_number) {
  model_data <- cam_day_model |>
    filter(
      class == class_number,
      is.finite(temp),
      !is.na(abun),
      n_systems > 0
    )

  glmmTMB(
    abun ~ scale(temp) + scale(I(temp^2)) + offset(log(n_systems)),
    data = model_data,
    family = nbinom2
  )
}

fm_bt <- fit_camera_model(4)
fm_bl <- fit_camera_model(3)

# Daily camera predictions
pred_daily <- tmax_orig |>
  dplyr::select(station, day, temp) |>
  mutate(
    station = as.character(station),
    day = as.Date(day),
    n_systems = 1
  ) |>
  filter(
    day >= as.Date("2010-05-01"),
    day <= as.Date("2020-08-17")
  )

pred_daily$pred_cam_bt <- predict(fm_bt, newdata = pred_daily, type = "response")
pred_daily$pred_cam_bl <- predict(fm_bl, newdata = pred_daily, type = "response")

site_station <- tibble(
  site = c("FBG", "GFH", "HAR", "SIP", "SST", "WAN"),
  station = c("2704", "3226", "4032", "2044", "4036", "3126")
)

pred_daily <- pred_daily |>
  left_join(site_station, by = "station")

tereno_bt <- tereno |>
  filter(GenSpec == "Bombus_terrestris")

tereno_bl <- tereno |>
  filter(GenSpec == "Bombus_lapidarius")

# Sampling-period predictions
aggregate_predictions <- function(trap_data, prediction_column) {
  trap_data |>
    rowwise() |>
    mutate(
      sum_cam_pred = sum(
        pred_daily[[prediction_column]][
          pred_daily$site == site &
            pred_daily$day >= StartDate &
            pred_daily$day <= EndDate
        ],
        na.rm = TRUE
      )
    ) |>
    ungroup()
}

aggregated_pred_bt <- aggregate_predictions(tereno_bt, "pred_cam_bt")
aggregated_pred_bl <- aggregate_predictions(tereno_bl, "pred_cam_bl")

save(
  aggregated_pred_bt,
  file = file.path(path_ana, "aggregated_pred_bt_mean_tmax.RData")
)
save(
  aggregated_pred_bl,
  file = file.path(path_ana, "aggregated_pred_bl_mean_tmax.RData")
)

# Annual validation models
aggregate_yearly <- function(data) {
  data |>
    group_by(site, YearValue) |>
    summarise(
      sum_pantrap = sum(sum_pantrap, na.rm = TRUE),
      sum_cam_pred = sum(sum_cam_pred, na.rm = TRUE),
      .groups = "drop"
    )
}

bt_yearly <- aggregate_yearly(aggregated_pred_bt)
bl_yearly <- aggregate_yearly(aggregated_pred_bl)

model_bt <- glmmTMB(
  sum_pantrap ~ sum_cam_pred + (1 | site),
  data = bt_yearly,
  family = nbinom2
)
model_bl <- glmmTMB(
  sum_pantrap ~ sum_cam_pred + (1 | site),
  data = bl_yearly,
  family = nbinom2
)

sim_bt <- simulateResiduals(model_bt)
plot(sim_bt)

sim_bl <- simulateResiduals(model_bl)
plot(sim_bl)

summary(model_bt)
summary(model_bl)

extract_glmm_performance <- function(model, predictor_name) {
  fe <- broom.mixed::tidy(
    model,
    effects = "fixed",
    conf.int = TRUE,
    exponentiate = FALSE
  ) %>%
    filter(term == predictor_name) %>%
    mutate(
      IRR      = exp(estimate),
      IRR_low  = exp(conf.low),
      IRR_high = exp(conf.high)
    ) %>%
    select(
      term,
      beta = estimate,
      conf.low,
      conf.high,
      IRR,
      IRR_low,
      IRR_high,
      p.value
    )
  
  r2 <- performance::r2_nakagawa(model)
  
  list(
    fixed_effect = fe,
    r2 = r2
  )
}

perf_bt <- extract_glmm_performance(
  model_bt,
  predictor_name = "sum_cam_pred"
)

perf_bt$fixed_effect
perf_bt$r2

perf_bl <- extract_glmm_performance(
  model_bl,
  predictor_name = "sum_cam_pred"
)

perf_bl$fixed_effect
perf_bl$r2

# Trend figures
plot_trends_z <- function(data, output_prefix) {
  plot_data <- data |>
    mutate(
      YearValue = as.integer(YearValue),
      site = factor(site),
      sum_pantrap = as.numeric(sum_pantrap),
      sum_cam_pred = as.numeric(sum_cam_pred)
    ) |>
    group_by(site) |>
    mutate(
      pantrap_z = as.numeric(scale(sum_pantrap)),
      camera_z = as.numeric(scale(sum_cam_pred))
    ) |>
    ungroup() |>
    pivot_longer(
      c(pantrap_z, camera_z),
      names_to = "series",
      values_to = "value"
    ) |>
    mutate(
      series = recode(
        series,
        pantrap_z = "Observed pan-trap catches (annual)",
        camera_z = "Predicted population size (annual)"
      ),
      series = factor(
        series,
        levels = c(
          "Observed pan-trap catches (annual)",
          "Predicted population size (annual)"
        )
      )
    )

  plot_object <- ggplot(
    plot_data,
    aes(
      x = YearValue,
      y = value,
      colour = series,
      linetype = series,
      group = series
    )
  ) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
    geom_line(linewidth = 0.7, lineend = "round") +
    geom_point(size = 1.4) +
    facet_wrap(~ site, ncol = 3) +
    scale_colour_manual(values = c(
      "Observed pan-trap catches (annual)" = "#0072B2",
      "Predicted population size (annual)" = "#D55E00"
    )) +
    scale_linetype_manual(values = c(
      "Observed pan-trap catches (annual)" = "solid",
      "Predicted population size (annual)" = "22"
    )) +
    scale_x_continuous(
      breaks = seq(
        min(plot_data$YearValue, na.rm = TRUE),
        max(plot_data$YearValue, na.rm = TRUE),
        by = 2
      )
    ) +
    labs(
      x = "Year",
      y = "Standardised signal (z-score)",
      colour = NULL,
      linetype = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11, margin = margin(b = 6)),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 9),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      panel.spacing = unit(8, "pt"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = 9),
      legend.key.width = unit(18, "pt")
    ) +
    guides(
      colour = guide_legend(nrow = 1, byrow = TRUE),
      linetype = guide_legend(nrow = 1, byrow = TRUE)
    )

  ggsave(
    file.path(path_res, paste0(output_prefix, "_trend_z.pdf")),
    plot_object,
    width = 180,
    height = 120,
    units = "mm",
    device = cairo_pdf
  )
  ggsave(
    file.path(path_res, paste0(output_prefix, "_trend_z.tiff")),
    plot_object,
    width = 180,
    height = 120,
    units = "mm",
    dpi = 600,
    compression = "lzw"
  )

  plot_object
}

# Effect figures
plot_effect_glmm <- function(data, output_prefix) {
  predictor_mean <- mean(data$sum_cam_pred, na.rm = TRUE)
  predictor_sd <- sd(data$sum_cam_pred, na.rm = TRUE)

  model_data <- data |>
    mutate(sum_cam_pred_sc = (sum_cam_pred - predictor_mean) / predictor_sd)

  model <- glmmTMB(
    sum_pantrap ~ sum_cam_pred_sc + (1 | site),
    data = model_data,
    family = nbinom2
  )

  new_data <- tibble(
    sum_cam_pred = seq(
      min(data$sum_cam_pred, na.rm = TRUE),
      max(data$sum_cam_pred, na.rm = TRUE),
      length.out = 200
    ),
    site = NA
  ) |>
    mutate(sum_cam_pred_sc = (sum_cam_pred - predictor_mean) / predictor_sd)

  prediction <- predict(
    model,
    newdata = new_data,
    type = "link",
    se.fit = TRUE,
    allow.new.levels = TRUE
  )

  new_data <- new_data |>
    mutate(
      fit_link = as.numeric(prediction$fit),
      se_link = as.numeric(prediction$se.fit),
      fit = exp(fit_link),
      lwr = exp(fit_link - 1.96 * se_link),
      upr = exp(fit_link + 1.96 * se_link)
    )

  plot_object <- ggplot() +
    geom_ribbon(
      data = new_data,
      aes(sum_cam_pred, ymin = lwr, ymax = upr),
      fill = "#56B4E9",
      alpha = 0.25
    ) +
    geom_line(
      data = new_data,
      aes(sum_cam_pred, fit),
      linewidth = 0.9,
      colour = "#0072B2",
      lineend = "round"
    ) +
    geom_point(
      data = data,
      aes(sum_cam_pred, sum_pantrap),
      shape = 21,
      size = 2,
      stroke = 0.3,
      fill = "#E69F00",
      colour = "black",
      alpha = 0.85
    ) +
    labs(
      x = "Predicted population size (annual)",
      y = "Observed pan-trap catches (annual)"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11, margin = margin(b = 6)),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9)
    )

  ggsave(
    file.path(path_res, paste0(output_prefix, "_effect.pdf")),
    plot_object,
    width = 180,
    height = 130,
    units = "mm",
    device = cairo_pdf
  )
  ggsave(
    file.path(path_res, paste0(output_prefix, "_effect.tiff")),
    plot_object,
    width = 180,
    height = 130,
    units = "mm",
    dpi = 600,
    compression = "lzw"
  )

  plot_object
}

# Validation figures
p_bt_trend <- plot_trends_z(bt_yearly, "Fig_BT")
p_bl_trend <- plot_trends_z(bl_yearly, "Fig_BL")
p_bt_eff <- plot_effect_glmm(bt_yearly, "Fig_BT")
p_bl_eff <- plot_effect_glmm(bl_yearly, "Fig_BL")

keep_facet_site <- function(plot_object, selected_site = "FBG") {
  plot_object %+% filter(plot_object$data, site == selected_site) +
    facet_wrap(~ site, ncol = 1)
}

p_bl_fbg <- keep_facet_site(p_bl_trend) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Year",
    y = "Standardised signal (z-score)"
  ) +
  scale_x_continuous(limits = c(2010, 2020), breaks = seq(2010, 2020, 2)) +
  theme(
    plot.title = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

p_bt_fbg <- keep_facet_site(p_bt_trend) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Year",
    y = "Standardised signal (z-score)"
  ) +
  scale_x_continuous(limits = c(2010, 2020), breaks = seq(2010, 2020, 2)) +
  theme(
    plot.title = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

p_bl_eff2 <- p_bl_eff +
  labs(
    title = NULL,
    subtitle = expression(italic("Bombus lapidarius")),
    x = "Predicted population size (annual)",
    y = "Observed pan-trap catches (annual)"
  ) +
  theme(plot.title = element_blank())

p_bt_eff2 <- p_bt_eff +
  labs(
    title = NULL,
    subtitle = expression(italic("Bombus terrestris")),
    x = "Predicted population size (annual)",
    y = "Observed pan-trap catches (annual)"
  ) +
  theme(plot.title = element_blank())

p_combo <- ((p_bl_eff2 | p_bl_fbg) / (p_bt_eff2 | p_bt_fbg)) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.margin = margin(6, 10, 6, 6),
    plot.subtitle = element_text(size = 10, face = "italic", margin = margin(b = 4))
  )

ggsave(
  file.path(path_res, "Fig_Bombus_validation_COMBINED.pdf"),
  p_combo,
  width = 180,
  height = 190,
  units = "mm",
  device = cairo_pdf
)
ggsave(
  file.path(path_res, "Fig_Bombus_validation_COMBINED.tiff"),
  p_combo,
  width = 180,
  height = 190,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

# Temporal skill table
site_temporal_skill <- function(data) {
  data |>
    mutate(site = factor(site)) |>
    group_by(site) |>
    mutate(
      pantrap_z = as.numeric(scale(sum_pantrap)),
      camera_z = as.numeric(scale(sum_cam_pred))
    ) |>
    summarise(
      n = sum(is.finite(pantrap_z) & is.finite(camera_z)),
      r = cor(pantrap_z, camera_z, use = "complete.obs", method = "spearman"),
      rmse_z = sqrt(mean((pantrap_z - camera_z)^2, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    mutate(
      r2 = r^2,
      z = atanh(r),
      se = 1 / sqrt(pmax(n - 3, 1)),
      r_lwr = tanh(z - 1.96 * se),
      r_upr = tanh(z + 1.96 * se)
    ) |>
    dplyr::select(site, n, r, r_lwr, r_upr, r2, rmse_z) |>
    arrange(desc(r))
}

tab_bt <- site_temporal_skill(bt_yearly)
tab_bl <- site_temporal_skill(bl_yearly)

tab_wide <- bind_rows(
  transmute(
    tab_bl,
    Site = as.character(site),
    Species = "Bombus lapidarius (ρ²)",
    rho2 = r2
  ),
  transmute(
    tab_bt,
    Site = as.character(site),
    Species = "Bombus terrestris (ρ²)",
    rho2 = r2
  )
) |>
  mutate(
    rho2 = if_else(
      abs(rho2) < 0.001,
      "<0.001",
      sprintf("%.3f", rho2)
    )
  ) |>
  pivot_wider(
    names_from = Species,
    values_from = rho2
  ) |>
  arrange(Site)

write.csv(
  tab_wide,
  file.path(path_res, "Table_site_temporal_skill_rho2.csv"),
  row.names = FALSE
)

# Site-level effect model
r2_long <- bind_rows(
  transmute(
    tab_bt,
    site,
    species = "Bombus terrestris",
    r2
  ),
  transmute(
    tab_bl,
    site,
    species = "Bombus lapidarius",
    r2
  )
)

habitat_vars <- c(
  "EUNIS_div_agg",
  "arable",
  "forest",
  "managed_grassland",
  "semi_natural",
  "urban",
  "crop_div",
  "crop_rich",
  "mass_fl_area",
  "mass_fl_pa",
  "EUNIS_rich_agg",
  "EUNIS_div",
  "EUNIS_rich"
)

land_by_site <- read.csv(
  file.path(path_data, "landuse.csv"),
  row.names = NULL
) |>
  mutate(
    site = str_extract(
      Name,
      "^(FBG|GFH|HAR|SIP|SST|WAN)"
    )
  ) |>
  filter(!is.na(site)) |>
  group_by(site) |>
  summarise(
    across(
      all_of(habitat_vars),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

# Habitat-variable clustering
habitat_data <- land_by_site |>
  ungroup() |>
  dplyr::select(
    all_of(habitat_vars)
  )

cor_matrix <- cor(
  habitat_data,
  method = "spearman",
  use = "pairwise.complete.obs"
)

var_cluster <- hclust(
  as.dist(1 - abs(cor_matrix)),
  method = "complete"
)

plot(
  var_cluster,
  hang = -1,
  xlab = "",
  sub = "",
  ylab = expression(1 - abs(rho))
)

abline(
  h = 0.3,
  lty = 2
)

variable_clusters <- cutree(
  var_cluster,
  h = 0.3
)

cluster_table <- tibble(
  variable = names(variable_clusters),
  cluster = as.integer(variable_clusters)
) |>
  arrange(cluster, variable)

cluster_table

write.csv(
  cluster_table,
  file.path(
    path_res,
    "Table_habitat_variable_clusters.csv"
  ),
  row.names = FALSE
)

model_data <- r2_long |>
  left_join(
    land_by_site,
    by = "site"
  ) |>
  mutate(
    species = factor(species)
  )

# Backward model selection
m_full <- glmmTMB(
  r2 ~
    crop_div +
    semi_natural +
    EUNIS_rich_agg +
    mass_fl_area +
    (1 | species),
  data = model_data,
  family = beta_family(link = "logit")
)

summary(m_full)
drop1(m_full, test = "Chisq")

m_1 <- update(
  m_full,
  . ~ . - EUNIS_rich_agg
)

summary(m_1)
drop1(m_1, test = "Chisq")

m_2 <- update(
  m_1,
  . ~ . - mass_fl_area
)

summary(m_2)
drop1(m_2, test = "Chisq")

m_3 <- update(
  m_2,
  . ~ . - semi_natural
)

summary(m_3)

m <- m_2

sim_m <- simulateResiduals(m, n = 1000)
plot(sim_m)

summary(m)

# Site-level effect figure
make_effect_plot <- function(model, data, focal_variable, x_label) {
  new_data <- tibble(
    crop_div = mean(data$crop_div, na.rm = TRUE),
    semi_natural = mean(data$semi_natural, na.rm = TRUE),
    species = levels(data$species)[1]
  ) |>
    slice(rep(1, 200))

  new_data[[focal_variable]] <- seq(
    min(data[[focal_variable]], na.rm = TRUE),
    max(data[[focal_variable]], na.rm = TRUE),
    length.out = 200
  )

  prediction <- predict(
    model,
    newdata = new_data,
    type = "link",
    se.fit = TRUE,
    re.form = NA
  )

  new_data <- new_data |>
    mutate(
      fit = plogis(prediction$fit),
      lwr = plogis(
        prediction$fit -
          1.96 * prediction$se.fit
      ),
      upr = plogis(
        prediction$fit +
          1.96 * prediction$se.fit
      )
    )

  ggplot() +
    geom_ribbon(
      data = new_data,
      aes(
        x = .data[[focal_variable]],
        ymin = lwr,
        ymax = upr
      ),
      fill = "grey70",
      alpha = 0.35
    ) +
    geom_line(
      data = new_data,
      aes(
        x = .data[[focal_variable]],
        y = fit
      ),
      linewidth = 0.9,
      colour = "#2C7BB6"
    ) +
    geom_point(
      data = data,
      aes(
        x = .data[[focal_variable]],
        y = r2
      ),
      shape = 21,
      size = 2.2,
      stroke = 0.4,
      fill = NA,
      colour = "black"
    ) +
    labs(
      x = x_label,
      y = expression(rho^2)
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_blank())
}

p_crop <- make_effect_plot(
  m,
  model_data,
  "crop_div",
  "Crop diversity [Shannon]"
)

p_semi <- make_effect_plot(
  m,
  model_data,
  "semi_natural",
  "Semi-natural habitat proportion"
)

p_effects <- p_crop + p_semi +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 13
    ),
    plot.tag.position = c(0.02, 0.98)
  )

ggsave(
  file.path(
    path_res,
    "Fig_rho2_effects_crop_seminatural.png"
  ),
  p_effects,
  width = 180,
  height = 90,
  units = "mm",
  dpi = 600
)

