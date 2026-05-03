#study I continuous variables analysis and plots 


# ---- 1. GROWTH GLMMs — ALL PLANTS ------------------------------------------
# ALL terms printed (no p-value filter)
# day centred; unique_id random intercept; cycle fixed block

growth_traits <- c("plant_height_cm", "total_tillers", "leaf_no_ms",
                   "dev_stage", "SPAD", "dry_biomass", "gwc")

cat("\n=== GROWTH TRAJECTORY GLMMs — ALL PLANTS ===\n")

growth_models <- map(growth_traits, function(trait) {
  df <- long %>%
    filter(!is.na(.data[[trait]])) %>%
    mutate(day_c = as.numeric(scale(day, center = TRUE, scale = FALSE)))
  
  frm <- as.formula(paste(
    trait,
    "~ cultivar * aphid_treatment * watering_treatment * day_c",
    "+ cycle + (1 | unique_id)"
  ))
  
  mod <- tryCatch(
    lmer(frm, data = df, REML = FALSE),
    error = function(e) {
      message("  Full model failed for ", trait, ": ", e$message)
      frm2 <- as.formula(paste(
        trait,
        "~ (cultivar + aphid_treatment + watering_treatment) * day_c",
        "+ cultivar:aphid_treatment + cultivar:watering_treatment",
        "+ aphid_treatment:watering_treatment + cycle + (1 | unique_id)"
      ))
      tryCatch(lmer(frm2, data = df, REML = FALSE), error = function(e2) NULL)
    }
  )
  
  if (!is.null(mod)) {
    cat("\n---", trait, "---\n")
    print(anova(mod, type = "III"))   # ALL terms printed
  }
  mod
})
names(growth_models) <- growth_traits


# ---- 2. GROWTH GLMMs — LAUREATE ONLY ---------------------------------------

cat("\n=== GROWTH TRAJECTORY GLMMs — LAUREATE ONLY ===\n")

growth_models_laur <- map(growth_traits, function(trait) {
  df <- long_laur %>%
    filter(!is.na(.data[[trait]])) %>%
    mutate(day_c = as.numeric(scale(day, center = TRUE, scale = FALSE)))
  
  frm <- as.formula(paste(
    trait,
    "~ aphid_treatment * watering_treatment * day_c",
    "+ cycle + (1 | unique_id)"
  ))
  
  mod <- tryCatch(
    lmer(frm, data = df, REML = FALSE),
    error = function(e) {
      message("  Laureate GLMM failed for ", trait, ": ", e$message); NULL
    }
  )
  
  if (!is.null(mod)) {
    cat("\n---", trait, "(laureate only) ---\n")
    print(anova(mod, type = "III"))
  }
  mod
})
names(growth_models_laur) <- growth_traits


# ---- 3. GROWTH GLMMs — CYCLE 2 ONLY ---------------------------------------

cat("\n=== GROWTH TRAJECTORY GLMMs — CYCLE 2 ONLY ===\n")

growth_models_2 <- map(growth_traits, function(trait) {
  df <- long_2 %>%
    filter(!is.na(.data[[trait]])) %>%
    mutate(day_c = as.numeric(scale(day, center = TRUE, scale = FALSE)))
  
  frm <- as.formula(paste(
    trait,
    "~ aphid_treatment * watering_treatment * cultivar * day_c",
    " + (1 | unique_id)"
  ))
  
  mod <- tryCatch(
    lmer(frm, data = df, REML = FALSE),
    error = function(e) {
      message("  cycle 2 GLMM failed for ", trait, ": ", e$message); NULL
    }
  )
  
  if (!is.null(mod)) {
    cat("\n---", trait, "(cycle 2 only) ---\n")
    print(anova(mod, type = "III"))
  }
  mod
})
names(growth_models_2) <- growth_traits


# ---- 2. GROWTH TRAJECTORY PLOTS -------------------------------------------
# Aesthetics:
#   colour    = watering_treatment (steelblue / tomato)
#   linetype  = cultivar (solid bere / dashed laureate)
#   shape     = aphid_treatment  (1 = open circle no-aphids; 17 = filled triangle aphids)
#   linewidth = 0.4  (thin)
#   error bars: alpha = 0.5, position_dodge
#   position  = dodge to separate overlapping groups

plot_growth_traj <- function(trait, ylab, data = long, show_facet = TRUE,
                             include_cultivar = TRUE) {
  summ <- data %>%
    filter(!is.na(.data[[trait]])) %>%
    { if (include_cultivar)
      group_by(., cycle, cultivar, watering_treatment, aphid_treatment, day)
      else
        group_by(., cycle, watering_treatment, aphid_treatment, day) } %>%
    summarise(mean_val = mean(.data[[trait]], na.rm = TRUE),
              se_val   = sem(.data[[trait]]),
              .groups  = "drop")
  
  if (include_cultivar) {
    summ <- summ %>%
      mutate(grp = interaction(cultivar, watering_treatment, aphid_treatment))
  } else {
    summ <- summ %>%
      mutate(grp = interaction(watering_treatment, aphid_treatment))
  }
  
  p <- ggplot(summ, aes(x = day, y = mean_val,
                        colour = watering_treatment,
                        shape  = aphid_treatment,
                        group  = grp)) +
    geom_line(aes(linetype = if (include_cultivar) cultivar else NULL),
              linewidth = 0.4, position = pd) +
    geom_point(size = 2.2, position = pd, fill = "white", stroke = 1.1) +
    geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                  linetype = "solid",
                  width = 0.6, position = pd, alpha = 0.45) +
    scale_colour_manual(
      values = trt_cols,
      name   = "Watering treatment",
      labels = c("drought" = "Droughted", "water" = "Watered")
    ) +
    scale_shape_manual(
      values = c("no aphids" = 1, "aphids" = 17),
      name   = "Aphid treatment",
      labels = c("no aphids" = "No Aphids", "aphids" = "Aphid Infested")
    ) +
    scale_y_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.10))) +
    labs(x = "Days from Aphid Inoculation", y = ylab) +
    theme_classic(base_size = 13) +
    theme(
      panel.grid       = element_blank(),
      axis.line        = element_line(colour = "black"),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold"),
      strip.placement  = "outside"
    )
  
  if (include_cultivar) {
    p <- p + scale_linetype_manual(
      values = c("bere" = "solid", "laureate" = "dashed"),
      name   = "Cultivar",
      labels = c("bere" = "Bere", "laureate" = "Laureate")
    )
  }
  
  if (show_facet) {
    p <- p + facet_wrap(
      ~ cycle,
      strip.position = "bottom",
      labeller = as_labeller(function(x) paste("Cycle", x))
    )
  }
  p
}

save_traj_plots <- function(trait, ylab, prefix = "traj_all",
                            include_cultivar = TRUE, data = long) {
  p_comb <- plot_growth_traj(trait, ylab, data = data,
                             show_facet = TRUE,
                             include_cultivar = include_cultivar)
  
  panel <- p_comb +
    plot_annotation(tag_levels = "A")
  
  ggsave(paste0(prefix, "_", gsub(" |/", "_", trait), ".jpeg"), panel,
         width = 14, height = 10, dpi = 300, units = "in")
}

# All-plant trajectory plots
walk2(
  growth_traits,
  c("Plant Height (cm)", "Total Tillers (n)", "Main Stem Leaf No.",
    "Dev. Stage (Zadoks)", "SPAD", "Dry Biomass (g)", "Gravimetric Water Content (%)"),
  ~ save_traj_plots(.x, .y, prefix = "traj_all", data = long)
)

