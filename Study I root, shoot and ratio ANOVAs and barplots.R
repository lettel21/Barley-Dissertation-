#root, shoot and ratio at harvest ANOVAs and bar plots 


# ---- 1. PERCENTAGE DIFFERENCE TABLE (APHID IMPACT) ------------------------
# Addresses: "was the % reduction in biomass greater for laureate than bere?"

cat("\n=== PERCENTAGE REDUCTION IN BIOMASS DUE TO APHIDS ===\n")

pct_table <- harvest_full %>%
  group_by(cultivar, aphid_treatment) %>%
  summarise(
    n              = sum(!is.na(dry_biomass)),
    mean_root      = mean(root,       na.rm = TRUE),
    mean_shoot     = mean(shoot,      na.rm = TRUE),
    mean_dry_bio   = mean(dry_biomass,na.rm = TRUE),
    mean_ratio     = mean(ratio,      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = aphid_treatment,
    values_from = c(n, mean_root, mean_shoot, mean_dry_bio, mean_ratio)
  ) %>%
  mutate(
    pct_red_root    = pct_reduction(`mean_root_no aphids`,    `mean_root_aphids`),
    pct_red_shoot   = pct_reduction(`mean_shoot_no aphids`,   `mean_shoot_aphids`),
    pct_red_biomass = pct_reduction(`mean_dry_bio_no aphids`, `mean_dry_bio_aphids`),
    d_biomass       = mapply(function(cv) {
      no_aph <- harvest_full$dry_biomass[harvest_full$cultivar == cv &
                                           harvest_full$aphid_treatment == "no aphids"]
      aph    <- harvest_full$dry_biomass[harvest_full$cultivar == cv &
                                           harvest_full$aphid_treatment == "aphids"]
      cohens_d(no_aph, aph)
    }, cultivar)
  )

cat("\nBiomass: mean (no aphids) | mean (aphids) | % reduction | Cohen's d\n")
pct_table %>%
  select(cultivar,
         `mean_dry_bio_no aphids`, `mean_dry_bio_aphids`,
         pct_red_biomass, d_biomass) %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  print()

cat("\nShoot mass % reduction by cultivar:\n")
pct_table %>% select(cultivar, `mean_shoot_no aphids`, `mean_shoot_aphids`, pct_red_shoot) %>%
  mutate(across(where(is.numeric), ~ round(., 2))) %>% print()

cat("\nRoot mass % reduction by cultivar:\n")
pct_table %>% select(cultivar, `mean_root_no aphids`, `mean_root_aphids`, pct_red_root) %>%
  mutate(across(where(is.numeric), ~ round(., 2))) %>% print()

# Also by watering treatment
cat("\n--- Biomass % reduction by cultivar × watering ---\n")
harvest_full %>%
  group_by(cultivar, watering_treatment, aphid_treatment) %>%
  summarise(mean_bio = mean(dry_biomass, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = aphid_treatment, values_from = mean_bio) %>%
  mutate(pct_red = pct_reduction(`no aphids`, aphids)) %>%
  mutate(across(where(is.numeric), ~ round(., 2))) %>%
  print()


# ---- 2. HARVEST ANOVAs — ALL PLANTS ----------------------------------------

all_harvest_traits <- c("root", "shoot", "ratio")

cat("\n=== HARVEST ANOVAs — ALL PLANTS ===\n")

harvest_models <- map(all_harvest_traits, function(trait) {
  df  <- harvest_full %>%
    filter(!is.na(.data[[trait]])) %>%
    droplevels()
  frm <- as.formula(
    paste(trait, "~ cultivar * aphid_treatment * watering_treatment + cycle")
  )
  mod <- aov(frm, data = df)
  eta <- partial_eta2(mod)
  
  cat("\n---", trait, "---\n")
  print(summary(mod))
  
  cat("Partial eta-squared (η²p):\n")
  print(eta)
  
  sw <- shapiro.test(residuals(mod))
  cat("Shapiro-Wilk: W =", round(sw$statistic, 3),
      ", p =", round(sw$p.value, 3), "\n")
  
  list(trait = trait, model = mod, eta2 = eta, shapiro = sw)
})
names(harvest_models) <- all_harvest_traits

# ---- 3. HARVEST ANOVAs - CYCLE 2 ONLY --------------------------------------
# Cycle term dropped; aphid × watering × cultivar model

cat("\n=== HARVEST ANOVAs — CYCLE 2 ONLY ===\n")

harvest_models_c2 <- map(all_harvest_traits, function(trait) {
  df  <- harvest_2 %>%
    filter(!is.na(.data[[trait]])) %>%
    droplevels()
  frm <- as.formula(
    paste(trait, "~ aphid_treatment * watering_treatment * cultivar")
  )
  mod <- aov(frm, data = df)
  eta <- partial_eta2(mod)
  
  cat("\n---", trait, "(cycle 2 only) ---\n")
  print(summary(mod))
  cat("Partial eta-squared (η²p):\n")
  print(eta)
  
  sw <- shapiro.test(residuals(mod))
  cat("Shapiro-Wilk: W =", round(sw$statistic, 3),
      ", p =", round(sw$p.value, 3), "\n")
  
  list(trait = trait, model = mod, eta2 = eta, shapiro = sw)
})
names(harvest_models_c2) <- all_harvest_traits


# ---- 4. HARVEST ANOVAs — LAUREATE ONLY -------------------------------------
# Cultivar term dropped; aphid × watering × cycle model

cat("\n=== HARVEST ANOVAs — LAUREATE ONLY ===\n")

harvest_models_laur <- map(all_harvest_traits, function(trait) {
  df  <- harvest_laur %>%
    filter(!is.na(.data[[trait]])) %>%
    droplevels()
  frm <- as.formula(
    paste(trait, "~ aphid_treatment * watering_treatment + cycle")
  )
  mod <- aov(frm, data = df)
  eta <- partial_eta2(mod)
  
  cat("\n---", trait, "(laureate only) ---\n")
  print(summary(mod))
  cat("Partial eta-squared (η²p):\n")
  print(eta)
  
  sw <- shapiro.test(residuals(mod))
  cat("Shapiro-Wilk: W =", round(sw$statistic, 3),
      ", p =", round(sw$p.value, 3), "\n")
  
  list(trait = trait, model = mod, eta2 = eta, shapiro = sw)
})
names(harvest_models_laur) <- all_harvest_traits


# ---- 5. SENSITIVITY: HARVEST ANOVAs — EXCLUDING FLAGGED DEAD PLANTS --------
# Four cycle 1 Bere aphid plants flagged any_dead = TRUE on day 23:
#   1_Left_F, 1_Right_F, 1_Left_H, 1_Right_H
# These are compared to the full model (section 6) as a sensitivity check.

dead_ids <- c("1_Left_F", "1_Right_F", "1_Left_H", "1_Right_H")

harvest_nodead   <- harvest_full %>% filter(!unique_id %in% dead_ids)
long_nodead      <- long         %>% filter(!unique_id %in% dead_ids)

cat("\n=== SENSITIVITY: HARVEST ANOVAs — DEAD PLANTS EXCLUDED ===\n")
cat("Excluded unique_ids:", paste(dead_ids, collapse = ", "), "\n")
cat("n remaining (harvest):", nrow(harvest_nodead), "\n\n")

harvest_models_nodead <- map(all_harvest_traits, function(trait) {
  df  <- harvest_nodead %>%
    filter(!is.na(.data[[trait]])) %>%
    droplevels()
  frm <- as.formula(
    paste(trait, "~ cultivar * aphid_treatment * watering_treatment + cycle")
  )
  mod <- aov(frm, data = df)
  eta <- partial_eta2(mod)
  
  cat("\n---", trait, "(dead plants excluded) ---\n")
  print(summary(mod))
  cat("Partial eta-squared (η²p):\n")
  print(eta)
  
  sw <- shapiro.test(residuals(mod))
  cat("Shapiro-Wilk: W =", round(sw$statistic, 3),
      ", p =", round(sw$p.value, 3), "\n")
  
  list(trait = trait, model = mod, eta2 = eta, shapiro = sw)
})
names(harvest_models_nodead) <- all_harvest_traits

# ---- 8. HARVEST BAR PLOTS --------------------------------------------------
# Three versions per trait: combined (both cycles), cycle 1, cycle 2
# Laureate-only versions also saved

plot_harvest_bar <- function(trait, ylab, data = harvest_full) {
  data %>%
    group_by(cultivar, watering_treatment, aphid_treatment) %>%
    summarise(mean_val = mean(.data[[trait]], na.rm = TRUE),
              se_val   = sem(.data[[trait]]),
              .groups  = "drop") %>%
    ggplot(aes(x = cultivar, y = mean_val, fill = watering_treatment)) +
    geom_col(position = position_dodge(0.75), width = 0.7, colour = "black") +
    geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                  position = position_dodge(0.75), width = 0.25, alpha = 0.6) +
    scale_fill_manual(values = trt_cols, name = "Watering treatment",
                      labels = c("drought" = "Droughted", "water" = "Watered")) +
    scale_x_discrete(labels = c("bere" = "Bere", "laureate" = "Laureate")) +
    scale_y_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.10))) +
    facet_wrap(~ aphid_treatment,
               strip.position = "bottom",
               labeller = labeller(aphid_treatment = c(
                 "aphids" = "Aphid Infested",
                 "no aphids" = "No Aphids"
               ))) +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = 11) +
    theme(
      strip.placement  = "outside",
      panel.grid       = element_blank(),
      axis.line        = element_line(colour = "black"),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold")
    )
}


# Function to save three versions (combined, c1, c2) for a bar plot
save_harvest_plots <- function(trait, ylab, prefix) {
  p_all <- plot_harvest_bar(trait, ylab, data = harvest_full) 
  p_c1  <- plot_harvest_bar(trait, ylab, data = harvest_full %>% filter(cycle == 1)) 
  p_c2  <- plot_harvest_bar(trait, ylab, data = harvest_full %>% filter(cycle == 2))
  panel <- (p_all | p_c1 | p_c2) +
    plot_annotation(tag_levels = "A")
  ggsave(paste0(prefix, "_", trait, ".jpeg"), panel,
         width = 16, height = 5, dpi = 300, units = "in")
}

# All plants harvest bar panels
save_harvest_plots("root",           "Dry Root Mass (g)",                        "harvest_bars_all")
save_harvest_plots("shoot",          "Dry Shoot Mass (g)",                       "harvest_bars_all")
save_harvest_plots("ratio",          "Root:Shoot Ratio",                     "harvest_bars_all")



