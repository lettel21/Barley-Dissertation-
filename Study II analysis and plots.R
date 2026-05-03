# GLASSHOUSE EXPERIMENT — FULL ANALYSIS SCRIPT


# =============================================================================
# SECTION 1: PACKAGES
# =============================================================================

install.packages(c("tidyverse", "janitor", "lme4", "lmerTest",
                   "car", "patchwork", "ggfortify"))

library(tidyverse)
library(janitor)
library(lme4)
library(lmerTest)   # adds p-values to lme4 output
library(car)        # Levene's test
library(patchwork)  # multi-panel figures
library(ggfortify)  # ggplot-style PCA biplots


# =============================================================================
# SECTION 2: LOAD & CLEAN DATA
# =============================================================================

# Read CSV — treat "na" strings and blanks as NA
GH_data <- read_csv("Glasshouse experiment .csv",
                    na = c("", "NA", "na"))

# Standardise column names (snake_case, removes trailing spaces)
GH_data <- clean_names(GH_data)

# Verify column names look correct
names(GH_data)

# Remove entirely blank rows
# (those rows have NA in the cultivar column, which is column 1)
GH_data <- GH_data %>%
  filter(!is.na(cultivar))

# Set factor types
GH_data <- GH_data %>%
  mutate(
    cultivar  = factor(cultivar,  levels = c("bere", "laureate")),
    treatment = factor(treatment, levels = c("drought", "water")),
    plant_id  = as.factor(plant_id)
  )

# Create a combined group factor for plotting
GH_data <- GH_data %>%
  mutate(group = factor(
    interaction(cultivar, treatment),
    levels = c("bere.drought", "bere.water", "laureate.drought", "laureate.water")
  ))

# Quick overview
glimpse(GH_data)
summary(GH_data)


# =============================================================================
# SECTION 3: PIVOT CONTINUOUS DATA TO LONG FORMAT
# =============================================================================

# Columns matching "_day" are the repeated-measures variables
# Pattern: {measure}_day{day_number} → e.g. mv_day8, gwc_day15, spad_day21

GH_long <- GH_data %>%
  pivot_longer(
    cols      = matches("_day"),          # selects all mv/gwc/spad day columns
    names_to  = c("measure", "day"),
    names_sep = "_day",                   # splits at "_day"
    values_to = "value"
  ) %>%
  mutate(day = as.numeric(day)) %>%
  pivot_wider(
    names_from  = measure,
    values_from = value
  )

# Resulting columns: cultivar, treatment, plant_id, group, day, mv, gwc, spad
# (plus harvest columns — these repeat per row but aren't used in this section)

glimpse(GH_long)

# =============================================================================
# SECTION 4: HARVEST DATA — ASSUMPTION CHECKS
# =============================================================================

# Harvest variables of interest
harvest_vars <- c("root_mass", "shoot_mass", "ratio",
                  "ms_height",
                  "primary_tillers", "secondary_tillers", "total_tillers",
                  "ms_nodes", "dev_stage")
# not currently included, "stems_with_nodes", "stems_w_1", "stems_w_2", "stems_w_3","stems_w_4","stems_w_5")

# Note on dev_stage: this is the Zadoks growth stage — an ordinal scale.
# Most plants are at stage 33, so variance is very low. Treat ANOVA results
# with caution; it is included here for completeness.

# Note on ms_nodes (and Stems_w data): count data. Normality unlikely to be perfectly met
# with small n, but ANOVA is reasonably robust to mild departures.

# Isolate harvest data — one row per plant, ensure numerics
GH_harvest <- GH_data %>%
  select(cultivar, treatment, plant_id, group,
         all_of(harvest_vars)) %>%
  mutate(across(all_of(harvest_vars), as.numeric))

# ---- 4a: Shapiro-Wilk normality test (per group, per variable) ----
# p > 0.05 = normality assumption met for that group

shapiro_results <- GH_harvest %>%
  group_by(cultivar, treatment) %>%
  summarise(
    across(
      all_of(harvest_vars),
      ~ if (sum(!is.na(.)) >= 3) shapiro.test(na.omit(.))$p.value else NA_real_,
      .names = "sw_{.col}"
    ),
    .groups = "drop"
  )

print(shapiro_results)

# ---- 4b: QQ plots (overall, one per variable) ----
par(mfrow = c(3, 3))
for (v in harvest_vars) {
  qqnorm(GH_harvest[[v]], main = paste("QQ —", v), na.action = na.exclude)
  qqline(GH_harvest[[v]], na.action = na.exclude)
}
par(mfrow = c(1, 1))

# ---- 4c: Levene's test — homogeneity of variance across groups ----
# p > 0.05 = variances are homogeneous (assumption met)

for (v in harvest_vars) {
  cat("\nLevene's test —", v, "\n")
  formula <- as.formula(paste(v, "~ cultivar * treatment"))
  print(leveneTest(formula, data = GH_harvest))
}


# =============================================================================
# SECTION 5: HARVEST DATA — TWO-WAY ANOVAs
# =============================================================================

# Helper function: runs ANOVA, prints diagnostic plots, and runs Tukey post-hoc
run_anova <- function(var_name, data) {
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("TWO-WAY ANOVA:", var_name, "\n")
  cat(strrep("=", 60), "\n", sep = "")
  
  formula <- as.formula(paste(var_name, "~ cultivar * treatment"))
  model   <- aov(formula, data = data)
  
  print(summary(model))
  
  # Diagnostic plots (residuals, QQ, scale-location, leverage)
  par(mfrow = c(2, 2))
  plot(model, main = var_name)
  par(mfrow = c(1, 1))
  
  # Tukey HSD post-hoc — useful even when not significant to see direction
  cat("\nTukey HSD:\n")
  print(TukeyHSD(model))
}

for (v in harvest_vars) {
  run_anova(v, GH_harvest)
}

library(tidyverse)
library(patchwork)

# =============================================================================
# SECTION 6: HARVEST DATA — BAR CHARTS WITH SEM ERROR BARS
# =============================================================================

# ── Shared aesthetics ──────────────────────────────────────────────────────────
treat_colours <- c("drought" = "tomato", "water" = "steelblue")
line_colours  <- c("drought" = "tomato", "water" = "steelblue")

harvest_vars <- c("root_mass", "shoot_mass", "ratio",
                  "ms_height", "primary_tillers", "secondary_tillers",
                  "total_tillers", "ms_nodes", "dev_stage")

# ── Summary stats for bar charts ───────────────────────────────────────────────
group_stats <- GH_harvest %>%
  group_by(cultivar, treatment) %>%
  summarise(
    across(all_of(harvest_vars),
           list(mean = ~ mean(., na.rm = TRUE),
                sem  = ~ sd(., na.rm = TRUE) / sqrt(sum(!is.na(.)))),
           .names = "{.fn}_{.col}"),
    .groups = "drop"
  ) %>%
  mutate(cultivar = recode(cultivar, "bere" = "Bere", "laureate" = "Laureate"))

# ── Bar chart function ─────────────────────────────────────────────────────────
make_bar <- function(mean_col, sem_col, y_label) {
  pd <- position_dodge(width = 0.7)
  ggplot(group_stats,
         aes(x = cultivar, y = .data[[mean_col]], fill = treatment)) +
    geom_bar(stat = "identity", colour = "black", width = 0.6, position = pd) +
    geom_errorbar(aes(ymin = .data[[mean_col]] - .data[[sem_col]],
                      ymax = .data[[mean_col]] + .data[[sem_col]]),
                  width = 0.25, linewidth = 0.6, position = pd) +
    scale_fill_manual(values = treat_colours,
                      name   = "Watering treatment",
                      labels = c("drought" = "Droughted", "water" = "Watered")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = y_label) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 9))
}

p_root   <- make_bar("mean_root_mass",         "sem_root_mass",         "Dry Root Mass (g)")
p_shoot  <- make_bar("mean_shoot_mass",        "sem_shoot_mass",        "Dry Shoot Mass (g)")
p_ratio  <- make_bar("mean_ratio",             "sem_ratio",             "Root:Shoot Ratio")
p_height <- make_bar("mean_ms_height",         "sem_ms_height",         "Plant Height (cm)")
p_ptill  <- make_bar("mean_primary_tillers",   "sem_primary_tillers",   "No. Primary Tillers")
p_stilt  <- make_bar("mean_secondary_tillers", "sem_secondary_tillers", "No. Secondary Tillers")

harvest_combined <- (p_root + p_shoot + p_ratio +
                       p_height + p_ptill + p_stilt) +
  plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "right",
        plot.tag = element_text(face = "bold", size = 12))

harvest_combined
ggsave("harvest_bar_charts.jpeg", harvest_combined,
       width = 12, height = 8, dpi = 300)

# =============================================================================
# SECTION 7: CONTINUOUS DATA — GLMMs
# =============================================================================
# Model structure: outcome ~ cultivar * treatment * day + (1 | plant_id)
# Random intercept for plant_id accounts for repeated measures on the same plant.
# If random effect variance ≈ 0, simplify to lm() — see notes below.
#
# Assumptions are checked AFTER fitting the model (residuals require a fitted model).
# =============================================================================

# ---- 7a: SPAD model ----

model_spad <- lmer(spad ~ cultivar * treatment * day + (1 | plant_id),
                   data = GH_long)
summary(model_spad)
anova(model_spad)

# Assumption checks
# 1. Normality of residuals
qqnorm(resid(model_spad), main = "QQ Plot — SPAD Residuals")
qqline(resid(model_spad))

# 2. Homogeneity of variance (residuals vs fitted)
plot(fitted(model_spad), resid(model_spad),
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted — SPAD")
abline(h = 0, lty = 2)

# 3. Flag large residuals
cat("Large SPAD residuals (|resid| > 8):\n")
print(which(abs(resid(model_spad)) > 8))

# 4. Normality of random effects
qqnorm(ranef(model_spad)$plant_id[, 1], main = "QQ Plot — SPAD Random Effects")
qqline(ranef(model_spad)$plant_id[, 1])


# ---- 7b: GWC model (calibrated soil moisture — preferred over mV) ----

model_gwc <- lmer(gwc ~ cultivar * treatment * day + (1 | plant_id),
                  data = GH_long)
summary(model_gwc)
anova(model_gwc)

# If the random effect variance for plant_id is 0 (or near 0), the model
# simplifies — uncomment the line below to use a standard linear model instead:
# model_gwc_simple <- lm(gwc ~ cultivar * treatment * day, data = GH_long)

# Assumption checks
qqnorm(resid(model_gwc), main = "QQ Plot — GWC Residuals")
qqline(resid(model_gwc))

plot(fitted(model_gwc), resid(model_gwc),
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted — GWC")
abline(h = 0, lty = 2)

cat("Large GWC residuals (|resid| > 8):\n")
print(which(abs(resid(model_gwc)) > 8))

qqnorm(ranef(model_gwc)$plant_id[, 1], main = "QQ Plot — GWC Random Effects")
qqline(ranef(model_gwc)$plant_id[, 1])

#simpler model, random effects 0 <- below correct?
model_gwc_simple <- lm(gwc ~ cultivar * treatment * day, data = GH_long)
summary(model_gwc_simple)
anova(model_gwc_simple)


# ---- 7c: mV model (raw sensor values — retained for reference only) ----
# GWC is the calibrated and preferred measure; mV is included for comparison.

model_mv <- lmer(m_v ~ cultivar * treatment * day + (1 | plant_id),
                 data = GH_long)
summary(model_mv)
anova(model_mv)

# If random effect variance ≈ 0, use simpler model:
# model_mv_simple <- lm(mv ~ cultivar * treatment * day, data = GH_long)


# =============================================================================
# SECTION 8: CONTINUOUS DATA — LINE GRAPHS
# =============================================================================

# ── Summary stats for line graphs ─────────────────────────────────────────────
sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

calc_line_stats <- function(data, var) {
  data %>%
    group_by(cultivar, treatment, day) %>%
    summarise(mean_val = mean(.data[[var]], na.rm = TRUE),
              sem_val  = sem(.data[[var]]),
              .groups  = "drop")
}

spad_stats <- calc_line_stats(GH_long, "spad")
gwc_stats  <- calc_line_stats(GH_long, "gwc")

# ── Line graph function ────────────────────────────────────────────────────────
make_line <- function(stats_df, y_label) {
  pd <- position_dodge(width = 0.3)
  ggplot(stats_df,
         aes(x = day, y = mean_val,
             colour   = treatment,
             linetype = cultivar,
             group    = interaction(cultivar, treatment))) +
    geom_line(linewidth = 0.9, position = pd) +
    geom_point(size = 2.2, position = pd) +
    geom_errorbar(aes(ymin = mean_val - sem_val,
                      ymax = mean_val + sem_val),
                  width = 0.6, alpha = 0.4, position = pd) +
    scale_colour_manual(values = line_colours,
                        name   = "Watering treatment",
                        labels = c("drought" = "Droughted", "water" = "Watered")) +
    scale_linetype_manual(values = c("bere" = "solid", "laureate" = "dashed"),
                          name   = "Cultivar",
                          labels = c("bere" = "Bere", "laureate" = "Laureate")) +
    scale_x_continuous(breaks = c(8, 15, 21, 28),
                       labels = c("Day 8", "Day 15", "Day 21", "Day 28")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Days From Treatment Commenced", y = y_label) +
    theme_classic(base_size = 11) +
    theme(legend.position = "right",
          axis.line = element_line(colour = "black"))
}

p_spad <- make_line(spad_stats, "SPAD")
p_gwc  <- make_line(gwc_stats,  "Gravimetric Water Content (%)")

p_spad
p_gwc

ggsave("SPAD_over_time.jpeg", p_spad, width = 8, height = 5, dpi = 300)
ggsave("gwc_over_time.jpeg",  p_gwc,  width = 8, height = 5, dpi = 300)

# =============================================================================
# SECTION 9: PCA
# =============================================================================
#

# Variables to include in PCA
pca_vars <- c("root_mass", "shoot_mass", "ratio",
              "ms_height",
              "primary_tillers", "secondary_tillers")

# Prepare input — complete cases only (plant 'a' removed by drop_na)
GH_pca_input <- GH_harvest %>%
  select(cultivar, treatment, plant_id, all_of(pca_vars)) %>%
  mutate(across(all_of(pca_vars), as.numeric)) %>%
  drop_na()

cat("Plants retained in PCA:", nrow(GH_pca_input), "\n")
cat("Plants removed (missing data):",
    setdiff(levels(GH_harvest$plant_id), GH_pca_input$plant_id), "\n")

# Run PCA
pca_result <- prcomp(GH_pca_input[, pca_vars], scale = TRUE, center = TRUE)

# Summary — variance explained per component
summary(pca_result)

# Variable loadings — which traits drive each PC
print(pca_result$rotation)

# ---- Scree plot ----
scree_df <- data.frame(
  PC       = paste0("PC", seq_along(pca_result$sdev)),
  Variance = (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
)
scree_df$PC <- factor(scree_df$PC, levels = scree_df$PC)

p_scree <- ggplot(scree_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "grey70", colour = "black") +
  geom_line(aes(group = 1), colour = "tomato", size = 1) +
  geom_point(colour = "tomato", size = 2) +
  labs(
    x     = "Principal Component",
    y     = "Variance Explained (%)",
    title = "Scree Plot — PCA of Harvest Traits"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

p_scree
ggsave("PCA_scree.jpeg", p_scree, width = 6, height = 4, dpi = 300)

# ---- Biplot — PC1 vs PC2 ----
p_biplot <- autoplot(
  pca_result,
  data                  = GH_pca_input,
  colour                = "treatment",
  shape                 = "cultivar",
  loadings              = TRUE,
  loadings.label        = TRUE,
  loadings.label.size   = 3,
  loadings.colour       = "black",
  loadings.label.colour = "black",
  size                  = 3
) +
  scale_colour_manual(
    values = c("drought" = "tomato", "water" = "steelblue"),
    labels = c("drought" = "Drought", "water" = "Water")
  ) +
  scale_shape_manual(
    values = c("bere" = 16, "laureate" = 17),
    labels = c("bere" = "Bere", "laureate" = "Laureate")
  ) +
  labs(
    colour = "Treatment",
    shape  = "Cultivar",
    title  = "PCA Biplot — Harvest Traits\n(PC1 vs PC2)"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

p_biplot
ggsave("PCA_biplot.jpeg", p_biplot, width = 7, height = 6, dpi = 300)
