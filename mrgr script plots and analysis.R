#mrgr plots and analysis 

# ---- 0. PACKAGES -----------------------------------------------------------
library(tidyverse)
library(janitor)
library(lme4)
library(lmerTest)   # Satterthwaite / KR df for lmer
library(car)        # Anova() Type II; leveneTest()
library(emmeans)    # estimated marginal means, pairwise contrasts
library(patchwork)


# ---- 1. LOAD & CLEAN -------------------------------------------------------

raw <- read_csv("mrgr3.csv", show_col_types = FALSE)

mrgr <- raw %>%
  # Drop Excel's trailing empty columns
  select(any_of(c("plant_id", "cycle", "cultivar",
                  "treatment", "treatment ",
                  "old_mass", "old_mass ",
                  "young_mass", "ln_old", "ln_young", "RGR"))) %>%
  rename_with(~ str_trim(.)) %>%
  mutate(across(where(is.character), str_trim)) %>%
  filter(!is.na(plant_id)) %>%
  mutate(
    cycle     = factor(cycle, levels = c("cycle 1", "cycle 2"),
                       labels = c("1", "2")),
    cultivar  = factor(cultivar, levels = c("Bere", "Laureate")),
    treatment = factor(treatment, levels = c("W", "D"),
                       labels = c("water", "drought")),
    # plant_id letters are reused across cycles — build a unique id
    unique_id = factor(paste(cycle, plant_id, sep = "_"))
  )

cat("=== Aphid-level rows per cell ===\n")
mrgr %>% count(cycle, cultivar, treatment) %>% print()

cat("\n=== Plants per cell ===\n")
mrgr %>%
  distinct(cycle, cultivar, treatment, unique_id) %>%
  count(cycle, cultivar, treatment) %>%
  print()

cat("\n=== Aphids per plant ===\n")
mrgr %>% count(cycle, cultivar, treatment, unique_id) %>% print(n = 30)

cat("\n=== Summary of RGR ===\n")
mrgr %>%
  group_by(cycle, cultivar, treatment) %>%
  summarise(n_aphids = n(),
            n_plants = n_distinct(unique_id),
            mean_RGR = mean(RGR, na.rm = TRUE),
            sd_RGR   = sd(RGR,  na.rm = TRUE),
            .groups  = "drop") %>%
  print()

glimpse(mrgr)


# ---- 2. HELPERS ------------------------------------------------------------

sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

trt_cols <- c("water"   = "steelblue",
              "drought" = "tomato")

# Safe Levene's test — returns NA row if cells have < 2 observations.
safe_levene <- function(formula, data) {
  tryCatch(leveneTest(formula, data = data),
           error = function(e) {
             message("Levene's test not computable: ", e$message); NULL
           })
}


# ---- 3. ASSUMPTION CHECKS (raw RGR) ---------------------------------------
# The paper's protocol runs these on the raw response, before model fitting.
# With mixed models we also check residuals below.

cat("\n\n=== ASSUMPTION CHECKS ON RAW RGR ===\n")

cat("\n-- Shapiro-Wilk (overall) --\n")
print(shapiro.test(mrgr$RGR))

cat("\n-- Shapiro-Wilk within cells (cultivar:treatment:cycle) --\n")
mrgr %>%
  group_by(cycle, cultivar, treatment) %>%
  summarise(n = n(),
            W = if (n() >= 3) shapiro.test(RGR)$statistic else NA_real_,
            p = if (n() >= 3) shapiro.test(RGR)$p.value   else NA_real_,
            .groups = "drop") %>%
  print()

cat("\n-- Levene's test (cultivar:treatment:cycle cells) --\n")
lev_all <- safe_levene(RGR ~ interaction(cultivar, treatment, cycle), data = mrgr)
print(lev_all)


# ---- 4. MODEL A — LAUREATE ONLY -------------------------------------------
# Aphid-level RGR with plant as random intercept; treatment and cycle as
# fixed effects.  Cycle is a fixed effect (only 2 levels — too few to
# estimate a random-effect variance component reliably).

cat("\n\n=== MODEL A: Laureate only ===\n")

laur <- mrgr %>% filter(cultivar == "Laureate") %>% droplevels()

cat("Plants per cell (Laureate only):\n")
laur %>% distinct(cycle, treatment, unique_id) %>%
  count(cycle, treatment) %>% print()

cat("\n-- Levene's test within Laureate (treatment:cycle cells) --\n")
print(safe_levene(RGR ~ interaction(treatment, cycle), data = laur))

# Full (interaction) and reduced (no interaction) models, ML-fit for LRT
mod_A_full_ml <- lmer(RGR ~ treatment * cycle + (1 | unique_id),
                      data = laur, REML = FALSE)
mod_A_red_ml  <- lmer(RGR ~ treatment + cycle + (1 | unique_id),
                      data = laur, REML = FALSE)

cat("\n-- LRT: treatment:cycle interaction --\n")
lrt_A <- anova(mod_A_red_ml, mod_A_full_ml)
print(lrt_A)

p_int_A <- lrt_A$`Pr(>Chisq)`[2]
cat(sprintf("\nInteraction p = %.3f — ", p_int_A))
if (!is.na(p_int_A) && p_int_A < 0.05) {
  cat("SIGNIFICANT: keeping full (interaction) model.\n")
  mod_A <- lmer(RGR ~ treatment * cycle + (1 | unique_id),
                data = laur, REML = TRUE)
  emm_A_spec <- ~ treatment | cycle
} else {
  cat("NON-SIGNIFICANT: adopting reduced (no-interaction) model.\n")
  mod_A <- lmer(RGR ~ treatment + cycle + (1 | unique_id),
                data = laur, REML = TRUE)
  emm_A_spec <- ~ treatment
}

cat("\n-- Fixed-effects ANOVA (Type II, Kenward-Roger) --\n")
print(Anova(mod_A, type = 2, test.statistic = "F"))

cat("\n-- Summary --\n")
print(summary(mod_A))

cat("\n-- Residual diagnostics --\n")
print(shapiro.test(residuals(mod_A)))

cat("\n-- Estimated marginal means --\n")
emm_A <- emmeans(mod_A, emm_A_spec)
print(emm_A)
print(pairs(emm_A, adjust = "tukey"))

# Always also print the cell-mean emmeans from the FULL model so the
# cycle-specific pattern is available descriptively even if the reduced
# model is the one used for inference.
cat("\n-- Descriptive cell means (from full model, for reference) --\n")
emm_A_full <- emmeans(lmer(RGR ~ treatment * cycle + (1 | unique_id),
                           data = laur, REML = TRUE),
                      ~ treatment | cycle)
print(emm_A_full)


# ---- 5. MODEL B — CYCLE 2 ONLY --------------------------------------------
# Bere only appears in cycle 2; the only clean Bere-vs-Laureate comparison
# is within cycle 2.  NOTE: cycle 2 Bere W has only 1 plant.

cat("\n\n=== MODEL B: Cycle 2 only ===\n")

c2 <- mrgr %>% filter(cycle == "2") %>% droplevels()

cat("Plants per cell (cycle 2):\n")
c2 %>% distinct(cultivar, treatment, unique_id) %>%
  count(cultivar, treatment) %>% print()

cat("\n-- Levene's test within cycle 2 (cultivar:treatment cells) --\n")
print(safe_levene(RGR ~ interaction(cultivar, treatment), data = c2))

mod_B_full_ml <- lmer(RGR ~ cultivar * treatment + (1 | unique_id),
                      data = c2, REML = FALSE)
mod_B_red_ml  <- lmer(RGR ~ cultivar + treatment + (1 | unique_id),
                      data = c2, REML = FALSE)

cat("\n-- LRT: cultivar:treatment interaction --\n")
lrt_B <- anova(mod_B_red_ml, mod_B_full_ml)
print(lrt_B)

p_int_B <- lrt_B$`Pr(>Chisq)`[2]
cat(sprintf("\nInteraction p = %.3f — ", p_int_B))
if (!is.na(p_int_B) && p_int_B < 0.05) {
  cat("SIGNIFICANT: keeping full (interaction) model.\n")
  mod_B <- lmer(RGR ~ cultivar * treatment + (1 | unique_id),
                data = c2, REML = TRUE)
  emm_B_spec <- ~ cultivar * treatment
} else {
  cat("NON-SIGNIFICANT: adopting reduced (no-interaction) model.\n")
  mod_B <- lmer(RGR ~ cultivar + treatment + (1 | unique_id),
                data = c2, REML = TRUE)
  emm_B_spec <- ~ cultivar + treatment
}

cat("\n-- Fixed-effects ANOVA (Type II, Kenward-Roger) --\n")
print(Anova(mod_B, type = 2, test.statistic = "F"))

cat("\n-- Summary --\n")
print(summary(mod_B))

cat("\n-- Residual diagnostics --\n")
print(shapiro.test(residuals(mod_B)))

cat("\n-- Estimated marginal means --\n")
emm_B <- emmeans(mod_B, emm_B_spec)
print(emm_B)
print(pairs(emm_B, adjust = "tukey"))


# ---- 6. MODEL C — COMBINED, EXPLORATORY ------------------------------------
# Bere is absent from cycle 1; cultivar and cycle are partially confounded.
# Fit for completeness only.

cat("\n\n=== MODEL C: Combined (exploratory, unbalanced) ===\n")

mod_C_full_ml <- lmer(RGR ~ cultivar * treatment + cycle + (1 | unique_id),
                      data = mrgr, REML = FALSE)
mod_C_red_ml  <- lmer(RGR ~ cultivar + treatment + cycle + (1 | unique_id),
                      data = mrgr, REML = FALSE)

cat("\n-- LRT: cultivar:treatment interaction --\n")
lrt_C <- anova(mod_C_red_ml, mod_C_full_ml)
print(lrt_C)

p_int_C <- lrt_C$`Pr(>Chisq)`[2]
cat(sprintf("\nInteraction p = %.3f — ", p_int_C))
if (!is.na(p_int_C) && p_int_C < 0.05) {
  cat("SIGNIFICANT: keeping full (interaction) model.\n")
  mod_C <- lmer(RGR ~ cultivar * treatment + cycle + (1 | unique_id),
                data = mrgr, REML = TRUE)
  emm_C_spec <- ~ cultivar * treatment
} else {
  cat("NON-SIGNIFICANT: adopting reduced (no-interaction) model.\n")
  mod_C <- lmer(RGR ~ cultivar + treatment + cycle + (1 | unique_id),
                data = mrgr, REML = TRUE)
  emm_C_spec <- ~ cultivar + treatment
}

cat("\n-- Fixed-effects ANOVA (Type II, Kenward-Roger) --\n")
print(Anova(mod_C, type = 2, test.statistic = "F"))

cat("\n-- Residual diagnostics --\n")
print(shapiro.test(residuals(mod_C)))

cat("\n-- Estimated marginal means --\n")
emm_C <- emmeans(mod_C, emm_C_spec)
print(emm_C)
print(pairs(emm_C, adjust = "tukey"))


# ---- 7. PLOT A — Laureate, cycle × treatment ------------------------------
# Requirements (per user):
#   - Laureate only, four bars: cycle 1 water, cycle 1 drought,
#     cycle 2 water, cycle 2 drought
#   - y axis starts at zero
#   - no grid behind the plot
#   - y axis label "MRGR"
#   - no n or plant-count labels above bars
#   - panel labelled "A"

summary_cell <- mrgr %>%
  group_by(cycle, cultivar, treatment) %>%
  summarise(mean_RGR = mean(RGR, na.rm = TRUE),
            se_RGR   = sem(RGR),
            n_aphids = n(),
            n_plants = n_distinct(unique_id),
            .groups  = "drop")

summary_laur <- summary_cell %>% filter(cultivar == "Laureate")

p_A <- ggplot(summary_laur,
              aes(x = cycle, y = mean_RGR, fill = treatment)) +
  geom_col(position = position_dodge(0.8), width = 0.7, colour = "black") +
  geom_errorbar(aes(ymin = mean_RGR - se_RGR,
                    ymax = mean_RGR + se_RGR),
                position = position_dodge(0.8), width = 0.25) +
  scale_fill_manual(values = trt_cols, name = "Watering") +
  scale_x_discrete(labels = c("1" = "Cycle 1", "2" = "Cycle 2")) +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.10))) +
  labs(x = NULL, y = "MRGR", tag = "A") +
  theme_classic(base_size = 12) +
  theme(
    panel.grid       = element_blank(),
    axis.line        = element_line(colour = "black"),
    plot.tag         = element_text(face = "bold", size = 16),
    plot.tag.position = c(0.02, 0.97),
    legend.position  = "right"
  )

ggsave("mrgr3_plotA_laureate.jpeg", p_A,
       width = 7, height = 5, dpi = 300, units = "in")


# ---- 8. OTHER FIGURES (retained for reference) ----------------------------

# 8a. Bar plot — full dataset (facet by cycle)
p_bar_all <- ggplot(summary_cell,
                    aes(x = cultivar, y = mean_RGR, fill = treatment)) +
  geom_col(position = position_dodge(0.8), width = 0.7, colour = "black") +
  geom_errorbar(aes(ymin = mean_RGR - se_RGR, ymax = mean_RGR + se_RGR),
                position = position_dodge(0.8), width = 0.25) +
  facet_wrap(~ cycle, labeller = label_both) +
  scale_fill_manual(values = trt_cols, name = "Watering") +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.10))) +
  labs(x = "Cultivar", y = "MRGR") +
  theme_classic(base_size = 12) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey90", colour = NA))

# 8b. Boxplot — full dataset
p_box_all <- ggplot(mrgr, aes(x = cultivar, y = RGR, fill = treatment)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6,
               outlier.shape = NA, alpha = 0.8) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width  = 0.8),
             size = 1.5, alpha = 0.7) +
  facet_wrap(~ cycle, labeller = label_both) +
  scale_fill_manual(values = trt_cols, name = "Watering") +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Cultivar", y = "MRGR") +
  theme_classic(base_size = 12) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey90", colour = NA))

# 8c. Boxplot — Laureate only
p_box_laur <- ggplot(laur, aes(x = cycle, y = RGR, fill = treatment)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6,
               outlier.shape = NA, alpha = 0.8) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width  = 0.8),
             size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = trt_cols, name = "Watering") +
  scale_x_discrete(labels = c("1" = "Cycle 1", "2" = "Cycle 2")) +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "MRGR") +
  theme_classic(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.line  = element_line(colour = "black"))

ggsave("mrgr3_barplot_all.jpeg", p_bar_all,
       width = 10, height = 5, dpi = 300, units = "in")
ggsave("mrgr3_boxplot_laureate.jpeg", p_box_laur,
       width = 7,  height = 5, dpi = 300, units = "in")
ggsave("mrgr3_boxplot_all.jpeg",      p_box_all,
       width = 10, height = 5, dpi = 300, units = "in")


