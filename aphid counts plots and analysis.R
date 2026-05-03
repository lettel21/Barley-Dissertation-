#aphid counts plots and analysis 

# ---- 1. APHID COUNTS -------------------------------------------------------

cat("\n=== APHID COUNT GLMM — ALL PLANTS (aphid treatment only) ===\n")

aphid_df <- long %>%
  filter(aphid_treatment == "aphids", !is.na(log_aphid_count)) %>%
  mutate(day_c = as.numeric(scale(day, center = TRUE, scale = FALSE)))

aphid_mod <- lmer(
  log_aphid_count ~ cultivar * watering_treatment * day_c +
    cycle + (1 | unique_id),
  data = aphid_df, REML = FALSE
)
cat("--- Aphid count model ---\n")
print(anova(aphid_mod, type = "III"))

cat("\n=== APHID COUNT GLMM — LAUREATE ONLY ===\n")
aphid_df_laur <- long_laur %>%
  filter(aphid_treatment == "aphids", !is.na(log_aphid_count)) %>%
  mutate(day_c = as.numeric(scale(day, center = TRUE, scale = FALSE)))

cat("\n=== APHID COUNT GLMM — cycle 2 ONLY ===\n")
aphid_df_2 <- long_2 %>%
  filter(aphid_treatment == "aphids", !is.na(log_aphid_count)) %>%
  mutate(day_c = as.numeric(scale(day, center = TRUE, scale = FALSE)))

aphid_mod_2 <- lmer(
  log_aphid_count ~ cultivar * watering_treatment * day_c + (1 | unique_id),
  data = aphid_df_2, REML = FALSE
)

aphid_mod_laur <- lmer(
  log_aphid_count ~ watering_treatment * day_c +
    cycle + (1 | unique_id),
  data = aphid_df_laur, REML = FALSE
)
cat("--- Aphid count model (laureate only) ---\n")
print(anova(aphid_mod_laur, type = "III"))
cat("--- Aphid count model (cycle 2 only) ---\n")
print(anova(aphid_mod_2, type = "III"))

# Residual checks (all plants model)
par(mfrow = c(1, 2))
qqnorm(residuals(aphid_mod)); qqline(residuals(aphid_mod), col = "red")
plot(fitted(aphid_mod), residuals(aphid_mod),
     xlab = "Fitted", ylab = "Residuals", main = "Aphid residuals")
par(mfrow = c(1, 1))

# Aphid plots — three versions
save_aphid_plots <- function(data, prefix, include_cultivar = TRUE) {
  summ_fn <- function(d) {
    if (include_cultivar) {
      d %>%
        group_by(cycle, cultivar, watering_treatment, day) %>%
        summarise(mean_log = mean(log_aphid_count, na.rm = TRUE),
                  se_log   = sem(log_aphid_count), .groups = "drop") %>%
        mutate(grp = interaction(cultivar, watering_treatment))
    } else {
      d %>%
        group_by(cycle, watering_treatment, day) %>%
        summarise(mean_log = mean(log_aphid_count, na.rm = TRUE),
                  se_log   = sem(log_aphid_count), .groups = "drop") %>%
        mutate(grp = watering_treatment)
    }
  }
  
  mk <- function(d, title_sfx) {
    s <- summ_fn(d %>% filter(!is.na(log_aphid_count)))
    p <- ggplot(s, aes(x = day, y = mean_log,
                       colour   = watering_treatment,
                       group    = grp)) +
      geom_line(linewidth = 0.4, position = pd) +
      geom_point(size = 2.2, position = pd, shape = 16, stroke = 1.1) +
      geom_errorbar(aes(ymin = mean_log - se_log, ymax = mean_log + se_log),
                    width = 0.6, position = pd, alpha = 0.45) +
      scale_colour_manual(values = trt_cols, name = "Watering") +
      scale_y_continuous(limits = c(0, NA),
                         expand = expansion(mult = c(0, 0.10))) +
      labs(x = "Day in insectary", y = "Log aphid count", title = title_sfx) +
      theme_classic(base_size = 11) +
      theme(panel.grid    = element_blank(),
            axis.line     = element_line(colour = "black"),
            strip.background = element_blank(),
            strip.text    = element_text(face = "bold"))
    if (include_cultivar) {
      p <- p + aes(linetype = cultivar) +
        scale_linetype_manual(values = c("bere" = "solid", "laureate" = "dashed"),
                              name = "Cultivar")
    }
    p
  }
  
  p_comb <- mk(data, "Both cycles") + facet_wrap(~ cycle, labeller = label_both)
  p_c1   <- mk(data %>% filter(cycle == 1), "Cycle 1")
  p_c2   <- mk(data %>% filter(cycle == 2), "Cycle 2")
  
  panel <- (p_comb / (p_c1 | p_c2)) +
    plot_annotation(tag_levels = "A",
                    caption = "Aphid-treatment plants only; points = mean ± SEM")
  ggsave(paste0(prefix, "_aphid_counts.jpeg"), panel,
         width = 12, height = 10, dpi = 300, units = "in")
}

save_aphid_plots(aphid_df,      "aphids_all",  include_cultivar = TRUE)
save_aphid_plots(aphid_df_laur, "aphids_laur", include_cultivar = FALSE)
