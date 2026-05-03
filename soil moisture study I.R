# ---- 1. SOIL MOISTURE -------------------------------------------------------
#all plants----
cat("\n=== SOIL MOISTURE (GWC) — ALL PLANTS ===\n")

moisture_df <- long %>%
  filter(!is.na(gwc)) %>%
  mutate(day_c = as.numeric(scale(day, center = TRUE, scale = FALSE)))

moisture_mod <- lmer(
  gwc ~ cultivar * watering_treatment * day_c + aphid_treatment +
    cycle + (1 | unique_id),
  data = moisture_df, REML = FALSE
)
cat("--- GWC (all plants) ---\n")
print(anova(moisture_mod, type = "III"))

# Laureate only -----

cat("\n=== SOIL MOISTURE (GWC) — LAUREATE ONLY ===\n")
moisture_df_laur <- long_laur %>%
  filter(!is.na(gwc)) %>%
  mutate(day_c = as.numeric(scale(day, center = TRUE, scale = FALSE)))

moisture_mod_laur <- lmer(
  gwc ~ watering_treatment * day_c + aphid_treatment +
    cycle + (1 | unique_id),
  data = moisture_df_laur, REML = FALSE
)
cat("--- GWC (laureate only) ---\n")
print(anova(moisture_mod_laur, type = "III"))

# cycle 2 only -----

cat("\n=== SOIL MOISTURE (GWC) — CYCLE 2 ONLY ===\n")
moisture_df_2 <- long_2 %>%
  filter(!is.na(gwc)) %>%
  mutate(day_c = as.numeric(scale(day, center = TRUE, scale = FALSE)))

moisture_mod_2 <- lmer(
  gwc ~ aphid_treatment * watering_treatment * cultivar * day_c +
   + (1 | unique_id),
  data = moisture_df_2, REML = FALSE
)


cat("--- GWC (cycle 2 only) ---\n")
print(anova(moisture_mod_2, type = "III"))

#plots ------------
# Soil moisture plots — three versions
save_moisture_plots <- function(data, prefix) {
  mk <- function(d, title_sfx) {
    d %>%
      filter(!is.na(gwc)) %>%
      group_by(cycle, watering_treatment, day) %>%
      summarise(mean_val = mean(gwc, na.rm = TRUE),
                se_val   = sem(gwc), .groups = "drop") %>%
      ggplot(aes(x = day, y = mean_val, colour = watering_treatment,
                 group = watering_treatment)) +
      geom_line(linewidth = 0.4, position = pd) +
      geom_point(size = 2, shape = 1, stroke = 1.1, position = pd) +
      geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                    width = 0.6, position = pd, alpha = 0.45) +
      scale_colour_manual(values = trt_cols, name = "Watering") +
      scale_y_continuous(limits = c(0, NA),
                         expand = expansion(mult = c(0, 0.10))) +
      labs(x = "Day in insectary", y = "GWC (%)", title = title_sfx) +
      theme_classic(base_size = 11) +
      theme(panel.grid    = element_blank(),
            axis.line     = element_line(colour = "black"),
            strip.background = element_blank(),
            strip.text    = element_text(face = "bold"))
  }
  
  p_comb <- mk(data, "Both cycles") + facet_wrap(~ cycle, labeller = label_both)
  p_c1   <- mk(data %>% filter(cycle == 1), "Cycle 1")
  p_c2   <- mk(data %>% filter(cycle == 2), "Cycle 2")
  
  panel <- (p_comb / (p_c1 | p_c2)) +
    plot_annotation(tag_levels = "A")
  ggsave(paste0(prefix, "_gwc.jpeg"), panel,
         width = 12, height = 10, dpi = 300, units = "in")
}

save_moisture_plots(long,      "moisture_all")
save_moisture_plots(long_laur, "moisture_laur")

