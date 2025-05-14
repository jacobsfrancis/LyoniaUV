# ---------------------------------------------------------------------
# LOAD LIBRARIES
# ---------------------------------------------------------------------
library(tidyverse)   # Data wrangling and plotting
library(zoo)         # Rolling means
library(lme4)        # Linear mixed models
library(emmeans)     # Estimated marginal means
library(multcomp)    # Group letters for emmeans

# ---------------------------------------------------------------------
# STEP 1: Import and Inspect Raw Spectral Data
# ---------------------------------------------------------------------
uv_data <- read_csv("~/Documents/ScienceProjects/2025/LyoniaUV/merged_uvvis_data.csv")

glimpse(uv_data)  # Preview structure

uv_data %>%
  count(species, `window/petal`, sort = TRUE)  # Sample counts by group

# Quick visual QC of a few random spectra
uv_data %>%
  filter(Filename %in% sample(unique(Filename), 3)) %>%
  ggplot(aes(Wavelength, Value, color = Filename)) +
  geom_line() +
  geom_vline(xintercept = c(280, 400), linetype = "dashed", color = "blue") +
  labs(title = "Raw Spectra with UV Range", x = "Wavelength (nm)", y = "Transmission")
# ---------------------------------------------------------------------
# STEP 2: Interpolate Artifact Regions and Trim Spectra
# ---------------------------------------------------------------------
# Visualize potential artifact regions for manual inspection
uv_data %>%
  filter(species %in% c("fruticosa", "lucida")) %>%
  ggplot(aes(Wavelength, Value, group = Filename)) +
  geom_line(alpha = 0.1) +
  geom_vline(xintercept = seq(400, 700, 10), color = "gray90", linetype = "dotted") +
  labs(title = "Raw Transmission Spectra", y = "Transmission", x = "Wavelength (nm)") +
  theme_minimal() +
  scale_y_continuous(limits = c(-10, 100)) +
  scale_x_continuous(limits = c(437,450))

# Interpolation helper function to clean sensor artifacts
interpolate_artifact <- function(x, low, high) {
  interpolated_chunks <- x %>%
    group_by(Filename) %>%
    group_modify(~ {
      non_artifact <- .x %>% filter(Wavelength < low | Wavelength > high)
      artifact <- .x %>% filter(Wavelength >= low & Wavelength <= high)
      if (nrow(artifact) > 0 && nrow(non_artifact) > 1) {
        interpolated <- approx(
          x = non_artifact$Wavelength,
          y = non_artifact$Value,
          xout = artifact$Wavelength,
          rule = 2
        )
        artifact$Value <- interpolated$y
      }
      bind_rows(non_artifact, artifact) %>% arrange(Wavelength)
    }) %>%
    ungroup()
  return(interpolated_chunks)
}

# Apply interpolation to specific artifact-prone regions also not3
# we are trimming from 280-700 here that is a choice we should make sure we
# justify in the paper
uv_trimmed <- uv_data %>%
  filter(Wavelength >= 280, Wavelength <= 800) %>%
  mutate(`window/petal` = case_when(
    str_to_lower(`window/petal`) == "w" ~ "window",
    str_to_lower(`window/petal`) == "p" ~ "petal",
    TRUE ~ as.character(`window/petal`)
  ))

uv_interpolated <- uv_trimmed %>%
  interpolate_artifact(415, 428) %>%
  interpolate_artifact(438,442) %>%
  interpolate_artifact(444, 454) %>%
  interpolate_artifact(481, 499) %>%
  interpolate_artifact(520, 550)

# Visualize interpolated regions
uv_compare_interp <- uv_trimmed %>%
  rename(Original = Value) %>%
  inner_join(uv_interpolated %>% dplyr::select(Filename, Wavelength, Interpolated = Value),
             by = c("Filename", "Wavelength")) %>%
  mutate(interpolated = !near(Original, Interpolated))


# Sample visual comparison of raw vs interpolated curves
example_files <- sample(unique(uv_compare_interp$Filename), 4)
uv_compare_interp %>%
  filter(Filename %in% example_files) %>%
  pivot_longer(cols = c(Original, Interpolated), names_to = "Type", values_to = "Value") %>%
  ggplot(aes(x = Wavelength, y = Value, color = Type)) +
  geom_line(size = 0.8) +
  geom_point(data = . %>% filter(interpolated & Type == "Interpolated"),
             aes(x = Wavelength, y = Value), color = "red", size = 1.5) +
  facet_wrap(~ Filename, scales = "free_y") +
  labs(title = "Comparison of Original vs Interpolated Spectra (Random Sample)",
       y = "Transmission", x = "Wavelength (nm)", color = "Curve Type") +
  theme_minimal()

# ---------------------------------------------------------------------
# STEP 3: Smooth Each Sample with Rolling Average
# ---------------------------------------------------------------------
uv_smoothed <- uv_interpolated %>%
  arrange(Filename, Wavelength) %>%
  group_by(Filename) %>%
  mutate(smoothed_value = rollmean(Value, k = 100, fill = NA, align = "center")) %>%
  ungroup()

# ---------------------------------------------------------------------
# STEP 4: Filter Out Poor-Quality Samples
# ---------------------------------------------------------------------
bad_samples <- uv_smoothed %>%
  group_by(Filename) %>%
  summarise(max_trans = max(smoothed_value, na.rm = TRUE)) %>%
  filter(max_trans > 125)

bad_neg <- uv_smoothed %>%
  group_by(Filename) %>%
  summarise(min_trans = min(smoothed_value, na.rm = TRUE)) %>%
  filter(min_trans < -10)

bad_all <- union(bad_samples$Filename, bad_neg$Filename)

uv_clean <- uv_smoothed %>%
  filter(!Filename %in% bad_all)

# ---------------------------------------------------------------------
# STEP 5: Calculate AUC for UV and IR Ranges
# ---------------------------------------------------------------------



uv_auc_summary <- uv_clean %>%
  filter(Wavelength >= 280, Wavelength <= 400) %>%
  group_by(Filename, species, `window/petal`, plant) %>%
  dplyr::summarise(
    AUC_uv = {
      df <- na.omit(cur_data())
      trapz(df$Wavelength, df$smoothed_value)
    },
    .groups = "drop"
  )

ir_auc_summary <- uv_clean %>%
  filter(Wavelength >= 700, Wavelength <= 800) %>%
  group_by(Filename, species, `window/petal`, plant) %>%
  summarise(
    AUC_ir = {
      df <- na.omit(cur_data())
      trapz(df$Wavelength, df$smoothed_value)
    },
    .groups = "drop"
  )

uv_clean %>%
  filter(Wavelength >= 280, Wavelength <= 400) %>%
  group_by(Filename) %>%
  summarise(n_na = sum(is.na(smoothed_value))) %>%
  filter(n_na > 0)

# ---------------------------------------------------------------------
# STEP 6: Model AUC Differences (UV and IR)
# ---------------------------------------------------------------------
lm_uv <- lm(AUC_UV ~ species * `window/petal` + plant, data = uv_clean_summary)
lm_ir <- lm(AUC_IR ~ species * `window/petal` + plant, data = uv_clean_summary)

# ---------------------------------------------------------------------
# STEP 7: Estimated Marginal Means and Tukey CLDs
# ---------------------------------------------------------------------
emm_uv <- emmeans(lm_uv, ~ species * `window/petal`)
cld_uv <- cld(emm_uv, Letters = letters, adjust = "tukey") %>%
  mutate(group_label = paste(species, `window/petal`, sep = "_"))

emm_ir <- emmeans(lm_ir, ~ species * `window/petal`)
cld_ir <- cld(emm_ir, Letters = letters, adjust = "tukey") %>%
  mutate(group_label = paste(species, `window/petal`, sep = "_"))

# ---------------------------------------------------------------------
# STEP 8: Plot AUC with Estimated Means and CLDs
# ---------------------------------------------------------------------
# UV AUC Plot
uv_plot <- ggplot() +
  geom_jitter(data = uv_clean_summary,
              aes(x = group_label, y = AUC_UV, color = `window/petal`),
              width = 0.15, alpha = 0.3, size = 1) +
  geom_point(data = cld_uv, aes(x = group_label, y = emmean), size = 3) +
  geom_errorbar(data = cld_uv, aes(x = group_label, ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_text(data = cld_uv, aes(x = group_label, y = upper.CL + 5, label = .group)) +
  theme_minimal() +
  labs(title = "UV AUC by Species and Tissue", y = "AUC (280–400 nm)", x = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# IR AUC Plot
ir_plot <- ggplot() +
  geom_jitter(data = uv_clean_summary,
              aes(x = group_label, y = AUC_IR, color = `window/petal`),
              width = 0.15, alpha = 0.3, size = 1) +
  geom_point(data = cld_ir, aes(x = group_label, y = emmean), size = 3) +
  geom_errorbar(data = cld_ir, aes(x = group_label, ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_text(data = cld_ir, aes(x = group_label, y = upper.CL + 5, label = .group)) +
  theme_minimal() +
  labs(title = "IR AUC by Species and Tissue", y = "AUC (700–1100 nm)", x = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(uv_plot)
print(ir_plot)

#### EVERYTHING BELOW HERE IS EXTRA RIGHT NOW#####

# ---------------------------------------------------------------------
# STEP 6: Summarize Mean UV Transmission (280–400 nm) per Sample
# ---------------------------------------------------------------------
uv_summary_stats <- uv_clean %>%
  filter(Wavelength >= 280, Wavelength <= 400) %>%
  group_by(Filename, species, `window/petal`, plant) %>%
  summarise(mean_uv_trans = mean(smoothed_value, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    plant = as.factor(plant),  # Treat as fixed factor
    group_label = paste(species, `window/petal`, sep = "_")
  )

# ---------------------------------------------------------------------
# STEP 7: Linear Model for Species × Tissue Effects
# ---------------------------------------------------------------------
lm_model <- lm(mean_uv_trans ~ species * `window/petal` + plant, data = uv_summary_stats)
summary(lm_model)
anova(lm_model)

# ---------------------------------------------------------------------
# STEP 8: Estimated Marginal Means + Tukey Group Letters
# ---------------------------------------------------------------------
emm <- emmeans(lm_model, ~ species * `window/petal`)
cld_df <- cld(emm, Letters = letters, adjust = "tukey") %>%
  mutate(group_label = paste(species, `window/petal`, sep = "_"))

# ---------------------------------------------------------------------
# STEP 9: Plot Estimated Means with Group Letters + Raw Points
# ---------------------------------------------------------------------
ggplot() +
  # Individual raw sample points
  geom_jitter(
    data = uv_summary_stats,
    aes(x = group_label, y = mean_uv_trans, color = `window/petal`),
    width = 0.15,
    alpha = 0.3,
    size = 1
  ) +
  
  # Model-estimated group means
  geom_point(
    data = cld_df,
    aes(x = group_label, y = emmean, color = `window/petal`),
    width = 0.6
  ) +
  
  # 95% CI error bars
  geom_errorbar(
    data = cld_df,
    aes(x = group_label, ymin = lower.CL, ymax = upper.CL, color = `window/petal`), width=0.2
  ) +
  
  # Group letters above bars
  geom_text(
    data = cld_df,
    aes(x = group_label, y = upper.CL + 2, label = .group),
    size = 6
  ) +
  
  labs(
    title = "UV Transmission by Species and Tissue Type",
    x = "Group (Species × Tissue)",
    y = "Estimated Mean UV Transmission (280–400 nm)",
    fill = "Tissue",
    color = "Tissue"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---------------------------------------------------------------------
# IR ANALYSIS (700–1100 nm)
# ---------------------------------------------------------------------

# Define IR range
IR_MIN <- 500
IR_MAX <- 1100

# STEP 1: Summarize mean IR transmittance per sample
uv_summary_stats_ir <- uv_clean %>%
  filter(Wavelength >= IR_MIN, Wavelength <= IR_MAX) %>%
  group_by(Filename, species, `window/petal`, plant) %>%
  summarise(mean_ir_trans = mean(smoothed_value, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    plant = as.factor(plant),  # treat plant as a fixed factor
    group_label = paste(species, `window/petal`, sep = "_")
  )

# STEP 2: Fit linear model with interaction
lm_model_ir <- lm(mean_ir_trans ~ species * `window/petal` + plant, data = uv_summary_stats_ir)

# View model results
summary(lm_model_ir)
anova(lm_model_ir)

# STEP 3: Estimated marginal means and Tukey groupings
library(emmeans)
library(multcomp)

emm_ir <- emmeans(lm_model_ir, ~ species * `window/petal`)
cld_ir <- cld(emm_ir, Letters = letters, adjust = "tukey") %>%
  mutate(group_label = paste(species, `window/petal`, sep = "_"))

# STEP 4: Final IR plot with points, means, error bars, and letters
ggplot() +
  # Individual samples (low alpha)
  geom_jitter(
    data = uv_summary_stats_ir,
    aes(x = group_label, y = mean_ir_trans, color = `window/petal`),
    width = 0.15,
    alpha = 0.3,
    size = 1
  ) +
  
  # Group means (bars)
  geom_col(
    data = cld_ir,
    aes(x = group_label, y = emmean, fill = `window/petal`),
    color = "black",
    width = 0.6
  ) +
  
  # 95% confidence intervals
  geom_errorbar(
    data = cld_ir,
    aes(x = group_label, ymin = lower.CL, ymax = upper.CL),
    width = 0.2
  ) +
  
  # Tukey group letters
  geom_text(
    data = cld_ir,
    aes(x = group_label, y = upper.CL + 2, label = .group),
    size = 5
  ) +
  
  labs(
    title = "IR Transmission by Species and Tissue Type",
    x = "Group (Species × Tissue)",
    y = "Estimated Mean IR Transmission (700–1100 nm)",
    fill = "Tissue",
    color = "Tissue"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#AUC


#PCA

spectral_matrix <- uv_clean_sep %>%
  filter(Wavelength >= 281.1, Wavelength <= 698) %>%
  dplyr::select(Filename, Wavelength, smoothed_value) %>%
  pivot_wider(names_from = Wavelength, values_from = smoothed_value)

glimpse(spectral_matrix)

spectral_values <- spectral_matrix %>% dplyr::select(-Filename)
pca_result <- prcomp(spectral_values, scale. = TRUE)
metadata <- uv_clean_sep %>%
  distinct(Filename, species, `window/petal`)

spectral_matrix <- spectral_matrix %>%
  mutate(filename = tolower(trimws(Filename)))

metadata <- metadata %>%
  mutate(filename = tolower(trimws(Filename)))
spectral_with_meta <- left_join(spectral_matrix, metadata, by = "filename")

library(ggfortify)

autoplot(pca_result,
         data = spectral_with_meta,
         colour = 'species',
         shape = `window/petal`,
         size = 3) +
  theme_minimal() +
  labs(title = "PCA of UV-Visible Spectra",
       color = "Species",
       shape = "Tissue")
