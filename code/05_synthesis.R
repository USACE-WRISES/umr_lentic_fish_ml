library(tidyverse)
library(tidymodels)
library(sf)
library(patchwork)
library(here)
library(grid) # For textGrob in Figure 2

source("code/config.R")
load(paste0(OUTPUT_PATH, "results/03_thresholds_diagnostics.RData"))
load(paste0(OUTPUT_PATH, "results/04_importance_pdp.RData"))

sp_guilds <- read_csv(here("data", "sp_guilds.csv")) %>%
  rename_with(~ gsub(" ", "_", .x)) %>%
  mutate(Common_Name = str_to_title(Common_Name))

results_all <- left_join(results_all, sp_guilds)
cv_preds_all <- left_join(cv_preds_all, sp_guilds)

results_tab <- select(results_all, Common_Name, Scientific_Name, "Species Code" = species,  cv_scheme, var_set, mtry, min_n, roc_auc) %>%
  mutate(cv_scheme = recode(cv_scheme, "spatial" = "Pool-based", "cluster" = "Lake-cluster", "random" = "Random"),
         var_set = recode(var_set, "blgl" = "Bluegill", "rest" = "Restoration"),
         roc_auc = signif(round(roc_auc, 2), 2)) %>%
  rename("Common Name" =Common_Name, "Scientific Name" = Scientific_Name, "CV Scheme"=cv_scheme, "Variable Set"=var_set, "Minimum N"=min_n, "AUC"=roc_auc)

write_csv(results_tab, "output/tabS2_results_tab.csv")

# ── Calibration plot ────────────────────────────────────────────────────

calibration_check <- cv_preds_all %>%
  filter(model_type == "spatial_rest") %>%
  group_by(species) %>%
  mutate(bin = ntile(.pred_1, 10)) %>% 
  group_by(species, Common_Name, bin) %>%
  summarise(
    mean_pred = mean(.pred_1),
    mean_obs  = mean(occ == "1"),
    n         = n(),
    .groups   = "drop"
  ) #%>%
  #mutate(Common_Name = str_replace_all(Common_Name, " ", "\n"))

ggplot(calibration_check, aes(x = mean_pred, y = mean_obs)) +
  geom_point() +
  geom_abline(linetype = "dashed") +
  geom_smooth(method = "loess", se = FALSE, colour = "steelblue") +
  #scale_size_continuous(range = c(1, 5)) +
  labs(x = "Mean predicted probability",
       y = "Observed fraction present") +
  facet_wrap(~Common_Name, 
             ncol = 4)+
  #coord_equal()+
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave("output/figures/figS3_cal_plots.jpeg", width = 7, height = 9)

# ── Threshold evaluation  ────────────────────────────────────────────────────
results_all %>%
  filter(model_type == "spatial_rest") %>%
ggplot(aes(spec, sens)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = Common_Name), size = 3) +
  labs(y = "Sensitivity", x = "Specificity") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank())
  
ggsave("output/figures/figS4_sens_spec.jpeg", width = 8, height = 7)

# ── Cross-CV comparison ────────────────────────────────────────────────────
cv_comparison <- results_all %>%
  filter(var_set == "rest") %>%
  dplyr::select(Common_Name, cv_scheme, roc_auc, tss, dist_to_01,
         sens, spec, observed_prevalence) %>%
  pivot_wider(names_from = cv_scheme,
              values_from = c(roc_auc, tss, dist_to_01, sens, spec)) %>%
  mutate(
    auc_drop_spatial = roc_auc_random - roc_auc_spatial,
    auc_drop_cluster = roc_auc_random - roc_auc_cluster
  )

median(cv_comparison$auc_drop_spatial)
range(cv_comparison$auc_drop_spatial)

median(cv_comparison$auc_drop_cluster)
range(cv_comparison$auc_drop_cluster)

# ── Variable set comparison ────────────────────────────────────────────────
varset_comparison <- results_all %>%
  filter(cv_scheme == "cluster") %>%
  dplyr::select(Common_Name, var_set, roc_auc) #%>%
  #pivot_wider(names_from = var_set, values_from = roc_auc) %>%
  #mutate(patch_improvement = rest - blgl)

# Order Common_Name based on 
species_order <- results_all %>% 
  filter(cv_scheme != "random") %>%
  group_by(Common_Name) %>%
  summarize(mean_auc = mean(roc_auc, na.rm = TRUE)) %>%
  arrange(mean_auc) %>%
  pull(Common_Name)

# Linking CV Observations
df <- results_all %>%
  filter(cv_scheme != "random") %>%
  mutate(
    species_f = factor(Common_Name, levels = species_order),
    y_base = as.numeric(species_f),
    y_offset = ifelse(var_set == "blgl", -0.2, 0.2),
    y = y_base + y_offset
  )

seg_df <- df %>%
  dplyr::select(species_f, var_set, cv_scheme, roc_auc, y) %>%
  tidyr::pivot_wider(
    names_from = cv_scheme,
    values_from = roc_auc
  )

ggplot() +
  geom_hline(yintercept = seq(1.5, length(species_order) - 0.5, by = 1),   # background separators
             color = "grey50", linewidth = 0.5) +
  geom_vline(aes(xintercept = .9), linetype = 2)+
  geom_vline(aes(xintercept = .8), linetype = 2)+
  geom_vline(aes(xintercept = .7), linetype = 2)+
  geom_segment(data = seg_df,   # segments
               aes(x = spatial,
                   xend = cluster,
                   y = y,
                   yend = y,
                   color = var_set),
    linewidth = 0.8, alpha = 0.7) +
  geom_point(data = df,
             aes(x = roc_auc, y = y, color = var_set, shape = cv_scheme),
             size = 3)+
  scale_y_continuous(breaks = seq_along(species_order), 
                     labels = species_order,
                     expand = c(0.015,0.015)) + # To deal with weird spacing at top and bottom of plot
  scale_color_manual(values = c("blgl" = "#1f78b4",
                                "rest" = "#b2df8a"),
                     labels = c("blgl" = "Bluegill",
                                "rest" = "Restoration")) +
  scale_shape_manual(values = c(16,17), labels = c("cluster" = "Cluster",
                                "spatial" = "Pool")) +
  labs(x = "AUC", y = "Species",
    color = "Variable set", shape = "CV scheme") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor.y = element_blank())+
  # 1. Turn off clipping so annotations can draw outside the panel boundary
  coord_cartesian(clip = "off") +
  # 2. Add text using viewport (npc) coordinates relative to the plot area
  annotation_custom(textGrob("Poor", x = 0.16, y = 1.03,  gp = gpar(fontface = "italic"))) +
  annotation_custom(textGrob("Useful", x = 0.57, y = 1.03, gp = gpar(fontface = "italic"))) +
  annotation_custom(textGrob("Excellent", x = 0.92, y = 1.03, gp = gpar(fontface = "italic"))) +
  # 3. Pad the outer margins so the text doesn't hit the image file edge
  theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10))

ggsave("output/figures/fig2_mod_comp.jpeg", width = 9, height = 6.5)

# Version for presentations
ggplot() +
  geom_hline(yintercept = seq(1.5, length(species_order) - 0.5, by = 1),   # background separators
             color = "grey50", linewidth = 0.5) +
  geom_vline(aes(xintercept = .9), linetype = 2)+
  geom_vline(aes(xintercept = .8), linetype = 2)+
  geom_vline(aes(xintercept = .7), linetype = 2)+
  # geom_segment(data = seg_df,   # segments
  #              aes(x = spatial,
  #                  xend = cluster,
  #                  y = y,
  #                  yend = y,
  #                  color = var_set),
  #              linewidth = 0.8, alpha = 0.7) +
  geom_point(data = filter(df, cv_scheme == "cluster"),
             aes(x = roc_auc, y = y, color = var_set, shape = var_set),
             size = 3)+
  scale_y_continuous(breaks = seq_along(species_order), 
                     labels = species_order,
                     expand = c(0.015,0.015)) + # To deal with weird spacing at top and bottom of plot
  scale_color_manual(values = c("blgl" = "#1f78b4",
                                "rest" = "#b2df8a"),
                     labels = c("blgl" = "Bluegill",
                                "rest" = "Restoration")) +
  scale_shape_manual(values = c(16,17), labels = c("blgl" = "Bluegill",
                                                   "rest" = "Restoration")) +
  labs(x = "AUC", y = "Species",
       color = "Variable set", shape = "Variable set") +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "bottom")

# Linking Var. Set Observations
df <- results_all %>%
  filter(cv_scheme != "random") %>%
  mutate(
    species_f = factor(Common_Name, levels = species_order),
    y_base = as.numeric(species_f),
    y_offset = ifelse(cv_scheme == "spatial", -0.2, 0.2),
    y = y_base + y_offset
  )

seg_df <- df %>%
  select(species_f, var_set, cv_scheme, roc_auc, y) %>%
  tidyr::pivot_wider(
    names_from = var_set,
    values_from = roc_auc
  )

ggplot() +
  geom_hline(yintercept = seq(1.5, length(species_order) - 0.5, by = 1),   # background separators
             color = "grey50", linewidth = 0.5) +
  geom_vline(aes(xintercept = .9), linetype = 2)+
  geom_vline(aes(xintercept = .8), linetype = 2)+
  geom_vline(aes(xintercept = .7), linetype = 2)+
  geom_segment(data = seg_df,   # segments
               aes(x = blgl,
                   xend = rest,
                   y = y,
                   yend = y,
                   color = cv_scheme),
               linewidth = 0.8, alpha = 0.7) +
  geom_point(data = df,
             aes(x = roc_auc, y = y, color = cv_scheme, shape = var_set),
             size = 3)+
  scale_y_continuous(breaks = seq_along(species_order), 
                     labels = species_order,
                     expand = c(0.015,0.015)) + # To deal with weird spacing at top and bottom of plot
  scale_color_manual(values = c("cluster" = "#1f78b4",
                                "spatial" = "#b2df8a"),
                     labels = c("Cluster", "Spatial")) +
  scale_shape_discrete(labels = c("Bluegill", "Restoration")) +
  labs(x = "AUC", y = "Common_Name",
       color = "CV Scheme", shape = "Variable Set") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor.y = element_blank())

# ── Cross-validation comparison ─────────────────────────────────
random_auc <- results_all %>%
  filter(cv_scheme == "random") %>%
  select(Common_Name, random_auc = roc_auc, var_set)
  
results_diff <- left_join(results_all, random_auc) %>%
  filter(!cv_scheme %in% c("random", "cluster")) %>%
  mutate(auc_drop = random_auc - roc_auc) %>%
  rename(spatial_auc = roc_auc) %>%
  mutate(
    auc_drop_z  = (auc_drop - mean(auc_drop)) / sd(auc_drop),
    flag = auc_drop_z > 1)

# ── Figures ────────────────────────────────────────────────────────────────

# Figure S5, Performance by prevalence
prev_auc <- results_all %>%
  filter(model_type %in% c("spatial_rest")) %>% select(observed_prevalence, roc_auc)

cor(prev_auc$observed_prevalence, prev_auc$roc_auc)

results_all %>%
  filter(model_type %in% c("spatial_rest")) %>%
  ggplot(aes(observed_prevalence, roc_auc))+
  geom_point()+
  ggrepel::geom_text_repel(aes(label = Common_Name), size = 3) +
  labs(x = "Observed prevalence", y = "AUC") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank())

ggsave("output/figures/figS5_prev_auc.jpeg", width = 8, height = 6.5)

# Figure S4, Cross-validation comparison
pool_rest <- results_all %>%
  filter(model_type %in% c("spatial_rest"))

random_rest <- results_all %>%
  filter(model_type %in% c("random_rest")) %>%
  select(species, random_roc_auc = roc_auc, var_set) %>% 
  right_join(pool_rest) %>%
  mutate(auc_drop = random_roc_auc - roc_auc) %>%
  rename(spatial_roc_auc = roc_auc) %>%
  select(species, spatial_roc_auc, random_roc_auc, auc_drop, var_set) %>% 
  mutate(
    auc_drop_z  = (auc_drop - mean(auc_drop)) / sd(auc_drop),
    flag = auc_drop_z > 1) 

pool_blgl <- results_all %>%
  filter(model_type %in% c("spatial_blgl"))

random_blgl <- results_all %>%
  filter(model_type %in% c("random_blgl")) %>%
  select(species, random_roc_auc = roc_auc, var_set) %>% 
  right_join(pool_blgl) %>%
  mutate(auc_drop = random_roc_auc - roc_auc) %>%
  rename(spatial_roc_auc = roc_auc) %>%
  select(species, spatial_roc_auc, random_roc_auc, auc_drop, var_set) %>% 
  mutate(
    auc_drop_z  = (auc_drop - mean(auc_drop)) / sd(auc_drop),
    flag = auc_drop_z > 1)

all_cv <- bind_rows(random_rest, random_blgl) %>%
  mutate(var_set = recode(var_set, 
                          "blgl" = "Bluegill variables",
                          "rest" = "Restoration variables")) %>% 
  left_join(sp_guilds)

spatial_names <- filter(all_cv, flag == "TRUE")

# Spatial comparisons
ggplot(all_cv, aes(x = spatial_roc_auc, y = random_roc_auc, label = Common_Name)) +
  geom_abline(linetype = "dashed") +
  geom_point(aes(colour = flag), size = 3) +
  ggrepel::geom_text_repel(data = spatial_names, size = 3) +
  scale_colour_manual(values = c("FALSE" = "grey40", "TRUE" = "firebrick"),
                      labels = c("Acceptable", "High spatial dependency"),
                      name   = NULL) +
  facet_wrap(~var_set)+
  labs(x = "Pool CV AUC", y = "Random CV AUC") +
  coord_equal(xlim = c(0.62, .95), ylim = c(0.62, .95)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

ggsave("output/figures/figS6_CV_comp.jpeg", width = 9, height = 7)

save(cv_comparison, varset_comparison, 
     file = paste0(OUTPUT_PATH, "results/05_synthesis.RData"))
cat("05 complete.\n")