library(tidyverse)
library(tidymodels)
library(sf)
library(patchwork)

source("code/config.R")
load(paste0(OUTPUT_PATH, "results/03_thresholds_diagnostics.RData"))
load(paste0(OUTPUT_PATH, "results/04_importance_pdp.RData"))

# Load fish data for supplemental maps
#fish_all_cv <- read.csv("data/cleaned/fish_all_cv.csv")

# ── Cross-CV comparison ────────────────────────────────────────────────────
cv_comparison <- results_all %>%
  filter(var_set == "rest") %>%
  select(species, cv_scheme, roc_auc, tss, dist_to_01,
         sens, spec, observed_prevalence) %>%
  pivot_wider(names_from = cv_scheme,
              values_from = c(roc_auc, tss, dist_to_01, sens, spec)) %>%
  mutate(
    auc_drop_spatial = roc_auc_random - roc_auc_spatial,
    auc_drop_cluster = roc_auc_random - roc_auc_cluster
  )

# ── Variable set comparison ────────────────────────────────────────────────
varset_comparison <- results_all %>%
  filter(cv_scheme == "cluster") %>%
  select(species, var_set, roc_auc) #%>%
  #pivot_wider(names_from = var_set, values_from = roc_auc) %>%
  #mutate(patch_improvement = rest - blgl)

# Order species based on 
species_order <- results_all %>% 
  filter(cv_scheme != "random") %>%
  group_by(species) %>%
  summarize(mean_auc = mean(roc_auc, na.rm = TRUE)) %>%
  arrange(mean_auc) %>%
  pull(species)

# Linking CV Observations
df <- results_all %>%
  filter(cv_scheme != "random") %>%
  mutate(
    species_f = factor(species, levels = species_order),
    y_base = as.numeric(species_f),
    y_offset = ifelse(var_set == "blgl", -0.2, 0.2),
    y = y_base + y_offset
  )

seg_df <- df %>%
  select(species_f, var_set, cv_scheme, roc_auc, y) %>%
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
                                "rest" = "#b2df8a")) +
  labs(x = "AUC", y = "Species",
    color = "Variable set", shape = "CV scheme") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor.y = element_blank())

# Linking Var. Set Observations
df <- results_all %>%
  filter(cv_scheme != "random") %>%
  mutate(
    species_f = factor(species, levels = species_order),
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
  labs(x = "AUC", y = "Species",
       color = "CV Scheme", shape = "Variable Set") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor.y = element_blank())


# ── Cross-validation comparison ─────────────────────────────────
random_auc <- results_all %>%
  filter(cv_scheme == "random") %>%
  select(species, random_auc = roc_auc, var_set)
  
results_diff <- left_join(results_all, random_auc) %>%
  filter(!cv_scheme %in% c("random", "cluster")) %>%
  mutate(auc_drop = random_auc - roc_auc) %>%
  rename(spatial_auc = roc_auc) %>%
  mutate(
    auc_drop_z  = (auc_drop - mean(auc_drop)) / sd(auc_drop),
    flag = auc_drop_z > 1)


# ── Species clustering by habitat affinity ─────────────────────────────────
cluster_matrix <- signed_importance %>%
  pivot_wider(id_cols = species, names_from = variable,
              values_from = signed_imp) %>%
  column_to_rownames("species") %>%
  replace(is.na(.), 0)

dist_mat   <- dist(cluster_matrix, method = "euclidean")
hclust_res <- hclust(dist_mat, method = "ward.D2")
species_clusters <- cutree(hclust_res, k = 4) %>%
  enframe(name = "species", value = "habitat_guild")

# ── Figures ────────────────────────────────────────────────────────────────
# (your existing plotting code for calibration panels,
#  importance heatmaps, PDP overlays, CV comparisons)

# Variable comparison plot

ggplot(varset_comparison)+
  #geom_point(aes(reorder(`Common Name`, roc_auc), roc_auc, color = vars), size = 3)+
  geom_point(aes(roc_auc, reorder(species, roc_auc), color = var_set), size = 3)+
  geom_vline(aes(xintercept = .9), linetype = 2)+
  geom_vline(aes(xintercept = .8), linetype = 2)+
  geom_vline(aes(xintercept = .7), linetype = 2)+
  labs(x = "", y = "Model Performance", color = "Model type")+
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom")













save(cv_comparison, varset_comparison, 
     species_clusters, hclust_res,
     file = paste0(OUTPUT_PATH, "results/05_synthesis.RData"))
cat("05 complete.\n")