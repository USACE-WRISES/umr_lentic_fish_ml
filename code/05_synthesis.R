source("config.R")
load(paste0(OUTPUT_PATH, "results/03_thresholds_diagnostics.RData"))
load(paste0(OUTPUT_PATH, "results/04_importance_pdp.RData"))

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
  filter(cv_scheme == "spatial") %>%
  select(species, var_set, roc_auc) %>%
  pivot_wider(names_from = var_set, values_from = roc_auc) %>%
  mutate(patch_improvement = all - rest)

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

save(cv_comparison, varset_comparison, 
     species_clusters, hclust_res,
     file = paste0(OUTPUT_PATH, "results/05_synthesis.RData"))
cat("05 complete.\n")