source("config.R")
load(paste0(OUTPUT_PATH, "models/01_folds.RData"))
load(paste0(OUTPUT_PATH, "results/03_thresholds_diagnostics.RData"))
load(paste0(OUTPUT_PATH, "results/04_importance_pdp.RData"))
load(paste0(OUTPUT_PATH, "results/05_synthesis.RData"))

# Pull thresholds for primary model type
primary_results <- results_all %>%
  filter(model_type == "spatial_rest")

best_thresh_list <- setNames(
  as.list(primary_results$optimal_threshold),
  primary_results$species
)

# Assemble everything the Shiny app needs
save(
  # Core prediction objects
  final_rf_list,
  best_thresh_list,
  sp_list,
  
  # Performance for display
  results_all,
  cv_comparison,
  
  # Diagnostics for display
  importance_df_all,
  pdp_all,
  signed_importance,
  species_clusters,
  
  file = paste0(OUTPUT_PATH, "shiny/shiny_objects.RData")
)

cat("Shiny objects saved. File size:",
    round(file.size(paste0(OUTPUT_PATH, "shiny/shiny_objects.RData")) / 1e6, 1),
    "MB\n")