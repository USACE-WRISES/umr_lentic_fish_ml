library(tidymodels)

source("code/config.R")
load(paste0(OUTPUT_PATH, "models/01_folds.RData"))
load(paste0(OUTPUT_PATH, "results/03_thresholds_diagnostics.RData"))

# Run importance and PDPs only for primary model type
PRIMARY_MODEL <- "spatial_rest"
load(paste0(OUTPUT_PATH, "models/02_", PRIMARY_MODEL, ".RData"))

importance_df_all <- tibble()
pdp_all           <- tibble()
final_rf_list     <- vector("list", length(sp_list))
names(final_rf_list) <- sp_list

for (i in seq_along(sp_list)) {
  sp         <- sp_list[i]
  fish_train <- cluster_folds_list[[sp]]$data %>%
    dplyr::select(all_of(VARS_REST))
  
  cat("Importance/PDP:", sp, format(Sys.time(), "%H:%M:%S"), "\n")
  
  best_params <- race_results[[i]]$tune_result %>%
    select_best(metric = "roc_auc")
  
  dummy_rec <- recipe(occ ~ ., data = fish_train) %>%
    update_role(cluster, new_role = "id") %>%
    step_impute_bag(all_predictors())
  
  rf_final_spec <- rand_forest(
    mtry = best_params$mtry, min_n = best_params$min_n,
    trees = N_TREES) %>%
    set_engine("ranger", probability = TRUE,
               importance = IMPORTANCE, num.threads = 1) %>%
    set_mode("classification")
  
  final_rf <- workflow() %>%
    add_recipe(dummy_rec) %>%
    add_model(rf_final_spec) %>%
    fit(data = fish_train)
  
  final_rf_list[[i]] <- final_rf
  rf_model           <- extract_fit_engine(final_rf)
  
  # Variable importance
  rf_importance <- ranger::importance(rf_model)
  importance_df <- tibble(variable  = names(rf_importance),
                          importance = rf_importance,
                          species   = sp)
  importance_df_all <- bind_rows(importance_df_all, importance_df)
  
  # PDPs with quantile grids
  top_vars <- names(sort(rf_importance, decreasing = TRUE))
  
  pdp_list_sp <- lapply(top_vars, function(v) {
    pd <- tryCatch({
      if (v %in% CONTINUOUS_VARS) {
        obs_vals  <- fish_train[[v]][!is.na(fish_train[[v]])]
        #grid_vals <- quantile(obs_vals, probs = seq(0.01, 0.99, length.out = 25))
        q <- quantile(obs_vals, probs = c(0.01, 0.99), na.rm = TRUE)
        grid_vals <- seq(q[1], q[2], length.out = 25)
        pred_grid <- setNames(data.frame(as.numeric(grid_vals)), v)
        pdp::partial(object = rf_model, pred.var = v, train = fish_train,
                     prob = TRUE, which.class = "1", pred.grid = pred_grid)
      } else {
        pdp::partial(object = rf_model, pred.var = v, train = fish_train,
                     prob = TRUE, which.class = "1", grid.resolution = 25)
      }
    }, error = function(e) {
      message("PDP failed: ", sp, " ", v, ": ", e$message); NULL
    })
    if (!is.null(pd)) { pd$variable <- v; pd$species <- sp }
    pd
  })
  
  pdp_all <- bind_rows(pdp_all, 
                       bind_rows(Filter(Negate(is.null), pdp_list_sp)))
}

# # Function to extract features from PDP plots
# extract_pdp_features <- function(pd_df, group_keys) {
#   var_name <- group_keys$variable  # get variable name from group keys, not from pd_df
#   
#   continuous_vars <- c("temp", "depth", "current", "do")
#   
#   if (!var_name %in% continuous_vars) {
#     return(tibble(
#       cor_linear     = NA_real_,
#       peak_location  = NA_real_,
#       response_range = max(pd_df$yhat) - min(pd_df$yhat)
#     ))
#   }
#   
#   x <- pd_df[[var_name]]
#   y <- pd_df$yhat
#   
#   # Guard against all-NA or single-value x
#   if (all(is.na(x)) || length(unique(x)) < 3) {
#     return(tibble(
#       cor_linear     = NA_real_,
#       peak_location  = NA_real_,
#       response_range = max(y) - min(y)
#     ))
#   }
#   
#   tibble(
#     cor_linear     = cor(x, y, method = "spearman"),
#     peak_location  = (x[which.max(y)] - min(x)) / (max(x) - min(x)),
#     response_range = max(y) - min(y)
#   )
# }
# 
# # Signed importance for clustering
# pdp_features <- pdp_all %>%
#   group_by(species, variable) %>%
#   group_modify(extract_pdp_features) %>%   # your existing function
#   ungroup()
# 
# signed_importance <- importance_df_all %>%
#   group_by(species) %>%
#   mutate(importance_scaled = importance / sum(importance)) %>%
#   ungroup() %>%
#   left_join(pdp_features %>% select(species, variable, cor_linear),
#             by = c("species","variable")) %>%
#   mutate(signed_imp = case_when(
#     variable %in% CONTINUOUS_VARS  ~ importance_scaled * sign(cor_linear),
#     variable %in% CATEGORICAL_VARS ~ importance_scaled,
#     TRUE                           ~ importance_scaled
#   ))

save(final_rf_list, importance_df_all, pdp_all,
     #pdp_features, signed_importance,
     file = paste0(OUTPUT_PATH, "results/04_importance_pdp.RData"))
cat("04 complete.\n")
