library(tidymodels)
library(finetune)
library(future)
library(doFuture)

options(future.globals.maxSize = 1.5 * 1024^3)
plan(multisession, workers = parallel::detectCores() - 1)

registerDoFuture()

source("code/config.R")
load(paste0(OUTPUT_PATH, "models/01_folds.RData"))

# ── Core tuning function ────────────────────────────────────────────────────
run_one_model_type <- function(sp_list, folds_list, var_set,
                               model_label, n_trees = N_TREES) {
  
  var_cols <- if (var_set == "rest") VARS_REST else VARS_BLGL
  
  race_results_list <- vector("list", length(sp_list))
  names(race_results_list) <- sp_list
  
  for (i in seq_along(sp_list)) {
    sp         <- sp_list[i]
    fish_train <- folds_list[[sp]]$data %>%
      dplyr::select(all_of(var_cols))
    fish_folds <- folds_list[[sp]]$folds
    
    cat(model_label, sp, format(Sys.time(), "%H:%M:%S"), "\n")
    
    dummy_rec <- recipe(occ ~ ., data = fish_train) %>%
      update_role(cluster, new_role = "id") %>%
      step_impute_bag(all_predictors())
    
    rf_spec <- rand_forest(mtry = tune(), min_n = tune(), 
                           trees = n_trees) %>%
      set_engine("ranger", num.threads = 1,
                 importance = "impurity", probability = TRUE) %>%
      set_mode("classification")
    
    rf_params <- parameters(finalize(mtry(), fish_train), min_n())
    base_grid <- grid_space_filling(rf_params, size = GRID_SIZE)
    
    my_metrics <- metric_set(roc_auc, brier_class)
    race_ctrl  <- control_race(save_pred = TRUE, 
                               parallel_over = "everything",
                               save_workflow = TRUE)
    wf <- workflow() %>% add_recipe(dummy_rec) %>% add_model(rf_spec)
    
    result <- tryCatch(
      { r <- tune_race_anova(wf, resamples = fish_folds,
                             grid = base_grid, metrics = my_metrics,
                             control = race_ctrl)
      list(tune_result = r, tune_method = "race_anova") },
      error = function(e) {
        message("race_anova failed, trying race_win_loss...")
        tryCatch(
          { r <- tune_race_win_loss(wf, resamples = fish_folds,
                                    grid = base_grid, metrics = my_metrics,
                                    control = race_ctrl)
          list(tune_result = r, tune_method = "race_win_loss") },
          error = function(e2) {
            message("race_win_loss failed, using tune_grid...")
            r <- tune_grid(wf, resamples = fish_folds,
                           grid = base_grid, metrics = my_metrics,
                           control = control_grid(save_pred = TRUE,
                                                  save_workflow = TRUE))
            list(tune_result = r, tune_method = "tune_grid")
          }
        )
      }
    )
    
    race_results_list[[i]] <- result
  }
  
  race_results_list
}

# ── Run all six model types ─────────────────────────────────────────────────
for (mt in names(MODEL_TYPES)[3:6]) {
  cfg        <- MODEL_TYPES[[mt]]
  folds_list <- switch(cfg$cv,
                       spatial = spatial_folds_list,
                       cluster = cluster_folds_list,
                       random  = random_folds_list
  )
  
  cat("\n========", mt, "========\n")
  race_results <- run_one_model_type(
    sp_list     = sp_list,
    folds_list  = folds_list,
    var_set     = cfg$vars,
    model_label = mt
  )
  
  save(race_results, sp_list,
       file = paste0(OUTPUT_PATH, "models/02_", mt, ".RData"))
  cat(mt, "saved.\n")
}
