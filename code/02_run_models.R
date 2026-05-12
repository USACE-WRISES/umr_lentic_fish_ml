source("code/config.R")

library(tidymodels)
library(finetune)
library(future)
library(furrr)

load(paste0(OUTPUT_PATH, "models/01_folds.RData"))

# ── Set up parallel processing ──────────────────────────────────────────────
options(future.globals.maxSize = 1.5 * 1024^3)

N_CORES <- parallel::detectCores() - 1  # leave one core free
plan(multisession, workers = N_CORES)

cat("Using", N_CORES, "cores for parallel processing\n")

# ── Core tuning function (unchanged) ────────────────────────────────────────
run_one_species <- function(sp, fish_train, fish_folds, var_cols,
                            model_label, n_trees = N_TREES) {
  
  # Subset to requested variables
  fish_train <- fish_train %>% dplyr::select(all_of(var_cols))
  print(table(fish_train$occ))
  print(str(fish_train$occ))
  
  dummy_rec <- recipe(occ ~ ., data = fish_train) %>%
    update_role(cluster, new_role = "id") %>%
    step_impute_bag(all_predictors())
  
  rf_spec <- rand_forest(mtry = tune(), min_n = tune(), 
                         trees = n_trees) %>%
    set_engine("ranger", num.threads = 1,  # IMPORTANT: 1 thread per worker
               importance = "impurity", probability = TRUE) %>%
    set_mode("classification")
  
  rf_params <- parameters(finalize(mtry(), fish_train), min_n())
  base_grid <- grid_space_filling(rf_params, size = GRID_SIZE)
  
  my_metrics <- metric_set(roc_auc, brier_class)
  race_ctrl  <- control_race(
    save_pred     = TRUE,
    parallel_over = "resamples",  # parallel within species
    save_workflow = TRUE
  )
  
  wf <- workflow() %>% add_recipe(dummy_rec) %>% add_model(rf_spec)
  
  result <- tryCatch(
    { r <- tune_race_anova(wf, resamples = fish_folds,
                           grid = base_grid, metrics = my_metrics,
                           control = race_ctrl)
    list(tune_result = r, tune_method = "race_anova") },
    error = function(e) {
      message(sp, ": race_anova failed, trying race_win_loss...")
      tryCatch(
        { r <- tune_race_win_loss(wf, resamples = fish_folds,
                                  grid = base_grid, metrics = my_metrics,
                                  control = race_ctrl)
        list(tune_result = r, tune_method = "race_win_loss") },
        error = function(e2) {
          message(sp, ": race_win_loss failed, using tune_grid...")
          r <- tune_grid(wf, resamples = fish_folds,
                         grid = base_grid, metrics = my_metrics,
                         control = control_grid(save_pred = TRUE,
                                                save_workflow = TRUE))
          list(tune_result = r, tune_method = "tune_grid")
        }
      )
    }
  )
  
  list(
    species     = sp,
    tune_result = result$tune_result,
    tune_method = result$tune_method
  )
}

# ── Run all six model types ─────────────────────────────────────────────────
for (mt in names(MODEL_TYPES)) {
  cfg        <- MODEL_TYPES[[mt]]
  folds_list <- switch(cfg$cv,
                       spatial = spatial_folds_list,
                       cluster = cluster_folds_list,
                       random  = random_folds_list
  )
  
  var_cols <- if (cfg$vars == "rest") VARS_REST else VARS_BLGL
  
  cat("\n========", mt, "========\n")
  cat("Starting", length(sp_list), "species at", 
      format(Sys.time(), "%H:%M:%S"), "\n")
  
  # Parallel map over species
  race_results <- future_map(
    sp_list,
    function(sp) {
      run_one_species(
        sp          = sp,
        fish_train  = folds_list[[sp]]$data,
        fish_folds  = folds_list[[sp]]$folds,
        var_cols    = var_cols,
        model_label = mt,
        n_trees     = N_TREES
      )
    },
    .options = furrr_options(seed = TRUE),  # reproducibility
    .progress = TRUE  # progress bar across species
  )
  
  # Convert list to named list by species
  names(race_results) <- sp_list
  
  save(race_results, sp_list,
       file = paste0(OUTPUT_PATH, "models/02_", mt, ".RData"))
  
  cat(mt, "complete at", format(Sys.time(), "%H:%M:%S"), "\n")
}

# Clean up parallel workers
plan(sequential)
