source("config.R")
load(paste0(OUTPUT_PATH, "models/01_folds.RData"))

# ── Threshold optimization function ────────────────────────────────────────
optimize_threshold <- function(cv_preds, method = THRESH_METHOD) {
  
  thresh_perf <- tibble(threshold = THRESH_SEQ) %>%
    mutate(
      sens = map_dbl(threshold, function(t) {
        preds <- cv_preds %>%
          mutate(pred = factor(ifelse(.pred_1 >= t, "1","0"),
                               levels = c("1","0")))
        yardstick::sens(preds, truth = occ, estimate = pred)$.estimate
      }),
      spec = map_dbl(threshold, function(t) {
        preds <- cv_preds %>%
          mutate(pred = factor(ifelse(.pred_1 >= t, "1","0"),
                               levels = c("1","0")))
        yardstick::spec(preds, truth = occ, estimate = pred)$.estimate
      })
    ) %>%
    mutate(
      tss          = sens + spec - 1,
      dist_to_01   = sqrt((1 - sens)^2 + (1 - spec)^2),
      sens_spec_diff = abs(sens - spec)
    )
  
  opt_col <- switch(method,
                    dist_to_01     = "dist_to_01",
                    tss            = "tss",
                    sens_spec_diff = "sens_spec_diff"
  )
  
  # Minimise for dist_to_01 and sens_spec_diff, maximise for tss
  if (method == "tss") {
    best <- thresh_perf %>%
      filter(.data[[opt_col]] == max(.data[[opt_col]], na.rm = TRUE))
  } else {
    best <- thresh_perf %>%
      filter(.data[[opt_col]] == min(.data[[opt_col]], na.rm = TRUE))
  }
  
  list(
    threshold     = median(best$threshold),
    threshold_min = min(best$threshold),
    threshold_max = max(best$threshold),
    thresh_perf   = thresh_perf
  )
}

# ── Process all model types ─────────────────────────────────────────────────
all_results  <- list()
all_cv_preds <- list()

for (mt in names(MODEL_TYPES)) {
  
  cat("\nProcessing thresholds:", mt, "\n")
  load(paste0(OUTPUT_PATH, "models/02_", mt, ".RData"))
  
  results_mt  <- list()
  cv_preds_mt <- list()
  
  for (i in seq_along(sp_list)) {
    sp         <- sp_list[i]
    tune_res   <- race_results[[i]]$tune_result
    tune_meth  <- race_results[[i]]$tune_method
    fish_train <- spatial_folds_list[[sp]]$data  # always use full training data
    
    best_params <- tune_res %>% select_best(metric = "roc_auc")
    cv_preds    <- tune_res %>% collect_predictions(parameters = best_params)
    
    # Threshold optimization
    thresh_res  <- optimize_threshold(cv_preds, method = THRESH_METHOD)
    best_thresh <- thresh_res$threshold
    
    cv_preds_final <- cv_preds %>%
      mutate(pred_opt = factor(ifelse(.pred_1 >= best_thresh, "1","0"),
                               levels = c("1","0")))
    
    # Tuning metrics
    tuning_metrics <- cv_preds %>%
      summarise(
        roc_auc = yardstick::roc_auc_vec(truth = occ, estimate = .pred_1),
        brier   = yardstick::brier_class_vec(truth = occ, estimate = .pred_1)
      )
    
    # Threshold metrics
    thresh_metrics <- cv_preds_final %>%
      summarise(
        sens         = yardstick::sens_vec(truth = occ, estimate = pred_opt),
        spec         = yardstick::spec_vec(truth = occ, estimate = pred_opt),
        bal_accuracy = yardstick::bal_accuracy_vec(truth = occ, estimate = pred_opt),
        mcc          = yardstick::mcc_vec(truth = occ, estimate = pred_opt),
        tss          = sens + spec - 1,
        dist_to_01   = sqrt((1 - sens)^2 + (1 - spec)^2)
      )
    
    # Calibration diagnostics
    cal_diag <- cv_preds %>%
      mutate(bin = ntile(.pred_1, 10)) %>%
      group_by(bin) %>%
      summarise(mean_pred = mean(.pred_1),
                mean_obs  = mean(occ == "1"),
                n         = n(), .groups = "drop") %>%
      filter(n >= 5) %>%
      summarise(
        cal_monotonicity = cor(mean_pred, mean_obs, method = "spearman"),
        n_bins           = n(),
        max_pred         = max(mean_pred),
        n_inversions     = sum(diff(mean_obs) < -0.05),
        max_drop         = ifelse(n() > 1, abs(min(diff(mean_obs))), 0)
      ) %>%
      mutate(
        cal_flag = case_when(
          n_bins <= 3 & max_pred < 0.40              ~ "compressed_range",
          cal_monotonicity >= 0.80 & n_inversions == 0 ~ "good",
          cal_monotonicity >= 0.60 & n_inversions <= 1 &
            max_drop < 0.10                           ~ "acceptable",
          n_inversions >= 2 | max_drop >= 0.10        ~ "non_monotone",
          TRUE                                        ~ "uncertain"
        )
      )
    
    obs_prev       <- mean(cv_preds$occ == "1")
    mean_pred_prob <- mean(cv_preds$.pred_1)
    pred_prev      <- mean(cv_preds_final$pred_opt == "1")
    
    results_mt[[sp]] <- bind_cols(tuning_metrics, thresh_metrics) %>%
      mutate(
        optimal_threshold    = best_thresh,
        threshold_min        = thresh_res$threshold_min,
        threshold_max        = thresh_res$threshold_max,
        observed_prevalence  = obs_prev,
        mean_predicted_prob  = mean_pred_prob,
        prop_classified_pres = pred_prev,
        species              = sp,
        model_type           = mt,
        cv_scheme            = MODEL_TYPES[[mt]]$cv,
        var_set              = MODEL_TYPES[[mt]]$vars,
        tune_method          = tune_meth,
        mtry                 = best_params$mtry,
        min_n                = best_params$min_n
      ) %>%
      bind_cols(cal_diag)
    
    cv_preds_mt[[sp]] <- cv_preds_final %>%
      mutate(species = sp, model_type = mt)
  }
  
  all_results[[mt]]  <- bind_rows(results_mt)
  all_cv_preds[[mt]] <- bind_rows(cv_preds_mt)
}

results_all  <- bind_rows(all_results)
cv_preds_all <- bind_rows(all_cv_preds)

save(results_all, cv_preds_all,
     file = paste0(OUTPUT_PATH, "results/03_thresholds_diagnostics.RData"))
cat("03 complete.\n")