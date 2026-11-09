# =============================================================================
# R/models.R
# Model fitting for the Maine municipal housing-permit analysis.
#
# Task 1 (regression): explain the municipal building-permit *rate* (permits
#   per 1,000 residents) from municipal characteristics -- income, home value,
#   population density, tenure and housing-stock composition -- with town size
#   and county as controls.
# Task 2 (classification): which towns permitted *any* multifamily (2+ unit)
#   housing, using the same characteristics.
#
# All models are compared on identical repeated 10-fold cross-validation folds.
# Tuned models are tuned inside each resample and then re-evaluated on the same
# folds with their selected hyperparameters, so the reported metrics are not the
# optimistic "best grid cell" values.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(tidymodels)
  library(ranger)
  library(xgboost)
})

tidymodels_prefer(quiet = TRUE)

reg_metrics <- metric_set(rmse, mae, rsq)
clf_metrics <- metric_set(roc_auc, accuracy, sens, spec)

# the covariate set shared by both tasks
covariates <- c("log_pop", "log_inc", "log_hv", "log_dens",
                "owner_share", "mf_stock_share", "med_year_built")

# -----------------------------------------------------------------------------
# Helper: tune (when needed) then fit_resamples on the same folds
# -----------------------------------------------------------------------------
evaluate_workflow <- function(wf, folds, metrics, grid = NULL, select_by = "rmse") {
  if (is.null(grid)) {
    res <- fit_resamples(wf, resamples = folds, metrics = metrics,
                         control = control_resamples(save_pred = TRUE))
    return(list(cv = res, tuned = NULL, wf = wf))
  }
  tuned <- tune_grid(wf, resamples = folds, grid = grid, metrics = metrics,
                     control = control_grid(save_pred = TRUE))
  best  <- select_best(tuned, metric = select_by)
  wf_fin <- finalize_workflow(wf, best)
  cv <- fit_resamples(wf_fin, resamples = folds, metrics = metrics,
                      control = control_resamples(save_pred = TRUE))
  list(cv = cv, tuned = tuned, wf = wf_fin, best = best)
}

summarise_cv <- function(cv, model_name) {
  collect_metrics(cv) |>
    select(.metric, mean, std_err, n) |>
    mutate(model = model_name, .before = 1)
}

# -----------------------------------------------------------------------------
# Helper: partial dependence of a fitted workflow on one numeric predictor
# (average prediction with that predictor set to a grid, everything else left
# at its observed value)
# -----------------------------------------------------------------------------
pdp_numeric <- function(fit, data, var, n = 60, prob = FALSE) {
  grid <- seq(min(data[[var]]), max(data[[var]]), length.out = n)
  map_dfr(grid, function(g) {
    nd <- data
    nd[[var]] <- g
    if (prob) {
      p <- predict(fit, new_data = nd, type = "prob")$.pred_Yes
    } else {
      p <- predict(fit, new_data = nd)$.pred
    }
    tibble(!!var := g, pred = mean(p))
  })
}

# =============================================================================
# Task 1: regression -- permits per 1,000 residents
# =============================================================================
run_regression_models <- function(d, seed = 2025, v = 10, repeats = 3) {

  d <- d |>
    filter(!is.na(log1p_rate)) |>
    mutate(county = factor(county))

  set.seed(seed)
  folds <- vfold_cv(d, v = v, repeats = repeats)

  fml_pop  <- as.formula(paste("log1p_rate ~ log_pop"))
  fml_geo  <- as.formula(paste("log1p_rate ~ log_pop + county"))
  fml_char <- as.formula(paste("log1p_rate ~", paste(covariates, collapse = " + ")))
  fml_full <- as.formula(paste("log1p_rate ~",
                               paste(c(covariates, "county"), collapse = " + ")))

  rec_pop  <- recipe(fml_pop,  data = d)
  rec_geo  <- recipe(fml_geo,  data = d) |> step_dummy(all_nominal_predictors())
  rec_char <- recipe(fml_char, data = d)
  rec_full <- recipe(fml_full, data = d) |> step_dummy(all_nominal_predictors())

  lin <- function() linear_reg() |> set_engine("lm") |> set_mode("regression")

  wf_pop  <- workflow() |> add_recipe(rec_pop)  |> add_model(lin())
  wf_geo  <- workflow() |> add_recipe(rec_geo)  |> add_model(lin())
  wf_char <- workflow() |> add_recipe(rec_char) |> add_model(lin())
  wf_full <- workflow() |> add_recipe(rec_full) |> add_model(lin())

  # penalised models need standardised predictors and the full design
  rec_pen <- rec_full |> step_normalize(all_numeric_predictors())
  wf_ridge <- workflow() |> add_recipe(rec_pen) |>
    add_model(linear_reg(penalty = tune(), mixture = 0) |>
                set_engine("glmnet") |> set_mode("regression"))
  wf_lasso <- workflow() |> add_recipe(rec_pen) |>
    add_model(linear_reg(penalty = tune(), mixture = 1) |>
                set_engine("glmnet") |> set_mode("regression"))

  wf_rf  <- workflow() |> add_recipe(rec_full) |>
    add_model(rand_forest(mtry = tune(), min_n = tune(), trees = 1000) |>
                set_engine("ranger", importance = "permutation") |> set_mode("regression"))
  wf_xgb <- workflow() |> add_recipe(rec_full) |>
    add_model(boost_tree(trees = 1000, tree_depth = tune(), learn_rate = tune(),
                         min_n = tune(), loss_reduction = 0) |>
                set_engine("xgboost") |> set_mode("regression"))

  grid_pen <- grid_regular(penalty(), levels = 20)
  # mtry range depends on the number of columns created by the dummy step
  n_pred <- ncol(juice(prep(rec_full)))
  grid_rf  <- grid_regular(mtry(range = c(2, floor(n_pred / 2))),
                           min_n(range = c(2, 12)), levels = 3)
  grid_xgb <- grid_regular(tree_depth(range = c(2, 6)),
                           learn_rate(range = c(-2.3, -1)),
                           min_n(range = c(5, 15)), levels = 3)

  res <- list(
    "Size only (OLS)"                = evaluate_workflow(wf_pop,  folds, reg_metrics),
    "Size + county (OLS)"            = evaluate_workflow(wf_geo,  folds, reg_metrics),
    "Characteristics (OLS)"          = evaluate_workflow(wf_char, folds, reg_metrics),
    "Characteristics + county (OLS)" = evaluate_workflow(wf_full, folds, reg_metrics),
    "Ridge"                          = evaluate_workflow(wf_ridge, folds, reg_metrics, grid_pen),
    "Lasso"                          = evaluate_workflow(wf_lasso, folds, reg_metrics, grid_pen),
    "Random forest"                  = evaluate_workflow(wf_rf,   folds, reg_metrics, grid_rf),
    "XGBoost"                        = evaluate_workflow(wf_xgb,  folds, reg_metrics, grid_xgb)
  )

  metrics <- map2_dfr(res, names(res), ~ summarise_cv(.x$cv, .y))
  preds   <- map2_dfr(res, names(res),
                      ~ collect_predictions(.x$cv) |> mutate(model = .y))

  # ---- final fits on the full sample ----------------------------------------
  fit_full <- wf_full |> fit(data = d)
  fit_rf   <- finalize_workflow(wf_rf,  select_best(res[["Random forest"]]$tuned, metric = "rmse")) |> fit(data = d)
  fit_xgb  <- finalize_workflow(wf_xgb, select_best(res[["XGBoost"]]$tuned, metric = "rmse")) |> fit(data = d)

  ols_tidy   <- fit_full |> extract_fit_parsnip() |> tidy(conf.int = TRUE, conf.level = 0.95)
  ols_glance <- fit_full |> extract_fit_parsnip() |> glance()

  # ---- standardised coefficients: which characteristics matter most? --------
  fit_std <- workflow() |>
    add_recipe(rec_full |> step_normalize(all_numeric_predictors())) |>
    add_model(lin()) |>
    fit(data = d)
  std_coef <- fit_std |> extract_fit_parsnip() |> tidy() |>
    filter(term != "(Intercept)") |>
    mutate(term = str_remove(term, "^county")) |>
    arrange(desc(abs(estimate)))

  # ---- collinearity of the numeric design -----------------------------------
  vif_tab <- tryCatch({
    lm_fit <- lm(fml_char, data = d)
    car::vif(lm_fit) |> as.data.frame() |> tibble::rownames_to_column("term") |>
      rename(vif = 2) |> arrange(desc(vif))
  }, error = function(e) NULL)

  # ---- permutation importance from the forest -------------------------------
  rf_imp <- fit_rf |> extract_fit_engine() |> ranger::importance() |>
    enframe(name = "term", value = "importance") |>
    mutate(term = str_remove(term, "^county")) |>
    group_by(term) |> summarise(importance = sum(importance), .groups = "drop") |>
    arrange(desc(importance))

  # ---- partial dependence: the three characteristics of interest -----------
  pdp <- bind_rows(
    pdp_numeric(fit_xgb, d, "log_dens") |> mutate(variable = "log_dens"),
    pdp_numeric(fit_xgb, d, "log_hv")   |> mutate(variable = "log_hv"),
    pdp_numeric(fit_xgb, d, "log_inc")  |> mutate(variable = "log_inc")
  ) |> mutate(model = "XGBoost")

  # OLS line for the same variables, holding everything else at the mean town
  pdp_ols <- function(var) {
    m <- lm(fml_full, data = d)
    grid <- seq(min(d[[var]]), max(d[[var]]), length.out = 60)
    nd <- d[rep(1, length(grid)), ]
    nd[[var]] <- grid
    tibble(!!var := grid, pred = predict(m, newdata = nd), variable = var,
           model = "OLS")
  }
  pdp_ols_all <- bind_rows(pdp_ols("log_dens"), pdp_ols("log_hv"), pdp_ols("log_inc"))

  # ---- robustness: strictly positive rates, log(rate) ----------------------
  d_pos <- d |> filter(permits_per_1000 > 0) |> mutate(log_rate = log(permits_per_1000))
  set.seed(seed)
  folds_pos <- vfold_cv(d_pos, v = v, repeats = repeats)
  rec_pos <- recipe(log_rate ~ ., data = d_pos |> select(log_rate, all_of(c(covariates, "county")))) |>
    step_dummy(all_nominal_predictors())
  wf_pos <- workflow() |> add_recipe(rec_pos) |> add_model(lin())
  cv_pos <- fit_resamples(wf_pos, resamples = folds_pos, metrics = reg_metrics,
                          control = control_resamples(save_pred = TRUE))
  fit_pos <- wf_pos |> fit(data = d_pos)

  hp <- map2_dfr(res, names(res), function(x, nm) {
    if (is.null(x$best)) return(tibble())
    x$best |> mutate(model = nm, .before = 1)
  })

  list(
    metrics = metrics,
    predictions = preds,
    ols_tidy = ols_tidy,
    ols_glance = ols_glance,
    std_coef = std_coef,
    vif = vif_tab,
    rf_importance = rf_imp,
    pdp = pdp,
    pdp_ols = pdp_ols_all,
    ols_fit = fit_full,
    positive = list(metrics = summarise_cv(cv_pos, "OLS (positive rates only)"),
                    tidy = fit_pos |> extract_fit_parsnip() |> tidy(conf.int = TRUE),
                    n = nrow(d_pos)),
    hyperparameters = hp,
    n = nrow(d)
  )
}

# =============================================================================
# Task 2: classification -- did the town permit any multifamily housing?
# =============================================================================
run_classification_models <- function(d, seed = 2025, v = 10, repeats = 3) {

  d <- d |> filter(!is.na(mf_any)) |> mutate(county = factor(county))

  set.seed(seed)
  folds <- vfold_cv(d, v = v, repeats = repeats, strata = mf_any)

  fml <- as.formula(paste("mf_any ~", paste(c(covariates, "county"), collapse = " + ")))
  base_rec <- recipe(fml, data = d) |> step_dummy(all_nominal_predictors())

  # size-only baseline, for a like-for-like comparison with the full design
  rec_pop <- recipe(mf_any ~ log_pop, data = d)
  wf_pop <- workflow() |> add_recipe(rec_pop) |>
    add_model(logistic_reg() |> set_engine("glm") |> set_mode("classification"))

  wf_logit <- workflow() |> add_recipe(base_rec) |>
    add_model(logistic_reg() |> set_engine("glm") |> set_mode("classification"))
  wf_pen <- workflow() |> add_recipe(base_rec |> step_normalize(all_numeric_predictors())) |>
    add_model(logistic_reg(penalty = tune(), mixture = 0) |>
                set_engine("glmnet") |> set_mode("classification"))

  n_pred <- ncol(juice(prep(base_rec)))
  grid_rf  <- grid_regular(mtry(range = c(2, floor(n_pred / 2))),
                           min_n(range = c(2, 12)), levels = 3)

  wf_rf <- workflow() |> add_recipe(base_rec) |>
    add_model(rand_forest(mtry = tune(), min_n = tune(), trees = 1000) |>
                set_engine("ranger", importance = "permutation") |> set_mode("classification"))
  wf_xgb <- workflow() |> add_recipe(base_rec) |>
    add_model(boost_tree(trees = 1000, tree_depth = tune(), learn_rate = tune(),
                         min_n = tune(), loss_reduction = 0) |>
                set_engine("xgboost") |> set_mode("classification"))

  res <- list(
    "Size only (logistic)" = evaluate_workflow(wf_pop, folds, clf_metrics),
    "Logistic regression" = evaluate_workflow(wf_logit, folds, clf_metrics),
    "Penalised logistic"  = evaluate_workflow(wf_pen, folds, clf_metrics,
                                              grid_regular(penalty(), levels = 20),
                                              select_by = "roc_auc"),
    "Random forest"       = evaluate_workflow(wf_rf, folds, clf_metrics, grid_rf,
                                              select_by = "roc_auc"),
    "XGBoost"             = evaluate_workflow(wf_xgb, folds, clf_metrics,
                                              grid_regular(tree_depth(range = c(2, 6)),
                                                           learn_rate(range = c(-2.3, -1)),
                                                           min_n(range = c(5, 15)), levels = 3),
                                              select_by = "roc_auc")
  )

  metrics <- map2_dfr(res, names(res), ~ summarise_cv(.x$cv, .y))
  preds   <- map2_dfr(res, names(res),
                      ~ collect_predictions(.x$cv) |> mutate(model = .y))

  fit_logit <- wf_logit |> fit(data = d)
  fit_xgb   <- finalize_workflow(wf_xgb, select_best(res[["XGBoost"]]$tuned, metric = "roc_auc")) |> fit(data = d)

  logit_tidy <- fit_logit |> extract_fit_parsnip() |> tidy(conf.int = TRUE, exponentiate = TRUE)

  # ROC curves from the pooled out-of-fold predictions
  roc <- preds |> group_by(model) |>
    group_modify(~ roc_curve(.x, truth = mf_any, .pred_Yes)) |> ungroup()

  # calibration of the best model, by bin of predicted probability
  best_model <- metrics |> filter(.metric == "roc_auc") |> slice_max(mean) |> pull(model)
  best_preds <- preds |> filter(model == best_model)
  calib <- best_preds |>
    mutate(bin = ntile(.pred_Yes, 6)) |>
    group_by(bin) |>
    summarise(n = n(), pred = mean(.pred_Yes),
              obs = mean(mf_any == "Yes"), .groups = "drop")

  thr_grid <- seq(0.05, 0.95, by = 0.01)
  thr <- map_dfr(thr_grid, function(th) {
    lab <- factor(if_else(best_preds$.pred_Yes >= th, "Yes", "No"), levels = c("No", "Yes"))
    tp <- sum(lab == "Yes" & best_preds$mf_any == "Yes")
    fp <- sum(lab == "Yes" & best_preds$mf_any == "No")
    fn <- sum(lab == "No"  & best_preds$mf_any == "Yes")
    tn <- sum(lab == "No"  & best_preds$mf_any == "No")
    tibble(threshold = th, sens = tp / (tp + fn), spec = tn / (tn + fp),
           youden = sens + spec - 1)
  })
  best_thr <- thr |> slice_max(youden, n = 1) |> pull(threshold)

  list(metrics = metrics, predictions = preds, logit_tidy = logit_tidy,
       roc = roc, calibration = calib, threshold = thr,
       best_model = best_model, best_threshold = best_thr,
       pdp = pdp_numeric(fit_xgb, d, "log_pop", prob = TRUE) |> mutate(model = "XGBoost"),
       n = nrow(d), n_pos = sum(d$mf_any == "Yes"), folds = folds,
       hyperparameters = map2_dfr(res, names(res),
                                  function(x, nm) if (is.null(x$best)) tibble() else x$best |> mutate(model = nm)))
}

# =============================================================================
# Caching wrapper: reuse fitted results across renders unless forced
# =============================================================================
get_model_results <- function(data, which = c("regression", "classification"),
                              dir = "results", force = FALSE) {
  which <- match.arg(which)
  f <- file.path(dir, paste0(which, ".rds"))
  if (!force && file.exists(f)) {
    message("loading cached ", which, " results from ", f)
    return(readRDS(f))
  }
  res <- if (which == "regression") run_regression_models(data) else run_classification_models(data)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(res, f)
  res
}
