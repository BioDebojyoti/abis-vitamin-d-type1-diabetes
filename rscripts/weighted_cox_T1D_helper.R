clean_terms <- function(x) {
  if (is.na(x)) {
    return(NA_character_)
  } else {
    return(gsub("_", " ", gsub(
      "_cat", "", gsub(
        "_dichotomous",
        "",
        gsub("_during_pregnancy|_any_pregnancy", " (pregnancy)", x)
      )
    )))
  }
}

add_colon <- function(x) {
  # Replace underscores with spaces first
  x <- gsub("_", " ", x)
  
  # Split at camel-case boundary (lower → upper) OR last space
  x <- gsub("([a-z])([A-Z])", "\\1: \\2", x)   # add colon before uppercase
  x <-  gsub(
    "_cat", "", gsub(
      "_dichotomous",
      "",
      gsub("_during_pregnancy|_any_pregnancy", " (pregnancy)",x)))
  
  return(x)
}

format_numeric <- function(x) {
  ifelse(
    x < 0.01,
    formatC(x, format = "e", digits = 2),
    formatC(x, format = "f", digits = 2)
  )
}


OR_tidy_base_model <- function(model, model_name, predictor="conc_25OHD", 
                               df_cox_full = if (exists("df_cox_full", envir = .GlobalEnv)) 
                                 get("df_cox_full", envir = .GlobalEnv) 
                               else 
                                 NULL, firth_model = FALSE, 
                               predictor_increment = "+ 1") {
  
  
  quantile_sub = ifelse(grepl("quantile|quartile", predictor) ||
                          is.factor(df_cox_full[[predictor]]),
                        "",
                        predictor_increment)
  
  if(firth_model){
    result <- model %>%
      coxphf::tidy.coxphf(exponentiate = TRUE, conf.int = TRUE)
  } else {
    result <- model %>%
      broom::tidy(exponentiate = TRUE, conf.int = TRUE)
  }
  
  return(
    result %>%
      dplyr::mutate(p.value = ifelse(
        p.value < 0.001,
        formatC(p.value, format = "e", digits = 2),
        formatC(p.value, format = "f", digits = 3)
      )) %>%
      dplyr::mutate(Model = model_name, .before = 1) %>%
      dplyr::mutate(`OR (95% CI)` = paste0(
        round(estimate, 2),
        " (",
        round(conf.low, 2),
        "-",
        round(conf.high, 2),
        ")"
      )) %>%
      dplyr::rename(Term = term, `P-value` = p.value) %>%
      dplyr::mutate(Level = gsub(predictor,quantile_sub, Term)) %>%
      dplyr::mutate(Term = predictor) %>%
      dplyr::select(Model, Term, Level, `OR (95% CI)`, `P-value`) %>%
      dplyr::mutate(Term = ifelse(
        row_number() > 1 &
          Term == lag(Term),
        NA_character_,
        gsub("_", " ", gsub(
          "_cat", "", gsub(
            "_dichotomous",
            "",
            gsub(
              "Birth_weight_cat",
              "Birth weight category",
              gsub("Preterm_dichotomous", "", Term)
            )
          )
        ))
      ),
      Model = ifelse(row_number() > 1 &
                       Model == lag(Model), NA_character_, Model))
  )
}

OR_tidy_model_specific_covariates <- function(model,
                                              predictor = "conc_25OHD",
                                              model_name = "Base + Child covariates",
                                              covariates_as_pattern = "",
                                              covariate_levels_as_pattern = "", 
                                              df_cox_full = if (exists("df_cox_full", envir = .GlobalEnv)) 
                                                get("df_cox_full", envir = .GlobalEnv) 
                                              else 
                                                NULL,
                                              firth_model = FALSE, 
                                              predictor_increment = "+ 1") {
  
  quantile_sub = ifelse(grepl("quantile|quartile", predictor) ||
                          is.factor(df_cox_full[[predictor]]),
                        "",
                        predictor_increment)
  if(firth_model){
    result <- model %>%
      coxphf::tidy.coxphf(exponentiate = TRUE, conf.int = TRUE)
  } else {
    result <- model %>%
      broom::tidy(exponentiate = TRUE, conf.int = TRUE)
  }
  return(
    result %>% 
      dplyr::mutate(p.value = ifelse(
        p.value < 0.001,
        formatC(p.value, format = "e", digits = 2),
        formatC(p.value, format = "f", digits = 2)
      )) %>%
      dplyr::mutate(Model = model_name, .before = 1) %>%
      dplyr::mutate(`OR (95% CI)` = paste0(
        round(estimate, 2),
        " (",
        round(conf.low, 2),
        "-",
        round(conf.high, 2),
        ")"
      )) %>%
      dplyr::rename(Term = term, `P-value` = p.value) %>%
      dplyr::select(Model, Term, `OR (95% CI)`, `P-value`) %>%
      dplyr::mutate(
        Level = gsub(covariates_as_pattern, "", gsub(predictor, quantile_sub, Term)),
        .before = 3
      ) %>%
      dplyr::mutate(Term = gsub(covariate_levels_as_pattern, "", Term)) %>%
      dplyr::mutate(
        Term = ifelse(
          row_number() > 1 &
            Term == lag(Term),
          NA_character_,
          gsub("_", " ", gsub(
            "_cat", "", gsub(
              "_dichotomous",
              "",
              gsub(
                "Birth_weight_cat",
                "Birth weight category",
                gsub("Preterm_dichotomous", "", Term)
              )
            )
          ))
        ),
        Model = ifelse(row_number() > 1 &
                         Model == lag(Model), NA_character_, Model)
      )
  )
}


tidy_base_model <- function(model, model_name, predictor="conc_25OHD", 
                            df_cox_full = if (exists("df_cox_full", envir = .GlobalEnv)) 
                              get("df_cox_full", envir = .GlobalEnv) 
                            else 
                              NULL, 
                            predictor_increment = "+ 1") {
  
  
  quantile_sub = ifelse(grepl("quantile|quartile", predictor) ||
                          is.factor(df_cox_full[[predictor]]),
                        "",
                        predictor_increment)
  
  s <- summary(model)
  
  result <- data.frame(
    term = rownames(s$coefficients),
    HR = format_numeric(s$coefficients[, "exp(coef)"]),
    lower_CI = format_numeric(s$conf.int[, "lower .95"]),
    upper_CI = format_numeric(s$conf.int[, "upper .95"]),
    p.value = format_numeric(s$coefficients[, "Pr(>|z|)"]),
    row.names = NULL
  )
  
  return(
    result %>%
      dplyr::mutate(Model = model_name, .before = 1) %>%
      dplyr::mutate(`HR (95% CI)` = paste0(HR, " (", lower_CI, "-", upper_CI, ")")) %>%
      dplyr::rename(Term = term, `P-value` = p.value) %>%
      dplyr::mutate(Level = gsub(predictor,quantile_sub, Term)) %>%
      dplyr::mutate(Term = predictor) %>%
      dplyr::select(Model, Term, Level, `HR (95% CI)`, `P-value`) %>%
      dplyr::mutate(
        Term = ifelse(
          row_number() > 1 &
            Term == lag(Term),
          NA_character_,
          gsub("_", " ", gsub(
            "_cat", "", gsub(
              "_dichotomous",
              "",
              gsub(
                "Birth_weight_cat",
                "Birth weight category",
                gsub("Preterm_dichotomous", "", Term)
              )
            )
          ))
        ),
        Model = ifelse(row_number() > 1 &
                         Model == lag(Model), NA_character_, Model)
      )
  )
}

tidy_model_specific_covariates <- function(model,
                                           predictor = "conc_25OHD",
                                           model_name = "Base + Child covariates",
                                           covariates_as_pattern = "",
                                           covariate_levels_as_pattern = "", 
                                           df_cox_full = if (exists("df_cox_full", envir = .GlobalEnv)) 
                                             get("df_cox_full", envir = .GlobalEnv) 
                                           else 
                                             NULL, 
                                           predictor_increment = "+ 1") {
  
  quantile_sub = ifelse(grepl("quantile|quartile", predictor) || is.factor(df_cox_full[[predictor]]), "", predictor_increment)
  
  s <- summary(model)
  
  result <- data.frame(
    term = rownames(s$coefficients),
    HR = format_numeric(s$coefficients[, "exp(coef)"]),
    lower_CI = format_numeric(s$conf.int[, "lower .95"]),
    upper_CI = format_numeric(s$conf.int[, "upper .95"]),
    p.value = format_numeric(s$coefficients[, "Pr(>|z|)"]),
    row.names = NULL
  )
  
  return(
    result %>% 
      dplyr::mutate(Model = model_name, .before = 1) %>%
      dplyr::mutate(`HR (95% CI)` = paste0(HR, " (", lower_CI, "-", upper_CI, ")")) %>%
      dplyr::rename(Term = term, `P-value` = p.value) %>%
      dplyr::select(Model, Term, `HR (95% CI)`, `P-value`) %>%
      dplyr::mutate(Level = gsub(covariates_as_pattern, "", gsub(predictor, quantile_sub, Term)),
                    .before = 3) %>%
      dplyr::mutate(Term = gsub(covariate_levels_as_pattern, "", Term)) %>%
      dplyr::mutate(
        Term = ifelse(
          row_number() > 1 &
            Term == lag(Term),
          NA_character_,
          gsub("_", " ", gsub(
            "_cat", "", gsub(
              "_dichotomous",
              "",
              gsub(
                "Birth_weight_cat",
                "Birth weight category",
                gsub("Preterm_dichotomous", "", Term)
              )
            )
          ))
        ),
        Model = ifelse(row_number() > 1 &
                         Model == lag(Model), NA_character_, Model)
      )
  )
}

run_stratified_cox_models <- function(data, 
                                      predictor, 
                                      strata_vars, 
                                      covariate_list = NULL,
                                      covariates_as_pattern = "",
                                      covariate_levels_as_pattern = "", 
                                      predictor_increment = "+ 1") {
  
  
  quantile_sub = ifelse(grepl("quantile|quartile", predictor) ||
                          is.factor(data[[predictor]]),
                        "",
                        predictor_increment )
  results <- list()
  models <- list()
  
  for (strata_var in strata_vars) {
    # Construct formula
    if(is.null(covariate_list)){
      fmla <- as.formula(paste0("Surv(exit_age, event = event) ~ ",predictor," + strata(", strata_var, ")"))
    } else {
      fmla <- as.formula(paste0("Surv(exit_age, event = event) ~ ",predictor," + strata(", strata_var, ") + ", paste0(covariate_list, collapse = " + ")))
    }
    
    # print(fmla)
    
    # Fit Cox model
    cox_model <- survival::coxph(fmla, data = data)
    
    # Tidy and format
    tidy_result <- broom::tidy(cox_model, exponentiate = TRUE, conf.int = TRUE) %>%
      dplyr::rename(Term = term,OR = estimate,`P-value` = p.value) %>%
      dplyr::mutate(`HR (95% CI)` = paste0(round(OR, 2), " (", round(conf.low, 2), "-", round(conf.high, 2), ")")) %>%
      dplyr::mutate(Strata = strata_var) %>%
      dplyr::select(Term, Strata, `HR (95% CI)`, `P-value`) %>%
      dplyr::mutate(across(where(is.numeric), ~ ifelse(
        .x < 0.001,
        formatC(.x, format = "e", digits = 2),
        formatC(.x, format = "f", digits = 2)
      ))) %>%
      dplyr::mutate(Level = gsub(covariates_as_pattern, "", gsub(predictor, quantile_sub, Term)), .before = 3) %>%
      dplyr::mutate(Term = gsub(covariate_levels_as_pattern, "", Term)) %>%
      dplyr::mutate(Term = ifelse(row_number() > 1 & Term == lag(Term), NA_character_, gsub("_", " ", gsub("_cat", "", gsub("_dichotomous",
                                                                                                                            "", gsub("Birth_weight_cat",  "Birth weight category", gsub("Preterm_dichotomous", "Preterm", Term)))))
      )) %>%
      dplyr::mutate(Strata = ifelse(
        row_number() > 1 & Strata == lag(Strata),
        NA_character_,
        gsub("_", " ", gsub(
          "_cat", "", gsub(
            "_dichotomous",
            "",
            gsub("_during_pregnancy|_any_pregnancy", " (pregnancy)", Strata)
          )
        ))
      ))
    
    
    
    # Store in list
    results[[strata_var]] <- tidy_result
    models[[strata_var]] <- cox_model
  }
  
  # Combine all results
  final_result <- bind_rows(results) 
  
  return(list(final_result,models))
}

first_letter_capital <- function(s){
  if(is.na(s)){
    return(NA_character_)
  } else {
    s_list <- stringr::str_split(s, "")[[1]]
    return(
      paste0(toupper(s_list[1]),paste0(s_list[c(2:length(s_list))], collapse = ""))
    )
  }
}

CI_significance <- function(s){
  
  m <- stringr::str_match(s, "\\(([^-]+)-([^\\)]+)\\)")
  low  <- as.numeric(m[,2])
  high <- as.numeric(m[,3])
  out <- (low < 1 & high < 1) | (low > 1 & high > 1)
  out
}

table2flex <- function(df){
  return(
    df %>% 
      dplyr::mutate(
        Term = sapply(Term, first_letter_capital)
      ) %>%
      flextable::flextable() %>% 
      flextable::autofit() %>% 
      flextable::set_table_properties(width = 1, layout = "autofit") %>% 
      flextable::color(
        i = ~ as.numeric(`P-value`) < 0.05,
        j = c("Term", "Level", "P-value"),
        color = "red")  %>% 
      flextable::color(
        i = ~ grepl("25OHD|DBP",Term),
        j = c("Term"),
        color = "blue") 
  )
}

make_comp_table <- function(m1, m2, label) {
  an <- anova(m1, m2, test = "LRT") # likelihood ratio test
  
  data.frame(
    Comparison = label,
    logLik_m1 = as.numeric(logLik(m1)),
    logLik_m2 = as.numeric(logLik(m2)),
    AIC_m1 = AIC(m1),
    AIC_m2 = AIC(m2),
    # BIC_m1 = BIC(m1),
    # BIC_m2 = BIC(m2),
    Chisq = an$Chisq[2],
    df = an$Df[2],
    P.value = an$"Pr(>|Chi|)"[2]
  ) %>% 
    dplyr::rename(
      `P-value` = P.value,
      `AIC \n(no-interaction model)` = AIC_m1,      
      `AIC \n(interaction model)` = AIC_m2,
      `logLik \n(no-interaction model)` = logLik_m1,      
      `logLik \n(interaction model)` = logLik_m2,
    ) %>%
    dplyr::mutate(across(where(is.numeric), ~ ifelse(
      abs(.x) < 0.001,
      formatC(.x, format = "e", digits = 2),
      formatC(.x, format = "f", digits = 2)
    )))
  
}

tidy_ggsurvplot <- function(ggsurvplot_obj, facet_grid_formula = NULL){
  
  if(!is.null(facet_grid_formula)){
    ggsurvplot_obj$plot <- ggsurvplot_obj$plot + facet_grid(as.formula(facet_grid_formula))
  }  
  
  ggsurvplot_obj$table$data$strata <- sapply(ggsurvplot_obj$table$data$strata, clean_terms)
  ggsurvplot_obj$table <- ggsurvplot_obj$table + 
    scale_y_discrete(labels = function(x) sapply(x, clean_terms)) +
    theme(legend.position = "none")
  ggsurvplot_obj$table <- ggsurvplot_obj$table + theme(axis.title.y = element_blank())
  
  ggsurvplot_obj$plot <- ggsurvplot_obj$plot + theme(legend.position = "none")
  
  ggsurvplot_obj
}

# Case–cohort weighted Cox (using time-to-event):
#   f <- n_subcohort / n_full_cohort
sampling_fractions <- function(full_cohort, sub_cohort){
  f <- nrow(sub_cohort)/ nrow(full_cohort)
  f_cases <- nrow(sub_cohort %>% dplyr::filter(event==1))/ nrow(full_cohort %>% dplyr::filter(event==1))
  
  print("Sampling fractions:")
  print("Sub cohort:")
  print(f)
  print("Subcohort cases:")
  print(f_cases)
  
  return(list(f,f_cases))  
}

spline_cox_analysis <- function(data, 
                                time, event, 
                                exposure, 
                                covariates = NULL, 
                                df_spline = 3) {
  # -----------------------------
  # 1. Build formula
  # -----------------------------
  covariate_part <- if (!is.null(covariates)) paste(covariates, collapse = " + ") else ""
  formula_text <- paste0("Surv(", time, ", ", event, ") ~ splines::ns(", exposure, ", df = ", df_spline, ")",
                         if (covariate_part != "") paste0(" + ", covariate_part) else "")
  formula <- as.formula(formula_text)
  
  print(formula)
  # -----------------------------
  # 2. Fit spline Cox model
  # -----------------------------
  cox_model <- coxph(formula, data = data)
  
  # -----------------------------
  # 3. Create prediction grid across exposure range
  # -----------------------------
  # grid <- data.frame(exposure_value = seq(min(data[[exposure]], na.rm = TRUE),
  #                                         max(data[[exposure]], na.rm = TRUE),
  #                                         length.out = 100))
  
  # Use 1st–99th percentile instead of absolute min/max to avoid extreme tails
  obs_range <- quantile(data[[exposure]], probs = c(0.01, 0.99), na.rm = TRUE)
  grid <- data.frame(exposure_value = seq(obs_range[1], obs_range[2], length.out = 100))
  
  if (!is.null(covariates)) {
    for (cov in covariates) {
      col <- data[[cov]]
      if (is.numeric(col)) {
        grid[[cov]] <- median(col, na.rm = TRUE)
      } else {
        grid[[cov]] <- levels(data[[cov]])[1]
      }
    }
  }
  names(grid)[1] <- exposure
  print(str(grid))
  
  # -----------------------------
  # 4. Predict log-HR and CI
  # -----------------------------
  lp <- predict(cox_model, newdata = grid, type = "lp", se.fit = TRUE)
  grid$HR <- exp(lp$fit)
  grid$HR_low <- exp(lp$fit - 1.96 * lp$se.fit)
  grid$HR_high <- exp(lp$fit + 1.96 * lp$se.fit)
  
  # -----------------------------
  # 5. Quartile-based table (nearest grid points)
  # -----------------------------
  quartiles <- quantile(data[[exposure]], probs = c(0.25, 0.5, 0.75, 1), na.rm = TRUE)
  attr(quartiles, "names") <- NULL
  
  table_quartiles <- sapply(quartiles, function(q) {
    idx <- which.min(abs(grid[[exposure]] - q))
    grid[idx, c(exposure, "HR", "HR_low", "HR_high")]
  }, simplify = FALSE) %>%
    bind_rows() %>%
    dplyr::mutate(Quartile = as.character(quartiles)) %>%
    dplyr::select(Quartile, everything())
  
  # -----------------------------
  # 6. Publication-ready plot
  # -----------------------------
  plot <- ggplot(grid, aes_string(x = exposure, y = "HR")) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = HR_low, ymax = HR_high), alpha = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(x = paste0(exposure, " (continuous)"),
         y = "Adjusted Hazard Ratio",
         title = paste0("Adjusted HR across ", exposure, " using splines")) +
    theme_minimal(base_size = 14)
  
  return(list(cox_model = cox_model,
              prediction_grid = grid,
              table_quartiles = table_quartiles,
              plot = plot))
}
