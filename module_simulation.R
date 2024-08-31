MODULE_SIMULATION <- modules::module({
  
  import("stats")
  import("dplyr")
  import("tidyr")
  
  MODULE_BAYES_FACTOR <- use("module_bayes_factor.R")$MODULE_BAYES_FACTOR
  
  N_SUBJ        = c(20,50,100)   # number of subjects
  N_INTERGROUP  = c (2, 5, 10)   # number of ingroup stimuli
  N_OUTGROUP    = c(2, 5, 10)  # number of outgroup stimuli
  RHO           = c(-0.7, -0.5, 0, 0.7, 0.5)   # correlation between intercept and slope
  N_REP = 20
  MODEL_TYPES = c("pois", "logit","guassian")
  #MODEL_TYPES = c("logit")
  
  CONTROL_GLM <- list(
    
    iter_EM = 20,
    
    tol1 = 1e-03,
    tol2 = 1e-04,
    tol3 = 1e-08,
    nAGQ=20,
    iter_qN_outer=30,
    max_coef_value = 5000
  )
  
  # The current version of GLMMadaptive, don't supper normal distribution  
  temp_fun <- gaussian()
  temp_fun$variance <-  NULL # this it no needed
  temp_fun$family <- "gaussian2" # if the name is not change, GLMMadaptive will stop and subgest to use lmer
  temp_fun$log_dens <- function(y, eta, mu_fun, phis, eta_zi) {
    
    
    scale <- pmax(exp(phis), .Machine$double.eps)
    
    m_y <- mu_fun(eta)
    
    out <- dnorm(y,m_y,scale,log=TRUE)
    
    
    attr(out, "mu_y") <- m_y
    
    return(out)
    
  }
  
  BETA_0_NORMAL <- 800
  BETA_1_NORMAL <- 50
  
  
  BETA_0_POIS <- 5
  BETA_1_POIS <- 4
  
  BETA_0_LOGIT <- 0.5
  BETA_1_LOGIT <- 1
  
  
  
  my_sim_data <- function(
    n_subj     = 100,   # number of subjects
    n_ingroup  =  25,   # number of ingroup stimuli
    n_outgroup =  25,   # number of outgroup stimuli
    omega_0    =  80,   # by-item random intercept standard deviation
    tau_0      = 100,   # by-subject random intercept standard deviation
    tau_1      =  40,   # by-subject random slope standard deviation
    rho        = 0.2,   # correlation between intercept and slope
    sigma      = 200,   # residual standard deviation
    model_type = "guassian", # type of model: gaussian, poisson, or logistic
    exclude_items = TRUE # whether to exclude item-specific effects
  ) {
    
    # Create a data frame for items with their categories and random intercepts
    items <- data.frame(
      item_id = seq_len(n_ingroup + n_outgroup),
      category = rep(c("ingroup", "outgroup"), c(n_ingroup, n_outgroup)),
      X_i = rep(c(-0.5, 0.5), c(n_ingroup, n_outgroup)),
      O_0i = rnorm(n = n_ingroup + n_outgroup, mean = 0, sd = omega_0)
    )
    
    # Create the variance-covariance matrix for the random effects
    cov_mx <- matrix(
      c(tau_0^2,             rho * tau_0 * tau_1,
        rho * tau_0 * tau_1, tau_1^2),
      nrow = 2, byrow = TRUE
    )
    
    # Generate random effects for subjects using multivariate normal distribution
    subjects <- data.frame(subj_id = seq_len(n_subj),
                           MASS::mvrnorm(n = n_subj,
                                         mu = c(T_0s = 0, T_1s = 0),
                                         Sigma = cov_mx))
    
    # Generate the model data based on the specified model type
    if (model_type == "guassian") {
      
      # Gaussian model
      m_data <- crossing(subjects, items) %>%
        mutate(e_si = rnorm(nrow(.), mean = 0, sd = sigma),
               RT = BETA_0_NORMAL + T_0s + ((!exclude_items) * O_0i) + 
                 (BETA_1_NORMAL + T_1s) * X_i + e_si) %>%
        select(subj_id, item_id, category, X_i, RT)
      
    } else if (model_type == "pois") {
      
      # Poisson model
      m_data <- crossing(subjects, items) %>%
        mutate(RT = rpois(nrow(.), 
                          exp((BETA_0_POIS + T_0s + ((!exclude_items) * O_0i) + 
                                 (BETA_1_POIS + T_1s) * X_i) / 30))) %>%
        select(subj_id, item_id, category, X_i, RT)
      
    } else {
      
      # Logistic model
      m_data <- crossing(subjects, items) %>%
        mutate(RT = rbinom(nrow(.), 1, 
                           plogis(BETA_0_LOGIT + T_0s + ((!exclude_items) * O_0i) + 
                                    (BETA_1_LOGIT + T_1s) * X_i))) %>%
        select(subj_id, item_id, category, X_i, RT)
    }
    
    # Return the generated data
    return(m_data)
  }
  
  
  
  # Define a function 'get_formula' that takes three arguments: 'response', 'predictors', and 'slope_var'
  get_formula <- function(response, predictors, slope_var) {
    
    # Filter predictors containing the substring "rand_" and remove "rand_" from their names
    predictors_rand <- predictors[predictors |> 
                                    stringr::str_detect("rand_")] |> stringr::str_remove("rand_")
    
    # Filter predictors containing the substring "fixed_" and remove "fixed_" from their names
    predictors_fixed <-  predictors[predictors |> 
                                      stringr::str_detect("fixed_")] |> stringr::str_remove("fixed_")
    
    # Add "1" to the list of fixed predictors and ensure uniqueness
    predictors_fixed <- c("1", predictors_fixed) |> unique()
    
    # Create the base formula for fixed effects by concatenating the fixed predictors with "+"
    base_formula <- paste0(predictors_fixed, collapse = "+")
    
    # Combine the response variable with the base formula for fixed effects
    base_formula <- paste0(response, " ~ ", base_formula)
    
    # Check if 'slope_var' is NA; if true, create random terms without slope variable
    if(is.na(slope_var))  {
      random_terms <- paste0("1 | ", paste0(predictors_rand, collapse = "+"))
    } else {
      # If 'slope_var' is provided, include it in the random terms
      random_terms <- paste0("1 + ", slope_var ," | ", paste0(predictors_rand, collapse = "+"))
    }
    
    # Create the random part of the formula
    ramdom_part <- paste0(" ~ ", random_terms)
    
    # Create a list containing the random and fixed parts of the formula as formulas
    formula_list <- list(
      "ramdom_part" = ramdom_part  |> as.formula(),
      "fixed_part"  = base_formula |> as.formula()
    )
    
    # Return the list of formulas
    return(formula_list)
  }
  
  
  # Define a function to fit a statistical model
  fit_model <- function(simul_data, m_formula, random_varaibles, model_type) {
    
    # Suppress warnings and messages to keep the output clean
    suppressWarnings({
      suppressMessages({ 
        # Try to fit the model
        fitted_model <- ({
          # If model type is "gaussian"
          if(model_type == "guassian") {
            # Fit a linear model as a temporary step
            model_temp <- lm(m_formula$fixed_part, data=simul_data)
            
            # Initialize values for the mixed model
            initial_values <- list(
              "betas" = coef(model_temp),
              "D" = diag(length(random_varaibles)+1),
              "phis" = 0
            )
            
            # Fit a mixed effects model with a Gaussian family
            model <- GLMMadaptive::mixed_model(
              fixed = m_formula$fixed_part, 
              random = m_formula$ramdom_part, 
              data = simul_data,
              n_phis = 1,
              initial_values = initial_values,
              family = temp_fun, 
              control = CONTROL_GLM
            )
          } else if(model_type == "pois") {
            # Fit a mixed effects model with a Poisson family
            model <- GLMMadaptive::mixed_model(
              fixed = m_formula$fixed_part, 
              random = m_formula$ramdom_part,
              data = simul_data,
              family = poisson(),
              control = CONTROL_GLM
            )
          } else {
            # Fit a mixed effects model with a binomial family
            model <- GLMMadaptive::mixed_model(
              fixed = m_formula$fixed_part, 
              random = m_formula$ramdom_part,
              data = simul_data,
              family = binomial(),
              control = CONTROL_GLM
            )
          }
          
          model # Return the fitted model
        })
      })
    })
    
    # browser()
    
    # Return NULL if model fitting resulted in an error
    if ("try-error" %in% class(fitted_model)) { return(NULL) }
    
    # Extract fixed effects from the fitted model and convert to list
    m_result <- fitted_model |> GLMMadaptive::fixef() |> as.list()
    
    # Extract the Hessian matrix from the fitted model
    hessian <- fitted_model$Hessian
    
    # Calculate determinants and log determinants of the Hessian matrix
    m_result$det <- hessian |> Matrix::det()
    m_result$det_abs <- m_result$det |> abs()
    m_result$log_det <- m_result$det_abs |> log()
    
    # Calculate the variance of the fitted model
    m_result$var <- as.numeric(fitted_model$phis)**2
    
    if(length(m_result$var) != 1) {
      m_result$var <- NA_real_
    }
    
    # Calculate AIC and BIC for model selection criteria
    m_result$aic <- AIC(fitted_model)
    m_result$bic <- BIC(fitted_model)
    
    # Calculate a custom BIC using the Diogenes module
    m_result$bic_dio <- MODULE_BAYES_FACTOR$get_diogenes_bic(fitted_model)
    
    # Initialize AICc and calculate if possible
    m_result$aicc <- NA_real_
    temp_aicc <- try({ MODULE_BAYES_FACTOR$AICc(fitted_model) })
    if (!("try-error" %in% class(temp_aicc))) { m_result$aicc <- temp_aicc }
    
    # Initialize and calculate custom BIC without parameters if possible
    m_result$diogenes_bic_without_params <- NA_real_
    temp_diogenes_bic_without_params <- try({ MODULE_BAYES_FACTOR$get_diogenes_bic_without_params(fitted_model) })
    if (!("try-error" %in% class(temp_diogenes_bic_without_params))) {
      m_result$diogenes_bic_without_params <- temp_diogenes_bic_without_params
    }
    
    # Initialize and calculate CAIF if possible
    m_result$caif <- NA_real_
    temp_caif <- try({ MODULE_BAYES_FACTOR$get_caif(fitted_model) })
    if (!("try-error" %in% class(temp_caif))) { m_result$caif <- temp_caif }
    
    # Initialize and calculate ICOMP if possible
    m_result$icomp <- NA_real_
    temp_icomp <- try({ MODULE_BAYES_FACTOR$get_icomp(fitted_model) })
    if (!("try-error" %in% class(temp_icomp))) { m_result$icomp <- temp_icomp }
    
    # Return the result list
    return(m_result)
  }
  
  
  
  # Define a function to run a single model on the provided data
  run_single_model <- function(data, m_vars_data) {
    
    # Initialize an empty list to store results
    result_list <- list()
    
    # Simulate data based on the provided parameters
    simul_data <- my_sim_data(
      n_subj = data$n_subj,
      n_ingroup = data$n_ingroup,
      n_outgroup = data$n_outgroup,
      omega_0 = data$omega_0,
      tau_0 = data$tau_0,
      tau_1 = data$tau_1,
      rho = data$rho,
      sigma = data$sigma,
      model_type = data$model_type,
      exclude_items = TRUE
    )
    
    # Select predictor variables, excluding the random slope
    m_predictors <- select(m_vars_data, -random_slope) |> colnames()
    model_response_var <- "RT"  # Define the response variable
    
    # Loop through each row in m_vars_data
    for(idx in 1:nrow(m_vars_data)) {
      
      try({
        # Get the current random slope variable
        m_slope_var <- pull(m_vars_data, random_slope)[idx]
        
        # Get the current predictors based on the logical selection
        temp_index <- select(m_vars_data, -random_slope)[idx, ] |> as.logical()
        m_current_predictors <- m_predictors[temp_index]
        
        # Generate the current formula for the model
        curren_formula <- get_formula(model_response_var, m_current_predictors, m_slope_var)
        
        # Count the number of random variables
        m_random_variables <- m_current_predictors[m_current_predictors |> 
                                                     stringr::str_detect("rand_")] |> length()
        
        # Fit the model using the current formula and random variables
        fitted_model <- fit_model(simul_data, curren_formula, m_random_variables, data$model_type)
        
        # If the model fitting resulted in NULL, skip to the next iteration
        if(is.null(fitted_model)) { next }
        
        # Append the fitted model results to the current list
        current_list <- base::append(fitted_model, m_vars_data[idx, ] |> as.list())
        current_list <- base::append(current_list, data)
        
        # Add the current list to the result list
        result_list[[length(result_list) + 1]] <- current_list
      })
    }
    
    # browser()
    
    # Create a directory to save the results
    base_dir_name <- "23_05_26"
    dir.create(base_dir_name, showWarnings = FALSE, recursive = TRUE)
    
    # Save the result list as an RDS file
    saveRDS(result_list, file.path(base_dir_name, stringr::str_glue("model_{data$grid_idx}_{data$rep_id}.Rds")))
    
    # Return the result list
    return(result_list)
  }
  
  # Define a function to return a grid of model variables
  return_grid <- function() {
    
    model_response_var <- "RT"  # Define the response variable
    
    # Define fixed predictors, excluding the response variable
    m_predictors_fixed <- c("X_i") |> setdiff(model_response_var)
    
    # Define random predictors, excluding the response variable
    m_predictors_rand <- c("subj_id") |> setdiff(model_response_var)
    
    # Create a grid of all combinations (Cartesian product) of fixed predictors
    m_vars_data_fixed <- do.call(data.table::CJ, replicate(length(m_predictors_fixed), 0:1, FALSE))
    colnames(m_vars_data_fixed) <- paste0("fixed_", m_predictors_fixed)
    
    # Create a grid of all combinations (Cartesian product) of random predictors
    m_vars_data_rand <- do.call(data.table::CJ, replicate(length(m_predictors_rand), 0:1, FALSE))
    colnames(m_vars_data_rand) <- paste0("rand_", m_predictors_rand)
    
    # Cross join the fixed and random predictors grids
    m_vars_data <- cross_join(m_vars_data_fixed, m_vars_data_rand) |> 
      filter(rand_subj_id != 0)  # Filter out rows where random subject ID is zero
    
    # Bind rows to include both random slopes and no random slopes
    m_vars_data <- bind_rows(
      m_vars_data |> mutate(random_slope = "X_i"),  # Add a column for random slope with "X_i"
      m_vars_data |> mutate(random_slope = NA_character_)  # Add a column for random slope with NA
    )
    
    # Return the resulting grid of variables
    return(m_vars_data)
  }
  
  rep_df <- function(df, n_reps, id_var) {
    
    # browser()
    result_list <- list()
    
    for(i in 1:n_reps) {
      result_list[[i]] <- df 
      result_list[[i]][[id_var]] <- i
    }
    
    result_list <- result_list |> bind_rows()
    
    return(result_list)
    
  }
  
  # Define a function to generate a simulation grid with various parameters
  get_simulation_grid <- function(
    n_subj     = N_SUBJ,         # Number of subjects
    n_ingroup  = N_INTERGROUP,   # Number of ingroup stimuli
    n_outgroup = N_OUTGROUP,     # Number of outgroup stimuli
    model_type = MODEL_TYPES,    # Distributions of the response variable
    omega_0    = 80,             # By-item random intercept standard deviation
    tau_0      = 100,            # By-subject random intercept standard deviation
    tau_1      = 40,             # By-subject random slope standard deviation
    rho        = RHO,            # Correlation between intercept and slope
    sigma      = 200,            # Residual standard deviation
    reps       = N_REP           # Number of repetitions
  ) {
    
    # If 'reps' is a single number, create a sequence from 1 to that number
    # if(length(reps) < 2) {
    #   reps <- 1:reps
    # }
    
    # Get the grid of variables
    m_vars_data <- return_grid()
    
    # browser()
    
    # Create a data frame with all combinations of the specified parameters
    m_grid <- expand.grid(
      n_subj = n_subj,
      n_ingroup = n_ingroup,
      n_outgroup = n_outgroup,
      omega_0 = omega_0,
      tau_0 = tau_0,
      tau_1 = tau_1,
      rho = rho,
      sigma = sigma,
      model_type = model_type
    ) |> 
      data.frame() |> 
      mutate(model_type = model_type |> as.character()) |>  # Convert model_type to character
      rep_df(reps, id_var="rep_id") |>                        # number of repetitions
      mutate(grid_idx = 1:n()) |>                           # Add a column for grid index
      sample_n(n()) |> 
      purrr::transpose()                                    # Transpose the data frame
    
    # Return a list containing the grid and the variable data
    return(list(
      "m_grid" = m_grid,
      "m_vars_data" = m_vars_data
    ))
  }
  
  
})