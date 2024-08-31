MODULE_BAYES_FACTOR <- modules::module({
  
  import("stats")
  import("dplyr")
  import("tidyr")
  import("purrr")
  
  
  # Define a function to compute the Hessian matrix for a linear model
  get_hessiano_lm <- function(lm_model) {
    
    # Get the covariance matrix of the linear model coefficients
    cov_matrix <- lm_model |> stats::vcov()
    
    # Calculate the inverse of the covariance matrix
    inv_matrix <- Matrix::solve(cov_matrix)
    
    # Get the summary of the linear model to extract sigma (residual standard error)
    temp_summary <- summary(lm_model)
    
    # Compute the Hessian matrix using the formula 2 * (sigma^2) * inverse(covariance matrix)
    result <- 2 * (temp_summary$sigma**2) * inv_matrix
    
    # Return the resulting Hessian matrix
    return(result)
  }
  
  
  # Define a function to compute the Diogenes BIC for a given model
  get_diogenes_bic <- function(model_1) {
    
    # Check if the model is a mixed model ("MixMod" class)
    if(("MixMod" %in% class(model_1))) {
      # Extract the Hessian matrix directly from the mixed model
      hessian_1 <- model_1$Hessian
    } else {
      # Compute the Hessian matrix for a linear model using the custom function
      hessian_1 <- get_hessiano_lm(model_1)
    }
    
    # Compute the determinant of the Hessian matrix and take the absolute value
    hessian_1_det <- hessian_1 |> Matrix::det() |> abs()
    
    # Compute the Diogenes BIC by adding the log of the determinant to the standard BIC
    bic_genes <- BIC(model_1) + log(hessian_1_det)
    
    # Return the computed Diogenes BIC
    return(bic_genes)
  }
  
  
  
  # Define a function to compute the Diogenes BIC without parameters for a given model
  get_diogenes_bic_without_params <- function(model_1) {
    
    # Check if the model is a mixed model ("MixMod" class)
    if(("MixMod" %in% class(model_1))) {
      # Extract the Hessian matrix directly from the mixed model
      hessian_1 <- model_1$Hessian
    } else {
      # Compute the Hessian matrix for a linear model using the custom function
      hessian_1 <- get_hessiano_lm(model_1)
    }
    
    # Get the number of parameters (columns in the Hessian matrix)
    n_params <- hessian_1 |> ncol()
    
    # Get the number of observations (rows in the model's data)
    n_obs <- model_1$data |> nrow() 
    
    # Compute the determinant of the Hessian matrix and take the absolute value
    hessian_1_det <- hessian_1 |> Matrix::det() |> abs()
    
    # Compute the Diogenes BIC without parameters
    bic_genes <- BIC(model_1) + log(hessian_1_det) - (n_params * log(n_obs))
    
    # Return the computed Diogenes BIC without parameters
    return(bic_genes)
  }
  
  # Define a function to compute the CAIF (Consistent AIC with Fisher Information) for a given model
  get_caif <- function(model_1) {
    
    # Check if the model is a mixed model ("MixMod" class)
    if(("MixMod" %in% class(model_1))) {
      # Extract the Hessian matrix directly from the mixed model
      hessian_1 <- model_1$Hessian
    } else {
      # Compute the Hessian matrix for a linear model using the custom function
      hessian_1 <- get_hessiano_lm(model_1)
    }
    
    # Get the number of parameters (columns in the Hessian matrix)
    n_params <- hessian_1 |> ncol()
    
    # Get the number of observations (rows in the model's data)
    n_obs <- model_1$data |> nrow() 
    
    # Compute the determinant of the Hessian matrix and take the absolute value
    hessian_1_det <- hessian_1 |> Matrix::det() |> abs()
    
    # Compute the CAIF by adding the standard BIC, twice the number of parameters, and the log of the Hessian determinant
    caif_value <- BIC(model_1) + (2 * n_params) + log(hessian_1_det)
    
    # Return the computed CAIF value
    return(caif_value)
  }
  
  # Define a function to compute the ICOMP (Information Complexity Criterion) for a given model
  get_icomp <- function(model_1) {
    
    # Check if the model is a mixed model ("MixMod" class)
    if(("MixMod" %in% class(model_1))) {
      # Extract the Hessian matrix directly from the mixed model
      hessian_1 <- model_1$Hessian
    } else {
      # Compute the Hessian matrix for a linear model using the custom function
      hessian_1 <- get_hessiano_lm(model_1)
    }
    
    # Get the number of parameters (columns in the Hessian matrix)
    n_params <- hessian_1 |> ncol()
    
    # Get the number of observations (rows in the model's data)
    n_obs <- model_1$data |> nrow() 
    
    # Compute the determinant of the Hessian matrix and take the absolute value
    hessian_1_det <- hessian_1 |> Matrix::det() |> abs()
    
    # Compute the trace of the inverse Hessian matrix and take the absolute value
    m_trace <- hessian_1 |> Matrix::solve() |> diag() |> abs() |> sum()
    
    # Compute the ICOMP value using the specified formula
    icomp_value <- BIC(model_1) - (n_params * log(n_obs)) + log(hessian_1_det) + (n_params * log(m_trace / n_params))
    
    # Return the computed ICOMP value
    return(icomp_value)
  }
  
  # This function is taken from AICcmodavg
  AICc <-
    function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1, ...){
      
      if(is.null(nobs)) {
        n <- length(mod$fitted)
      } else {n <- nobs}
      
      LL <- logLik(mod)[1]
      K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
      
      if(c.hat == 1) {
        if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))} else{AICc <- -2*LL+2*K}
      }
      if(c.hat > 1 && c.hat <= 4) {
        K <- K+1
        if(second.ord==TRUE) {
          AICc <- (-2*LL/c.hat)+2*K*(n/(n-K-1))
          ##adjust parameter count to include estimation of dispersion parameter
        } else{
          AICc <- (-2*LL/c.hat)+2*K}
      }
      
      if(c.hat > 4) stop("High overdispersion and model fit is questionable\n")
      if(c.hat < 1) stop("You should set \'c.hat\' to 1 if < 1, but values << 1 might also indicate lack of fit\n")
      
      ##check if negative binomial and add 1 to K for estimation of theta if glm( ) was used
      if(!is.na(charmatch(x="Negative Binomial", table=family(mod)$family))) {
        if(!identical(class(mod)[1], "negbin")) { #if not negbin, add + 1 because k of negbin was estimated glm.convert( ) screws up logLik
          K <- K+1
          if(second.ord == TRUE) {
            AICc <- -2*LL+2*K*(n/(n-K-1))
          } else {
            AICc <- -2*LL+2*K
          }
        }
        if(c.hat != 1) stop("You should not use the c.hat argument with the negative binomial")
      }
      ##add 1 for theta parameter in negative binomial fit with glm( )
      
      ##check if gamma and add 1 to K for estimation of shape parameter if glm( ) was used
      if(identical(family(mod)$family, "Gamma") && c.hat > 1) stop("You should not use the c.hat argument with the gamma")
      
      ##an extra condition must be added to avoid adding a parameter for theta with negative binomial when glm.nb( ) is fit which estimates the correct number of parameters
      if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
      AICc
    }
  
  
  
  
  
})