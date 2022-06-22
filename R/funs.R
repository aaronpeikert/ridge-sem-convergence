## function to specify correlation matrix depending on various variables

matrices <-
  function(n_pred,
           n_noise,
           Rsquare,
           coll,
           a,
           trait_t = 0.9,
           trait_m = 0.4,
           method_m = 0.8,
           n_trait = 3,
           n_method = 3,
           ...) {
    # number of predictors need to be a multiple of 3 so they can be equally divided into the sizes small, middle and large
    if (n_pred %% 3 != 0) {
      stop("Stop it! n_pred needs to be a multiple of 3")
    }
    ## constructing lambda matrix
    u <-
      c(rep(c(trait_t, rep(
        0.0, times = (n_pred + n_noise + 1)
      )), times = n_trait))
    v <-
      c(rep(c(trait_m, method_m, rep(
        0.0, times = (n_pred + n_noise)
      )), times = n_method))
    w <- c(u, v)
    lam <-
      matrix(
        w,
        nrow = n_trait + n_method,
        ncol = (n_pred + n_noise + 2),
        byrow = T
      )
    for (i in 1:(n_pred + n_noise)) {
      c <- c(a, rep(0.0, times = (n_pred + n_noise + 1)))
      c[i + 2] <-
        sqrt(1 - a ^ 2) ## communalities always need to add up to 1
      lam <- rbind(lam, c)
    }
    ## constructing phi matrix
    
    ########## execute this separately after fitting model to compare estimated standardized betas with the ones calculated here
    R <-
      diag(x = 1,
           nrow = n_pred / 3,
           ncol = n_pred / 3) # correlation matrix for predictors
    R <- replace(R, R == 0, coll)
    
    # first order correlation between each predictor and outcome given Rsquare and correlation among predictors
    r_out_small = sqrt(((1 / 9) * Rsquare) / (sum(solve(R))))
    r_out_middle = sqrt(((3 / 9) * Rsquare) / (sum(solve(R))))
    r_out_large = sqrt(((5 / 9) * Rsquare) / (sum(solve(R))))
    
    # standardized beta coeffcients
    beta_small = (1 / 9) * Rsquare / ((n_pred / 3) * r_out_small) # partial correlation between predictor and outcome
    beta_middle = (3 / 9) * Rsquare / ((n_pred / 3) * r_out_middle)
    beta_large = (5 / 9) * Rsquare / ((n_pred / 3) * r_out_large)
    #########
    
    phi <- diag(1, nrow = n_pred + 2, ncol = n_pred + 2)
    
    phi[2, 3:(2 + (n_pred / 3))] <- r_out_small
    phi[2, (3 + (n_pred / 3)):(2 + (2 * (n_pred / 3)))] <-
      r_out_middle
    phi[2, (3 + (2 * (n_pred / 3))):(2 + (3 * (n_pred / 3)))] <-
      r_out_large
    
    phi[3:(2 + (n_pred / 3)), 2] <- r_out_small
    phi[(3 + (n_pred / 3)):(2 + (2 * (n_pred / 3))), 2] <-
      r_out_middle
    phi[(3 + (2 * (n_pred / 3))):(2 + (3 * (n_pred / 3))), 2] <-
      r_out_large
    
    ## auxiliary matrix in order to set collinearity so that predictors with different sizes do not correlate
    G <- diag(x = 1, nrow = n_pred, ncol = n_pred)
    
    G[1:(n_pred / 3), 1:(n_pred / 3)] <-
      # only predictors of the same size are allowed to correlate
      replace(G[1:(n_pred / 3), 1:(n_pred / 3)],
              G[1:(n_pred / 3), 1:(n_pred / 3)] == 0, coll)
    
    G[((n_pred / 3) + 1):((n_pred / 3) * 2), ((n_pred / 3) + 1):((n_pred /
                                                                    3) * 2)] <-
      replace(G[((n_pred / 3) + 1):((n_pred / 3) * 2), ((n_pred / 3) + 1):((n_pred /
                                                                              3) * 2)],
              G[((n_pred / 3) + 1):((n_pred / 3) * 2), ((n_pred / 3) + 1):((n_pred /
                                                                              3) * 2)] == 0, coll)
    
    G[(((n_pred / 3) * 2) + 1):n_pred, (((n_pred / 3) * 2) + 1):n_pred] <-
      replace(G[(((n_pred / 3) * 2) + 1):n_pred, (((n_pred / 3) * 2) + 1):n_pred],
              G[(((n_pred / 3) * 2) + 1):n_pred, (((n_pred / 3) * 2) + 1):n_pred] == 0, coll)
    
    ## preliminary phi matrix for informative predictors
    phi[3:(n_pred + 2), 3:(n_pred + 2)] <- G
    
    ## final phi matrix
    H <-
      diag(x = 1,
           nrow = n_pred + n_noise + 2,
           ncol = n_pred + n_noise + 2)
    H[1:(n_pred + 2), 1:(n_pred + 2)] <- phi
    phi <- H
    
    ## constructing theta matrix
    
    theta <-
      diag(
        1,
        nrow = (n_trait + n_method + n_pred + n_noise),
        ncol = (n_trait + n_method + n_pred + n_noise)
      )
    
    com <-
      apply(
        lam,
        1,
        FUN = function(x)
          sum(x ^ 2)
      ) # calculating item communalities
    theta <-
      replace(theta, theta == 1, 1 - com) # adding variance to complement communalities to a value of 1
    
    ## calculating covariance (i.e. correlation) matrix
    
    cor = lam %*% phi %*% t(lam) + theta
    
    ## adding row and colnames corresponding to lavaan model
    rownames(cor) <-
      c(paste0("t", 1:n_trait),
        paste0("m", 1:n_method),
        paste0("x", 1:(n_pred + n_noise)))
    colnames(cor) <-
      c(paste0("t", 1:n_trait),
        paste0("m", 1:n_method),
        paste0("x", 1:(n_pred + n_noise)))
    round(cor, 2)
    
    return(cor)
  }


gen_data <- function(n_pred,
                     n_noise,
                     Rsquare,
                     coll,
                     a,
                     trait_t = 0.9,
                     trait_m = 0.4,
                     method_m = 0.8,
                     n_trait = 3,
                     n_method = 3,
                     n_sample = 100) {
  cormat <- matrices(
    n_pred = n_pred,
    n_noise = n_noise,
    Rsquare = Rsquare,
    coll = coll,
    a = a,
    trait_t = trait_t,
    trait_m = trait_m,
    method_m = method_m,
    n_trait = n_trait,
    n_method = n_method
  )
  as.data.frame(MASS::mvrnorm(n_sample, mu = c(rep(
    0, times = (6 + n_pred + n_noise)
  )), Sigma = cormat))
}

remove_spaces <- function(x){
  lines <- purrr::simplify(stringr::str_split(x, "\n"), "character")
  stringr::str_c(stringr::str_squish(lines), collapse = "\n")
}

lavaan_model <-
  function(n_pred, n_noise, n_allpred = n_pred + n_noise, ...) {
    # dots get discarded
    out <- paste(
      "### Trait variable ###
        T =~ 1*t1 + t2 + t3 +
        m1 + m2 + m3

        ### Method factor ###
        M =~ 1*m1 + m2 + m3
        T ~~ 0*M

        ### Regressing covariates on Trait ###
        T =~" ,
      paste0("x", 1:n_allpred, collapse = " + "),
      
      "\n\n### Residualizing predictors ###\n",
      paste0("c", 1:n_allpred, " =~ NA*", "x", 1:n_allpred, "\n", collapse = ""),
      "\n### restricting latent predictor variances ###\n",
      paste0("c", 1:n_allpred, " ~~ 1*", "c", 1:n_allpred, "\n", collapse = ""),
      "\n### restricting manifest predictor error variances ###\n",
      paste0("x", 1:n_allpred, " ~~ 0*", "x", 1:n_allpred, "\n", collapse = ""),
      "\n### Restricting correlations between Trait and latent predictors ###\nT ~~",
      paste0("0*c", 1:n_allpred, collapse = " + "),
      "\n\n### Regress method factor on residualized covariates ###\nM ~",
      paste0("c", 1:n_allpred, collapse = " + ")
    )
    remove_spaces(out)
  }

repeat_simulation <- function(n_repeat, data_grid, model_list){
  furrr::future_map_dfr(seq_len(n_repeat), function(discard)simulation(data_grid, model_list),
                        .options = furrr::furrr_options(seed = TRUE))
}

simulation <- function(data_grid, model_list){
  results <- simulation_(data_grid, model_list)
  tidyr::unnest(dplyr::mutate(data_grid, results = results), results)
}

simulation_ <- function(data_grid, model_list){
  furrr::future_pmap(data_grid, fit_models, models = model_list,
                     .options = furrr::furrr_options(seed = TRUE))
}

fit_models <- function(models, ...){
  lavaan_model <- lavaan_model(...)
  d <- gen_data(...)
  purrr::map_dfr(models, fit_model, lavaan_model, data = d, .id = "model")
}

fit_model <- function(model, lavaan_model, data){
  if(is.function(model$fun))fit_fun <- model$fun
  else fit_fun <- rlang::eval_tidy(rlang::parse_expr(model$fun))
  model$args$data <- data
  model$args$lavaan_model <- lavaan_model
  exec(fit_fun, !!!model$args)
}
