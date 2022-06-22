data_grid <- tidyr::expand_grid(
  n_pred = 9,
  trait_t = 0.9, # loadings for manifest reference trait variables on trait
  trait_m = 0.4, # loadings for manifest method variables on trait
  method_m = 0.8, # loadings for manifest method variables on method
  n_trait = 3, # number of manifest trait variables
  n_method = 3, # number of manifest method variabels
  a = 0.2, # predictor loading on trait
  n_noise = c(0, 10, 25),
  Rsq = c(0.1, 0.2, 0.4, 0.6),
  coll = c(0, 0.2, 0.5, 0.9), # convergence issues are most pronounced for a value of 9
  n_sample = c(100, 250, 500, 1000)
)
