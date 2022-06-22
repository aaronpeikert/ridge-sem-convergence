source(here::here("R", "funs.R"))
source(here::here("R", "data_grid.R"))

library(tidyverse)

string2partable <- function(model, cov){
  data <-
    MASS::mvrnorm(1e3,
                  mu = c(rep(0, times = ncol(cov))),
                  Sigma = cov,
                  empirical = TRUE)
  model_object <- lavaan::sem(model = model, sample.cov = cov, sample.nobs = 1000)
  lavaan::parametertable(model_object)
}

lpartable2jpartable <- function(lpartable){
  lv <- unique(dplyr::pull(dplyr::filter(lpartable, op == "=~"), "lhs"))
  jpartable <- mutate(lpartable,
                      lhs_copy = ifelse(op == "=~", rhs, lhs),
                      rhs_copy = ifelse(op == "=~", lhs, rhs),
                      op = ifelse(op == "=~", "~", op)) %>% 
    select(-lhs, -rhs) %>% 
    select(id, lhs = lhs_copy, op, rhs = rhs_copy, everything()) %>% 
    mutate(latent = lhs %in% lv)
  select(jpartable,
      from = lhs,
      parameter_type = op,
      to = rhs,
      free,
      start,
      label,
      estimate = est,
      latent
    )
}

true_covariances <- data_grid %>%
  select(-n_sample) %>% 
  unique()

true_covariances <- true_covariances %>% mutate(
    file = pmap(true_covariances, function(...)
      map2_chr(
        names(true_covariances), list(...), str_c
      )) %>%
      map(str_c, collapse = "") %>%
      map_chr(digest::digest),
    cov_file = file %>% 
      str_c("cov_files/", ., ".csv"),
    cov = true_covariances %>%
      pmap(matrices),
    model = map2_chr(n_pred, n_noise, ~lavaan_model(.x, .y)),
    lparameter_table = map2(model, cov, string2partable),
    jparameter_table = map(lparameter_table, lpartable2jpartable),
    par_file = file %>%
      str_c("par_files/", ., ".csv")
  )

write_mat <- function(mat, mat_file, ...){
  write_csv(as.data.frame(mat), mat_file)
}

fs::dir_create(here::here("cov_files"))
true_covariances %>%
  select(mat = cov, mat_file = cov_file) %>%
  pwalk(write_mat)

fs::dir_create(here::here("par_files"))
true_covariances %>%
  select(mat = jparameter_table, mat_file = par_file) %>%
  pwalk(write_mat)

data_grid <- data_grid %>% 
  left_join(true_covariances) %>% 
  select(-cov, -jparameter_table, -lparameter_table)

readr::write_csv(data_grid, "data_grid.csv")
capture.output(sessionInfo(), file = "RsessionInfo.txt")
