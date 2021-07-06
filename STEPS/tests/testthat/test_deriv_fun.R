# "deriv_fun"

library(STEPS)




test_that("test get_number_IDs returns right number", {
  function_output <- .get_number_ids(data.frame("ID" = c(1, 2), "other" = c(3, 4)))
  expect_equal(2, function_output)

})


test_that("test get_number_of_pairs returns right number", {
  source(file.path('tests_utils.R'))
  config <- STEPS_configuration(n = 3)
  X_table <- x_sim_fun(config$IDs_tab, config$difficulties_per_ID, pX = 2)
  ordered_sets <- SIM_sets_fun(data = X_table, obs_cov = config$obs_cov)

  function_output <- .get_number_of_pairs(ordered_sets$setB)
  expect_equal(1, function_output)


})

test_that("test averaging repeated measures returns all samples and right calculation", {
  source(file.path('tests_utils.R'))
  n_subjects <- 3
  n_levels <- 3
  config <- STEPS_configuration(n = n_subjects)
  X_table <- x_sim_fun(config$IDs_tab, config$difficulties_per_ID, pX = 2)
  ordered_sets <- SIM_sets_fun(data = X_table, obs_cov = config$obs_cov)

  function_output0.5 <- .average_repeated_measures(c(rep(0, n_subjects * n_levels), rep(1, n_subjects * n_levels)), ordered_sets$setR)
  function_output0 <- .average_repeated_measures(c(rep(-1, n_subjects * n_levels), rep(1, n_subjects * n_levels)), ordered_sets$setR)

  expect_equal(length(function_output0.5), n_subjects * n_levels)
  expect_true(all(function_output0.5 == rep(0.5, n_subjects * n_levels)))


  expect_true(all(function_output0 == rep(0, n_subjects * n_levels)))
})


test_that("test global options of datatable auto index is False", {
  expect_false(getOption("datatable.auto.index"))
})


test_that("test calculate deriv of B, W and R terms returns what expected", {
  source(file.path('tests_utils.R'))
  n_subjects <- 3
  n_levels <- 3
  n_repetitions <- 2
  config <- STEPS_configuration(n = n_subjects)
  X_table <- x_sim_fun(config$IDs_tab, config$difficulties_per_ID, pX = 2)
  ordered_sets <- SIM_sets_fun(data = X_table, obs_cov = config$obs_cov)

  Fx <- rep(0, n_subjects * n_levels * n_repetitions)
  n.id <- .get_number_ids(ordered_sets$setR, ID_col = "ID")
  n_B <- .get_number_of_pairs(ordered_sets$setB)
  Fx.mean <- .average_repeated_measures(Fx, ordered_sets$setR)


  b_deriv <- .calculate_deriv_b_term(ordered_sets$setB, ordered_sets$setR, Fx.mean)
  expect_true(all((b_deriv %>% filter(ID == 1))$b.term == 1 / (n_repetitions * n_levels) ) )
  expect_true(all((b_deriv %>% filter(ID == 3))$b.term == -1 / (n_repetitions * n_levels) ) )
  expect_true(all((b_deriv %>% filter(ID == 2))$b.term == 0) )

  w_deriv <- .calculate_deriv_w_term(ordered_sets$setW, ordered_sets$setR, Fx.mean)
  expect_true(all((w_deriv %>% filter(ID == 1, difficulty1 == 0))$w.term == 2 / (n_repetitions * n_levels) ) )
  expect_true(all((w_deriv %>% filter(ID == 1, difficulty1 == 2))$w.term == -2 / (n_repetitions * n_levels) ) )
  expect_true(all((w_deriv %>% filter(ID == 1, difficulty1 == 1))$w.term == 0) )

  r_deriv <- .calculate_deriv_r_term(ordered_sets$setR, Fx, Fx.mean)
  expect_true(all(r_deriv$r.term == 0 ) )

  full_deriv <- deriv_fun(Fx = Fx,
                          setB = ordered_sets$setB, setW = ordered_sets$setW, setR = ordered_sets$setR,
                          lambdaB = 1, lambdaW = 1, lambdaR = 1, h = 0.5)
  expect_true(all(full_deriv == ((1 / n_B) * b_deriv$b.term) + ((1 / n.id) * w_deriv$w.term) + ((1 / n.id) * r_deriv$r.term)))

  Fx_zero_and_one <- c(rep(0, n_subjects * n_levels), rep(1, n_subjects * n_levels))
  Fx.mean <- .average_repeated_measures(Fx_zero_and_one, ordered_sets$setR)
  r_deriv <- .calculate_deriv_r_term(ordered_sets$setR, Fx_zero_and_one, Fx.mean)
  expect_true(all((r_deriv %>% filter(R == 1))$r.term == -1 / (n_levels * n_repetitions ) ))
  expect_true(all((r_deriv %>% filter(R == 2))$r.term == 1 / (n_levels * n_repetitions ) ))

})

