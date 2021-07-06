# "ord_pairs_fun"
library(STEPS)

test_that("ord pairs fun return right order", {
  pair <- c(1, 2)
  sort_by <- c(1, 2)
  expect_identical(c(1,2), .ord_pairs_fun(pair = pair, sort.by = sort_by))

  pair <- c(1, 2)
  sort_by <- c(2, 1)
  expect_identical(c(2, 1), .ord_pairs_fun(pair = pair, sort.by = sort_by))

  pair <- c(1, 2)
  sort_by <- c(1, 1)
  expect_identical(NULL, .ord_pairs_fun(pair = pair, sort.by = sort_by))
})



