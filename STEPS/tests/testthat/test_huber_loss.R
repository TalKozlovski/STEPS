context("huber_loss and grad huber_loss")
library(STEPS)

test_that("huber_loss returns expected values", {
  expect_identical(0, huber_loss(-1, h = 0.5))
  expect_identical(0, huber_loss(-0.5, h = 0.5))
  expect_identical(0.125, huber_loss(0, h = 0.5))
  expect_identical(0.5, huber_loss(0.5, h = 0.5))
  expect_identical(1, huber_loss(1, h = 0.5))



  expect_identical(c(0, 0.125, 0.5, 1), huber_loss(a = c(-1, 0, 0.5, 1), h = 0.5))
})




test_that("grad_huber_loss returns expected values", {
  expect_identical(0, grad_huber_loss(-1, h = 0.5))
  expect_identical(0, grad_huber_loss(-0.5, h = 0.5))
  expect_identical(0.5, grad_huber_loss(0, h = 0.5))
  expect_identical(1, grad_huber_loss(0.5, h = 0.5))
  expect_identical(1, grad_huber_loss(1, h = 0.5))



  expect_identical(c(0, 0, 0.5, 0.75, 1, 1), grad_huber_loss(a = c(-1, -0.5, 0, 0.25,  0.5, 1), h = 0.5))
})



