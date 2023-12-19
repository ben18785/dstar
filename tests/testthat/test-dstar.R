test_that("create_kde works", {
  # works with bkde (default)
  n <- 100000
  x <- stats::rnorm(n)
  f <- create_kde(x)
  f_1 <- f(1)
  f_0 <- f(0)
  expect_true(f_0 > f_1)

  # works with ash
  f <- create_kde(x, method = "ash1")
  f_1 <- f(1)
  f_0 <- f(0)
  expect_true(f_0 > f_1)

  # allows different bandwidths
  f <- create_kde(x, bandwidth = 0.1)
  f_1 <- f(1)
  f_0 <- f(0)
  expect_true(f_0 > f_1)

  # fails when different method called
  expect_error(create_kde(x, method = "test"))
})

test_that("probability_comparison works ok", {
  n <- 3
  probs1 <- rep(0.8, n)
  probs2 <- rep(0.2, n)
  expect_equal(probability_comparison(probs1, probs2), 1)

  probs1 <- c(0.8, NA, 0.8)
  expect_equal(probability_comparison(probs1, probs2), 2 / 3)

  probs2 <- c(0.2, NA, 0.2)
  expect_equal(probability_comparison(probs1, probs2), 1)
})

test_that("dstar throws errors when incorrect arguments given", {
  n <- 100000
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  density_prior <- data.frame(x = x1)
  expect_error(dstar(x1, x2, density_prior))
  expect_error(dstar(x1, x2, density_posterior = density_prior))
})

test_that("dstar works ok when using kde", {
  n <- 1000
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  d <- dstar(x1, x2)
  expect_true(d < 0.1)
  expect_true(d >= 0)

  x2 <- stats::rnorm(n, 10, 1)
  d <- dstar(x1, x2)
  expect_true(d > 0.9)
  expect_true(d <= 1)

  x2 <- stats::rnorm(n, 0, 10)
  d <- dstar(x1, x2)
  expect_true(d > 0.5)
  expect_true(d <= 1)

  # test runs with different kde
  d <- dstar(x1, x2, method = "ash1")
  expect_true(d > 0.5)
  expect_true(d <= 1)

  # test runs with different bandwidths
  d <- dstar(x1, x2,
    method = "ash1",
    bandwidth_prior = 0.1,
    bandwidth_posterior = 0.1
  )
  expect_true(d > 0.5)
  expect_true(d <= 1)

  # test runs with different training percentages
  d <- dstar(x1, x2, training_percent = 0.6)
  expect_true(d > 0.5)
  expect_true(d <= 1)
})

test_that("dstar works ok with long tailed distributions", {
  n <- 1000
  x1 <- stats::rnorm(n)
  x2 <- stats::rt(n, df = 0.1)
  d <- dstar(x1, x2)
  expect_true(d > 0.1)
  expect_true(d <= 1)
})

test_that("dstar works when supplying analytical densities", {
  n <- 1000
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)

  # supplying only prior density
  d <- dstar(x1, x2, density_prior = stats::dnorm)
  expect_true(d < 0.1)
  expect_true(d >= 0)

  # supplying only posterior density
  d <- dstar(x1, x2, density_posterior = stats::dnorm)
  expect_true(d < 0.1)
  expect_true(d >= 0)

  # supplying both
  d <- dstar(x1, x2,
    density_prior = stats::dnorm,
    density_posterior = stats::dnorm
  )
  expect_true(d < 0.1)
  expect_true(d >= 0)
})
