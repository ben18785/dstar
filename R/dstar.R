#' Replicates ash1 but without print statement
#'
#' @param bins number of bins
#'
#' @return list which is useable by kde methods
ash2 <- function(bins) {
  m <- 5 # not sure on what these mean
  kopt <- c(2, 2) # not sure on what these mean
  nc <- bins$nc
  ab <- bins$ab
  nbin <- length(nc)
  r <- .Fortran("ash1", as.integer(m), as.integer(nc), as.integer(nbin),
    as.double(ab), as.integer(kopt),
    t = double(nbin), f = double(nbin),
    double(m), ier = integer(1), PACKAGE = "ash"
  )
  list(x = r$t, y = r$f, m = m, ab = ab, kopt = kopt, ier = r$ier)
}

#' Create kde give a set of draws
#'
#' @param x a set of univariate draws
#' @param method a kde method (either bkde or ash1)
#' @param bandwidth bandwidth of kde
#' @param ... other arguments passed to kde method
#'
#' @return a function which when evaluated returns estimated density at a point
create_kde <- function(x, method = "bkde", bandwidth = NULL, ...) {
  if (!method %in% c("bkde", "ash1")) {
    stop("Method must either be bkde or ash1.")
  }

  bandwidth1 <- bandwidth
  if (method == "bkde") {
    if (is.null(bandwidth)) {
      a <- KernSmooth::bkde(x, ...)
    } else {
      a <- KernSmooth::bkde(x, bandwidth = bandwidth1, ...)
    }
  } else if (method == "ash1") {
    if (is.null(bandwidth)) {
      a <- ash2(ash::bin1(x), ...)
    } else {
      a <- ash2(ash::bin1(x, nbin = round(bandwidth1 * 50)), ...)
    }
  }

  f <- stats::approxfun(a$x, a$y, method = "constant")

  kde <- function(x) {
    f_val <- f(x)
    dplyr::if_else(is.na(f_val), 0, f_val)
  }
  kde
}



#' Compares probability densities across one reference model to another
#'
#' @param prob1 Set of probability density values from one model (this is reference set)
#' @param prob2 Set of probability density values from another model (this is comparison set)
#'
#' @return the fraction of cases where prob1 > prob2
probability_comparison <- function(prob1, prob2) {
  k <- 0
  n <- length(prob1)
  for (i in seq_along(prob1)) {
    if ((!is.na(prob1[i])) & (!is.na(prob2[i]))) {
      k <- k + dplyr::if_else(prob1[i] > prob2[i], 1, 0)
    } else if (!is.na(prob1[i])) {
      k <- k + 1
    } else if (is.na(prob1[i]) & is.na(prob2[i])) {
      n <- n - 1
    }
  }

  k / n
}

#' Calculates the d* measure of information gain from after seeing the data in
#' Bayesian inference
#'
#' @param draws_prior a set of univariate draws of a parameter from the prior
#' @param draws_posterior a set of univariate draws of a parameter from the posterior
#' @param density_prior the density representing the univariate prior on the parameter, if known
#' @param density_posterior the density representing the univariate posterior on the parameter,
#' if known
#' @param method the method for performing kernel density estimation (either bkde or ash1)
#' @param training_percent the fraction of draws to use for training set (defaults to 0.7)
#' with the remainder comprising the testing set
#' @param bandwidth_prior the bandwidth to use if/when fitting a kernel density estimator to the draws in draws_prior
#' @param bandwidth_posterior the bandwidth to use if/when fitting a kernel density estimator to the draws in draws_posterior
#' @param ... other arguments to pass to create_kde
#'
#' @return a d* value, which is between 0 and 1. A value of 0 indicates the prior
#' and posterior and identical; a value of 1 indicates they don't overlap at all
#' @export
#'
#' @examples
#' # shows how to use function supposing prior is theta ~ normal(0, 1) and posterior is
#' # theta|data ~ normal(1, 1)
#' n <- 1000 # number of prior and posterior draws you generate
#' draws_prior <- rnorm(n, 0, 1)
#' draws_posterior <- rnorm(n, 1, 1)
#' dstar(draws_prior, draws_posterior)
dstar <- function(draws_prior, draws_posterior,
                  density_prior = NULL, density_posterior = NULL,
                  method = "bkde",
                  training_percent = 0.7,
                  bandwidth_prior = NULL,
                  bandwidth_posterior = NULL,
                  ...) {
  # check if density_prior or density_posterior are defined (i.e. a probability density)
  # if so check they must be functions
  if (!is.null(density_prior)) {
    isfunc1 <- is.function(density_prior)
  } else {
    isfunc1 <- NULL
  }
  if (!is.null(density_posterior)) {
    isfunc2 <- is.function(density_posterior)
  } else {
    isfunc2 <- NULL
  }

  if (!is.null(isfunc1)) {
    if (!isfunc1) {
      stop("density_prior must be a function or NULL type.")
    }
  }
  if (!is.null(isfunc2)) {
    if (!isfunc2) {
      stop("density_posterior must be a function or NULL type.")
    }
  }
  if (is.null(isfunc1)) {
    isfunc1 <- FALSE
  }
  if (is.null(isfunc2)) {
    isfunc2 <- FALSE
  }

  if (isfunc1 && isfunc2) {
    # use whole sets given
    prob11 <- purrr::map_dbl(draws_prior, density_prior)
    prob12 <- purrr::map_dbl(draws_prior, density_posterior)
    prob21 <- purrr::map_dbl(draws_posterior, density_prior)
    prob22 <- purrr::map_dbl(draws_posterior, density_posterior)
  } else {
    # training and testing sets
    idxs <- seq_along(draws_prior)
    training_ids <- sample(idxs,
      size = round(training_percent * length(draws_prior))
    )
    testing_ids <- setdiff(idxs, training_ids)
    draws_prior_train <- draws_prior[training_ids]
    draws_posterior_train <- draws_posterior[training_ids]
    draws_prior_test <- draws_prior[testing_ids]
    draws_posterior_test <- draws_posterior[testing_ids]
    if (isfunc1) {
      density_posterior <- create_kde(draws_posterior_train,
        method,
        bandwidth = bandwidth_posterior,
        ...
      )
    } else if (isfunc2) {
      density_prior <- create_kde(draws_prior_train,
        method,
        bandwidth = bandwidth_prior,
        ...
      )
    } else {
      density_prior <- create_kde(draws_prior_train,
        method,
        bandwidth = bandwidth_prior,
        ...
      )
      density_posterior <- create_kde(draws_posterior_train,
        method,
        bandwidth = bandwidth_posterior,
        ...
      )
    }
    prob11 <- purrr::map_dbl(draws_prior_test, density_prior)
    prob12 <- purrr::map_dbl(draws_prior_test, density_posterior)
    prob21 <- purrr::map_dbl(draws_posterior_test, density_prior)
    prob22 <- purrr::map_dbl(draws_posterior_test, density_posterior)
  }
  # calculate accuracy
  acc1 <- probability_comparison(prob11, prob12)
  acc2 <- probability_comparison(prob22, prob21)
  acc_overall <- 0.5 * (acc1 + acc2)

  # calculate score
  score <- 2 * (acc_overall - 0.5)
  score <- dplyr::if_else(score < 0, 0, score)
  score
}
