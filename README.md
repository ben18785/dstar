
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dstar

<!-- badges: start -->
<!-- badges: end -->

d\* is (stochastic) statistic that quantifies the distance between two
distributions using draws from each

## Installation

You can install the development version of dstar from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ben18785/dstar")
```

## Example

If the two distributions are very similar d\* is close to 0.

``` r
library(dstar)

n <- 1000
draws_1 <- stats::rnorm(n, 0, 1)
draws_2 <- stats::rnorm(n, 0.1, 1)
dstar(draws_1, draws_2)
#> [1] 0.003333333
```

If the two distributions are different d\* is close to 1.

``` r
draws_1 <- stats::rnorm(n, 0, 1)
draws_2 <- stats::rnorm(n, 5, 1)
dstar(draws_1, draws_2)
#> [1] 0.9833333
```
