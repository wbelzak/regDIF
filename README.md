
<!-- README.md is generated from README.Rmd. -->

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/wbelzak/regDIF.svg?branch=master)](https://travis-ci.org/wbelzak/regDIF)
[![Codecov test
coverage](https://codecov.io/gh/wbelzak/regDIF/branch/master/graph/badge.svg)](https://codecov.io/gh/wbelzak/regDIF?branch=master)
<!-- badges: end -->

# regDIF: Regularized Differential Item Functioning

An R package that performs regularization of differential item
functioning (DIF) parameters in item response theory (IRT) models using
a penalized expectation-maximization algorithm.

## Version 1.0.0 Features

  - Handles multiple continuous and categorical DIF covariates.
  - Supports binary and ordinal item responses.
  - Includes LASSO, ridge, MCP, and elastic net penalty functions for
    regularization.

## Installation

To get the current development version from Github:

``` r
# install.packages("devtools")
devtools::install_github("wbelzak/regDIF")
```

## Getting Started

A simulated data example with 6 item responses (binary) and 3 background
variables (gender, age, study) is used to demonstrate the `regDIF`
package.

``` r
library(regDIF)
head(ida)
#>   item1 item2 item3 item4 item5 item6 age gender study
#> 1     0     0     0     0     0     0  -2     -1    -1
#> 2     0     0     0     0     0     0   0     -1    -1
#> 3     0     0     0     0     0     0   3     -1    -1
#> 4     0     1     1     1     1     1   1     -1    -1
#> 5     0     0     0     0     0     0  -2     -1    -1
#> 6     1     0     0     0     0     0   1     -1    -1
```

First, the item responses and predictor values are separately specified.

``` r
item.data <- ida[, 1:6]
pred.data <- ida[, 7:9]
```

Second, the `regDIF()` function fits a sequence of 10 tuning parameter
values.

``` r
fit <- regDIF(item.data, pred.data, num.tau = 10)
```

Finally, the DIF results are shown.

``` r
summary(fit)
#> Call:
#> regDIF(item.data = item.data, pred.data = pred.data, num.tau = 10)
#> 
#> Optimal model (out of 10):
#>          tau          bic 
#>    0.1934603 4074.7270000 
#> 
#> Non-zero DIF effects:
#>    item4.int.age    item5.int.age item5.int.gender  item5.int.study 
#>            0.263           -0.247           -0.656            0.349
```

## Questions

Please send any questions to <wbelzak@gmail.com>.
