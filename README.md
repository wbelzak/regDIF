
<!-- README.md is generated from README.Rmd. -->

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/wbelzak/regDIF.svg?branch=master)](https://travis-ci.org/wbelzak/regDIF)
[![Codecov test
coverage](https://codecov.io/gh/wbelzak/regDIF/branch/master/graph/badge.svg)](https://codecov.io/gh/wbelzak/regDIF?branch=master)
<!-- badges: end -->

# regDIF: Regularized Differential Item Functioning

This R package performs regularization of differential item functioning
(DIF) parameters in item response theory (IRT) models using a penalized
expectation-maximization algorithm.

## Version 1.1.0 Features

regDIF can:

  - Handle multiple continuous and categorical DIF covariates;
  - Support binary, ordinal, and continuous item responses;
  - Use LASSO, ridge, MCP, elastic net, and group penalty functions for
    regularization;
  - Support parallel processing to improve model estimation speed;
  - Allow for proxy data to be used in place of estimating latent
    variable scores, which leads to much faster estimation speed.

## Installation

To get the current released version from CRAN:

``` r
install.packages("regDIF")
```

To get the current development version from Github:

``` r
# install.packages("devtools")
devtools::install_github("wbelzak/regDIF")
```

## Getting Started

A simulated data example with 6 item responses (binary) and 3 background
variables (gender, age, study) is available in the regDIF package:

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

First, the item responses and predictor values are separately specified:

``` r
item.data <- ida[, 1:6]
pred.data <- ida[, 7:9]
```

Second, the `regDIF()` function fits a sequence of 10 tuning parameter
values using a penalized EM algorithm, which assumes a normal latent
variable affects all item responses:

``` r
fit <- regDIF(item.data, pred.data, num.tau = 10)
```

The DIF results are shown below:

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

When estimation speed is slow, proxy data may be used in place of latent
score estimation:

``` r
fit_proxy <- regDIF(item.data, pred.data, prox.data = rowSums(item.data))
```

``` r
summary(fit_proxy)
#> Call:
#> regDIF(item.data = item.data, pred.data = pred.data, prox.data = rowSums(item.data))
#> 
#> Optimal model (out of 100):
#>          tau          bic 
#>    0.2241118 3550.1030000 
#> 
#> Non-zero DIF effects:
#> item3.int.gender    item4.int.age    item5.int.age item5.int.gender 
#>            0.121            0.339           -0.192           -0.589 
#>  item5.int.study item6.int.gender item2.slp.gender  item3.slp.study 
#>            0.519            0.035            0.144           -0.062 
#> item5.slp.gender 
#>           -0.126
```

Other penalty functions (besides LASSO) may also be used. For instance,
the elastic net penalty uses a second tuning parameter, `alpha`, to vary
the ratio of LASSO to ridge penalties:

``` r
fit_proxy_net <- regDIF(item.data, pred.data, prox.data = rowSums(item.data), alpha = .5)
```

``` r
summary(fit_proxy_net)
#> Call:
#> regDIF(item.data = item.data, pred.data = pred.data, prox.data = rowSums(item.data), 
#>     alpha = 0.5)
#> 
#> Optimal model (out of 100):
#>          tau          bic 
#>    0.4377184 3562.5710000 
#> 
#> Non-zero DIF effects:
#> item3.int.gender    item4.int.age    item5.int.age item5.int.gender 
#>            0.104            0.277           -0.209           -0.466 
#>  item5.int.study item6.int.gender item2.slp.gender  item3.slp.study 
#>            0.369            0.033            0.119           -0.054 
#> item5.slp.gender 
#>           -0.099
```

## Questions

Please send any questions to <wbelzak@gmail.com>.
