
<!-- README.md is generated from README.Rmd. -->
<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/wbelzak/regDIF.svg?branch=master)](https://travis-ci.org/wbelzak/regDIF)
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
- Allow for proxy data to be used in place of estimating latent variable
  scores, which leads to much faster estimation speed.

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
#>    0.1753246 4081.6941000 
#> 
#> Non-zero DIF effects:
#>    item4.int.age    item5.int.age item5.int.gender  item5.int.study 
#>           0.2153          -0.0897          -0.5717           0.6018 
#>  item4.slp.study item5.slp.gender 
#>          -0.0936          -0.1764
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
#>    0.2766486 3540.8070000 
#> 
#> Non-zero DIF effects:
#> item3.int.gender    item4.int.age item5.int.gender  item5.int.study 
#>           0.0955           0.2200          -0.5118           0.7040 
#> item2.slp.gender  item4.slp.study item5.slp.gender 
#>           0.1102          -0.1413          -0.1384
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
#>    0.5685967 3563.7495000 
#> 
#> Non-zero DIF effects:
#> item3.int.gender    item4.int.age    item5.int.age item5.int.gender 
#>           0.0681           0.1672          -0.0939          -0.3463 
#>  item5.int.study item2.slp.gender  item4.slp.study item5.slp.gender 
#>           0.4346           0.0778          -0.1172          -0.1379
```

## Questions

Please send any questions to <wbelzak@gmail.com>.
