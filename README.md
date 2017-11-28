
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`tstmle`

[![Travis-CI Build
Status](https://travis-ci.org/podTockom/tstmle.svg?branch=master)](https://travis-ci.org/podTockom/tstmle)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/podTockom/tstmle?branch=master&svg=true)](https://ci.appveyor.com/project/podTockom/tstmle)
[![Coverage
Status](https://img.shields.io/codecov/c/github/podTockom/tstmle/master.svg)](https://codecov.io/github/podTockom/tstmle?branch=master)
[![Project Status: WIP - Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

> Estimation and Inference for Causal Effects with Single Time Series

**Authors:**

## What’s `tstmle`?

The `tstmle` package implements targeted maximum likelihood estimation
of causal effects based on the observation of a single time series
(i.e., interrupted time series analysis). The approach considers the
case where at each time point in a time series, a covariate vector,
treatment, and outcome are observed in chronological order. A family of
causal effects are defined as the impacts of stochastic interventions on
a subset of the treatment nodes on a future outcome. This general
formulation of the statistical estimation problem subsumes many other
important estimation problems, including but not limited to classical
time series models, group sequential adaptive designs, and even
independent and identically distributed data where the summary measure
of the past is simply the empty set. For details on the theoretical
underpinnings of the approach, the interested reader may consider
consulting (**???**). For a general introduction to the targeted
learning methodology and statistical causal inference, please consult
van der Laan and Rose (2011) and van der Laan and Rose (2017).

-----

## Installation

You can install a stable release of `tstmle` from GitHub via
[`devtools`](https://www.rstudio.com/products/rpackages/devtools/) with:

``` r
devtools::install_github("podTockom/tstmle")
```

<!--

In the future, the package will be available from
[CRAN](https://cran.r-project.org/) and can be installed via


```r
install.packages("tstmle")
```

-->

-----

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/podTockom/tstmle/issues).

-----

## Example

To illustrate how `tstmle` may be used to ascertain the effect of an
intervention on a single time series, consider the following example:

``` r
# Forthcoming...check back later.
```

-----

## License

© 2017

The contents of this repository are distributed under the MIT license.
See below for details:

    The MIT License (MIT)
    
    Copyright (c) 2017-2018
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

-----

## References

van der Laan, Mark J, and Sherri Rose. 2011. *Targeted Learning: Causal
Inference for Observational and Experimental Data*. Springer Science &
Business Media.

———. 2017. *Targeted Learning in Data Science: Causal Inference for
Complex Longitudinal Studies*. Springer Science & Business Media.
