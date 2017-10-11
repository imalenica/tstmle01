
<!-- README.md is generated from README.Rmd. Please edit that file -->
R/`tstmle`
==========

\[\[Travis-CI Build Status\] \[\[AppVeyor Build Status\] \[\[Coverage Status\]\[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)\](<http://opensource.org/licenses/MIT>)

> Causal Effect based on a Single Time Series

**Authors:**

Introduction
------------

This package implements targeted maximum likelihood estimation of causal effects based on observing a single time series. We consider the case where at each time within a time series we observe in chronological order a covariate vector, treatment, and an outcome. We define a family of causal effects in terms of stochastic interventions on a subset of the treatment nodes on a future outcome. This general formulation of the statistical estimation problem includes many important estimation problems including classical time series models, group sequential adaptive designs, and even independent and identically distributed data when the summary measure of the past is an empty set.

Installation
------------

You can install a stable release of `tstmle` from GitHub via [`devtools`](https://www.rstudio.com/products/rpackages/devtools/) with:

``` r
devtools::install_github("podTockom/tstmle")
```

In the future, the package will be available from [CRAN](https://cran.r-project.org/) via

``` r
install.packages("tstmle")
```

Issues
------

If you encounter any bugs or have any specific feature requests, please [file an issue](https://github.com/podTockom/tstmle/issues).

Example
-------

License
-------

© 2017

The contents of this repository are distributed under the MIT license. See below for details:

    The MIT License (MIT)

    Copyright (c) 2016-2017 

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
