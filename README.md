
R/`tstmle01`
============

[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip) [![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

<!-- README.md is generated from README.Rmd. Please edit that file -->
> Estimation and Inference for Causal Effects with Single Binary Time Series

What's `tstmle01`?
------------------

The `tstmle01` package implements the targeted maximum likelihood estimation (TMLE) of the marginal causal effect based on the observation of a single binary time-series (van der Laan and Rose 2017). Current implementation supports iterative TMLE; check back for one-step and online TMLE (van der Laan and Gruber 2016), (van der Laan and Lendle 2014).

For the (faster) package implementing a more general methodology and other time-series based target parameters, see [tstmle](https://github.com/podTockom/tstmle/). In particular, [tstmle](https://github.com/podTockom/tstmle/) implements the data-dependent, *C*<sub>*o*(*t*)</sub>-specific, causal effect and the adaptive design for learning the optimal treatment rule within a single time series (van der Laan and Malenica 2018). Here, initial estimation is based on the [sl3](https://github.com/jeremyrcoyle/sl3) package, which constructs ensemble models with proven optimality properties for time-series data (Malenica and van der Laan 2018).

We emphasize that this general formulation of the statistical estimation problem subsumes many other important estimation problems, including but not limited to classical time series models, group sequential adaptive designs, and even independent and identically distributed data when the summary measure of the past is simply the empty set.

------------------------------------------------------------------------

Installation
------------

You can install a stable release of `tstmle01` from GitHub via [`devtools`](https://www.rstudio.com/products/rpackages/devtools/) with:

``` r
devtools::install_github("podTockom/tstmle01")
```

------------------------------------------------------------------------

Issues
------

If you encounter any bugs or have any specific feature requests, please [file an issue](https://github.com/podTockom/tstmle01/issues).

------------------------------------------------------------------------

Example
-------

To illustrate how to ascertain the effect of an intervention on a single binary time series, consider the following example:

``` r
#Simulated data based on the simcausal package:
load(data)

#Estimate of the expected value of the outcome at time 5, under intervention on Anode 3:
res<-tstmle01(data,freqY=3,freqA=3,freqW=3,t=5,Anode=3,intervention1=1)
res$psi
```

------------------------------------------------------------------------

License
-------

© 2017 [Ivana Malenica](https://github.com/podTockom)

The contents of this repository are distributed under the MIT license. See below for details:

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

------------------------------------------------------------------------

References
----------

Malenica, Ivana, and Mark J van der Laan. 2018. “Oracle Inequality for Cross-Validation Estimator Selector for Dependent Time-Ordered Experiments.”

van der Laan, Mark J, and Susan Gruber. 2016. “One-Step Targeted Minimum Loss-Based Estimation Based on Universal Least Favorable One-Dimensional Submodels.” Working Paper 347. U.C. Berkeley Division of Biostatistics Working Paper Series.

van der Laan, Mark J, and Samuel D Lendle. 2014. “Online Targeted Learning.” Working Paper 330. U.C. Berkeley Division of Biostatistics Working Paper Series.

van der Laan, Mark J, and Ivana Malenica. 2018. “Robust Estimation of Data-Dependent Causal Effects Based on Observing a Single Time-Series.”

van der Laan, Mark J, and Sherri Rose. 2017. *Targeted Learning in Data Science: Causal Inference for Complex Longitudinal Studies*. Springer Science & Business Media.
