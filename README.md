
Deep Time-Series Learning
=========================

<!-- README.md is generated from README.Rmd. Please edit that file -->
> Estimation and Inference for Causal Effects with Single Binary Time Series

This repo provides basic implementation of the targeted maximum likelihood estimation of causal effects based on the observation of a single binary time series. The approach considers the case where at each time point in a time series, a covariate vector, treatment, and outcome are observed in chronological order. A family of causal effects are defined as the impacts of stochastic interventions on a subset of the treatment nodes on a future outcome. This general formulation of the statistical estimation problem subsumes many other important estimation problems, including but not limited to classical time series models, group sequential adaptive designs, and even independent and identically distributed data where the summary measure of the past is simply the empty set. For a general introduction to the targeted learning methodology and statistical causal inference, please consult van der Laan and Rose (2011) and van der Laan and Rose (2017). For details on the theoretical underpinnings of the approach, the interested reader may consider consulting van der Laan (2017) to start. Actual R package implementing a more general methodology and other time-series based target parameters, see [tstmle](https://github.com/podTockom/tstmle/).

------------------------------------------------------------------------

Example
-------

To illustrate how to ascertain the effect of an intervention on a single binary time series, consider the following example:

``` r
# Forthcoming...check back later.
```

------------------------------------------------------------------------

References
----------

van der Laan, Mark J. 2017. “Online Targeted Maximum Likelihood Estimation of Causal Effects Based on Observing a Single Time Series.” In *Targeted Learning in Data Science: Causal Inference for Complex Longitudinal Studies*, 263–97. Springer Science & Business Media.

van der Laan, Mark J, and Sherri Rose. 2011. *Targeted Learning: Causal Inference for Observational and Experimental Data*. Springer Science & Business Media.

———. 2017. *Targeted Learning in Data Science: Causal Inference for Complex Longitudinal Studies*. Springer Science & Business Media.
