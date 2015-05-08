<!-- README.md is generated from README.Rmd. Please edit that file -->
brms
====

The <b>brms</b> package provides an interface to fit bayesian generalized linear mixed models using Stan, which is a C++ package for obtaining Bayesian inference using the No-U-turn sampler (see <http://mc-stan.org/>). The formula syntax is very similar to that of the package lme4 to provide a familiar and simple interface for performing regression analyses.

How to use brms
===============

``` r
library(brms)
```

As a simple example, we use poisson regression to model the seizure counts in epileptic patients to investigate whether the treatment (represented by variable Trt\_c) can reduce the seizure counts. A random intercept is incorporated to account for the variance between patients.

``` r
fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient), data = epilepsy, family = "poisson")
```

If rstan is not installed, brm will return the Stan model, the required data, and the parameters of interest, which are the important prerequisites to fit the model in Stan. If rstan is installed, the model is fitted automatically and the results (i.e. posterior samples) can be investigated.

<!-- 

```r
print(fit) 
#> Inference for Stan model: model.
#> 2 chains, each with iter=2000; warmup=500; thin=1; 
#> post-warmup draws per chain=1500, total post-warmup draws=3000.
#> 
#>                         mean se_mean   sd    2.5%     25%     50%     75%   97.5% n_eff Rhat
#> b_Intercept             1.61    0.00 0.08    1.46    1.56    1.61    1.66    1.76  3000    1
#> b_log_Age_c             0.48    0.01 0.38   -0.30    0.23    0.49    0.73    1.22   665    1
#> b_log_Base4_c           1.06    0.00 0.11    0.84    0.99    1.06    1.14    1.28   766    1
#> b_Trt_c                -0.34    0.01 0.16   -0.66   -0.45   -0.34   -0.23   -0.01   817    1
#> b_log_Base4_c__Trt_c    0.33    0.01 0.22   -0.10    0.18    0.33    0.48    0.76   821    1
#> sd_patient_Intercept    0.55    0.00 0.07    0.43    0.50    0.55    0.59    0.69  3000    1
#> lp__                 3197.39    0.18 5.96 3185.11 3193.39 3197.60 3201.53 3208.74  1053    1
#> 
#> Samples were drawn using NUTS(diag_e) at Fri May 08 09:56:03 2015.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

```r
brm.plot(fit) 
```

<img src="README-plot-1.png" title="" alt="" style="display: block; margin: auto;" />
-->
How to install brms
===================

``` r
install.packages("brms")
```

Without having rstan installed, the function brm will return the Stan model, the required data, and the parameters of interest. To allow brm to fit the model automatically, the package rstan has to be installed manually, as it is not on CRAN, yet. However, the developers of Stan and rstan are currently working on a version to be uploaded on CRAN. In the meantime, instructions on how to install rstan can be found at <https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started>.

<!-- Before you will be able to actually fit bayesian models with brms, the package rstan has to be installed manually, as it is not on CRAN, yet. First, you need a C++ compiler. See https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started#prerequisites for instructions on how to get one. Second, install rstan by running the following R code (the number behind 'j' in the first line corresponds to the number of cores to use for the installation). This may take a few minutes and you should restart R after the installation.


```r
Sys.setenv(MAKEFLAGS = "-j1") 
source('http://mc-stan.org/rstan/install.R', echo = TRUE, max.deparse.length = 2000)
install_rstan()
```
-->
