<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build Status](https://travis-ci.org/paul-buerkner/brms.svg?branch=master)](https://travis-ci.org/paul-buerkner/brms) [![Coverage Status](https://codecov.io/github/paul-buerkner/brms/coverage.svg?branch=master)](https://codecov.io/github/paul-buerkner/brms?branch=master) [![CRAN Version](http://www.r-pkg.org/badges/version/brms)](https://cran.r-project.org/package=brms)

brms
====

The **brms** package provides an interface to fit Bayesian generalized (non-)linear multilevel models using Stan, which is a C++ package for performing full Bayesian inference (see <http://mc-stan.org/>). The formula syntax is very similar to that of the package lme4 to provide a familiar and simple interface for performing regression analyses. A wide range of distributions and link functions are supported, allowing users to fit -- among others -- linear, robust linear, count data, survival, response times, ordinal, zero-inflated, hurdle, and even self-defined mixture models all in a multilevel context. Further modeling options include non-linear and smooth terms, auto-correlation structures, censored data, meta-analytic standard errors, and quite a few more. In addition, all parameters of the response distribution can be predicted in order to perform distributional regression. Prior specifications are flexible and explicitly encourage users to apply prior distributions that actually reflect their beliefs. Model fit can easily be assessed and compared with posterior predictive checks and leave-one-out cross-validation.

<!--

-->
How to use brms
===============

``` r
library(brms)
```

As a simple example, we use poisson regression to model the seizure counts in epileptic patients to investigate whether the treatment (represented by variable Trt\_c) can reduce the seizure counts. Two group-level intercepts are incorporated to account for the variance between patients as well as for the residual variance.

``` r
fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient) + (1|obs), 
           data = epilepsy, family = "poisson")
#> Compiling the C++ model
#> Start sampling
```

The results (i.e. posterior samples) can be investigated using

``` r
summary(fit, waic = TRUE) 
#>  Family: poisson(log) 
#> Formula: count ~ log_Age_c + log_Base4_c * Trt_c + (1 | patient) + (1 | obs) 
#>    Data: epilepsy (Number of observations: 236) 
#> Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1; 
#>          total post-warmup samples = 4000
#>     ICs: LOO = Not computed; WAIC = 1145.42
#>  
#> Group-Level Effects: 
#> ~obs (Number of levels: 236) 
#>               Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
#> sd(Intercept)     0.37      0.04     0.29     0.46       1239    1
#> 
#> ~patient (Number of levels: 59) 
#>               Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
#> sd(Intercept)     0.51      0.07     0.38     0.66       1182    1
#> 
#> Population-Level Effects: 
#>                   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
#> Intercept             1.56      0.08     1.40     1.72       1208    1
#> log_Age_c             0.48      0.37    -0.26     1.21        760    1
#> log_Base4_c           1.06      0.11     0.85     1.28       1202    1
#> Trt_c                -0.33      0.16    -0.63    -0.03       1203    1
#> log_Base4_c:Trt_c     0.35      0.22    -0.08     0.78       1212    1
#> 
#> Samples were drawn using sampling(NUTS). For each parameter, Eff.Sample 
#> is a crude measure of effective sample size, and Rhat is the potential 
#> scale reduction factor on split chains (at convergence, Rhat = 1).
```

On the top of the output, some general information on the model is given, such as family, formula, number of iterations and chains, as well as the WAIC, which is an information criterion for Bayesian models. Next, group-level effects are displayed seperately for each grouping factor in terms of standard deviations and (in case of more than one group-level effect per grouping factor; not displayed here) correlations between group-level effects. On the bottom of the output, population-level effects are displayed. If incorporated, autocorrelation effects and family specific parameters (e.g., the residual standard deviation 'sigma' in normal models) are also given.

In general, every parameter is summarized using the mean ('Estimate') and the standard deviation ('Est.Error') of the posterior distribution as well as two-sided 95% credible intervals ('l-95% CI' and 'u-95% CI') based on quantiles. The last two values ('Eff.Sample' and 'Rhat') provide information on how well the algorithm could estimate the posterior distribution of this parameter. If 'Rhat' is considerably greater than 1, the algorithm has not yet converged and it is necessary to run more iterations and / or set stronger priors.

To visually investigate the chains as well as the posterior distributions, you can use

``` r
plot(fit) 
```

An even more detailed investigation can be achieved by applying the shinystan package:

``` r
launch_shiny(fit) 
```

There are several methods to compute and visualize model predictions. Suppose that we want to predict responses (i.e. seizure counts) of a person in the treatment group (`Trt_c = 0.5`) and in the control group (`Trt_c = -0.5`) with average age and average number of previous seizures. Than we can use

``` r
newdata <- data.frame(Trt_c = c(0.5, -0.5), log_Age_c = 0, log_Base4_c = 0)
predict(fit, newdata = newdata, allow_new_levels = TRUE, probs = c(0.05, 0.95))
#>   Estimate Est.Error 5%ile 95%ile
#> 1  4.94375  4.108930     0     13
#> 2  6.83875  5.450842     1     17
```

We need to set `allow_new_levels = TRUE` because we want to predict responses of a person that was not present in the data used to fit the model. While the `predict` method returns predictions of the responses, the `fitted` method returns predictions of the regression line.

``` r
fitted(fit, newdata = newdata, allow_new_levels = TRUE, probs = c(0.05, 0.95))
#>   Estimate Est.Error    5%ile   95%ile
#> 1  5.00792  3.568603 1.408698 11.79776
#> 2  6.95871  4.881711 1.969220 16.04973
```

Both methods return the same etimate (up to random error), while the latter has smaller variance, because the uncertainty in the regression line is smaller than the uncertainty in each response. If we want to predict values of the original data, we can just leave the `newdata` argument empty.

A related feature is the computation and visualization of marginal effects, which can help in better understanding the influence of the predictors on the response.

``` r
plot(marginal_effects(fit, probs = c(0.05, 0.95)))
```

For a complete list of methods to apply on **brms** models see

``` r
methods(class = "brmsfit") 
#>  [1] add_ic                  as.array                as.data.frame          
#>  [4] as.matrix               as.mcmc                 coef                   
#>  [7] control_params          expose_functions        family                 
#> [10] fitted                  fixef                   formula                
#> [13] hypothesis              launch_shiny            log_lik                
#> [16] log_posterior           logLik                  loo                    
#> [19] LOO                     loo_linpred             loo_predict            
#> [22] loo_predictive_interval marginal_effects        marginal_smooths       
#> [25] model.frame             neff_ratio              ngrps                  
#> [28] nobs                    nsamples                nuts_params            
#> [31] pairs                   parnames                plot                   
#> [34] posterior_predict       posterior_samples       pp_check               
#> [37] pp_mixture              predict                 predictive_error       
#> [40] print                   prior_samples           prior_summary          
#> [43] ranef                   residuals               rhat                   
#> [46] stancode                standata                stanplot               
#> [49] summary                 update                  VarCorr                
#> [52] vcov                    waic                    WAIC                   
#> see '?methods' for accessing help and source code
```

Details on formula syntax, families and link functions, as well as prior distributions can be found on the help page of the brm function:

``` r
help("brm") 
```

More instructions on how to use **brms** are given in the package's main vignette.

``` r
vignette("brms_overview") 
```

FAQ
===

How do I install brms?
----------------------

To install the latest release version from CRAN use

``` r
install.packages("brms")
```

The current developmental version can be downloaded from github via

``` r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("paul-buerkner/brms", dependencies = TRUE)
```

Because **brms** is based on Stan, a C++ compiler is required. The program Rtools (available on <https://cran.r-project.org/bin/windows/Rtools/>) comes with a C++ compiler for Windows. On Mac, you should install Xcode. For further instructions on how to get the compilers running, see the prerequisites section on <https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started>.

What is the best way to ask a question or propose a new feature?
----------------------------------------------------------------

Questions can be asked in the google group [brms-users](https://groups.google.com/forum/#!forum/brms-users). To propose a new feature or report a bug, please open an issue on [github](https://github.com/paul-buerkner/brms). Of course, you can always write me an email (<paul.buerkner@gmail.com>).

How can I extract the generated Stan code?
------------------------------------------

If you have already fitted a model, just apply the `stancode` method on the fitted model object. If you just want to generate the Stan code without any model fitting, use the `make_stancode` function.

Can I avoid compiling models?
-----------------------------

When you fit your model for the first time with **brms**, there is currently no way to avoid compilation. However, if you have already fitted your model and want to run it again, for instance with more samples, you can do this without recompilation by using the `update` method. For more details see

``` r
help("update.brmsfit")
```

How can I specify non-linear or distributional models?
------------------------------------------------------

Specification of non-linear or distributional models requires multiple formulae. In **brms**, the function `brmsformula` (or short `bf`) is used to combine all formulae into one object, which can then be passed to the `formula` argument of `brm`. More help is given in

``` r
help("brmsformula")
```

For a detailed discussion of some examples see

``` r
vignette("brms_nonlinear")
```

``` r
vignette("brms_distreg")
```

What is the difference between brms and rstanarm?
-------------------------------------------------

**rstanarm** is an R package similar to **brms** that also allows to fit regression models using **Stan** for the backend estimation. Contrary to **brms**, **rstanarm** comes with precompiled code to save the compilation time (and the need for a C++ compiler) when fitting a model. However, as **brms** generates its **Stan** code on the fly, it offers much more flexibility in model specification than **rstanarm**. Also, multilevel models are currently fitted a bit more efficiently in **brms**. For a detailed comparison of **brms** with other common R packages implementing multilevel models, see

``` r
vignette("brms_overview")
```
