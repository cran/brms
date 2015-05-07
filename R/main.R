#' Fit Bayesian Generalized Linear Mixed Models
#' 
#' Fit a Bayesian generalized linear mixed model using Stan
#' 
#' @param formula An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. 
#'   The details of model specification are given under 'Details'.
#' @param data An optional data frame, list or environment  (or object coercible by \code{as.data.frame} to a data frame) containing 
#'  the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}, 
#'  typically the environment from which \code{brm} is called.
#' @param family A vector of one or two character strings. The first string indicates the distribution of the dependent variable (the 'family'). Currently, the following families are supported:
#'  \code{"gaussian"}, \code{"student"}, \code{"cauchy"}, \code{"poisson"}, \code{"binomial"}, \code{"categorical"}, \code{"gamma"}, \code{"exponential"}, 
#'  \code{"weibull"}, \code{"cumulative"}, \code{"cratio"}, \code{"sratio"}, and \code{"acat"}.
#'  The second string indicates the link function, which must be supported by the distribution of the dependent variable. 
#'  If not specified, default link functions are used. Further information is provided under 'Details'.
#' @param prior A named list of character strings specifing the prior distributions of the parameters. Further information
#'  is provided under 'Details'.
#' @param partial A one sided formula of the form \code{~ partial.effects} specifing the predictors that can vary between categories in non-cumulative ordinal models
#'  (i.e. in families \code{"cratio"}, \code{"sratio"}, or \code{"acat"}).
#' @param threshold A character string indicating the type of thresholds (i.e. intercepts) used in an ordinal model.
#'  \code{"flexible"} provides the standard unstructured thresholds and \code{"equidistant"} restricts the distance between consecutive thresholds to the same value.
#' @param post.pred A flag to indicate if posterior predictives of the dependent variables should be generated
#' @param fit An instance of S4 class \code{stanfit} derived from a previous fit; defaults to \code{NA}. If fit is not \code{NA}, the compiled model associated with the fitted
#'  result is re-used. Make sure that all other arguments of \code{brm} are the same as at the time at which the \code{stanfit} object was generated.
#' @param n.chains Number of Markov chains (default: 2)
#' @param n.iter Number of total iterations per chain (including burnin; default: 2000)
#' @param n.warmup A positive integer specifying number of warmup (aka burnin) iterations. This also specifies the number of iterations used for stepsize adaptation, 
#'   so warmup samples should not be used for inference. The number of warmup should not be larger than \code{n.iter} and the default is 500.
#' @param n.thin Thinning rate. Must be a positive integer. Set \code{n.thin > 1} to save memory and computation time if \code{n.iter} is large. Default is 1, that is no thinning.
#' @param n.cluster	Number of clusters to use to run parallel chains. Default is 1. 
#' @param pars Either "auto", "all", or a character vector containing the names of the parameters to be saved. If \code{"auto"}, the function \code{\link[brms:brm.pars]{brm.pars}} chooses
#'   the observed parameters.
#' @param inits A list with n.chains elements; each element of the list is itself a list of starting values for the model, or a function creating (possibly random) initial values. 
#'   If inits is \code{NULL} (the default), Stan will generate initial values for parameters. Other options supported by Stan are also possible.
#' @param save.model Either \code{NULL} or a character string. In the latter case, the model code is
#'   saved in a file with its name specified by \code{save.model} in the current working directory.
#' @param seed Positive integer. Used by \code{set.seed} to make results reproducable.  
#' @param engine A character string, either \code{"stan"} (the default) or \code{"jags"}. Specifies which program should be used to fit the model. 
#'  Note that \code{jags} is currently implemented for testing purposes only, does not allow full functionality and is not supported or documented.
#' @param ... Further arguments to be passed to Stan.
#' 
#' @return An object of class \code{stanfit}, which contains the posterior samples. If rstan is not installed,
#'  a named list containing the Stan model, the required data and the parameters of interest is returned instead.
#' 
#' @details Fit a generalized linear mixed model, which incorporates both fixed-effects parameters and random effects in a linear predictor  
#'   via full bayesian inference using Stan. \cr
#'   
#'   \bold{Formula syntax}
#'   
#'   The \code{formula} argument accepts formulas of the following syntax: 
#'   
#'   \code{response | addition ~ fixed + (random | group)} 
#'   
#'   Multiple grouping factors each with multiple random effects are possible. With the exception of \code{addition}, 
#'   this is basically \code{lme4} syntax.
#'   The optional argument \code{addition} has different meanings depending on the \code{family} argument. 
#'   
#'   For families \code{gaussian}, \code{student}, and \code{cauchy} 
#'   it may be a variable specifying the standard errors of the observation, thus allowing to perform meta-analysis. 
#'   Suppose that the variable \code{yi} contains the effect sizes from the studies and \code{sei} the 
#'   corresponding standard errors. Then, fixed and random effects meta-analyses can be conducted
#'   using the formulae \code{yi | sei ~ 1} and \code{yi | sei ~ 1 + (1|study)}, respectively, where 
#'   \code{study} is a variable uniquely identifying every study.
#'   If desired, meta-regressen can be performed via \code{yi | sei ~ 1 + mod1 + mod2 + (1|study)} 
#'   or \code{yi | sei ~ 1 + mod1 + mod2 + (1 + mod1 + mod2|study)}, where
#'   \code{mod1} and \code{mod2} represent moderator variables.
#'   
#'   For family \code{binomial}, addition may be a variable indicating the number of trials 
#'   underlying each observation. In \code{lme4} syntax, we may write for instance 
#'   \code{cbind(success, trials - success)}, which is equivalent
#'   to \code{success | trials} in \code{brms} syntax. If the number of trials
#'   is constant across all observation (say \code{10}), we may also write \code{success | 10}. 
#'   
#'   For family \code{categorical} and all ordinal families, \code{addition} specifies the number of 
#'   categories for each observation, either with a variable name or a single number.
#'   
#'   For families \code{gamma}, \code{exponential}, and \code{weibull}, \code{addition} may contain 
#'   a logical variable (or a variable than can be coerced to logical) indicating
#'   if the response variable is left censored (corresponding to \code{TRUE}) or not censored 
#'   (corresponding to \code{FALSE}). \cr
#' 
#'   \bold{Families and link functions}
#'   
#'   Family \code{gaussian} with \code{identity} link leads to linear regression. Families \code{student}, and \code{cauchy}
#'   with \code{identity} link leads to robust linear regression that is less influenced by outliers. 
#'   Family \code{poisson} with \code{log} link leads to poisson regression for count data. 
#'   Family \code{binomial} with \code{logit} link leads to logistic regression and family \code{categorical} to
#'   multi-logistic regression when there are more than two possible outcomes.
#'   Families \code{cumulative}, \code{cratio} ('contiuation ratio'), \code{sratio} ('stopping ratio'), 
#'   and \code{acat} ('adjacent category') leads to ordinal regression. Families \code{gamma}, \code{weibull}, and \code{exponential}
#'   can be used (among others) for survival regression when combined with the \code{log} link.
#'   
#'   In the following, we list all possible links for each family.
#'   The families \code{gaussian}, \code{student}, and \code{cauchy} accept the links (as names) \code{identity}, \code{log}, and \code{inverse};
#'   the \code{poisson} family the links \code{log}, \code{identity}, and \code{sqrt}; 
#'   families \code{binomial}, \code{cumulative}, \code{cratio}, \code{sratio}, and \code{acat} the links \code{logit}, \code{probit}, \code{probit_approx}, and \code{cloglog};
#'   family  \code{categorical} the link \code{logit}; families \code{gamma}, \code{weibull}, and \code{exponential} the links \code{log}, \code{identity}, and \code{inverse}. 
#'   The first link mentioned for each family is the default. \cr    
#'   
#'   
#'   
#'   \bold{Prior distributions}
#'   
#'   Below, we describe the usage of the \code{prior} argument and list some common prior distributions 
#'   for parameters in \code{brms} models. 
#'   A complete overview on possible prior distributions is given in the Stan Reference Manual available at 
#'   \url{http://mc-stan.org/}.
#'   
#'   \code{brm} performs no checks if the priors are written in correct Stan language.
#'   Instead, Stan will check their correctness when the model is parsed to C++ and returns an error if they are not.
#'   Currently, there are four types of parameters in \code{brms} models, 
#'   for which the user can specify prior distributions. \cr
#'   
#'   1. Fixed effects 
#'   
#'   Every fixed (and partial) effect has its corresponding regression parameter. These parameters are named as
#'   \code{b_(fixed)}, where \code{(fixed)} represents the name of the corresponding fixed effect. 
#'   Suppose, for instance, that \code{y} is predicted by \code{x1} and \code{x2} 
#'   (i.e. \code{y ~ x1+x2} in formula syntax). 
#'   Then, \code{x1} and \code{x2} have regression parameters \code{b_x1} and \code{b_x2} respectively. 
#'   The default prior for fixed effects parameters is an improper flat prior over the reals. 
#'   Other common options are normal priors or uniform priors over a finite interval.
#'   If we want to have a normal prior with mean 0 and standard deviation 5 for \code{b_x1}, 
#'   and a uniform prior between -10 and 10 for \code{b_x2},
#'   we can specify this via \cr
#'   \code{prior = list(b_x1 = "normal(0,5)", b_x2 = "uniform(-10,10)")}. 
#'   To put the same prior (e.g. a normal prior) on all fixed effects at once, 
#'   we may write as a shortcut \code{prior = } \cr \code{list(b = "normal(0,5)")}. In addition, this
#'   leads to faster sampling in Stan, because priors can be vectorized. \cr
#'   
#'   2. Standard deviations of random effects
#'   
#'   Each random effect of each grouping factor has a standard deviation named
#'   \code{sd_(group)_(random)}. Consider, for instance, the formula \code{y ~ x1+x2+(1+x1|z)}.
#'   We see that the intercept as well as \code{x1} are random effects nested in the grouping factor \code{z}. 
#'   The corresponding standard deviation parameters are named as \code{sd_z_Intercept} and \code{sd_z_x1} respectively. 
#'   These parameters are restriced to be non-negative and, by default, 
#'   have a half cauchy prior with 'mean' 0 and 'standard deviation' 5. 
#'   We could make this explicit by writing \code{prior = list(sd = "cauchy(0,5)")}. 
#'   One common alternative is a uniform prior over a positive interval. \cr
#'   
#'   3. Correlations of random effects 
#'   
#'   If there is more than one random effect per grouping factor, the correlations between those random
#'   effects have to be estimated. 
#'   However, in \code{brms} models, the corresponding correlation matrix \eqn{C} does not have prior itself. 
#'   Instead, a prior is defined for the cholesky factor \eqn{L} of \eqn{C}. They are related through the equation
#'     \deqn{L * L' = C.} 
#'   The prior \code{"lkj_corr_cholesky(eta)"} with \code{eta > 0} is essentially the only prior for 
#'   cholesky factors of correlation matrices.
#'   If \code{eta = 1} (the default) all correlations matrices are equally likely a priori. If \code{eta > 1}, 
#'   extreme correlations become less likely, 
#'   whereas \code{0 < eta < 1} results in higher probabilities for extreme correlations. 
#'   The cholesky factors in \code{brms} models are named as 
#'   \code{L_(group)}, (e.g., \code{L_z} if \code{z} is the grouping factor). \cr
#'   
#'   4. Parameters for specific families 
#'   
#'   Some families need additional parameters to be estimated. 
#'   Families \code{gaussian}, \code{student}, and \code{cauchy} need the parameter \code{sigma} 
#'   to account for the standard deviation of the response variable around the regression line
#'   (not to be confused with the standard deviations of random effects). 
#'   By default, \code{sigma} has an improper flat prior over the positiv reals. 
#'   Furthermore, family \code{student} needs the parameter \code{nu} representing 
#'   the degrees of freedom of students t distribution. 
#'   By default, \code{nu} has prior \code{"uniform(1,60)"}. 
#'   Families \code{gamma} and \code{weibull} need the parameter \code{shape} 
#'   that has a \code{"gamma(0.01,0.01)"} prior by default. \cr
#'   
#'   \bold{Parameters of interest}
#'   
#'   If \code{pars = "auto"} (the default) only certain parameters are returned by \code{brm}. 
#'   These are the fixed (and partial) regression parameters, 
#'   the random effects standard deviations and correlations, as well as parameters specific to certain 
#'   families (see also section 'Prior distributions'). By default, the random effects themselves (named as \code{r_(group)}) are not returned, 
#'   and one has to use \code{pars = "reffects"} to add them to the set of returned parameters.
#'   When \code{post.pred = }\code{TRUE}, posterior predictive samples are also included.
#' 
#' @examples
#' \dontrun{ 
#' ### Poisson Regression for the number of seizures in epileptic patients
#' ### using half cauchy priors for standard deviations of random effects 
#' fit_e <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient) + (1|visit), 
#'            data = epilepsy, family = c("poisson", "log"), 
#'            prior = list(sd = "cauchy(0,2.5)"))
#' brm.plot(fit_e)   
#' print(fit_e) 
#' 
#' ### Ordinal regression (with family 'sratio') modeling patient's rating 
#' ### of inhaler instructions using normal priors for fixed effects parameters
#' fit_i <- brm(rating ~ treat + period + carry, data = inhaler, 
#'               family = "sratio", prior = prior = list(b = "normal(0,5)"))
#' brm.plot(fit_i)
#' print(fit_i)    
#' 
#' ### Surivival Regression (with family 'weibull') modeling time between 
#' ### first and second recurrence of an infection in kidney patients
#' ### time | cens indicates which values in variable time are left censored
#' fit_k <- brm(time | cens ~ age + sex + disease, data = kidney, family = "weibull")
#' brm.plot(fit_k) 
#' print(fit_k)             
#' }
#' 
#' @import parallel
#' @import Rcpp
#' @export 
brm <- function(formula, data = NULL, family = c("gaussian", "identity"), prior = list(),
                partial = NULL, threshold = "flexible", post.pred = FALSE,  fit = NA, 
                n.chains = 2, n.iter = 2000, n.warmup = 500, n.thin = 1, n.cluster = 1,
                pars = "auto", inits = "random", save.model = NULL, 
                seed = 12345, engine = "stan", ...) {
  dots <- list(...)  
  if (n.chains %% n.cluster != 0) stop("n.chains must be a multiple of n.cluster")
  if (engine %in% c("stan","jags")) stan <- engine == "stan"
  else stop("engine must be either stan or jags")
  if (!is.element(threshold,c("flexible","equidistant"))) 
    stop("threshold must be either flexible or equidistant")
  
  set.seed(seed)
  supl.data <- brm.data(formula, data = data, family = family, prior = prior, 
                        partial = partial, engine = engine, ...) 
  if (is.function(inits) | (is.character(inits) & !is.element(inits, c("random", "0")))) 
    inits <- replicate(n.chains, do.call(inits, list()), simplify = FALSE)
  if (pars %in% c("auto","reffects")) 
    pars <- brm.pars(formula, data, family = family[1], partial = partial, threshold = threshold,
      reffects = pars == "reffects", engine = engine, post.pred = post.pred)
  else if (pars == "all") pars <- NA
  else if (!is.vector(pars)) 
    stop("Argument pars must be either 'auto', 'reffects', 'all', or a vector of parameter names")
  
  if (stan) {
    if (is(fit, "stanfit")) 
      warning(paste("When passing a stanfit object to argument 'fit', make sure that all other arguments are the same", 
              "as at the time at which the stanfit object was generated.")) 
    model <- stan.model(formula = formula, data = data, family = family[1], link = brm.link(family), 
                        prior = prior, partial = partial, post.pred = post.pred, 
                        threshold = threshold, save.model = save.model)
    if (!requireNamespace("rstan", quietly = TRUE)) {
      warning(paste("Package rstan is not installed yet so that the model cannot be fitted.
        Returning the Stan model, the required data, and the parameters of interest instead.
        Please see https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
        for instructions on how to install rstan.")) 
      return(list(model = model, data = supl.data, pars = pars))
    }

    if (n.cluster > 1) {
      if (is.character(inits) | is.numeric(inits)) inits <- rep(inits, n.chains)
      cl <- makeCluster(n.cluster)
      clusterEvalQ(cl, require(rstan))
      clusterExport(cl = cl, c("fit", "supl.data", "pars", "inits", "n.iter",
                               "n.warmup", "n.thin"), envir = environment())
      fit <- rstan::stan(model_code = model, data = supl.data, chains = 0, fit = fit, ...)
      sflist <- parLapply(cl, 1:n.chains, fun = function(i)  
        rstan::stan(fit = fit, data = supl.data, iter = n.iter, pars = pars, init = inits[i],
                    warmup = n.warmup, thin = n.thin, chains = 1, chain_id = i))
      fit <- rstan::sflist2stanfit(sflist)
      stopCluster(cl)
    } 
    else fit <- rstan::stan(model_code = model, data = supl.data, pars = pars, init = inits, 
                           iter = n.iter, chains = n.chains, warmup = n.warmup, thin = n.thin, 
                           fit = fit, ...)
  } 
  else {
    warning("Engine 'jags' is currently implemented for testing purposes only and we do not support its usage.")
    if (is.character(inits) | is.numeric(inits)) 
      inits <- replicate(n.chains, bugs.inits(formula = formula, data = data, family = family[1], 
        partial = partial, threshold = threshold, engine = engine, range = dots$range), simplify = FALSE)
    model <- brm.bugs(formula = formula, data = data, family = family[1], link = brm.link(family), 
                      prior = prior, partial = partial, threshold = threshold, 
                      post.pred = post.pred, save.model = save.model)
    if (!requireNamespace("R2jags", quietly = TRUE)) {
      warning(paste0("Package 'R2jags' is not installed yet so that the model cannot be fitted.
        Returning the Bugs model, the required data, and the parameters of interest instead."))
      return(list(model = model, data = supl.data, pars = pars))
    }  
    fit <- suppressWarnings(R2jags::jags(supl.data, inits, pars, textConnection(model), 
             n.chains = n.chains, n.iter = n.iter, n.burnin = n.warmup, n.thin = n.thin))
  }
  fit
}