params <-
list(EVAL = TRUE)

## ----SETTINGS-knitr, include=FALSE------------------------------------------------------
stopifnot(require(knitr))
options(width = 90)
knit_hooks$set(pngquant = knitr::hook_pngquant)
opts_chunk$set(
  comment = NA,
  message = FALSE,
  warning = FALSE,
  eval = if (isTRUE(exists("params"))) params$EVAL else FALSE,
  dev = "ragg_png",
  dpi = 72,
  fig.retina = 1.5,
  fig.asp = 0.8,
  fig.width = 5,
  out.width = "60%",
  fig.align = "center",
  pngquant = "--speed=1 --quality=50"
)
library(brms)
ggplot2::theme_set(theme_default())

## ---------------------------------------------------------------------------------------
data("nhanes", package = "mice")
head(nhanes)

## ---------------------------------------------------------------------------------------
library(mice)
m <- 5
imp <- mice(nhanes, m = m, print = FALSE)

## ----results = 'hide', message = FALSE--------------------------------------------------
fit_imp1 <- brm_multiple(bmi ~ age*chl, data = imp, chains = 2)

## ---------------------------------------------------------------------------------------
summary(fit_imp1)

## ---------------------------------------------------------------------------------------
plot(fit_imp1, variable = "^b", regex = TRUE)

## ---------------------------------------------------------------------------------------
library(posterior)
draws <- as_draws_array(fit_imp1)
# every dataset has nc = 2 chains in this example
nc <- nchains(fit_imp1) / m
draws_per_dat <- lapply(1:m, 
  \(i) subset_draws(draws, chain = ((i-1)*nc+1):(i*nc))
)
lapply(draws_per_dat, summarise_draws, default_convergence_measures())

## ---------------------------------------------------------------------------------------
conditional_effects(fit_imp1, "age:chl")

## ----results = 'hide', message = FALSE--------------------------------------------------
bform <- bf(bmi | mi() ~ age * mi(chl)) +
  bf(chl | mi() ~ age) + set_rescor(FALSE)
fit_imp2 <- brm(bform, data = nhanes)

## ---------------------------------------------------------------------------------------
summary(fit_imp2)
conditional_effects(fit_imp2, "age:chl", resp = "bmi")

## ---------------------------------------------------------------------------------------
nhanes$se <- rexp(nrow(nhanes), 2)

## ----results = 'hide', message = FALSE, eval = FALSE------------------------------------
# bform <- bf(bmi | mi() ~ age * mi(chl)) +
#   bf(chl | mi(se) ~ age) + set_rescor(FALSE)
# fit_imp3 <- brm(bform, data = nhanes)

