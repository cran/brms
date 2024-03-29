% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{epilepsy}
\alias{epilepsy}
\title{Epileptic seizure counts}
\format{
A data frame of 236 observations containing information
  on the following 9 variables.
\describe{
 \item{Age}{The age of the patients in years}
 \item{Base}{The seizure count at 8-weeks baseline}
 \item{Trt}{Either \code{0} or \code{1} indicating
   if the patient received anti-convulsant therapy}
 \item{patient}{The patient number}
 \item{visit}{The session number from \code{1} (first visit)
   to \code{4} (last visit)}
 \item{count}{The seizure count between two visits}
 \item{obs}{The observation number, that is
   a unique identifier for each observation}
 \item{zAge}{Standardized \code{Age}}
 \item{zBase}{Standardized \code{Base}}
}
}
\source{
Thall, P. F., & Vail, S. C. (1990).
 Some covariance models for longitudinal count data with overdispersion.
 \emph{Biometrics, 46(2)}, 657-671. \cr

Breslow, N. E., & Clayton, D. G. (1993).
 Approximate inference in generalized linear mixed models.
 \emph{Journal of the American Statistical Association}, 88(421), 9-25.
}
\usage{
epilepsy
}
\description{
Breslow and Clayton (1993) analyze data initially
  provided by Thall and Vail (1990) concerning
  seizure counts in a randomized trial of anti-convulsant
  therapy in epilepsy. Covariates are treatment,
  8-week baseline seizure counts, and age of the patients in years.
}
\examples{
\dontrun{
## poisson regression without random effects.
fit1 <- brm(count ~ zAge + zBase * Trt,
            data = epilepsy, family = poisson())
summary(fit1)
plot(fit1)

## poisson regression with varying intercepts of patients
## as well as normal priors for overall effects parameters.
fit2 <- brm(count ~ zAge + zBase * Trt + (1|patient),
            data = epilepsy, family = poisson(),
            prior = set_prior("normal(0,5)"))
summary(fit2)
plot(fit2)
}

}
\keyword{datasets}
