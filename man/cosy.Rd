% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formula-ac.R
\name{cosy}
\alias{cosy}
\title{Set up COSY correlation structures}
\usage{
cosy(time = NA, gr = NA)
}
\arguments{
\item{time}{An optional time variable specifying the time ordering
of the observations. By default, the existing order of the observations
in the data is used.}

\item{gr}{An optional grouping variable. If specified, the correlation
structure is assumed to apply only to observations within the same grouping
level.}
}
\value{
An object of class \code{'cosy_term'}, which is a list
  of arguments to be interpreted by the formula
  parsing functions of \pkg{brms}.
}
\description{
Set up a compounds symmetry (COSY) term in \pkg{brms}. The function does
not evaluate its arguments -- it exists purely to help set up a model with
COSY terms.
}
\examples{
\dontrun{
data("lh")
lh <- as.data.frame(lh)
fit <- brm(x ~ cosy(), data = lh)
summary(fit)
}

}
\seealso{
\code{\link{autocor-terms}}
}
