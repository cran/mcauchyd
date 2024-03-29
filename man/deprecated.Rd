% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{deprecated}
\alias{deprecated}
\alias{mvdcd}
\alias{mvrcd}
\alias{plotmvcd}
\alias{contourmvcd}
\title{Deprecated functions}
\usage{
mvdcd(x, mu, Sigma, tol = 1e-06)

mvrcd(n, mu, Sigma, tol = 1e-06)

plotmvcd(
  mu,
  Sigma,
  xlim = c(mu[1] + c(-10, 10) * Sigma[1, 1]),
  ylim = c(mu[2] + c(-10, 10) * Sigma[2, 2]),
  n = 101,
  xvals = NULL,
  yvals = NULL,
  xlab = "x",
  ylab = "y",
  zlab = "f(x,y)",
  col = "gray",
  tol = 1e-06,
  ...
)

contourmvcd(
  mu,
  Sigma,
  beta,
  xlim = c(mu[1] + c(-10, 10) * Sigma[1, 1]),
  ylim = c(mu[2] + c(-10, 10) * Sigma[2, 2]),
  zlim = NULL,
  npt = 30,
  nx = npt,
  ny = npt,
  main = "Multivariate generalised Gaussian density",
  sub = NULL,
  nlevels = 10,
  levels = pretty(zlim, nlevels),
  tol = 1e-06,
  ...
)
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}

Use \code{\link[=dmcd]{dmcd()}} instead of \code{mvdcd()}.

Use \code{\link[=rmcd]{rmcd()}} instead of \code{mvrcd()}.

Use \code{\link[=plotmcd]{plotmcd()}} instead of \code{plotmvcd()}.
Use \code{\link[=contourmcd]{contourmcd()}} instead of \code{contourmvcd()}.

Use \code{\link[=contourmcd]{contourmcd()}} instead of \code{contourmvcd()}.
}
\keyword{internal}
