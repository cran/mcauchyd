% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotmcd.R
\name{plotmcd}
\alias{plotmcd}
\title{Plot of the Bivariate Cauchy Density}
\usage{
plotmcd(mu, Sigma, xlim = c(mu[1] + c(-10, 10)*Sigma[1, 1]),
                ylim = c(mu[2] + c(-10, 10)*Sigma[2, 2]), n = 101,
                xvals = NULL, yvals = NULL, xlab = "x", ylab = "y",
                zlab = "f(x,y)", col = "gray", tol = 1e-6, ...)
}
\arguments{
\item{mu}{length 2 numeric vector.}

\item{Sigma}{symmetric, positive-definite square matrix of order 2. The scatter matrix.}

\item{xlim, ylim}{x-and y- limits.}

\item{n}{A one or two element vector giving the number of steps in the x and y grid, passed to \code{\link[rgl]{plot3d.function}}.}

\item{xvals, yvals}{The values at which to evaluate \code{x} and \code{y}. If used, \code{xlim} and/or \code{ylim} are ignored.}

\item{xlab, ylab, zlab}{The axis labels.}

\item{col}{The color to use for the plot. See \code{\link[rgl]{plot3d.function}}.}

\item{tol}{tolerance (relative to largest variance) for numerical lack of positive-definiteness in Sigma, for the estimation of the density. see \code{\link{dmcd}}.}

\item{...}{Additional arguments to pass to \code{\link[rgl]{plot3d.function}}.}
}
\value{
Returns invisibly the probability density function.
}
\description{
Plots the probability density of the multivariate Cauchy distribution with 2 variables
with location parameter \code{mu} and scatter matrix \code{Sigma}.
}
\examples{
mu <- c(1, 4)
Sigma <- matrix(c(0.8, 0.2, 0.2, 0.2), nrow = 2)
plotmcd(mu, Sigma)

}
\references{
N. Bouhlel, D. Rousseau, A Generic Formula and Some Special Cases for the Kullback–Leibler Divergence between Central Multivariate Cauchy Distributions.
Entropy, 24, 838, July 2022.
\doi{10.3390/e24060838}
}
\seealso{
\code{\link{dmcd}}: probability density of a multivariate Cauchy density

\code{\link{contourmcd}}: contour plot of a bivariate Cauchy density.

\code{\link[rgl]{plot3d.function}}: plot a function of two variables.
}
\author{
Pierre Santagostini, Nizar Bouhlel
}
