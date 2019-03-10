\name{essHist-package}
\alias{essHist-package}
\alias{essHist}
\docType{package}
\title{
\packageTitle{essHist}
}
\description{
\packageDescription{essHist}
}
\details{

% The DESCRIPTION file:
% \packageDESCRIPTION{essHist}
% \packageIndices{essHist}
\tabular{ll}{
Package: \tab essHist\cr
Type: \tab Package\cr
Version: \tab 1.2.0\cr
Date: \tab 2019-03-07\cr
License: \tab GPL-3\cr
}
Index:
\preformatted{
essHistogram    Compute the essential histogram
checkHistogram  Check any estimator by the multiscale confidence set

genIntv         Generate the system of intervals
msQuantile      Simulate the quantile of multiscale statistics

dmixnorm        Compute density function of Gaussian mixtures
pmixnorm        Compute distribution function of Gaussian mixtures
rmixnorm        Generate random number of Gaussian mixtures
paramExample    Output detailed parameters for some famous examples
}
}
\author{
\packageAuthor{essHist}

Maintainer: \packageMaintainer{essHist}
}
\references{
Li, H., Munk, A., Sieling, H., and Walther, G. (2016). The essential histogram. arXiv:1612.07216}
\keyword{package}
\keyword{distribution}
\keyword{nonparametric}
% \seealso{
% \code{\link[stepR]{smuceR}, \link[stepR]{jsmurf}}
% }
\examples{
# Simulate data
set.seed(123)
n = 300
y = rnorm(n)

# Compute the essential histogram
eh = essHistogram(y, plot = FALSE)

# Plot results
#     compute oracle density
x  = seq(min(y), max(y), length.out = n)
od = dnorm(x)
#     compare with orcle density
plot(x, od, type = "l", xlab = NA, ylab = NA, col = "red")
lines(eh)
legend("topright", c("Oracle density", "Essential histogram"), 
       lty = c(1,1), col = c("red", "black"))

###########################################################

# Evaluate other method e.g. R default histogram function
# Data: mixture of Gaussians 1/3 N(0,0.5) + 1/3 N(5,1) + 1/3 N(15,2)
set.seed(123)
n = 300
y = rmixnorm(n, mean = c(0, 5, 15), sd = c(0.5, 1, 2))

# Oracle density
sy = sort(y)
ho = dmixnorm(sy, mean = c(0, 5, 15), sd = c(0.5, 1, 2))

# R default histogram
h  = hist(y, plot = FALSE)

# Check R default histogram to local multiscale constriants
b = checkHistogram(h, y)
lines(sy, ho, col = "red")
legend("topright", c("R-Histogram", "Truth"), col = c("black", "red"), lty = c(1,1))
}


