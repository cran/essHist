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
Version: \tab 1.2.2\cr
Date: \tab 2019-05-10\cr
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
type = 'skewed_unimodal'
n = 500
y = rmixnorm(n, type = type)

# Compute the essential histogram
eh = essHistogram(y, plot = FALSE)

# Plot results
#     compute oracle density
x  = sort(y)
od = dmixnorm(x, type = type)
#     compare with orcle density
plot(x, od, type = "l", xlab = NA, ylab = NA, col = "red", main = type)
lines(eh)
legend("topleft", c("Oracle density", "Essential histogram"), 
       lty = c(1,1), col = c("red", "black"))



##### Evaluate other method
set.seed(123)
# Data: mixture of Gaussians "harp"
n = 500
y = rmixnorm(n, type = 'harp')

# Oracle density
x = sort(y)
ho = dmixnorm(x, type = 'harp')

# R default histogram
h  = hist(y, plot = FALSE)

# Check R default histogram to local multiscale constraints
b = checkHistogram(h, y, ylim=c(-0.1,0.16))
lines(x, ho, col = "red")
rug(x, col = 'blue')
legend("topright", c("R-Histogram", "Truth"), col = c("black", "red"), lty = c(1,1))
}


