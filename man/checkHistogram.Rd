\name{checkHistogram}
\alias{checkHistogram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check any estimator by means of the multiscale confidence set
}
\description{
Give the locations (i.e. intervals) where the multiscale constraint is violated, and the change-points that are removable.
}
\usage{
checkHistogram(h, y, alpha = 0.1, q = NA, plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{
a numeric vector specifying values of a histogram at sample points or a \code{hitogram} class object (i.e. the return value of \code{hist}).
}
  \item{y}{
a numeric vector containing the data.
}
  \item{alpha}{
significance level; if \code{q} is missing, \code{q} is chosen as the (1-\code{alpha})-quantile of the null distribution of the multiscale statistic via Monte Carlo simulation, see also \code{\link{msQuantile}}.
}
  \item{q}{
threshold of the multiscale constraint.
}
  \item{plot}{
logical. If \code{TRUE}, the input estimator is potted,together with evaluation information. More precisely, at the very bottom, intervals where local constaints are violated are plotted. In the middle short vertical lines that indicate possibly removable change-points are drawn above a light blue horizontal line. Right below the light blue line, it plots a horizontal gray scale strap, the darkness of which reflects the number of violation intervals covering a given location, as a summary of violation information.
}
}
\value{
A data frame provides the intervals where the corresponding local side constraint is violated; an empty data frame if there is no violation. The data frame constains the following four columns
\item{leftIndex}{left index of an interval}
\item{rightIndex}{right index of an interval}
\item{leftEnd}{left end point of an interval}
\item{rightEnd}{right end point of an interval}
}
\references{
Li, H., Munk, A., Sieling, H., and Walther, G. (2016).The essential histogram. arXiv:1612.07216.
}
\note{
Only the intervals on which the input histogram is constant will be checked!
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{msQuantile}}, \code{\link{essHistogram}}
}
\examples{
set.seed(123)
# Data: mixture of Gaussians 1/3 N(0,0.5) + 1/3 N(5,1) + 1/3 N(15,2)
n  = 500
y  = rmixnorm(n, mean = c(0, 5, 15), sd = c(0.5, 1, 2))

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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}% use one of  RShowDoc("KEYWORDS")
