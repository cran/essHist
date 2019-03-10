\name{checkHistogram}
\alias{checkHistogram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check any histogram estimator by means of the multiscale confidence set
}
\description{
Give the locations (i.e. intervals) where the multiscale constraint is violated, and the change-points that are removable.
}
\usage{
checkHistogram(h, y, alpha = 0.1, q = NA, intv = genIntv(length(y)), 
               plot = TRUE, verbose = TRUE, xlim = range(y), 
               ylim = NULL, xlab = "", ylab = "", yaxt = "n", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{
a numeric vector specifying values of a histogram at sample points or a \code{hitogram} class object (i.e. the return value of \code{\link{hist}}).
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
\item{intv}{
  a data frame provides the system of intervals on which the multiscale statistic is defined. The data frame constains the following two columns   
  
  \code{left} left index of an interval
  
  \code{right} right index of an interval
  
  By default, it is set to the sparse interval system proposed by Rivera and Walther (2013), see also Li et al. (2016).
}
  \item{plot}{
logical. If \code{TRUE}, the input estimator is potted, together with evaluation information. More precisely, at the very bottom, intervals where local constaints are violated are plotted. In the middle short vertical lines that indicate possibly removable change-points are drawn above a light blue horizontal line. Right below the light blue line, it plots a horizontal gray scale strap, the darkness of which reflects the number of violation intervals covering a given location, as a summary of violation information.
}
\item{verbose}{
logical. If \code{TRUE} (default) it prints some details about the computation; otherwise nothing is printed.   
}
\item{xlim, ylim}{
numeric vectors of length 2 (default \code{xlim = range(y), ylim = NULL}): see \code{\link{plot}}.
}
\item{xlab}{
a title for the \code{x} axis (default empty string): see \code{\link{title}} and \code{\link{plot}}.
}
\item{ylab}{
a title for the \code{y} axis (default empty string): see \code{\link{title}} and \code{\link{plot}}.
}
\item{yaxt}{
A character which specifies the \code{y} axis type (default \code{"n"}): see \code{\link{par}}. 
}
\item{...}{
further arguments and \code{\link{graphical parameters}} passed to \code{\link{plot}} (if \code{plot} = \code{TRUE}).
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
Li, H., Munk, A., Sieling, H., and Walther, G. (2016). The essential histogram. arXiv:1612.07216.
}

\note{
The argument \code{intv} is internally adjusted via function \code{.validInterval} to ensure it is reasonable, and contains no empty intervals in case of tied observations. Only the intervals on which the input histogram is constant will be checked! 
}


\seealso{
\code{\link{essHistogram}},
\code{\link{genIntv}},
\code{\link{msQuantile}}
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
