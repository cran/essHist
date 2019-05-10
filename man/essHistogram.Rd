\name{Essential Histogram}
\alias{essHistogram}
\title{
The Essential Histogram
}
\description{
Compute the essential histogram via (pruned) dynamic programming.
}
\usage{
essHistogram(x, alpha = 0.5, q = NULL, intv = NULL, plot = TRUE, 
              mode = ifelse(anyDuplicated(x),"Gen","Con"), 
              xname = deparse(substitute(x)), ...)
}
\arguments{
  \item{x}{
a numeric vector containing the data.
}
  \item{alpha}{
  significance level; default as 0.5. One should set \code{alpha = 0.1} or even smaller if   confidence statements have to be made, while one can set \code{alpha = 0.9} if the goal is to explore the data for potential features with tolerance to false positives. The default value is only a trade-off.
}
  \item{q}{
  threshold value; by default, \code{q} is chosen as the (1-\code{alpha})-quantile of the null distribution of the multiscale statistic via Monte Carlo simulation, see also \code{\link{msQuantile}}.
}
\item{intv}{
  a data frame provides the system of intervals on which the multiscale statistic is defined. The data frame constains the following two columns   
  
  \code{left} left index of an interval
  
  \code{right} right index of an interval
  
  By default, it is set to the sparse interval system proposed by Rivera and Walther (2013), see also Li et al. (2016).
}
  \item{plot}{
  logical. If \code{TRUE} (default), a histogram is plotted, otherwise a list of breaks and counts is returned. In the latter case, a warning is used if (typically graphical) arguments are specified that only apply to the \code{plot = TRUE} case.
}
  \item{mode}{
  \code{"Con"} for continuous distribution functions
  
  \code{"Gen"} for general (possibly with discontinuous) distribution functions
  
  By default, \code{"Con"} is chosen if there is no tied observations; otherwise, \code{"Gen"} is chosen; see Li et al. (2016) for further details.
  }
  \item{xname}{
  a character string with the actual \code{x} argument name.
}
  \item{...}{
  further arguments and \code{\link{graphical parameters}} passed to \code{\link{plot.histogram}} and thence to \code{\link{title}} and \code{\link{axis}} (if \code{plot = TRUE}).
}
}
 \value{
  An object of class "\code{histogram}", which is of the same class as returned by function \code{\link{hist}}. 
}

\details{
  The essential histogram is defined as the histogram with least blocks within the multiscale constraint. The one with highest likelihood is picked if there are more than one solutions. The essential histogram involves only one parameter \code{q}, the threshold of the multiscale constraint. Such a parameter can be chosen by means of the significance level \code{alpha}, which leads to nature statistical significance statements for the multiscale constraint. The computational complexity is often linear in terms of sample size, although the worst complexity bound is quadratic up to a log-factor in case of the sparse interval system. See Li et al. (2016) for further details.
}
\references{
Li, H., Munk, A., Sieling, H., and Walther, G. (2016). The essential histogram. arXiv:1612.07216.

Rivera, C., & Walther, G. (2013). Optimal detection of a jump in the intensity of a {P}oisson process or in a density with likelihood ratio statistics. Scand. J. Stat. 40, 752--769.
}
% \author{
%%  ~~who you are~~
% }
\note{
The argument \code{intv} is internally adjusted to ensure it contains no empty intervals, especially in case of tied observations. The first block of the returned histogram is a closed interval, and the rest blocks are left open right closed intervals. All the printing messages can be disabled by calling \code{\link{suppressMessages}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{checkHistogram}},
\code{\link{genIntv}}, 
\code{\link{hist}}, 
\code{\link{msQuantile}}
}

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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}% use one of  RShowDoc("KEYWORDS")
\keyword{nonparametric}% __ONLY ONE__ keyword per line
