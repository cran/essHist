\name{essHistogram}
\alias{essHistogram}
\title{
The Essential Histogram
}
\description{
Compute the essential histogram via dynamic programming.
}
\usage{
essHistogram(x, alpha = 0.5, q = NA, plot = TRUE, xname = deparse(substitute(x)), ...)
}
\arguments{
  \item{x}{
a numeric vector containing the data.
}
  \item{alpha}{
significance level; if \code{q} is missing, \code{q} is chosen as the (1-\code{alpha})-quantile of the null distribution of the multiscale statistic via Monte Carlo simulation, see also \code{\link{msQuantile}}.
}
  \item{q}{
threshold value.
}
  \item{plot}{
logical. If \code{TRUE} (default), a histogram is plotted, otherwise a list of breaks and counts is returned. In the latter case, a warning is used if (typically graphical) arguments are specified that only apply to the plot = \code{TRUE} case.
}
  \item{xname}{
a character string with the actual \code{x} argument name.
}
  \item{...}{
further arguments and \code{\link{graphical parameters}} passed to \code{\link{plot.histogram}} and thence to \code{\link{title}} and \code{\link{axis}} (if \code{plot} = \code{TRUE}).
}
}
 \value{
 An object of class "\code{histogram}", which is the same as \code{\link{hist}}. 
}

\details{
The essential histogram is defined as the histogram with least blocks within the multiscale constraint. The one with highest likelihood is picked if there are more than one solutions. The essential histogram involves only one parameter \code{q}, the threshold of the multiscale constraint. Such a parameter can be chosen by means of the significance level \code{alpha}, which leads to nature statistical significance statements for the multiscale constraint. See Li et al. (2016) for further details.
}
\references{
Li, H., Munk, A., Sieling, H., and Walther, G. (2016).The essential histogram. arXiv:1612.07216.}
% \author{
%%  ~~who you are~~
% }
% \note{
%%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{msQuantile}}, \code{\link{hist}}
}
\examples{
# simulate data
set.seed(123)
n = 300
x = rnorm(n)

# compute the essential histogram
eh = essHistogram(x, xname = "Gauss")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}% use one of  RShowDoc("KEYWORDS")
\keyword{nonparametric}% __ONLY ONE__ keyword per line
