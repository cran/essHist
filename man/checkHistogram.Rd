\name{checkHistogram}
\alias{checkHistogram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check any histogram estimator by means of the multiscale confidence set
}
\description{
Provide the locations, i.e., intervals, where features are potentially missing (a.k.a. false negatives), and the break-points that are potentially redundant (a.k.a. false positives), by means of the multiscale confidence set. 
}
\usage{
checkHistogram(h, x, alpha = 0.1, q = NULL, intv = NULL, 
                mode = ifelse(anyDuplicated(x),"Gen","Con"), 
                plot = TRUE, xlim = NULL, ylim = NULL, 
                xlab = "", ylab = "", yaxt = "n", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{
a numeric vector specifying values of a histogram at sample points; or a \code{hitogram} class object (i.e. the return value of \code{\link{hist}}).
}
  \item{x}{
a numeric vector containing the data.
}
  \item{alpha}{
significance level, default as 0.1, see also \code{\link{essHistogram}}.
}
  \item{q}{
threshold of the multiscale constraint; by default, \code{q} is chosen as the (1-\code{alpha})-quantile of the null distribution of the multiscale statistic via Monte Carlo simulation, see also \code{\link{msQuantile}}.
}
\item{intv}{
  a data frame provides the system of intervals on which the multiscale statistic is defined. The data frame constains the following two columns   
  
  \code{left} left index of an interval
  
  \code{right} right index of an interval
  
  By default, it is set to the sparse interval system proposed by Rivera and Walther (2013), see also Li et al. (2016).
}
\item{mode}{
  \code{"Con"} for continuous distribution functions
  
  \code{"Gen"} for general (possibly with discontinuous) distribution functions
  
  By default, \code{"Con"} is chosen if there is no tied observations; otherwise, \code{"Gen"} is chosen; see Li et al. (2016) for further details.
  }
  \item{plot}{
logical. If \code{TRUE}, the input estimator is potted, together with evaluation information. More precisely, at the very bottom, intervals where local constaints are violated are plotted. In the middle short vertical lines that indicate possibly removable change-points are drawn above a light blue horizontal line. Right below the light blue line, it plots a horizontal gray scale strap, the darkness of which reflects the number of violation intervals covering a given location, as a summary of violation information.
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
A list consists of one data frame, and one numeric vector: 

\item{violatedIntervals}{A data frame provides the intervals where the corresponding local side constraint is violated; an empty data frame if there is no violation. It constains the following four columns

\code{leftIndex} left index of an interval

\code{rightIndex} right index of an interval

\code{leftEnd} left end point of an interval

\code{rightEnd} right end point of an interval

An empty \code{data.frame} is returned if there is no violation. 
}

\item{removableBreakpoints}{A numeric vector contains all removable breakpoints, with zero length if there is no removable breakpoint. }
}

\references{
Li, H., Munk, A., Sieling, H., and Walther, G. (2016). The essential histogram. arXiv:1612.07216.
}

\details{
This function presents a visualization: the upper part plots the given histogram; in the middle part short vertical lines mark all removable break-points; in the lower part intervals of violation are shown, and a graybar below the middle horizontal line (blue) sumarizes such violations with the darkness scaling with the number of violation intervals covering a location. See Examples below and Li et al. (2016) for further details. 
}

\note{
The argument \code{intv} is internally adjusted ensure it contains no empty intervals in case of tied observations. Only the intervals on which the input histogram is constant will be checked! All the printing messages can be disabled by calling \code{\link{suppressMessages}}.
}


\seealso{
\code{\link{essHistogram}},
\code{\link{genIntv}},
\code{\link{msQuantile}}
}
\examples{
set.seed(123)
# Data: mixture of Gaussians "harp"
n = 500
y = rmixnorm(n, type = 'harp')

# Oracle density
x = sort(y)
ho = dmixnorm(x, type = 'harp')

# R default histogram
h  = hist(y, plot = FALSE)

# Check R default histogram to local multiscale constriants
b = checkHistogram(h, y, ylim=c(-0.1,0.16))
lines(x, ho, col = "red")
rug(x, col = 'blue')
legend("topright", c("R-Histogram", "Truth"), col = c("black", "red"), lty = c(1,1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}% use one of  RShowDoc("KEYWORDS")
