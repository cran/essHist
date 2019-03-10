\name{Generate Intervals}
\alias{genIntv}
\title{
Generate the system of intervals
}
\description{
Generate the system of intervals on which the multiscale statistic is defined, see Li et al. (2016).
}
\usage{
genIntv(n, type = c("Sparse", "Full"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
  number of observations.
}
  \item{type}{
  type of interval system. \code{type = "Sparse"} (default) is the sparse system proposed by Rivera and Walther (2013), see also Li et al. (2016). \code{type = "Full"} is the system of all possible intervals with end-index ranging from \code{1} to \code{n}.   
}
}
\value{
A data frame provides the system of intervals, and consists two columns
\item{left}{left index of an interval}
\item{right}{right index of an interval}
}
\references{
Li, H., Munk, A., Sieling, H., and Walther, G. (2016). The essential histogram. arXiv:1612.07216.


Rivera, C., & Walther, G. (2013). Optimal detection of a jump in the intensity of a {P}oisson process or in a density with likelihood ratio statistics. Scand. J. Stat. 40, 752--769.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{checkHistogram}},
\code{\link{essHistogram}},
\code{\link{msQuantile}}
}
\examples{
n    = 5
intv = genIntv(n,"Full")
print(intv)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}% use one of  RShowDoc("KEYWORDS")
