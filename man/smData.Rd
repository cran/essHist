\name{smData}
\alias{smData}
\alias{smoothData}
\alias{resolveTies}
\alias{breakTies}
\title{
Smooth discrete data
}
\description{
Resolve equal values (a.k.a., tied values) of observations by randomizing locally according to a Gaussian distribution with a rather small variance. 
}
\usage{
smData(x, sw = NA)
}
\arguments{
  \item{x}{
  a numeric vector containing the data.
}
  \item{sw}{
  a positive number; it specifies the spread width of the randomization procedure; its default value is from the minimal gap between two different values of observations.
}
}
\details{
The essential histogram (Li et al, 2016) is designed based on the assumption that the underlying distribution function is continuous. Such assumption is natural as it guarantees the existence of density with respect to the Lebesgue measure. However, in pratice, one also faces discrete distributions, whose distribution function is piece-wise constant, thus discontinuous. The function \code{smData} implements a simple idea of adapting the essential histogram to discrete data: more precisely, the Dirac delta density is approximated by a thin Guassian density, and the resulted approximation has continuous distribution. 

The function \code{smData} is automatically called, when \code{essHistogram} is called. Note that \code{smData} only sorts the observations \code{x} if there is no tied values. 
}
\value{
A vector of length \code{length(x)} is returned, i.e., modified observations with no tied values, and ordered increasingly. 
}
\references{
Li, H., Munk, A., Sieling, H., and Walther, G. (2016). The essential histogram. arXiv:1612.07216.}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{essHistogram}}
}
\examples{
# generate Poisson data (discrete)
set.seed(123)
n      = 100 # number of observations
lambda = 5
x.dis  = rpois(n, lambda)

# smooth discrete data
x.sm   = smData(x.dis)

# compute the essential histogram
eh = essHistogram(x.sm, xname = "Poisson")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
