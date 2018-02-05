\name{msQuantile}
\alias{msQuantile}
\title{
Quantile of the multiscale statistics
}
\description{
Simulate quantiles of the multiscale statistics under the null hypothesis.
}
\usage{
msQuantile(n, alpha = c(0.1), nsim = 5000)
}
\arguments{
  \item{n}{
  number of observations
}
  \item{alpha}{
  significance level; the (1-\code{alpha})-quantile of the null distribution of the multiscale statistic via Monte Carlo simulation
}
  \item{nsim}{
numer of Monte Carlo simulations
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A scalar value with attributes \code{names} storing (1 - \code{alpha}). See Li et al. (2016) for further details.
}
\references{
Li, H., Munk, A., Sieling, H., and Walther, G. (2016).The essential histogram. arXiv:1612.07216.}
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
n    = 100 # number of observations
nsim = 100 # number of simulations

alpha = c(0.1, 0.9) # significance level
q     = msQuantile(n, alpha, nsim)

print(q)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{nonparametric}
