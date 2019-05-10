\name{Multiscale Quantiles}
\alias{msQuantile}
\title{
Quantile of the multiscale statistics
}
\description{
Simulate quantiles of the multiscale statistics under any continuous distribution function.
}
\usage{
msQuantile(n, alpha = c(0.5), nsim = 5e3, is.sim = (n < 1e4), 
            intv = genIntv(n), mode = c("Con", "Gen"), ...)
}
\arguments{
  \item{n}{
  number of observations.
  }
  \item{alpha}{
  significance level; default as 0.5, see also \code{\link{essHistogram}}. Like \code{\link{quantile}}, it can also be a vector. 
}
  \item{nsim}{
  numer of Monte Carlo simulations.
} 
  \item{is.sim}{
  logical. If \code{TRUE} (default if \code{n} < 10,000) the quantile is determined via Monte Carlo simulations, which might take a long time; otherwise (default if \code{n} >= 10,000) it uses the quantile with \code{n} = 10,000, which has been precomputed and stored.
}
  \item{intv}{
  a data frame provides the system of intervals on which the multiscale statistic is defined. The data frame constains the following two columns   
  
  \code{left} left index of an interval
  
  \code{right} right index of an interval
  
  By default, it is set to the sparse interval system proposed by Rivera and Walther (2013), see \code{\link{genIntv}} and also Li et al. (2016).
}
  \item{mode}{
  \code{"Con"} for continuous distribution functions (default)
  
  \code{"Gen"} for general (possibly with discontinuous) distribution functions
  
  See Li et al. (2016) for further details.
}
  \item{...}{
  further arguments passed to function \code{\link{quantile}}. 
}
}
\details{
  Empirically, it turns out that the quantile of the multiscale statistic converges fast to that of the limit distribution as the number of observations \code{n} increases. Thus, for the sake of computational efficiency, the quantile with \code{n} = 10,000 are used by default for that with \code{n} > 10,000, which has already been precomputed and stored. Of course, for arbitrary sample size \code{n}, one can always simulate the quantile by setting  \code{is.sim = TRUE}, and use the precomputed value by setting \code{is.sim = FALSE}. For a given sample size \code{n}, simulations are once computed, and then automatically recorded in the R memory for later usage. For memory efficiency, only the last simulation is stored.  
}
\value{
  A vector of length \code{length(alpha)} is returned, the same structure as returned by funtion \code{\link{quantile}} with option \code{\link{names} = FALSE}; The values are the (1-\code{alpha})-quantile(s) of the null distribution of the multiscale statistic via Monte Carlo simulation, corresponding to (1-\code{alpha})-confidence level(s). See Li et al. (2016) for further details.
}
\references{
Li, H., Munk, A., Sieling, H., and Walther, G. (2016). The essential histogram. arXiv:1612.07216.


Rivera, C., & Walther, G. (2013). Optimal detection of a jump in the intensity of a {P}oisson process or in a density with likelihood ratio statistics. Scand. J. Stat. 40, 752--769.
}

%\author{
%%  ~~who you are~~
%}
\note{
All the printing messages can be disabled by calling \code{\link{suppressMessages}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{checkHistogram}},
\code{\link{essHistogram}},
\code{\link{genIntv}}
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
