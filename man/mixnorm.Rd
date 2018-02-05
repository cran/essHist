\name{Mixed normals}
\alias{rmixnorm}
\alias{pmixnorm}
\alias{dmixnorm}
\alias{paramExample}
\title{
The mixture of normal distributions
}
\description{
Density, distribution function and random generation for the mixture of normals with each component specified by \code{mean} and \code{sd}, and mixture weights by \code{prob}. \code{paramExample} gives detailed parameters for some examples specified by \code{type}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dmixnorm(x, mean, sd, prob = rep(1/length(mean),length(mean)), type = NULL, ...)
pmixnorm(x, mean, sd, prob = rep(1/length(mean),length(mean)), type = NULL, ...)
rmixnorm(n, mean, sd, prob = rep(1/length(mean),length(mean)), type = NULL)
paramExample(type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  vector of locations.
  } 
  \item{n}{
  integer; number of observations.
  } 
  \item{mean}{
  vector of means for each mixture component. 
  } 
  \item{sd}{
  vector of standard deviations for each mixture component. 
  }
  \item{prob}{
  vector of prior probability for each mixture component (i.e. mixture weights). 
  }
  \item{type}{
  a (case insensitive) character string of example name; It includes examples from Marron & Wand (1992): "MW1", ..., "MW15", or equivalently "guass", "skewed_unimodal", "strong_skewed", "kurtotic_unimodal", "outlier", "bimodal", "separated_bimodal", "skewed_bimodal", "trimodal", "claw", "double_claw", "asymmetric_claw", "asymmetric_double_claw", "smooth_comb", "discrete_comb"; It also includes "harp" example from Li et al. (2016).
  }
  \item{...}{
  further arguments passed to \code{\link{dnorm}} and \code{\link{pnorm}}.
  }
}
\details{
Users have to either provide \code{mean}, \code{sd} and optionally \code{prob}; or \code{type}. In case of providing \code{type}, the values of \code{mean}, \code{sd} and \code{prob} are ignored.  

If \code{prob} is not specified it assumes the default value of equal weights. Each component is computed via \code{\link{dnorm}}, \code{\link{pnorm}} and \code{\link{rnorm}}.
}
\value{
\code{dmixnorm} gives the density, \code{pmixnorm} gives the distribution function, and \code{rmixnorm} generates random deviates.

The length of the result is determined by \code{n} for \code{rmixnorm}, and is the length of \code{x} for \code{dmixnorm} and \code{pmixnorm}.

\code{paramExample} gives a data frame with components \code{mean}, \code{sd} and \code{prob}. 
}
\references{
%% ~put references to the literature/web site here ~
Li, H., Munk, A., Sieling, H., and Walther, G. (2016).The essential histogram. arXiv:1612.07216.

Marron, J. S., & Wand, M. P. (1992). Exact mean integrated squred error. The Annals of Statistics, 20(2), 712--736.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\link{Normal} for standard normal distributions; \link{Distributions} for other standard distributions.  
}
\examples{
## Example claw 
type = "claw" # or equivalently "MW10"
#     generate random numbers
n = 500
Y = rmixnorm(n, type = type)
#     compute the density
x = seq(min(Y), max(Y), length.out = n)
f = dmixnorm(x, type = type)
#     compute the distribution
F = pmixnorm(x, type = type)
#     plots
op = par(mfrow = c(1,2))
plot(x, f, type = "l", main = "Claw Density")
points(Y, rep(0,n))
plot(x, F, type = "l", main = "Claw Distribution")
points(Y, rep(0,n))
par(op)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}% use one of  RShowDoc("KEYWORDS")
\keyword{distribution}% __ONLY ONE__ keyword per line
