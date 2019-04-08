###########################################################
##########  Mixture of Gaussians  #########################
###########################################################
# random number generator
rmixnorm <- function(n, mean = c(0), sd = rep(1,length(mean)), 
                     prob = rep(1,length(mean)), type = NULL) {
  if (!(length(n) == 1 && is.finite(n) && n == round(n) && n > 0))
    stop("Invalid input 'n'!")
  if (is.null(type)) {
    if (!is.numeric(mean))
      stop("'mean' must be numeric!")
    mean = mean[is.finite(mean)]
    nc   = length(mean)
    if (nc == 0)
      stop("Input 'mean' must be finite!")
    if (!is.numeric(sd))
      stop("'sd' must be numeric!")
    sd = sd[is.finite(sd)]
    if (any(sd < 0) || length(sd) == 0)
      stop("'sd' must be positive and finite!")
    if (!is.numeric(prob))
      stop("'prob' must be numeric!")
    prob = prob[is.finite(prob) & prob >= 0]
    if (length(prob) == 0 || sum(prob) == 0)
      stop("'prob' must be positive!")
    prob = prob/sum(prob)
    if (nc != length(sd) || length(prob) != nc) 
      stop("Input 'mean', 'sd', and 'prob' must have the same length!")
  } else {
    param = paramExample(type)
    mean  = param$mean
    sd    = param$sd
    prob  = param$prob
    nc    = length(mean)
    # ret   = rmixnorm(n, mean = mean, sd = sd, prob = prob)
  }
  beta = sample(1:nc, size = n, replace = TRUE, prob = prob)
  rnorm(n, mean = mean[beta], sd = sd[beta])
}
# density function
dmixnorm <- function(x, mean = c(0), sd = rep(1,length(mean)), 
                     prob = rep(1,length(mean)), type = NULL, ...) {
  if (!is.numeric(x))
    stop("'x' must be numeric!")
  x = x[is.finite(x)]
  if (length(x) == 0)
    stop("Invalid length(x)!")
  if (is.null(type)) {
    if (!is.numeric(mean))
      stop("'mean' must be numeric!")
    mean = mean[is.finite(mean)]
    nc   = length(mean)
    if (nc == 0)
      stop("Input 'mean' must be finite!")
    if (!is.numeric(sd))
      stop("'sd' must be numeric!")
    sd = sd[is.finite(sd)]
    if (any(sd < 0) || length(sd) == 0)
      stop("'sd' must be positive and finite!")
    if (!is.numeric(prob))
      stop("'prob' must be numeric!")
    prob = prob[is.finite(prob) & prob >= 0]
    if (length(prob) == 0 || sum(prob) == 0)
      stop("'prob' must be positive!")
    prob = prob/sum(prob)
    if (nc != length(sd) || length(prob) != nc) 
      stop("Input 'mean', 'sd', and 'prob' must have the same length!")
  } else {
    param = paramExample(type)
    mean  = param$mean
    sd    = param$sd
    prob  = param$prob
    nc    = length(mean)
    # ret   = dmixnorm(x, mean = mean, sd = sd, prob = prob, ...)
  }
  n   = length(x)
  ret = rep(0, n)
  for (i in 1:nc) 
    ret = ret + prob[i]*dnorm(x, mean = mean[i], sd = sd[i], ...)
  ret
}
# distribution function
pmixnorm <- function(x, mean = c(0), sd = rep(1,length(mean)), 
                     prob = rep(1,length(mean)), type = NULL, ...) {
  if (!is.numeric(x))
    stop("'x' must be numeric!")
  x = x[is.finite(x)]
  if (length(x) == 0)
    stop("Invalid length(x)!")
  if (is.null(type)) {
    if (!is.numeric(mean))
      stop("'mean' must be numeric!")
    mean = mean[is.finite(mean)]
    nc   = length(mean)
    if (nc == 0)
      stop("Input 'mean' must be finite!")
    if (!is.numeric(sd))
      stop("'sd' must be numeric!")
    sd = sd[is.finite(sd)]
    if (any(sd < 0) || length(sd) == 0)
      stop("'sd' must be positive and finite!")
    if (!is.numeric(prob))
      stop("'prob' must be numeric!")
    prob = prob[is.finite(prob) & prob >= 0]
    if (length(prob) == 0 || sum(prob) == 0)
      stop("'prob' must be positive!")
    prob = prob/sum(prob)
    if (nc != length(sd) || length(prob) != nc) 
      stop("Input 'mean', 'sd', and 'prob' must have the same length!")
  } else {
    param = paramExample(type)
    mean  = param$mean
    sd    = param$sd
    prob  = param$prob
    nc    = length(mean)
    # ret   = pmixnorm(x, mean = mean, sd = sd, prob = prob, ...)
  }
  ret = rep(0, length(x))
  for (i in 1:nc) 
    ret = ret + prob[i]*pnorm(x, mean = mean[i], sd = sd[i], ...)
  ret
}
# parameters for some examples: (Marron & Wand '92), harp
paramExample <- function(type) {
  supportTypeS = c("mw1", "mw2", "mw3", "mw4", "mw5", "mw6", "mw7", "mw8", 
                   "mw9", "mw10", "mw11", "mw12", "mw13", "mw14", "mw15",
                   "gauss", "skewed_unimodal", "strong_skewed",  "claw",
                   "kurtotic_unimodal", "outlier", "bimodal", "trimodal",
                   "separated_bimodal", "skewed_bimodal", "double_claw",
                   "asymmetric_claw", "asymmetric_double_claw", "harp", 
                   "smooth_comb", "discrete_comb")
  type = match.arg(tolower(type), supportTypeS)
  if (type %in% c("mw1", "guass")) { 
    mean = 0
    sd   = 1
    prob = 1
  } else if (type %in% c("mw2","skewed_unimodal")) {
    mean = c(0, 1/2, 13/12)
    sd   = c(1, 2/3, 5/9)
    prob = c(1/5, 1/5, 3/5)
  } else if (type %in% c("mw3","strong_skewed")) { 
    mean = 3*((2/3)^(0:7) - 1)
    sd   = (2/3)^(0:7)
    prob = rep(1/8, 8)
  } else if (type %in% c("mw4","kurtotic_unimodal")) { 
    mean = c(0, 0)
    sd   = c(1, 1/10)
    prob = c(2/3, 1/3)
  } else if (type %in% c("mw5","outlier")) { 
    mean = c(0, 0)
    sd   = c(1, 1/10)
    prob = c(1/10, 9/10)
  } else if (type %in% c("mw6","bimodal")) {
    mean = c(-1, 1)
    sd   = c(2/3, 2/3)
    prob = c(1/2, 1/2)
  } else if (type %in% c("mw7","separated_bimodal")) { 
    mean = c(-3/2, 3/2)
    sd   = c(1/2, 1/2)
    prob = c(1/2, 1/2)
  } else if (type %in% c("mw8","skewed_bimodal")) { 
    mean = c(0, 3/2)
    sd   = c(1, 1/3)
    prob = c(3/4, 1/4)
  } else if (type %in% c("mw9","trimodal")) {  
    mean = c(-6/5, 6/5, 0)
    sd   = c(3/5, 3/5, 1/4)
    prob = c(9/20, 9/20, 1/10)
  } else if (type %in% c("mw10","claw")) { 
    mean = c(0, (0:4)/2-1)
    sd   = c(1, rep(1/10, 5))
    prob = c(1/2, rep(1/10, 5))
  } else if (type %in% c("mw11","double_claw")) {
    mean = c(-1, 1, (-3:3)/2)
    sd   = c(2/3, 2/3, rep(1/100, 7))
    prob = rep(c(49/100, 1/350), c(2, 7))
  } else if (type %in% c("mw12","asymmetric_claw")) { 
    mean = c(0, (-2:2)+1/2)
    sd   = c(1, 2^(2:-2)/10)
    prob = c(1/2, 2^(3:-1)/31)
  } else if (type %in% c("mw13","asymmetric_double_claw")) {  
    mean = c(2*(0:1)-1, -(1:3)/2, (1:3)/2)
    sd   = rep(c(2/3, 1/100, 7/100), c(2, 3, 3))
    prob = rep(c(46/100, 1/300, 7/300), c(2, 3, 3))
  } else if (tolower(type) %in% c("mw14","smooth_comb")) {
    mean = (65-96*(1/2)^(0:5))/21
    sd   = 32/63/2^(0:5)
    prob = 2^(5:0)/63
  } else if (tolower(type) %in% c("mw15","discrete_comb")) {  
    mean = c((12*(0:2)-15), 2*(8:10))/7
    sd   = rep(c(2/7,1/21), c(3,3))
    prob = rep(c(2/7,1/21), c(3,3))
  } else if (tolower(type) == 'harp') {
    mean = c(0,   5, 15, 30, 60)
    sd   = c(0.5, 1,  2,  4,  8)
    prob = rep(0.2, 5)
  } else {
    stop("Unknown input 'type'!") 
  }
  data.frame("mean" = mean, "sd" = sd, "prob" = prob) 
}
