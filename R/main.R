# Get the quantile of the multiscale statistic via Monte Carlo simulations
msQuantile <- function(n, alpha = c(0.1), nsim = 5e3, verbose = TRUE, is.sim = (n < 1e4), ...)
{
  # ##### R implementation
  # if (exists('lastSimulation') && lastSimulation$n == n && lastSimulation$nsim == nsim) {
  #   # print("Use the last simulation!")
  #   stat = lastSimulation$stat
  # } else {
  #   # print("Start Monte-Carlo simulation ...")
  #   LR <- function(x, len, pen, n)
  #   {
  #     sqrt(2)*sqrt(pmax(len*log(len/n/x) + (n-len) * log((1-len/n)/(1-x)), 0)) - pen
  #   }
  #
  #   # compute set of intervals (Rivera & Walther '13)
  #   l   = 2:floor(log(n/log(n))/log(2))
  #   m_l = n * 2^(-l)
  #   d_l = ceiling(m_l/(6*sqrt(l)))
  #
  #   nintv = 0
  #   for (i in 1:length(d_l)) {
  #     for (len in unique(ceiling(m_l[i]/d_l[i]):floor(2*m_l[i]/d_l[i]))) {
  #       nintv = nintv + floor((n-1)/d_l[i]) - len + 1
  #     }
  #   }
  #   Left  = numeric(nintv) # vector of left bounds of intervals
  #   Right = numeric(nintv) # vector of right values of intervals
  #
  #   # print(sprintf("nintv = %d", nintv))
  #
  #   cnt = 0
  #   for (i in 1:length(d_l)) {
  #     lens = unique(ceiling(m_l[i]/d_l[i]):floor(2*m_l[i]/d_l[i]))*d_l[i]
  #     gridsize = d_l[i]
  #
  #     for (start in 0:floor((n-1)/gridsize)) {
  #       for (len in lens) {
  #         if ((1+start*gridsize + len )<= n) {
  #           cnt        = 1+cnt
  #           Left[cnt]  = 1+start*gridsize
  #           Right[cnt] = 1+start*gridsize + len
  #         }
  #       }
  #     }
  #   }
  #
  #   # print(sprintf("cnt = %d", cnt))
  #
  #   Len   = Right-Left+1
  #   Pen   = sqrt(2*(1+log(n^2/Len/(n-Len))))
  #
  #   # simulate quantile
  #   stat = numeric(nsim)
  #   for (t in 1:nsim) {
  #     U       = sort(runif(n,0,1))
  #     stat[t] = max(na.omit(LR(U[Right]-U[Left], Right-Left+1, Pen, n)))
  #   }
  #
  #   # store simulations
  #   lastSimulation <<- list(n=n, nsim=nsim, stat=stat)
  # }
  
  if (is.sim) {
    ##### fast C++ implementation
    stat = .msQuantile(n, nsim, verbose)
  } else {
    path = system.file("extdata", package = "essHist")
    stat = readRDS(file.path(path, "simn1e4.rds"))
  }
  
  quantile(na.omit(stat), 1-alpha, ...)
}


# Get the local bounds for multiscale constraint
.msBounds <- function (y, q = unname(q))
{
  y = sort(y)
  n = length(y)

  ##### compute set of intervals
  l   = 2:floor(log(n/log(n))/log(2))
  m_l = n * 2^(-l)
  d_l = ceiling(m_l/(6*sqrt(l)))

  nintv = 0
  for (i in 1:length(d_l)) {
    for (len in unique(ceiling(m_l[i]/d_l[i]):floor(2*m_l[i]/d_l[i]))) {
      nintv = nintv + floor((n-1)/d_l[i]) - len + 1
    }
  }
  Left  = numeric(nintv) # vector of left bounds of intervals
  Right = numeric(nintv) # vector of right values of intervals
  Lower = rep(NA, nintv) # lower bounds for average density on intervals
  Upper = rep(NA, nintv) # upper bounds for average density on intervals

  cnt = 0
  for (i in 1:length(d_l)){
    lens = unique(ceiling(m_l[i]/d_l[i]):floor(2*m_l[i]/d_l[i]))*d_l[i]
    gridsize = d_l[i]
    for (start in 0:floor((n-1)/gridsize)) {
      for (len in lens) {
        if ((1+start*gridsize + len )<= n) {
          cnt        = 1 + cnt
          Left[cnt]  = 1 + start*gridsize
          Right[cnt] = 1 + start*gridsize + len
        }
      }
    }
  }
  Len   = Right-Left+1

  #### Compute upper and lower bounds via quasi-newton (see package stepR for further details)
  # tol   = 1e-9
  maxIt = 10

  len = unique(Len)
  pen = sqrt(2*(1+log(n^2/len/(n-len)))) # vector of length-dependent penalty

  a     = len * log(len/n) + (n-len) * log(1 - len/n) - 1/2 * (q + pen)^2
  lower = (len/n) * 0.99
  z     = atanh( 2 * lower - 1 )
  # zprob = atanh( 2 * (len/n) - 1 )
  for(i in 1:maxIt) {
    #     print(lower)
    #     if(all(abs( len * log(lower) + (n-len) * log(1 - lower) - a ) < tol)) {
    #       i <- NA
    #       break
    #     }
    a2    = -0.5 * len * ( 1 - tanh(z)^2 )
    a1    = len * (1 - tanh(z)) - (n-len) * (1 + tanh(z))
    a0    = a - len * log(lower) - (n-len) * log(1 - lower)
    p2    = a1 / a2 / 2
    root  = pmax(p2^2 + a0 / a2, 0)
    z     = z + ifelse(root > 0 & len != 0, -p2 - sqrt(root), a0 / a1)
    lower = (tanh(z) + 1) / 2
  }
  upper = 1 - ( 1 - (len/n) ) * 0.9
  z     = atanh( 2 * upper - 1 )
  for(i in 1:maxIt) {
    # if(all(abs( len * log(upper) + (n-len) * log(1 - upper) - a ) < tol)) {
    #   i <- NA
    #   break
    # }
    a2    = -0.5 * len * ( 1 - tanh(z)^2 )
    a1    = len * (1 - tanh(z)) - (n-len) * (1 + tanh(z))
    a0    = a - len * log(upper) - (n-len) * log(1 - upper)
    p2    = a1 / a2 / 2
    root  = pmax(p2^2 + a0 / a2, 0)
    z     = z + ifelse(root > 0 & len != 0, -p2 + sqrt(root), a0 / a1)
    upper = (tanh(z) + 1) / 2
  }

  if (any(is.na(lower))) {
    # print(sprintf("Lower bounds: Intervals of length %d are ignored (n = %d)!", len[is.na(lower)], n))
    lower[is.na(lower)] = -Inf
  }
  if (any(is.na(upper))) {
    # print(sprintf("Upper bounds: Intervals of length %d are ignored (n = %d)!", len[is.na(upper)], n))
    upper[is.na(upper)] = Inf
  }

  # update Bounds
  for(i in 1:length(len)) {
    Lower[Len==len[i]] = lower[i]
    Upper[Len==len[i]] = upper[i]
  }

  ########### matrix with all informations
  eLen = y[Right] - y[Left]
  bnd  = data.frame(li = Left, ri = Right, lower = Lower/eLen, upper = Upper/eLen)
  bnd  = bnd[order(bnd$li, bnd$ri),]
  st   = cumsum(sapply(tapply(bnd$li, ordered(bnd$li, levels = 1:max(bnd$ri)), identity), length))
  st   = c(0, st[-length(st)]) # C-style
  st[is.na(tapply(bnd$li, ordered(bnd$li, levels = 1:max(bnd$ri)), length))] = NA
  si   = c(st[!is.na(st)], nrow(bnd)) + 1
  feas = sapply(1:(length(si)-1), function(i) with(bnd[si[i]:(si[i+1]-1),], {wi <- ri == li[1]; if(any(wi)) max(lower[wi]) <= min(upper[wi]) else TRUE}))
  list(bounds = bnd, start = as.integer(st),  feasible = all(feas))
}


# Compute the essential histogram
# essHistogram <- function(y, alpha = 0.1, q = NA, confband = FALSE)
essHistogram <- function(x, alpha = 0.5, q = NA, plot = TRUE, verbose =TRUE, xname = deparse(substitute(x)), ...)
{
  y = sort(x)
  n = length(y)
  # Tackle duplicated data samples: slightly shift the samples by 'eps'
  dupId = duplicated(y)
  if (any(dupId)) {
    nDup = diff(c(which(!dupId), n+1))
    nDup = nDup[nDup > 1] - 1
    eps  = max(1e-12, 0.1*min(diff(y[!dupId]))/max(nDup))
    y[dupId] = y[dupId] + eps*sequence(nDup)
  }

  if (is.na(q)) q = msQuantile(n = n, alpha = alpha, verbose = verbose)
  bnd = .msBounds(y, q = q)
  ##### fast C++ implementation
  if (verbose) cat("Dynamic programming ...")
  eh = .boundedHistogram(cumsum(y), bnd$start, bnd$bounds$ri-1, bnd$bounds$lower, bnd$bounds$upper)
  if (verbose) cat(" ... end!\n")
  
  nSeg       = length(eh$value)
  breakPoint = c(3*y[1]-y[2], y[eh$rightIndex[-nSeg]]+y[eh$rightIndex[-nSeg]+1], 3*y[n]-y[n-1])/2
  ret = list("breaks"   = breakPoint,
             "counts"   = round(eh$value*diff(breakPoint)*n), 
             "density"  = eh$value,
             "mids"     = (breakPoint[1:nSeg]+breakPoint[2:(nSeg+1)])/2,
             "xname"    = xname,
             "equidist" = FALSE) # (var(diff(breakPoint)) <= 1e-16)
  class(ret) = "histogram"
  if (plot) plot(ret, ...)
  ret
}

.preSimulation <- function() {
  scan()
}
