# Get the quantile of the multiscale statistic via Monte Carlo simulations
msQuantile <- function(n, alpha = c(0.1), nsim = 5e3, is.sim = (n < 1e4), intv = genIntv(n),
                       verbose = TRUE, ...) {
  if (is.sim) {
    ##### fast C++ implementation
    stat = .msQuantile(as.integer(intv$left-1), as.integer(intv$right-1), n, nsim, verbose)
  } else {
    if (verbose) { cat('Load precomputed quantiles ... ') }
    path = system.file("extdata", package = "essHist")
    stat = readRDS(file.path(path, "simn1e4.rds"))
    if (verbose) { cat('... end!\n') }
  }
  quantile(na.omit(stat), 1-alpha, ...)
}

# # Debug: for simulation of Tn
# dms_stat <- function(n, nsim) {
#   intv = genIntv(n)
#   .msQuantile(as.integer(intv$left-1), as.integer(intv$right-1), n, nsim, TRUE)
# }

# Generate system of intervals of form (a, b] with a < b
genIntv <- function (n, type = c('Sparse', 'Full')) {
  type = match.arg(type)
  if (type == 'Sparse') {
    l     = 2:floor(log(n/log(n))/log(2))
    m_l   = n * 2^(-l)
    d_l   = ceiling(m_l/(6*sqrt(l)))
    nintv = 0
    for (i in 1:length(d_l)) {
      for (len in unique(ceiling(m_l[i]/d_l[i]):floor(2*m_l[i]/d_l[i]))) {
        nintv = nintv + floor((n-1)/d_l[i]) - len + 1
      }
    }
    left  = numeric(nintv) # vector of left bounds of intervals
    right = numeric(nintv) # vector of right values of intervals
    cnt   = 0
    for (i in 1:length(d_l)){
      lens = unique(ceiling(m_l[i]/d_l[i]):floor(2*m_l[i]/d_l[i]))*d_l[i]
      gridsize = d_l[i]
      for (start in 0:floor((n-1)/gridsize)) {
        for (len in lens) {
          if ((1+start*gridsize + len )<= n) {
            cnt        = 1 + cnt
            left[cnt]  = 1 + start*gridsize
            right[cnt] = 1 + start*gridsize + len
          }
        }
      }
    }
  } else if (type == 'Full') {
    left  = rep(1:(n-1),(n-1):1)
    right = unlist(lapply(1:(n-1), function(x) (x+1):n))
  } 
  data.frame('left'=left,'right'=right)
}

# ensure 'intv' are valid
.validInterval <- function (intv, y) {
  y = sort(y)
  n = length(y)
  ind.valid = (intv$left < intv$right) & (intv$left >= 1) & (intv$right <= n)
  left  = intv$left[ind.valid]
  right = intv$right[ind.valid]
  # adjust for discrete data 
  dupId = duplicated(rev(y))
  if (any(dupId)) {
    # # method A: remove all improper intervals
    # ind.valid  = (left %in% left.valid) & (right %in% right.valid)
    # left  = left[ind.valid]
    # right = right[ind.valid]
    
    # # method B: modify improper intervals
    # right.valid = c((n:1)[!dupId],1)    # decreasing
    # nintv = length(left)
    # #     move right end towards right
    # right.r = sapply(1:nintv, function (x) tail(right.valid[right.valid>=right[x]],n=1))
    # #     move left end towards right
    # left.r  = sapply(1:nintv, function (x) tail(right.valid[right.valid>=left[x]],n=1))
    # #     move left end towards left
    # right.valid = rev(right.valid)
    # left.l  = sapply(1:nintv, function (x) tail(right.valid[right.valid<=left[x]],n=1))
    # #     move right end towards left
    # right.l = sapply(1:nintv, function (x) tail(right.valid[right.valid<=right[x]],n=1))
    # #     combine to make intervals
    # right = rep(c(right.l, right.r),2)
    # left  = c(rep(left.l,2),rep(left.r,2))
    # ind.valid = left < right
    # left  = left[ind.valid]
    # right = right[ind.valid]
    
    # method C: a faster way of modifying improper intervals
    right.valid = c((n:1)[!dupId])
    left.valid  = c(right.valid[-1],1)
    ind.valid   = (left %in% left.valid) & (right %in% right.valid)
    if (!all(ind.valid)) {
      itEnd  = unique(data.frame('leE'=y[intv$left[!ind.valid]],'riE'=y[intv$right[!ind.valid]]))
      nitE   = length(itEnd$leE)
      # move left endpoint leftward
      leMd.l = sapply(1:nitE,function (x) max(which.max(y >= itEnd$leE[x])-1,1))
      # move left endpoint rightward
      leMd.r = sapply(1:nitE, function (x) tail((1:n)[y <= itEnd$leE[x]],n=1))
      # move right endpoint leftward
      riMd.l = sapply(1:nitE,function (x) max(which.max(y >= itEnd$riE[x])-1,1))
      # move right endpoint rightward
      riMd.r = sapply(1:nitE, function (x) tail((1:n)[y <= itEnd$riE[x]],n=1))
      riMd = rep(c(riMd.l, riMd.r),2)
      leMd = c(rep(leMd.l,2),rep(leMd.r,2))
      idMd = leMd < riMd
      leMd = leMd[idMd]
      riMd = riMd[idMd]
    }
    left  = c(left[ind.valid],leMd)
    right = c(right[ind.valid],riMd)
  }
  # remove repetitive intervals
  unique(data.frame('left'=left, 'right'=right))
}

# Precompute the bounds on every interval
.intv2Bounds <- function (intv, y, q) {
  y = sort(y)
  n = length(y)
  q = unname(q)
  
  nintv = nrow(intv)
  Left  = intv$left
  Right = intv$right
  
  od = order(Left, Right)
  Left  = Left[od]
  Right = Right[od]
  Len   = Right-Left+(Left == 1) # Right-Left+1
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
  for (i in 1:maxIt) {
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
  
  nintv = length(Len)
  Lower = rep(NA, nintv) # lower bounds for average density on intervals
  Upper = rep(NA, nintv) # upper bounds for average density on intervals
  # update Bounds
  for (i in 1:length(len)) {
    Lower[Len==len[i]] = lower[i]
    Upper[Len==len[i]] = upper[i]
  }
  
  # remove redudant constraints
  eLen = y[Right] - y[Left]
  vInd  = is.finite(Lower) & is.finite(Upper) & (eLen != 0)
  Lower = Lower[vInd]
  Upper = Upper[vInd]
  Left  = Left[vInd]
  Right = Right[vInd]
  eLen  = eLen[vInd]
  
  ########### matrix with all informations
  bnd  = data.frame(li = Left, ri = Right, lower = Lower/eLen, upper = Upper/eLen)
  remove('Lower', 'Upper', 'Left', 'Right', 'eLen', 'Len', 'lower', 'upper', 'len', 'vInd')
  bnd  = bnd[order(bnd$li, bnd$ri),]
  # st   = cumsum(sapply(tapply(bnd$li, ordered(bnd$li, levels = 1:max(bnd$ri)), identity), length))
  # st   = c(0, st[-length(st)]) # C-style
  # st[is.na(tapply(bnd$li, ordered(bnd$li, levels = 1:max(bnd$ri)), length))] = NA
  st = rep(NA,n)
  st[bnd$li[1]] = 0
  aux = which(diff(bnd$li) != 0)
  st[bnd$li[aux+1]] = aux
  si   = c(st[!is.na(st)], nrow(bnd)) + 1
  feas = sapply(1:(length(si)-1), function(i) with(bnd[si[i]:(si[i+1]-1),], {wi <- ri == li[1]; if(any(wi)) max(lower[wi]) <= min(upper[wi]) else TRUE}))
  list(bounds = bnd, start = as.integer(st),  feasible = all(feas))
}

# Compute the essential histogram
essHistogram <- function(x, alpha = 0.5, q = NA, intv = genIntv(length(x)), plot = TRUE, 
                         verbose = TRUE, xname = deparse(substitute(x)), ...) {
  y = sort(x)
  n = length(y)
  intv = .validInterval(intv, y)
  # print(sprintf('Number of intervals: %d', length(intv$left)))
  if (is.na(q)) q = msQuantile(n, alpha, intv = intv, verbose = verbose)
  if (verbose) cat("Dynamic programming ...") 
  bnd = .intv2Bounds(intv, y, q)
  ##### fast C++ implementation
  eh = .boundedHistogram(y, as.integer(bnd$start), as.integer(bnd$bounds$ri-1), 
                         bnd$bounds$lower, bnd$bounds$upper)
  if (verbose) cat(" ... end!\n")
  
  nSeg       = length(eh$value)
  # breakPoint = c(3*y[1]-y[2], y[eh$rightIndex[-nSeg]]+y[eh$rightIndex[-nSeg]+1], 3*y[n]-y[n-1])/2
  breakPoint = c(y[1],y[eh$rightIndex])
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

# .preSimulation <- function() {
#   scan()
# }
