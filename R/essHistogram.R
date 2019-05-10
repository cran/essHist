# Compute the essential histogram
essHistogram <- function(x, alpha = 0.5, q = NULL, intv = NULL, 
                         plot = TRUE, mode = ifelse(anyDuplicated(x),"Gen","Con"), 
                         xname = deparse(substitute(x)), ...) { 
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  y = sort(x[is.finite(x)])
  n = as.integer(length(y))
  if (is.na(n)) 
    stop("invalid length(x)") 
  if (length(unique(y)) < 2) {# call R-built hist for 1 data point
    message(sprintf('Only %d distinct finite observations, call "graphics::hist()" instead!',length(unique(y))))
    return (hist(x,freq=FALSE,plot=plot,...))
  }
  if (is.null(intv))
    intv = genIntv(n)
  if (!('left' %in% names(intv) && 'right' %in% names(intv)))
    stop("'intv' must have two fields 'left' and 'right'!")
  intv = intv[(intv$left<intv$right)&(intv$left>=1)&(intv$right<=n),]
  if (nrow(intv) < 1) 
    stop("No valid intervals in 'intv'!")
  if (is.null(q)) 
    q = msQuantile(n, alpha, intv = intv, mode = mode)
  else if (!missing(alpha) || !missing(mode) || !missing(intv))
    warning("Use input 'q' and ignore 'alpha', 'intv' or 'mode'!")
  if (length(q) > 1) {
    q = q[1]
    warning("Length of 'q' or 'alpha' > 1: use only the first value, and ignore the others!") 
  }
  if (!is.finite(q)) 
    stop("'q' must be a finite value!")
  minQ = .minThreshold(n,intv)
  if (q < minQ)
    stop(sprintf("Empty constraint set: threshold 'q' should be >= %g!", minQ))
  message(sprintf('The threshold is %g', q))
  message("Dynamic programming ...", appendLF=FALSE) 
  # parse discrete data
  intv   = .validInterval(intv, y)
  cumcnt = c(0,which(diff(y)!=0),n) # cumulative counts
  y      = unique(y)                # unique data values
  # y      = c(2*y[1]-y[2], y)
  y = c(y[1],(y[1]+y[2])/2,y[-1]) # to include the first observation 
                                  # by allowing break at the middle point 
                                  # of the first two distinct observations
  nintv = nrow(intv)
  intv$left  = sapply(1:nintv, function (x) { if (intv$left[x] == 1) 1 else which(intv$left[x] == cumcnt)})
  intv$right = sapply(1:nintv, function (x) which(intv$right[x] == cumcnt))
  if (length(intv$left)!=nintv || length(intv$right)!=nintv)
    stop('There are some invalid intervals')
  bnd = .intv2Bounds(intv, y, cumcnt, q)
  ##### fast C++ implementation
  eh = .boundedHistogram(y, cumcnt, as.integer(bnd$start), as.integer(bnd$bounds$ri-1), 
                         bnd$bounds$lower, bnd$bounds$upper)
  message(" ... end!")
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
  if (plot) 
    plot(ret, ...)
  ret
}

# Precompute the bounds on every interval
.intv2Bounds <- function (intv, y, cumcnt, q) {
  if (is.unsorted(y))
    stop("'y' must be sorted!")
  n = cumcnt[length(cumcnt)]
  # q = unname(q)
  
  nintv = nrow(intv)
  Left  = intv$left
  Right = intv$right
  
  od    = order(Left, Right)
  Left  = Left[od]
  Right = Right[od]
  Len   = cumcnt[Right]-cumcnt[Left] # Right-Left+1
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
  
  if (any(is.na(lower)))
    lower[is.na(lower)] = -Inf
  if (any(is.na(upper))) 
    upper[is.na(upper)] = Inf
  
  # nintv = length(Len)
  Lower = rep(NA, nintv) # lower bounds for average density on intervals
  Upper = rep(NA, nintv) # upper bounds for average density on intervals
  # update Bounds
  for (i in 1:length(len)) {
    Lower[Len==len[i]] = lower[i]
    Upper[Len==len[i]] = upper[i]
  }
  
  # remove redudant constraints
  eLen = y[Right] - y[Left]
  vInd  = (eLen != 0) # & is.finite(Lower) & is.finite(Upper) 
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
  st = rep(NA,length(y))
  st[bnd$li[1]] = 0
  aux = which(diff(bnd$li) != 0)
  st[bnd$li[aux+1]] = aux
  si   = c(st[!is.na(st)], nrow(bnd)) + 1
  feas = sapply(1:(length(si)-1), 
                function(i) with(bnd[si[i]:(si[i+1]-1),], 
                                 {wi <- ri == li[1]; if(any(wi)) max(lower[wi]) <= min(upper[wi]) else TRUE}))
  list(bounds = bnd, start = as.integer(st),  feasible = all(feas))
}