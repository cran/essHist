checkHistogram <- function(h, x, alpha = 0.1, q = NULL, intv = NULL, 
                           mode = ifelse(anyDuplicated(x),"Gen","Con"),
                           plot = TRUE, xlim = NULL, ylim = NULL, 
                           xlab = "", ylab = "", yaxt = "n", ...) { 
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  sy = sort(x[is.finite(x)])
  n  = as.integer(length(sy))
  if (is.na(n) || n < 2) 
    stop("length(x) must be > 1!") 
  if ("histogram" %in% class(h)) 
    h = .hist2vec(h, sy)#$h
  if (!is.numeric(h))
    stop("'h' constains non-numeric values!")
  h = h[is.finite(h)]
  if (n != length(h)) 
    stop("length of 'x' does not match that of 'h'!")
  # prepare data
  ris = c(which(diff(h) != 0), n)
  lis = c(1, which(diff(h) != 0)) # lis = c(1, which(diff(h) != 0) + 1)
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
  else if (!missing(alpha))
    warning("Use input 'q' and ignore 'alpha'!")
  if (length(q) > 1) {
    q = q[1]
    warning("Length of 'q' or 'alpha' > 1: use only the first value, and ignore the others!") 
  }
  if (!is.finite(q)) 
    stop("'q' must be a finite value!")
  minQ = .minThreshold(n,intv)
  if (q < minQ)
    stop(sprintf("Empty constraint set: threshold 'q' should be >= %g!", minQ))
  # parse discrete data
  intv = .validInterval(intv, sy)
  bnds = .intv2Bounds(intv, sy, c(0,2:n), q)$bounds

  # check feasibility for each local constraint
  message("Find violated intervals ...", appendLF=FALSE)
  nintv     = nrow(bnds)
  indexIntv = rep(FALSE, nintv)
  for (i in 1:length(lis)) 
    indexIntv = indexIntv | ((bnds$li >= lis[i]) & (bnds$ri <= ris[i]))
  indexIntv[indexIntv] = (h[bnds$ri[indexIntv]] > bnds$upper[indexIntv]) | (h[bnds$ri[indexIntv]] < bnds$lower[indexIntv])
  if (any(indexIntv)) 
    vioIntvs = data.frame(leftIndex  = bnds$li[indexIntv],
                          rightIndex = bnds$ri[indexIntv],
                          leftEnd    = sy[bnds$li[indexIntv]],
                          rightEnd   = sy[bnds$ri[indexIntv]])
  else
    vioIntvs = data.frame(NULL)
  message(" ... end!")
  
  # check whether jumps are removable
  message("Find removable jumps ...", appendLF=FALSE)
  jmpflag = rep(0, length(lis)-1)
  if (length(lis) > 1) {
    for (i in 1:(length(lis)-1)) {
      for (j in (i+1):length(ris)){
        indexIntv = (bnds$li >= lis[i]) & (bnds$ri <= ris[j])
        height = (ris[j]-ris[i]+(ris[i]==1))/n/(sy[ris[j]] - sy[lis[i]])
        if (any(indexIntv) && (min(bnds$upper[indexIntv]) < height || height < max(bnds$lower[indexIntv]))) 
          break
        else 
          jmpflag[i:(j-1)] = jmpflag[i:(j-1)] + 1
      }
    }
  }
  message(" ... end!")
    
  # visualization
  if (plot == TRUE) {
    #     graphical parameter
    good_ratio = 0.618
    if (is.null(xlim))
      xlim = range(sy)
    if (is.null(ylim)) {
      ymax = 1.8*max(h)
      ymin = -good_ratio * ymax * ((nrow(vioIntvs) > 0) * 0.9 + 0.1)
      ylim = c(ymin,ymax)
    } else {
      ymin = ylim[1]
      ymax = ylim[2]
    }
    #     draw estimator
    plot(NULL, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, yaxt = yaxt, ...)
    yat = par()$yaxp
    yat = seq(yat[1], yat[2], length.out = yat[3]+1)
    yat = unique(c(0, yat[yat >= 0]))
    axis(2, at = yat)
    lines(sy, h, type = "s")
    #     draw removable jumps
    if (any(jmpflag != 0)) {
      jmpflag = 1 - jmpflag/max(jmpflag)
      for (i in 1:length(jmpflag)) {
        if (jmpflag[i] != 1)
          segments(sy[lis[i+1]], -0.1 * ymax * good_ratio, sy[lis[i+1]], -0.05 * ymax * good_ratio)
        #, col = rgb(jmpflag[i], jmpflag[i], jmpflag[i])
      }
      rmBrk = sy[lis[(jmpflag != 1)+1]]
    } else
      rmBrk = numeric(0)
    abline(h = -0.1 * ymax * good_ratio, col = "lightblue")
    if (nrow(vioIntvs) > 0) {
      #     draw a strap
      strap = rep(0, n)
      for (i in 1:nrow(vioIntvs)) 
        strap[vioIntvs$leftIndex[i]:vioIntvs$rightIndex[i]] = strap[vioIntvs$leftIndex[i]:vioIntvs$rightIndex[i]] + 1
      ymid      = c((sy[1:(n-1)]+sy[2:n])/2, sy[n])
      strap = 1 - strap / max(strap)
      for (i in 1:n) 
        rect(ymid[i], 0.23*ymin, ymid[i+1], 0.18*ymin, col = rgb(strap[i], strap[i], strap[i]), border = NA)
      #     draw intervals
      for (i in 1:nrow(vioIntvs)) {
        z = runif(1) # z = (i %% nlev) / nlev
        lines(sy[vioIntvs$leftIndex[i]:vioIntvs$rightIndex[i]],
              rep(z*ymin*0.69+ymin*0.31, vioIntvs$rightIndex[i]-vioIntvs$leftIndex[i]+1),
              col = rgb(0.9, 0.9, 0.9))
      }
      mtext("Intervals of violation", side = 2, line = 1, at = 0.7*ymin)
    }
  }

  # return value
  list("violatedIntervals" = vioIntvs, "removableBreakpoints" = rmBrk)
}

.hist2vec <- function(h, x) {
  if (length(h$mids) != length(h$breaks)-1 || length(h$mids) != length(h$density)
      || length(h$mids) != length(h$counts)) {
    stop("'h' is invalid as a histogram object!")
  }
  n = sum(h$counts)
  # if (missing(x)) 
  #   x = seq(min(h$breaks), max(h$breaks), length.out = n)
  # else 
  #   x = sort(x)
  if (length(x) != n) 
    stop("length of 'x' does not match that of 'h'!")
  dval = numeric(n)
  dval[x >= h$breaks[1] & x <= h$breaks[2]] = h$density[1]
  if (length(h$breaks) > 2) {
    for (i in 2:(length(h$breaks)-1)) 
      dval[x > h$breaks[i] & x <= h$breaks[i+1]] = h$density[i]
  }
  # data.frame("x" = x, "h" = dval)
  dval
}
