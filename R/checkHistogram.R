checkHistogram <- function(h, y, alpha = 0.1, q = NA, plot = TRUE)
{
  if ("histogram" %in% class(h)) h = .hist2vec(h, y)$h
  # prepare data
  if (length(h) != length(y)) 
    stop("length of data (", length(y), ") and length of estimator (", length(h),") does not match!")
  n   = length(y)
  y   = sort(y)
  ris = c(which(diff(h) != 0), n)
  lis = c(1, which(diff(h) != 0) + 1)
  if (is.na(q)) q = msQuantile(n, alpha)
  bnds = .msBounds(y, q)$bounds

  vioIntvs = data.frame() # violated intervals

  # check feasibility for each local constraint
  # print("Find violated intervals")
  nintv     = nrow(bnds)
  indexIntv = rep(FALSE, nintv)
  for (i in 1:length(lis)) 
    indexIntv = indexIntv | ((bnds$li >= lis[i]) & (bnds$ri <= ris[i]))
  indexIntv[indexIntv] = (h[bnds$li[indexIntv]] > bnds$upper[indexIntv]) | (h[bnds$li[indexIntv]] < bnds$lower[indexIntv])
  if (any(indexIntv)) 
    vioIntvs = rbind(vioIntvs, data.frame(leftIndex  = bnds$li[indexIntv],
                                          rightIndex = bnds$ri[indexIntv],
                                          leftEnd    = y[bnds$li[indexIntv]],
                                          rightEnd   = y[bnds$ri[indexIntv]]))

  # check whether jumps are removable
  # print("Find removable jumps")
  jmpflag = rep(0, length(lis)-1)
  for (i in 1:(length(lis)-1)) {
    for (j in (i+1):length(ris)){
      indexIntv = (bnds$li >= lis[i]) & (bnds$ri <= ris[j])
      if (any(indexIntv) && (min(bnds$upper[indexIntv]) < max(bnds$lower[indexIntv]))) {
        break
      } else {
        jmpflag[i:(j-1)] = jmpflag[i:(j-1)] + 1
      }
    }
  }

  # visualization
  if (plot == TRUE) {
    #     graphical parameter
    ymax = 1.8*max(h)
    good_ratio = 0.618
    ymin = -good_ratio * ymax * ((nrow(vioIntvs) > 0) * 0.9 + 0.1)
    #     draw estimator
    plot(NULL, xlim = range(y), ylim = c(ymin,ymax), xlab = "", ylab = "", yaxt = "n")
    yat = par()$yaxp
    yat = seq(yat[1], yat[2], length.out = yat[3]+1)
    yat = unique(c(0, yat[yat >= 0]))
    axis(2, at = yat)
    lines(y, h, type = "s")
    #     draw removable jumps
    if (any(jmpflag != 0)) {
      jmpflag = 1 - jmpflag/max(jmpflag)
      for (i in 1:length(jmpflag)) {
        if (jmpflag[i] != 1)
        segments(y[lis[i+1]], -0.1 * ymax * good_ratio,
                 y[lis[i+1]], -0.05 * ymax * good_ratio,
                 col = rgb(jmpflag[i], jmpflag[i], jmpflag[i]))
      }
      abline(h = -0.1 * ymax * good_ratio, col = "lightblue")
    }
    if (nrow(vioIntvs) > 0) {
      #     draw a strap
      strap = rep(0, n)
      for (i in 1:nrow(vioIntvs)) {
        strap[vioIntvs$leftIndex[i]:vioIntvs$rightIndex[i]] =
          strap[vioIntvs$leftIndex[i]:vioIntvs$rightIndex[i]] + 1
      }
      ymid      = c((y[1:(n-1)]+y[2:n])/2, y[n])
      strap = 1 - strap / max(strap)
      for (i in 1:n) {
        rect(ymid[i], 0.23*ymin, ymid[i+1], 0.18*ymin,
             col = rgb(strap[i], strap[i], strap[i]), border = NA)
      }
      #     draw intervals
      for (i in 1:nrow(vioIntvs)) {
        z = runif(1) # z = (i %% nlev) / nlev
        lines(y[vioIntvs$leftIndex[i]:vioIntvs$rightIndex[i]],
              rep(z*ymin*0.69+ymin*0.31, vioIntvs$rightIndex[i]-vioIntvs$leftIndex[i]+1),
              col = rgb(0.9, 0.9, 0.9))
      }
    }
    mtext("Intervals of violation", side = 2, line = 1, at = 0.7*ymin)
  }

  # return value
  vioIntvs
}

.hist2vec <- function(h, x) {
  if (missing(x)) 
    x = seq(min(h$breaks), max(h$breaks), length.out = sum(h$counts))
  else
    x = sort(x)
  n    = length(h)
  dval = numeric(n)
  dval[x >= h$breaks[1] & x <= h$breaks[2]] = h$density[1]
  if (length(h$breaks) > 2) {
    for (i in 2:(length(h$breaks)-1)) 
      dval[x > h$breaks[i] & x <= h$breaks[i+1]] = h$density[i]
  }
  data.frame("x" = x, "h" = dval)
}
