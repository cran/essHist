# Generate system of intervals of form (a, b] with a < b
genIntv <- function (n, type = c('Sparse', 'Full')) {
  if (!(is.finite(n) && n > 1 && n == round(n) && length(n) == 1)) 
    stop("'n' should be an integer >= 2")
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
  unique(data.frame('left'=left,'right'=right))
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
    ind.valid   = (left %in% left.valid) & (right %in% right.valid)# & (y[left] < y[right])
    if (!all(ind.valid)) {
      itEnd  = unique(data.frame('leE'=y[left[!ind.valid]],'riE'=y[right[!ind.valid]]))
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
      idMd = (leMd < riMd)# & (y[leMd] < y[riMd])
      leMd = leMd[idMd]
      riMd = riMd[idMd]
      left  = c(left[ind.valid],leMd)
      right = c(right[ind.valid],riMd)
    }
    # left  = left[ind.valid]
    # right = right[right.valid]
  }
  # remove repetitive intervals
  unique(data.frame('left'=left, 'right'=right))
}