.pkgglobalenv <- new.env(parent = emptyenv())
# Get the quantile of the multiscale statistic via Monte Carlo simulations
msQuantile <- function(n, alpha = c(0.5), nsim = 5e3, is.sim = (n < 1e4), intv = NULL, ...) {
  if (!(length(n) == 1 && is.finite(n) && n > 1 && n == round(n))) 
    stop("'n' must be a single integer >= 2")
  if (!is.numeric(alpha) || any(alpha < 0) || any(alpha > 1) || length(alpha) < 1)
    stop("'alpha' (significance level) must be in [0, 1]!")
  if (is.null(intv))
    intv = genIntv(n)
  if (!('left' %in% names(intv) && 'right' %in% names(intv)) || length(intv$left) != length(intv$right))
    stop("'intv' must have two fields 'left' and 'right', with equal length!")
  # validate intervals
  intv$left  = as.integer(intv$left)
  intv$right = as.integer(intv$right)
  vid = intv$left < intv$right & intv$left > 0 & intv$right <= n
  if (!any(vid)) 
    stop("No valid intervals in 'intv'!")
  intv$left  = intv$left[vid]
  intv$right = intv$right[vid]
  intv       = intv[order(intv$left, intv$right),]
  
  if (is.sim) {
    # stat = .simQuantile(n,nsim,intv)
    if (exists('last_simulation',envir=.pkgglobalenv) == TRUE && 
        .pkgglobalenv$last_simulation$n == n && .pkgglobalenv$last_simulation$nsim >= nsim && 
        nrow(.pkgglobalenv$last_simulation$intv) == nrow(intv) && 
        all(.pkgglobalenv$last_simulation$intv == intv)) {
      message("Quantile is retrieved form the previous simulation!")
      stat = .pkgglobalenv$last_simulation$stat
    } else {
      message("Quantile simulation might take a while ... ", appendLF = FALSE)
      ##### fast C++ implementation
      stat = .msQuantile(as.integer(intv$left-1), as.integer(intv$right-1), n, nsim)
      assign("last_simulation", list('n'=n,'nsim'=nsim,'intv'=intv,'stat'=stat), 
             envir=.pkgglobalenv)
      message("... end!")
    }
  } else {
    message('Load precomputed quantiles ... ', appendLF=FALSE) 
    path = system.file("extdata", package = "essHist")
    stat = readRDS(file.path(path, "simn1e4.rds"))
    message('... end!') 
  }
  quantile(na.omit(stat), 1-alpha, names = FALSE, ...)
}

# # Debug: for simulation of Tn
# dms_stat <- function(n, nsim) {
#   intv = genIntv(n)
#   .msQuantile(as.integer(intv$left-1), as.integer(intv$right-1), n, nsim, TRUE)
# }

# Compute the minimal threshold
.minThreshold <- function (n, intv) {
  theta = unique(intv$right - intv$left + (intv$left == 1))/n
  max(-sqrt(2 - 2*log (theta*(1-theta))))
} 
