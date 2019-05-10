.pkgglobalenv <- new.env(parent = emptyenv())
# Get the quantile of the multiscale statistic via Monte Carlo simulations
msQuantile <- function(n, alpha = c(0.5), nsim = 5e3, is.sim = (n < 1e4), intv = genIntv(n),
                       mode = c("Con", "Gen"), ...) {
  mode <- match.arg(mode)
  if (!(length(n) == 1 && is.finite(n) && n > 1 && n == round(n))) 
    stop("'n' must be a single integer >= 2")
  if (!is.numeric(alpha) || any(alpha < 0) || any(alpha > 1) || length(alpha) < 1)
    stop("'alpha' (significance level) must be in [0, 1]!")
  if (is.sim) {
    if (!('left' %in% names(intv) && 'right' %in% names(intv)) || length(intv$left) != length(intv$right))
      stop("'intv' must have two fields 'left' and 'right', with equal length!")
    # validate intervals
    intv$left  = as.integer(intv$left)
    intv$right = as.integer(intv$right)
    vid = intv$left < intv$right & intv$left > 0 & intv$right <= n
    if (mode == "Gen")
      vid = vid & intv$left + 1 < intv$right
    if (!any(vid)) 
      stop("No valid intervals in 'intv'!")
    intv$left  = intv$left[vid]
    intv$right = intv$right[vid]
    intv       = intv[order(intv$left, intv$right),]
    # stat = .simQuantile(n,nsim,intv)
    if (exists('last_simulation',envir=.pkgglobalenv) == TRUE && 
        .pkgglobalenv$last_simulation$n == n && 
        .pkgglobalenv$last_simulation$nsim >= nsim && 
        nrow(.pkgglobalenv$last_simulation$intv) == nrow(intv) && 
        all(.pkgglobalenv$last_simulation$intv == intv) && 
        .pkgglobalenv$last_simulation$mode == mode) {
      message("Quantile is retrieved form the previous simulation!")
      stat = .pkgglobalenv$last_simulation$stat
    } else {
      message("Quantile simulation might take a while ... ", appendLF = FALSE)
      ##### fast C++ implementation
      stat = .msQuantile(as.integer(intv$left-1), as.integer(intv$right-1), n, nsim, (mode=="Gen"))
      assign("last_simulation", list('n'=n,'nsim'=nsim,'intv'=intv,'stat'=stat,'mode'=mode), 
             envir=.pkgglobalenv)
      message("... end!")
    }
  } else {
    message('Load precomputed quantiles ... ', appendLF=FALSE) 
    path = system.file("extdata", package = "essHist")
    fileName = switch(mode, Con="simn1e4c.rds", Gen="simn1e4g.rds")
    stat = readRDS(file.path(path, fileName))
    message('... end!') 
  }
  quantile(na.omit(stat), 1-alpha, names = FALSE, ...)
}

# # Debug: for simulation of Tn
# dms_stat <- function(n, nsim=5e3, isGen=TRUE) {
#   intv = genIntv(n)
#   .msQuantile(as.integer(intv$left-1), as.integer(intv$right-1), n, nsim, isGen)
# }

# # Codes for pre-simulation of quantiles
# n    = 1e4
# nsim = 1e5
# stat = dms_stat(n, nsim, FALSE)
# saveRDS(stat,file = "simn1e4c.rds")
# stat = dms_stat(n, nsim, TRUE)
# saveRDS(stat,file = "simn1e4g.rds")


# Compute the minimal threshold
.minThreshold <- function (n, intv) {
  theta = unique(intv$right - intv$left + (intv$left == 1))/n
  max(-sqrt(2 - 2*log (theta*(1-theta))))
} 
