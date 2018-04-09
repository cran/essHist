# Smooth discrete data
smData <- function(x, sw = NA) {
  ret   = sort(x)
  dupId = duplicated(ret)
  if (any(dupId)) {
    if (is.na(sw)) sw = max(1e-12, 0.02*min(diff(ret[!dupId])))
    ret[dupId] = ret[dupId] + sw*rnorm(sum(dupId))
  } 
  sort(ret)
}