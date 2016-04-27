# Zero-truncated poisson distribution

rpoisson0 <- function(n, lambda){
  U <- runif(n)   # the uniform sample
  t <- -log(1 - U*(1 - exp(-lambda))) # the "first" event-times
  T1 <- (lambda - t)   # the set of (T-t)
  
  X <- rpois(n, T1) + 1 # the final truncated Poisson sample
  return(X)
}