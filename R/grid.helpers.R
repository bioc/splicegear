# grid.numeric2npc <- function(x, xlim=NULL, lower.blank=0, upper.blank=0) {
#   if (is.null(xlim)) {
#     offset <- min(x, na.rm=TRUE)
#     scale <- max(x, na.rm=TRUE)
#   } else {
#     offset <- min(xlim)
#     scale <- max(xlim, na.rm=TRUE)
#   }
#   scale <- scale - offset
#   offset <- offset + lower.blank
#   scale <- scale + upper.blank
#   x <- (x - offset) / scale
#   return(x)
# }

grid.make.numeric2npc <- function(x=NULL, xlim=NULL, lower.blank=0, upper.blank=0) {
  if (! is.null(x)) {
    offset <- min(x, na.rm=TRUE)
    scale <- max(x, na.rm=TRUE)
  } else if (! is.null(xlim)) {
    offset <- min(xlim)
    scale <- max(xlim, na.rm=TRUE)
  } else {
    stop("Only 'x' or 'xlim' are expected to be defined !")
  }
  scale <- scale - offset
  offset <- offset + lower.blank
  scale <- scale + upper.blank
  
  f <- function(x) {
    x <- (x - offset) / scale
    return(x)
  }
  
  return(f)
}

# grid.expand.gp <- function(n, parlist=list(), ...) {
#   parlist <- c(parlist, substitute(list(...)))
#   lapply(parlist, function(x) rep(x, length=n))
# }

grid.expand.gp <- function(n, parlist=list()) {
  gp <- vector("list", length=n)
  parlist <- lapply(parlist, function(x) rep(x, length=n))
  for (i in seq(1, n, length=n)) {
    gp[[i]] <- do.call("gpar", lapply(parlist, function(x) x[i]))
  }
  return(gp)
}
