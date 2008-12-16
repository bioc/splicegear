as.data.frame.SpliceSites <- function(x, row.names=NA, optional=NA, ...) {

    if (! is(x, "SpliceSites"))
        stop("the argument should inherit class 'SpliceSites'")

  nr.typeII <- length(x@spsiteIIpos)

  pdat.I <- x@spsiteIpos.pData@pData
  rownames(pdat.I) <- seq(1, nrow(pdat.I), length=nrow(pdat.I))
  pdat.II <- x@spsiteIIpos.pData@pData
  rownames(pdat.II) <- seq(nrow(pdat.I)+1, nrow(pdat.I) + nr.typeII, length=nr.typeII)

  if (nrow(pdat.I) == 0 && nrow(pdat.II) == 0) {
    m <- list()
  } else if (nrow(pdat.I) == 0) {
    m <- pdat.II
  } else if (nrow(pdat.II) == 0) {
    m <- pdat.I
  } else {
    m <- merge(pdat.I, pdat.II, by.x = c(0, seq(along=pdat.I)), by.y = c(0, seq(along=pdat.II)), all = TRUE)
  }

  rv <- do.call(data.frame, c(list(begin = c(x@spsiteIpos[, 1], x@spsiteIIpos),
                                     end = c(x@spsiteIpos[, 2], rep(NA, nr.typeII)),
                                     lapply(m, function(x) {x})
                                     ))
                )


  return(rv)
}
