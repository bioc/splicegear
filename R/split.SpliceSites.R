split.SpliceSites <- function(x, f=list(typeI=NA, typeII=NA), drop=NULL, ...) {

  get.factor <- function(type.as, f) {
    
    dt <- switch(type.as,
                 typeI = pData(x@spsiteIpos.pData),
                 typeII = pData(x@spsiteIIpos.pData)
                 )

    f <- factor(numeric(0))
    if (is.character(f) && length(f) == 1 && nrow(dt) > 0) {
      if (! f %in% colnames(dt))
        stop(paste(paste("parameter 'f=", f, sep=""),
             "' is not one of ", colnames(dt), collapse=", "))
      else
        f <- dt[f]
    }
    return(f)
  }
  
  get.pdata <- function(type.as, f) {
    if (length(f) == 0) {
      dt <- list()
    } else {
      dt <- switch(type.as,
                   typeI = cbind(x@spsiteIpos, pData(x@spsiteIpos.pData)),
                   typeII = cbind(x@spsiteIIpos, pData(x@spsiteIIpos.pData))
                   )
      
      dt <- split(dt, f)
    }
    
    return(dt)
  }
  
  if (! is.factor(f$typeI)) {
    if (is.na(f$typeI))
      f$typeI <- get.factor("typeI", "site")
    else
      f$typeI <- get.factor("typeI", f$typeI)
  }
  
  if (! is.factor(f$typeII)) {
    if (is.na(f$typeII))
      f$typeII <- get.factor("typeII", "site")
    else
      f$typeII <- get.factor("typeII", f$typeII)
  }
  
  
  r.I <- get.pdata("typeI", f$typeI)
  r.II <- get.pdata("typeII", f$typeII)
  
  r <- vector("list", length=length(r.I) + length(r.II))

  for (i in seq(along=r.I)) {
    r[[i]] <- new("SpliceSites", seq = x@seq, seq.length = x@seq.length,
                  spsiteIpos = as.matrix(r.I[[i]][, 1:2]),
                  spsiteIpos.pData = new("AnnotatedDataFrame", pData=r.I[[i]][, -c(1,2)],
                    varLabels=as.list(names(r.I[[i]])[-c(1,2)]))
                  )
  }

  for (i in seq(along=r.II)) {
    r[[i+length(r.I)]] <- new("SpliceSites", seq = x@seq, seq.length = x@seq.length,
                              spsiteIIpos = as.integer(as.vector(r.II[[i]][, 1])),
                              spsiteIIpos.pData = new("AnnotatedDataFrame", pData=r.II[[i]][, -c(1)],
                                varLabels=as.list(names(r.II[[i]])[-c(1)]))
                              )
  }
  
  return(r)
}
