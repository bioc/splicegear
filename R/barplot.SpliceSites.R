barplot.SpliceSites <- function(height, type.as=c("typeI", "typeII", "all"),
                                info="tissue", ...) {

  type.as <- match.arg(type.as)
  dt <- switch(type.as,
               typeI = pData(height@spsiteIpos.pData),
               typeII = pData(height@spsiteIIpos.pData),
               all = rbind(pData(height@spsiteIpos.pData), pData(height@spsiteIIpos.pData))
               )

  if (nrow(dt)>0 && ! info %in% colnames(dt))
    stop(paste(paste("parameter 'info=", info,"' is not one of ", sep=""),
               colnames(dt), collapse=", "))

  dt <- dt[[info]]

  if (is.null(dt))
    ct <- 0
  else
    ct <- table(dt)

  r <- barplot(ct, ...)

  invisible(r)
  
}
