.First.lib <- function(libname, pkgname, where) {

  message <- TRUE
  
  if (missing(where)) {
    where <- match(paste("package:", pkgname, sep=""), search())
    if(is.na(where)) {
      warning(paste("Not a package name: ",pkgname))
      return()
        }
    where <- pos.to.env(where)
  }

  if (message) {
    cat("splicegear: development version.\n")
    cat("Make sure you have the latest devel version of the package Biobase.\n")
  }
  
  require(methods, quietly=TRUE) || stop("The package 'methods' is required !")
  require(Biobase, quietly=TRUE) || stop("The package 'Biobase' is required !")

  .initSpliceSitesMethods(where=where)
  .initSpliceExprSetMethods(where=where)
  .initProbes(where=where)
  ##where <- match(paste("package:", pkgname, sep=""), search())
  ##where <- as.environment(paste("package:", pkgname, sep=""))

  ##cacheMetaData(as.environment(where))
}
