.First.lib <- function(libname, pkgname, where) {

  message <- TRUE

  if (missing(where)) {
    where <- match(paste("package:", pkgname, sep=""), search())
    if(is.na(where)) {
      warning(paste("Not a package name: ", pkgname))
      return()
    }
    where <- pos.to.env(where)
  }

  require(methods, quietly=TRUE) || stop("The package 'methods' is required !")
  require(Biobase, quietly=TRUE) || stop("The package 'Biobase' is required !")
  require(grid, quietly=TRUE) || stop("The package 'grid' is required !")

  if (message) {
    cat("splicegear loaded.\n")
    if ("1.3.25" > packageDescription("Biobase")[c("Version")]) {
      cat("Please source the function 'as.data.frame.exprSet' available at\n",
          "http://www.cbs.dtu.dk/laurent/download/splicegear/")
    }
  }

  if(.Platform$OS.type == "windows" && require(Biobase) && interactive()
        && .Platform$GUI ==  "Rgui"){
        addPDF2Vig("splicegear")
    }

}
