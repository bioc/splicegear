##.initSpliceSitesGenomicMethods <- function(where=where) {

    setClass("SpliceSitesGenomic",
             representation(variants = "list"),
             contains="SpliceSites")
##             where=where)

    setMethod("plot", signature(x="SpliceSitesGenomic", y="missing"),
              function(x, ...) {
                plot.SpliceSitesGenomic(x, ...)
              })
##, where = where)
    
##  }


##a <- new("SpliceSitesGenomic", seq.length=as.integer(10), spsiteIpos=matrix(c(1, 3, 5, 2, 3.5, 8), nc=2),
##         variants=list(c(1,2,3), c(2,3)))
