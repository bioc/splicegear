require(Biobase, quietly=TRUE) || stop("Could not load the package 'Biobase'")

  
setClass("SpliceSites",
         representation(
                        seq = "character",   # the genomic sequence (if available)
                        seq.length = "integer", 
                        spsiteIpos = "matrix", # a two-columns matrix (window for type I splice site),
                          spsiteIIpos = "integer", # position for type II splice site
                        spsiteIIIpos = "matrix", # a two-columns matrix (window for type III splice site),
                        spsiteIpos.pData = "phenoData",
                        spsiteIIpos.pData = "phenoData",
                        spsiteIIIpos.pData = "phenoData"
                        ))


  setMethod("show", signature(object = "SpliceSites"),
            function(object) {
              cat("Alternative splicing sites (SpliceSites):\n")
              cat("\tseq is ", object@seq.length, " bp long", sep="")
              if (object@seq == "")
                cat(" (warning: sequence not included).\n")
              else
                cat(".\n")
              cat("\t", nrow(object@spsiteIpos), " type I splice site(s)\n", sep="")
              show(object@spsiteIpos.pData)
              cat("\t", length(object@spsiteIIpos), " type II splice site(s)\n", sep="")
              show(object@spsiteIIpos.pData)
              cat("\t", nrow(object@spsiteIIIpos), " type III splice site(s)\n", sep="")
            })


if( !isGeneric("plot") )
  setGeneric("plot", function(x, y, ...)
             standardGeneric("plot"))


setMethod("plot",
          signature(x="SpliceSites", y="missing"),
          function(x, ...) {
              plot.SpliceSites(x, ...)
            })


setMethod("initialize", "SpliceSites",
          function(.Object, 
                   seq = "", seq.length = as.integer(-1),
                   spsiteIpos = matrix(0, 0, 0),
                   spsiteIIpos = integer(0),
                   spsiteIIIpos = matrix(0, 0, 0),
                   spsiteIpos.pData = new("phenoData"),
                   spsiteIIpos.pData = new("phenoData"),
                   spsiteIIIpos.pData = new("phenoData"),
                   ...)
          {
            if (seq == "" && seq.length == -1)
              stop("'seq' or 'seq.length' must be defined.")
            .Object <- callNextMethod()
            .Object@seq = seq
            .Object@seq.length = seq.length
            .Object@spsiteIpos = spsiteIpos
            .Object@spsiteIIpos = spsiteIIpos
            .Object@spsiteIIIpos = spsiteIIIpos
            .Object@spsiteIpos.pData = spsiteIpos.pData
            .Object@spsiteIIpos.pData = spsiteIIpos.pData
            .Object@spsiteIIIpos.pData = spsiteIIIpos.pData
            return(.Object)
          })


if( !isGeneric("grid.plot") )
  setGeneric("grid.plot", function(x, y, ...)
             standardGeneric("grid.plot"))

setMethod("grid.plot",
          signature(x="SpliceSites", y="missing"),
          function(x, ...) {
              grid.plot.SpliceSites(x, ...)
            })
