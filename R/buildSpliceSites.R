queryPALSdb <- function(query, disp = c("data", "browser"),
                        field = c("keyword", "ug.id", "gb.id",
                          "human.cytoband", "mouse.cytoband", "cluster_count"),
                        species = c("human", "mouse"),
                        e.value = "1e-1",
                        ident.threshold = c("90% 50b", "95% 50b", "90% 45b"),
                        verbose = FALSE) {

  require(annotate, quietly=TRUE) || stop()
  
  ##url.base <- "http://140.129.151.131/cgi-bin/palsdb/big_xml.cgi"
  url.base <- "http://140.129.151.155/~laurent/cgi-bin/palsdb/big_xml.cgi"

  disp <- match.arg(disp)
  field <- match.arg(field)
  species <- match.arg(species)
  ident.threshold <- match.arg(ident.threshold)
  
  format <- switch(disp,
                   data="xml",
                   browser="html")

  field <- switch(field,
                  keyword="A",
                  ug.id="B",
                  gb.id="C",
                  human.cytoband="E",
                  mouse.cytoband="F")

  query.tag <- "keyword"
  
  ident.threshold <- switch(ident.threshold,
                            "95% 45b"="A",
                            "90% 50b"="B",
                            "95% 50b"="C",
                            "90% 45b"="D"
                            )

  url <- paste(url.base,
               paste(paste("format", format, sep="="),
                     paste(query.tag, query, sep="="),
                     paste("field", field, sep="="),
                     paste("species", species, sep="="),
                     paste("evalue", e.value, sep="="),
                     paste("constraint", ident.threshold, sep="="),
                     "submit=Submit", ## ?!
                     sep="&"),
               sep="?")
  if (verbose) {
    print(url)
  }
  if (disp == "data") {
    require(XML, quietly=TRUE) || stop("Library XML required !")
    return(.handleXML(url))
    ##return(paste(readLines(url(url)), collapse=""))
  }
  else {
    browseURL(url)
  }
}

buildSpliceSites <- function(xml, verbose=TRUE) {

  get.pData <- function(this.result, what) {
    leaf <- this.result[["Lib-info"]][[what]]
    if (is.null(leaf)) {
      r <- NA
    } else {
      r <- as.character(xmlChildren(leaf)$text)[6]
    }
    return(r)
  }
  
  resultQuery <- xml$doc$children["PALSdbResultQuery"][[1]]

  entries.i <- which(names(resultQuery) == "PALSdbEntry")

  n.EST.hit <- sum(unlist(lapply(resultQuery[entries.i],
                                 function (x) length(names(x)) - 1)))
  
  spsites.list <- vector("list", length=length(entries.i))
  names(spsites.list) <- rep(NULL, length=length(entries.i))
  ##spsites.list <- vector("list", length=n.EST.hit)
  ##names(spsites.list) <- rep(NULL, n.EST.hit)

  ##i.list <- 1

  if (verbose)
    cat(length(entries.i), " entrie(s) in the set.\n")
  
  for (i in seq(along=entries.i)) {
    e.i <- entries.i[i]

    entries.k <- which(names(resultQuery[[e.i]]) == "AS-info")
    ug.id <- resultQuery[[e.i]][["Head-info"]][["ug-cluster-id"]]
    ug.id <- as.character(xmlChildren(ug.id)$text)[6]
    seq.len <- resultQuery[[e.i]][["Head-info"]][["ref-len"]]
    seq.len <- as.integer(as.character(xmlChildren(seq.len)$text)[6])

    n.ASinfo <- sum(unlist(lapply(resultQuery[[e.i]][entries.k],
                                 function(x) sum(names(x) == "AS-seq-info"))))
    spsiteIpos <- vector("list", length = n.ASinfo)
    spsiteIIpos <- vector("list", length = n.ASinfo)
    pData.tissue <- vector("character", length = n.ASinfo)
    pData.histology <- vector("character", length = n.ASinfo)
    pData.cellline <- vector("character", length = n.ASinfo)
    pData.site <- vector("character", length = n.ASinfo)
    j.offset <- 0
    
    if (verbose)
      cat(" entrie ", i, " has ", n.ASinfo, " element(s).\n")
  
    for (k in seq(along=entries.k)) {
      e.k <- entries.k[k]
      
      i.ASinfo <- names(resultQuery[[e.i]][[e.k]]) == "AS-seq-info"

      entries.j <- which(i.ASinfo)
      
      if (verbose)
        cat("  sub-entrie ", k, " has ", length(entries.j), " element(s).\n")

      for (j in seq(along=entries.j)) {
        ##cat ("j=", j, "(j.offset=", j.offset, ")\n")
        j.offset <- j.offset + 1
        e.j <- entries.j[j]
        this.result <- resultQuery[[e.i]][[e.k]]
        typeAS <- this.result[["Alter-info"]][["AS"]][["AS-type"]]
        typeAS <- as.integer(as.character(xmlChildren(typeAS)$text)[6])
        
        if (typeAS == 1) {
          pos1 <- this.result[["Alter-info"]][["AS"]][[2]]
          pos1 <- as.character(xmlChildren(pos1)$text)[6]
          pos2 <- this.result[["Alter-info"]][["AS"]][[3]]
          pos2 <- as.character(xmlChildren(pos2)$text)[6]
          spsiteIpos[[j.offset]] <- c(as.integer(pos1),
                               as.integer(pos2))
        }
        if (typeAS == 2) {
          pos2 <- resultQuery[[e.i]][[e.k]][["Alter-info"]][["AS"]][[2]]
          pos2 <- as.character(xmlChildren(pos2)$text)[6]
          spsiteIIpos[[j.offset]] <- c(as.integer(pos2))
        }
        pData.tissue[j.offset] <- get.pData(this.result[[e.j]], "Hit-tissue")
        pData.histology[j.offset] <- get.pData(this.result[[e.j]], "Hit-histology")
        pData.site[j.offset] <- k
      }

    }
    spsiteIpos.i <- ! unlist(lapply(spsiteIpos, is.null))
    
    if (sum(spsiteIpos.i) == 0)
      spsiteIpos <- matrix(0, 0, 0)
    else
      spsiteIpos <- matrix(unlist(spsiteIpos[spsiteIpos.i]), nc=2, byrow=TRUE)
    
    spsiteIIpos.i <- ! unlist(lapply(spsiteIIpos, is.null))
    
    if (sum(spsiteIIpos.i) == 0)
      spsiteIIpos <- integer(0)
    else
      spsiteIIpos <- unlist(spsiteIIpos[spsiteIIpos.i])
    
    ##spsiteIpos.pData <- new("phenoData", pData=data.frame(tissue = spsiteIpos.pData.tissue[spsiteIpos.i]))
    spsiteIpos.pData <- new("phenoData")
    spsiteIpos.pData@pData <-
      data.frame(tissue = pData.tissue[spsiteIpos.i],
                 histology = pData.histology[spsiteIpos.i],
                 site = pData.site[spsiteIpos.i])
    spsiteIIpos.pData <- new("phenoData")
    spsiteIIpos.pData@pData <-
      data.frame(tissue = pData.tissue[spsiteIIpos.i],
                 histology = pData.histology[spsiteIIpos.i],
                 site = pData.site[spsiteIIpos.i])
    
    spsites.list[[i]] <- new("SpliceSites", seq.len = seq.len,
                             spsiteIpos = spsiteIpos,
                             spsiteIIpos = spsiteIIpos,
                             spsiteIpos.pData = spsiteIpos.pData,
                             spsiteIIpos.pData = spsiteIIpos.pData
                             )
    
    names(spsites.list)[i] <- ug.id
    ##i.list <- i.list + 1
  }
  
  return(spsites.list)
  
}


# group.pdata <- function(spsites, type.as) {

#   if (type.as == 1) {
#     u <- unique(spsites@spsiteIpos)
#     pdat <- pData(spsites@spsiteIpos.pData)
#     pos <- spsites@spsiteIpos
#   } else if (type.as ==2) {
#     pos <- matrix(spsites@spsiteIIpos, nc = 1)
#     u <- unique(spsites@spsiteIIpos)
#     pdat <- pData(spsites@spsiteIIpos.pData)
#   }
  
#   u.n <- nrow(u)
  
#   r <- vector("list", length=u.n)
#   names(r) <- apply(u, 1, paste, collapse="..")

#   for (r.i in seq(1, u.n, length=u.n)) {
#     i <- apply(pos, 1, function(x) all(x == u[r.i, ]))
#     r[[r.i]] <- pdat[i, ]
#   }

#   return(r)
# }
