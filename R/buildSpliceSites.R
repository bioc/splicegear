buildSpliceSites <- function(xml, verbose=TRUE) {

  get.pData <- function(this.result, what) {
    leaf <- this.result[[what]]
    if (is.null(leaf)) {
      r <- NA
    } else {
      r <- as.character(xmlChildren(leaf)$text)[6]
    }
    return(r)
  }
  
  resultQuery <- xml$doc$children[["ResultQuery"]]
  
  entries.i <- which(names(resultQuery) == "Entry")

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

    entries.k <- which(names(resultQuery[[e.i]]) == "Alt-Splice")
    ug.id <- resultQuery[[e.i]][["Reference-sequence"]][["ug-cluster-id"]]
    ug.id <- as.character(xmlChildren(ug.id)$text)[6]
    seq.len <- resultQuery[[e.i]][["Reference-sequence"]][["ref-len"]]
    seq.len <- as.integer(as.character(xmlChildren(seq.len)$text)[6])

    n.ASinfo <- sum(unlist(lapply(resultQuery[[e.i]][entries.k],
                                 function(x) sum(names(x) == "Hit-info"))))
    spsiteIpos <- vector("list", length = n.ASinfo)
    spsiteIIpos <- vector("list", length = n.ASinfo)
    pData.tissue <- vector("character", length = n.ASinfo)
    pData.histology <- vector("character", length = n.ASinfo)
    pData.cellline <- vector("character", length = n.ASinfo)
    pData.site <- vector("character", length = n.ASinfo)
    ##pData.other <- vector("list", length=n.ASinfo)
    j.offset <- 0
    
    if (verbose)
      cat(" entrie ", i, " has ", n.ASinfo, " element(s).\n")
  
    for (k in seq(along=entries.k)) {
      e.k <- entries.k[k]
      
      i.ASinfo <- names(resultQuery[[e.i]][[e.k]]) == "Hit-info"

      entries.j <- which(i.ASinfo)

      this.result <- resultQuery[[e.i]][[e.k]]
      typeAS <- xmlAttrs(this.result)["Type"]
      typeAS <- as.integer(typeAS)

      if (typeAS == 1) {
        pos.i <- which(names(this.result[["Site-info"]]) == "Pos")
        pos1 <- this.result[["Site-info"]][[pos.i[1]]]
        pos1 <- as.character(xmlChildren(pos1)$text)[6]
        pos2 <- this.result[["Site-info"]][[pos.i[2]]]
        pos2 <- as.character(xmlChildren(pos2)$text)[6]
      }
      if (typeAS == 2) {
        pos2 <- resultQuery[[e.i]][[e.k]][["Site-info"]][["Pos"]]
        pos2 <- as.character(xmlChildren(pos2)$text)[6]
      }
      
      if (verbose)
        cat("  sub-entrie ", k, " has ", length(entries.j), " element(s).\n")

      for (j in seq(along=entries.j)) {
        ##cat ("j=", j, "(j.offset=", j.offset, ")\n")
        j.offset <- j.offset + 1
        e.j <- entries.j[j]
        if (typeAS == 1) {
          spsiteIpos[[j.offset]] <- c(as.integer(pos1),
                                      as.integer(pos2))
        }
        if (typeAS == 2) {
          spsiteIIpos[[j.offset]] <- c(as.integer(pos2))
        }
        pData.tissue[j.offset] <- get.pData(this.result[[e.j]], "Hit-tissue")
        pData.histology[j.offset] <- get.pData(this.result[[e.j]], "Hit-histology")
        
        ##pData.other[[j.offset]] <- other.pdata
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
    
    spsites.list[[i]] <- new("SpliceSites", seq.length = seq.len,
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
