queryPALSdb <- function(query, disp = c("data", "browser"),
                        field=c("keyword", "ug_id", "gb_id", "cluster_count")) {
  
  url.base <- "http://palsdb.ym.edu.tw/cgi-bin?"
  
  url <- paste(url.base, field, sep="&")
  
  if (disp == "data") {
    return(.handleXML(query))
  }
  else {
    browseURL(query)
  }
}

buildSpliceSites <- function(xml) {

  get.pData <- function(this.result, what) {
    leaf <- this.result[["AS-seq-info"]][["Lib-info"]][[what]]
    if (is.null(leaf)) {
      r <- NA
    } else {
      r <- as.character(xmlChildren(leaf)$text)[6]
    }
    return(r)
  }
  
  resultQuery <- xml$doc$children["PALSdbResultQuery"][[1]]

  entries.i <- which(names(resultQuery) == "PALSdbEntry")
  
  spsites.list <- vector("list", length=length(entries.i))
  names(spsites.list) <- rep(NULL, length(entries.i))
  
  for (i in entries.i) {
    ug.id <- resultQuery[[i]]["Head-info"][[1]][["ug-cluster-id"]]
    ug.id <- as.character(xmlChildren(ug.id)$text)[6]
    seq.len <- resultQuery[[i]][["Head-info"]][["ref-len"]]
    seq.len <- as.integer(as.character(xmlChildren(seq.len)$text)[6])
    
    i.ASinfo <- names(resultQuery[[i]]) == "AS-info"
    
    spsiteIpos <- vector("list", length=sum(i.ASinfo))
    spsiteIIpos <- vector("list", length=sum(i.ASinfo))
    spsiteIpos.pData.tissue <- vector("character", length=sum(i.ASinfo))
    spsiteIpos.pData.histology <- vector("character", length=sum(i.ASinfo))
    spsiteIpos.pData.cellline <- vector("character", length=sum(i.ASinfo))
    spsiteIIpos.pData.tissue <- vector("character", length=sum(i.ASinfo))
    
    for (j in which(i.ASinfo)) {

      this.result <- resultQuery[[i]][[j]]
      typeAS <- this.result[["Alter-info"]][["AS"]][["AS-type"]]
      typeAS <- as.integer(as.character(xmlChildren(typeAS)$text)[6])
      
      if (typeAS == 1) {
        pos1 <- this.result[["Alter-info"]][["AS"]][[2]]
        pos1 <- as.character(xmlChildren(pos1)$text)[6]
        pos2 <- this.result[["Alter-info"]][["AS"]][[3]]
        pos2 <- as.character(xmlChildren(pos2)$text)[6]
        spsiteIpos[[j]] <- c(as.integer(pos1),
                             as.integer(pos2))
        spsiteIpos.pData.tissue[j] <- get.pData(this.result, "Hit-tissue")
        spsiteIpos.pData.tissue[j] <- get.pData(this.result, "Hit-histology")
      }
      if (typeAS == 2) {
        pos2 <- resultQuery[[i]][[j]][["Alter-info"]][["AS"]][[2]]
        pos2 <- as.character(xmlChildren(pos2)$text)[6]
        spsiteIIpos[[j]] <- c(as.integer(pos2))
        spsiteIIpos.pData.tissue[j] <- get.pData(this.result, "Hit-tissue")
        spsiteIIpos.pData.tissue[j] <- get.pData(this.result, "Hit-histology")
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
    spsiteIpos.pData@pData <- data.frame(tissue = spsiteIpos.pData.tissue[spsiteIpos.i])
    
    spsites.list[[i]] <- new("SpliceSites", seq.len = seq.len,
                             spsiteIpos = spsiteIpos, spsiteIIpos = spsiteIIpos,
                             spsiteIpos.pData = spsiteIpos.pData
                             )
    
    names(spsites.list)[i] <- ug.id
  }
  
  return(spsites.list)
  
}
