queryPALSdb <- function(query, disp = c("data", "browser"),
                        field = c("keyword", "ug.id", "gb.id",
                          "human.cytoband", "mouse.cytoband", "cluster_count"),
                        species = c("human", "mouse"),
                        e.value = "1e-1",
                        ident.threshold = c("90% 50b", "95% 50b", "90% 45b"),
                        verbose = FALSE) {

  require(annotate, quietly=TRUE) || stop()
  
  ##url.base <- "http://palsdb.ym.edu.tw/cgi-bin?"
  url.base <- "http://140.129.151.155/~cedric/cgi-bin/big_xml.cgi"

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
  
  for (i in seq(along=entries.i)) {
    e.i <- entries.i[i]
    ug.id <- resultQuery[[e.i]]["Head-info"][[1]][["ug-cluster-id"]]
    ug.id <- as.character(xmlChildren(ug.id)$text)[6]
    seq.len <- resultQuery[[e.i]][["Head-info"]][["ref-len"]]
    seq.len <- as.integer(as.character(xmlChildren(seq.len)$text)[6])
    
    i.ASinfo <- names(resultQuery[[e.i]]) == "AS-info"
    
    spsiteIpos <- vector("list", length=sum(i.ASinfo))
    spsiteIIpos <- vector("list", length=sum(i.ASinfo))
    pData.tissue <- vector("character", length=sum(i.ASinfo))
    pData.histology <- vector("character", length=sum(i.ASinfo))
    pData.cellline <- vector("character", length=sum(i.ASinfo))
        
    for (j in which(i.ASinfo)) {

      this.result <- resultQuery[[e.i]][[j]]
      typeAS <- this.result[["Alter-info"]][["AS"]][["AS-type"]]
      typeAS <- as.integer(as.character(xmlChildren(typeAS)$text)[6])
      
      if (typeAS == 1) {
        pos1 <- this.result[["Alter-info"]][["AS"]][[2]]
        pos1 <- as.character(xmlChildren(pos1)$text)[6]
        pos2 <- this.result[["Alter-info"]][["AS"]][[3]]
        pos2 <- as.character(xmlChildren(pos2)$text)[6]
        spsiteIpos[[j]] <- c(as.integer(pos1),
                             as.integer(pos2))
      }
      if (typeAS == 2) {
        pos2 <- resultQuery[[e.i]][[j]][["Alter-info"]][["AS"]][[2]]
        pos2 <- as.character(xmlChildren(pos2)$text)[6]
        spsiteIIpos[[j]] <- c(as.integer(pos2))
      }
      pData.tissue[j] <- get.pData(this.result, "Hit-tissue")
      pData.histology[j] <- get.pData(this.result, "Hit-histology")
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
                 histology = pData.histology[spsiteIpos.i])
     spsiteIIpos.pData <- new("phenoData")
    spsiteIIpos.pData@pData <-
      data.frame(tissue = pData.tissue[spsiteIIpos.i],
                 histology = pData.histology[spsiteIIpos.i])
    
    spsites.list[[i]] <- new("SpliceSites", seq.len = seq.len,
                             spsiteIpos = spsiteIpos,
                             spsiteIIpos = spsiteIIpos,
                             spsiteIpos.pData = spsiteIpos.pData,
                             spsiteIIpos.pData = spsiteIIpos.pData
                             )
    
    names(spsites.list)[i] <- ug.id
  }
  
  return(spsites.list)
  
}
