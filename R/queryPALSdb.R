getPALSdbURL <- function(query, disp = c("data", "browser"),
                        field = c("keyword", "ug.id", "gb.id",
                          "human.cytoband", "mouse.cytoband", "cluster_count"),
                        species = c("human", "mouse"),
                        e.value = "1e-1",
                        ident.threshold = c("90% 50b", "95% 50b", "90% 45b")) {
  
  url.base <- "http://palsdb.ym.edu.tw/cgi-bin/palsdb/big_xml.cgi"
  ##url.base <- "http://140.129.151.155/~laurent/cgi-bin/palsdb/big_xml.cgi"

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
}

queryPALSdb <- function(query, disp = c("data", "browser"),
                        field = c("keyword", "ug.id", "gb.id",
                          "human.cytoband", "mouse.cytoband", "cluster_count"),
                        species = c("human", "mouse"),
                        e.value = "1e-1",
                        ident.threshold = c("90% 50b", "95% 50b", "90% 45b")) {

  disp <- match.arg(disp)
  
  url <- getPALSdbURL(query=query, disp = disp,
                      field = field,
                      species = species,
                      e.value = e.value,
                      ident.threshold = ident.threshold)
  
  if (disp == "data") {
    return(.handleXML(url))
    ##return(paste(readLines(url(url)), collapse=""))
  }
  else {
    browseURL(url)
  }
}

