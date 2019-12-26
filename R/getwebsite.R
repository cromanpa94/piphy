#' This is not my function. I borrowed from somewhere

getwebsites<- function( websites ){
  gsub("sql_getcluster.cgi","sql_getcluster_fasta.cgi"  , websites)
}
