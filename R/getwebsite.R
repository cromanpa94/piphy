#' This is not my function. I borrowed from somewhere

getwebsites<- function( websites ){
  url_2 <- websites
  doc_2 <- htmlParse(url_2)
  links_2 <- xpathSApply(doc_2, "//a/@href")
  todo<-grep("getalign",links_2, value=T )
  matches <- regmatches(todo, gregexpr("[[:digit:]]+", todo))
  num<-as.numeric(unlist(matches))
  newweb<-paste0("http://sirloinpope.com/cgi-bin/sql_getcluster_fasta.cgi?format=gi&db=194&ti=",num[1], "&cl=",num[2],"&ntype=1" )
}
