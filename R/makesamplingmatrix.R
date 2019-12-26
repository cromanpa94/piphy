#' Make a sampling matrix for the csv files in the current working directory

makesamplingmatrix<-function(folder=".", name="SM.csv"){
  options(warn=-1)
  files = list.files(folder, "*.csv", full.names = TRUE)
  cls<-lapply(files, function(x) read.csv(x, stringsAsFactors = FALSE))
  sampling_matrix <- Reduce(function(...) merge(..., by="name", all=TRUE), cls)
  names(sampling_matrix)<-c("Species", paste0("Cluster_",gsub("[^\\d]+", "", files, perl=TRUE)))
  write.csv(sampling_matrix, name)
  return(sampling_matrix)
  options(warn=-0)

}
