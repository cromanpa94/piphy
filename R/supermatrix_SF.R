supermatrix_SF<-function (missing = "-", prefix = "concatenated", save = T, threshold=0,
                        pattern="Cluster") {

  if(is.null(threshold)){
    file.names <- list.files(pattern=".fasta")
    DNA <- list()
    for (i in 1:length(file.names)) {
      print(paste("Reading alignment", i))
      DNA[[i]] <- ape::read.dna(file = file.names[i], format = "f",
                           as.character = T)
    }

  }else{
    file.names <- list.files(pattern=".fasta")
    file.names <- grep(pattern, file.names, value = T)

    TAL<- sapply(1:length(file.names), function(x){
      str<-file.names[x]
      sf1<-as.numeric(unlist(regmatches(str,gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*",str, perl=TRUE))))[2]
      sf1>threshold
    })
    file.names<-file.names[TAL==T]

    DNA <- list()
    for (i in 1:length(file.names)) {
      print(paste("Reading alignment", i))
      DNA[[i]] <- read.dna(file = file.names[i], format = "f",
                           as.character = T)
    }

  }


  total.length <- 0
  for (i in 1:length(DNA)) {
    total.length <- total.length + ncol(DNA[[i]])
  }
  taxa <- vector()
  for (i in 1:length(DNA)) {
    counter <- length(taxa) + 1
    for (j in 1:nrow(DNA[[i]])) {
      taxa <- c(taxa, row.names(DNA[[i]])[!row.names(DNA[[i]]) %in%
                                            taxa])
    }
  }
  seqmatrix <- matrix(as.character(missing), length(taxa),
                      total.length)
  rownames(seqmatrix) <- taxa
  partitions <- as.data.frame(matrix(, length(DNA), 3))
  colnames(partitions) <- c("part", "start", "stop")
  c.col <- 0
  print("Creating supermatrix")
  for (i in 1:length(DNA)) {
    print(paste("Processing alignment", i))
    gene <- DNA[[i]]
    print(i)
    for (j in 1:nrow(gene)) {
      c.row <- which(rownames(seqmatrix) == rownames(gene)[j])
      seqmatrix[c.row, (c.col + 1):(c.col + ncol(gene))] <- gene[j,
                                                                 ]
    }
    partitions[i, 1:3] <- c(file.names[i], c.col + 1, c.col +
                              ncol(gene))
    c.col <- c.col + ncol(gene)
  }
  results <- list()
  results[[1]] <- partitions
  results[[2]] <- seqmatrix
  if (save == T) {
    print("saving files")
    write.dna(seqmatrix, file = paste(prefix, threshold,".fasta",
                                      sep = ""), format = "f")
    write.csv(partitions, row.names = F, file = paste(prefix,threshold,
                                                      ".partitions.csv", sep = ""))
  }
  return(results)
}
