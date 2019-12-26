#' Retrieve orthologous clusters from genbank for a given set of ingroup and
#' outgroup taxa. This function also performs multiple sequence alignment
#' and runs Aliscore on the resulting sequences.
#'
#' @param ingroup The name of a target clade
#' @param MSA Whether or not perform sequence alignment on the retrieved clusters
#' @param ALI Whether or nor run Aliscore on the aligned sequences
#'
#'
#' @example
#' orthoclusters_genbank("centrolene", "cochranella", MSA=T, ALI=T)
#'
#' @export

orthoclusters_genbank <- function(ingroup,
                               outgroup,
                               MSA = F,
                               ALI = F) {
  ##Delete files
  if (file.exists("sequences.fasta")) {
    unlink("sequences.fasta", recursive = T)
  } else{
  }
  if (length(grep("cluster_*", list.files("."), value = T)) == 0) {

  } else {
    unlink(grep("cluster_*", list.files("."), value = T))
  }
  if(ALI == T & MSA ==F ){stop("Aliscore is only available when clusters are aligned" )}
  if (length(setdiff("msa", rownames(installed.packages()))) > 0 & MSA==T ) {
    stop("Please install msa first")
  }
  if (length(setdiff("ips", rownames(installed.packages()))) > 0 & ALI==T ) {
    stop("Please install ips first")
  }
  if (is.null(ingroup)) {
    stop("Please select an ingroup")
  }
  if (class(ingroup) != "character") {
    stop("Please use a character vector for the ingroup")
  }
  if (file.exists("Unaligned")) {
    stop("Please work in a new working directory or delete the existing files")
  }


  ##GetWD
  mainDir <-  getwd()
  options(warn = -1)
  dir.create(file.path(mainDir, "Unaligned"))
  options(warn = -0)
  setwd(file.path(mainDir, "Unaligned"))

  taxa <- c(ingroup, outgroup)

  cat("Downloading sequences")

  lapply(1:length(taxa), function(x) {
    tr <-
      entrez_search(
        db = "nuccore",
        term = paste0(taxa[x], "[ORGN]"),
        use_history = TRUE
      )
    for (seq_start in seq(1, tr$count, tr$count / 10)) {
      recs <- entrez_fetch(
        db = "nuccore",
        web_history = tr$web_history,
        rettype = "fasta",
        retstart = seq_start
      )
      cat(recs, file = "sequences.fasta", append = TRUE)
      cat(round(seq_start + 49), "sequences downloaded\r")
    }

  })

  seqs <- ape::read.dna("sequences.fasta", "fasta")

  cat("Removing sequences for unconfirmed species")


  seqs <- seqs[!duplicated(names(seqs))]
  seqs <- seqs[-grep("sp.|cf.|aff.|UNVERIFIED:", names(seqs))]

  write.dna(seqs, "sequences.fasta", format = "fasta")

  cat("Blasting sequences")


  args <-
    paste0("-in ",
           getwd(),
           "/sequences.fasta",
           " -input_type fasta -dbtype nucl")
  system2(
    command = 'makeblastdb',
    args = args,
    wait = FALSE,
    stdout = TRUE
  )


  que <-
    paste0(
      "-outfmt 6 -query sequences.fasta -out test.txt -db ",
      getwd(),
      "/sequences.fasta",
      " -num_threads 4"
    )
  system2(
    command = 'blastn',
    args = que,
    wait = FALSE,
    stdout = TRUE
  )


  BLAST1 <- read.blast(file = "test.txt")
  sim_mat <- simMatrix(BLAST1)

  res.d <- sim2dist(sim_mat)
  hc <- hclust(res.d)
  hc2 <- cutree(hc, h = 0.999999)

  ##Write clusters
  cat(length(unique(hc2)), "clusters found")

  seqs <- ape::read.dna("sequences.fasta", "fasta")
  lapply(1:length(unique(hc2)), function(x) {
    names <- sapply(strsplit(names(seqs), " "), head, 1)
    seqs2 <- seqs[which(names %in% names(which(hc2 == x)))]

    cat(
      file = paste0("cluster_", x, ".fasta"),
      paste(
        paste0(">", names(seqs2)),
        sapply(as.character(seqs2), paste, collapse =
                 ""),
        sep = "\n"
      ),
      sep = "\n"
    )

  })

  ##Sampling

  seqs <- list.files(".", pattern = "cluster")

  for (i in 1:length(seqs)) {
    sequences <- readDNAStringSet(seqs[i])
    names(sequences) <- gsub(" ", "_", names(sequences))
    sequence_names <- cbind.data.frame(name = word(names(sequences), 2, 3, fixed("_")),
                                       AN = word(names(sequences), 1, 1, fixed("_")))

    names(sequences) <- sequence_names$name

    summary_unique_sequences <-
      sequence_names[!duplicated(sequence_names$name), ]

    sequences <- sequences[!duplicated(sequence_names$name)]

    anys <- grep("_sp|_aff", summary_unique_sequences$name)
    if (length(anys) > 0) {
      summary_unique_sequences <- summary_unique_sequences[-anys, ]
      sequences <- sequences[!anys]
    }

    if (length(sequences) > 2) {
      write.dna(
        as.DNAbin(sequences),
        file = paste0(ingroup, "_", "Cluster" , i, ".fasta"),
        format = "fasta"
      )

      write.csv(
        summary_unique_sequences[c("name", "AN")],
        file = paste0(ingroup, "_", "Cluster" , i, ".csv") ,
        row.names = F
      )
    }
    plot.progress(i / length(seqs))

  }


  #Delete all files in the current directory

  file.remove(
    list.files(
      ".",
      pattern = "cluster|sequences",
      include.dirs = F,
      full.names = T,
      recursive = T
    )
  )

  sm <-
    makesamplingmatrix(".", name = "Sampling_matrix_unaligned.csv")

  setwd(mainDir)

  if (MSA == TRUE) {
    alignments<-align_piphy(mainDir = mainDir, clade=ingroup)
    write.csv(sm, "Sampling_matrix_aligned.csv")
    if (ALI == TRUE) {
      setwd(mainDir)
      sm<-ali_piphy(mainDir, alignments, sm=sm, clade=ingroup)
      write.csv(sm, "Sampling_matrix_alignedAli.csv")
      setwd(mainDir)
      unlink(list.dirs()[grep("Aliscore_", list.dirs(), ignore.case = T)], recursive =
               TRUE)
    }
  }
}
