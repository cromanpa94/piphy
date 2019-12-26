#' This function downloads all the  ortholog clusters from the Phylota
#' pipeline for any given clade in the database. Additional subrutines are available
#' to align all downloaded sequences and perform aaliscore.
#'
#' @param clade The name of a valid taxonomic group with sequences in Phylota
#' @param MSA Whether or not perform sequence alignment on the retrieved clusters
#' @param ALI Whether or nor run Aliscore on the aligned sequences
#'
#'
#' @example
#' get_phylota("centrolene", MSA = T, ALI =T)
#'
#' @export

get_phylota <- function(clade,
                        MSA = F,
                        ALI = F) {
  if (is.null(clade)) {
    stop("Please select a Clade")
  }
  if (class(clade) != "character") {
    stop("Please use a character vector for Clade")
  }
  if (file.exists("Unaligned")) {
    stop("Please work in a new working directory or delete the existing files")
  }
  if(ALI == T & MSA ==F ){stop("Aliscore is only available when clusters are aligned" )}
  if (length(setdiff("msa", rownames(installed.packages()))) > 0 & MSA==T ) {
    stop("Please install msa first")
  }
  if (length(setdiff("ips", rownames(installed.packages()))) > 0 & ALI==T ) {
    stop("Please install ips first")
  }

  ##GetWD
  mainDir <-  getwd()
  setwd(mainDir)
  ti <- get_uid(sciname = clade)[1]
  ##Get Main website
  url <-
    paste0(
      "http://sirloinpope.com/cgi-bin/sql_getclusterset.cgi?ti=",
      ti,
      "&ntype=1&piflag=1&dflag=0&db=194"
    )
  doc <- htmlParse(url)
  links <- xpathSApply(doc, "//a/@href")
  clusters <- grep("getcluster", links, value = T)
  clusters <- clusters[grep(ti, clusters)]
  cat("Retrieving clusters")
  seqs <- pblapply(clusters, getwebsites)
  options(warn = -1)
  dir.create(file.path(mainDir, "Unaligned"))
  options(warn = -0)

  setwd(file.path(mainDir, "Unaligned"))

  cat("Downloading clusters")
  for (i in 1:length(seqs)) {
    sequences <- ReadFasta(gsub("format=gi", "format=all", seqs[i]$href))
    sequences$name <- gsub(" ", "_", sequences$name)
    sequences$gi <- word(sequences$name, 2, 2, fixed("|"))
    sequences$name <- word(sequences$name, 1, 2, fixed("_"))
    sequences$name <-
      sapply(strsplit(sequences$name, "|", fixed = T), `[`, 3)

    summary_unique_sequences <-
      sequences[!duplicated(sequences$name),]
    anys <- grep("_sp|_aff", summary_unique_sequences$name)
    if (length(anys) > 0) {
      summary_unique_sequences <- summary_unique_sequences[-anys,]
    }

    seqRFLP::write.fasta(
      dataframe2fas(summary_unique_sequences[c("name", "sequence")]),
      file = paste0(clade, "_", "Cluster" , i, ".fasta")
    )
    write.csv(
      summary_unique_sequences[c("name", "gi")],
      file = paste0(clade, "_", "Cluster" , i, ".csv") ,
      row.names = F
    )
    plot.progress(i / length(seqs))

  }

  cat("\nDone!!")

  ##Make sampling matrix------
  sm <- makesamplingmatrix(".", name = "Sampling_matrix_unaligned.csv")
  setwd(mainDir)

  ##Do MSA--------
  if (MSA == TRUE) {
    alignments<- align_piphy(mainDir = mainDir, clade=clade)
    write.csv(sm, "Sampling_matrix_aligned.csv")
    if (ALI == TRUE) {
      setwd(mainDir)
      sm<-ali_piphy(mainDir,alignments, sm=sm, clade=clade)
      write.csv(sm, "Sampling_matrix_alignedAli.csv")
      setwd(mainDir)
      unlink(list.dirs()[grep("Aliscore_", list.dirs(), ignore.case = T)], recursive=TRUE)

    }
  }
  setwd(mainDir)
  cat("\nDone!!")
}
