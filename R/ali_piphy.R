ali_piphy<-function(mainDir){

if ("Aliscore_v.2.0" %in% list.files() == FALSE) {
  download.file("http://software.zfmk.de/ALISCORE_v2.0.zip",
                'Aliscore_v.2.0.zip')
  unzip("Aliscore_v.2.0.zip")
  unzip("ALISCORE_v2.0/Aliscore_v.2.0.zip")
} else{
}
options(warn = -1)
dir.create(file.path(mainDir, "Aliscore"))
options(warn = 0)

setwd(file.path(mainDir, "Aliscore"))

aln_aliscored <- list()
for (i in 1:length(alignments)) {
  aln_aliscored[[i]] <-
    aliscore(
      alignments[[i]],
      gaps = "ambiguous",
      w = 3,
      exec  = paste0(mainDir, "/Aliscore_v.2.0/Aliscore.02.2.pl")
    )
  write.dna(
    aln_aliscored[[i]],
    paste0(
      clade,
      "_",
      "AliScore",
      "_Alignment",
      "_Cluster_" ,
      i,
      ".fasta"
    ),
    format = "fasta"
  )
  plot.progress(i / length(alignments))
  cat("\nCurated alignments are under ALISCORE folder")

}

'%nin%' = Negate('%in%')

for (i in 2:ncol(sm)) {
  target <- sm[!is.na(sm[, i]) , 1]
  in_original <- labels(aln_aliscored[[i - 1]]) %nin% target
  if (any(in_original)) {
    sm[which(sm[, 1] == labels(aln_aliscored[[i - 1]])[in_original]), i] <-
      NA
  }
}

}
