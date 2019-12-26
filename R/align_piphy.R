align_piphy<-function(mainDir){
cat("\nalignment in process. Please be patient\n")

setwd(file.path(mainDir, "unaligned"))

temp = list.files(pattern = "*Cluster")

##Chose method: "Muscle", "ClustalOmega", "ClustalW"

align_PHYLOTA <- function(file.names) {
  mySequences <- readDNAStringSet(file.names)
  myFirstAlignment <- msa(mySequences, "ClustalOmega")
  aln <-
    as.DNAbin(msaConvert(myFirstAlignment, type = "phangorn::phyDat"))
}

alignments <- pblapply(temp, align_PHYLOTA)

options(warn = -1)
dir.create(file.path(mainDir, "Aligned"))
options(warn = -0)

setwd(file.path(mainDir, "Aligned"))

for (i in 1:length(alignments)) {
  write.dna(alignments[[i]],
            paste0(clade, "_", "Alignment", "Cluster" , i, ".fasta"),
            format = "fasta")
}

}
