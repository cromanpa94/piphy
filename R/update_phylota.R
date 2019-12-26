update_phylota <-
  function(ingroup,
           outgroup,
           genes = NULL,
           MSA = FALSE,
           ALI = FALSE,
           correct_db = TRUE,
           delete_all = TRUE) {
    ##New version
    options(warn = -1)

    clade <- ingroup
    fn <- "Unaligned"
    fn2 <- "Aligned"
    fn2.1 <- "Aliscore"
    fn3 <- "phylota.sampling.csv"
    fn4 <- "included.species.csv"
    if (file.exists(fn))
      unlink(fn, recursive = T)
    if (file.exists(fn2))
      unlink(fn2, recursive = T)
    if (file.exists(fn2.1))
      unlink(fn2.1, recursive = T)
    if (file.exists(fn3))
      unlink(fn3, recursive = T)
    if (file.exists(fn4))
      unlink(fn4, recursive = T)
    if (length(grep("Molecular_*", list.files("."), value = T)) == 0) {
    } else {
      unlink(grep("Molecular_*", list.files("."), value = T))
    }

    cat("\n Get PhyLota clusters \n")

    piphy::get_phylota(ingroup, MSA = MSA, ALI = ALI)
    mainDirect <-  getwd()
    mainDir<-mainDirect

    subwd <- paste0(mainDirect, "/Unaligned/")

    ##Check for new species in genbank

    cat("Checking for novel species in genbank \n")

    spp_genbank <-
      taxize::downstream(ingroup, db = "ncbi", downto = 'species')[[1]]

    spp_sampled <-read.csv(list.files(
      path = subwd,
      pattern = c("Sampling_matrix"),
      full.names = T
    ))
    new_spp <-
      setdiff(gsub(" ", "_", spp_genbank$childtaxa_name), spp_sampled$Species)

    if (length(new_spp) == 0) {
      if (file.exists(fn))
        unlink(fn, recursive = T)

      stop("Nothing to be added")
    } else{
      cat("Tranforming each cluster into a Blast DB \n")

      ##First make DB for each cluster
      files <- list.files(path = subwd, pattern = c("Cluster"))
      files <- Filter(function(x)
        grepl("\\.fasta$", x), files)

      subwd<- if(Sys.info()[['sysname']]=="Windows") {
        gsub("\\", "/", subwd, fixed=T) } else {subwd}

      for (i in 1:length(files)) {
        args <-
          paste0("-in ", subwd, files[i], " -input_type fasta -dbtype nucl")
        system2(
          command = 'makeblastdb',
          args = args,
          wait = FALSE,
          stdout = TRUE
        )
        Sys.sleep(1)
        print(i)
      }



      cat(
        length(new_spp),
        " species can be added to the ",
        dim(spp_sampled)[1],
        "that are already sampled in PhyloTa \n"
      )

        sampled_gene_names <- genes


      ##Then run each sequence against each cluster

      ##Start!
      ##
      new_spp2 <-   if (is.null(outgroup) == T) {
        new_spp
      } else{
        new_spp <- c(new_spp, outgroup)
      }

      data <-
        scrape.genbank(gsub("_", " ", new_spp), sampled_gene_names)



      ###Now, I download each sequence and test where it fits better

      cat("\nFitting the new sequences into the existing clusters \n")

      acce_new <- as.character(na.omit(unlist(data[, -1])))

      if (length(acce_new) == 0) {
        stop("No new sequences found")
      }

      n_clust <- length(list.files(path = subwd, pattern = c("csv")))
      species_included <- list()
      for (i in 1:length(acce_new)) {
        sequence <- read.GenBank(acce_new[i])
        Acc_spp <- names(sequence)
        names(sequence) <- attr(sequence, "species")

        for (j in 1:n_clust) {
          write.dna(sequence, "sequence.fasta", format = "fasta")
          que <-
            paste0("-outfmt 6 -query sequence.fasta -out test.txt -db ",
                   subwd,
                   files[j])
          system2(
            command = 'blastn',
            args = que,
            wait = FALSE,
            stdout = TRUE
          )

          if (file.info("test.txt")$size > 0) {
            cluster <- ape::read.dna(paste0(subwd, files[j]), format = "fasta")
            cluster[length(cluster) + 1] <- sequence
            names(cluster)[length(cluster)] <- names(sequence)
            write.dna(cluster, paste0(subwd, files[j]), format = "fasta")
            species_included[[i]] <-
              data.frame(
                Included_species = names(sequence),
                AN = Acc_spp,
                Cluster = files[j]
              )
            cat("*******Sequence matched!******* \n")
            break ##If sequence is matched, we shold stop.
          } else {
            cat("****Still working \n")
          }

        }

        cat("**Done with sequence",
            i,
            "of",
            length(acce_new),
            "****\n")

      }
      df <- do.call("rbind", species_included)

      if (file.exists("sequence.fasta"))
        unlink("sequence.fasta")
      if (file.exists("test.txt"))
        unlink("test.txt")


    }

    setwd(mainDirect)
    cat("\n Retrieving accession numbers \n")

    ##Need to include the new seqs to the cluster matrix



    ##Do MSA--------
    ##Read fasta files
    if(MSA==TRUE){

      cat("\nalignment in process. Please be patient\n")

      setwd(file.path(mainDir, "Unaligned"))

      temp = list.files(pattern="*.fasta")

      ##Chose method: "Muscle", "ClustalOmega", "ClustalW"

      align_PHYLOTA<-function(file.names){
        mySequences <- readDNAStringSet(file.names)
        myFirstAlignment <- msa(mySequences, "ClustalOmega")
        aln <- as.DNAbin(msaConvert(myFirstAlignment, type="phangorn::phyDat"))
      }

      alignments<-pblapply(temp,align_PHYLOTA)

      options(warn=-1)
      dir.create(file.path(mainDir, "Aligned"))
      options(warn=-0)

      setwd(file.path(mainDir, "Aligned"))

      for (i in 1: length(alignments)){
        write.dna(alignments[[i]], paste0(clade, "_", "Alignment","Cluster" ,i, ".fasta"), format="fasta")
      }
      cat("\nAligned and unaligned sequences are in separated folders within your working directory")

      if (ALI==TRUE){
        ####For Aliscore------

        setwd(mainDir)

        if( "Aliscore_v.2.0" %in% list.files() == FALSE){
          download.file("http://software.zfmk.de/ALISCORE_v2.0.zip",'Aliscore_v.2.0.zip')
          unzip("Aliscore_v.2.0.zip")
          unzip("ALISCORE_v2.0/Aliscore_v.2.0.zip")
        } else
          options(warn=-1)
        dir.create(file.path(mainDir, "Aliscore"))
        options(warn=0)

        setwd(file.path(mainDir, "Aliscore"))

        plot.progress <- function(percent) {
          plot(c(0,100), c(0,1), type='n', xlab='', ylab='', yaxt='n')
          rect(0, 0.1, percent*100, 0.9, col='blue')
          title(paste('Progress: ', round(percent*100,2), '%', sep=''))
        }


        aln_aliscored<-list()
        for (i in 1: length(alignments)){
          aln_aliscored[[i]]<-aliscore(alignments[[i]], gaps = "ambiguous", w = 3, path = paste0(mainDir, "/Aliscore_v.2.0"))
          write.dna(aln_aliscored[[i]], paste0(clade, "_", "AliScore", "_Alignment","_Cluster_" ,i, ".fasta"), format="fasta")
          plot.progress(i/length(alignments))
          cat("\nCurated alignments are under ALISCORE folder")

        }
        setwd(mainDirect)
      } else {}
      setwd(mainDirect)
    } else {}

    setwd(mainDir)
    setwd(mainDirect)
    ####Remove subspecies for the last time

    dfs <-
      lapply(list.files(
        pattern = c(".csv"),
        full.names = T
      ), read.csv)

    files_names<-list.files(
      pattern = c(".csv"),
      full.names = T
    )

    dfs_co<-lapply(1:length(dfs), function(x){
     vec<- gsub(" ", "_",as.character(dfs[[x]][,grep("species",colnames(dfs[[x]]),ignore.case=TRUE)]))
      num1<-which(sapply(strsplit(vec, "_"), length)==3)
      df_corr<- if(length(num1)==0){ dfs[[x]]}else{dfs[[x]][-num1,]}
      return(df_corr)
    })


    lapply(1:length(dfs_co), function(x){
      write.csv(dfs_co[[x]], file=files_names[x])
    })


    alns <-
      lapply(list.files(
        pattern = c(".fasta"),
        full.names = T, recursive = T
      ), ape::read.dna, "fasta")

    names_alns <-list.files(
      pattern = c(".fasta"),
      full.names = T, recursive = T
    )


    lapply(1:length(alns), function(x){
      aln_ex<-alns[[x]]
      num1<-which(sapply(strsplit(labels(aln_ex), "_"), length)==3)
      if(length(num1)==0){ }else{
        newaln<- if(is.null(dim(aln_ex)) ==F){
          aln_ex[-(num1),]
        }else{
          aln_ex[-(num1)]
        }
        write.dna(newaln, names_alns[x], format = "fasta")
      }
    })

    ####Re organize datasets

    dfs <-
      lapply(list.files(
        pattern = c(".csv"),
        full.names = T
      ), read.csv)

    files_names<-list.files(
      pattern = c(".csv"),
      full.names = T
    )

    species_correct<- unique(gsub(" ", "_",as.character(dfs[[2]][,grep("species",colnames(dfs[[x]]),ignore.case=TRUE)])))

    dfs_co<-lapply(1:length(dfs), function(x){
      spps<- gsub(" ", "_",as.character(dfs[[x]][,grep("species",colnames(dfs[[x]]),ignore.case=TRUE)]))
      new_df<-dfs[[x]][!duplicated(spps),]
      spps<- gsub(" ", "_",as.character(new_df[,grep("species",colnames(new_df),ignore.case=TRUE)]))
      new_df<-new_df[which(spps %in% species_correct),]
      return(new_df)
    })


    lapply(1:length(dfs_co), function(x){
      write.csv(dfs_co[[x]], file=files_names[x])
    })


    alns <-
      lapply(list.files(
        pattern = c(".fasta"),
        full.names = T, recursive = T
      ), ape::read.dna, "fasta")

    names_alns <-list.files(
      pattern = c(".fasta"),
      full.names = T, recursive = T
    )


    lapply(1:length(alns), function(x){
      aln_ex<-alns[[x]]
      num1<-which(sapply(strsplit(labels(aln_ex), "_"), length)==3)
      num2<-which(! labels(aln_ex) %in% species_correct)

      num3<-c(num1,num2)

      if(length(num3)==0){ }else{
        newaln<- if(is.null(dim(aln_ex)) ==F){
          aln_ex[-(num3),]
        }else{
          aln_ex[-(num3)]
        }
        write.dna(newaln, names_alns[x], format = "fasta")
      }
    })


    cat("\nDone!!")

    options(warn = 0)
  }
