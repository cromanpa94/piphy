
#' update_phylota(ingroup = "centrolene", outgroup = "cochranella",
#'               genes = c("12S","RAG1", "ND1","16S"), MSA=T, ALI=T)

update_phylota <-
  function(ingroup,
           outgroup,
           genes = NULL,
           MSA = FALSE,
           ALI = FALSE) {
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

    piphy::get_phylota(ingroup, MSA = F, ALI = F)
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

      acce_new <- unique(as.character(na.omit(unlist(data[, -1]))))

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
            if( names(sequence) %in% names(cluster)  ){}else{
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
            }
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

    cat("\n Retrieving accession numbers \n")

    ##Need to include the new seqs to the cluster matrix


    addedseqs<-reshape(df, idvar = "Included_species", timevar = "Cluster", direction = "wide")
    colnames(addedseqs)<-c("Species", paste0("Cluster_",as.numeric(gsub("[^\\d]+", "", colnames(addedseqs)[-1], perl=TRUE))))


    final_sampling<-as.data.frame(rbindlist(list(spp_sampled, addedseqs), fill = TRUE))[,-1]


    write.csv( final_sampling,"Unaligned/Sampling_matrix_unaligned.csv" )

    setwd(mainDirect)

    file.remove(
      list.files(
        ".",
        pattern = ".fasta.",
        include.dirs = F,
        full.names = T,
        recursive = T
      )
    )

    if (MSA == TRUE) {
     alignments<- align_piphy(mainDir = mainDir, clade=ingroup)
      write.csv(final_sampling, "Sampling_matrix_aligned.csv")
      setwd(mainDir)
      if (ALI == TRUE) {
        setwd(mainDir)
        sm<-ali_piphy(mainDir, alignments, sm=final_sampling, clade=ingroup)
        write.csv(sm, "Sampling_matrix_aliscore.csv")

        setwd(mainDir)
        unlink(list.dirs()[grep("Aliscore_", list.dirs(), ignore.case = T)], recursive =
                 TRUE)
      }
    }

    cat("\nDone!!")

    options(warn = 0)
  }
