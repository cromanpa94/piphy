# Genbank scraping functions
# Copyright (C) 2015 Anusha Beer <anbeer29@gmail.com>
#
# Permission to use, copy, modify, and/or distribute this software for
# any purpose with or without fee is hereby granted, provided that the
# above copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
library("reutils")
# Search Genbank for search term e.g. species name. Returns a dataframe
# of all available (nucleotide) sequences on Genbank.
#
# Give it a string e.g. species name. Returns a dataframe. Each row is a
# result and each column is an attribute of the results. Some of the imprtant
# attributes are "Title", "OSLT". Title contains a description of the
# sequence. OSLT contains acession number.
grab.results <- function (term) {
  # Search for the given term on nuccore. This gives us a list of
  # record IDs.
  ids <- esearch(term, db="nuccore")
  # Grab summaries for the given record IDs, as a sort-of data frame.
  sum <- esummary(ids, db="nuccore")
  data <- content(sum, as="parsed")
  # For some reason, this parser gives us lists of lists instead of a
  # proper data frame (which should be lists of vectors). Return a
  # fixed-up version.
  data.frame(lapply(data, as.character), stringsAsFactors=FALSE)
}
#######################################################################
# Takes a dataframe of Genbank results and a keyword string and returns a
# string containing an accession number. This is the accession number of
# the first result in the table that had the keyword in its Title field.
first.of.type <- function (results, type) {
  # Filter rows whose title contains the given text
  filtered <- results[grepl(type, results$Title, fixed=TRUE),]
  # Return the first accession number
  filtered$OSLT[1]
}
#######################################################################
# Given a list of species and genes, return a dataframe of accession numbers.
#
# Takes string vectors as arguments and returns a dataframe. The rows of the
# dataframe are species. There is a column for each gene. These columns are
# filled with accession numbers. The species names are put in a column called
# "Species".
scrape.genbank <- function (species, genes) {
  # Create dataset skeleton
  data <- data.frame(Species=species)
  for (i in genes) {
    data[,i] <- rep(NA, length(species))
  }
  # Look up data for each species
  for (i in 1:length(species)) {
    n <- species[i]
    print(sprintf("Looking up %s (%d/%d)...", n, i, length(species)))
    tryCatch({
      r <- grab.results(paste(n, "[Primary Organism]"))
      for (g in genes) {
        data[i,g] <- first.of.type(r, g)
      }
    }, error = function(e) { })
  }
  data
}

