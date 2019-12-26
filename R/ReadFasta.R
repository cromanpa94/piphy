ReadFasta<-function(file) {
  options(warn=-1)
  fasta<-readLines(file)
  options(warn=0)
  fasta[1]<-gsub("<html><pre>", "", fasta[1])
  fasan<-fasta
  fasan_2<-fasan[-length(fasta)[1]]
  # Identify header lines
  ind<-grep(">", fasan_2)
  # Identify the sequence lines
  s<-data.frame(ind=ind, from=ind+1, to=c((ind-1)[-1], length(fasan_2)))
  # Process sequence lines
  seqs<-rep(NA, length(ind))
  for(i in 1:length(ind)) {
    seqs[i]<-paste(fasan_2[s$from[i]:s$to[i]], collapse="")
  }
  # Create a data frame
  DF<-data.frame(name=gsub(">", "", fasan_2[ind]), sequence=seqs)
  # Return the data frame as a result object from the function
  return(DF)
}
