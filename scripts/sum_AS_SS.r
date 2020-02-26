args = commandArgs(trailingOnly=TRUE)

#Define a function for combining two or more sequence tables:  
sumSequenceTables <- function(table1, table2, ..., orderBy = "abundance") {
  # Combine passed tables into a list
  tables <- list(table1, table2)
  tables <- c(tables, list(...))
  # Validate tables
  if(!(all(sapply(tables, dada2:::is.sequence.table)))) {
    stop("At least two valid sequence tables, and no invalid objects, are expected.")
  }
  sample.names <- rownames(tables[[1]])
  for(i in seq(2, length(tables))) {
    sample.names <- c(sample.names, rownames(tables[[i]]))
  }
  seqs <- unique(c(sapply(tables, colnames), recursive=TRUE))
  sams <- unique(sample.names)
  # Make merged table
  rval <- matrix(0L, nrow=length(sams), ncol=length(seqs))
  rownames(rval) <- sams
  colnames(rval) <- seqs
  for(tab in tables) {
    rval[rownames(tab), colnames(tab)] <- rval[rownames(tab), colnames(tab)] + tab
  }
  # Order columns
  if(!is.null(orderBy)) {
    if(orderBy == "abundance") {
      rval <- rval[,order(colSums(rval), decreasing=TRUE),drop=FALSE]
    } else if(orderBy == "nsamples") {
      rval <- rval[,order(colSums(rval>0), decreasing=TRUE),drop=FALSE]
    }
  }
  rval
}

###Sum sense and antisense sequence tables
stAS <- file.path(args[1],"seqtab_AS_RDS")
stnsAS <- file.path(args[1],"seqtab.nochim_AS_RDS")
stSS <- file.path(args[1],"seqtab_SS_RDS")
stnsSS <- file.path(args[1],"seqtab.nochim_SS_RDS")
seqtab.nochim_AS <- readRDS(stnsAS)
seqtab.nochim_SS <- readRDS(stnsSS)
seqtab_AS <- readRDS(stAS)
seqtab_SS <- readRDS(stSS)
sumtable <- sumSequenceTables(seqtab_SS,seqtab_AS) 
nochim_sumtable <- sumSequenceTables(seqtab.nochim_SS,seqtab.nochim_AS)
st <- file.path(args[1],"seqtab_RDS")
stns <- file.path(args[1],"seqtab.nochim_RDS")
saveRDS(sumtable,st)
saveRDS(nochim_sumtable,stns)
