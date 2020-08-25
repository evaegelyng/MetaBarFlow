args = commandArgs(trailingOnly=TRUE)

#Define a function for combining two or more tables, collapsing samples with similar names:  
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

##Merge Libraries
fileList = dir(path=args[1], pattern=NULL, all.files=FALSE,full.names=TRUE)
sumtable_all <- readRDS( paste(fileList[1], "/seqtab_RDS", sep='') )
nochim_sumtable_all <- readRDS( paste(fileList[1], "/seqtab.nochim_RDS", sep='') )
for (f in fileList[-1]){
  sumtable <- readRDS( paste(f, "/seqtab_RDS", sep='') )
  sumtable_nochim <- readRDS( paste(f, "/seqtab.nochim_RDS", sep='') )
  sumtable_all <- sumSequenceTables(sumtable_all,sumtable)
  nochim_sumtable_all <- sumSequenceTables(nochim_sumtable_all,sumtable_nochim)
}

stBoth <- file.path(args[2],"seqtab_Both")
stnsBoth <- file.path(args[2],"seqtab.nochim_Both")
saveRDS(sumtable_all,stBoth)
saveRDS(nochim_sumtable_all,stnsBoth)

#Transpose table, assign names, extract sequences and saving table, for further processing:
trans_nochim_sumtable <- as.data.frame(t(nochim_sumtable_all))
#Get DNA sequences
sequences <- row.names(trans_nochim_sumtable)
#Assign new rownames
row.names(trans_nochim_sumtable) <- paste0("seq",seq.int(nrow(trans_nochim_sumtable)))
tbname <- file.path(args[2],"DADA2_nochim.table")
{write.table(trans_nochim_sumtable,tbname,sep="\t",col.names = NA, quote=FALSE)}
#Extract OTUs (sequences)
sinkname <- file.path(args[2],"DADA2_nochim.otus")
{
  sink(sinkname)
  for (seqX in seq.int(nrow(trans_nochim_sumtable))) {
    header <- paste0(">","seq",seqX,"\n")
    cat(header)
    seqq <- paste0(sequences[seqX],"\n")
    cat(seqq)
  }
  sink()
}

#Define function to extract sequences sample-wise
extrSamDADA2 <- function(my_table) {
  out_path <- file.path(args[2],"DADA2_extracted_samples_nochim")
  if(!file_test("-d", out_path)) dir.create(out_path)
  for (sampleX in seq(1:dim(my_table)[1])){
    sinkname <- file.path(out_path, paste0(rownames(my_table)[sampleX],".fas"))
    {
      sink(sinkname)
      for (seqX in seq(1:dim(my_table)[2])) {
        if (my_table[sampleX,seqX] > 0) {
          header <- paste0(">",rownames(my_table)[sampleX],";size=",my_table[sampleX,seqX],";","\n")
          cat(header)
          seqq <- paste0(colnames(my_table)[seqX],"\n")
          cat(seqq)
        }
      }
      sink()
    }
  }
}

#Extract samplewise sequences from the non-chimera table using the above function:
extrSamDADA2(nochim_sumtable_all)



#Transpose table, assign names, extract sequences and saving table, for further processing:
trans_raw_sumtable <- as.data.frame(t(sumtable_all))
#Get DNA sequences
sequences <- row.names(trans_raw_sumtable)
#Assign new rownames
row.names(trans_raw_sumtable) <- paste0("seq",seq.int(nrow(trans_raw_sumtable)))
tbname <- file.path(args[2],"DADA2_raw.table")
{write.table(trans_raw_sumtable,tbname,sep="\t",col.names = NA, quote=FALSE)}
#Extract OTUs (sequences)
sinkname <- file.path(args[2],"DADA2_raw.otus")
{
  sink(sinkname)
  for (seqX in seq.int(nrow(trans_raw_sumtable))) {
    header <- paste0(">","seq",seqX,"\n")
    cat(header)
    seqq <- paste0(sequences[seqX],"\n")
    cat(seqq)
  }
  sink()
}

#Define function to extract sequences sample-wise
extrSamDADA2 <- function(my_table) {
  out_path <- file.path(args[2],"DADA2_extracted_samples_raw")
  if(!file_test("-d", out_path)) dir.create(out_path)
  for (sampleX in seq(1:dim(my_table)[1])){
    sinkname <- file.path(out_path, paste0(rownames(my_table)[sampleX],".fas"))
    {
      sink(sinkname)
      for (seqX in seq(1:dim(my_table)[2])) {
        if (my_table[sampleX,seqX] > 0) {
          header <- paste0(">",rownames(my_table)[sampleX],";size=",my_table[sampleX,seqX],";","\n")
          cat(header)
          seqq <- paste0(colnames(my_table)[seqX],"\n")
          cat(seqq)
        }
      }
      sink()
    }
  }
}

#Extract samplewise sequences from the non-chimera table using the above function:
extrSamDADA2(sumtable_all)
