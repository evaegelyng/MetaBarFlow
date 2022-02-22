# Merging the resulting tables of the "sense" and the "antisense" analyses
# Authors: This script is based on the script C_Processing_with_DADA2.Rmd by Tobias G. Frøslev (see Frøslev et al. 2017). The custom function sumSequenceTables has been replaced with the standard mergeSequenceTables function, which now has a "sum" setting. 

args = commandArgs(trailingOnly=TRUE)

###Sum sense and antisense sequence tables
stAS <- file.path(args[1],"seqtab_AS_RDS")
stnsAS <- file.path(args[1],"seqtab.nochim_AS_RDS")
stSS <- file.path(args[1],"seqtab_SS_RDS")
stnsSS <- file.path(args[1],"seqtab.nochim_SS_RDS")
seqtab.nochim_AS <- readRDS(stnsAS)
seqtab.nochim_SS <- readRDS(stnsSS)
seqtab_AS <- readRDS(stAS)
seqtab_SS <- readRDS(stSS)
sumtable <- mergeSequenceTables(seqtab_SS,seqtab_AS, repeats="sum")
nochim_sumtable <- mergeSequenceTables(seqtab.nochim_SS,seqtab.nochim_AS, repeats="sum")
st <- file.path(args[1],"seqtab_RDS")
stns <- file.path(args[1],"seqtab.nochim_RDS")
saveRDS(sumtable,st)
saveRDS(nochim_sumtable,stns)
