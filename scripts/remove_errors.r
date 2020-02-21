args = commandArgs(trailingOnly=TRUE)

library(dada2)

main_path <- args[1]

start.time <- Sys.time()

#Processing the set of files containing the forward primer in the R1 reads (the sense reads):
filt_path <- file.path(main_path, "DADA2_SS/filtered/matched") 
fns <- list.files(filt_path)
fastqs <- fns[grepl(".fastq.gz", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("_F_", fastqs)]
fnRs <- fastqs[grepl("_R_", fastqs)]
filtFs <- file.path(filt_path, fnFs)
filtRs <- file.path(filt_path, fnRs)
filtFs <- filtFs[sapply(filtFs, file.size) > 1]
filtRs <- filtRs[sapply(filtRs, file.size) > 1]
SSsample.names <- sapply(strsplit(basename(filtFs), "_F_"), `[`, 1)
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- SSsample.names
names(derepRs) <- SSsample.names
SSdadaFs <- dada(derepFs, err=errF, multithread = TRUE)
SSdadaRs <- dada(derepRs, err=errR, multithread = TRUE)
SSmergers <- mergePairs(SSdadaFs, derepFs, SSdadaRs, derepRs, verbose=TRUE,minOverlap = 5)
seqtab_SS <- makeSequenceTable(SSmergers[names(SSmergers)])
seqtab.nochim_SS <- removeBimeraDenovo(seqtab_SS, verbose=TRUE)
stSS <- file.path(args[1],"seqtab_SS_RDS")
stnsSS <- file.path(args[1],"seqtab.nochim_SS_RDS")
saveRDS(seqtab_SS,stSS)
saveRDS(seqtab.nochim_SS,stnsSS)

#Then DADA2 processing of "the antisense" reads:
filt_path <- file.path(main_path, "DADA2_AS/filtered/matched") 
fns <- list.files(filt_path)
fastqs <- fns[grepl(".fastq.gz", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("_F_", fastqs)]
fnRs <- fastqs[grepl("_R_", fastqs)]
filtFs <- file.path(filt_path, fnFs)
filtRs <- file.path(filt_path, fnRs)
filtFs <- filtFs[sapply(filtFs, file.size) > 1]
filtRs <- filtRs[sapply(filtRs, file.size) > 1]
ASsample.names <- sapply(strsplit(basename(filtFs), "_F_"), `[`, 1)
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- ASsample.names
names(derepRs) <- ASsample.names
ASdadaFs <- dada(derepFs, err=errF, multithread = TRUE)
ASdadaRs <- dada(derepRs, err=errR, multithread = TRUE)
ASmergers <- mergePairs(ASdadaFs, derepFs, ASdadaRs, derepRs, verbose=TRUE,minOverlap = 5)
seqtab_AS <- makeSequenceTable(ASmergers[names(ASmergers)])
seqtab.nochim_AS <- removeBimeraDenovo(seqtab_AS, verbose=TRUE)
stAS <- file.path(args[1],"seqtab_AS_RDS")
stnsAS <- file.path(args[1],"seqtab.nochim_AS_RDS")
saveRDS(seqtab_AS,stAS)
saveRDS(seqtab.nochim_AS,stnsAS)
