args = commandArgs(trailingOnly=TRUE)
######################################################################################################################################################

# Run the taxonomy analysis of your metabarcoding data, after running it through DADA2 and BLASTN:

######################################################################################################################################################

print(args[1])
print(args[2])

# Read the completed blast results into a table
IDtable <- read.csv(file = args[1], sep='\t', header=F, as.is=TRUE)

# Read the possible problematic TaxIDs as a table
MergedTaxIDs<-read.table("~/eDNA/faststorage/blastdb/nt/taxdump/MergedTaxIDs", header=TRUE)

# Add header information
names(IDtable) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qlen","qcovs","sgi","sseq","ssciname","staxid")

# Extract only those rows where the qcovs score is 100
IDtable <- IDtable[IDtable$qcovs==100,]

# Optionally make a smaller test table first, to test if the script will run properly, e.g. just the first 100 rows:
# IDtable <- IDtable[1:100,]

# Load required packages
library(taxizedb) 
library(dplyr)
library(tidyr)

# Define the four taxonomy script functions (where FunctionX is a wrapper that runs 1, 2, 3 in one go):

# FunctionX
# Wrapper function using the three main functions - each step can also be done manually
assign_taxonomy <- function(table,upper_margin=0.5,lower_margin=2, remove = c("")) {
  pf <- prefilter(table, upper_margin, lower_margin, remove)
  
  # Replace old tax ids with new ones if they have a match in mergedtaxids dataframe #CKF 
  pf$OldTaxID<-pf$staxid # Making an extra column on pf, same as staxid but with a new name so it can match with the mergedtaxid dataframe
  pf_intermediate<-pf %>%
  dplyr::left_join(MergedTaxIDs, by=c("OldTaxID"))  # now we can join with MergedTaxIDs using the oldtaxid column
  
  # Print which IDs are being changed 
  pf_changed<-pf_intermediate %>%
  dplyr::filter(!is.na(NewTaxID))
  
  if (nrow(pf_changed)>0) {
  
  for (i in 1:nrow(pf_changed)) {
  pf_sub<-pf_changed[i,]
  print(paste("Warning:", pf_sub$qseqid, "has an outdated taxid. Overwriting outdated taxid:", pf_sub$OldTaxID, "to new taxid:", pf_sub$NewTaxID, sep=" "))
  }
  
 }  
  
  pf_new<-pf_intermediate %>%
  dplyr::mutate(staxid=ifelse(is.na(NewTaxID), staxid, NewTaxID)) %>% # If there is no match, then staxid stays the same, otherwise it takes the id in NewTaxID
  dplyr::select(-c(OldTaxID, NewTaxID)) # Now remove the two columns added to the dataframe, no longer needed 
  
  
  ##
  
  gc <- get_classification(pf_new)
  cf <- evaluate_classification(gc[[1]])     # AGR # gc[2] is the wrong taxid not classified 
  result <- list(classified_table=cf$taxonon_table, all_classifications=cf$all_taxa_table, all_classifications_summed=cf$all_taxa_table_summed,upper=upper_margin, lower=lower_margin, removed=remove)
  if (length(gc[[2]]) != 0) {
    print(paste0("Taxids not found in the classification: ", gc[2])) # AGR - Print the list of not matched taxids  
  }
  return(result) 
}

# Function1
# Filter data OTU-wise according to upper and lower margin set, and taxa to exclude
prefilter <- function(IDtable, upper_margin=0.5, lower_margin=2, remove = c("uncultured", "environmental")) {
  new_IDtable <- IDtable[0,] # prepare filtered matchlist
  IDtable <- IDtable[!IDtable$staxid == "N/A",]
  ids <- names(table(IDtable$qseqid))
  i=1
  o=length(ids)
  for (name in ids) {
    test <- IDtable[which(IDtable$qseqid == name),] # select all lines for a query
    if (nchar(remove[1])>0) {
      test2 <- test
      for (rm in 1:length(remove)) {
        test2 <- test2[!grepl(remove[rm], test2$ssciname,ignore.case = TRUE),]
      }
      if (nrow(test2) > 1) {test <- test2}
    }
    max <- max(test$pident)
    upper <- max-upper_margin
    lower <- max-lower_margin
    test <- test[which(test$pident >= lower),] # select all lines for a query
    test$margin <- "lower"
    test[test$pident >= upper,"margin"] <- "upper"
    
    new_IDtable = rbind(new_IDtable,test) # add this row to the filtered IDtable
    i=i+1
  }
  return(new_IDtable)
}

# Function2
# Get full taxonomic path for all hits within the upper limit of each OTU. Identical species are only queried once.

get_classification <- function(IDtable2) {
  require(taxizedb)
  all_staxids <- names(table(IDtable2$staxid[IDtable2$margin=="upper"])) # get all taxids for table
  all_classifications <- list() # prepare list for taxize output
  o=length(all_staxids) # number of taxids
  
  Start_from <- 1 # change if loop needs to be restarted due to time-out
  
  wrong_taxid_matches <- c()
  remove_entries <- c()
  
  #Get ncbi classification of each entry
  for (cl in Start_from:o) { # the taxize command "classification" can be run on the all_staxids vector in one line, but often there is
    #a timeout command, therefor this loop workaround.
    print(paste0("step 1 of 3: processing: ", cl , " of ", o , " taxids")) # make a progressline (indicating the index the loops needs to be
    #restarted from if it quits)
    tax_match <- classification(all_staxids[cl], db = "ncbi")   # AGR 
    if (is.na(tax_match) == TRUE) {
      wrong_taxid_matches <- c(wrong_taxid_matches,all_staxids[cl])
      remove_entries <- c(remove_entries,cl)
    } else {
      all_classifications[cl] <- tax_match
    }                                                           # AGR 
  }
  
  # In case there are still bad taxIDs, we delete all sequences with the bad tax ID (Adrian + Mads solve 15-01-2020)
  if (length(wrong_taxid_matches) != 0) {
    all_classifications <- all_classifications[-remove_entries]
  }
  
  #Construct a taxonomic path from each classification
  output <- data.frame(staxid=character(),kingdom=character(), phylum=character(),class=character(),order=character(),family=character(),genus=character(),species=character(), stringsAsFactors=FALSE)
  totalnames <- length(all_staxids) - length(wrong_taxid_matches)
  
  ## - 1 ## this is if you have one NA, would run on test file, but not necessarily on other files with more NAs
  for (curpart in seq(1:totalnames)) {
    print(paste0("step 2 of 3: progress: ", round(((curpart/totalnames) * 100),0) ,"%")) # make a progress line
    currenttaxon <- all_classifications[curpart][[1]]
    if (nchar(currenttaxon[1]) > 0) {
      spec <- all_staxids[curpart]
      output[curpart,"kingdom"] <- currenttaxon[which(currenttaxon$rank == "kingdom"),"name"][1]
      output[curpart,"phylum"] <- currenttaxon[which(currenttaxon$rank == "phylum"),"name"][1]
      output[curpart,"class"] <- currenttaxon[which(currenttaxon$rank == "class"),"name"][1]
      output[curpart,"order"] <- currenttaxon[which(currenttaxon$rank == "order"),"name"][1]
      output[curpart,"family"] <- currenttaxon[which(currenttaxon$rank == "family"),"name"][1]
      output[curpart,"genus"] <- currenttaxon[which(currenttaxon$rank == "genus"),"name"][1]
      output[curpart,"species"] <- currenttaxon[which(currenttaxon$rank == "species"),"name"][1]
      output[curpart,"staxid"] <-  spec # add that row to the filtered IDtable
    }
  }
  taxonomic_info <- merge(IDtable2,output,by = "staxid", all=TRUE)
  taxonomic_info$species[is.na(taxonomic_info$species)] <- taxonomic_info$ssciname[is.na(taxonomic_info$species)]
  return(list(taxonomic_info,wrong_taxid_matches))
}

# Function3
# Function for evaluating the taxonomic assignment of each OTU. All hits within the upper margin are used in the evaluation weighted by their evalue, so that suboptimal matches have a lower weight. All hits within the lower margin are put into the output (but not used for evaluating classification)
evaluate_classification <- function(classified) {
  require(tidyr)
  require(dplyr)
  ids <- names(table(classified$qseqid))
  i <- 1
  for (name in ids) {
    print(paste0("last step: progress: ", round(((i/length(ids)) * 100),0) ,"%")) # make a progressline
    test <- classified[which(classified$qseqid == name),]
    test2 <- test %>% filter(margin == "upper")
    test2$score <- 100*(1/test2$evalue)/sum(1/test2$evalue)  # HER BEREGSES SCOREN FOR ALLE MATCHES PER OTU
    test4 <- test2 %>% filter(margin == "upper") %>%
      dplyr::select(margin,qseqid,sgi,sseq,staxid,pident,score,qcovs,kingdom,phylum,class,order,family,genus,species) %>% 
      group_by(qseqid,kingdom, phylum,class,order,family,genus,species) %>% 
      mutate(species_score=sum(score)) %>% 
      group_by(qseqid,kingdom, phylum,class,order,family,genus) %>% 
      mutate(genus_score=sum(score)) %>%
      group_by(qseqid,kingdom, phylum,class,order,family) %>% 
      mutate(family_score=sum(score))%>%
      group_by(qseqid,kingdom, phylum,class,order) %>% 
      mutate(order_score=sum(score)) %>%
      group_by(qseqid,kingdom, phylum,class) %>% 
      mutate(class_score=sum(score)) %>%
      group_by(qseqid,kingdom, phylum) %>% 
      mutate(phylum_score=sum(score)) %>%
      group_by(qseqid,kingdom) %>% 
      mutate(kingdom_score=sum(score)) %>% ungroup() %>%
      arrange(-kingdom_score,-phylum_score,-class_score,-order_score,-family_score,-genus_score,-species_score)
    test3 <- test4 %>% slice(1)
    test5 <- test4 %>% distinct(qseqid,sgi,sseq,pident,qcovs,kingdom,phylum,class,order,family,genus,species,sseq,kingdom_score,phylum_score,class_score,order_score,family_score,genus_score,species_score) 
    string1 <- test %>% dplyr::select(species,pident) %>% 
      distinct(species,pident) %>% arrange(-pident) %>% t()
    string2 <- toString(unlist(string1))
    test3$alternatives <- string2
    if (i == 1){result <- test3} else {
      result <- rbind(result,test3)
    }
    if (i == 1){result2 <- test2} else {
      result2 <- rbind(result2,test2)
    }
    if (i == 1){result3 <- test5} else {
      result3 <- rbind(result3,test5)
    }
    i=i+1
  }
  total_result <- list(taxonon_table = result, all_taxa_table=result2, all_taxa_table_summed=result3)
  return(total_result)
}


# Classify your OTUs by running the wrapper (functionX) "assign_taxonomy" like this. Consider whether it makes sense to remove specific hits or taxa from the evaluation (see explanation below):
my_classified_result <- assign_taxonomy(IDtable, upper_margin = 0.5, lower_margin = 2, remove = c("uncultured", "environmental"))

#Write the result to a table
write.table(my_classified_result$classified_table, args[2], sep = "\t", quote = F, row.names = F)

# Explanation to input
#   INPUT.blasthits is the blast-results
#   upper_margin is the margin used for suboptimal hits used for classification - e.g. a margin of 0.5 means that hits of 100% to 99.5% is used og 95% to 94.5%, if the best hit is 100% or 95% respectively.
#   lower_margin: hits down to this margin from the best hit are shown in the output as alternative possibilities, but not used for taxonomic evaluation.
#   remove: a vector of taxa to exclude from the evaluation. Could be e.g. remove = c("uncultured","environmental") to exclude hits with no precise annotation, or names of species known not to be in the study area.
#
# Explanation to output
#   the output is a list with
#   $classified_table : the table with all OTUs classified. One row per OTU
#        this table contains the estimated best classification at all taxonomic levels, based on a weighted score (of the evalue) of the hits in the upper_margin, 
#        each taxonomic level gets a score indicating the agreement on the selected classification at that level..
#        also a string of alternatives and their matches (%) this string includes hits from both upper and lower margin
#   $all_classifications: this is the table used to make the classified_table. It contains all hits above lower_magrin for all OTUs and their classifications (only upper_margin).
#   ...and the input parameters
