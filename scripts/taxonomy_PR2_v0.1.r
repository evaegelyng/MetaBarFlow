#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Libraries
library(dada2)

# Function to read a multiline FASTA file properly into R (*Generated using AI*. Avoids including more R-packages in the environment that could otherwise handle multiline fastafiles)
read_fasta <- function(file_path) {
  # Step 1: Read the entire file
  lines <- readLines(file_path)
  
  # Initialize variables
  headers <- c()  # to store the sequence headers
  sequences <- c()  # to store the sequences
  
  # Temporary variable to hold sequence data
  current_sequence <- ""
  
  # Step 2: Loop through each line of the file
  for (line in lines) {
    if (substr(line, 1, 1) == ">") {
      # If the line starts with ">", it's a header
      if (current_sequence != "") {
        # If there's a current sequence, save it
        sequences <- c(sequences, current_sequence)
        current_sequence <- ""  # reset for the next sequence
      }
      headers <- c(headers, line)  # store the header
    } else {
      # If the line is a sequence, append it to the current sequence
      current_sequence <- paste0(current_sequence, line)
    }
  }
  
  # Add the last sequence (after the loop ends)
  if (current_sequence != "") {
    sequences <- c(sequences, current_sequence)
  }
  
  # Step 3: Return a list with headers and sequences
  return(list(headers = headers, sequences = sequences))
}

# dada2 tutorial for assigning taxonomy
# https://benjjneb.github.io/dada2/assign.html

# Load sequences
sequences_fasta <- read_fasta(file = args[1])

# Split sequences and sequence names, remove ">" from sequence name
seqs <- sequences_fasta$sequences
seq_names <- sequences_fasta$headers
seq_names <- gsub(">", "", seq_names)

# Set seed for reproducibility
set.seed(100)

# Assign taxonomy using the RDP classifier incorporated in DADA2's "assignTaxonomy"
taxa <- assignTaxonomy(seqs, 
                       refFasta = "/YOUR_PATH/pr2_version_5.0.0_SSU_dada2.fasta.gz", # wget https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_dada2.fasta.gz
                       minBoot = 80,
                       multithread = 1, 
                       tryRC = T, 
                       outputBootstraps = T)

# Collate the sequence names, bootstrap values and taxonomic annotations into a table and write output
df <- data.frame(seq_names, taxa$boot, taxa$tax)
write.table(x = df, file = args[2], quote = F, sep = "\t", row.names = F, col.names = c("Seq","Boot_Domain","Boot_Supergroup","Boot_Division","Boot_Subdivision","Boot_Class","Boot_Order","Boot_Family","Boot_Genus","Boot_Species","Domain","Supergroup","Division","Subdivision","Class","Order","Family","Genus","Species"))
