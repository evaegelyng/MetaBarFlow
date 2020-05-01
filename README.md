### This workflow is designed for processing fastq files from Illumina sequencing of metabarcoding libraries. I.e. sequencing libraries built on pools of amplicons from complex samples (e.g. eDNA samples or bulk samples with DNA from several taxa), which have each been PCR-amplified using a unique combination of oligonucleotide tags on the primers. The workflow operates with amplicon sequence variants (ASVs) throughout, i.e. sequences are not collapsed into OTUs. The workflow includes demultiplexing, quality and error filtering, BLAST searching against GenBank's nt database, and taxonomic classification of the ASVs based on the BLAST hits. The main outputs of this workflow are an ASV table (which ASVs are found in which samples) and the taxonomic classification of these ASVs. Enjoy!

#### Make overall directories

  mkdir -p backup/data tmp results

#### Copy the scripts folder from the Github repository to the backup folder. 

#### Put the workflow.py file in the main directory, and the conda environment description in the backup folder

#### In the backup folder, add a readme file with explanations about the project. Ideally, put an appropriate readme file in the scripts and data folders as well

  touch backup/README.txt 

#### Add symbolic links to files and folders in backup for ease of use

  ln -s backup/scripts/ scripts
  
  ln -s backup/data/ data
  
  ln -s backup/environment.yml environment.yml
  
  ln -s backup/README.txt README.txt

#### Make a folder with a subfolder for each person involved in the (data analysis of) the project

  mkdir people

  mkdir "name1"

  mkdir "name2"

#### Duplicate root folder hierarchy in each person's folder
  
#### Create a conda environment based on the description file

  conda env create -f environment.yml

#### Install R package taxizedb. In the shell, run

  git config --global http.sslCAInfo /etc/ssl/certs/ca-bundle.crt # Network goes into a proxy, we need to give the certificate of the proxy to git        
  
  unset https_proxy # Unable the proxy anyway
  
  unset http_proxy

#### In R (still on the cluster, type R to start), run 

  devtools::install_github("ropensci/taxizedb")
  
  library("taxizedb")
  
  db_download_ncbi()

#### The taxizedb NCBI database should be updated regularly to keep up to date with the GenBank nt database (there seems to be some lag in the taxizedb online database) 

#### In the root data folder, download the raw sequencing data using wget directly in the terminal:

  nohup wget --continue --quiet URL &
 
#### Unzip the tar.gz or .tar file 

  tar -xzf filename.tar.gz
  
or
 
  tar -xvf filename.tar

#### Check md5 sum for the fastq.gz files in each library to make sure they are intact

md5sum -c MD5.txt
  
#### Unzip the fastq.gz files

  gunzip *gz
  
#### In each library data folder, make a file named tags.txt containing the sample names and corresponding forward and reverse tag sequences (see an example in data folder of this repository). Remember to put the library/PCR replicate number at the end of each sample name (e.g. "SAMPLE1_1" for library 1, "SAMPLE1_2" for library 2 and so on), so that PCR replicates will be kept separate when the data from the different libraries are merged. You can start by making the file for library 1 in excel, then save as tab separated file, transfer to the server, change from Windows to UNIX format (important!) and then use this file as a template for the remaining libraries (just replace replicate number in the sample names). Ensure that there is one (1) empty line at the end of each file. 

#### In each library data folder, make a file named batchfileDADA2.list containing the fastq file names, the primer sequences, and the minimum length required for a read (unaligned, i.e. forward or reverse read) after trimming of primers and tags (see example in data folder of this repository). If your primers contain inosine bases ("I"), these need to be replaced with "N", as the software does not recognize "I". Again, you can start by making the file for the first library, then copy it to the other libraries like this:

  echo lib1 lib2 | xargs -n 1 cp batchfileDADA2.list
  
#### Then change the fastq file names and the tag file name to the appropriate ones. Make sure that the file is in UNIX format. Again, ensure that there is one (1) empty line at the end of the file.

#### If appropriate, change the minimum length requirement in the match_pairs.r script. Check whether it would be appropriate to change any of the options set for the blastn command. In the taxonomy.r script, consider whether you want to keep hits to "uncultured" and "environmental" sequences and if so, adjust the "remove" parameter to change this. Also consider whether the upper and lower margins should be adjusted (see explanation in the script).    

#### In the workflow file, replace the project name and the path to the raw data with your own name and path. If appropriate, change the length and quality requirements provided to the sickle command. 

#### Activate the environment
  
  conda activate havblitz
  
#### Run gwf workflow from main folder

gwf run

#### Check the status of the workflow or a specific part using e.g.

gwf status taxonomy_96

#### Increase no. of cores, memory requirements and/or time limits if needed, or decrease if you need less resources. 

#### The outputs from this workflow that you will normally use for further analyses are primarily the ASV table of which unique sequences are in which samples (DADA2_nochim.table in tmp) and the taxonomic classification of these ASVs (classified.txt in tmp/taxonomy). Further analyses can be done on your local computer in R.

#### NB! Marie Lund at Microbiology recommends to do dada2 cleaning before demultiplexing! Perhaps this should be implemented
