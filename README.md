### This workflow is designed for processing fastq files from Illumina sequencing of metabarcoding libraries. I.e. sequencing libraries built on pools of amplicons from complex samples (e.g. eDNA samples or bulk samples with DNA from several taxa), which have each been PCR-amplified using a unique combination of oligonucleotide tags on the primers. The workflow operates with amplicon sequence variants (ASVs) throughout, i.e. sequences are not collapsed into OTUs. The workflow includes demultiplexing, quality and error filtering, BLAST searching against a local sequence database, and taxonomic classification of the ASVs based on the BLAST hits. The main outputs of this workflow are an ASV table (which ASVs are found in which samples) and the taxonomic classification of these ASVs. Enjoy!

#### Make overall directories

  `mkdir -p backup/data tmp results`

#### Copy the scripts folder from the Github repository to the backup folder. 

#### Put the appropriate workflow file (for COI data: workflow_bold_nt.py - otherwise workflow.py) in the main directory, and the conda environment description in the backup folder

#### In the backup folder, add a readme file with explanations about the project. Ideally, put an appropriate readme file in the scripts and data folders as well

  `touch backup/README.txt` 

#### Add symbolic links to files and folders in backup for ease of use

```
  ln -s backup/scripts/ scripts
  
  ln -s backup/data/ data
  
  ln -s backup/environment.yml environment.yml
  
  ln -s backup/README.txt README.txt
```
 
#### Create a conda environment based on the description file

  `conda env create --name projectname -f environment.yml`
  
#### If you cannot create the environment based on the description file (updated packages may cause problems), create your own environment, beginning with the packages that are directly called in the scripts (cutadapt, sickle, taxizedb etc.). If you still have trouble once you have installed these and their dependencies, check the list of packages specified in the file Clean_210922.yml. This file contains only the packages that were common between currently working environments of Mads, Martin and Eva. It can be helpful to use mamba to install packages, as it is faster than conda. 

#### The taxizedb NCBI database should be updated regularly to keep up to date with the GenBank nt database (there seems to be some lag in the taxizedb online database) 

```
  library("taxizedb")
  
  db_download_ncbi()
```

#### In the root data folder, download the raw sequencing data using wget and the csv file from Novogene ("Export link" on the data website):

  `wget -c --progress=dot:giga --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -i YOUR_CSV.csv`
 
#### Unzip the tar.gz or .tar file 

  `tar -xzf filename.tar.gz`
  
or
 
  `tar -xvf filename.tar`
  
#### After unzipping, remember to move the zipped data folder to an independent location, such as a portable hard drive for backup. This is to avoid unnecessary use of expensive server backup storage and to have a local backup, which is independent of the server, and faster to retrieve. 

#### Remember to check the sequencing report to get an overview of the quality and amount of data. If this is not satisfactory, consider asking for resequencing.

#### Check md5 sum for the fastq.gz files in each library to make sure they are intact

  `md5sum -c MD5.txt`
  
#### Unzip the fastq.gz files in each library

  `gunzip *gz`

or if you have many libraries, run the following for the entire raw data folder

  ``for i in `find . -name "*.gz"`; do gunzip $i; done &``  

#### Use the software fastqc (base environment) to further inspect the quality of the data:

  `sbatch --account eDNA YOUR_PATH/scripts/fastqc.sh`
   
#### In each library data folder, make a tab separated file named tags.txt containing the sample names and corresponding forward and reverse tag sequences (see an example in data folder of this repository). Remember to put the library number/PCR replicate number at the end of each sample name (e.g. "SAMPLE1_1" for library 1, "SAMPLE1_2" for library 2 and so on. Check that none of the sample names themselves contain these endings, e.g. "SAMPLE_1"). This way, PCR replicates will be kept separate when the data from the different libraries are merged. You can start by making the file for library 1 in excel, transfer to the server, change from Windows to UNIX format (important!) and then use this file as a template for the remaining libraries (just replace replicate number in the sample names). 

#### The script create_batch.sh can be used to make a file (batchfileDADA2.list) in each library data folder containing the fastq file names, the primer sequences, and the minimum length required for a read (unaligned, i.e. forward or reverse read) after trimming of primers and tags. Replace the primer sequences and length specified in the script with those appropriate for your own project. If your primers contain inosine bases ("I"), these need to be replaced with "N", as the software does not recognize "I". 

#### If appropriate, change the minimum length requirement in the match_pairs.r script. Check whether it would be appropriate to change any of the options set for the blastn command. In the taxonomy.r script, consider whether you for instance want to keep hits to "uncultured" and "environmental" sequences and if so, adjust the "remove" parameter to change this. Also consider whether the upper and lower margins should be adjusted (see explanation in the script).    

#### In the workflow file, replace the project name and the path to the raw data with your own. If appropriate, change the length and quality requirements provided to the sickle command. 

#### Activate the environment
  
  `conda activate YOUR_ENV`
  
#### Run gwf workflow from main folder (if using workflow_bold_nt.py, you need to specify this filename after the command)

  `gwf run`

#### If you get an error mentioning the backend, try to reconfigure this to slurm

  `gwf config set backend slurm`

#### Check the status of the workflow using 

  `gwf status` 

#### By adding the name of a specific target after the above command, you can see the status of this target.

#### As the function splitting your fasta file of OTUs before BLASTing may output a smaller number of files than the 99 files specified (it seems the software has certain thresholds for the number of sequences that can go in each file), double-check in the .stderr log file that the number of sequences of the separate files add up to the total sequence number.

#### Increase no. of cores, memory requirements and/or time limits if needed, or decrease if you need less resources. You can check your realized resource use for a target using the package gwf-utilization:

```
   conda install -c micknudsen gwf-utilization
   
   gwf utilization
```

#### The outputs from this workflow that you will normally use for further analyses are primarily the ASV table of which unique sequences are in which samples (DADA2_nochim.table) and the taxonomic classification of these ASVs (classified.txt). Further analyses can be done on your local computer in R.

#### Remember to backup your raw data, metadata, scripts and conda environment(s), and final outputs!
