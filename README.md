### This workflow allows efficient parallel processing of fastq files from Illumina sequencing of DNA metabarcoding libraries. I.e. sequencing libraries built on pools of amplicons from complex samples (e.g. eDNA samples or bulk samples with DNA from several taxa), which have each been PCR-amplified using a unique combination of oligonucleotide tags on the primers. The workflow operates with amplicon sequence variants (ASVs) throughout, i.e. sequences are not collapsed into OTUs. The workflow includes demultiplexing, quality and error filtering, BLAST searching, and taxonomic classification. Databases used for BLAST searching and taxonomic classification are downloaded locally, facilitating the analysis of large numbers of ASVs. The main outputs of the workflow are an ASV table (which ASVs are found in which samples) and the taxonomic classifications of these ASVs, based on a Last Common Ancestor (LCA) approach. Overlaps in sequence similarity to query sequences are used to determine which BLAST hits to include when assigning an LCA (Sigsgaard et al. 2020). Importantly, the automatic classifications always need to be checked carefully and manually, using general knowledge of taxonomy, local species occurrences, synonyms etc., in order to produce a more realistic final list of taxa. Enjoy!

#### Make overall directories

  `mkdir -p backup/data tmp results`

#### Copy the scripts folder and the conda environment description from the Github repository to the backup folder. 

#### Put the workflow.py file in the main directory.

#### In the backup folder, add a readme file with explanations about the project. Ideally, put an appropriate readme file in the scripts and data folders as well

  `touch backup/README.txt` 

#### Add symbolic links in the main folder to files and folders in backup

```
  ln -s backup/scripts/ scripts
  
  ln -s backup/data/ data
  
  ln -s backup/environment.yml environment.yml
  
  ln -s backup/README.txt README.txt
```
 
#### Create a conda environment based on the description file

`conda env create --name projectname -f environment.yml`
  
#### In Jensen et al., the package taxizedb was installed with devtools, as it was not yet on conda. If you want to reproduce exactly the pipeline in Jensen et al., see footnote for installation details. Otherwise, load your conda environment, and install the package taxizedb:

```  
  conda activate projectname
  
  conda install -c conda-forge r-taxizedb
```

#### If you cannot create the environment based on the description file (updated packages may cause problems), create your own environment, beginning with the packages that are directly called in the scripts (cutadapt, sickle, taxizedb etc.). It can be helpful to use mamba to install packages, as it is faster than conda.

   `conda install -c conda-forge mamba`
   
#### The taxizedb NCBI database should be updated regularly to keep up to date with the GenBank nt database (there seems to be some lag in the taxizedb online database) 

```
  R
  
  library("taxizedb")
  
  db_download_ncbi()
```

#### In the root data folder, download the raw sequencing data. If you have been provided with a csv file with links to the files you can run the following command to download all files at once:

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

#### Use the software fastqc to further inspect the quality of the data. If you do not have it in your base environment, install to your project environment:

```
  mamba install -c bioconda fastqc 
  
  sbatch YOUR_PATH/scripts/fastqc.sh
```
   
#### In each library data folder, make a tab separated file named tags.txt containing the sample names and corresponding forward and reverse tag sequences (see an example in data folder of this repository). Remember to put the library number/PCR replicate number at the end of each sample name (e.g. "SAMPLE1_1" for library 1, "SAMPLE1_2" for library 2 and so on. Check that none of the sample names themselves contain these endings, e.g. "SAMPLE_1"). This way, PCR replicates will be kept separate when the data from the different libraries are merged. You can start by making the file for library 1 in excel, transfer to the server, and then use this file as a template for the remaining libraries (just replace replicate number in the sample names). Note that these txt-files should be in UNIX format (not Windows - can be checked e.g. using notepad++). In some cases, it is also necessary to add an empty line at the end of each tags file, and to remove the tab separator ("\t") in all instances of "line.split("\t")" in the workflow.

#### The script create_batch.sh can be used to make a file (batchfileDADA2.list) in each library data folder containing the fastq file names, the primer sequences, and the minimum length required for a read (unaligned, i.e. forward or reverse read) after trimming of primers and tags. Replace the primer sequences and length specified in the script with those appropriate for your own project. If your primers contain inosine bases ("I"), these need to be replaced with "N", as the software does not recognize "I". 

#### If appropriate, change the minimum length requirement in the match_pairs.r script. Check whether it would be appropriate to change any of the options set for the blastn command, and add your own database path. A widely used blast database is the NCBI GenBank "nt" database, which can be downloaded accordingly (inside a fitting directory, and with a stable internet connection):

`update_blastdb.pl nt --timeout 500`

#### Update taxonomy for blast folder

```
update_blastdb.pl taxdb
tar -zvxf taxdb.tar.gz
```

#### For all "nt.XX.tar.gz" files, run the following:

`tar -zvxf nt.XX.tar.gz`

#### If using the NCBI taxonomy, create the table MergedTaxIDs to translate old, deprecated taxids to the corresponding new taxid. This file should be updated every so often to account for newly merged taxIDs. The file is generated from the merged.dmp file in the taxdb folder downloaded from NCBI:

#### Download and unzip the most recent taxdump folder (e.g.)

```
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/new_taxdump_2021-12-01.zip

unzip new_taxdump_2021-12-01.zip
```

#### Create the MergedTaxIDs file from the merged.dmp file

```
echo -e OldTaxID'\t'NewTaxID > MergedTaxIDs

less merged.dmp | cut -f1,3 >> MergedTaxIDs
```

#### In the taxonomy.r script, add your own path to the MergeTaxIDs table. Also, consider whether you for instance want to keep hits to "uncultured" and "environmental" sequences and if so, adjust the "remove" parameter to change this. Also consider whether the upper and lower margins should be adjusted (see explanation in the script).    

#### In the workflow file, replace the project name and the path to the raw data with your own. If appropriate, change the length and quality requirements provided to the sickle command. 

#### Create an API key for NCBI and paste it into the taxonomy.r script (replace "YOUR_KEY")

#### Run the gwf workflow from the main folder

  `gwf run`

#### If you get an error mentioning the backend, try to reconfigure this to slurm

  `gwf config set backend slurm`

#### Check the status of the workflow using 

  `gwf status` 

#### By adding the name of a specific target after the above command, you can see the status of this target. E.g:

 `gwf status demultiplex*` 

#### As the function splitting your fasta file of OTUs before BLAST searching may output a smaller number of files than the 99 files specified (it seems the software has a minimum threshold for the number of sequences that can go in each file), double-check in the .stderr log file that the number of sequences of the separate files add up to the total sequence number. Note that a hidden folder named ".gwf/logs" is where you will find your log files. Because the number of input files for BLAST searching is unknown until the split function has run, the remaining targets of the workflow can only be started once splitting is complete. To start the remaining targets, just use:

 `gwf run`

#### Increase no. of cores, memory requirements and/or time limits if needed, or decrease if you need less resources. You can check your realized resource use for a target using the package gwf-utilization:

```
   conda install -c micknudsen gwf-utilization
   
   gwf utilization
```

#### The outputs from this workflow that you will normally use for further analyses are primarily the ASV table of which unique sequences are in which samples (DADA2_nochim.table) and the taxonomic classification of these ASVs (classified.txt). Further analyses can be done on your local computer in R.

#### Remember to backup your raw data, metadata, scripts and conda environment(s), and final outputs!

#### Footnote on taxizedb installation with devtools:

```
  git config --global http.sslCAInfo /etc/ssl/certs/ca-bundle.crt # Network goes into a proxy, we need to give the certificate of the proxy to git        
  
  unset https_proxy
  
  unset http_proxy

  R 

  devtools::install_github("ropensci/taxizedb")
```

### Key Contributors

#### This pipeline was developed in the eDNA research group at the Department of Biology, Aarhus University, by:

#### [Eva Egelyng Sigsgaard](https://github.com/evaegelyng): Lead developer and maintainer since 2019

#### Philip Francis Thomsen (Principal Investigator): Scientific input, especially on taxonomic classification

#### [Adrián Gómez Repollés](https://github.com/adriangeerre): Developer, diverse contributions 

### Suggested Citation

#### Please link to this GitHub repository and refer to the publication: Jensen et al. Short-term temporal variation of coastal marine eDNA. Environmental DNA XX (2022).

### Acknowledgements

#### The scripts called by the gwf workflow were mainly written by [Tobias G. Frøslev](https://github.com/tobiasgf) (see Frøslev et al. 2017), and are to a large extent based on the DADA2 package (Callahan et al. 2016). Thanks to [Dan Søndergaard](https://github.com/dansondergaard) for help getting started with the [gwf workflow tool](https://docs.gwf.app/), and a special mention to [Samuele Soraggi](https://github.com/SamueleSoraggi) for assistance with Python scripting and troubleshooting. We thank GenomeDK at the Bioinformatic Research Center (BiRC), Aarhus University, for providing computational resources. This work was supported by the The Velux Foundations (grant 21517), the Carlsberg Foundation (grant CF18-0949) and The Faculty of Natural Sciences, Aarhus University (grant 27744).

### Citations (please see Jensen et al. for references to all software packages)

#### Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: high-resolution sample inference from Illumina amplicon data. Nature methods, 13(7), 581-583.

#### Frøslev, T. G., Kjøller, R., Bruun, H. H., Ejrnæs, R., Brunbjerg, A. K., Pietroni, C., & Hansen, A. J. (2017). Algorithm for post-clustering curation of DNA amplicon data yields reliable biodiversity estimates. Nature Communications, 8(1), 1188.

#### Sigsgaard, E. E., Olsen, K., Hansen, M. D., Hansen, O. L. P., Høye, T. T., Svenning, J. C., & Thomsen, P. F. (2020). Environmental DNA metabarcoding of cow dung reveals taxonomic and functional diversity of invertebrate assemblages. Molecular ecology, 30(13), 3374-3389.

### Questions

#### If you have questions or issues, please email Eva Egelyng Sigsgaard (eva.sigsgaard@bio.au.dk) or Mads R. Jensen (mrj@bio.au.dk), or leave a comment on this repository.
