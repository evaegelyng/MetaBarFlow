#### Make overall directories

  mkdir -p backup/scripts backup/data tmp results

#### Put the workflow.py file in the main directory, and the conda environment description in the backup folder

#### In the backup folder, add a readme file (and remember to put explanations about the project in it!)

  touch backup/environment.yml backup/README.txt

#### Add symbolic links to files and folders in backup for ease of use

  ln -s backup/scripts/ scripts
  
  ln -s backup/data/ data
  
  ln -s backup/templates.py templates.py
  
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

#### In R, run 

  devtools::install_github("ropensci/taxizedb")
  
  library("taxizedb")
  
  db_download_ncbi()

#### The taxizedb NCBI database should be updated regularly to keep up to date with the GenBank nt database (as much as possible - there seems to be some lag in the taxizedb online database) 

#### In the root data folder, download the raw sequencing data using wget directly in the terminal, e.g. using a screen or with:

  nohup wget --continue --quiet URL &
 
#### Unzip the tar.gz file 

  tar -xzf filename.tar.gz
  
#### Unzip the fastq.gz files

  gunzip *gz
  
#### In each library data folder, make a file containing the sample names and corresponding forward and reverse tag sequences. Remember to put the library/PCR replicate number at the end of each sample name (e.g. "SAMPLE1_1" for library 1, "SAMPLE1_2" for library 2 and so on), so that PCR replicates will be kept separate when the data from the different libraries are merged. You can start by making the file for library 1 in excel (Havblitz Bib1.xlsx), then save as tab separated file, transfer to the server, change from Windows to UNIX format (important!) and then use this file as a template for the remaining libraries (just replace replicate number in the sample names). Ensure that there is one (1) empty line at the end of each file. 

#### In each library data folder, make a file containing the fastq file names, the primer sequences, and the minimum length required for a read (unaligned, i.e. forward or reverse read) after trimming of primers and tags (see "batchfileDADA2.list" files in HavBlitz project). Again, you can start by making the file for the first library, then copy it to the other libraries like this:

  echo lib1 lib2 | xargs -n 1 cp batchfileDADA2.list
  
#### Then change the fastq file names and the tag file name to the appropriate ones. Make sure that the file is in UNIX format. Again, ensure that there is one (1) empty line at the end of the file.  

#### Add all scripts to the scripts folder. 

#### Activate your environment
  
  conda activate myproject
  
#### Run gwf workflow from main folder

gwf run

#### NB! Marie Lund at Microbiology recommends to do dada2 cleaning before demultiplexing!