from gwf import Workflow
import os, sys
import math
from glob import glob

project_name = "havblitz"

#gwf = Workflow(defaults={"account": "edna"}) #Should give higher job priority, but does not work - remind Dan
gwf = Workflow()

#Demultiplex
libraries = [x for x in glob("data/X201SC19122834-Z01-F001/raw_data/*") if os.path.isdir(x)]

for library_root in libraries:
    library_id = os.path.basename(library_root)
    input_files = glob(library_root + "/*.fq")
    
    with open(os.path.join(library_root, "tags.txt")) as tags_file:
        for line in tags_file:
            output_files = []
            tag_id, fseq, rseq = line.split()

            output_files.append("tmp/{}/DADA2_AS/{}_R1.fastq".format(library_id, tag_id))
            output_files.append("tmp/{}/DADA2_AS/{}_R2.fastq".format(library_id, tag_id))
            output_files.append("tmp/{}/DADA2_SS/{}_R1.fastq".format(library_id, tag_id))
            output_files.append("tmp/{}/DADA2_SS/{}_R2.fastq".format(library_id, tag_id))

            gwf.target(
                name="demultiplex_{}_{}_{}".format(project_name,library_id, tag_id),
                inputs=input_files,
                outputs=output_files,
                cores=1,
                memory="4g",
                walltime="1:00:00",
            ) << """
                mkdir -p tmp/{library_id}
                ./scripts/demultiplex.sh {library_root} tmp/{library_id} {tag_id} {tag_fseq} {tag_rseq}
            """.format(library_root=library_root, library_id=library_id, tag_id=tag_id, tag_fseq=fseq, tag_rseq=rseq) 
            
#Quality trimming of reads
for library_root in libraries:
    library_id = os.path.basename(library_root)
    
    with open(os.path.join(library_root, "tags.txt")) as tags_file:
        for line in tags_file:
            input_files = []
            tag_id, fseq, rseq = line.split()
    
            input_files.append("tmp/{}/DADA2_AS/{}_R1.fastq".format(library_id, tag_id))
            input_files.append("tmp/{}/DADA2_AS/{}_R2.fastq".format(library_id, tag_id))
            input_files.append("tmp/{}/DADA2_SS/{}_R1.fastq".format(library_id, tag_id))
            input_files.append("tmp/{}/DADA2_SS/{}_R2.fastq".format(library_id, tag_id))
    
            output_files = []
            tag_id, fseq, rseq = line.split()

            output_files.append("tmp/{}/DADA2_AS/filtered/{}_F_filtered.fastq".format(library_id, tag_id))
            output_files.append("tmp/{}/DADA2_AS/filtered/{}_R_filtered.fastq".format(library_id, tag_id))
            output_files.append("tmp/{}/DADA2_SS/filtered/{}_F_filtered.fastq".format(library_id, tag_id))
            output_files.append("tmp/{}/DADA2_SS/filtered/{}_R_filtered.fastq".format(library_id, tag_id))
            
            folderAS="tmp/{}/DADA2_AS/filtered".format(library_id)
            folderSS="tmp/{}/DADA2_SS/filtered".format(library_id)            
            
            gwf.target(
                name="sickle_{}_{}_{}".format(project_name,library_id, tag_id),
                inputs=input_files,
                outputs=output_files,
                cores=1,
                memory="4g",
                walltime="1:00:00",
            ) << """
                mkdir -p {folderAS}
                mkdir -p {folderSS}
                sickle se -l 50 -q 28 -x -t sanger -f {inputASF} -o {outputASF}
                sickle se -l 50 -q 28 -x -t sanger -f {inputASR} -o {outputASR}
                sickle se -l 50 -q 28 -x -t sanger -f {inputSSF} -o {outputSSF}
                sickle se -l 50 -q 28 -x -t sanger -f {inputSSR} -o {outputSSR} 
            """.format(folderAS=folderAS,folderSS=folderSS,inputASF=input_files[0],inputASR=input_files[1],inputSSF=input_files[2],inputSSR=input_files[3],outputASF=output_files[0],outputASR=output_files[1],outputSSF=output_files[2],outputSSR=output_files[3])
       
#Matching paired reads
for library_root in libraries:
    library_id = os.path.basename(library_root)
    
    with open(os.path.join(library_root, "tags.txt")) as tags_file:
        for line in tags_file:
            input_files = []
            tag_id, fseq, rseq = line.split()
    
            input_files.append("tmp/{}/DADA2_AS/filtered/{}_F_filtered.fastq".format(library_id, tag_id))
            input_files.append("tmp/{}/DADA2_AS/filtered/{}_R_filtered.fastq".format(library_id, tag_id))
            input_files.append("tmp/{}/DADA2_SS/filtered/{}_F_filtered.fastq".format(library_id, tag_id))
            input_files.append("tmp/{}/DADA2_SS/filtered/{}_R_filtered.fastq".format(library_id, tag_id))
    
            output_files = []
            tag_id, fseq, rseq = line.split()
            
            output_files.append("tmp/{}/DADA2_AS/filtered/matched/{}_F_matched.fastq.gz".format(library_id, tag_id))
            output_files.append("tmp/{}/DADA2_AS/filtered/matched/{}_R_matched.fastq.gz".format(library_id, tag_id))
            output_files.append("tmp/{}/DADA2_SS/filtered/matched/{}_F_matched.fastq.gz".format(library_id, tag_id))
            output_files.append("tmp/{}/DADA2_SS/filtered/matched/{}_R_matched.fastq.gz".format(library_id, tag_id))
                        
            folderAS="tmp/{}/DADA2_AS/filtered/matched".format(library_id)
            folderSS="tmp/{}/DADA2_SS/filtered/matched".format(library_id)            
            
            gwf.target(
                name="match_{}_{}_{}".format(project_name,library_id, tag_id),
                inputs=input_files,
                outputs=output_files,
                cores=1,
                memory="4g",
                walltime="1:00:00",
            ) << """
                mkdir -p {folderAS}
                mkdir -p {folderSS}
                Rscript ./scripts/match_pairs.r {inputASF},{inputASR},{inputSSF},{inputSSR} {outputASF},{outputASR},{outputSSF},{outputSSR} 
                 if grep -q "removed" ".gwf/logs/match_{project_name}_{library_id}_{tag_id}.stderr"
                 then
                  echo "" > "tmp/{library_id}/DADA2_AS/filtered/matched/{tag_id}_F_matched.fastq.gz"
                  echo "" > "tmp/{library_id}/DADA2_AS/filtered/matched/{tag_id}_R_matched.fastq.gz"
                  echo "" > "tmp/{library_id}/DADA2_SS/filtered/matched/{tag_id}_F_matched.fastq.gz"
                  echo "" > "tmp/{library_id}/DADA2_SS/filtered/matched/{tag_id}_R_matched.fastq.gz"
                 fi
            """.format(folderAS=folderAS,folderSS=folderSS,inputASF=input_files[0],inputASR=input_files[1],inputSSF=input_files[2],inputSSR=input_files[3],outputASF=output_files[0],outputASR=output_files[1],outputSSF=output_files[2],outputSSR=output_files[3],project_name=project_name,library_id=library_id,tag_id=tag_id)         

#Removing likely erroneous sequences 
for library_root in libraries:
    library_id = os.path.basename(library_root)
    
    with open(os.path.join(library_root, "tags.txt")) as tags_file:
        for line in tags_file:
            input_files = []
            tag_id, fseq, rseq = line.split()

            input_files.append("tmp/{}/DADA2_AS/filtered/matched/{}_F_matched.fastq.gz".format(library_id, tag_id))
            input_files.append("tmp/{}/DADA2_AS/filtered/matched/{}_R_matched.fastq.gz".format(library_id, tag_id))
            input_files.append("tmp/{}/DADA2_SS/filtered/matched/{}_F_matched.fastq.gz".format(library_id, tag_id))
            input_files.append("tmp/{}/DADA2_SS/filtered/matched/{}_R_matched.fastq.gz".format(library_id, tag_id))
            
    output_files = []
    
    output_files.append("tmp/{}/seqtab_AS_RDS".format(library_id))
    output_files.append("tmp/{}/seqtab.nochim_AS_RDS".format(library_id))
    output_files.append("tmp/{}/seqtab_SS_RDS".format(library_id))
    output_files.append("tmp/{}/seqtab.nochim_SS_RDS".format(library_id))
                                
    gwf.target(
      name="remove_errors_{}_{}".format(project_name,library_id),
      inputs=input_files,
      outputs=output_files,
      cores=4,
      memory="4g",
      walltime="1:00:00",
    ) << """
      Rscript ./scripts/remove_errors.r tmp/{library_id}
      """.format(library_id=library_id)         
            
#Summing sense and antisense sequence tables
for library_root in libraries:
    library_id = os.path.basename(library_root)
    input_files = []
    
    input_files.append("tmp/{}/seqtab_AS_RDS".format(library_id))
    input_files.append("tmp/{}/seqtab.nochim_AS_RDS".format(library_id))
    input_files.append("tmp/{}/seqtab_SS_RDS".format(library_id))
    input_files.append("tmp/{}/seqtab.nochim_SS_RDS".format(library_id))
    
    output_files = []
    
    output_files.append("tmp/{}/seqtab_RDS".format(library_id))
    output_files.append("tmp/{}/seqtab.nochim_RDS".format(library_id))
    
    gwf.target(
      name="sum_AS_SS_{}_{}".format(project_name,library_id),
      inputs=input_files,
      outputs=output_files,
      cores=1,
      memory="2g",
      walltime="1:00:00",
    ) << """
      Rscript ./scripts/sum_AS_SS.r tmp/{library_id}
      """.format(library_id=library_id)    

#Summing sequence tables of all libraries
input_files = []
for library_root in libraries:
    library_id = os.path.basename(library_root)   
    input_files.append("tmp/{}/seqtab_RDS".format(library_id))
    input_files.append("tmp/{}/seqtab.nochim_RDS".format(library_id))
 
output_files = []
    
output_files.append("tmp/seqtab_Both")
output_files.append("tmp/seqtab.nochim_Both")
output_files.append("tmp/DADA2_raw.table")
output_files.append("tmp/DADA2_raw.otus")
output_files.append("tmp/DADA2_nochim.table")
output_files.append("tmp/DADA2_nochim.otus")
    
gwf.target(
   name="sum_libraries_{}".format(project_name),
   inputs=input_files,
   outputs=output_files,
   cores=1,
   memory="2g",
   walltime="1:00:00",
 ) << """
   Rscript ./scripts/sum_libraries.r tmp/
   """.format(library_id=library_id)                

#BLAST search and taxonomic assignment

###Split fasta file (the nochim one with chimeras removed) into K parts
def splitter(inputFile, K=100):
    inputs = [inputFile]
    outFiles = [inputFile + '.split/'+'DADA2_nochim.part_'+'{:0>3d}'.format(i)+'.fasta' for i in range(1,K+1)]
    outIndex = inputFile + '.seqkit.fai'
    outputs = outFiles + [outIndex]
    options = {
        'cores': 1,
        'memory': '2g',
        'walltime': '1:00:00'
    }
    spec = '''
    seqkit split {inputFile} -p {K} -2
    '''.format(inputFile=inputFile, K=K)
    return inputs, outputs, options, spec

#####blast a single k-th file
def blaster(fileName, k, outFolder):
    inputFasta = fileName+'.split/'+'DADA2_nochim.part_'+'{:0>3d}'.format(k)+'.fasta'
    inputs = [inputFasta]
    outBlast = outFolder + '/blast.' + str(k) + '.blasthits'
    outLog = outFolder + '/blast.' + str(k) + '.txt'
    outputs = [
      outBlast,
      outLog
    ]
    options = {
        'cores': 4,
        'memory': '32g',
        'walltime': '4:00:00'
    }
    spec = '''
    export BLASTDB=/faststorage/project/eDNA/blastdb/nt/
    mkdir -p {out}/blast

    echo "RUNNING THREAD {k} BLAST"
    blastn -db /faststorage/project/eDNA/blastdb/nt/nt -max_target_seqs 500 -num_threads 4 -outfmt "6 std qlen qcovs sgi sseq sscinames staxid" -out {outBlast} -qcov_hsp_perc 90 -perc_identity 80 -query {inputFasta}
    echo "hello" > {outLog}
    echo "DONE THREAD {k}"
    '''.format(out=outFolder, k=k, inputFasta=inputFasta, outBlast=outBlast, outLog=outLog)
    return inputs, outputs, options, spec

def taxonomy(fileName, taxonomyFolder, blastFolder, k):
    inputFile = blastFolder + '/blast.' + str(k) + '.blasthits'
    inputs = [inputFile , blastFolder + '/blast.' + str(k) + '.txt']
    outputFile = taxonomyFolder + '/taxonomy.' + str(k) + '.txt'
    outputs = [outputFile]
    options = {
        'cores': 1,
        'memory': '32g',
        'walltime': '4:00:00'
    }
    
    spec = '''
    mkdir -p {taxonomyFolder}
    Rscript scripts/taxonomy.r {inputFile} {outputFile}
    '''.format(taxonomyFolder=taxonomyFolder, inputFile=inputFile, outputFile=outputFile) 
    return inputs, outputs, options, spec

K=100
inputName = 'tmp/DADA2_nochim.otus'

gwf.target_from_template( 'split', splitter(inputFile=inputName, K=K) )
                                                                
for k in range(1,K+1):
  gwf.target_from_template( 'blaster_{}'.format(k), blaster(fileName=inputName, k=k, outFolder='tmp/blast') )
  gwf.target_from_template( 'taxonomy_{}'.format(k), taxonomy(fileName=inputName, taxonomyFolder='tmp/taxonomy', blastFolder='tmp/blast', k=k) )
