
<img width="524" height="252" alt="xpression_collector_logo" src="https://github.com/user-attachments/assets/ba8bb939-0c83-4724-a82f-524acdb64910" />

# Xpression_collector - an end-to-end pipeline for RNA-seq data processing

Xpression_collector is a simple Python-based pipeline for fetching RNA-seq reads from SRA and processing them to get the expression TPM files for further downstream analyses. It can also work with locally fetched SRA folders already existing in the user's system, circumventing the fetching step as needed. Following are the salient features of the pipeline:

-> Optional BUSCO analysis of the PEP file obtained from the input CDS file

-> Reattempting capability to handle network disruptions during SRA prefetch and fasterq-dump

-> Batch processing of SRA files from SRA and iterative removal of SRA files after kallisto quantification of every batch to free-up storage space

-> Data processing from locally available SRA folders 

-> Gene expression quantification by pseudopalignment using kallisto

-> Removal of invalid RNA-seq samples

-> Optional removal of isoforms

<h2>📢 News & Updates</h2>

-> **v0.2 is out and is the latest stable release!**

-> **v0.2 updates include**:

  -> buscolineage flag option entry has been simplified
  
  -> SRA file fetching and kallisto runs take place in an interleaved manner improving processing time
  
  -> Better error handling and network issue tolerant fall backs in SRA file fetching
  
  -> You can now merge already exisiting TPM/ repr_TPM files with the ones produced in the current run of Xpression_collector

## Workflow

<img width="5676" height="6723" alt="Xpression_collector" src="https://github.com/user-attachments/assets/0978c426-ba5a-463e-aa2b-5365d2cb852d" />


## Installation and dependencies

### (1) Manual installation

```
git clone https://github.com/PuckerLab/Xpression_collect
```
**Mandatory dependencies**
- **Tools** - kallisto, sratoolkit
- **Python libraries** - pandas (v2.3.1 or greater), matplotlib (v3.10.5 or greater)

**Optional dependencies**
- **Tools** - BUSCO 
- **Python libraries** - pyfiglet (v1.0.2 or greater), rich (v15.0.0 or greater)

### (2) Docker installation

```
docker pull shakunthalan/xpression_collector:latest
```

It is recommended to run the docker image as a user and not root:

```
 docker run --rm -u $(id -u) -v /path/to/data:/path/to/data shakunthalan/xpression_collector:latest
```

**Note:** If you are using a docker image of the tool, you need not specify the dependencies' full paths while running the tool; The dependencies are built-in in the docker image and that makes it simple to use the tool across systems without the need for local installations.

### (3) Installation in a conda environment

This method of installation installs all the dependencies in a conda environment using the environment.yml file in this repository

```
git clone https://github.com/PuckerLab/Xpression_collect

cd Xpression_collect

conda env create -f environment.yml

conda activate xpression_collector
```

## Running the pipeline

```
Usage:

python3 Xpression_collector.py --cds <CDS_FILE>
                               --sra <TXT_FILE_WITH_SRA_ACCESSIONS_LIST>
                               --out <OUTPUT_DIR>
MANDATORY:

--cds                  STR     Full path to CDS file

Either provide a list of SRA accessions to fetch or path to already available main folder with subfolders named according to the samples/ accessions

--sra                  STR     Full path to TXT file with one SRA accession per line


--out                  STR     Full path to output directory


OPTIONAL:

--sample_name          STR     Name of the species; default is sample

--min_sra_file_size    INT    Minimum file size cutoff in MB to check prefetched SRA file sizes to cath sralite files; default cutoff is 1MB

--annotation_qc        STR     yes or no for BUSCO-based QC of the PEP file; default is yes

--threads              STR     Total number of cores for running the pipeline; default is 4

--parallel_prefetch    STR     Number of SRA accessions to be prefetched parallely; default is 2

--attempts             STR     Number of attempts at prefetching and accession from SRA; default is just 1 attempt

--wait                 STR     Base wait time for sleep in case of network delays for prefetch; Increases exponentially
                               with a base of 2 for each reattempt; default is 20 seconds

--remove_isoforms      STR   Optional step to remove isoforms; yes or no; default is yes

--merge_tpms           STR   Full path to config file containing the full paths to the filtered_tpm, and/or repr_filtered_tpms to be merged

--min                  INT   Minimum percentage expression of top 100 genes

--max                  INT   Maximum percentage expression of top 100 genes

--black                STR   SRA IDs to be removed or blacklisted in a TXT file with one SRA accession ID per line

--busco                STR    Full path to BUSCO; Specify busco_docker if BUSCO is installed via docker

--busco_lineage        STR    Full path to the config file to specify the BUSCO lineage> <Tab separated TXT file with sample name
                              (should be the same as sample_name) in the first column and BUSCO lineage in the second column

--busco_version        STR    Version of BUSCO> default is v6.0.0

--container_version    STR    Container version of the BUSCO docker image

--docker_host_path     STR    Host path to be mounted on for running the docker image

-docker_container_path STR    Container path of the docker image

--organism_type        STR    eukaryote or prokaryote for BUSCO analysis; default is eukaryote

--kallisto             STR    Full path to the kallisto executable

--prefetch             STR    Full path to the prefetch executable

--fasterq-dump         STR    Full path to the fasterq-dump executable

--uncompressed         STR     Provide this flag if your read files are uncompressed

--clean_up             STR    Clean up the TMP folder after a successful run; default is yes

```
## More details

-> The pipeline has the ability to restart after disruptions. Simply repeat the command you used before the disruption

-> The TMP folder gets cleaned up after a successful run of the pipeline. You can retain it in case you want to take a deeper look

-> --parallel_prefetch flag determines the number of parallel prefetch processes at a time

  The optimal number of parallel prefetches cannot be defined precisely as it is influenced by the network bandwidth and the disk space available to store  .sra files

  But some performance optimizations have shown that a prefetch process uses anywhere between 30-45% of a CPU

  Based on this the pipeline runs two parallel prefetches by default with a dynamic queueing mechanism
  
  In the dynamic queueing and pushing design, once an accession's prefetch is done, it is pushed to the next step of fasterq-dump+pigz 
  and kallisto which function serially since both of these are CPU-bound tasks showing maximum performance when utilizing all available 
  cores at a time. The vacant spot in the prefetch block is now taken by the next accession in line making it a combinatorial scheme 
  utilizing parallelization and serialization as necessary

  Hence it is recommended to use the default parameter for --parallel_prefetch flag unless you are sure about abundance of storage space 
  and excellent network bandwidth

-> Alternative isoforms removal needs the GFF file to be given along with the CDS file and it is important for the CDS FASTA headers and the transcript identifiers in the GFF file to match to process without errors
  
-> Detailed guidelines on preparing the config file for --gff_config flag can be found in https://github.com/ShakNat/DupyliCate

-> In case the headers in your CDS file are formatted differently when compared to the transcript identifiers in the GFF file supplied, you can use a helper called Fasta_fix.py to tackle this issue. The helper script and its detailed usage can be found in https://github.com/ShakNat/DupyliCate





