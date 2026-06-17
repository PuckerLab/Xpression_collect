
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

-> **v0.15 is out and is the latest stable release!**

-> **v0.15 updates include**:

  -> buscolineage flag option entry has been simplified
  
  -> SRA file fetching and kallisto runs take place in an interleaved manner saving processing time as batch size and total sample size increase
  
  -> Better error handling and network issue tolerant fall backs in SRA file fetching
  
  -> You can now merge already exisiting TPM/ repr_TPM files with the ones produced in the current run of Xpression_collector

## Workflow

<img width="5673" height="6723" alt="Xpression_collector" src="https://github.com/user-attachments/assets/59649a82-1531-4ba5-accf-47936399fefd" />

## Installation and dependencies

### (1) Manual installation

```
git clone https://github.com/PuckerLab/Xpression_collect
```
**Mandatory dependencies**
- **Tools** - kallisto
- **Python libraries** - pandas (v2.3.1 or greater), matplotlib (v3.10.5 or greater)

**Optional dependencies**
- **Tools** - BUSCO, sratoolkit, BLAST, MAFFT
- **Python libraries** - pyfiglet (v1.0.2 or greater)

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
                               --sra [<TXT_FILE_WITH_SRA_ACCESSIONS_LIST> | --readfiles <FOLDER_WITH_SRA_ACCESSION_SUBFOLDERS>]
                               --out <OUTPUT_DIR>
MANDATORY:

--cds                  STR     Full path to CDS file

Either provide a list of SRA accessions to fetch or path to already available main folder with subfolders named according to the samples/ accessions

--sra                  STR     Full path to TXT file with one SRA accession per line

--readfiles            STR     Full path to folder with SRA accession subfolders each containing FASTQ files

--out                  STR     Full path to output directory


OPTIONAL:

--sample_name          STR     Name of the species; default is sample

--uncompressed         STR     Provide this flag if your read files are uncompressed

--annotation_qc        STR     yes or no for BUSCO-based QC of the PEP file; default is yes

--threads              STR     Total number of cores for running the pipeline; default is 4

--batch_size           STR     Number of SRA accessions to be fetched per batch; default is 1

--attempts             STR     Number of attempts at prefetching and accession from SRA; default is just 1 attempt

--wait                 STR     Base wait time for sleep in case of network delays for prefetch; Increases exponentially
                               with a base of 2 for each reattempt; default is 20 seconds

--remove_isoforms      STR   Optional step to remove isoforms; yes or no; default is yes

--min                  INT   Minimum percentage expression of top 100 genes

--max                  INT   Maximum percentage expression of top 100 genes

--black                STR   SRA IDs to be removed or blacklisted in a TXT file with one SRA accession ID per line

--scorecut             INT   BLAST bit score cutoff for isoform purging; default is 100

--simcut               INT   BLAST similarity cutoff for isoform purging; default is 99.0

--lencut               INT   Length cutoff for isoform purging; default is 100

--snvcut               INT   Number of single nucleotide variants allowed between two nucleotide sequences to group them as isoforms or not; default is 5

--blast                STR   Full path to BLAST aligner

--eval                 STR   evalue cutoff for self BLAST used in isoform purging; default is 1e-10

--mafft                STR    Full path to MAFFT

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

```


# Xpression_collect
Repository with scripts for data fetching from SRA and other downstream analyses to obtain comprehensive gene expression database

<p align="center">
<img width="600" alt="xpression_collect drawio(1)" src="https://github.com/user-attachments/assets/aa647d2b-a4db-416d-868c-c257ce05a163" />
</p>

**fetch_sras_run_kallisto.py**: script to fetch SRA files in batches, process kallisto for the batch, 
obtain counts, tpm files for the batch, and delete the SRA batch for the next run

**kallisto_pipeline3.py**: wrapper script to run kallisto 

**merge_kallisto_ouput3.py**: script to merge kallisto output files

```
python3 fetch_sras_run_kallisto.py --sra /vol/data/A_thaliana_SRA.txt \
--kallisto /vol/data/scripts/kallisto_pipeline3.py \
--kallisto_merge /vol/data/scripts/merge_kallisto_output3.py \
--cds /vol/data/datasets/A_thaliana.cds.fasta \
--threads 20 --batch_size 10 \
--out /vol/data/results/Athaliana
```

**combine_tpm_robust.py**: script to combine multiple TPM files into a single TPM file. 
Individual TPM files must end with the extension TPM.txt; model command to use the script is:

```
python3 combine_tpm_robust.py -i /vol/data/results/Athaliana \
-o /vol/data/results/Athaliana/Athaliana_exp
```

 **combine_counts_robust.py**: script to combine multiple counts file into a single counts file. 
 Individual counts file must end with the extension Counts.txt; model command to use the script is:

```
python3 combine_counts_robust.py -i /vol/data/results/Athaliana/Athaliana_exp \
-o /vol/data/results/Athaliana/Athaliana_exp
```

**filter_RNAseq_samples.py**: script to QC filter the expression files for genomic sequences

```
python3 filter_RNAseq_samples.py --tpms /vol/data/results/Athaliana_exp/Combined_tpm.tsv \
--counts /vol/data/results/Athaliana/Athaliana_exp/Combined_counts.tsv \
--out /vol/data/results/Athaliana/Athaliana_exp/Arabidopsis_thaliana.tpms.tsv
```

**isoform_purge_wrapper.py**: script to obtain CDS file, PEP file, and expression file with only primary transcripts

**isoform_purger.py**: script ot obtain CDS file with primary transcript

**transeq.py**: script to translate CDS sequences to PEP

```
python3 isoform_purge_wrapper.py --cds /vol/data/datasets/A_thaliana.cds.fasta \
--exp /vol/data/results/Athaliana/Athaliana_exp/Arabidopsis_thaliana.tpms.tsv \
--purger /vol/data/scripts/isoform_purger.py \
--transeq /vol/data/scripts/transeq.py \
--out /vol/data/results/Athaliana/Athaliana_exp
```


