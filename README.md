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


