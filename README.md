# Xpression_collect
Repository with scripts for data fetching from SRA and other downstream analyses to obtain comprehensive gene expression database

**kallisto_pipeline3.py**: wrapper script to run kallisto 

**merge_kallisto_ouput3.py**: script to merge kallisto output files

**fetch_sras_run_kallisto.py**: script to fetch SRA files in batches, process kallisto for the batch, obtain counts, tpm files for the batch, and delete the SRA batch for the next run

**filter_RNAseq_samples.py**: script to QC filter the expression files for genomic sequences

**isoform_purger.py**: script ot obtain CDS file with primary transcript

**transeq.py**: script to translate CDS sequences to PEP

**isoform_purger.py**: script to obtain CDS with primary transcript file, PEP with primary transscript file, and expression file with only primary transcripts


