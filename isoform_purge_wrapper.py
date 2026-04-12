# Shakunthala Natarajan #
## bug reports contact: s64snata@uni-bonn.de ##
### v0.1 ###

__usage__ = """
					python3 isoform_purge_wrapper.py
					--cds <FULL_PATH_TO_CDS_FILE_WITH_ISOFORMS_WITH_THE_EXTENSION .cds.fasta>
					--exp <FULL_PATH_TO_TPM_OR_COUNTS_EXPRESSION_TABLE>
					--purger <FULL_PATH_TO_isoform_purger.py_INCLUDING_SCRIPT_NAME>
					--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
					"""

# --- begin imports --- #
import os, sys, glob, re
from pathlib import Path
import subprocess
import Pandas as pd
# --- end of imports --- #

def keep_primary_transcript_cds_pep(cds_file, purge_script, transeq_script,outdir):
	species = ((os.path.basename(cds_file)).split(".cds"))[0]
	cmd = f'python3 {purge_script} --in {cds_file} --out {outdir}'
	p = subprocess.Popen(args=cmd, shell=True)
	p.communicate()
	cds_file_purged = (outdir.glob("*.fasta"))
	cmd = f'mv {cds_file_purged} {species}_primary_transcript.cds.fasta'
	p = subprocess.Popen(args=cmd, shell=True)
	p.communicate()
	pepfile_purged = os.path.join(outdir, f'{species}_primary_transcript.pep.fasta')
	cds_file_purged_renamed = (outdir.glob("*.cds.fasta"))
	cmd = f'python3 {transeq_script} --in {cds_file_purged_renamed} --out {pepfile_purged}'
	p = subprocess.Popen(args=cmd, shell=True)
	p.communicate()
	return species, cds_file_purged_renamed

def keep_primary_transcript_exp(species, primary_transcript_cds_file, exp_file):
	primary_transcripts = []
	with open (primary_transcript_cds_file, 'r') as f:
		for line in f:
			if line.startswith(">"):
				primary_transcripts.append(line[1:].strip())
	df = pd.read_csv(exp_file)
	filtered_df = df[df['gene'].isin(primary_transcripts)]
	if '.tpm' in exp_file:
		filtered_df.to_csv(f'{species}_primary_transcript.tpm.tsv', sep="\t", index=False)
	elif '.count' in exp_file:
		filtered_df.to_csv(f'{species}_primary_transcript.counts.tsv', sep="\t", index=False)

def main (arguments):
	cds_file = arguments[arguments.index('--cds')+1]
	purge_script = arguments[arguments.index('--purger')+1]
	transeq_script = arguments[arguments.index('--transeq')+1]
	exp_file = arguments[arguments.index('--exp')+1]
	outdir = arguments[arguments.index('--out')+1]
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	species, primary_transcript_cds_file = keep_primary_transcript_cds_pep(cds_file, purge_script, transeq_script,outdir)
	keep_primary_transcript_exp(species, primary_transcript_cds_file, exp_file)


if '--exp' in sys.argv and '--out' in sys.argv and '--cds' in sys.argv and '--purger' in sys.argv and '--transeq' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )