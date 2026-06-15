### Shakunthala Natarajan ###
### bug reports: s64snata@uni-bonn.de ###

__version__=0.15
__usage__="""
			python3 Xpression_collector.py
			--cds <Full path to CDS file>
			--sample_name <Name of the sample you are analyzing>
			--sra <Full path to TXT file with one SRA accession per line>
			--readfiles <Full path to folder with SRA accession subfolders each containing FASTQ files>
			--uncompressed <Provide this flag if your read files are uncompressed>
			--annotation_qc <yes or no for BUSCO-based QC of the PEP file> default is yes
			--threads <Total number of cores for running the pipeline> default is 4
			--batch_size <Number od SRA accessions to be fetched per batch> default is 1
			--attempts <Number of attempts at prefetching and accession from SRA> default is just 1 attempt
			--wait <Base wait time for sleep in case of network delays for prefetch; Increases exponentially with a base of 2 for each reattempt> default is 20 seconds
			--min_sra_file_size <Minimum file size cutoff in MB to check prefetched SRA file sizes to cath sralite files that might cause downstream errors> default cutoff is 1MB
			--remove_isoforms <Optional step to remove isoforms>< yes or no> default is yes
			--merge_tpms <Full path to config file containing the full paths to the filtered_tpm, and/or repr_filtered_tpms to be merged>
						 <First column will have paths to filtered_tpm files one per line; Second column will have paths to filtered_repr_tpm files one per line>
			--min <Minimum percentage expression of top 100 genes>
			--max <Maximum percentage expression of top 100 genes>
			--black <SRA IDs to be removed or blacklisted in a TXT file with one SRA accession ID per line>
			--scorecut <BLAST bit score cutoff for isoform purging> default is 100
			--simcut <BLAST similarity cutoff for isoform purging> default is 99.0
			--lencut <Length cutoff for isoform purging> default is 100
			--snvcut <Number of single nucleotide variants allowed between two nucleotide sequences to group them as isoforms or not> default is 5
			--blast <Full path to BLAST aligner>
			--eval <evalue cutoff for self BLAST used in isoform purging> default is 1e-10
			--mafft <Full path to MAFFT>
			--busco <Full path to BUSCO> <Specify busco_docker if BUSCO is installed via docker>
			--busco_lineage <Specify the BUSCO lineage> default is auto
			--busco_version <Version of BUSCO> default is v6.0.0
			--container_version <Container version of the BUSCO docker image>
			--docker_host_path <Host path to be mounted on for running the docker image>
			-docker_container_path <Container path of the docker image>
			--organism_type <eukaryote or prokaryote> default is eukaryote
			--kallisto <Full path to the kallisto executable>
			--prefetch <Full path to the prefetch executable>
			--fasterq-dump <Full path to the fasterq-dump executable>
			--out <Full path to output directory>
			"""

### --- start imports --- ###
import os,sys,glob, re, time
import gzip
import subprocess
import tempfile
import copy
import logging
import shutil
import pandas as pd
from functools import reduce
from pathlib import Path
from collections import Counter, OrderedDict
import matplotlib.pyplot as plt
from threading import Thread
from queue import Queue
from operator import itemgetter
try:
	from pyfiglet import Figlet
	f = Figlet(font='doom',width=200)
	print(f.renderText('Xpression collector'))
except ImportError:
	pass
### --- end of imports --- ###

fetch_queue = Queue()
barrier = None

#global definition of dictionary with problematic characters for dendropy processing, busco processing and their respective placeholders
replacements = {
		"'": "§quo",
		":": "§col",
		",": "§com",
		"(": "§lbr",
		")": "§rbr",
		";": "§sco",
		"[": "§lsbr",
		"]": "§rsbr",
		"\\": "§bsl",
		"\"": "§dquo",
		"/": "§fsl"
		}

def load_multiple_fasta_file(fasta_file):
	"""Load all sequences from a (possibly wrapped) FASTA file into a dict."""
	content = {}
	header = None
	seq_chunks = []
	with open(fasta_file, "r") as f:
		for line in f:
			if not line.strip():
				continue
			if line.startswith(">"):
				if header is not None:
					content[header] = "".join(seq_chunks)
				header = line[1:].strip().split()[0]  # trim after first whitespace
				seq_chunks = []
			else:
				seq_chunks.append(line.strip())
		if header is not None:
			content[header] = "".join(seq_chunks)
	return content

def translate(seq, genetic_code, unknown="X", internal_stop_to_x=False, keep_terminal_stop=True):
	"""Translate DNA to protein.
	   - Unknown/ambiguous codons -> 'X'
	   - Stop codons per table -> '*'
	   - Optionally convert internal '*' to 'X' (keeps final '*' if present).
	"""
	seq = seq.upper().replace("U", "T")
	pep = []
	n = len(seq) // 3
	for i in range(n):
		codon = seq[i*3:i*3+3]
		aa = genetic_code.get(codon)
		if aa is None:
			aa = unknown
		pep.append(aa)
	pep = "".join(pep)

	if internal_stop_to_x and "*" in pep:
		if keep_terminal_stop and pep.endswith("*"):
			pep = pep[:-1].replace("*", "X") + "*"
		else:
			pep = pep.replace("*", "X")
	return pep

def translate_file(logger, in_fa, out_fa, genetic_code, internal_stop_to_x=False):
	seqs = load_multiple_fasta_file(in_fa)
	internal_stop_count = 0
	with open(out_fa, "w") as out:
		for header, nt in seqs.items():
			pep = translate(nt, genetic_code, unknown="X",
							internal_stop_to_x=internal_stop_to_x,
							keep_terminal_stop=True)
			if "*" in pep[:-1]:
				internal_stop_count += 1
			out.write(f">{header}\n{pep}\n")
	logger.info(f" {os.path.basename(in_fa)} -> {os.path.basename(out_fa)} | "
		  f"seqs: {len(seqs)} | internal-stops (pre-fix): {internal_stop_count}")

def gather_inputs(spec):
	"""Return list of input files. If directory, grab common FASTA extensions."""
	if os.path.isdir(spec):
		exts = ("*.fa", "*.fasta", "*.fna", "*.fas")
		files = []
		for ext in exts:
			files.extend(glob.glob(os.path.join(spec, ext)))
		files.sort()
		return files
	else:
		return [spec]

def make_output_path(out_spec, in_file):
	"""Decide output path: if out_spec is a dir or endswith '/', write there; else treat as file path."""
	if out_spec.endswith(os.sep) or os.path.isdir(out_spec):
		os.makedirs(out_spec, exist_ok=True)
		base = os.path.basename(in_file)
		root = os.path.splitext(base)[0]
		return os.path.join(out_spec, f"{root}.pep.fa")
	else:
		# Single input -> exact output file path
		# Multi-input with a single-file out_spec is not supported
		return out_spec

#functions to check for presence of working BUSCO installation in the user's system
def check_native_busco(user_path=None):
	# distinguish between a bare command name and an explicit file path
	if user_path and os.path.isabs(user_path):  # only treat as file path if absolute
		if os.path.isfile(user_path) and os.access(user_path, os.X_OK):
			try:
				subprocess.run([user_path, "--version"], check=True,
							   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
				return user_path
			except Exception:
				return None
		else:
			return None

	# bare command name — search PATH
	busco_path = shutil.which(user_path or "busco")
	if busco_path:
		try:
			subprocess.run([busco_path, "--version"], check=True,
						   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
			return busco_path
		except Exception:
			return None

	return None


def check_docker_busco():
	# CHeck docker installation
	if shutil.which("docker") is None:
		return None

	# Check for presence of 'busco' in the name of the local image?
	try:
		images = subprocess.check_output(["docker", "images", "--format", "{{.Repository}}:{{.Tag}}"]).decode().splitlines()
		busco_images = [img for img in images if "busco" in img.lower()]
		if not busco_images:
			return None
	except Exception:
		return None

	# Try running 'busco --version' inside one of them
	for image in busco_images:
		try:
			subprocess.run(
				["docker", "run", "--rm", image, "busco", "--version"],
				check=True,
				stdout=subprocess.DEVNULL,
				stderr=subprocess.DEVNULL
			)
			return image  # BUSCO is usable via this Docker image
		except Exception:
			continue

	return None  # BUSCO not available via Docker

def detect_busco(busco_user_path):

	if busco_user_path == 'busco_docker':
		docker_busco = check_docker_busco()
		if docker_busco:
			return docker_busco

	else:
		#try native busco
		native = check_native_busco(busco_user_path)
		if native:
			return native

	return None

#function to predownload BUSCO datasets
def pre_download_databases(logger, busco_path, org_type, busco_dir, busco_db_dir):
	"""Pre-download BUSCO databases to avoid concurrent download issues"""
	logger.info("Pre-downloading BUSCO databases...")

	# Ensure the database directory exists
	os.makedirs(busco_db_dir, exist_ok=True)

	dummy_fasta = os.path.join(busco_dir, "dummy.faa")
	with open(dummy_fasta, 'w') as f:
		f.write(">dummy\nMVKIILFVGLLFSSVTYGC\n")

	if org_type == 'eukaryote':
		cmd = f"{busco_path} -i {dummy_fasta} -m proteins --auto-lineage-euk -q -o dummy_download --out_path {busco_dir} --download_path {busco_db_dir}"
	elif org_type == 'prokaryote':
		cmd = f"{busco_path} -i {dummy_fasta} -m proteins --auto-lineage-prok -q -o dummy_download --out_path {busco_dir} --download_path {busco_db_dir}"

	logger.info(f"Downloading BUSCO databases to: {busco_db_dir}")
	# Redirect stderr to devnull to suppress busco error messages for the above dummy dataset
	with open(os.devnull, 'w') as devnull:
		p = subprocess.Popen(args=cmd, shell=True, stderr=devnull, stdout=devnull)
		p.communicate()

	# Cleanup
	os.remove(dummy_fasta)
	dummy_dir = os.path.join(busco_dir, "dummy_download")
	if os.path.exists(dummy_dir):
		shutil.rmtree(dummy_dir)

	logger.info("Database pre-download completed.")

#function to clean gene names for BUSCO analysis
def clean_headers_for_busco(gene_id):
	for char, replacement in replacements.items():
		safe_gene_id = gene_id.replace(char, replacement)
	return safe_gene_id

#function to run BUSCO
def run_busco(logger, orgname, fasta_file, busco_path, busco_dir, busco_dir_final, busco_threads_per_organism, org_type,host_cache_dir, busco_db_dir, buscolineage):
	#code to create a temp CDS fasta file to clean headers and use this clean protein file for busco analysis - this temp file will be deleted after the busco run completes
	tmp_fasta = os.path.join(busco_dir, f"{orgname}_temp.cds.fasta")
	try:
		# Process and write headers + sequences
		with open(fasta_file, "r") as infile, open(tmp_fasta, "w") as outfile:
			for line in infile:
				if line.startswith(">"):
					new_header = clean_headers_for_busco(line.strip()[1:])  # remove ">"
					outfile.write(f">{new_header}\n")
				else:
					outfile.write(line)
		if org_type == 'eukaryote':
			if buscolineage=='auto':
				cmd = busco_path + ' -i ' + tmp_fasta + ' -m proteins --auto-lineage-euk -q -o ' + orgname + ' -c ' + busco_threads_per_organism + ' --out_path ' + busco_dir + ' --download_path ' + busco_db_dir
			else:
				lineage=buscolineage
				cmd = busco_path + ' -i ' + tmp_fasta + ' -m proteins -l '+ lineage + ' -q -o ' + orgname + ' -c ' + busco_threads_per_organism + ' --out_path ' + busco_dir + ' --download_path ' + busco_db_dir
		elif org_type == 'prokaryote':
			if buscolineage=='auto':
				cmd = busco_path + ' -i ' + tmp_fasta + ' -m proteins --auto-lineage-prok -q -o ' + orgname + ' -c ' + busco_threads_per_organism + ' --out_path ' + busco_dir + ' --download_path ' + busco_db_dir
			else:
				lineage=buscolineage
				cmd = busco_path + ' -i ' + tmp_fasta + ' -m proteins -l ' + lineage + ' -q -o ' + orgname + ' -c ' + busco_threads_per_organism + ' --out_path ' + busco_dir + ' --download_path ' + busco_db_dir
		p = subprocess.Popen(args=cmd, shell=True)
		p.communicate()
	finally:
		#remove the temporary pep fasta file
		if os.path.exists(tmp_fasta):
			os.remove(tmp_fasta)
	species_dir = Path(busco_dir) / orgname
	# Get all run_* folders (ignore auto_lineage and subdirs of auto_lineage)
	run_dirs = [
		p for p in species_dir.iterdir()
		if p.is_dir() and p.name.startswith("run_") and "auto_lineage" not in p.name
	]

	# Decide which run folder to use
	feff_result = str(orgname) + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA'
	if not run_dirs:
		logger.error(f"BUSCO error for {orgname}")
		selected_dir = ''

	elif len(run_dirs) == 1:
		selected_dir = run_dirs[0]

	else:
		# Prefer a folder that does NOT include "eukaryota"
		if org_type == 'eukaryote':
			filtered = [d for d in run_dirs if "eukaryota" not in d.name.lower()]
			if not filtered:
				logger.info(f"Multiple run_* folders found for {orgname}, but all are eukaryota-based.")
				selected_dir = run_dirs[0] # use the first eukaryota run folder
			else:
				selected_dir = filtered[0]  # Use the first non-eukaryota run folder
		elif org_type == 'prokaryote':
			filtered = [d for d in run_dirs if "prokaryota" not in d.name.lower()]
			if not filtered:
				logger.info(f"Multiple run_* folders found for {orgname}, but all are prokaryota-based..")
				selected_dir = run_dirs[0] # use the first prokaryota run folder
			else:
				selected_dir = filtered[0]  # Use the first non-prokaryota run folder

	#Find full_table.tsv inside selected folder (ignore subfolders of auto_lineage if any)
	if selected_dir:
		full_tables = list(selected_dir.rglob("full_table.tsv"))
	else:
		full_tables = []

	if not full_tables:
		logger.info(f"No 'full_table.tsv' found for {orgname}")

	else:
		# Get the path to the full_table.tsv
		out_path = full_tables[0]
		full_table = str(out_path)
		complete_genes = []
		all_entries = []
		busco_list = []
		with open(full_table, 'r') as f:
			for line in f:
				if not line.startswith('#'):
					columns = line.strip().split('\t')
					if columns and  len(columns) >= 3:  # Ensure the line is not empty
						busco_list.append(columns[0])  # First column contains the BUSCO ID
						busco_id = columns[0]
						status = columns[1]
						gene_id = columns[2]

						all_entries.append({
							'busco_id': busco_id,
							'status': status,
							'gene_id': gene_id
						})
		# Count frequencies
		busco_id_counts = Counter([entry['busco_id'] for entry in all_entries])
		status_counts = Counter([entry['status'] for entry in all_entries])

		# Extract single-copy genes (Strategy 2)
		complete_from_single = []
		complete_from_multi = []

		for entry in all_entries:
			if entry['status'] == 'Complete' and entry['gene_id']:
				if busco_id_counts[entry['busco_id']] == 1:
					complete_from_single.append(entry['gene_id'])
				else:
					complete_from_multi.append(entry['gene_id'])

		single_copy_busco_genes = complete_from_single

		gene_counts = Counter(busco_list)  # Count occurrences of each gene
		frequency_counts = Counter(gene_counts.values())  # Count occurrences of each frequency

		x_values = sorted(frequency_counts.keys())  # Unique occurrence frequencies
		y_values = [frequency_counts[x] for x in x_values]  # Count of genes occurring those many times

		#code to extract C, D percentages from busco log file
		busco_log_path = next(Path(species_dir).rglob("busco.log"))
		busco_log = str(busco_log_path)
		with open (busco_log,'r') as f:
			contents = f.read()
			# Find all matches of the C/S/D/F/M pattern to take the C, D values of the specific BUSCO lineage instead of the general one that is presented first
			matches = re.findall(r"C:(\d+(?:\.\d+)?)%\[S:\d+(?:\.\d+)?%,D:(\d+(?:\.\d+)?)%", contents)
		if not matches:
			logger.warning(f"No C/D values found in BUSCO log. BUSCO run might be incomplete for {orgname}")
		else:
			# If multiple matches, use the last one (usually the lineage-specific one)
			C_value = float(matches[-1][0])
			D_value = float(matches[-1][1])
			plt.bar(x_values, y_values, color='blue', alpha=0.7)
			plt.xlabel("Occurrence Frequency")
			plt.ylabel("Number of BUSCO Genes")
			plt.title("Gene Occurrence Frequency Distribution")
			plt.xticks(x_values)  # Ensure all x-axis values are labeled
			plt.grid(axis='y', linestyle='--', alpha=0.7)
			plot_file = os.path.join(busco_dir_final, orgname+'_gene_frequency_distribution.png')
			plt.savefig(plot_file)
			plt.close()
			total_genes = sum(frequency_counts.values())  # Total number of BUSCO genes
			Feff_numerator = sum(freq * count for freq, count in frequency_counts.items())  # Σ (C_i * N_i)
			Feff = Feff_numerator / total_genes if total_genes > 0 else 1  # Avoid division by zero
			feff_result = str(orgname)+'\t'+str(Feff)+'\t'+str(C_value)+'\t'+str(D_value)
	return feff_result

def load_IDs(filename):
	IDs = []
	if not filename.lower().endswith('.gz'):
		with open(filename, "r") as f:
			line = f.readline()
			while line:
				ID = line.strip()
				if len(ID) > 3:
					if "\t" in line:
						tmp = line.strip().split('\t')
						for each in tmp:
							if len(each) > 3:
								IDs.append(each)
					else:
						IDs.append(ID)
				line = f.readline()
	else:
		with gzip.open(filename, "rt") as f:
			line = f.readline()
			while line:
				ID = line.strip()
				if len(ID) > 3:
					if "\t" in line:
						tmp = line.strip().split('\t')
						for each in tmp:
							if len(each) > 3:
								IDs.append(each)
					else:
						IDs.append(ID)
				line = f.readline()
	return IDs

#functions to run kallisto
def get_data_for_jobs_to_run(readfile_status, logger, read_file_folders, final_output_folder, index_file, tmp_cluster_folder):
	"""! @brief collect all infos to run jobs """

	jobs_to_do = []
	for folder in read_file_folders:
		ID = folder.split('/')[-1]
		status = True

		# --- get read file --- #
		PE_status = True
		SRA = False
		if readfile_status is not None:
			read_file1 = folder + "/" + ID + "_R1_001.fastq.gz"
		else:
			read_file1 = folder + "/" + ID + "_R1_001.fastq"
		if not os.path.isfile(read_file1):
			# print "ERROR: file missing - " + read_file1
			PE_status = False
			if not PE_status:
				if readfile_status is not None:
					read_file1 = folder + "/" + ID + "_pass_1.fastq.gz"
				else:
					read_file1 = folder + "/" + ID + "_pass_1.fastq"
				if os.path.isfile(read_file1):
					PE_status = True
					SRA = True
					if readfile_status is not None:
						read_file2 = folder + "/" + ID + "_pass_2.fastq.gz"
					else:
						read_file2 = folder + "/" + ID + "_pass_2.fastq"
					if not os.path.isfile(read_file2):
						# print "ERROR: file missing - " + read_file2
						PE_status = False
				else:
					if readfile_status is not None:
						read_file1 = folder + "/" + ID + "_1.fastq.gz"
					else:
						read_file1 = folder + "/" + ID + "_1.fastq"
					if os.path.isfile(read_file1):
						PE_status = True
						SRA = True
						if readfile_status is not None:
							read_file2 = folder + "/" + ID + "_2.fastq.gz"
						else:
							read_file2 = folder + "/" + ID + "_2.fastq"
						if not os.path.isfile(read_file2):
							# print "ERROR: file missing - " + read_file2
							PE_status = False
					else:
						if readfile_status is not None:
							read_file1 = folder + "/" + ID + "_R1.fq.gz"
						else:
							read_file1 = folder + "/" + ID + "_R1.fq"
						if os.path.isfile(read_file1):
							PE_status = True
							SRA = True
							if readfile_status is not None:
								read_file2 = folder + "/" + ID + "_R2.fq.gz"
							else:
								read_file2 = folder + "/" + ID + "_R2.fq"
							if not os.path.isfile(read_file2):
								# print "ERROR: file missing - " + read_file2
								PE_status = False
		if not SRA:
			if readfile_status is not None:
				read_file2 = folder + "/" + ID + "_R2_001.fastq.gz"
			else:
				read_file2 = folder + "/" + ID + "_R2_001.fastq"
			if not os.path.isfile(read_file2):
				# print "ERROR: file missing - " + read_file2
				PE_status = False
		if not PE_status:
			if readfile_status is not None:
				read_file1 = folder + "/" + ID + "_R1_001.fastq.gz"
			else:
				read_file1 = folder + "/" + ID + "_R1_001.fastq"
			if not os.path.isfile(read_file1):
				if readfile_status is not None:
					read_file1 = folder + "/" + ID + ".fastq.gz"
				else:
					read_file1 = folder + "/" + ID + ".fastq"
				if not os.path.isfile(read_file1):
					if readfile_status is not None:
						read_file1 = folder + "/" + ID + "_1.fastq.gz"
					else:
						read_file1 = folder + "/" + ID + "_1.fastq"
					if not os.path.isfile(read_file1):
						if readfile_status is not None:
							read_file1 = folder + "/" + ID + "_R1.fq.gz"
						else:
							read_file1 = folder + "/" + ID + "_R1.fq"
						if not os.path.isfile(read_file1):
							if readfile_status is not None:
								read_file1 = folder + "/" + ID + "_pass.fastq.gz"
							else:
								read_file1 = folder + "/" + ID + "_pass.fastq"
							if not os.path.isfile(read_file1):
								status = False
			read_file2 = False

		# --- get reference for quantification --- #
		output_dir = os.path.join(tmp_cluster_folder, f'{ID}/')
		if not os.path.exists(output_dir):
			os.makedirs(output_dir)
		tmp_result_file = output_dir + "abundance.tsv"
		try:
			timestr = time.strftime("%Y_%m_%d_")
		except ModuleNotFoundError:
			timestr = ""
		final_result_file = final_output_folder + timestr + ID + ".tsv"
		if os.path.isfile(final_result_file):
			status = False
		if status:
			jobs_to_do.append(
				{'r1': read_file1, 'r2': read_file2, 'out': output_dir, 'index': index_file, 'tmp': tmp_result_file,
				 'fin': final_result_file, "ID": ID})
	return jobs_to_do

def job_executer(logger, jobs_to_run, kallisto, threads):
	"""! @brief run all jobs in list """

	for idx, job in enumerate(jobs_to_run):
		logger.info("running job " + str(idx + 1) + "/" + str(len(jobs_to_run)) + " - " + job["ID"] + "\n")

		if job['r2']:
			cmd2 = " ".join([kallisto, "quant", "--index=" + job['index'], "--output-dir=" + job['out'],
							 "--threads " + str(threads), job['r1'], job['r2']])
		else:
			cmd2 = " ".join(
				[kallisto, "quant", "--index=" + job['index'], "--single -l 200 -s 100", "--output-dir=" + job['out'],
				 "--threads " + str(threads), job['r1']])
		p = subprocess.Popen(args=cmd2, shell=True)
		p.communicate()

		p = subprocess.Popen(args="cp " + job["tmp"] + " " + job["fin"], shell=True)
		p.communicate()

#functions to merge the TPM. Counts files per batch
def load_counttable(counttable):
	"""! @brief load data from counttable """

	counts = {}
	tpms = {}
	with open(counttable, "r") as f:
		f.readline()  # remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			counts.update({parts[0]: float(parts[3])})
			tpms.update({parts[0]: float(parts[4])})
			line = f.readline()
	return counts, tpms


def generate_mapping_table(logger, gff_file):
	"""! @brief generate transcript to gene mapping table """

	transcript2gene = {}
	with open(gff_file, "r") as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] in ["mRNA", "transcript"]:
					try:
						ID = parts[-1].split(';')[0].split('=')[1]
						parent = parts[-1].split('arent=')[1]
						if ";" in parent:
							parent = parent.split(';')[0]
						transcript2gene.update({ID: parent})
					except:
						logger.error(line)
			line = f.readline()
	return transcript2gene


def map_counts_to_genes(logger, transcript2gene, counts):
	"""! @brief map transcript counts to parent genes """

	error_collector = []
	gene_counts = {}
	for key in counts.keys():
		try:
			gene_counts[transcript2gene[key]] += counts[key]
		except KeyError:
			try:
				gene_counts.update({transcript2gene[key]: counts[key]})
			except KeyError:
				error_collector.append(key)
				gene_counts.update({key: counts[key]})
	if len(error_collector) > 0:
		logger.info("number of unmapped transcripts: " + str(len(error_collector)) + "\n")
	return gene_counts


def generate_output_file(output_file, data):
	"""! @brief generate output file for given data dictionary """
	try:
		timestr = time.strftime("%Y_%m_%d")
	except ModuleNotFoundError:
		timestr = ""
	samples = list(sorted(list(data.keys())))

	with open(output_file, "w") as out:
		out.write("\t".join([timestr] + samples) + '\n')
		for gene in list(sorted(list(data.values())[0].keys())):
			new_line = [gene]
			for sample in samples:
				new_line.append(data[sample][gene])
			out.write("\t".join(map(str, new_line)) + '\n')

#function to read tpm/ counts TXT files for combining them
def read_tpm_file(filename):
	"""Read TPM file into a dictionary with gene names as keys"""
	with open(filename, 'r') as f:
		lines = f.readlines()

	# First row: date + SRA accessionso
	header = lines[0].rstrip('\n').split('\t')
	date_info = header[0]
	sra_accessions = header[1:]  # Skip first column (date)

	# Rest: gene_name + expression values
	gene_data = OrderedDict()
	for line in lines[1:]:
		cols = line.rstrip('\n').split('\t')
		gene_name = cols[0]
		expression_values = cols[1:]
		gene_data[gene_name] = expression_values

	return date_info, sra_accessions, gene_data

#functions to load and read expression files for the RNA-seq sample removal QC step
def load_all_TPMs(exp_file):
	"""! @brief load all values from given TPM file """

	data = {}
	genes = []
	with open(exp_file, "r") as f:
		headers = f.readline().strip()
		if "\t" in headers:
			headers = headers.split('\t')
		else:
			headers = [headers]
		if headers[0] == "gene":
			headers = headers[1:]
		elif headers[0] == exp_file.split('/')[-1][:10]:
			headers = headers[1:]
		for header in headers:
			data.update({header: []})
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			genes.append(parts[0])
			for idx, val in enumerate(parts[1:]):
				data[headers[idx]].append(float(val))
			line = f.readline()
	return data, genes

def load_black_IDs(black_list_file):
	"""! @brief load IDs from given black list """

	black_list = {}
	with open(black_list_file, "r") as f:
		line = f.readline()
		while line:
			if len(line) > 3:
				black_list.update({line.strip(): None})
			line = f.readline()
	return black_list

#functions to remove isoforms
def load_fasta(fasta_file):
	"""! @brief load FASTA alignment into dictionary	"""

	sequences = {}
	with open(fasta_file) as f:
		header = f.readline()[1:].strip()
		if " " in header:
			header = header.split(' ')[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
				sequences.update({header: "".join(seq)})
				header = line.strip()[1:]
				if " " in header:
					header = header.split(' ')[0]
				seq = []
			else:
				seq.append(line.strip())
			line = f.readline()
		sequences.update({header: "".join(seq)})
	return sequences


def load_hits_per_bait(blast_result_file, scorecut, simcut, lencut):
	"""! @brief load BLAST hits per bait
	@return dictionary with baits as keys and lists of hits as values
	"""

	hits = {}
	with open(blast_result_file, "r") as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if parts[0] != parts[1]:
				if float(parts[-1]) > scorecut:  # BLAST hits are filtered with user defined criteria
					if float(parts[2]) > simcut:
						if int(parts[3]) > lencut:
							try:
								if parts[-1] not in hits[parts[0]]:
									hits[parts[0]].append(parts[1])  # add to existing dictionary entry
							except KeyError:
								hits.update({parts[0]: [parts[1]]})  # generate new entry in dictionary
			line = f.readline()
	return hits


def load_aln_similarity(alignment, snv_cutoff):
	"""! @brief load alignment similarity """

	results = []
	seqIDs = list(alignment.keys())
	for idx1, ID1 in enumerate(seqIDs):
		for idx2, ID2 in enumerate(seqIDs[1:]):
			if idx2 >= idx1:
				aln_seq1 = alignment[ID1]
				aln_seq2 = alignment[ID2]
				snv_counter = 0
				indel_counter = 0
				for i1, nt in enumerate(aln_seq1):
					if nt != "-" and aln_seq2[i1] != "-":
						if nt != aln_seq2[i1]:
							snv_counter += 1
					else:
						indel_counter += 1
				if snv_counter <= snv_cutoff:
					results.append({'ID1': ID1, 'ID2': ID2, 'snvs': snv_counter, 'indels': indel_counter})
	return results


def identify_isoforms(alignment_results, snv_cutoff):
	"""! @brief identify isoforms """

	isoform_groups = []
	for group in alignment_results:
		tmp_groups = []
		for entry in group:
			if entry['snvs'] <= snv_cutoff:
				status = False  # genes not included yet
				for idx, g in enumerate(tmp_groups):
					if entry['ID1'] in g:
						tmp_groups[idx].append(entry['ID2'])
						status = True
						break
					elif entry['ID2'] in g:
						tmp_groups[idx].append(entry['ID1'])
						status = True
						break
				if not status:
					tmp_groups.append([entry['ID1'], entry['ID2']])
		for g in tmp_groups:  # go through all tmp groups to break up too nested structure
			if len(g) > 0:  # only collect non empty lists
				isoform_groups.append(g)
	return isoform_groups


def identify_repr_isoform_per_group(logger, isoforms, seqs):
	"""! @brief identify representative isoform per group """

	repr_seq_collection = {}
	repr_id_collection = {}
	for group in isoforms:
		#debug lines
		if any('TC00001' in id for id in group):
			print(f'isoform group is {group}')
		tmp_seqs = []
		for ID in group:
			tmp_seqs.append({'ID': ID, 'seq': seqs[ID], 'len': len(seqs[ID])})
		sorted_list = sorted(tmp_seqs, key=itemgetter('len'))
		repr_seq_collection.update({sorted_list[-1]['ID']: sorted_list[-1]['seq']})
		ID_values=[]
		for element in sorted_list[:-1]:
			ID_values.append(element['ID'])
		repr_id_collection.update({sorted_list[-1]['ID']: ID_values})#make a dictionary where key is the repr isoform and the value is a list of all isoforms in the group represented by this isoform
	return repr_seq_collection, repr_id_collection

#function to produce TPM file without isoforms
def keep_primary_transcript_exp(repr_ids, repr_tpm_file, repr_counts_file, primary_transcript_cds_file, exp_file):
	primary_transcripts = []
	with open (primary_transcript_cds_file, 'r') as f:
		for line in f:
			if line.startswith(">"):
				primary_transcripts.append(line[1:].strip())

	df = pd.read_csv(exp_file, sep='\t', index_col='gene')
	original_transcript_order = df.index.tolist()#obtaining the transcript order from the original file

	result_rows = []
	all_mapped_transcripts = set()#to collect transcripts from the isoform representative collection
	for retained_transcript, purged_transcripts in repr_ids.items():#collect transcripts to be purged for summing up their expression values onto the retained primary transcript
		transcripts_to_sum = [retained_transcript]+purged_transcripts
		all_mapped_transcripts.update(transcripts_to_sum)

		summed_up_row = df.loc[transcripts_to_sum].sum()
		summed_up_row.name = retained_transcript
		result_rows.append(summed_up_row)

	df_summed = pd.DataFrame(result_rows)
	df_summed.index.name = 'gene'

	df_unmapped = df[~df.index.isin(all_mapped_transcripts)]#dataframe with transcripts not found in the repr_ids dic mapping

	df_total = pd.concat([df_summed, df_unmapped])
	df_total = df_total.reset_index()
	#preserving original transcript order in the isoform purged file
	retained_order = [t for t in original_transcript_order if t in df_total['gene'].values]
	df_total = df_total.set_index('gene').loc[retained_order].reset_index()

	if '.tpm' in exp_file:
		df_total.to_csv(repr_tpm_file, sep="\t", index=False)
	elif '.count' in exp_file:
		df_total.to_csv(repr_counts_file, sep="\t", index=False)

	return repr_tpm_file

#function to fetch SRA files
def fetch_worker (attempts,sradir, accession,prefetch_command, minimum_sra_file_size_threshold, fasterq_dump, cores, logger, completed_accessions_file, base_wait, failed_accessions_file):
	for attempt in range(attempts):
		# Create subfolder for this accession
		acc_dir = os.path.join(sradir, accession)
		if os.path.exists(acc_dir):
			shutil.rmtree(acc_dir)
		os.makedirs(acc_dir, exist_ok=True)
		prefetched_file = os.path.join(acc_dir, f"{accession}.sra")
		try:
			cmd = f"{prefetch_command} --max-size 200G {accession} -O {sradir}"
			prefetch_result = subprocess.run(cmd, shell=True)
			if prefetch_result.returncode != 0:  # if prefetch fails due to issues lik network disruption
				raise RuntimeError(f"prefetch for {accession} failed with code {prefetch_result.returncode}")
			if not os.path.exists(prefetched_file):  # if prefetch did not fetch an SRA file in the first place
				raise RuntimeError(f"prefetch for {accession} did not fetch a .sra file")
			filesize = (os.path.getsize(prefetched_file))/(1024 ** 2)#converting the file size returned in bytes by getsize to Mb
			if filesize < minimum_sra_file_size_threshold:  # in case prefetch gets the sralite files
				raise RuntimeError(
					f"File size of the prefetched file for {accession} is smaller than the threshold size of {minimum_sra_file_size_threshold}")

			# fasterq-dump
			cmd = f"{fasterq_dump} --split-3 --outdir {acc_dir} --skip-technical --threads {cores} {prefetched_file}"
			fasterq_dump_result = subprocess.run(cmd, shell=True)
			if fasterq_dump_result.returncode != 0:
				raise RuntimeError(
					f"fasterq-dump failed for {accession} with error code {fasterq_dump_result.returncode}")

			# pigz (generalized)
			fq_patterns = [
				os.path.join(acc_dir, f"{accession}*.fastq"),
				os.path.join(acc_dir, f"{accession}*.fq"),
			]
			fq_files = []
			for pattern in fq_patterns:
				fq_files.extend(glob.glob(pattern))
			pigz_result = None
			if fq_files:
				cmd = f"pigz -p {cores} " + " ".join(fq_files)
				pigz_result = subprocess.run(cmd, shell=True)
				if pigz_result.returncode != 0:
					raise RuntimeError(f"pigz failed with code {pigz_result.returncode}")
				# After pigz, check what was created
				logger.info(f"Files in {acc_dir}:")
				all_valid = False
				for f in os.listdir(acc_dir):
					logger.info(f"  {f}")
			gz_files = glob.glob(os.path.join(acc_dir, "*.gz"))
			all_valid = gz_files and all(os.path.getsize(f) > 0 for f in gz_files)  # checks for all fetched gzipped file sizes and if all are non-empty
			if not all_valid:
				raise RuntimeError(f"gz files missing or empty for {accession}")
			# if fetching is successful mark as completed accession and break out of the retry loop
			with open(completed_accessions_file, 'a') as out:
				out.write(f'{accession}\n')
				out.flush()
				os.fsync(out.fileno())
			successfull_accession = accession
			fetch_queue.put(successfull_accession)
			break
		except RuntimeError as e:
			logger.warning(f"Attempt {attempt + 1} failed for {accession}: {e}")
			# clean up partial files before retrying
			if os.path.exists(acc_dir):
				shutil.rmtree(acc_dir)
				os.makedirs(acc_dir, exist_ok=True)
			if attempt < attempts - 1:
				wait = base_wait * (2 ** attempt)
				logger.info(f"Retrying {accession} in {wait}s")
				time.sleep(wait)
			else:
				logger.error(f"All attempts exhausted for {accession}, marking as failed")
				with open(failed_accessions_file, 'a') as out:
					out.write(f'{accession}\n')
					out.flush()
					os.fsync(out.fileno())

#function to perform kallisto quantification as soon as SRA file is fetched
def kallisto_worker(index_file, sradir,readfile_status, logger,kallistodir,kallisto,cds_file,cores,tmpdir):
	while True:
		accession = fetch_queue.get()
		if accession is barrier:
			break
		# Run Kallisto
		# --- load data --- #
		acc_dir = os.path.join(sradir, accession)
		single_read_file_folders = [acc_dir]
		logger.info("Number of FASTQ file folders detected: " + str(len(single_read_file_folders)) + "\n")

		# --- prepare jobs to run --- #
		jobs_to_run = get_data_for_jobs_to_run(readfile_status, logger, single_read_file_folders, kallistodir, index_file,tmpdir)
		logger.info("Number of jobs to run: " + str(len(jobs_to_run)) + "\n")

		# --- run jobs --- #
		job_executer(logger, jobs_to_run, kallisto, cores)
		fetch_queue.task_done()

#function to check if tsv files for merge have the same genes in the same order in the first column
def check_gene_consistency(file_paths):
	gene_lists = {}
	for f in file_paths:
		genes = pd.read_csv(f, sep='\t', usecols=['gene'])['gene'].tolist()
		gene_lists[f] = genes

	base_file = file_paths[0]
	base_genes = gene_lists[base_file]

	errors = []
	for f in file_paths[1:]:
		other_genes = gene_lists[f]

		if base_genes == other_genes:
			continue  # identical names and order, no issue

		# collect errors
		if set(base_genes) == set(other_genes):
			# Same genes but wrong order — find first mismatch position
			mismatches = [(i, base_genes[i], other_genes[i])
						  for i in range(len(base_genes))
						  if base_genes[i] != other_genes[i]]
			errors.append(
				f"\nFile: {f}"
				f"\nGene names match but order differs."
				f"\nFirst mismatch at row {mismatches[0][0]+1}: "
				f"Expected '{mismatches[0][1]}', found '{mismatches[0][2]}'"
				f" ({len(mismatches)} mismatched positions total)"
			)
		else:
			# Different gene sets entirely
			only_in_base = set(base_genes) - set(other_genes)
			only_in_other = set(other_genes) - set(base_genes)
			errors.append(
				f"\nFile: {f}"
				+ (f"\nGenes only in {base_file}: {only_in_base}" if only_in_base else "")
				+ (f"\nGenes only in {f}: {only_in_other}" if only_in_other else "")
			)

	if errors:
		return False, "Gene mismatch detected:" + "".join(errors)
	return True, None

def merge_expression_tsvs(base_file, additional_files):
	all_files = [base_file] + additional_files

	# --- Validate before merging ---
	is_consistent, error_msg = check_gene_consistency(all_files)
	if not is_consistent:
		raise ValueError(error_msg)

	dfs = [pd.read_csv(f, sep='\t', index_col='gene') for f in all_files]
	merged = reduce(lambda left, right: left.join(right, how='inner'), dfs)
	return merged.reset_index()

def main(arguments):
	if '--sample_name' in arguments:
		orgname = arguments[arguments.index('--sample_name')+1]
	else:
		orgname = 'sample'
	if  '--sra' in arguments:#list of SRA accessions to be fetched from NCBI SRA
		todo_sras = arguments[arguments.index('--sra')+1]
	elif '--readfiles' in arguments:#Full path to folder with RNA-seq read files
		todo_sras = arguments[arguments.index('--readfiles')+1]
	if '--fastq_pattern' in arguments:
		pattern_names = arguments[arguments.index('--fastq_pattern')+1]#specify fastq read file name pattern for the paired end files separated by commas without spaces like - _pass_1,_pass_2_
		pattern_names_list = pattern_names.split(',')
	else:
		pattern_names_list = ["_pass_1", "_pass_2"]

	cdsfile = arguments[arguments.index('--cds')+1]#full path to CDS file

	if '--annotation_qc' in arguments:# yes or no for BUSCO-based QC of the PEP file
		qc = arguments[arguments.index('--annotation_qc')+1]
	else:
		qc = 'yes'

	# Option for user to give full path to busco
	if '--busco' in arguments:
		busco_path = arguments[arguments.index('--busco') + 1]  # busco_docker or default busco
	else:
		busco_path = 'busco'

	# Option for user to specify BUSCO lineage per input organism
	if '--busco_lineage' in arguments:
		buscolineage = arguments[arguments.index('--busco_lineage') + 1]
	else:
		buscolineage = 'auto'

	if '--busco_version' in arguments:
		busco_version = arguments[arguments.index('--busco_version') + 1]
	else:
		busco_version = "v6.0.0"  # the most recent version of BUSCO at the time of writing this script

	if '--container_version' in arguments:
		container_version = arguments[arguments.index('--container_version') + 1]
	else:
		container_version = 'cv1'

	if '--docker_host_path' in arguments:
		host_path = arguments[arguments.index('--docker_host_path') + 1]
	else:
		host_path = "NA"

	if '--docker_container_path' in arguments:
		container_path = arguments[arguments.index('--docker_container_path') + 1]
	else:
		container_path = "NA"

	if '--organism_type' in arguments:
		org_type = arguments[arguments.index('--organism_type')+1]
	else:
		org_type = 'eukaryote'

	if '--kallisto' in arguments:#full path to kallisto
		kallisto = arguments[arguments.index('--kallisto')+1]
	else:
		kallisto = 'kallisto'
	if '--fasterq-dump' in arguments:#full path to fasterq-dump from NCBI SRA toolkit
		fasterq_dump = arguments[arguments.index('--fasterq-dump')+1]
	else:
		fasterq_dump = 'fasterq-dump'

	if '--prefetch' in arguments:
		prefetch_command = arguments[arguments.index('--prefetch')+1]
	else:
		prefetch_command = "prefetch"
	cds_file = arguments[arguments.index('--cds')+1]
	if '--threads' in arguments:
		cores = int(arguments[arguments.index('--threads')+1])
	else:
		cores = 4
	if '--batch_size' in arguments:
		batch_size=int(arguments[arguments.index('--batch_size')+1])
	else:
		batch_size = 1

	if '--remove_isoforms' in arguments:
		remove_isoforms = arguments[arguments.index('--remove_isoforms')+1]
	else:
		remove_isoforms = 'yes'
	#flags to remove SRA samples after qc check based on the contribution of top 100 genes to the gene expression matrix
	if '--min' in arguments:
		min_cutoff = int(arguments[arguments.index('--min') + 1])
	else:
		min_cutoff = 10  # value in percent
	if '--max' in arguments:
		max_cutoff = int(arguments[arguments.index('--max') + 1])
	else:
		max_cutoff = 80

	if '--mincounts' in arguments:
		min_counts = int(arguments[arguments.index('--mincounts') + 1])
	else:
		min_counts = 1000000

	if '--black' in arguments:
		black_list_file = arguments[arguments.index('--black') + 1]
		black_list = load_black_IDs(black_list_file)
	else:
		black_list = {}

	#flags to remove isoforms
	if '--scorecut' in arguments:
		scorecut = int(arguments[arguments.index('--scorecut') + 1])
	else:
		scorecut = 100

	if '--simcut' in arguments:
		simcut = float(arguments[arguments.index('--simcut') + 1])
	else:
		simcut = 99.0

	if '--lencut' in arguments:
		lencut = int(arguments[arguments.index('--lencut') + 1])
	else:
		lencut = 100

	if '--snvcut' in arguments:
		snv_cutoff = int(arguments[arguments.index('--snvcut') + 1])
	else:
		snv_cutoff = 5
	outdir = arguments[arguments.index('--out') + 1]
	if outdir[-1] != "/":
		outdir += "/"
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	logfile = os.path.join(outdir, 'Xpression_collector.log')

	# Create a logger
	logger = logging.getLogger("Xpression_collector_logger")
	logger.setLevel(logging.DEBUG)  # Capture all levels of logs

	# Create a file handler to write logs to a file
	file_handler = logging.FileHandler(logfile)
	file_handler.setLevel(logging.DEBUG)  # Write all logs to file

	# Create a stream handler to print logs to console
	console_handler = logging.StreamHandler()
	console_handler.setLevel(logging.DEBUG)  # Show only INFO and above in console

	# Create a common log format
	formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
	file_handler.setFormatter(formatter)
	console_handler.setFormatter(formatter)

	# Step 5: Add both handlers to the logger
	logger.addHandler(file_handler)
	logger.addHandler(console_handler)

	if '--sra' in arguments:
		sradir=os.path.join(outdir,'SRA')
		if not os.path.exists(sradir):
			os.makedirs(sradir)
	tmpdir = os.path.join(outdir, 'TMP')
	if not os.path.exists(tmpdir):
		os.makedirs(tmpdir)
		if tmpdir[-1] != '/':
			tmpdir += "/"

	#Option for user to give full path to blast
	#only BLAST offers blastn DIAMOND does not offer blastn option, it offers just protein vs protein or nucleotide vs protein alignments
	if '--blast' in arguments:
		aligner = arguments[arguments.index('--blast')+1]
		aligner = os.path.join(aligner, '')#to ensure that / is automatically added after the path is mentioned by the user if it is not already there
		aligner = aligner.replace(os.sep, '/')#to ensure adding only forward slash for ubuntu and windows systems
		makeblastdb = os.path.join(aligner + 'makeblastdb')
		blastn = os.path.join(aligner + 'blastn')
	else:
		makeblastdb = 'makeblastdb'
		blastn = 'blastn'

	#Option for user to give full path to diamond
	if '--diamond' in arguments:
		aligner = arguments[arguments.index('--diamond')+1]
		aligner = os.path.join(aligner, '')#to ensure that / is automatically added after the path is mentioned by the user if it is not already there
		aligner = aligner.replace(os.sep, '/')#to ensure adding only forward slash for ubuntu and windows systems
		diamond = os.path.join(aligner + 'diamond')
	else:
		diamond = 'diamond'

	#evalue for self local alignment in the isoform purger steps
	if '--eval' in arguments:
		eval = arguments[arguments.index('--eval')+1]
	else:
		eval = '1e-10'

	#Option for user to give full path to MAFFT
	if '--mafft' in arguments:
		mafft = arguments[arguments.index('--mafft')+1]#full path to mafft including the mafft.bat file
	else:
		mafft = 'mafft' #defaults to v7.526 of MAFFT - the most recent version of MAFFT while developing this script
	kallistofinaldir = os.path.join(outdir, 'Kallisto_run_final')
	if not os.path.exists(kallistofinaldir):
		os.makedirs(kallistofinaldir)

	#flag to obtain number of attempts at prefetching and accession from SRA; Default is just 1 attempt
	if '--attempts' in arguments:
		attempts = int(arguments[arguments.index('--attempts')+1])
	else:
		attempts = 1

	#flag to obtain the base wait time for sleep in case of network delays for prefetch
	if '--wait' in arguments:
		base_wait = arguments[arguments.index('--wait')+1]
	else:
		base_wait = 20 #20 seconds is the default base wait time

	if '--uncompressed' in arguments:
		readfile_status = None
	else:
		readfile_status = 'gz'

	if '--min_sra_file_size' in arguments:#enter file size in Mb; default is 1 Mb
		minimum_sra_file_size_threshold = int(arguments[arguments.index('--min_sra_file_size')+1])
	else:
		minimum_sra_file_size_threshold = 1

	merge_filtered_tpm= []
	merge_filtered_repr_tpm = []
	if '--merge_tpms' in arguments:#full path to config file containing the full paths to the filtered_tpm, and/or repr_filtered_tpms to be merged; first column will have paths to filtered_tpm files one per line; second column will have paths to filtered_repr_tpm files one per line
		tpm_config = arguments[arguments.index('--merge_tpms')+1]
		with open (tpm_config, 'r') as f:
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				merge_filtered_tpm.append(parts[0])
				if len(parts)>1:
					merge_filtered_repr_tpm.append(parts[1])
				line=f.readline()


	logger.info("Welcome to Xpression_collector!")

	#code block to obtain PEP file from CDS file

	input_spec = cdsfile
	output_spec = os.path.join(outdir,f'{orgname}.pep.fasta')
	internal_stop_to_x = ('--internal-stop-to-x' in arguments)

	# Standard code table (nuclear, 1)
	genetic_code = {
		'CTT':'L','ATG':'M','AAG':'K','AAA':'K','ATC':'I','AAC':'N','ATA':'I','AGG':'R',
		'CCT':'P','ACT':'T','AGC':'S','ACA':'T','AGA':'R','CAT':'H','AAT':'N','ATT':'I',
		'CTG':'L','CTA':'L','CTC':'L','CAC':'H','ACG':'T','CCG':'P','AGT':'S','CAG':'Q',
		'CAA':'Q','CCC':'P','TAG':'*','TAT':'Y','GGT':'G','TGT':'C','CGA':'R','CCA':'P',
		'TCT':'S','GAT':'D','CGG':'R','TTT':'F','TGC':'C','GGG':'G','TGA':'*','GGA':'G',
		'TGG':'W','GGC':'G','TAC':'Y','GAG':'E','TCG':'S','TTA':'L','GAC':'D','TCC':'S',
		'GAA':'E','TCA':'S','GCA':'A','GTA':'V','GCC':'A','GTC':'V','GCG':'A','GTG':'V',
		'TTC':'F','GTT':'V','GCT':'A','ACC':'T','TTG':'L','CGT':'R','TAA':'*','CGC':'R'
	}

	inputs = gather_inputs(input_spec)

	# If multiple inputs but output_spec is a file path, abort to avoid overwriting
	if len(inputs) > 1 and (not output_spec.endswith(os.sep) and not os.path.isdir(output_spec)):
		sys.exit("[ERROR] For multiple input files, --out must be a directory or end with '/'.")
	if output_spec.endswith(os.sep):
		os.makedirs(output_spec, exist_ok=True)

	for in_file in inputs:
		pep_file = make_output_path(output_spec, in_file)
		os.makedirs(os.path.dirname(pep_file), exist_ok=True)
		if not os.path.exists(pep_file):
			logger.info("Translating CDS file to PEP file")
			translate_file(logger, in_file, pep_file, genetic_code, internal_stop_to_x=internal_stop_to_x)

	# code to do BUSCO based QC
	if qc == 'yes':
		busco_dir_final = os.path.join(outdir, 'BUSCO_QC')
		if not os.path.exists(busco_dir_final):
			os.makedirs(busco_dir_final)
		busco_dir = os.path.join(outdir, 'BUSCO_DB')
		if not os.path.exists(busco_dir):
			os.makedirs(busco_dir)
		host_cache_dir = os.path.join(os.path.dirname(busco_dir), "BUSCO_cache")
		os.makedirs(host_cache_dir, exist_ok=True)
		busco_cmd = detect_busco(busco_path)
		if busco_cmd:
			busco_qc_result = os.path.join(busco_dir_final, 'BUSCO_QC.tsv')
			busco_db_dir = os.path.join(busco_dir, 'BUSCO_databases')
			busco_path_generated = "NA"
			container_workdir = busco_dir
			container_cache_dir = host_cache_dir
			if not os.path.exists(busco_qc_result):
				logger.info("Attempting BUSCO QC of PEP file generated")
				if busco_path == 'busco_docker':
					container_db_dir = busco_db_dir
					docker_image = f"ezlabgva/busco:{busco_version}_{container_version}"
					# writing bash script for executing busco installed with docker
					bash_script = f"""#!/bin/bash
							# BUSCO Docker wrapper script with proper cleanup
	
							# Container name based on process ID and timestamp for uniqueness
							CONTAINER_NAME="busco_$$_$(date +%s)_$(shuf -i 1000-9999 -n 1)"
	
							# Cleanup function
							cleanup() {{
								echo "Cleaning up Docker container: $CONTAINER_NAME" >&2
								docker stop "$CONTAINER_NAME" 2>/dev/null || true
								docker rm "$CONTAINER_NAME" 2>/dev/null || true
								exit $1
							}}
	
							# Set up signal handlers for proper cleanup
							trap 'cleanup 130' INT    # Ctrl+C (SIGINT)
							trap 'cleanup 143' TERM   # Termination (SIGTERM)
							trap 'cleanup 1' EXIT     # Any exit
							trap 'cleanup 1' ERR      # Any error
	
							# Check if Docker is available
							if ! command -v docker &>/dev/null; then
								echo "Error: Docker is not installed or not in PATH." >&2
								exit 1
							fi
	
							# Run BUSCO in Docker with automatic cleanup
							docker run --rm \\
								--name "$CONTAINER_NAME" \\
								-u $(id -u) \\
								-v "{host_path}:{container_path}" \\
								-v "{host_cache_dir}:{container_cache_dir}" \\
								-v "{busco_db_dir}:{container_db_dir}" \\
								-e XDG_CONFIG_HOME="{container_cache_dir}" \\
								-w "{container_workdir}" \\
								{docker_image} \\
								busco "$@"
	
							# Capture exit code and exit cleanly
							EXIT_CODE=$?
							exit $EXIT_CODE
							"""
					output_filename = os.path.join(busco_dir_final, "run_busco_docker.sh")
					with open(output_filename, "w") as f:
						f.write(bash_script)
					os.chmod(output_filename, 0o755)
					busco_path_generated = output_filename
					# Pre-download for Docker
					if not os.path.isdir(busco_db_dir):
						pre_download_databases(logger, busco_path_generated, org_type, busco_dir, busco_db_dir)
				else:
					if not os.path.isdir(busco_db_dir):
						# Pre-download for normal BUSCO
						pre_download_databases(logger, busco_path, org_type, busco_dir, busco_db_dir)
				logger.info("starting BUSCO run")
				ploidy_results = []
				busco_single_copy_lists = {}
				if busco_path != 'busco_docker':
					busco_result = run_busco(logger, orgname, str(pep_file), busco_path, busco_dir, busco_dir_final,str(cores), org_type, host_cache_dir,busco_db_dir, buscolineage)
				elif busco_path == 'busco_docker':
					busco_result = run_busco(logger, orgname, str(pep_file), busco_path_generated, busco_dir,busco_dir_final, str(cores), org_type,container_cache_dir, container_db_dir, buscolineage)
				with open(busco_qc_result, 'w') as out:
					out.write('Organism' + '\t' + 'Pseudo ploidy number' + '\t' + 'BUSCO Completeness (%)' + '\t' + 'BUSCO Duplication (%)' + '\n')
					out.write(str(busco_result) + '\n')
				logger.info(f"BUSCO QC check completed.")
		else:
			logger.warning("BUSCO not found. No QC check possible.")

	#code block to record completed accessions to tackle internet and network disruption interruptions
	completed_accessions = set()
	completed_accessions_file = os.path.join(outdir,'Fetch_completed_accession.txt')
	failed_accessions_file = os.path.join(outdir,'Fetch_failed_accession.txt')
	if os.path.exists(completed_accessions_file):
		with open (completed_accessions_file, 'r') as f:
			completed_accessions = set(line.strip() for line in f)
	# Clear failed accessions file at start of each run
	with open(failed_accessions_file, 'w') as f:
		pass
	if '--sra' in arguments:
		sra_accessions = [acc for acc in load_IDs(todo_sras) if acc not in completed_accessions]
	elif '--readfiles' in arguments:
		todo_sras_path = Path(todo_sras)
		sra_accessions = [f for f in todo_sras_path.iterdir() if f.is_dir() and f.name not in completed_accessions]#list of immediate subfolders with FASTQ read files
	batch_counter = 1
	while sra_accessions:
		batch = sra_accessions[:batch_size]  # take first N
		sra_accessions = sra_accessions[batch_size:]  # consume N
		logger.info(f"Processing batch {batch_counter} ({len(batch)} SRAs)")

		# Create Kallisto dir
		kallistodir = os.path.join(outdir, f'Kallisto_run_{batch_counter}')
		os.makedirs(kallistodir, exist_ok=True)
		if kallistodir[-1] != "/":
			kallistodir += "/"
		#create kallisto index file
		index_file = os.path.join(tmpdir, "index")
		if not os.path.isfile(index_file):
			logger.info("Starting Kallisto indexing")
			cmd1 = " ".join([kallisto, "index", "--index=" + index_file, "--make-unique", cds_file])
			p = subprocess.Popen(args=cmd1, shell=True)
			p.communicate()

		if '--sra' in arguments:
			# start kallisto consumer thread first
			kt = Thread(target=kallisto_worker, args=(index_file,sradir,readfile_status, logger,kallistodir,kallisto,cds_file,cores,tmpdir))
			kt.start()
			for accession in batch:
				fetch_worker(attempts, sradir, accession, prefetch_command, minimum_sra_file_size_threshold,fasterq_dump, cores, logger, completed_accessions_file, base_wait, failed_accessions_file)
			# signal consumer that no more accessions are coming
			fetch_queue.put(barrier)
			kt.join()

		elif '--readfiles' in arguments:
			#creating a temp folder containing symlinks to the subfolders in the specific batch
			tmp_dir_obj = tempfile.TemporaryDirectory(dir=outdir)#creating tmp_dir to store symlinks in the output folder specified by the user
			tmp_path = Path(tmp_dir_obj.name)
			for subfolder in batch:
				symlink = tmp_path / subfolder.name
				symlink.symlink_to(subfolder.resolve())
			sradir = tmp_path
			# Run Kallisto
			# --- load data --- #
			single_read_file_folders = [x[0] for x in os.walk(sradir, followlinks=True)][1:]
			logger.info("Number of FASTQ file folders detected: " + str(len(single_read_file_folders)) + "\n")

			# --- prepare jobs to run --- #
			index_file = os.path.join(tmpdir, "index")
			jobs_to_run = get_data_for_jobs_to_run(readfile_status, logger, single_read_file_folders, kallistodir, index_file, tmpdir)
			logger.info("Number of jobs to run: " + str(len(jobs_to_run)) + "\n")

			# --- generate index --- #
			if not os.path.isfile(index_file):
				logger.info("Starting Kallisto indexing and quantification per batch")
				cmd1 = " ".join([kallisto, "index", "--index=" + index_file, "--make-unique", cds_file])
				p = subprocess.Popen(args=cmd1, shell=True)
				p.communicate()

			# --- run jobs --- #
			job_executer(logger, jobs_to_run, kallisto, cores)

		# Merge the TPM, Counts numbers of SRA samples in a batch into a single TPM, Counts file respectively
		tpmfile = os.path.join(kallistofinaldir, f"TPM_{batch_counter}.txt")
		countsfile = os.path.join(kallistofinaldir, f"Counts_{batch_counter}.txt")

		counttables = glob.glob(os.path.join(kallistodir, "*.tsv"))
		count_data = {}
		tpm_data = {}
		transcript2gene = {}
		for filename in counttables:
			ID = filename.split('/')[-1].split('.')[0]
			if "_" in ID:  # only take ID if datetime string was included in file name
				ID = ID.split("_")[-1]
			counts, tpms = load_counttable(filename)
			# TPM are available and could be processed in the same way
			gene_counts = map_counts_to_genes(logger, transcript2gene, counts)
			gene_tpms = map_counts_to_genes(logger, transcript2gene, tpms)
			count_data.update({ID: gene_counts})
			tpm_data.update({ID: gene_tpms})

		if countsfile:
			generate_output_file(countsfile, count_data)
		if tpmfile:
			generate_output_file(tpmfile, tpm_data)

		# Clean SRA, TMP folders for next batch
		if '--sra' in arguments:
			shutil.rmtree(sradir)
			os.makedirs(sradir)
			for subdir in glob.glob(os.path.join(tmpdir, "SRR*")):
				shutil.rmtree(subdir)
		elif '--readfiles' in arguments:#removing the temporary folder with symlinks in case user wants to use own RNA-seq data read files
			tmp_dir_obj.cleanup()
		batch_counter += 1

	#code lines to merge the TPM/ Counts files
	# Get all TPM files from input folder
	patterns = ['*TPM*.txt', '*Counts*.txt']
	for pattern in patterns:
		search_pattern = os.path.join(kallistofinaldir, pattern)
		files = sorted(glob.glob(search_pattern))

		if not files:
			logger.error(f"Error: No TPM/ counts files found in '{kallistofinaldir}'!")
			logger.info(f"Searched for pattern: {search_pattern}")
			sys.exit(1)

		# Write combined output
		if pattern == '*TPM*.txt':
			unfiltered_output_file = os.path.join(tmpdir, f'{orgname}_unfiltered.tpms.tsv')
			if os.path.exists(unfiltered_output_file):
				continue
			else:
				logger.info("Merging TPM files from all batches")
		elif pattern == '*Counts*.txt':
			unfiltered_output_file = os.path.join(tmpdir, f'{orgname}_unfiltered.counts.tsv')
			if os.path.exists(unfiltered_output_file):
				break
			else:
				logger.info("Merging counts files from all batches")
		logger.info(f"\nWriting combined file to {unfiltered_output_file}...")

		logger.info(f"Input folder: {kallistofinaldir}")
		logger.info(f"Output folder: {outdir}")
		logger.info(f"Found {len(files)} files matching the pattern {pattern}")
		logger.info(f"Processing: {', '.join([os.path.basename(f) for f in files[:3]])}{'...' if len(files) > 3 else ''}\n")

		# Process first file to establish gene order
		logger.info(f"[1/{len(files)}] Reading {os.path.basename(files[0])} (establishing gene order)...")
		date1, sra_list, master_gene_data = read_tpm_file(files[0])
		all_sra_accessions = sra_list.copy()
		gene_order = list(master_gene_data.keys())  # Preserve gene order from first file

		logger.info(f"  - Found {len(gene_order)} genes")
		logger.info(f"  - Found {len(sra_list)} SRA samples")

		# Process remaining files
		for idx, tpm_file in enumerate(files[1:], start=2):
			logger.info(f"[{idx}/{len(files)}] Processing {os.path.basename(tpm_file)}...")
			date_info, sra_list, gene_data = read_tpm_file(tpm_file)
			# Validation checks
			if set(gene_data.keys()) != set(gene_order):
				missing_in_new = set(gene_order) - set(gene_data.keys())
				extra_in_new = set(gene_data.keys()) - set(gene_order)

				if missing_in_new:
					logger.warning(f"  WARNING: {len(missing_in_new)} genes missing in {os.path.basename(tpm_file)}")
					logger.warning(f"    First 5 missing: {list(missing_in_new)[:5]}")
				if extra_in_new:
					logger.warning(f"  WARNING: {len(extra_in_new)} extra genes in {os.path.basename(tpm_file)}")
					logger.warning(f"    First 5 extra: {list(extra_in_new)[:5]}")

				# Handle missing genes by adding empty values
				for gene in gene_order:
					if gene not in gene_data:
						gene_data[gene] = ['NA'] * len(sra_list)

			# Add expression values to master data (in correct gene order)
			for gene in gene_order:
				if gene in gene_data:
					master_gene_data[gene].extend(gene_data[gene])
				else:
					# Gene missing in this file - add NAs
					master_gene_data[gene].extend(['NA'] * len(sra_list))

			all_sra_accessions.extend(sra_list)
			logger.info(f"  - Added {len(sra_list)} SRA samples")

		with open(unfiltered_output_file, 'w') as out:
			# Write header row
			out.write('gene\t' + '\t'.join(all_sra_accessions) + '\n')

			# Write gene expression data (in original gene order)
			for gene in gene_order:
				out.write(gene + '\t' + '\t'.join(master_gene_data[gene]) + '\n')

		logger.info(f"\n Success!")
		logger.info(f"  - Total genes: {len(gene_order)}")
		logger.info(f"  - Total SRA samples before filtering: {len(all_sra_accessions)}")
		logger.info(f"  - Unfiltered output file: {unfiltered_output_file}")

	#code block to filter RNA-seq samples
	filtered_tpm_file = os.path.join(outdir, f"{orgname}.tpms.tsv")
	if not os.path.exists(filtered_tpm_file):
		logger.info("Filtering the TPM expression file")
		for f in os.listdir(tmpdir):
			if 'tpms' in f:
				unfiltered_tpm_file = os.path.join(tmpdir, f)
			elif 'counts' in f:
				unfiltered_counts_file = os.path.join(tmpdir, f)
		# --- run analysis of all data in folder/file --- #
		doc_file = os.path.join(tmpdir,f'{orgname}_qc.doc')
		valid_samples = []
		with open(doc_file, "w") as out:
			out.write("SampleName\tPercentageOfTop100\tPercentageOfTop500\tPercentageOfTop1000\n")
			TPM_data, genes = load_all_TPMs(unfiltered_tpm_file)
			count_data, genes = load_all_TPMs(unfiltered_counts_file)
			for key in sorted(list(TPM_data.keys())):
				new_line = [key]
				selection = sorted(TPM_data[key])
				counts = sum(count_data[key])  # calculate counts per library
				if counts >= min_counts:  # check for sufficient library size
					try:  # check for ID presence on black list
						black_list[key]
						new_line.append("ID on black list")
						out.write("\t".join(list(map(str, new_line))) + "\n")
					except KeyError:
						try:
							val = 100.0 * sum(selection[-100:]) / sum(selection)
						except ZeroDivisionError:
							val = 0
						new_line.append(val)
						if min_cutoff < val < max_cutoff:
							valid_samples.append(key)
						if len(selection) > 500 and val > 0:
							new_line.append(100.0 * sum(selection[-500:]) / sum(selection))
						else:
							new_line.append("n/a")
						if len(selection) > 1000 and val > 0:
							new_line.append(100.0 * sum(selection[-1000:]) / sum(selection))
						else:
							new_line.append("n/a")
						out.write("\t".join(list(map(str, new_line))) + "\n")
				else:
					new_line.append("insufficient counts: " + str(counts))
					out.write("\t".join(list(map(str, new_line))) + "\n")

		logger.info("number of valid sample: " + str(len(valid_samples)))
		logger.info("number of invalid sample: " + str(len(TPM_data.keys()) - len(valid_samples)))

		# --- generate output file --- #
		if len(valid_samples) > 0:
			with open(filtered_tpm_file, "w") as out:
				out.write("gene\t" + "\t".join(valid_samples) + "\n")
				for idx, gene in enumerate(genes):
					new_line = [gene]
					for sample in valid_samples:
						new_line.append(TPM_data[sample][idx])
					out.write("\t".join(list(map(str, new_line))) + "\n")
		else:
			logger.error("WARNING: no valid samples in data set!")
		# --- generate figure --- #
		fig_file = os.path.join(tmpdir, f'{orgname}_qc.pdf')
		values = []
		with open(doc_file, "r") as f:
			f.readline()  # remove header
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				try:
					values.append(float(parts[1]))
				except ValueError:
					pass
				line = f.readline()

		values = [x for x in values if str(x) != 'nan']

		fig, ax = plt.subplots()

		ax.hist(values, bins=100, color="green")
		ax.set_xlabel("Percentage of expression on top100 genes")
		ax.set_ylabel("Number of analyzed samples")

		fig.savefig(fig_file)
		if remove_isoforms == 'no':
			logger.info(f'Xpression_collector pipeline completed successfully!!!')

	#optional removal of isoforms
	isoform_reduced_cds_file = outdir + f"{orgname}_repr.cds.fasta"
	repr_output_spec = os.path.join(outdir, f'{orgname}_repr.pep.fasta')
	repr_tpm_file = os.path.join(outdir, f'{orgname}_repr.tpms.tsv')
	repr_counts_file = os.path.join(outdir, f'{orgname}_repr.counts.tsv')
	if remove_isoforms == 'yes':
		if os.path.exists(isoform_reduced_cds_file) and os.path.exists(repr_output_spec) and os.path.exists(repr_tpm_file) and os.path.exists(repr_counts_file):
			pass
		else:
			logger.info("Removing alternative isoforms")
			# --- clean file --- #
			clean_file = os.path.join(tmpdir, "clean.fasta")
			seqs = load_fasta(cdsfile)
			with open(clean_file, "w") as out:
				for key in list(seqs.keys()):
					out.write('>' + key + "\n" + seqs[key] + "\n")

			# --- run BLAST/DIAMOND vs self --- #
			dbname = os.path.join(tmpdir, "blastdb")
			blast_result_file = os.path.join(tmpdir, "blast_results.txt")
			cmd = f"{makeblastdb} -in " + clean_file + " -out " + dbname + " -dbtype nucl"
			p = subprocess.Popen(args=cmd, shell=True)
			p.communicate()

			cmd = f"{blastn} -query " + clean_file + " -db " + dbname + " -out " + blast_result_file + " -outfmt 6 -evalue " + str(eval) + " -num_threads " + str(cores)
			p = subprocess.Popen(args=cmd, shell=True)
			p.communicate()
			"""
			elif aligner_tool == 'diamond':
				cmd = diamond + ' makedb --in ' + clean_file + ' --db ' + dbname + ' --quiet '
				p = subprocess.Popen(args=cmd, shell=True)
				p.communicate()

				cmd = diamond + ' blastn --db ' + dbname + ' --evalue ' + str(eval) + ' --query ' + clean_file + ' --out ' + blast_result_file + ' --outfmt 6 --quiet --ultra-sensitive --threads ' + str(cores)
				p = subprocess.Popen(args=cmd, shell=True)
				p.communicate()
			"""

			# --- identify good BLAST hits around sequences --- #
			blast_hits = load_hits_per_bait(blast_result_file, scorecut, simcut, lencut)
			black_list = {}
			groups_to_analyze = []
			for key in list(blast_hits.keys()):
				try:
					black_list[key]  # check if sequence ID is already assigned to group
				except KeyError:
					group = blast_hits[key] + [key]  # combine key with all good BLAST hits
					if len(group) > 1:
						groups_to_analyze.append(group)
					for each in group:
						black_list.update({each: None})
			logger.info("number of groups to analyze: " + str(len(groups_to_analyze)))
			single_fasta_folder = os.path.join(tmpdir, "single_fasta_input/")
			if not os.path.exists(single_fasta_folder):
				os.makedirs(single_fasta_folder)
			for idx, group in enumerate(groups_to_analyze):
				output_file_name = single_fasta_folder + str(idx + 1) + ".fasta"
				if not os.path.isfile(output_file_name):
					with open(output_file_name, "w") as out:
						for ID in group:
							out.write('>' + ID + "\n" + seqs[ID] + "\n")

			# --- construct global alignment --- #
			fasta_files = glob.glob(single_fasta_folder + "*.fasta")
			for fasta in fasta_files:
				if not os.path.isfile(fasta + ".aln"):
					cmd = f"{mafft} " + fasta + " > " + fasta + ".aln 2> " + fasta + ".doc.txt"
					p = subprocess.Popen(args=cmd, shell=True)
					p.communicate()

			# --- calculate similarity matrix --- #
			aln_fasta_files = glob.glob(single_fasta_folder + "*.fasta.aln")
			logger.info("number of aligned FASTA files: " + str(len(aln_fasta_files)))
			alignment_results = []
			for aln_fasta in aln_fasta_files:
				alignment = load_fasta(aln_fasta)
				sim_per_aln = load_aln_similarity(alignment, snv_cutoff)
				alignment_results.append(sim_per_aln)

			# --- classify sequences as isoforms/paralogs --- #
			isoforms = identify_isoforms(alignment_results, snv_cutoff)
			logger.info("number of isoform groups: " + str(len(isoforms)))
			repr_isoforms, repr_ids = identify_repr_isoform_per_group(logger, isoforms, seqs)
			blacklist = {}
			for ID in [x for sublist in isoforms for x in sublist]:
				blacklist.update({ID: None})

			# --- write isoform-reduced output file --- #
			with open(isoform_reduced_cds_file, "w") as out:
				for key in list(seqs.keys()):
					try:
						blacklist[key]
					except KeyError:
						out.write('>' + key + "\n" + seqs[key] + "\n")
				for key in list(repr_isoforms.keys()):
					out.write('>' + key + "\n" + seqs[key] + "\n")

			#code block to produce PEP file without isoforms
			repr_input_spec = isoform_reduced_cds_file
			internal_stop_to_x = ('--internal-stop-to-x' in arguments)

			# Standard code table (nuclear, 1)
			genetic_code = {
				'CTT': 'L', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R',
				'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'ACA': 'T', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I',
				'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'AGT': 'S', 'CAG': 'Q',
				'CAA': 'Q', 'CCC': 'P', 'TAG': '*', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CCA': 'P',
				'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G',
				'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S', 'TTA': 'L', 'GAC': 'D', 'TCC': 'S',
				'GAA': 'E', 'TCA': 'S', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V',
				'TTC': 'F', 'GTT': 'V', 'GCT': 'A', 'ACC': 'T', 'TTG': 'L', 'CGT': 'R', 'TAA': '*', 'CGC': 'R'
			}

			inputs = gather_inputs(repr_input_spec)

			# If multiple inputs but output_spec is a file path, abort to avoid overwriting
			if len(inputs) > 1 and (not repr_output_spec.endswith(os.sep) and not os.path.isdir(repr_output_spec)):
				sys.exit("[ERROR] For multiple input files, --out must be a directory or end with '/'.")
			if repr_output_spec.endswith(os.sep):
				os.makedirs(repr_output_spec, exist_ok=True)

			for in_file in inputs:
				pep_file = make_output_path(repr_output_spec, in_file)
				os.makedirs(os.path.dirname(pep_file), exist_ok=True)
				translate_file(logger, in_file, pep_file, genetic_code, internal_stop_to_x=internal_stop_to_x)

			#code block to produce TPM file without alternative isoforms
			repr_tpm_filtered = keep_primary_transcript_exp(repr_ids, repr_tpm_file, repr_counts_file, isoform_reduced_cds_file, filtered_tpm_file)
			logger.info(f'Xpression_collector pipeline completed successfully!!!')

	# code block to merge filtered TPM file produced in this run with already existing filtered TPM files
	try:
		timestr = time.strftime("%Y_%m_%d_")
		merged_filtered_tpm_file = os.path.join(outdir, f'{timestr}_{orgname}_merged.tpms.tsv')
		merged_filtered_repr_tpm_file = os.path.join(outdir, f'{timestr}_{orgname}_merged_repr.tpms.tsv')
	except ModuleNotFoundError:
		merged_filtered_tpm_file = os.path.join(outdir, f'{orgname}_merged.tpms.tsv')
		merged_filtered_repr_tpm_file = os.path.join(outdir, f'{orgname}_merged_repr.tpms.tsv')

	if merge_filtered_tpm:
		try:
			merged_df = merge_expression_tsvs(filtered_tpm_file, merge_filtered_tpm)
			merged_df.to_csv(merged_filtered_tpm_file, sep='\t', index=False)
			print(f"Merge successful: {merged_df.shape[0]} genes x {merged_df.shape[1]} columns")
		except ValueError as e:
			print(f"Aborting merge:\n{e}")
	if merge_filtered_repr_tpm and remove_isoforms == 'yes':
		try:
			merged_df = merge_expression_tsvs(repr_tpm_filtered, merge_filtered_repr_tpm)
			merged_df.to_csv(merged_filtered_repr_tpm_file, sep='\t', index=False)
			print(f"Merge successful: {merged_df.shape[0]} genes x {merged_df.shape[1]} columns")
		except ValueError as e:
			print(f"Aborting merge:\n{e}")

if '--sra' in sys.argv and '--cds' in sys.argv and '--out' in sys.argv:
	main(sys.argv)
elif '--readfiles' in sys.argv and '--cds' in sys.argv and '--out' in sys.argv:
	main(sys.argv)
else:
	sys.exit(__usage__)
