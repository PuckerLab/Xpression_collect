### Shakunthala Natarajan ###
### bug reports: s64snata@uni-bonn.de ###

__version__=0.2
__usage__="""
			python3 Xpression_collector.py
			--cds <Full path to CDS file>
			--sample_name <Name of the sample you are analyzing>
			--sra <Full path to TXT file with one SRA accession per line>
			--min_sra_file_size <Minimum file size cutoff in MB to check prefetched SRA file sizes to cath sralite files that might cause downstream errors> default cutoff is 1MB
			--gff <Full path to GFF file if isoform removal needs to be performed>
			--gff_config <Full path to GFF config file specifying the different GFF attributes - child_attribute, child_parent_linker, parent_attribute> default attributes are ID, Parent, ID respectively
			--uncompressed <Provide this flag if your read files are uncompressed>
			--annotation_qc <yes or no for BUSCO-based QC of the PEP file> default is yes
			--threads <Total number of cores for running the pipeline> default is 4
			--parallel_prefetch <Number of SRA accessions to be prefetched parallely> default is 2
			--attempts <Number of attempts at prefetching and accession from SRA> default is just 1 attempt
			--wait <Base wait time for sleep in case of network delays for prefetch; Increases exponentially with a base of 2 for each reattempt> default is 20 seconds
			--remove_isoforms <Optional step to remove isoforms>< yes or no> default is yes
			--merge_tpms <Full path to config file containing the full paths to the filtered_tpm, and/or repr_filtered_tpms to be merged>
						 <First column will have paths to filtered_tpm files one per line; Second column will have paths to filtered_repr_tpm files one per line>
			--min <Minimum percentage expression of top 100 genes>
			--max <Maximum percentage expression of top 100 genes>
			--black <SRA IDs to be removed or blacklisted in a TXT file with one SRA accession ID per line>
			--busco <Full path to BUSCO> <Specify busco_docker if BUSCO is installed via docker>
			--busco_lineage <Specify the BUSCO lineage> default is auto
			--busco_version <Version of BUSCO> default is v6.0.0
			--container_version <Container version of the BUSCO docker image>
			--docker_host_path <Host path to be mounted on for running the docker image>
			--docker_container_path <Container path of the docker image>
			--organism_type <eukaryote or prokaryote> default is eukaryote
			--kallisto <Full path to the kallisto executable>
			--prefetch <Full path to the prefetch executable>
			--fasterq-dump <Full path to the fasterq-dump executable>
			--out <Full path to output directory>
			--clean_up <Clean up the TMP folder after a successful run> default is yes
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
import threading
from threading import Thread, Lock
from queue import Queue
from operator import itemgetter
from concurrent.futures import ProcessPoolExecutor, as_completed, TimeoutError
try:
	from pyfiglet import Figlet
	f = Figlet(font='doom',width=200)
	print(f.renderText('Xpression collector'))
except ImportError:
	pass

try:
	from rich.live import Live
	from rich.table import Table
	from rich.console import Console
	RICH_AVAILABLE = True
except ImportError:
	RICH_AVAILABLE = False

### --- end of imports --- ###

kallisto_queue = Queue()
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
#functions to provide rich display of progress using rich library components
status = {}
status_lock = Lock()  # prevents simultaneous writes from multiple threads

def update_status(accession, stage):
	with status_lock:
		status[accession] = stage

def build_dashboard():
	table = Table(title="Xpression_collector Live Status")
	table.add_column("Accession")
	table.add_column("Stage")

	# --- per accession rows ---
	for accession, stage in status.items():
		color = {
			"prefetching": "cyan",
			"prefetch failed": "red",
			"fasterq-dump & pigz": "yellow",
			"fasterq-dump & pigz failed": "red",
			"kallisto": "magenta",
			"kallisto failed": "red",
			"completed": "green"  # ← add this
		}.get(stage, "white")
		table.add_row(accession, f"[{color}]{stage}[/{color}]")

	# --- compute summary counts from status dictionary ---
	fully_completed = sum(1 for stage in status.values() if stage == "completed")
	failed_count = sum(1 for stage in status.values() if "failed" in stage)
	in_progress = sum(1 for stage in status.values() if stage != "completed" and "failed" not in stage)

	# --- summary row ---
	table.add_row("---", "---")
	table.add_row(
		"[bold]SUMMARY[/bold]",
		f"[green]completed:{fully_completed}[/green]  "
		f"[red]failed:{failed_count}[/red]  "
		f"[cyan]in progress:{in_progress}[/cyan]"
	)
	return table
def load_multiple_fasta_file(fasta_file):
	"""Load all sequences from a (possibly wrapped) FASTA file into a dict."""
	content = {}
	header = None
	seq_chunks = []
	#uncompressed FASTA file
	if fasta_file[-2:].lower() != 'gz':
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
	else:
		with gzip.open(fasta_file, "rt") as f:
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
		final_result_file = os.path.join(final_output_folder, f'{timestr}{ID}.tsv')
		if os.path.isfile(final_result_file):
			status = False
		if status:
			jobs_to_do.append(
				{'r1': read_file1, 'r2': read_file2, 'out': output_dir, 'index': index_file, 'tmp': tmp_result_file,
				 'fin': final_result_file, "ID": ID})
	return jobs_to_do

def job_executer(accession, logger, jobs_to_run, kallisto, threads, kallisto_completed_accessions,failed_accessions_file,subprocess_log):
	"""! @brief run all jobs in list """
	try:
		for idx, job in enumerate(jobs_to_run):
			logger.info("running job " + str(idx + 1) + "/" + str(len(jobs_to_run)) + " - " + job["ID"] + "\n")

			if job['r2']:
				cmd2 = " ".join([kallisto, "quant", "--index=" + job['index'], "--output-dir=" + job['out'], "--threads " + str(threads), job['r1'], job['r2']])
			else:
				cmd2 = " ".join([kallisto, "quant", "--index=" + job['index'], "--single -l 200 -s 100", "--output-dir=" + job['out'], "--threads " + str(threads), job['r1']])
			result = subprocess.run(cmd2, shell=True, capture_output=True, text=True)
			output = (result.stdout + result.stderr).replace('\r', '\n')
			subprocess_log.write(output)
			subprocess_log.flush()

			p = subprocess.Popen(args="cp " + job["tmp"] + " " + job["fin"], shell=True, stdout=subprocess_log, stderr=subprocess_log)
			p.communicate()
		if p.returncode !=0:
			with open(failed_accessions_file, 'a') as out:
				out.write(f'{accession}\n')
				out.flush()
				os.fsync(out.fileno())
			raise RuntimeError(f"Kallisto quantification did not complete successfully for {accession}")
		else:
			with open(kallisto_completed_accessions, 'a') as out:
				out.write(f'{accession}\n')
				out.flush()
				os.fsync(out.fileno())
			update_status(accession, "completed")
	except RuntimeError as e:
		update_status(accession, "kallisto failed")
		logger.error(f"Kallisto quantification did not complete successfully for {accession}")

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

	samples = list(sorted(list(data.keys())))

	with open(output_file, "w") as out:
		out.write("\t".join(['gene'] + samples) + '\n')
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

#function to produce TPM file without isoforms
def keep_primary_transcript_exp(repr_ids, repr_tpm_file, primary_transcript_cds_file, exp_file):
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

	return repr_tpm_file

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

#function for removing alternate transcripts from the peptide FASTA file
def isoform_clean(gff3_input_file, cds_dict, no_trans_cds, child_attribute, child_parent_linker):
	repr_ids={}
	no_gene_no_parent = False
	has_gene = False
	has_parent =False
	if gff3_input_file[-2:].lower() != 'gz':
		with open(gff3_input_file, "r") as f:
			gff_lines = f.readlines()
			# checking if gene feature is present in the GFF file
			has_gene = any(line.split('\t')[2].upper() == 'GENE' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)

	else:
		with gzip.open(gff3_input_file, "rt") as f:
			gff_lines = f.readlines()
			# checking if gene feature is present in the GFF file
			has_gene = any(line.split('\t')[2].upper() == 'GENE' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
	if has_mrna:
		coding_feature = 'MRNA'
	else:
		has_transcript = any(line.split('\t')[2].upper() == 'TRANSCRIPT' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
		if has_transcript:
			coding_feature = 'TRANSCRIPT'
		else:
			has_cds = any(line.split('\t')[2].upper() == 'CDS' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
			if has_cds:
				coding_feature = 'CDS'
			else:
				coding_feature = 'EXON'

	nogene_noparent_counter = 0
	if gff3_input_file[-2:].lower() != 'gz':  # uncompressed gff file
		transcripts_per_gene = {}
		with open(gff3_input_file, "r") as f:
			line = f.readline()
			while line:
				if line[0] != "#":
					no_gene_no_parent = False  # Reset for each line
					parts = line.strip().split('\t')
					if len(parts) > 2:
						if parts[2].upper() == coding_feature:
							partsnew = parts[-1].strip().split(';')
							# Check if any attribute starts with 'Parent='
							has_parent = any(attr.startswith(str(child_parent_linker) + '=') for attr in partsnew)
							if has_gene and has_parent:
								nogene_noparent_counter += 1
								for each in partsnew:
									pattern_par = r'^' + re.escape(child_parent_linker + '=') + r'.*$'
									if re.match(pattern_par, each):
										partsnew1 = str(each).replace(str(child_parent_linker) + '=', "")
							for every in partsnew:
								pattern_ID = r'^' + re.escape(child_attribute + '=') + r'.*$'
								if re.match(pattern_ID, every):
									partsnew0 = str(every).replace(str(child_attribute) + '=', "")
							try:
								transcripts_per_gene[partsnew1].append(partsnew0)
							except KeyError:
								transcripts_per_gene.update({partsnew1: [partsnew0]})
				line = f.readline()

	else:#compressed gff file
		transcripts_per_gene = {}
		with gzip.open(gff3_input_file, "rt") as f:
			line = f.readline()
			while line:
				if line[0] != "#":
					no_gene_no_parent = False  # Reset for each line
					parts = line.strip().split('\t')
					if len(parts) > 2:

						if parts[2].upper() == coding_feature:
							partsnew = parts[-1].strip().split(';')
							# Check if any attribute starts with 'Parent='
							has_parent = any(attr.startswith(str(child_parent_linker)+'=') for attr in partsnew)
							if has_gene and has_parent:
								nogene_noparent_counter += 1
								for each in partsnew:
									pattern_par = r'^' + re.escape(child_parent_linker + '=') + r'.*$'
									if re.match(pattern_par, each):
										partsnew1 = str(each).replace(str(child_parent_linker)+'=', "")
							for every in partsnew:
								pattern_ID = r'^' + re.escape(child_attribute + '=') + r'.*$'
								if re.match(pattern_ID, every):
									partsnew0 = str(every).replace(str(child_attribute)+'=', "")
							try:
								transcripts_per_gene[partsnew1].append(partsnew0)
							except KeyError:
								transcripts_per_gene.update({partsnew1: [partsnew0]})
				line = f.readline()

	gene_names = list(transcripts_per_gene.keys())
	with open(no_trans_cds, "w") as out:
		for gene in gene_names:
			trans_length = []
			isoform_list = []
			for trans in transcripts_per_gene[gene]:
				if trans in cds_dict:
					if len(transcripts_per_gene[gene]) < 2:
						out.write('>' + str(trans) + '\n' + str(cds_dict[trans]) + '\n')
					else:
						trans_length.append((trans, cds_dict[trans]))
						isoform_list.append(trans)
			if trans_length:
				best_trans, seq = max(trans_length, key=lambda x: len(x[1]))
				out.write('>' + str(best_trans) + "\n" + str(seq) + "\n")
				isoform_list.remove(best_trans)
				repr_ids[best_trans] = isoform_list.copy()

	return repr_ids

# function to perform kallisto quantification as soon as SRA file is fetched
def fasterqdump_kallisto_worker(fasterqpigz_completed_accessions, kallisto_completed_accessions, index_file, sradir, readfile_status, logger, kallistodir, kallisto, cores,tmpdir, fasterq_dump, attempts, base_wait, failed_accessions_file,fasterqpigz_completed_accessions_file,kallisto_completed_accessions_file,subprocess_log):
	while True:
		accession = kallisto_queue.get()
		if accession is barrier:
			break
		acc_dir = os.path.join(sradir, accession)
		prefetched_file = os.path.join(acc_dir, f"{accession}.sra")
		# Run fasterq-dump and Kallisto
		if accession not in fasterqpigz_completed_accessions:
			update_status(accession, "fasterq-dump & pigz")
			for attempt in range(attempts):
				try:
					# fasterq-dump
					cmd = f"{fasterq_dump} --split-3 --outdir {acc_dir} --skip-technical --threads {cores} {prefetched_file}"
					fasterq_dump_result = subprocess.run(cmd, shell=True, stdout=subprocess_log, stderr=subprocess_log)
					if fasterq_dump_result.returncode != 0:
						raise RuntimeError(f"fasterq-dump failed for {accession} with error code {fasterq_dump_result.returncode}")

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
						pigz_result = subprocess.run(cmd, shell=True, stdout=subprocess_log, stderr=subprocess_log)
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
					# if fetching is successful mark as completed accession, remove the prefetched file (.sra) and break out of the retry loop
					rm_cmd = f'rm {prefetched_file}'
					p = subprocess.Popen(args=rm_cmd, shell=True)
					p.communicate()
					with open(fasterqpigz_completed_accessions_file, 'a') as out:
						out.write(f'{accession}\n')
						out.flush()
						os.fsync(out.fileno())
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
						update_status(accession, "faster-dump & pigz failed")
						logger.error(f"All attempts exhausted for {accession}, marking as failed")
						with open(failed_accessions_file, 'a') as out:
							out.write(f'{accession}\n')
							out.flush()
							os.fsync(out.fileno())

		if accession not in kallisto_completed_accessions:
			update_status(accession, "kallisto")
			# --- load data --- #
			acc_dir = os.path.join(sradir, accession)
			single_read_file_folders = [acc_dir]
			logger.info("Number of FASTQ file folders detected: " + str(len(single_read_file_folders)) + "\n")

			# --- prepare jobs to run --- #
			jobs_to_run = get_data_for_jobs_to_run(readfile_status, logger, single_read_file_folders, kallistodir,index_file, tmpdir)
			logger.info("Number of jobs to run: " + str(len(jobs_to_run)) + "\n")

			# --- run jobs --- #
			job_executer(accession, logger, jobs_to_run, kallisto, cores, kallisto_completed_accessions_file,failed_accessions_file,subprocess_log)
			#remove the accession SRA folder after kallisto is completed for that accession
			shutil.rmtree(acc_dir)
		kallisto_queue.task_done()

#function to control fetching of SRA files through parallelized prefetch
def parallel_prefetch(prefetch_completed_accessions, accession, attempts,sradir, prefetch_command, minimum_sra_file_size_threshold, base_wait, failed_accessions_file, prefetch_completed_accessions_file, logger,subprocess_log):
	if accession not in prefetch_completed_accessions:
		update_status(accession, "prefetch")
		for attempt in range(attempts):
			# Create subfolder for this accession
			acc_dir = os.path.join(sradir, accession)
			if os.path.exists(acc_dir):
				shutil.rmtree(acc_dir)
			os.makedirs(acc_dir, exist_ok=True)
			prefetched_file = os.path.join(acc_dir, f"{accession}.sra")
			try:
				cmd = f"{prefetch_command} --max-size 200G {accession} -O {sradir}"
				prefetch_result = subprocess.run(cmd, shell=True, stdout=subprocess_log, stderr=subprocess_log)
				if prefetch_result.returncode != 0:  # if prefetch fails due to issues lik network disruption
					raise RuntimeError(f"prefetch for {accession} failed with code {prefetch_result.returncode}")
				if not os.path.exists(prefetched_file):  # if prefetch did not fetch an SRA file in the first place
					raise RuntimeError(f"prefetch for {accession} did not fetch a .sra file")
				filesize = (os.path.getsize(prefetched_file)) / (1024 ** 2)  # converting the file size returned in bytes by getsize to Mb
				if filesize < minimum_sra_file_size_threshold:  # in case prefetch gets the sralite files
					raise RuntimeError(f"File size of the prefetched file for {accession} is smaller than the threshold size of {minimum_sra_file_size_threshold}")
				kallisto_queue.put(accession)
				with open(prefetch_completed_accessions_file, 'a') as out:
					out.write(f'{accession}\n')
					out.flush()
					os.fsync(out.fileno())
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
					update_status(accession, "prefetch failed")
					logger.error(f"All attempts exhausted for {accession}, marking as failed")
					with open(failed_accessions_file, 'a') as out:
						out.write(f'{accession}\n')
						out.flush()
						os.fsync(out.fileno())

def main(arguments):
	if '--sample_name' in arguments:
		orgname = arguments[arguments.index('--sample_name')+1]
	else:
		orgname = 'sample'
	if  '--sra' in arguments:#list of SRA accessions to be fetched from NCBI SRA
		todo_sras = arguments[arguments.index('--sra')+1]
	if '--gff' in arguments:
		gff_file=arguments[arguments.index('--gff')+1]
	# gff file config params
	if '--gff_config' in arguments:
		gff_config_file = arguments[arguments.index('--gff_config') + 1]
		with open(gff_config_file, 'r') as f:
			for line in f:
				parts = line.strip().split()
				child_attribute = parts[0]
				child_parent_linker = parts[1]
				parent_attribute = parts[2]
	else:
		child_attribute = 'ID'
		child_parent_linker = 'Parent'
		parent_attribute = 'ID'
	if '--fastq_pattern' in arguments:
		pattern_names = arguments[arguments.index('--fastq_pattern')+1]#specify fastq read file name pattern for the paired end files separated by commas without spaces like - _pass_1,_pass_2_
		pattern_names_list = pattern_names.split(',')
	else:
		pattern_names_list = ["_pass_1", "_pass_2"]

	cds_file = arguments[arguments.index('--cds')+1]#full path to CDS file

	if '--annotation_qc' in arguments:# yes or no for BUSCO-based QC of the PEP file
		qc = arguments[arguments.index('--annotation_qc')+1]
	else:
		qc = 'yes'

	# Option for user to give full path to busco
	if '--busco' in arguments:
		busco_path = arguments[arguments.index('--busco') + 1]  # busco_docker or default busco
	else:
		busco_path = "/tools/xpcollect_env/bin/busco"

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
		org_type = arguments[arguments.index('--organism_type') + 1]
	else:
		org_type = 'eukaryote'

	if '--kallisto' in arguments:  # full path to kallisto
		kallisto = arguments[arguments.index('--kallisto') + 1]
	else:
		kallisto = "/tools/xpcollect_env/bin/kallisto"

	if '--fasterq-dump' in arguments:  # full path to fasterq-dump from NCBI SRA toolkit
		fasterq_dump = arguments[arguments.index('--fasterq-dump') + 1]
	else:
		fasterq_dump = '/tools/sratoolkit/bin/fasterq-dump'

	if '--prefetch' in arguments:
		prefetch_command = arguments[arguments.index('--prefetch') + 1]
	else:
		prefetch_command = "/tools/sratoolkit/bin/prefetch"

	if '--threads' in arguments:
		cores = int(arguments[arguments.index('--threads')+1])
	else:
		cores = 4
	if '--parallel_prefetch' in arguments:
		batch_size=int(arguments[arguments.index('--parallel_prefetch')+1])
	else:
		batch_size = 2

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

	outdir = arguments[arguments.index('--out') + 1]
	if outdir[-1] != "/":
		outdir += "/"
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	if '--clean_up' in arguments:
		clean_up = arguments[arguments.index('--clean_up')+1]
	else:
		clean_up = 'yes'

	logfile = os.path.join(outdir, 'Xpression_collector.log')
	subprocess_log = open(logfile, 'a')

	# Create a logger
	logger = logging.getLogger("Xpression_collector_logger")
	logger.setLevel(logging.DEBUG)  # Capture all levels of logs

	# Create a file handler to write logs to a file
	file_handler = logging.FileHandler(logfile)
	file_handler.setLevel(logging.DEBUG)  # Write all logs to file

	# Create a common log format
	formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
	file_handler.setFormatter(formatter)


	# Step 5: Add both handlers to the logger
	logger.addHandler(file_handler)
	#logging to console as well only if rich library is not found otherwise just log to file
	if not RICH_AVAILABLE:
		console_handler = logging.StreamHandler()
		console_handler.setLevel(logging.DEBUG)
		console_handler.setFormatter(formatter)
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

	input_spec = cds_file
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
	prefetch_completed_accessions_file = os.path.join(tmpdir,'prefetch_completed_accession.txt')
	fasterqpigz_completed_accessions_file = os.path.join(tmpdir,'fasterq-dump_pigz_completed_accession.txt')
	kallisto_completed_accessions_file = os.path.join(tmpdir, 'kallisto_quant_completed_accessions.txt')

	failed_accessions_file = os.path.join(outdir,'failed_accessions.txt')
	kallistodir=os.path.join(tmpdir,'Kallisto_abundances')
	if not os.path.exists(kallistodir):
		os.makedirs(kallistodir)

	prefetch_completed_accessions = set()
	fasterqpigz_completed_accessions = set()
	kallisto_completed_accessions = set()

	if os.path.exists(prefetch_completed_accessions_file):
		with open (prefetch_completed_accessions_file, 'r') as f:
			prefetch_completed_accessions = set(line.strip() for line in f)

	if os.path.exists(fasterqpigz_completed_accessions_file):
		with open (fasterqpigz_completed_accessions_file, 'r') as f:
			fasterqpigz_completed_accessions = set(line.strip() for line in f)

	if os.path.exists(kallisto_completed_accessions_file):
		with open (kallisto_completed_accessions_file, 'r') as f:
			kallisto_completed_accessions = set(line.strip() for line in f)


	# Clear failed accessions file at start of each run
	with open(failed_accessions_file, 'w') as f:
		pass
	sra_accessions = [acc for acc in load_IDs(todo_sras)]

	# create kallisto index file
	index_file = os.path.join(tmpdir, "index")
	if not os.path.isfile(index_file):
		logger.info("Starting Kallisto indexing")
		cmd1 = " ".join([kallisto, "index", "--index=" + index_file, "--make-unique", cds_file])
		p = subprocess.Popen(args=cmd1, shell=True, stdout=subprocess_log, stderr=subprocess_log)
		p.communicate()

	if RICH_AVAILABLE:
		command = Live(build_dashboard(), refresh_per_second=2)
	else:
		from contextlib import nullcontext
		command = nullcontext(None)

	with command as live_display:
		fk_thread = Thread(target=fasterqdump_kallisto_worker, args=(fasterqpigz_completed_accessions, kallisto_completed_accessions, index_file, sradir, readfile_status, logger, kallistodir, kallisto, cores,tmpdir, fasterq_dump, attempts, base_wait, failed_accessions_file,fasterqpigz_completed_accessions_file,kallisto_completed_accessions_file,subprocess_log))
		fk_thread.start()

		# keep active prefetch threads = batch dynamically
		active_prefetch_threads = []
		for accession in sra_accessions:
			# wait if batch prefetches are already running
			while len([t for t in active_prefetch_threads if t.is_alive()]) >= batch_size:
				live_display.update(build_dashboard())  #refresh while waiting
				time.sleep(5)  # check every 5 seconds

			# clean up finished threads
			active_prefetch_threads = [t for t in active_prefetch_threads if t.is_alive()]
			# start new prefetch thread for this accession
			t = Thread(target=parallel_prefetch, args=(prefetch_completed_accessions, accession, attempts,sradir, prefetch_command, minimum_sra_file_size_threshold, base_wait, failed_accessions_file, prefetch_completed_accessions_file, logger,subprocess_log))
			t.start()
			active_prefetch_threads.append(t)
			live_display.update(build_dashboard())  # refresh after starting new thread

		# wait for all remaining prefetch threads to finish
		for t in active_prefetch_threads:
			t.join()
			live_display.update(build_dashboard())  # refresh as each prefetch finishes

		# signal fasterq+kallisto consumer that no more accessions are coming
		kallisto_queue.put(barrier)
		fk_thread.join()  # wait for last Kallisto to finish before merge step
		live_display.update(build_dashboard())  # final refresh

	subprocess_log.close()#close the subprocess_log file handle since there are no subprocesses to be logged after this step

	# Merge the TPM, Counts of SRA samples into a single TPM, Counts file respectively
	tpmfile = os.path.join(kallistofinaldir, f'{orgname}_unfiltered.tpms.tsv')
	countsfile = os.path.join(kallistofinaldir, f'{orgname}_unfiltered.counts.tsv')

	if not (os.path.exists(tpmfile) and os.path.exists(countsfile)):
		logger.info(f'Merging TPM files and count files of all the samples')
		if RICH_AVAILABLE:
			sys.stdout.write(f'Merging TPM files and count files of all the samples\n')
			sys.stdout.flush()
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

	#code block to filter RNA-seq samples
	filtered_tpm_file = os.path.join(outdir, f"{orgname}.tpms.tsv")
	if not os.path.exists(filtered_tpm_file):
		logger.info("Filtering the TPM expression file")
		if RICH_AVAILABLE:
			sys.stdout.write(f'Filtering the TPM expression file\n')
			sys.stdout.flush()
		# --- run analysis of all data in folder/file --- #
		doc_file = os.path.join(tmpdir,f'{orgname}_qc.doc')
		valid_samples = []
		with open(doc_file, "w") as out:
			out.write("SampleName\tPercentageOfTop100\tPercentageOfTop500\tPercentageOfTop1000\n")
			TPM_data, genes = load_all_TPMs(tpmfile)
			count_data, genes = load_all_TPMs(countsfile)
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
		if RICH_AVAILABLE:
			sys.stdout.write(f"number of valid sample: {str(len(valid_samples))}\n")
			sys.stdout.flush()
			sys.stdout.write(f"number of invalid sample: {str(len(TPM_data.keys()) - len(valid_samples))}\n")
			sys.stdout.flush()

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
			if RICH_AVAILABLE:
				sys.stdout.write(f'WARNING: no valid samples in data set!\n')
				sys.stdout.flush()
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
			if RICH_AVAILABLE:
				sys.stdout.write(f'Xpression_collector pipeline completed successfully!!!\n')
				sys.stdout.flush()

	#optional removal of isoforms
	isoform_reduced_cds_file = outdir + f"{orgname}_repr.cds.fasta"
	repr_output_spec = os.path.join(outdir, f'{orgname}_repr.pep.fasta')
	repr_tpm_file = os.path.join(outdir, f'{orgname}_repr.tpms.tsv')
	if remove_isoforms == 'yes':
		if os.path.exists(isoform_reduced_cds_file) and os.path.exists(repr_output_spec) and os.path.exists(repr_tpm_file):
			pass
		else:
			logger.info("Removing alternative isoforms")
			if RICH_AVAILABLE:
				sys.stdout.write(f"Removing alternative isoforms\n")
				sys.stdout.flush()
			cds_dict = load_multiple_fasta_file(cds_file)
			repr_ids = isoform_clean(gff_file, cds_dict, isoform_reduced_cds_file, child_attribute, child_parent_linker)

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
			repr_tpm_filtered = keep_primary_transcript_exp(repr_ids, repr_tpm_file, isoform_reduced_cds_file, filtered_tpm_file)
			logger.info(f'Xpression_collector pipeline completed successfully!!!')
			if RICH_AVAILABLE:
				sys.stdout.write(f"Xpression_collector pipeline completed successfully!!!\n")
				sys.stdout.flush()

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
			logger.info(f"Merge of currently produced TPM file and user supplied TPM file successful: {merged_df.shape[0]} genes x {merged_df.shape[1]} columns")
			if RICH_AVAILABLE:
				sys.stdout.write(f"Merge of currently produced TPM file and user supplied TPM file successful: {merged_df.shape[0]} genes x {merged_df.shape[1]} columns")
				sys.stdout.flush()
			logger.info(f'Xpression_collector pipeline completed successfully!!!')
			if RICH_AVAILABLE:
				sys.stdout.write(f"Xpression_collector pipeline completed successfully!!!\n")
				sys.stdout.flush()
		except ValueError as e:
			logger.error(f"Aborting merge:\n{e}")
			if RICH_AVAILABLE:
				sys.stdout.write(f"Aborting merge:\n{e}")
				sys.stdout.flush()
	if merge_filtered_repr_tpm and remove_isoforms == 'yes':
		try:
			merged_df = merge_expression_tsvs(repr_tpm_filtered, merge_filtered_repr_tpm)
			merged_df.to_csv(merged_filtered_repr_tpm_file, sep='\t', index=False)
			logger.info(f"Merge of currently produced TPM file and user supplied TPM file successful: {merged_df.shape[0]} genes x {merged_df.shape[1]} columns")
			if RICH_AVAILABLE:
				sys.stdout.write(f"Merge of currently produced TPM file and user supplied TPM file successful: {merged_df.shape[0]} genes x {merged_df.shape[1]} columns\n")
				sys.stdout.flush()
			logger.info(f'Xpression_collector pipeline completed successfully!!!')
			if RICH_AVAILABLE:
				sys.stdout.write(f"Xpression_collector pipeline completed successfully!!!\n")
				sys.stdout.flush()
		except ValueError as e:
			logger.error(f"Aborting merge:\n{e}")
			if RICH_AVAILABLE:
				sys.stdout.write(f"Aborting merge:\n{e}\n")
				sys.stdout.flush()
	if clean_up=='yes':
		shutil.rmtree(tmpdir)

if '--sra' in sys.argv and '--cds' in sys.argv and '--out' in sys.argv:
	main(sys.argv)
else:
	sys.exit(__usage__)
