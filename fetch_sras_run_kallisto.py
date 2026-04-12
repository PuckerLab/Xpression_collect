__version__=0.1
__usage__="""
			python3 fetch_sras_run_kallisto.py
			--sra <full path to folder of count table files or a single count table file with sra ids as column names>
			--kallisto <full path to kallisto pipeline3 script including the script name>
			--kallisto_merge <full path to kallisto merge script including the script name>
			--out <full path to output directory>
			Additional instructions:
			--threads <no. of cores for fasterq-dump per accession and kallisto run per batch>
			--batch_size <no. of sras to be processed in one batch>
			"""

### --- start imports --- ###
import os,sys,glob
import gzip
import subprocess
import copy
import shutil
### --- end of imports --- ###
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

def main(arguments):
	todo_sras = arguments[arguments.index('--sra')+1]
	kallisto_run = arguments[arguments.index('--kallisto')+1] #full path to kallisto pipeline3 script including the script name
	kallisto_merge = arguments[arguments.index('--kallisto_merge')+1]#full path to kallisto merge script including the script name
	cds_file = arguments[arguments.index('--cds')+1]
	if '--threads' in arguments:
		cores = int(arguments[arguments.index('--threads')+1])
	else:
		cores = 4
	if '--batch_size' in arguments:
		batch_size=int(arguments[arguments.index('--batch_size')+1])
	else:
		batch_size = 500
	outdir = arguments[arguments.index('--out') + 1]
	if outdir[-1] != "/":
		outdir += "/"
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	sradir=os.path.join(outdir,'SRA')
	if not os.path.exists(sradir):
		os.makedirs(sradir)
	tmpdir = os.path.join(outdir, 'TMP')
	if not os.path.exists(tmpdir):
		os.makedirs(tmpdir)
	kallistofinaldir = os.path.join(outdir, 'Kallisto_run_final')
	if not os.path.exists(kallistofinaldir):
		os.makedirs(kallistofinaldir)
	sra_accessions = load_IDs(todo_sras)
	batch_counter = 1

	while sra_accessions:
		batch = sra_accessions[:batch_size]  # take first N
		sra_accessions = sra_accessions[batch_size:]  # consume N

		print(f"Processing batch {batch_counter} ({len(batch)} SRAs)")

		# prefetch
		for accession in batch:
			cmd = f"prefetch --max-size 200G {accession} -O {sradir}"
			subprocess.run(cmd, shell=True)
			# fasterq-dump + pigz
			# Create subfolder for this accession
			acc_dir = os.path.join(sradir, accession)
			os.makedirs(acc_dir, exist_ok=True)
			prefetched_file = os.path.join(acc_dir,f"{accession}.sra")
			# fasterq-dump
			cmd = f"fasterq-dump --split-3 --outdir {acc_dir} --skip-technical --threads {cores} {prefetched_file}"
			subprocess.run(cmd, shell=True)

			# pigz (generalized)
			fq_patterns = [
				os.path.join(acc_dir, f"{accession}*.fastq"),
				os.path.join(acc_dir, f"{accession}*.fq"),
			]
			fq_files = []
			for pattern in fq_patterns:
				fq_files.extend(glob.glob(pattern))
			if fq_files:
				cmd = f"pigz -p {cores} " + " ".join(fq_files)
				subprocess.run(cmd, shell=True)
			# After pigz, check what was created
			print(f"Files in {acc_dir}:")
			for f in os.listdir(acc_dir):
				print(f"  {f}")
			if os.path.exists(prefetched_file):
				os.remove(prefetched_file)

		# Create Kallisto dir
		kallistodir = os.path.join(outdir, f'Kallisto_run_{batch_counter}')
		os.makedirs(kallistodir, exist_ok=True)

		# Run Kallisto
		cmd = f"python3 {kallisto_run} --cds {cds_file} --reads {sradir} --tmp {tmpdir} --out {kallistodir} --cpus {cores}"
		print(cmd)
		subprocess.run(cmd, shell=True)


		# Merge
		tpmfile = os.path.join(kallistofinaldir, f"TPM_{batch_counter}.txt")
		countsfile = os.path.join(kallistofinaldir, f"Counts_{batch_counter}.txt")

		cmd = f"python3 {kallisto_merge} --in {kallistodir} --tpms {tpmfile} --counts {countsfile}"
		subprocess.run(cmd, shell=True)

		# Clean SRA, TMP folders for next batch
		shutil.rmtree(sradir)
		os.makedirs(sradir)
		for subdir in glob.glob(os.path.join(tmpdir, "SRR*")):
			shutil.rmtree(subdir)
		batch_counter += 1


if '--sra' in sys.argv and '--kallisto' in sys.argv and '--kallisto_merge' in sys.argv and '--cds' in sys.argv and '--out' in sys.argv:
	main(sys.argv)
else:
	sys.exit(__usage__)
