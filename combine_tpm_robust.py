### Shakunthala Natarajan ###

# !/usr/bin/env python3
import glob
import sys
import os
import argparse
from collections import OrderedDict


def read_tpm_file(filename):
	"""Read TPM file into a dictionary with gene names as keys"""
	with open(filename, 'r') as f:
		lines = f.readlines()

	# First row: date + SRA accessions
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


def main():
	# Set up argument parser
	parser = argparse.ArgumentParser(
		description='Combine multiple TPM files into a single file',
		formatter_class=argparse.RawDescriptionHelpFormatter,
		epilog='''
Examples:
  %(prog)s -i /path/to/tpm_files -o /path/to/output
  %(prog)s --input ./data --output ./results
        '''
	)
	parser.add_argument('-i', '--input', required=True,
						help='Path to folder containing TPM.txt files')
	parser.add_argument('-o', '--output', required=True,
						help='Path to folder where combined file will be saved')

	args = parser.parse_args()

	# Validate input folder
	if not os.path.isdir(args.input):
		print(f"Error: Input folder '{args.input}' does not exist!")
		sys.exit(1)

	# Create output folder if it doesn't exist
	if not os.path.exists(args.output):
		print(f"Creating output folder: {args.output}")
		os.makedirs(args.output)

	# Get all TPM files from input folder
	search_pattern = os.path.join(args.input, '*TPM*.txt')
	tpm_files = sorted(glob.glob(search_pattern))

	if not tpm_files:
		print(f"Error: No TPM files found in '{args.input}'!")
		print(f"Searched for pattern: {search_pattern}")
		sys.exit(1)

	print(f"Input folder: {args.input}")
	print(f"Output folder: {args.output}")
	print(f"Found {len(tpm_files)} TPM files")
	print(
		f"Processing: {', '.join([os.path.basename(f) for f in tpm_files[:3]])}{'...' if len(tpm_files) > 3 else ''}\n")

	# Process first file to establish gene order
	print(f"[1/{len(tpm_files)}] Reading {os.path.basename(tpm_files[0])} (establishing gene order)...")
	date1, sra_list, master_gene_data = read_tpm_file(tpm_files[0])
	all_sra_accessions = sra_list.copy()
	gene_order = list(master_gene_data.keys())  # Preserve gene order from first file

	print(f"  - Found {len(gene_order)} genes")
	print(f"  - Found {len(sra_list)} SRA samples")

	# Process remaining files
	for idx, tpm_file in enumerate(tpm_files[1:], start=2):
		print(f"[{idx}/{len(tpm_files)}] Processing {os.path.basename(tpm_file)}...")

		date_info, sra_list, gene_data = read_tpm_file(tpm_file)

		# Validation checks
		if set(gene_data.keys()) != set(gene_order):
			missing_in_new = set(gene_order) - set(gene_data.keys())
			extra_in_new = set(gene_data.keys()) - set(gene_order)

			if missing_in_new:
				print(f"  WARNING: {len(missing_in_new)} genes missing in {os.path.basename(tpm_file)}")
				print(f"    First 5 missing: {list(missing_in_new)[:5]}")
			if extra_in_new:
				print(f"  WARNING: {len(extra_in_new)} extra genes in {os.path.basename(tpm_file)}")
				print(f"    First 5 extra: {list(extra_in_new)[:5]}")

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
		print(f"  - Added {len(sra_list)} SRA samples")

	# Write combined output
	output_file = os.path.join(args.output, 'combined_TPM.tsv')
	print(f"\nWriting combined file to {output_file}...")

	with open(output_file, 'w') as out:
		# Write header row
		out.write('gene\t' + '\t'.join(all_sra_accessions) + '\n')

		# Write gene expression data (in original gene order)
		for gene in gene_order:
			out.write(gene + '\t' + '\t'.join(master_gene_data[gene]) + '\n')

	print(f"\n✓ Success!")
	print(f"  - Total genes: {len(gene_order)}")
	print(f"  - Total SRA samples: {len(all_sra_accessions)}")
	print(f"  - Output file: {output_file}")


if __name__ == "__main__":
	main()
