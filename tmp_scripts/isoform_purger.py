### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.21 ###

__usage__ = """
					python3 isoform_purger.py
					--in <INPUT_FILE>
					--out <OUTPUT_FOLDER>
					
					optional:
					--scorecut <SCORE_CUTOFF, INT>
					--simcut <SIMILARITY_CUTOFF,FLOAT>
					--lencut <LENGTH_CUTOFF,INT>
					--snvcut <SNV_CUTOFF,INT>
					
					bug reports and feature requests: b.pucker@tu-bs.de					
					"""

import os, glob, sys, subprocess
from operator import itemgetter

# --- end of imports --- #


def load_fasta( fasta_file ):
	"""! @brief load FASTA alignment into dictionary	"""
	
	sequences = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		if " " in header:
			header = header.split(' ')[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					if " " in header:
						header = header.split(' ')[0]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def load_hits_per_bait( blast_result_file, scorecut, simcut, lencut ):
	"""! @brief load BLAST hits per bait
	@return dictionary with baits as keys and lists of hits as values
	"""
	
	hits = {}
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if parts[0] != parts[1]:
				if float( parts[-1] ) > scorecut:	#BLAST hits are filtered with user defined criteria
					if float( parts[2] ) > simcut:
						if int( parts[3] ) > lencut:
							try:
								if parts[-1] not in hits[ parts[0] ]:
									hits[ parts[0] ].append( parts[1] )	#add to existing dictionary entry
							except KeyError:
								hits.update( { parts[0]: [ parts[1] ] } )	#generate new entry in dictionary
			line = f.readline()
	return hits


def load_aln_similarity( alignment, snv_cutoff ):
	"""! @brief load alignment similarity """
	
	results = []
	seqIDs = list( alignment.keys() )
	for idx1, ID1 in enumerate( seqIDs ):
		for idx2, ID2 in enumerate( seqIDs[1:] ):
			if idx2 >= idx1:
				aln_seq1 = alignment[ ID1 ]
				aln_seq2 = alignment[ ID2 ]
				snv_counter = 0
				indel_counter = 0
				for i1, nt in enumerate( aln_seq1 ):
					if nt != "-" and aln_seq2[ i1 ] != "-":
						if nt != aln_seq2[ i1 ]:
							snv_counter += 1
					else:
						indel_counter += 1
				if snv_counter <= snv_cutoff:
					results.append( { 'ID1': ID1, 'ID2': ID2, 'snvs': snv_counter, 'indels': indel_counter } )
	return results


def identify_isoforms( alignment_results, snv_cutoff ):
	"""! @brief identify isoforms """
	
	isoform_groups = []
	for group in alignment_results:
		tmp_groups = []
		for entry in group:
			if entry['snvs'] <= snv_cutoff:
				status = False	#genes not included yet
				for idx, g in enumerate( tmp_groups ):
					if entry['ID1'] in g:
						tmp_groups[ idx ].append( entry['ID2'] )
						status = True
						break
					elif entry['ID2'] in g:
						tmp_groups[ idx ].append( entry['ID1'] )
						status = True
						break
				if not status:
					tmp_groups.append( [entry['ID1'], entry['ID2'] ] )
		for g in tmp_groups:	#go through all tmp groups to break up too nested structure
			if len(g) > 0:	#only collect non empty lists
				isoform_groups.append( g )
	return isoform_groups


def identify_repr_isoform_per_group( isoforms, seqs ):
	"""! @brief identify representative isoform per group """
	
	repr_seq_collection = {}
	for group in isoforms:
		print( group )
		tmp_seqs = []
		for ID in group:
			tmp_seqs.append( { 'ID': ID, 'seq': seqs[ ID ], 'len': len( seqs[ ID ] ) } )
		sorted_list = sorted( tmp_seqs, key=itemgetter('len') )
		repr_seq_collection.update( { sorted_list[-1]['ID']: sorted_list[-1]['seq'] } )
	return repr_seq_collection


def main( arguments ):
	"""! @brief running everything """
	
	input_file = arguments[ arguments.index('--in')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	
	if '--scorecut' in arguments:
		scorecut = int( arguments[ arguments.index('--scorecut')+1 ] )
	else:
		scorecut = 100
	
	if '--simcut' in arguments:
		simcut = float( arguments[ arguments.index('--simcut')+1 ] )
	else:
		simcut = 99.0
	
	if '--lencut' in arguments:
		lencut = int( arguments[ arguments.index('--lencut')+1 ] )
	else:
		lencut = 100
	
	if '--snvcut' in arguments:
		snv_cutoff = int( arguments[ arguments.index('--snvcut')+1 ] )
	else:
		snv_cutoff = 5
	
	if output_folder[-1] != "/":
		output_folder += "/"
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	# --- clean file --- #
	clean_file = output_folder + "clean.fasta"
	seqs = load_fasta( input_file )
	with open( clean_file, "w" ) as out:
		for key in list( seqs.keys() ):
			out.write( '>' + key + "\n" + seqs[key] + "\n" )
	
	
	# --- run BLAST vs self --- #
	dbname = output_folder + "blastdb"
	blast_result_file = output_folder + "blast_results.txt"
	if not os.path.isfile( blast_result_file ):
		cmd = "makeblastdb -in " + clean_file + " -out " + dbname + " -dbtype nucl"
		p = subprocess.Popen( args= cmd, shell=True )
		p.communicate()
		
		cmd = "blastn -query " + clean_file + " -db " + dbname + " -out " + blast_result_file + " -outfmt 6 -evalue 0.00000000001 -num_threads 8"
		p = subprocess.Popen( args= cmd, shell=True )
		p.communicate()
	
	
	# --- identify good BLAST hits around sequences --- #
	blast_hits = load_hits_per_bait( blast_result_file, scorecut, simcut, lencut )
	black_list = {}
	groups_to_analyze = []
	for key in list( blast_hits.keys() ):
		try:
			black_list[ key ]	#check if sequence ID is already assigned to group
		except KeyError:
			group = blast_hits[ key ] + [ key ]	#combine key with all good BLAST hits
			if len( group ) > 1:
				groups_to_analyze.append( group )
			for each in group:
				black_list.update( { each: None } )
	print( "number of groups to analyze: " + str( len( groups_to_analyze ) ) )
	single_fasta_folder = output_folder + "single_fasta_input/"	
	if not os.path.exists( single_fasta_folder ):
		os.makedirs( single_fasta_folder )
	for idx, group in enumerate( groups_to_analyze ):
		output_file_name = single_fasta_folder + str( idx+1 ) + ".fasta"
		if not os.path.isfile( output_file_name ):
			with open( output_file_name, "w" ) as out:
				for ID in group:
					out.write( '>' + ID + "\n" + seqs[ ID ] + "\n" )
	
	# --- construct global alignment --- #
	fasta_files = glob.glob( single_fasta_folder + "*.fasta" )
	for fasta in fasta_files:
		if not os.path.isfile( fasta + ".aln" ):
			cmd = "mafft " + fasta + " > " + fasta + ".aln 2> " + fasta + ".doc.txt"
			p = subprocess.Popen( args= cmd, shell=True )
			p.communicate()
	
	# --- calculate similarity matrix --- #
	aln_fasta_files = glob.glob( single_fasta_folder + "*.fasta.aln" )
	print( "number of aligned FASTA files: " + str( len( aln_fasta_files ) ) )
	alignment_results = []
	for aln_fasta in aln_fasta_files:
		alignment = load_fasta( aln_fasta )
		sim_per_aln = load_aln_similarity( alignment, snv_cutoff )
		alignment_results.append( sim_per_aln )
		
	
	# --- classify sequences as isoforms/paralogs --- #
	isoforms = identify_isoforms( alignment_results, snv_cutoff )
	print( "number of isoform groups: " + str( len( isoforms ) ) )
	repr_isoforms = identify_repr_isoform_per_group( isoforms, seqs )
	blacklist = {}
	for ID in [ x for sublist in isoforms for x in sublist ]:
		blacklist.update( { ID: None } )
	
	# --- write isoform-reduced output file --- #
	output_file = output_folder + "isoforms_reduced.fasta"
	with open( output_file, "w" ) as out:
		for key in list( seqs.keys() ):
			try:
				blacklist[ key ]
			except KeyError:
				out.write( '>' + key + "\n" + seqs[ key ] + "\n" )
		for key in list( repr_isoforms.keys() ):
			out.write( '>' + key + "\n" + seqs[ key ] + "\n" )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
