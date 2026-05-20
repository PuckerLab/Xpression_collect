### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.4 ###

__usage__ = """
						python3 filter_RNAseq_samples.py
						--tpms <TPM_FILE>
						--counts <COUNT_FILE>
						--out <OUTPUT_FILE>
						--min <MIN_PERCENT_EXPRESSION_ON_TOP100>
						--max <MAX_PERCENT_EXPRESSION_ON_TOP100>
						
						optional:
						--mincounts <MIN_READ_NUMBER>[1000000]
						--black <ID_BLACK_LIST>
						"""

import os, sys, glob
import matplotlib.pyplot as plt

# --- end of imports --- #


def load_all_TPMs( exp_file ):
	"""! @brief load all values from given TPM file """
	
	data = {}
	genes = []
	with open( exp_file, "r" ) as f:
		headers = f.readline().strip()
		if "\t" in headers:
			headers = headers.split('\t')
		else:
			headers = [ headers ]
		if headers[0] == "gene":
			headers = headers[1:]
		elif headers[0] == exp_file.split('/')[-1][:10]:
			headers = headers[1:]
		for header in headers:
			data.update( { header: [] } )
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			genes.append( parts[0] )
			for idx, val in enumerate( parts[1:] ):
				data[ headers[ idx ] ].append( float( val ) )
			line = f.readline()
	return data, genes


def  load_black_IDs( black_list_file ):
	"""! @brief load IDs from given black list """
	
	black_list = {}
	with open( black_list_file, "r" ) as f:
		line = f.readline()
		while line:
			if len( line ) > 3:
				black_list.update( { line.strip(): None } )
			line = f.readline()	
	return black_list


def main( arguments ):
	"""! @brief run everything """
	
	tpm_file = arguments[ arguments.index('--tpms')+1 ]
	count_file = arguments[ arguments.index('--counts')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	
	if '--min' in arguments:
		min_cutoff = int( arguments[ arguments.index('--min')+1 ] )
	else:
		min_cutoff = 10	#value in percent
	if '--max' in arguments:
		max_cutoff = int( arguments[ arguments.index('--max')+1 ] )
	else:
		max_cutoff = 80
	
	if '--mincounts' in arguments:
		min_counts = int( arguments[ arguments.index('--mincounts')+1 ] )
	else:
		min_counts = 1000000
	
	if '--black' in arguments:
		black_list_file = arguments[ arguments.index('--black')+1 ]
		black_list = load_black_IDs( black_list_file )
	else:
		black_list = {}
	
	# --- run analysis of all data in folder/file --- #
	doc_file = output_file + ".doc"
	valid_samples = []
	with open( doc_file, "w" ) as out:
		out.write( "SampleName\tPercentageOfTop100\tPercentageOfTop500\tPercentageOfTop1000\n" )
		TPM_data, genes = load_all_TPMs( tpm_file )
		count_data, genes = load_all_TPMs( count_file )
		for key in sorted( list( TPM_data.keys() ) ):
			new_line = [ key ]
			selection = sorted( TPM_data[ key ] )
			counts = sum( count_data[ key ] )	#calculate counts per library
			if counts >= min_counts:	#check for sufficient library size
				try:	#check for ID presence on black list
					black_list[ key ]
					new_line.append( "ID on black list" )
					out.write( "\t".join( list( map( str, new_line ) ) ) + "\n" )
				except KeyError:
					try:
						val = 100.0 * sum( selection[-100:] ) / sum( selection )
					except ZeroDivisionError:
						val = 0
					new_line.append( val )
					if min_cutoff < val < max_cutoff:
						valid_samples.append( key )
					if len( selection ) > 500 and val > 0:
						new_line.append( 100.0 * sum( selection[-500:] ) / sum( selection ) )
					else:
						new_line.append( "n/a" )
					if len( selection ) > 1000 and val > 0:
						new_line.append( 100.0 * sum( selection[-1000:] ) / sum( selection ) )
					else:
						new_line.append( "n/a" )
					out.write( "\t".join( list( map( str, new_line ) ) ) + "\n" )
			else:
				new_line.append( "insufficient counts: " + str( counts ) )
				out.write( "\t".join( list( map( str, new_line ) ) ) + "\n" )
	
	print ( "number of valid sample: " + str( len( valid_samples ) ) )
	print ( "number of invalid sample: " + str( len( TPM_data.keys() ) - len( valid_samples ) ) )
	
	# --- generate output file --- #
	if len( valid_samples ) > 0:
		with open( output_file, "w" ) as out:
			out.write( "gene\t" + "\t".join( valid_samples )+ "\n" )
			for idx, gene in enumerate( genes ):
				new_line = [ gene ]
				for sample in valid_samples:
					new_line.append( TPM_data[ sample ][ idx ] )
				out.write( "\t".join( list( map( str, new_line ) ) ) + "\n" )
	else:
		print ( "WARNING: no valid samples in data set!" )
	# --- generate figure --- #
	fig_file = output_file + ".pdf"
	values = []
	with open( doc_file, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				values.append( float( parts[1] ) )
			except ValueError:
				pass
			line = f.readline()
	
	values = [ x for x in values if str(x) != 'nan' ]
	
	fig, ax = plt.subplots()
	
	ax.hist( values, bins=100, color="green" )
	ax.set_xlabel( "Percentage of expression on top100 genes" )
	ax.set_ylabel( "Number of analyzed samples" )
	
	fig.savefig( fig_file )


if '--tpms' in sys.argv and '--counts' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
