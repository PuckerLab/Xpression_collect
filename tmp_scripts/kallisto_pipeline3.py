### Boas Pucker ###
### Milan Borchert ###
### b.pucker@tu-bs.de ###
### v0.37 ###

__usage__ = """
					python3 kallisto_pipeline3.py
					--cds <CDS_FILE>
					--reads <FASTQ_FILE_FOLDER>
					--tmp <TMP_FOLDER>
					--out <FINAL_OUTPUT_FOLDER>
					
					optional:
					--kallisto <FULL_PATH_TO_KALLISTO>[kallisto]
					--cpus <NUMBER_OF_CPUS_TO_USE>[10]
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, sys, glob, time, subprocess


# --- end of imports --- #

def get_data_for_jobs_to_run( read_file_folders, final_output_folder, index_file, tmp_cluster_folder ):
	"""! @brief collect all infos to run jobs """
	
	jobs_to_do = []
	for folder in read_file_folders:
		ID = folder.split('/')[-1]
		status = True
		
		# --- get read file --- #
		PE_status = True
		SRA = False
		read_file1 = folder + "/" + ID + "_R1_001.fastq.gz"
		if not os.path.isfile( read_file1 ):
			#print "ERROR: file missing - " + read_file1
			PE_status = False
			if not PE_status:
				read_file1 = folder + "/" + ID + "_pass_1.fastq.gz"
				if os.path.isfile( read_file1 ):
					PE_status = True
					SRA = True
					read_file2 = folder + "/" + ID + "_pass_2.fastq.gz"
					if not os.path.isfile( read_file2 ):
						#print "ERROR: file missing - " + read_file2
						PE_status = False
				else:
					read_file1 = folder + "/" + ID + "_1.fastq.gz"
					if os.path.isfile( read_file1 ):
						PE_status = True
						SRA = True
						read_file2 = folder + "/" + ID + "_2.fastq.gz"
						if not os.path.isfile( read_file2 ):
							#print "ERROR: file missing - " + read_file2
							PE_status = False
					else:
						read_file1 = folder + "/" + ID + "_R1.fq.gz"
						if os.path.isfile( read_file1 ):
							PE_status = True
							SRA = True
							read_file2 = folder + "/" + ID + "_R2.fq.gz"
							if not os.path.isfile( read_file2 ):
								#print "ERROR: file missing - " + read_file2
								PE_status = False
		if not SRA:
			read_file2 = folder + "/" + ID + "_R2_001.fastq.gz"
			if not os.path.isfile( read_file2 ):
				#print "ERROR: file missing - " + read_file2
				PE_status = False
		if not PE_status:
			read_file1 = folder + "/" + ID + "_R1_001.fastq.gz"
			if not os.path.isfile( read_file1 ):
				read_file1 = folder + "/" + ID + ".fastq.gz"
				if not os.path.isfile( read_file1 ):
					read_file1 = folder + "/" + ID + "_1.fastq.gz"
					if not os.path.isfile( read_file1 ):
						read_file1 = folder + "/" + ID + "_R1.fq.gz"
						if not os.path.isfile( read_file1 ):
							read_file1 = folder + "/" + ID + "_pass.fastq.gz"
							if not os.path.isfile( read_file1 ):
								status = False
			read_file2 = False
		
		# --- get reference for quantification --- #
		output_dir = tmp_cluster_folder + ID + "/"
		if not os.path.exists( output_dir ):
			os.makedirs( output_dir ) 
		tmp_result_file = output_dir + "abundance.tsv"
		try:
		    timestr = time.strftime("%Y_%m_%d_")
		except ModuleNotFoundError:
		    timestr = ""
		final_result_file = final_output_folder + timestr + ID + ".tsv"
		if os.path.isfile( final_result_file ):
			status = False
		if status:
			jobs_to_do.append( { 'r1': read_file1, 'r2': read_file2, 'out': output_dir, 'index': index_file, 'tmp': tmp_result_file, 'fin': final_result_file, "ID": ID } )
	return jobs_to_do


def job_executer( jobs_to_run, kallisto, threads ):
	"""! @brief run all jobs in list """
	
	for idx, job in enumerate( jobs_to_run ):
		sys.stdout.write( "running job " + str( idx+1 ) + "/" + str( len( jobs_to_run ) ) + " - " + job["ID"] + "\n" )
		sys.stdout.flush()
		
		if job['r2']:
			cmd2 = " ".join( [ kallisto, "quant", "--index="+job['index'], "--output-dir="+job['out'], "--threads "+str( threads ), job['r1'], job['r2'] ] )
		else:
			cmd2 = " ".join( [ kallisto, "quant", "--index="+job['index'], "--single -l 200 -s 100", "--output-dir="+job['out'], "--threads "+str( threads ), job['r1'] ] )
		p = subprocess.Popen( args= cmd2, shell=True )
		p.communicate()

		p = subprocess.Popen( args= "cp " + job["tmp"] + " " + job["fin"] , shell=True )
		p.communicate()


def main( arguments ):
	"""! @brief run everything """

	cds_file = arguments[ arguments.index( '--cds' )+1 ]
	read_file_folder = arguments[ arguments.index( '--reads' )+1 ]
	tmp_cluster_folder = arguments[ arguments.index( '--tmp' )+1 ]
	final_output_folder = arguments[ arguments.index( '--out' )+1 ]
	
	if '--kallisto' in arguments:
		kallisto = arguments[ arguments.index( '--kallisto' )+1 ]
	else:
		kallisto = "kallisto"
	
	if '--cpus' in arguments:
		try:
			threads = int( arguments[ arguments.index( '--cpus' )+1 ] )
		except:
			threads = 10
	else:
		threads = 10
	
	if final_output_folder[-1] != "/":
		final_output_folder += "/"
	if tmp_cluster_folder[-1] != "/":
		tmp_cluster_folder += "/"
	
	if not os.path.exists( tmp_cluster_folder ):
		os.makedirs( tmp_cluster_folder )
	if not os.path.exists( final_output_folder ):
		os.makedirs( final_output_folder )
	
	# --- load data --- #
	single_read_file_folders = [ x[0] for x in os.walk( read_file_folder ) ][1:]
	sys.stdout.write( "number of FASTQ file folders detected: " + str( len( single_read_file_folders ) ) + "\n" )
	sys.stdout.flush()
	
	# --- prepare jobs to run --- #
	index_file = tmp_cluster_folder + "index"
	jobs_to_run = get_data_for_jobs_to_run( single_read_file_folders, final_output_folder, index_file, tmp_cluster_folder )
	sys.stdout.write( "number of jobs to run: " + str( len( jobs_to_run ) ) + "\n" )
	sys.stdout.flush()
	
	# --- generate index --- #
	if not os.path.isfile( index_file ):
		cmd1 = " ".join( [ kallisto, "index", "--index="+index_file, "--make-unique", cds_file ] )
		p = subprocess.Popen( args= cmd1, shell=True )
		p.communicate()
	
	# --- run jobs --- #
	job_executer( jobs_to_run, kallisto, threads )


if '--cds' in sys.argv and '--reads' in sys.argv and '--tmp' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
