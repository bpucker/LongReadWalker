### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.22 ###
### Long Read Walker (LRW) ###

#USE FILTLONG to reduce coverage to 5x ???

__cite__ = """Pucker, 2021: https://github.com/bpucker/LongReadWalker """

__usage__ = """
					python3 LRW.py
					--reads <READS_IN_FASTA_OR_FASTQ>
					--seed <SEED_SEQUENCE_IN_FASTA>
					--out <OUTPUT_FOLDER>
					
					optional:
					--rounds <NUMBER_OF_EXTENSION_ROUNDS>[5]
					--block <BLOCK_SIZE>[2000]
					--direction <DIRECTION_TO_EXTEND (up|down)>[down]
					--simcut <MIN_BLAST_HIT_SIMILARITY>[80]
					--lencut <MIN_BLAST_HIT_LENGTH>[500]
					--threads <NUMBER_OF_THREADS_FOR_BLAST>[8]
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import sys, os, re, gzip
from operator import itemgetter

# --- end of imports --- #



def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
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


def revcomp( seq ):
	"""! @brief construct reverse complement of sequence """
	
	new_seq = []
	
	bases = { 'a':'t', 't':'a', 'c':'g', 'g':'c' }
	for nt in seq.lower():
		try:
			new_seq.append( bases[nt] )
		except:
			new_seq.append( 'n' )
	return ''.join( new_seq[::-1] ).upper()


def merge_HSPs( input_hits, block_size ):
	"""! @brief merge all HSPs per subject """
	
	# --- filter all available HSPs --- #
	final = []
	for hsp in sorted( input_hits, key=itemgetter('score') )[::-1]:
		if len( final ) == 0:
			final.append( hsp )
		else:
			if hsp['orientation'] == final[0]['orientation']:	#exclude HSPs with different orientation
				# --- exclude query overlap with any of the existing hsps --- #
				exclude = False
				for each in final:
					if hsp['qstart'] < each['qend']:
						if hsp['qend'] > each['qstart']:
							exclude = True	#exclude this hsp
					if min( [ hsp['sstart'], hsp['send'] ] ) < max( [ each['sstart'], each['send'] ] ):
						if max( [ hsp['sstart'], hsp['send'] ] ) > min( [ each['sstart'], each['send'] ] ):
							exclude = True
					avg_ref_pos = ( final[0]['sstart']+final[0]['send'] ) / 2.0
					avg_pos = ( hsp['sstart'] + hsp['send'] ) / 2.0
					#check distance between center of best HSP and current HSP; use block size minus HSP sizes as cutoff
					if abs( avg_ref_pos - avg_pos ) > ( 1.2 * block_size - 0.5 * abs( hsp['send'] - hsp['sstart'] ) - 0.5 * abs( final[0]['sstart'] - final[0]['send'] ) ):
						exclude = True	#exclude HPS which are too far apart
				if not exclude:
					final.append( hsp )
	
	# --- prepare final infos --- #
	qstart = final[0]['qstart']
	qend = final[0]['qend']
	sstart = final[0]['sstart']
	send = final[0]['send']
	scores = []
	sims = []
	lens = []
	for each in final:
		if each['qstart'] < qstart:
			qstart = 0 + each['qstart']
		if each['qend'] > qend:
			qend = 0 + each['qstart']
		if each['orientation']:	#hits on forward strand
			if each['sstart'] < sstart:
				sstart = 0 + each['sstart']
			if each['send'] > send:
				send = 0 + each['send']
		else:	#hits on reverse strand
			if each['sstart'] > sstart:
				sstart = 0 + each['sstart']
			if each['send'] < send:
				send = 0 + each['send']
		scores.append( each['score'] )
		sims.append( each['sim']*each['len'] )
		lens.append( each['len'] )
	
	final_info = { 	'ID': final[0]['ID'],
							'qstart': qstart,
							'qend': qend,
							'sstart': sstart,
							'send': send,
							'score': sum( scores ),
							'sim': sum( sims ) / sum( lens ),
							'len': sum( lens ),
							'orientation': final[0]['orientation']
						}
	return final_info


def analyze_BLAST_results( blast_result_file, block_size, all_reads, direction, sim_cut, len_cut, read_black_list ):
	"""! @brief find best BLAST hit which allows to extend sequence """
	
	best_hits = {}
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if parts[0] != parts[1]:
				if int( parts[8] ) < int( parts[9] ):
					orientation = True
				else:
					orientation = False
				try:
					best_hits[ parts[1] ].append( { 	'ID': parts[1],
																		'qstart': int( parts[6] ),
																		'qend': int( parts[7] ),
																		'sstart': int( parts[8] ),
																		'send': int( parts[9] ),
																		'score': float( parts[-1] ),
																		'sim': float( parts[2] ),
																		'len': int( parts[3] ),
																		'orientation': orientation
																	} )
				except KeyError:
					best_hits.update( { parts[1]: [ { 	'ID': parts[1],
																			'qstart': int( parts[6] ),
																			'qend': int( parts[7] ),
																			'sstart': int( parts[8] ),
																			'send': int( parts[9] ),
																			'score': float( parts[-1] ),
																			'sim': float( parts[2] ),
																			'len': int( parts[3] ),
																			'orientation': orientation
																		} ] } )
			line = f.readline()
	
	# --- merging all HSPs per subject --- #
	merged_best_hits = []
	for key in list( best_hits.keys() ):
		merged = merge_HSPs( best_hits[ key ], block_size )
		if merged['sim'] > sim_cut:
			if merged['len'] > len_cut:
				try:
					read_black_list[ merged['ID'] ]
				except KeyError:
					merged_best_hits.append( merged )
	
	# --- sort BLAST hits and return best hit --- #
	sorted_hits = sorted( merged_best_hits, key=itemgetter('score') )[::-1]	#sorted by decreasing score
	for i in range( 10 ):
		try:
			hit = sorted_hits[ i ]
			sys.stdout.write( hit['ID'][:10] + " - " + str( hit['sim'] ) + " - " + str( hit['len'] ) + " - " + str( hit['sstart'] ) + " - " + str( hit['send'] ) + " - " + str( hit['score'] ) +"\n" )
		except:
			pass
	sys.stdout.write( "\n" )
	sys.stdout.flush()
	for hit in sorted_hits:
		if direction == "down":	#continue extension towards downstream
			if hit['sstart'] < hit['send']:	#forward hit
				if len( all_reads[ hit['ID'] ] ) - hit['send'] > block_size:
					return { 'ID': hit['ID'], 'start': len( all_reads[ hit['ID'] ] )-block_size-1, 'end': len( all_reads[ hit['ID'] ] )-1, 'orientation': True }, { 'ID': hit['ID'], 'start': hit['send'], 'end': len( all_reads[ hit['ID'] ] )-1, 'orientation': True }
			else:	#reverse hit
				if hit['send'] > block_size:
					return { 'ID': hit['ID'], 'start': 1, 'end': block_size, 'orientation': False }, { 'ID': hit['ID'], 'start': 1, 'end': hit['send'], 'orientation': False }
			#return: ID, start, end, orientation (T=keep like this; F=revcomp)
		else:	#continue extension towards upstream
			if hit['sstart'] < hit['send']:	#forward hit
				if hit['sstart'] > block_size:
					return { 'ID': hit['ID'], 'start': 1, 'end': block_size, 'orientation': True }, { 'ID': hit['ID'], 'start': 1, 'end': hit['sstart'], 'orientation': True }
			else:	#reverse hit
				if len( all_reads[ hit['ID'] ] )-hit['send'] > block_size:
					return { 'ID': hit['ID'], 'start': len( all_reads[ hit['ID'] ] )-block_size-1, 'end': len( all_reads[ hit['ID'] ] )-1, 'orientation': False }, { 'ID': hit['ID'], 'start': hit['send'], 'end': len( all_reads[ hit['ID'] ] )-1, 'orientation': False }
	sys.exit( "ERROR: no suitable read for extension detected" )


def get_read_file_type( read_file ):
	"""! @brief get read file type """
	
	if read_file.split('.')[-1] in [ "gz", "gzip", "GZ", "GZIP" ]:	#compressed read input file
		with gzip.open( read_file, "rb" ) as f:
			line = f.readline().decode('ascii')
			if line[0] == ">":
				return "fasta"
			elif line[0] == "@":
				return "fastq"
			else:
				print( line )
				sys.exit( "ERROR: read input file type not recognized" )
	else:
		with open( read_file, "r" ) as f:
			line = f.readline()
			if line[0] == ">":
				return "fasta"
			elif line[0] == "@":
				return "fastq"
			else:
				sys.exit( "ERROR: read input file type not recognized" )


def load_sequences_fastq( read_file ):
	"""! @brief load sequences from FASTQ file """
	
	sequences = {}
	if read_file.split('.')[-1] in [ "gz", "gzip", "GZ", "GZIP" ]:	#compressed read input file
		with gzip.open( read_file, "rb" ) as f:
			line = f.readline().decode('ascii')
			while line:
				header = line.strip()
				if " " in header:
					header = header.split(' ')[0]
				seq = f.readline().strip()
				sequences.update( { header: seq } )
				f.readline()	#useless row
				f.readline()	#quality row
				line = f.readline().decode('ascii')
	else:	#uncompressed read input file
		with open( read_file, "r" ) as f:
			line = f.readline()
			while line:
				header = line.strip()
				if " " in header:
					header = header.split(' ')[0]
				seq = f.readline().strip()
				sequences.update( { header: seq } )
				f.readline()	#useless row
				f.readline()	#quality row
				line = f.readline()
	return sequences


def fastq2fasta( read_file, output_folder ):
	"""! @brief convert FASTQ to FASTA file """
	
	fasta_file = output_folder + "reads.fasta"
	if read_file.split('.')[-1] in [ "gz", "gzip", "GZ", "GZIP" ]:	#compressed read input file
		with open( fasta_file, "w" ) as out:
			with gzip.open( read_file, "rb" ) as f:
				line = f.readline().decode('ascii')
				while line:
					out.write( '>' + line + f.readline().decode('ascii') )
					f.readline()	#useless line
					f.readline()	#sequence quality
					line = f.readline().decode('ascii')
	else:
		with open( fasta_file, "w" ) as out:
			with open( read_file, "r" ) as f:
				line = f.readline()
				while line:
					out.write( '>' + line + f.readline() )
					f.readline()	#useless line
					f.readline()	#sequence quality
					line = f.readline()
	return fasta_file


def main( arguments ):
	"""! @brief run everything """
	
	read_file = arguments[ arguments.index( '--reads' )+1 ]
	start_seq_file = arguments[ arguments.index( '--seed' )+1 ]
	output_folder = arguments[ arguments.index( '--out' )+1 ]
	
	if '--rounds' in arguments:
		total_rounds = int( arguments[ arguments.index( '--rounds' )+1 ] )
	else:
		total_rounds = 5	#number of extension rounds to perform
	
	if "--block" in arguments:
		block_size = int( arguments[ arguments.index( '--block' )+1 ] )
	else:
		block_size = 2000	#start seq needs to have this length
	
	if '--direction' in arguments:
		direction = arguments[ arguments.index( '--direction' )+1 ] 
	else:
		direction = "down"	#up|down
	
	if "--simcut" in arguments:
		sim_cut = int( arguments[ arguments.index( '--simcut' )+1 ] )
	else:
		sim_cut = 80	#minimal alignment similarity for BLAST hits to be considered
	
	if "--lencut" in arguments:
		len_cut = int( arguments[ arguments.index( '--lencut' )+1 ] )
	else:
		len_cut = 500	#minimal alignment length for BLAST hits to be considered
	
	if "--threads" in arguments:
		threads = int( arguments[ arguments.index( '--threads' )+1 ] )
	else:
		threads = 8

	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )

	blastdb = output_folder + "blastdb"
	
	read_file_type = get_read_file_type( read_file )
	if read_file_type == "fastq":
			read_file = fastq2fasta( read_file, output_folder )
			read_file_type = "fasta"
	
	if read_file.split('.')[-1] in [ "gz", "gzip", "GZ", "GZIP" ]:
		os.popen( "cp " + read_file + " " + output_folder + "reads.fasta.gz" )
		os.popen( "gunzip " + output_folder + "reads.fasta.gz" )
		read_file = output_folder + "reads.fasta"
	
	os.popen( "makeblastdb -in " + read_file + " -out " + blastdb + " -dbtype nucl" )

	all_reads = load_sequences( read_file )
	blast_results_file = output_folder + "results.txt"

	seed_file = output_folder + "0.fasta"
	os.popen( "cp " + start_seq_file + " " + seed_file )

	result_seq_file = output_folder + "resulting_sequence.fasta"
	doc_file = output_folder + "doc.txt"

	start_seq = load_sequences( seed_file ).values()[0]
	current_round = 1
	cum_len = len( start_seq )
	read_black_list = {}
	fin_seq_collection = []	#collect sequence blocks if running in upstream mode ("up")
	with open( result_seq_file, "w" ) as fin_seq_out:
		with open( doc_file, "w" ) as doc_out:
			if direction == "down":
				fin_seq_out.write( ">result\n" + start_seq )
			else:
				fin_seq_out.write( ">result\n" )
				fin_seq_collection.append( start_seq )
			while current_round < total_rounds:
				# --- run BLAST search --- #
				blast_result_file = output_folder + str( current_round ) + "_blast_hits.txt"
				os.popen( "blastn -query " + seed_file + " -db " + blastdb + " -out " + blast_result_file + " -outfmt 6 -evalue 0.001 -num_threads " + str( threads ) )
				
				# --- find best hit in different (!) read and collect infos for sequence extraction --- #
				next_read_info, seq_to_add_info = analyze_BLAST_results( blast_result_file, block_size, all_reads, direction, sim_cut, len_cut, read_black_list )
				cum_len += abs( seq_to_add_info['start'] - seq_to_add_info['end'] )
				read_black_list.update( { seq_to_add_info['ID']: None } )
				doc_out.write( "\t".join( map( str, [ seq_to_add_info['ID'], seq_to_add_info['start'], seq_to_add_info['end'], str( cum_len ), seq_to_add_info['orientation'] ] ) ) + "\n" )
				
				# --- add new sequence to final seq file --- #
				if seq_to_add_info['orientation']:
					seq = all_reads[ seq_to_add_info['ID'] ][ seq_to_add_info['start']:seq_to_add_info['end'] ]
					if direction == "down":
						fin_seq_out.write( seq )
					else:
						fin_seq_collection.append( seq )
				else:
					seq = revcomp( all_reads[ seq_to_add_info['ID'] ][ seq_to_add_info['start']:seq_to_add_info['end'] ] )
					if direction == "down":
						fin_seq_out.write( seq )
					else:
						fin_seq_collection.append( seq )
				
				# --- take X kb from end of  best hit --- #
				seed_file = output_folder + str( current_round ) + ".fasta"
				with open( seed_file, "w" ) as out:
					seq = all_reads[ next_read_info['ID'] ][ next_read_info['start']:next_read_info['end'] ]
					if next_read_info['orientation']:
						out.write( '>' + next_read_info['ID'] + "\n" + seq + "\n" )
					else:
						out.write( '>' + next_read_info['ID'] + "\n" + revcomp( seq ) + "\n" )
				
				current_round += 1
		if direction == "up":
			fin_seq_out.write( "".join( fin_seq_collection[::-1] ) )


if '--reads' in sys.argv and '--seed' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
