import re, os
import numpy as np
import multiprocessing as mp
from itertools import islice
from collections import defaultdict

nm_regex = re.compile('NM:i:[0-9]*')

def parse_lines(lines):
	sam = np.zeros((len(lines),3),dtype=int)
	reads = []
	c = 0
	for l in lines:
		fields = l.split('\t')
		if fields[2]!="*":
			cig = {'I':0,'D':0,'M':0}
			for val,ci in re.findall('(\d+)([IDM])',fields[5]): cig[ci] += int(val)
			sam[c,0] = refids_lookup.get(refid_regex.search(fields[2]).group()[slice:],-1) # refid code from lookup (gi or av) or -1
			sam[c,1] = int(fields[3]) # POS
			sam[c,2] = (cig['M'] + cig['I']) - (int(nm_regex.search(l).group()[5:]) + cig['D'] + cig['I']) #MATCHES -> LEN - (NM + INDELS)
			reads.append(fields[0])
			c = c + 1
	return sam[:c,:],reads

def parse_lines_extended(lines):
	sam = np.zeros((len(lines),3),dtype=int)
	reads = []
	c = 0
	for l in lines:
		fields = l.split('\t')
		if fields[2]!="*":
			cig = {'I':0,'D':0,'M':0,'=':0,'X':0}
			for val,ci in re.findall('(\d+)([IDM=X])',fields[5]): cig[ci] += int(val)
			sam[c,0] = refids_lookup.get(refid_regex.search(fields[2]).group()[slice:],-1)  # refid code from lookup (gi or av) or -1
			sam[c,1] = int(fields[3]) # POS
			sam[c,2] = (cig['X'] + cig['='] + cig['M'] + cig['I']) - (cig['X'] + cig['D'] + cig['I']) #MATCHES -> LEN - (NM + INDELS)
			reads.append(fields[0])
			c = c + 1
	return sam[:c,:],reads

def parse_sam(sam_file,sam_format,rl,reference_mode,threads):
	
	global refids_lookup
	global refid_regex
	global slice
	
	# Global refids lookup
	refids_lookup = rl
	
	if reference_mode=="gi":
		refid_regex = re.compile('gi\|[0-9]*')
		slice = 3
	else:
		refid_regex = re.compile('[A-Z0-9_\.]*') #without ">" from fasta regex
		slice = 0
	
	refs = []
	reads = defaultdict(lambda: len(reads))
	pool = mp.Pool(threads)

	file = open(sam_file, 'r')
	while True:
		last_pos = file.tell()
		l = file.readline()
		fields = l.split('\t')
		if fields[0]=='@SQ': # list of references on the database
			# example fields:
			# 0: @SQ
			# 1: SN:NC_010571.1
			# 2: LN:5957605
			# refs: list of lists, each sublist is a pair of
			# 1. refids_lookup.get(accession)->dudes internal id for the reference accession
			# 2. LN attribute
			refs.append([refids_lookup.get(refid_regex.search(fields[1][3:]).group()[slice:], -1), int(fields[2][3:])])
		elif l[0]=='@': continue # other headers
		else: 
			first_aln_line = l
			break
	file.seek(last_pos) #return to the last iterated line (first sam result)
	
	n_lines = 1000
	res = []	
	
	# Check for NM flag
	if sam_format=="nm" and not nm_regex.search(first_aln_line):
		print("\n -- Warning: NM flag not found. Switching to extended CIGAR mode (-i 'ex')")
		sam_format = "ex"
	
	if sam_format=="nm":
		while True:		
			next_n_lines = list(islice(file, n_lines)) # select chunk of lines
			if not next_n_lines: break
			res.append(pool.apply_async(parse_lines,args=(next_n_lines,))) # start the process in the pool
	elif sam_format=="ex":
		while True:		
			next_n_lines = list(islice(file, n_lines)) # select chunk of lines
			if not next_n_lines: break
			res.append(pool.apply_async(parse_lines_extended,args=(next_n_lines,))) # start the process in the pool

	sam_all = []
	for r in res:
		# assign unique id for each read (reads defaultdict)
		if r.get()[0].any():
			rd = np.array([reads[read] for read in r.get()[1]]) # READID
			sm = r.get()[0] # sam data
			sam_all.append(np.insert(sm,sm[0,:].size,rd,axis=1)) # append readid column in the data array and add to the results

	pool.close()
	pool.join()	
	file.close()
	
	return np.concatenate(sam_all),np.array(refs)