import re, os
import numpy as np
import multiprocessing as mp
from itertools import islice
from collections import defaultdict

nm_regex = re.compile(r'(?<=NM:i:)[^\t]+')
xa_regex = re.compile(r'(?<=XA:Z:)[^\t]+')

def parse_lines(lines):
	sam = []
	reads = []
	for l in lines:
		fields = l.split('\t')
		if not int(fields[1]) & 4 : #bitwise flag check (sam flag for read unmapped (0x4))
			if fields[2]=="*": continue # no reference sequence

			if fields[9]=="*": # no sequences on sam file
				seq_len = sum(map(int,re.findall("[0-9]+(?=M)",fields[5]))) # sequence len. from cigar string M's
			else: 
				seq_len = len(fields[9]) # sequence len. directly from sequence
			
			# Edit distance from NM option field (0 if not found)
			try:
			    edit_dis = int(nm_regex.search('\t'.join(fields[11:])).group(0))
			except:
				edit_dis = 0

			sam.append([refids_lookup.get(refid_regex.search(fields[2]).group(),-1), # refid code from lookup (gi or av) or -1
						int(fields[3]), # POS
						seq_len-edit_dis]) # read score
			reads.append(fields[0])

			# look for secondary matches on XA option field (e.g. XA:Z:NC_010364.1,1119958,1120108,-,3;NZ_CP011246.1,1191943,1192093,+,2;)
			if sec_aln:
				sec_aln_sam = parse_sec_aln('\t'.join(fields[11:]))
				sam.extend(sec_aln_sam) # add to current sam
				reads.extend([fields[0]]*len(sec_aln_sam)) # add the read identifier * the number of secondary alignments

	return np.array(sam), reads

def parse_lines_extended(lines):
	sam = []
	reads = []
	for l in lines:
		fields = l.split('\t')
		if not int(fields[1]) & 4 : #bitwise flag check (sam flag for read unmapped (0x4))
			if fields[2]=="*": continue # no reference sequence
			cig = {'I':0,'D':0,'M':0,'=':0,'X':0,'S':0}
			for val,ci in re.findall('(\d+)([IDM=XS])',fields[5]): cig[ci] += int(val)
			sam.append([refids_lookup.get(refid_regex.search(fields[2]).group(),-1),  # refid code from lookup (gi or av) or -1
						int(fields[3]), # POS
						(cig['M'] + cig['I'] + cig['S'] + cig['='] + cig['X']) - (cig['X'] + cig['D'])])
			reads.append(fields[0])
			# look for secondary matches on XA option field (e.g. XA:Z:NC_010364.1,1119958,1120108,-,3;NZ_CP011246.1,1191943,1192093,+,2;)
			if sec_aln:
				sec_aln_sam = parse_sec_aln('\t'.join(fields[11:]))
				sam.extend(sec_aln_sam) # add to current sam
				reads.extend([fields[0]]*len(sec_aln_sam)) # add the read identifier * the number of secondary alignments

	return np.array(sam), reads

def parse_sec_aln(opt_fields):
	sec_aln_sam = []
	xa = xa_regex.search(opt_fields)
	if xa:
		for sec_match in xa.group(0).split(";")[:-1]: # [:-1] -> always an empty or \n
			ref,begin,end,_,edit_dis = sec_match.split(",")
			sec_aln_sam.append([refids_lookup.get(refid_regex.search(ref).group(),-1),
						int(begin),
						(int(end) - int(begin)) - int(edit_dis)])
	return sec_aln_sam

def parse_sam(sam_file, sam_format, sa, rl, reference_mode, threads):
	
	global sec_aln
	global refids_lookup
	global refid_regex
	
	# Global refids lookup
	refids_lookup = rl

	# Global secondary alignment flag
	sec_aln = sa
	
	# Regular expression based on reference_mode
	refid_regex = re.compile('(?<=gi\|)[0-9]+') if reference_mode=="gi" else re.compile('[A-Z0-9_\.]*')

	refs = []
	reads = defaultdict(lambda: len(reads))
	pool = mp.Pool(threads)

	file = open(sam_file, 'r')
	while True:
		last_pos = file.tell()
		l = file.readline()
		fields = l.split('\t')
		if fields[0]=='@SQ': # list of references on the database
			refs.append([refids_lookup.get(refid_regex.search(fields[1][3:]).group(),-1),int(fields[2][3:])])
		elif l[0]=='@': continue # other headers
		else: 
			first_aln_line = l
			break
	file.seek(last_pos) #return to the last iterated line (first sam result)
	
	n_lines = 50000
	res = []
	
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