#!/usr/bin/python3
# The MIT License (MIT)
# 
# Copyright (c) 2015 - Vitor C. Piro - PiroV@rki.de - vitorpiro@gmail.com
# Robert Koch-Institut, Germany
# All rights reserved.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import numpy as np
import argparse, subprocess
from dudes.Ranks import Ranks
from dudes.parse_names import parse_names
from collections import defaultdict
import multiprocessing as mp
import time, sys
import pandas as pd

def main():
	
	version = 'v0.08'
	
	total_tx = time.time()
	
	global nodes
		
	parser = argparse.ArgumentParser(description='DUDesDB')
	parser.add_argument('-m', metavar='<reference_mode>', dest="reference_mode", default="av", help="'gi' uses the GI as the identifier (For headers like: >gi|158333233|ref|NC_009925.1|) [NCBI is phasing out sequence GI numbers in September 2016]. 'av' uses the accession.version as the identifier (for headers like: >NC_013791.2).	Default: 'av'")
	parser.add_argument('-f', required=True, metavar='<fasta_file>', dest="fasta_file", help="Reference fasta file for header extraction [only]. Should be the same file used to generate the index database for the SAM file. Each sequence header should contain a identifier as defined in the reference mode.")
	parser.add_argument('-n', required=True, metavar='<nodes_file>', dest="nodes_file", help="nodes.dmp file [from NCBI taxonomy database ftp://ftp.ncbi.nih.gov/pub/taxonomy/]")
	parser.add_argument('-g', required=True, metavar='<ref2tax_files>', dest="ref2tax_files", nargs="*", help="reference id to taxid file(s): 'gi_taxid_nucl.dmp[.gz]' --> 'gi' mode, '*.accession2taxid[.gz]' --> 'av' mode [from NCBI taxonomy database ftp://ftp.ncbi.nih.gov/pub/taxonomy/]")
	parser.add_argument('-a', metavar='<names_file>', dest="names_file", help="names.dmp file [from NCBI taxonomy database ftp://ftp.ncbi.nih.gov/pub/taxonomy/]")
	parser.add_argument('-o', metavar='<output_prefix>', dest="output_prefix", default="dudesdb", help="Output prefix. Default: dudesdb")
	parser.add_argument('-t', metavar='<threads>', dest="threads", type=int, default=1, help="# of threads. Default: 1")
	parser.add_argument('-v', action='version', version='%(prog)s ' + version)
	args = parser.parse_args()
		
	if args.reference_mode not in ['gi','av']:
		print("Reference mode must be 'gi' or 'av'")
		return
	
	# Load fixed ranks
	fixed_ranks = Ranks.ranks
	
	# Load reference identifiers from fasta file (in a set for fast lookup)
	# refids = [identifier]
	sys.stdout.write("Extracting reference identifiers (%s) ..." % args.fasta_file)
	tx = time.time()
	try: # LINUX - faster
		if args.reference_mode=="gi":
			cmd = 'grep -o "^>gi|[0-9]*" ' + args.fasta_file + ' | sed "s/>gi|//g"'
		else:
			cmd = 'grep -o "^>[A-Z0-9_\.]*" ' + args.fasta_file + ' | sed "s/>//g"'
		out = subprocess.getoutput(cmd)
		refids = set(l for l in out.split('\n') if l)
	except: # general
		import re
		if args.reference_mode=="gi":
			regex = re.compile('gi\|[0-9]*')
			slice = 3
		else:
			regex = re.compile('>[A-Z0-9_\.]*')
			slice = 1	
		refids = set()
		for f in open(args.fasta_file,'r'):
			if f[0]==">":
				r = regex.search(f)
				if r: 
					if len(r.group())>1: refids.add(r.group()[slice:])

	if not refids:
		print("\n\tNo references found. Check the reference mode (-m) and your reference headers.")
		return
	else:
		print("\n\t%d unique references found" % len(refids))
	sys.stdout.write("Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")
		
	# Load nodes.dmp
	# nodes = [taxid,parent_taxid,rankid]
	sys.stdout.write("Loading nodes (%s) ..." % args.nodes_file)
	tx = time.time()
	ranks = defaultdict(lambda: len(ranks))
	nodes = [] 
	with open(args.nodes_file,'r') as fnodes:
		for line in fnodes:
			fields = line.split('\t|\t',3)
			nodes.append([int(fields[0]),int(fields[1]),ranks[fields[2]]])
	nodes = np.array(nodes)
	sys.stdout.write(" Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")

	# Load refid 2 taxid file
	# Verify if the entry is being used in the refids before and just output relevant rows
	# refid_taxid = [refid,taxid]
	sys.stdout.write("Loading taxids (%s) ..." % ",".join(args.ref2tax_files))
	tx = time.time()
	refids_lookup = defaultdict(lambda: len(refids_lookup)) #Auto-increase id dict
	refid_taxid = []
	if args.reference_mode=="gi":
		for file in args.ref2tax_files:
			refid_taxid.append(np.array(pd.read_csv(file, compression='gzip' if file.endswith(".gz") else None, sep='\t', header=None, dtype=int, converters={0:lambda x: refids_lookup[x] if x in refids else np.nan}).dropna(how='any'),dtype=int))
	else:
		for file in args.ref2tax_files:
			refid_taxid.append(np.array(pd.read_csv(file, compression='gzip' if file.endswith(".gz") else None, sep='\t', header=None, skiprows=1, usecols=[1,2], converters={1:lambda x: refids_lookup[x] if x in refids else np.nan}).dropna(how='any'),dtype=int))
	# Concatenate files together
	refid_taxid = np.concatenate(refid_taxid, axis=0)
	sys.stdout.write(" Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")
	
	# Output differences and verify the taxids on nodes file
	sys.stdout.write("Parsing nodes and taxids ...")
	tx = time.time()
	unique_refs_used = np.unique(refid_taxid[:,0]).size # found on the ref2tax files - should be already unique, just to make sure
	if not unique_refs_used:
		print("\n\tNo matches on taxonomy found. Check your reference to taxonomy file[s] (-g): ", ",".join(args.ref2tax_files))
		return
	
	seq_without_refid = len(refids) - unique_refs_used
	sys.stdout.write("\n\tIgnoring %d (out of %d) sequences without entries on %s\n" % ((seq_without_refid, len(refids), ",".join(args.ref2tax_files))))
	if seq_without_refid: 
		refids_nparr = np.array(list(refids))
		print("\n".join(refids_nparr[~np.in1d(refids_nparr, list(refids_lookup.keys()))]))
		del refids_nparr
		
	refid_with_valid_taxid = np.in1d(refid_taxid[:,1],nodes[:,0])
	seq_without_taxid = unique_refs_used - sum(refid_with_valid_taxid)
	sys.stdout.write("\tIgnoring %d (out of %d) sequences without taxid on %s\n" % (seq_without_taxid, unique_refs_used, args.nodes_file))
	if seq_without_taxid:
		refids_lookup_rev = {v: k for k, v in refids_lookup.items()}
		print("\n".join([refids_lookup_rev[r] + "\t" + str(t) for r,t in refid_taxid[~refid_with_valid_taxid]]))

		# filter from lookup (rev) and overwrite current lookup
		for r in refid_taxid[~refid_with_valid_taxid,0]:
			del refids_lookup_rev[r]
		refids_lookup = {v: k for k, v in refids_lookup_rev.items()}
		del refids_lookup_rev
		# filter out entries without taxid matches on nodes.dmp
		refid_taxid = refid_taxid[refid_with_valid_taxid]
		
	sys.stdout.write("Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")
	
	unique_refs_used_filtered = np.unique(refid_taxid[:,0]).size
	if not unique_refs_used_filtered:
		print("\n\tNo matches on nodes found. Check your file (-n): ", args.nodes_file)
		return
		
	# ------- MakeDB -----------
	sys.stdout.write("Creating database ...")
	tx = time.time()
	# Map rankids from nodes (all ranks) with the id from the fixed ranks
	fixed_ranks_id = dict()
	for i,r in enumerate(fixed_ranks):
		if r=="strain": #Special case - strain on my list - no rank on the nodes
			fixed_ranks_id[ranks["no rank"]] = i
		else:
			fixed_ranks_id[ranks[r]] = i

	# Generate all possible paths
	pool = mp.Pool(args.threads)
	paths = defaultdict(list)
	res = {}
	for taxid in np.unique(refid_taxid[:,1]):
		res[taxid] = pool.apply_async(generatePath,args=(taxid,fixed_ranks_id,))
	for taxid,r in res.items():
		paths[taxid] = r.get()
		
	# Generate data structure
	no_rank_after = fixed_ranks.index('species')
	no_rank_id = fixed_ranks.index('strain')
	refid_nodes = []
	res = []
	for refid,taxid in refid_taxid: 
		res.append(pool.apply_async(generateDB,args=(paths[taxid],refid,no_rank_after,no_rank_id,)))
	for r in res:
		refid_nodes.extend(r.get())
	pool.close()
	pool.join()	
	sys.stdout.write(" Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")

	# ------- Names -----------
	sys.stdout.write("Parsing names ...")
	tx = time.time()
	if args.names_file:
		unique_taxids = set(np.unique([entry[0] for path in list(paths.values()) for entry in path]))
		names = parse_names(args.names_file,unique_taxids)
	else:
		names = {}
	sys.stdout.write(" Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")
	
	# Save database
	sys.stdout.write("Saving database %s ..." % (args.output_prefix + ".npz"))
	tx = time.time()
	np.savez(args.output_prefix, version=version, reference_mode=args.reference_mode, refid_nodes=refid_nodes, refids_lookup=dict(refids_lookup), names=names)
	sys.stdout.write(" Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")
	
	sys.stdout.write("\nTotal elapsed time: " + str(time.time() - total_tx) + " seconds\n")
	
def generateDB(ptx,refid,no_rank_after,no_rank_id):
	refid_no = []
	parent_taxid = 1
	no_rank = False
	no_rank_count = 0
	for taxid,rankid in ptx[::-1]:
		if (no_rank==False and rankid!=no_rank_id) or no_rank: 
			refid_no.append([refid,taxid,parent_taxid,rankid+no_rank_count])
			parent_taxid = taxid
			if rankid==no_rank_after: no_rank = True # Only allow no_rank (strain) after some certain rank
			if rankid==no_rank_id: no_rank_count += 1
	return refid_no
	
def generatePath(taxid,fixed_ranks_id):
	path = []
	while taxid!=1:
		taxid, parent_taxid, rankid = nodes[nodes[:,0]==taxid,[0,1,2]]
		if rankid in list(fixed_ranks_id.keys()):
			path.append((taxid,fixed_ranks_id[rankid])) #use the rankid from fixed ranks
		taxid = parent_taxid
	return path
	
if __name__ == "__main__":
	main()
