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
from Ranks import Ranks
from parse_names import parse_names
from collections import defaultdict
import multiprocessing as mp
import time, sys
import shelve

def main():
	
	version = 'v0.07'
	
	total_tx = time.time()
	
	global nodes
		
	parser = argparse.ArgumentParser(description='DUDesDB')
	parser.add_argument('-m', metavar='<reference_mode>', dest="reference_mode", default="gi", help="'gi' uses the GI as the identifier (For headers like: >gi|158333233|ref|NC_009925.1|) [NCBI is phasing out sequence GI numbers in September 2016]. 'av' uses the accession.version as the identifier (for headers like: >NC_013791.2).	Default: 'gi'")
	parser.add_argument('-f', required=True, metavar='<fasta_file>', dest="fasta_file", help="Reference fasta file for header extraction [only]. Should be the same file used to generate the index database for the SAM file. Each sequence header should contain a identifier as defined in the reference mode.")
	parser.add_argument('-n', required=True, metavar='<nodes_file>', dest="nodes_file", help="nodes.dmp file [from NCBI taxonomy database ftp://ftp.ncbi.nih.gov/pub/taxonomy/]")
	parser.add_argument('-g', required=True, metavar='<ref2tax_file>', dest="ref2tax_file", help="Reference identifier (GI or accesion.version) to taxonomic ID file: 'gi_taxid_nucl.dmp' --> 'gi' mode, 'accession2taxid/nucl_gb.dump' --> 'av' mode [from NCBI taxonomy database ftp://ftp.ncbi.nih.gov/pub/taxonomy/]")
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
			cmd = 'grep -o "gi|[0-9]*" ' + args.fasta_file + ' | sed "s/gi|//g"'
		else:
			cmd = 'grep -o ">[A-Z0-9_\.]*" ' + args.fasta_file + ' | sed "s/>//g"'
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
	sys.stdout.write(" Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")
	
	if not refids:
		print("No references found. Check the reference mode (-m) and the patter from your fasta file")
		return
		
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
	sys.stdout.write("Loading taxids (%s) ..." % args.ref2tax_file)
	tx = time.time()
	refids_lookup = defaultdict(lambda: len(refids_lookup)) #Auto-increase id dict
	try: # pandas - faster
		import pandas as pd
		if args.reference_mode=="gi":
			refid_taxid = np.array(pd.read_csv(args.ref2tax_file, sep='\t', header=None, dtype=int, converters={0:lambda x: refids_lookup[x] if x in refids else np.nan}).dropna(how='any'),dtype=int)
		else:
			refid_taxid = np.array(pd.read_csv(args.ref2tax_file, sep='\t', header=None, skiprows=1, usecols=[1,2], converters={1:lambda x: refids_lookup[x] if x in refids else np.nan}).dropna(how='any'),dtype=int)
	except: # general
		refid_taxid = []
		if args.reference_mode=="gi":
			for f in open(args.ref2tax_file):
				fields = f.split('\t')
				if fields[0] in refids:
					refid_taxid.append([refids_lookup[fields[0]],int(fields[0])])
		else:
			with open(args.ref2tax_file,"r") as file:
				next(file)
				for f in file:
					fields = f.split('\t')
					if fields[1] in refids:
						refid_taxid.append([refids_lookup[fields[1]],int(fields[2])])
		refid_taxid = np.array(refid_taxid,dtype=int)
	sys.stdout.write(" Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")
	
	# Output differences and verify the taxids on nodes file
	sys.stdout.write("Parsing nodes and taxids ...")
	tx = time.time()
	unique_refs_used = np.unique(refid_taxid[:,0]).size
	sys.stdout.write("\n\t%d references without entries on %s" % ((len(refids) - unique_refs_used,args.ref2tax_file)))
	#print(np.setdiff1d(refids,refids_lookup.keys()))
	sub_refid_taxid = refid_taxid[np.in1d(refid_taxid[:,1],nodes[:,0])]
	sys.stdout.write("\n\t%d references without entries on %s\n" % (unique_refs_used - sub_refid_taxid[:,0].size,args.nodes_file))
	#print(refid_taxid[~np.in1d(refid_taxid[:,1],nodes[:,0])])
	sys.stdout.write(" Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")
	
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
	for taxid in np.unique(sub_refid_taxid[:,1]):
		res[taxid] = pool.apply_async(generatePath,args=(taxid,fixed_ranks_id,))
	for taxid,r in res.items():
		paths[taxid] = r.get()
		
	# Generate data structure
	no_rank_after = fixed_ranks.index('species')
	no_rank_id = fixed_ranks.index('strain')
	refid_nodes = []
	res = []
	for refid,taxid in sub_refid_taxid: 
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
