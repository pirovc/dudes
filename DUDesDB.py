#!/usr/bin/env python3
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
from dudes.Ranks import ranks as FIXED_RANKS
from dudes.parse_names import parse_names
from collections import defaultdict
import multiprocessing as mp
import time, sys
import pandas as pd
from dudes import VERSION


def main():
    total_tx = time.time()

    global nodes

    args = parse_args()

    refids = get_reference_identifiers(args.fasta_files, args.reference_mode)

    ranks = load_nodes(args.nodes_file)

    # Load refid 2 taxid file
    # Verify if the entry is being used in the refids before and just output relevant rows
    # refid_taxid = [refid,taxid]
    sys.stdout.write("Loading taxids (%s) ..." % ",".join(args.ref2tax_files))
    tx = time.time()
    refid_taxid, refids_lookup = load_refid2taxid_files(args.ref2tax_files, args.reference_mode, refids)
    sys.stdout.write(" Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")

    # Output differences and verify the taxids on nodes file
    sys.stdout.write("Parsing nodes and taxids ...")
    tx = time.time()
    refid_taxid, refids_lookup = remove_refs_without_taxids(refid_taxid, refids_lookup, nodes, args.ref2tax_files,
                                                            args.nodes_file, refids)
    sys.stdout.write("Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")

    unique_refs_used_filtered = np.unique(refid_taxid[:, 0]).size
    if not unique_refs_used_filtered:
        print("\n\tNo matches on nodes found. Check your file (-n): ", args.nodes_file)
        sys.exit(1)

    # ------- MakeDB -----------
    sys.stdout.write("Creating database ...")
    tx = time.time()
    paths, refid_nodes = generate_paths_and_database(args, FIXED_RANKS, ranks, refid_taxid)
    sys.stdout.write(" Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")

    # ------- Names -----------
    sys.stdout.write("Parsing names ...")
    tx = time.time()
    if args.names_file:
        unique_taxids = set(np.unique([entry[0] for path in list(paths.values()) for entry in path]))
        names = parse_names(args.names_file, unique_taxids)
    else:
        names = {}
    sys.stdout.write(" Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")

    save_database(args.output_prefix, args.reference_mode, names, refid_nodes, refids_lookup)

    sys.stdout.write("\nTotal elapsed time: " + str(time.time() - total_tx) + " seconds\n")


def generate_paths_and_database(args, fixed_ranks, ranks, refid_taxid):
    fixed_ranks_id = get_nodes_rank_to_fixed_rank_id_map(fixed_ranks, ranks)
    # Generate all possible taxonomy paths
    pool = mp.Pool(args.threads)
    paths = generate_all_possible_taxonomy_paths(fixed_ranks_id, pool, refid_taxid)
    refid_nodes = generate_database(fixed_ranks, paths, pool, refid_taxid)
    pool.close()
    pool.join()
    return paths, refid_nodes


def remove_refs_without_taxids(refid_taxid, refids_lookup, nodes, ref2tax_filenames, nodes_filename, refids):
    unique_refs_used = np.unique(
        refid_taxid[:, 0]).size  # found on the ref2tax files - should be already unique, just to make sure
    if not unique_refs_used:
        print("\n\tNo matches on taxonomy found. Check your reference to taxonomy file[s] (-g): ",
              ",".join(ref2tax_filenames))
        sys.exit(1)
    seq_without_refid = len(refids) - unique_refs_used
    sys.stdout.write("\n\tIgnoring %d (out of %d) sequences without entries on %s\n" % (
        (seq_without_refid, len(refids), ",".join(ref2tax_filenames))))
    if seq_without_refid:
        refids_nparr = np.array(list(refids))
        print("\n".join(refids_nparr[~np.in1d(refids_nparr, list(refids_lookup.keys()))]))
        del refids_nparr
    # 												  taxids from nodes file
    refid_with_valid_taxid = np.in1d(refid_taxid[:, 1], nodes[:, 0])
    seq_without_taxid = unique_refs_used - sum(refid_with_valid_taxid)
    sys.stdout.write("\tIgnoring %d (out of %d) sequences without taxid on %s\n" % (
        seq_without_taxid, unique_refs_used, nodes_filename))
    if seq_without_taxid:
        refids_lookup_rev = {v: k for k, v in refids_lookup.items()}
        print("\n".join([refids_lookup_rev[r] + "\t" + str(t) for r, t in refid_taxid[~refid_with_valid_taxid]]))

        # filter from lookup (rev) and overwrite current lookup
        for r in refid_taxid[~refid_with_valid_taxid, 0]:
            del refids_lookup_rev[r]
        refids_lookup = {v: k for k, v in refids_lookup_rev.items()}
        del refids_lookup_rev
        # filter out entries without taxid matches on nodes.dmp
        refid_taxid = refid_taxid[refid_with_valid_taxid]
    return refid_taxid, refids_lookup


def get_nodes_rank_to_fixed_rank_id_map(fixed_ranks, ranks):
    # Map rankids from nodes (all ranks) with the id from the fixed ranks
    fixed_ranks_id = dict()
    for i, r in enumerate(fixed_ranks):
        if r == "strain":  # Special case - strain on my list - no rank on the nodes
            fixed_ranks_id[ranks["no rank"]] = i
        else:
            fixed_ranks_id[ranks[r]] = i
    return fixed_ranks_id


def save_database(output_prefix, reference_mode, names, refid_nodes, refids_lookup):
    # Save database
    sys.stdout.write("Saving database %s ..." % (output_prefix + ".npz"))
    tx = time.time()
    np.savez(output_prefix, version=VERSION, reference_mode=reference_mode, refid_nodes=refid_nodes,
             refids_lookup=dict(refids_lookup), names=names)
    sys.stdout.write(" Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")


def generate_database(fixed_ranks, paths, pool, refid_taxid):
    # Generate data structure
    no_rank_after = fixed_ranks.index('species')
    no_rank_id = fixed_ranks.index('strain')
    refid_nodes = []
    res = []
    for refid, taxid in refid_taxid:
        res.append(
            pool.apply_async(generate_database_for_single_path, args=(paths[taxid], refid, no_rank_after, no_rank_id,)))
    for r in res:
        refid_nodes.extend(r.get())
    return refid_nodes


def generate_database_for_single_path(ptx, refid, no_rank_after, no_rank_id):
    """
    generate database for a single path

    :param ptx: lineage of refid, list of tuples, each tuple is a pair of taxid and it's index in fixed ranks, list is
    sorted from child to parent
    :param refid: id of reference accession
    :param no_rank_after: index id of most specific rank name in fixed ranks, after which the next rank is assigned
    "no rank" in taxonomy nodes dmp
    :param no_rank_id: index id of rank name "no rank" in fixed ranks
    :return: list of lists, one list per lineage taxid, ordered from parent to child,
    each list contains: query reference id, taxid from ptx, parent taxid of taxid from ptx, sum of rankid and count of
    no rank ids in the parent taxons between current taxid and taxid of no_rank_after id.
    occurrences of no_rank_id in the path before no_rank_after are skipped
    """
    refid_no = []
    parent_taxid = 1
    allow_no_rank = False
    no_rank_count = 0
    for taxid, rankid in ptx[::-1]:
        if (allow_no_rank == False and rankid != no_rank_id) or allow_no_rank:
            refid_no.append([refid, taxid, parent_taxid, rankid + no_rank_count])
            parent_taxid = taxid
            if rankid == no_rank_after: allow_no_rank = True  # Only allow no_rank (strain) after some certain rank
            if rankid == no_rank_id: no_rank_count += 1  # only reachable after no_rank_after was reached
    return refid_no


def generate_all_possible_taxonomy_paths(fixed_ranks_id, pool, refid_taxid):
    paths = defaultdict(list)
    res = {}
    for taxid in np.unique(refid_taxid[:, 1]):
        res[taxid] = pool.apply_async(generatePath, args=(taxid, fixed_ranks_id,))
    for taxid, r in res.items():
        paths[taxid] = r.get()
    return paths


def load_refid2taxid_files(refid2taxid_files, ref_mode, refids):
    refids_lookup = defaultdict(lambda: len(refids_lookup))  # Auto-increase id dict
    refid_taxid = []
    if ref_mode == "gi":
        for file in refid2taxid_files:
            refid_taxid.append(np.array(pd.read_csv(
                file, compression='gzip' if file.endswith(".gz") else None, sep='\t',
                header=None, dtype=int,
                converters={0: lambda x: refids_lookup[x] if x in refids else np.nan}
            ).dropna(how='any'), dtype=int))
    elif ref_mode == "av":
        for file in refid2taxid_files:
            refid_taxid.append(np.array(
                pd.read_csv(file, compression='gzip' if file.endswith(".gz") else None, sep='\t', header=None,
                            skiprows=1, usecols=[1, 2],
                            converters={1: lambda x: refids_lookup[x] if x in refids else np.nan}).dropna(how='any'),
                dtype=int))
    else:    # ref_mode == "up"
        for file in refid2taxid_files:
            refid_taxid.append(np.array(
                pd.read_csv(file, compression='gzip' if file.endswith(".gz") else None, sep='\t', header=None,
                            usecols=[0, 12],
                            converters={0: lambda x: refids_lookup[x] if x in refids else np.nan}).dropna(how='any'),
                dtype=int))
    # Concatenate files together
    refid_taxid = np.concatenate(refid_taxid, axis=0)
    return refid_taxid, refids_lookup


def load_nodes(nodes_file):
    global nodes
    # Load nodes.dmp
    # nodes = [taxid,parent_taxid,rankid]
    sys.stdout.write("Loading nodes (%s) ..." % nodes_file)
    tx = time.time()
    # ranks: key: rank name, value: unique ID
    ranks = defaultdict(lambda: len(ranks))
    nodes = []
    with open(nodes_file, 'r') as fnodes:
        for line in fnodes:
            fields = line.split('\t|\t', 3)
            #            taxid           parent_taxid,  rank name
            nodes.append([int(fields[0]), int(fields[1]), ranks[fields[2]]])
    nodes = np.array(nodes)
    sys.stdout.write(" Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")
    return ranks


def get_reference_identifiers(fasta_files, ref_mode):
    # Load reference identifiers from fasta file (in a set for fast lookup)
    # refids = {identifier}
    sys.stdout.write("Extracting reference identifiers (%s) ..." % ",".join(fasta_files))
    tx = time.time()
    # try:  # shell - faster
    if ref_mode == "gi":
        cmd = r'zgrep -h -o "^>gi|[0-9]*" ' + " ".join(fasta_files) + ' | sed "s/>gi|//g"'
    elif ref_mode == "av":
        cmd = r'zgrep -h -o "^>[A-Z0-9_\.]*" ' + " ".join(fasta_files) + ' | sed "s/>//g"'
    elif ref_mode == "up":
        cmd = r'zgrep -hoE "^>..\|[^|\s]+" ' + " ".join(fasta_files) + ' | sed "s/>..|//g"'
    return_code, out = subprocess.getstatusoutput("set -o pipefail; " + cmd)
    if return_code == 0:
        refids = set(l for l in out.split('\n') if l)
    else:    # python
        import re, gzip
        if ref_mode == "gi":
            regex = re.compile(r'gi\|[0-9]*')
            slice = 3
        elif ref_mode == "av":
            regex = re.compile(r'>[A-Z0-9_\.]*')
            slice = 1
        elif ref_mode == "up":
            regex = re.compile(r'(?:sp|tr)\|[A-Z0-9_.]+')
            slice = 3
        refids = set()
        for file in fasta_files:
            f = gzip.open(file, 'rt') if file.endswith(".gz") else open(file, 'r')
            for line in f:
                if line[0] == ">":
                    r = regex.search(line)
                    if r:
                        if len(r.group()) > 1: refids.add(r.group()[slice:])
    if not refids:
        print("\n\tNo references found. Check the reference mode (-m) and your reference headers.")
        sys.exit(1)
    else:
        print("\n\t%d unique references found" % len(refids))
    sys.stdout.write("Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")
    return refids


def parse_args():
    parser = argparse.ArgumentParser(prog='DUDesDB.py')
    parser.add_argument('-m', metavar='<reference_mode>', dest="reference_mode", default="av",
                        choices=['gi', 'av', 'up'],
                        help="""'gi' uses the GI as the identifier (For headers like: 
                        >gi|158333233|ref|NC_009925.1|) [NCBI is phasing out sequence GI numbers in September 2016]. 
                        'av' uses the accession.version as the identifier (for headers like: >NC_013791.2). 'up' uses 
                        the uniprot accession as identifier (for headers like: >sp|Q197F8|... Default: 'av'""")
    parser.add_argument('-f', required=True, metavar='<fasta_files>', dest="fasta_files", nargs="*",
                        help="Reference fasta file(s) for header extraction only, plain or gzipped -  the same file "
                             "used to generate the read mapping index. Each sequence header '>' should contain a "
                             "identifier as defined in the reference mode.")
    parser.add_argument('-g', required=True, metavar='<ref2tax_files>', dest="ref2tax_files", nargs="*",
                        help="reference id to taxid file(s): "
                             "'gi_taxid_nucl.dmp[.gz]' --> 'gi' mode, "
                             "'*.accession2taxid[.gz]' --> 'av' mode "
                             "[from NCBI taxonomy database ftp://ftp.ncbi.nih.gov/pub/taxonomy/]"
                             "'idmapping_selected.tab[.gz]' --> 'up' mode "
                             "[from https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz")
    parser.add_argument('-n', required=True, metavar='<nodes_file>', dest="nodes_file",
                        help="nodes.dmp file [from NCBI taxonomy database ftp://ftp.ncbi.nih.gov/pub/taxonomy/]")
    parser.add_argument('-a', metavar='<names_file>', dest="names_file",
                        help="names.dmp file [from NCBI taxonomy database ftp://ftp.ncbi.nih.gov/pub/taxonomy/]")
    parser.add_argument('-o', metavar='<output_prefix>', dest="output_prefix", default="dudesdb",
                        help="Output prefix. Default: dudesdb")
    parser.add_argument('-t', metavar='<threads>', dest="threads", type=int, default=1, help="# of threads. Default: 1")
    parser.add_argument('-v', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def generatePath(taxid, fixed_ranks_id):
    """
	generate the lineage path of taxids, only return taxids of fixed ranks.

	:param taxid: query taxid
	:param fixed_ranks_id: dict, key: id of first rank name occurrence rank in nodes_file, value: rank position in
	fixed ranks
	:return: list of tuples, each tuple is a pair of taxid and it's fixed rank index, list is sorted from child to
	parent
	"""
    path = []
    while taxid != 1:
        taxid, parent_taxid, rankid = nodes[nodes[:, 0] == taxid, [0, 1, 2]]
        if rankid in list(fixed_ranks_id.keys()):
            path.append((taxid, fixed_ranks_id[rankid]))  # use the rankid from fixed ranks
        taxid = parent_taxid
    return path


if __name__ == "__main__":
    main()
