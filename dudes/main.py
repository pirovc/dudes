import argparse
import multiprocessing as mp
import os
import sys
import time
from collections import defaultdict

import numpy as np

from dudes import VERSION
from dudes.Bins import Bins
from dudes.Ident import Ident
from dudes.Names import Names
from dudes.Ranks import Ranks
from dudes.Refs import Refs
from dudes.Rep import Rep
from dudes.SMap import SMap
from dudes.TTree import TTree
from dudes.Util import group_max, printDebug, getIndexRank, getNameRank
from dudes.parse_sam import parse_sam

np.set_printoptions(
    threshold=10000,
    suppress=True,
    formatter={"float": "{: 0.6f}".format},
    linewidth=1000,
)


def main():
    version = VERSION

    total_tx = time.time()

    global DEBUG
    global DEBUG_PLOTS_DIR
    global threads
    global fixed_ranks
    global thr_alpha
    global bin_cutoff
    global permutations

    global total_matches_filter
    global min_reference_matches_percent
    global min_reference_matches_number
    global min_group_size
    global bin_size

    parser = argparse.ArgumentParser(prog="DUDes.py")
    parser.add_argument(
        "-s",
        required=True,
        metavar="<sam_file>",
        dest="sam_file",
        help="Alignment/mapping file in SAM format. DUDes does not depend on any specific read mapper, but it requires header information (@SQ SN:gi|556555098|ref|NC_022650.1| LN:55956) and mismatch information (check -i)",
    )
    parser.add_argument(
        "-d",
        required=True,
        metavar="<database_file>",
        dest="database_file",
        help="Database file (output from DUDesDB [.npz])",
    )
    parser.add_argument(
        "-i",
        metavar="<sam_format>",
        dest="sam_format",
        default="nm",
        help="SAM file format ['nm': sam file with standard cigar string plus NM flag (NM:i:[0-9]*) for mismatches count | 'ex': just the extended cigar string]. Default: 'nm'",
    )
    parser.add_argument(
        "-t",
        metavar="<threads>",
        dest="threads",
        type=int,
        default=1,
        help="# of threads. Default: 1",
    )
    parser.add_argument(
        "-x",
        metavar="<taxid_start>",
        dest="taxid_start",
        type=int,
        default=1,
        help="Taxonomic Id used to start the analysis (1 = root). Default: 1",
    )
    parser.add_argument(
        "-m",
        metavar="<max_read_matches>",
        dest="max_read_matches",
        type=float,
        default=0,
        help="Keep reads up to this number/percentile of matches (0: off / 0-1: percentile / >=1: match count). Default: 0",
    )
    parser.add_argument(
        "-a",
        metavar="<min_reference_matches>",
        dest="min_reference_matches",
        type=float,
        default=0.001,
        help="Minimum number/percentage of supporting matches to consider the reference (0: off / 0-1: percentage / >=1: read number). Default: 0.001",
    )
    parser.add_argument(
        "-l",
        metavar="<last_rank>",
        dest="last_rank",
        default="species",
        help="Last considered rank [" + ",".join(Ranks.ranks) + "]. Default: 'species'",
    )
    parser.add_argument(
        "-b",
        metavar="<bin_size>",
        dest="bin_size",
        type=float,
        default=0.25,
        help="Bin size (0-1: percentile from the lengths of all references in the database / >=1: bp). Default: 0.25",
    )
    parser.add_argument(
        "-o",
        metavar="<output_prefix>",
        dest="output_prefix",
        default="",
        help="Output prefix. Default: STDOUT",
    )
    parser.add_argument(
        "--debug", action="store_true", help="print debug info to STDERR"
    )
    parser.add_argument(
        "--debug_plots_dir",
        default="",
        help="path to directory for writing debug plots to.",
    )
    parser.add_argument("-v", action="version", version="%(prog)s " + version)
    args = parser.parse_args()

    # Constants
    threads = args.threads
    fixed_ranks = Ranks.ranks
    taxid_start = args.taxid_start
    thr_alpha = 0.05
    bin_cutoff = 50
    permutations = 1000
    min_group_size = 5

    DEBUG = args.debug
    DEBUG_PLOTS_DIR = args.debug_plots_dir

    sys.stdout.write("- - - - - - - - - - - - - - - - - - - - -\n")
    sys.stdout.write("|\t\tDUDes %s\t\t|\n" % version)
    sys.stdout.write("- - - - - - - - - - - - - - - - - - - - -\n")
    sys.stdout.write("Output prefix = %s\n" % args.output_prefix)
    sys.stdout.write("SAM (format) = %s (%s)\n" % (args.sam_file, args.sam_format))
    sys.stdout.write("Database = %s\n" % args.database_file)
    sys.stdout.write("Threads = %d\n" % args.threads)
    sys.stdout.write(
        "TaxID Start/Last rank = %d/%s\n" % (args.taxid_start, args.last_rank)
    )
    sys.stdout.write("Max. Read Matches = %f\n" % args.max_read_matches)
    sys.stdout.write("Min. Ref. Matches = %f\n" % args.min_reference_matches)
    sys.stdout.write("Bin size = %f\n" % args.bin_size)
    sys.stdout.write("- - - - - - - - - - - - - - - - - - - - -\n")

    # Start pool before all to not account for memory several times (on some linux systems)
    if threads > 1:
        pool = mp.Pool(threads)
    else:
        pool = []

    # Calculate min. matches
    if args.min_reference_matches < 1:
        min_reference_matches_percent = args.min_reference_matches
        min_reference_matches_number = 1
    else:
        min_reference_matches_percent = 0
        min_reference_matches_number = int(args.min_reference_matches)

    # Load database file
    sys.stdout.write("Loading database file ...")
    tx = time.time()
    npzfile = np.load(args.database_file, allow_pickle=True)
    reference_mode = npzfile["reference_mode"].item()
    tax = np.array(npzfile["refid_nodes"])
    refids_lookup = npzfile["refids_lookup"].item()
    names = npzfile["names"].item()
    db_version = npzfile["version"].item()
    sys.stdout.write("\nReference mode = %s" % reference_mode)
    sys.stdout.write("\nDatabase version = %s\n" % db_version)
    sys.stdout.write(" Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")

    # Load sam file
    sys.stdout.write("Loading sam file ...")
    tx = time.time()
    if reference_mode == "up":
        pep_npzfile = np.load(args.sam_file, allow_pickle=True)
        sam = pep_npzfile["reads"]
        ref = pep_npzfile["reference_lengths"]
    else:
        sam, ref = parse_sam(
            args.sam_file, args.sam_format, refids_lookup, reference_mode, threads
        )
    sys.stdout.write(" Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")

    # Create objects
    smap = SMap(sam)
    refs = Refs(ref)
    ttree = TTree(tax)
    ident = Ident()
    names = Names(names)

    # Total matches before filtering
    total_matches_all = float(smap.getSize())

    # Filter nodes by last rank (don't limit if strain is selected)
    if args.last_rank in fixed_ranks and args.last_rank != "strain":
        fixed_ranks = fixed_ranks[: fixed_ranks.index(args.last_rank) + 1]
        ttree = ttree.getSubSet(
            ttree.getCol("RankID") <= fixed_ranks.index(args.last_rank)
        )

    # Filter references/matches without taxonomic entries (and the taxonomic entries based only on used references)
    sys.stdout.write("Filtering references/matches ...")
    tx = time.time()
    smap, refs, ttree = filterRefTax(smap, refs, ttree)
    sys.stdout.write(" Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")

    # Filter max matches (-m)
    if args.max_read_matches > 0:
        sys.stdout.write("Filtering max. read matches ...")
        tx = time.time()
        smap = filterMaxMatches(smap, args.max_read_matches, total_matches_all)
        sys.stdout.write(" Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")

    # Total matches after filtering
    total_matches_filter = float(smap.getSize())

    # Identify leaf nodes (for multiple testing correction)
    sys.stdout.write("Identifying leaf nodes ...")
    tx = time.time()
    ttree_order = np.lexsort((ttree.getCol("RankID"), ttree.getCol("RefID")))
    ttree_index = group_max(ttree.getCol("RefID"), ttree.getCol("RankID"), ttree_order)
    ttree.setLeafs(ttree.getSubSet(ttree_order).getSubSet(ttree_index))
    sys.stdout.write(" Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")

    # Calculate bin size (at least 75% of the reference sequences lens should be higher than the bin_size)
    if args.bin_size < 1:
        bin_size = int(
            np.round(np.percentile(refs.getCol("SeqLen"), args.bin_size * 100))
        )
        sys.stdout.write("Calculated bin size = %d\n" % bin_size)
    else:
        bin_size = int(args.bin_size)

    # Count iterations
    iter = 0

    # Main while
    while True:
        sys.stdout.write("Performing iteration " + str(iter + 1) + " ...")
        sys.stdout.flush()
        tx = time.time()
        debug_plot_dir_n = os.path.join(DEBUG_PLOTS_DIR, f"iter_{iter+1}")
        if DEBUG_PLOTS_DIR:
            os.makedirs(debug_plot_dir_n, exist_ok=True)
        iter_ident = Ident()
        # Save species matches to strain search
        species_matches = dict()
        # save indirect matches together (because they can overlap) to remove at once on the end of the loop
        removeMatches = np.zeros(smap.getSize(), dtype=bool)

        # When starting from any node but the root, add a temporary node with the chosen taxid_start
        if taxid_start > 1:
            iter_ident.add(
                Ident(
                    [
                        iter,
                        taxid_start,
                        1,
                        np.unique(
                            ttree.getSubSet(
                                ttree.getCol("TaxID") == taxid_start
                            ).getCol("RankID")
                        ),
                        0,
                        0,
                        0,
                        0,
                    ]
                )
            )

        if args.last_rank == "strain":
            stop_rank = "species"
        else:
            stop_rank = args.last_rank

        printDebug(
            DEBUG,
            "--------------------------------------------------------------------",
        )
        printDebug(DEBUG, "#### Iteration %d -> %s" % (iter, stop_rank))
        printDebug(DEBUG, iter_ident.ident)

        # First loop until last rank
        iter_ident.add(
            treeIter(
                pool, iter, smap, ttree, refs, taxid_start, stop_rank, debug_plot_dir_n
            )
        )

        if iter_ident.getSize():
            # Loop on current leafs of the identification
            for temp_leaf in iter_ident.getLeafs(1):
                directMatches = findDirectMatches(smap, ttree, temp_leaf["TaxID"])
                if directMatches.any():
                    indirectMatches = findIndirectMatches(smap, directMatches)
                    iter_ident.setMatchScoreSum(
                        iter,
                        temp_leaf["TaxID"],
                        np.sum(smap.getSubSet(directMatches).getCol("MatchScore")),
                    )
                    printDebug(
                        DEBUG,
                        "Leaf: %s - %d direct matches - %d indirect matches"
                        % (
                            temp_leaf["TaxID"],
                            np.sum(directMatches),
                            np.sum(indirectMatches),
                        ),
                    )
                    removeMatches = removeMatches + directMatches + indirectMatches
                    if (
                        temp_leaf["RankID"] == getIndexRank("species")
                        and args.last_rank == "strain"
                    ):
                        species_matches[temp_leaf["TaxID"]] = directMatches

            # Strain search
            if species_matches and args.last_rank == "strain":
                for species_taxid, dMatches in species_matches.items():
                    printDebug(
                        DEBUG, "#### Strain search - species taxid %d" % species_taxid
                    )
                    strains, totalMatchScoreSum = strainIdent(
                        pool,
                        iter,
                        smap.getSubSet(dMatches),
                        ttree,
                        refs,
                        species_taxid,
                        debug_plot_dir_n,
                    )
                    if strains.getSize():
                        # Reset MatchScoreSum from species removing the total used in the strains
                        iter_ident.setMatchScoreSum(
                            iter,
                            species_taxid,
                            iter_ident.getSubSet(
                                np.logical_and(
                                    iter_ident.getCol("Iter") == iter,
                                    iter_ident.getCol("TaxID") == species_taxid,
                                )
                            ).getCol("MatchScoreSum")
                            - totalMatchScoreSum,
                        )
                        iter_ident.add(strains)  # Add strains to iteration
        else:
            sys.stdout.write(" No more possible identifications\n")
            break

        sys.stdout.write(" Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")
        if taxid_start > 1 and iter_ident.getSize() == 1:
            break
        ident.add(iter_ident)
        smap = smap.getSubSet(~removeMatches)
        if smap.getSize() == 0:
            break
        iter += 1
    # end while True

    if ident.getSize():
        # Calculate abundances
        total_abundance_norm, ident = calcAbundance(ident, ttree, refs, taxid_start)

        printDebug(DEBUG, ident.ident)
        printDebug(
            DEBUG,
            "\n".join(
                ident.printTree(
                    taxid_start, total_matches_filter, total_abundance_norm, names
                )
            ),
        )

        # Representatives from the identification for output
        rep = Rep()
        rep.setRepresentatives(ident, args.last_rank)

        # Strain output
        if args.last_rank == "strain":
            if args.output_prefix:
                sys.stdout = open(args.output_prefix + "_strains.out", "w")
            printStrain(
                rep.getSubSet(rep.getCol("RankID") >= getIndexRank("species")),
                total_abundance_norm,
                names,
            )
            if args.output_prefix:
                sys.stdout.close()

        # Filter only best strain
        if args.last_rank == "strain":
            o, i = group_max(rep.getCol("ParentTaxID"), rep.getCol("Abundance"))
            rep_order = rep.getSubSet(o)
            rep_filter = rep_order.getSubSet(
                np.logical_or(rep_order.getCol("RankID") <= getIndexRank("species"), i)
            )
            rep_filter_sort = rep_filter.getSubSet(
                np.lexsort(
                    (-rep_filter.getCol("Abundance"), rep_filter.getCol("RankID"))
                )
            )
        else:
            rep_filter_sort = rep

        printDebug(DEBUG, rep.rep)
        printDebug(DEBUG, rep_filter_sort.rep)

        # Default output
        if args.output_prefix:
            sys.stdout = open(args.output_prefix + ".out", "w")
        printBioBoxes(
            rep_filter_sort,
            total_abundance_norm,
            names,
            args.sam_file,
            args.database_file,
        )
        if args.output_prefix:
            sys.stdout.close()
    else:
        sys.stdout.write(" No identifications")

    # Close threads
    if threads > 1:
        pool.close()
        pool.join()

    printDebug(DEBUG, args)

    if args.output_prefix:
        sys.stdout = sys.__stdout__
    sys.stdout.write("\n- - - - - - - - - - - - - - - - - - - - -\n")
    sys.stdout.write(
        "Total elapsed time: " + str(time.time() - total_tx) + " seconds\n"
    )


def strainIdent(pool, iter, smap_species, ttree, refs, taxid_species, debug_plot_dir):
    strains = Ident()
    totalMatchScoreSum = 0
    while True:
        iter_strain = treeIter(
            pool,
            iter,
            smap_species,
            ttree,
            refs,
            taxid_species,
            "strain",
            debug_plot_dir,
        )
        printDebug(DEBUG, iter_strain.ident)
        if iter_strain.getSize():
            # save indirect matches together (because they can overlap)
            removeMatches = np.zeros(smap_species.getSize(), dtype=bool)
            # loop on the identified strains
            # sort it by rankid desc. so higher nodes are calculated first (to account for the "intra"-strains)
            for temp_leaf_strain in iter_strain.getSubSet(
                np.argsort(iter_strain.getCol("RankID"))[::-1]
            ):
                directMatches = findDirectMatches(
                    smap_species, ttree, temp_leaf_strain["TaxID"]
                )
                indirectMatches = findIndirectMatches(smap_species, directMatches)
                mssum = sum(smap_species.getSubSet(directMatches).getCol("MatchScore"))
                # remove matches that were directly used in child nodes (calculated first)
                child_mssum = sum(
                    iter_strain.getSubSet(
                        iter_strain.getCol("ParentTaxID") == temp_leaf_strain["TaxID"]
                    ).getCol("MatchScoreSum")
                )
                iter_strain.setMatchScoreSum(
                    iter, temp_leaf_strain["TaxID"], mssum - child_mssum
                )
                removeMatches = removeMatches + directMatches + indirectMatches
                printDebug(
                    DEBUG,
                    "Leaf: %s - %d direct matches - %d indirect matches"
                    % (
                        temp_leaf_strain["TaxID"],
                        np.sum(directMatches),
                        np.sum(indirectMatches),
                    ),
                )

            strains.add(iter_strain)  # Add strains to iteration
            smap_species = smap_species.getSubSet(~removeMatches)
            totalMatchScoreSum += sum(iter_strain.getCol("MatchScoreSum"))
            if smap_species.getSize() == 0:
                break
        else:
            break

    return strains, totalMatchScoreSum


def treeIter(pool, iter, smap, ttree, refs, taxid_start, stop_rank, debug_plot_dir):
    """

    :param pool: multiprocessing Pool or empty list
    :param iter: current iteration number
    :param smap: SMap object
    :param ttree: TTree object
    :param refs: Refs object
    :param taxid_start: taxid to start from
    :param stop_rank: most specific rank to check
    :param debug_plot_dir: directory to write the debug plots to
    :return: Ident object
    """

    iter_ident = Ident()
    stack = [taxid_start]

    if taxid_start == 1:
        sub_smap = smap
    else:
        sub_smap = smap.getSubSet(findDirectMatches(smap, ttree, taxid_start))

    if sub_smap.getSize():
        bins = Bins(sub_smap, bin_size)

        while stack:
            parent_taxid = stack.pop()

            # Segment children by rank
            children_ranks = np.unique(
                ttree.getSubSet(ttree.getCol("ParentTaxID") == parent_taxid).getCol(
                    "RankID"
                )
            )

            if children_ranks.size == 0:  # No children, reached a leaf node
                continue
            elif children_ranks.size == 1:  # Only one rankid
                rankid = int(children_ranks)
            else:  # If there are more than one rankid among the children, choose one with higher matchscore sum
                rank_matches = 0
                for cr in children_ranks:
                    # sum binscores from this rankid (faster than get from the smap)
                    m = np.sum(
                        bins.getBins(
                            ttree.getSubSet(
                                np.logical_and(
                                    ttree.getCol("ParentTaxID") == parent_taxid,
                                    ttree.getCol("RankID") == cr,
                                )
                            ).getCol("RefID")
                        )[0]
                    )
                    if m >= rank_matches:
                        rankid = cr
                        rank_matches = m

            # Exit on last rank
            if rankid > getIndexRank(stop_rank) and stop_rank != "strain":
                continue

            # Select a subset of results (children of the parent node with the same rank)
            sub_ttree = ttree.getSubSet(
                np.logical_and(
                    ttree.getCol("ParentTaxID") == parent_taxid,
                    ttree.getCol("RankID") == rankid,
                )
            )

            # nodes tested on this rank
            nodes = np.unique(sub_ttree.getCol("TaxID"))

            # Creates groups in advance to minimize overhead
            groups = dict()
            matches = dict()
            for n in nodes:
                groups[n], matches[n] = bins.getBins(
                    sub_ttree.getSubSet(sub_ttree.getCol("TaxID") == n).getCol("RefID")
                )

            # Sort nodes by number of matches - USUALLY the ones with more matches are going to be choosen, so evaluating them first allow to skip the following nodes
            nodes = np.array(
                [
                    n[0]
                    for n in sorted(
                        list(matches.items()), key=lambda x: x[1], reverse=True
                    )
                ]
            )

            # initialize variables pvalues and critical values
            pvs = np.ones((nodes.size, nodes.size), float)
            cvs = np.zeros(nodes.size, float)

            # Control if the matrix has real rejections (real compared nodes or just filtered by parameters)
            real_rej = False

            # CRITICAL VALUES
            for pi, c_taxid in enumerate(nodes):
                if nodes.size > 1:
                    # Search how many leaf nodes are on the specific taxid
                    m_hypothesis = np.unique(
                        ttree.getLeafs()
                        .getSubSet(
                            np.in1d(
                                ttree.getLeafs().getCol("RefID"),
                                sub_ttree.getSubSet(
                                    sub_ttree.getCol("TaxID") == c_taxid
                                ).getCol("RefID"),
                            )
                        )
                        .getCol("TaxID")
                    ).size
                    # Bonferonni and meinshausen correction
                    cvs[pi] = (thr_alpha / float(len(nodes) - 1)) * (
                        m_hypothesis / float(ttree.getLeafs().getSize())
                    )
                else:
                    cvs[pi] = 1

            # Rejections/Identifications
            rjs = np.zeros((nodes.size, nodes.size), int)
            # 0 = not yet evaluated - diagonal
            # -1 = rejected/identified (p-val<=cv)
            # 1 = not rejected/not identified (p-val>cv)
            # Max number of possible rejections
            max_rej = (len(nodes) * 2) - 1

            for n in range(len(nodes)):

                # Store the instances of the pool
                pool_res = {}

                # Check based on the already calculated values if it is still possible to identify the node, otherwise skip it
                temp_rej = np.sum(rjs[n, :] == 1) + np.sum(rjs[:, n] == -1)
                if temp_rej > max_rej:
                    continue

                # iterate over all possible pairs of current node index n and all node indices (from nodes9
                for i, j in [(n, i) for i in range(len(nodes))] + [
                    (i, n) for i in range(len(nodes))
                ]:
                    if not rjs[
                        i, j
                    ]:  # Do not re-evaluate (Evaluate if entry was 0 (not yet evaluated))
                        if i == j and nodes.size > 1:
                            pvs[i, j] = 0
                            continue
                        elif (
                            (matches[nodes[i]] / total_matches_filter)
                            < min_reference_matches_percent
                            or matches[nodes[i]] < min_reference_matches_number
                            or groups[nodes[i]].size < min_group_size
                        ):
                            # Not enough matches or groups on the control group - Skip
                            pvs[i, j] = 1
                            rjs[i, j] = 1
                            continue
                        elif (
                            (matches[nodes[j]] / total_matches_filter)
                            < min_reference_matches_percent
                            or matches[nodes[j]] < min_reference_matches_number
                            or groups[nodes[j]].size < min_group_size
                        ):
                            # Not enough matches on the treatment group. This node is identified because the control has enough matches.
                            real_rej = True
                            pvs[i, j] = 0
                            rjs[i, j] = -1
                            continue
                        elif nodes.size == 1:
                            # No treatment group (only child with this rank on the branch). This node is identified because the control has enough matches.
                            real_rej = True
                            pvs[i, :] = 0
                            rjs[i, :] = -1
                            continue
                        else:
                            real_rej = True
                            if threads > 1:
                                pool_res[(i, j)] = pool.apply_async(
                                    perm_pval,
                                    args=(
                                        groups[nodes[i]],
                                        groups[nodes[j]],
                                        cvs[i],
                                        nodes[i],
                                        nodes[j],
                                        debug_plot_dir,
                                    ),
                                )
                            else:
                                pv = perm_pval(
                                    groups[nodes[i]],
                                    groups[nodes[j]],
                                    cvs[i],
                                    nodes[i],
                                    nodes[j],
                                    debug_plot_dir,
                                )
                                pvs[i, j] = pv  # store p-value
                                if pv > cvs[i]:
                                    rjs[i, j] = 1
                                else:
                                    rjs[i, j] = -1
                # end for i,j in [(n,i) for i in range(len(nodes))]+[(i,n) for i in range(len(nodes))]

                if threads > 1:
                    if pool_res:
                        for (i, j), res in list(pool_res.items()):
                            pv = res.get()
                            pvs[i, j] = pv  # store p-value
                            if pv > cvs[i]:
                                rjs[i, j] = 1
                            else:
                                rjs[i, j] = -1

                # After evaluating the node check if rejections are lower than the max and update
                temp_rej = np.sum(rjs[n, :] == 1) + np.sum(rjs[:, n] == -1)
                if temp_rej < max_rej:
                    max_rej = temp_rej

            # end for n in range(len(nodes))

            # If there are real rejections on the matrix
            if real_rej:
                # Rejections
                s_rows = np.sum(
                    rjs == 1, axis=1
                )  # Sum of nodes that could not be rejected against the others
                s_cols = np.sum(
                    rjs == -1, axis=0
                )  # Sum of nodes that could not reject the others
                comb_sum = s_rows + s_cols  # Combined sum of the rejections
                min_node = np.min(comb_sum)  # Set minimum as the best node(s)
                id_rej_rows = np.where(comb_sum == min_node)[0]

                printDebug(DEBUG, "-- PARENT TAXID %d" % parent_taxid)
                if DEBUG:
                    for i, pv in enumerate(pvs):
                        printDebug(
                            DEBUG,
                            ("\t%d\t%s %d %f")
                            % (
                                nodes[i],
                                pv,
                                np.sum(rjs[i, :] == 1) + np.sum(rjs[:, i] == -1),
                                cvs[i],
                            ),
                        )
                    printDebug(DEBUG, rjs)
                    printDebug(
                        DEBUG,
                        "\tIdentified nodes: %s "
                        % [nodes[idrej] for idrej in id_rej_rows],
                    )

                # Add rejections -> identification
                if id_rej_rows.size > 0:
                    for r in id_rej_rows:
                        taxid = nodes[r]
                        iter_ident.add(
                            Ident(
                                [
                                    iter,
                                    taxid,
                                    parent_taxid,
                                    rankid,
                                    matches[taxid],
                                    0,
                                    0,
                                    0,
                                ],
                                {(iter, taxid): pvs[r]},
                            )
                        )
                        stack.append(taxid)
        # end if real_rej
    # end while stack
    return iter_ident


def perm_pval(c, t, cv, taxon_id_control, taxon_id_treatment, debug_plot_dir):
    """
    calculate permutation based p-value of control group versus treatment group.

    :param c: array of control group bin scores
    :param t: array of treatment group bin scores
    :param cv: critical value
    :param taxon_id_control: taxon id of the control
    :param taxon_id_treatment: taxon id of the treatment
    :param debug_plot_dir: path to directory to save debug plots in
    :return:
    """
    # Normalize the groups (0-1)
    # mx = np.max(np.concatenate([c,t]))
    norm_c = np.sort(c)[::-1]  # /float(mx)
    norm_t = np.sort(t)[::-1]  # /float(mx)

    # Percentile based only on non-zero values
    ###cutoff = sum(norm_c>=np.percentile(norm_c[norm_c>0],100-bin_cutoff))
    cutoff = np.ceil(norm_c[norm_c > 0].size * (float(bin_cutoff) / 100))

    # If the number of groups is greater than the treatment, limit the size on the treatment group (only non-zeros)
    if norm_t[norm_t > 0].size < cutoff:
        cutoff = norm_t[norm_t > 0].size

    # Subset control and treatment groups based on the cutoff
    norm_c_sub = norm_c[: int(cutoff)]
    norm_t_sub = norm_t[: int(cutoff)]

    # Permutation test
    diff_obs = np.mean(norm_c_sub) - np.mean(norm_t_sub)
    combined = np.concatenate([norm_c_sub, norm_t_sub])

    diff_random = np.zeros(permutations)
    np.random.seed()  # re-seed for the threading
    for n in range(permutations):
        np.random.shuffle(combined)
        diff_random[n] = np.mean(combined[: norm_c_sub.size]) - np.mean(
            combined[norm_c_sub.size :]
        )

    pv = sum(np.greater_equal(diff_random, diff_obs)) / float(permutations)

    if DEBUG_PLOTS_DIR:
        plot_bin_scores(
            norm_c,
            norm_t,
            cutoff,
            diff_random,
            diff_obs,
            cv,
            taxon_id_control,
            taxon_id_treatment,
            debug_plot_dir,
        )

    return pv


def plot_bin_scores(
    norm_c,
    norm_t,
    cutoff,
    diff_random,
    diff_obs,
    cv,
    taxon_id_control,
    taxon_id_treatment,
    plot_dir,
):
    """
    Plot figure with two graphs:
    1. norm_c in blue, norm_t in red, cutoff -1 as vertical line
    2. histogrmm of diff_random, with diff_obs and cv as vertical lines

    :param norm_c:
    :param norm_t:
    :param cutoff:
    :param diff_random:
    :param diff_obs:
    :param cv:
    :param taxon_id_control:
    :param taxon_id_treatment:
    :param plot_dir:
    :return:
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(nrows=1, ncols=2, sharex=False, sharey=False)
    ax[0].plot(norm_c, lw=2.5, color="blue", label="A")
    ax[0].plot(norm_t, lw=2.5, color="red", label="B")
    ax[0].legend(fontsize=20)
    ax[0].axvline(cutoff - 1, color="black", lw=3)
    # ax[0].text(cutoff-1+5, 100, int(cutoff), fontsize=20)
    # ax[0].set_xticks([])
    ax[0].set_yticks([])
    ax[0].set_xlabel("bins", fontsize=20)
    ax[0].set_ylabel("bin score", fontsize=20)
    ax[1].hist(diff_random, 50, alpha=0.8, color="gray")
    ax[1].axvline(diff_obs, color="green", lw=3)
    ax[1].axvline(1 - np.percentile(diff_random, cv), color="red", lw=3, ls="--")
    ax[1].set_xticks([])
    ax[1].set_yticks([])
    xlim = ax[1].get_xlim()
    ax[1].set_xlim(xlim[0] - 100, xlim[1] + 100)
    fig.tight_layout()
    # build filename
    plt_fname = os.path.join(
        plot_dir, f"{taxon_id_control}_vs_{taxon_id_treatment}.svg"
    )
    fig.savefig(plt_fname)
    plt.close(fig)


def findDirectMatches(smap, ttree, taxid):
    """
    get indices of smap that are direct matches of taxid

    :param smap:
    :param ttree:
    :param taxid:
    :return:
    """

    # ttree - References from the identified taxid
    rej_ref_ttree = ttree.getSubSet(np.in1d(ttree.getCol("TaxID"), taxid)).getCol(
        "RefID"
    )
    # INDEX smap - References from the identified taxid --- DIRECT MATCHES
    direct_matches_smap_idx = np.in1d(smap.getCol("RefID"), rej_ref_ttree)
    return direct_matches_smap_idx


def findIndirectMatches(smap, direct_matches_smap_idx):
    # INDEX smap - All reads that have at least one match against one of the identified references
    rej_reads_smap_idx = np.in1d(
        smap.getCol("ReadID"), smap.getSubSet(direct_matches_smap_idx).getCol("ReadID")
    )

    # Get best match from the direct matches (because sometimes the same read can match several references)
    direct_matches = smap.getSubSet(direct_matches_smap_idx)
    order_dm, index_dm = group_max(
        direct_matches.getCol("ReadID"), direct_matches.getCol("MatchScore")
    )
    max_match_score_direct_matches = defaultdict(
        int,
        list(
            zip(
                np.unique(direct_matches.getCol("ReadID")),
                direct_matches.getSubSet(order_dm)
                .getSubSet(index_dm)
                .getCol("MatchScore"),
            )
        ),
    )

    # Matches from the reads used in the identification but assigned to other references
    indirect_matches_smap_idx = np.logical_and(
        rej_reads_smap_idx, ~direct_matches_smap_idx
    )
    indirect_matches = smap.getSubSet(indirect_matches_smap_idx)
    # Loop through all smap elements but just analyse the ones that are marked as indirect (non-intuitive but faster)
    for ui, udm in enumerate(indirect_matches_smap_idx):
        if udm:
            u = indirect_matches.next()
            if u["MatchScore"] > max_match_score_direct_matches[u["ReadID"]]:
                indirect_matches_smap_idx[ui] = False

    return indirect_matches_smap_idx


def filterRefTax(smap, refs, ttree):
    """
    Remove all references from smap and refs that were not found in refids_lookup (RefID == -1).
    Filter ttree to contain only those RefID entries which are also in the filtered refs.

    :param smap:
    :param refs:
    :param ttree:
    :return:
    """

    # References marked with -1 are the ones not found in the refids_lookup (on parse_sam)
    rem_refs = refs.getCol("RefID") == -1
    rem_matches = smap.getCol("RefID") == -1

    sys.stdout.write(
        "\n\t%d (%.2f%%) sequences without taxonomic entries on the database removed"
        % (np.sum(rem_refs), (np.sum(rem_refs) / float(refs.getSize())) * 100)
    )
    sys.stdout.write(
        "\n\t  %d (%.2f%%) matches removed\n"
        % (np.sum(rem_matches), (np.sum(rem_matches) / float(smap.getSize())) * 100)
    )
    # Filter matches
    fsmap = smap.getSubSet(~rem_matches)
    # Filter references
    frefs = refs.getSubSet(~rem_refs)
    # Filter tree
    fttree = ttree.getSubSet(np.in1d(ttree.getCol("RefID"), frefs.getCol("RefID")))

    return fsmap, frefs, fttree


def filterMaxMatches(smap, max_read_matches, total_matches_all):
    """
    Remove matches with more than this number/percentile of matches (0: off / 0-1: percentile / >=1: match count)

    :param smap:
    :param max_read_matches:
    :param total_matches_all:
    :return:
    """
    if max_read_matches < 1:
        mm = int(np.percentile(list(match_counter.values()), max_read_matches * 100))
    else:
        mm = int(max_read_matches)

    match_counter = np.bincount(smap.getCol("ReadID"))
    # ReadIDs with more than mm counts
    rm_reads = np.where(match_counter > mm)[0]
    if rm_reads.any():
        fsmap = smap.getSubSet(~np.in1d(smap.getCol("ReadID"), rm_reads))
        matches_rem = smap.getSize() - fsmap.getSize()
        sys.stdout.write(
            "\n\t%d (%.2f%%) reads with more than %d matches found"
            % (
                len(rm_reads),
                (len(rm_reads) / float(np.unique(smap.getCol("ReadID")).size)) * 100,
                mm,
            )
        )
        sys.stdout.write(
            "\n\t  %d (%.2f%%) matches removed\n"
            % (matches_rem, (matches_rem / total_matches_all) * 100)
        )
        return fsmap
    else:
        sys.stdout.write("\n\tNo reads with more than %d matches found\n" % mm)
        return smap


def calcAbundance(ident: Ident, ttree: TTree, refs: Refs, taxid_start):
    """
    Calculate abundance of identifications.

    Calculates individual and cumulative abundance of identifications.
    Calculates individual abundance for each identifications with a MatchScoreSum > 0.
    Individual abundance is an int and calculated by dividing the identified taxons "MatchScoreSum" by the length of all reference
    sequences belonging to that taxon, multiplied by 10^9 and rounded.
    Writes individual abundance into Ident object column "Abundance".
    Calculates cumulative abundance summing the individual abundance of an identification and all its child taxon IDs.
    Writes cumulative abundance into Ident object column "CumulativeAbundance".

    total_abundance_norm: sum of cumulative abundances of identifications with taxid == taxid_start
    special case: if taxid_start == 1, additionaly sums identifications with parent taxid == 1

    Args:
        ident: Ident object
        ttree: TTree object, containing refid_nodes
        refs: Refs object, containing numeric RefID and its sequence length
        taxid_start: taxonomic Id used to start the analysis

    Returns:
        tuple(total_abundance_norm, Ident)
    """
    # Single abundance
    for id_ in ident.getSubSet(ident.getCol("MatchScoreSum") > 0):
        lens = sum(
            # get array of sequence length for all references belonging to identified taxa
            refs.getSubSet(
                np.in1d(
                    refs.getCol("RefID"),
                    # get refids from ttree for taxids in ids
                    ttree.getSubSet(ttree.getCol("TaxID") == id_["TaxID"]).getCol(
                        "RefID"
                    ),
                )
            ).getCol("SeqLen")
        )
        abundance = id_["MatchScoreSum"] / float(lens)
        ident.setAbundance(
            id_["Iter"], id_["TaxID"], int(np.round(abundance * 10**9))
        )

    # Cumulative abundance
    total_abundance_norm = 0
    for iter in range(ident.getIter() + 1):
        ident_iter = ident.getSubSet(ident.getCol("Iter") == iter)
        for id_ in ident_iter:
            cumulative_abundance = 0
            stack = [id_["TaxID"]]
            while stack:
                node = stack.pop()
                if node:
                    # subset of identifications with taxid == node
                    node_d = ident_iter.getSubSet(ident_iter.getCol("TaxID") == node)
                    # Sum all the single abundances (MatchScoreSum > 0)
                    if node_d.getCol("Abundance") > 0:
                        cumulative_abundance += node_d.getCol("Abundance")
                    children = ident_iter.getSubSet(
                        ident_iter.getCol("ParentTaxID") == node_d.getCol("TaxID")
                    ).getCol("TaxID")
                    if children.any():
                        stack.extend(children)
            # Set abundance
            ident.setCumulativeAbundance(
                id_["Iter"], id_["TaxID"], cumulative_abundance
            )
            # Save total matches for the start node (to normalize)
            if id_["TaxID"] == taxid_start or (
                taxid_start == 1 and id_["ParentTaxID"] == taxid_start
            ):
                total_abundance_norm = total_abundance_norm + cumulative_abundance
    return total_abundance_norm, ident


def printStrain(rep, total_abundance_norm, names):
    # Strains
    for r in rep:
        if r["RankID"] == getIndexRank("species"):
            print(
                (
                    "%s\t%d\t%s\t%.6f"
                    % (
                        getNameRank(r["RankID"]),
                        r["TaxID"],
                        names.getName(r["TaxID"]),
                        np.round(
                            (r["Abundance"] / float(total_abundance_norm)) * 100, 5
                        ),
                    )
                )
            )
            for r_strains in rep.getSubSet(rep.getCol("ParentTaxID") == r["TaxID"]):
                print(
                    (
                        "\t%s\t%d\t%s\t%.6f"
                        % (
                            getNameRank(r_strains["RankID"]),
                            r_strains["TaxID"],
                            names.getName(r_strains["TaxID"]),
                            np.round(
                                (r_strains["Abundance"] / float(total_abundance_norm))
                                * 100,
                                5,
                            ),
                        )
                    )
                )


def printBioBoxes(rep, total_abundance_norm, names, sam_file, database_file):
    # BioBoxes output
    print("# Taxonomic Profiling Output")
    print(("@SampleID:%s" % sam_file))
    print("@Version:0.9.3")
    print(("@Ranks:%s" % "".join([r + "|" for r in fixed_ranks])[:-1]))
    print(("@TaxonomyID:%s" % database_file))
    print("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE")
    for r in rep:
        taxidpath = str(r["TaxID"])
        namepath = names.getName(r["TaxID"])
        taxid = r["TaxID"]
        parent_taxid = r["ParentTaxID"]
        for rank in reversed(list(range(r["RankID"]))):
            f = rep.getSubSet(
                np.logical_and(
                    rep.getCol("TaxID") == parent_taxid, rep.getCol("RankID") == rank
                )
            ).getCol("ParentTaxID")
            if f:
                taxidpath = str(int(parent_taxid)) + "|" + taxidpath
                namepath = names.getName(int(parent_taxid)) + "|" + namepath
                parent_taxid = f
            else:
                taxidpath = "|" + taxidpath
                namepath = "|" + namepath
            if parent_taxid == 1:
                break

        print(
            (
                "%d\t%s\t%s\t%s\t%.6f"
                % (
                    r["TaxID"],
                    getNameRank(r["RankID"]),
                    taxidpath,
                    namepath,
                    np.round((r["Abundance"] / float(total_abundance_norm)) * 100, 5),
                )
            )
        )