import re
import numpy as np
import multiprocessing as mp
from itertools import islice
from collections import defaultdict

from dudes.calculate_match_score import calculate_match_score

NM_REGEX = re.compile("NM:i:[0-9]*")


def parse_lines(lines):
    """

    :param lines:
    :return:
        numpy array, one row for each line that has a valid reference (not '*'), three columns
            col 0: int, dudes internal reference id, -1 if not in refids_lookup
            col 1: int, 1-based-leftmost-mapping-position of the read
            col 2: int, match score, calculated as follows:
              (matches + insertions) - (Edit distance to the reference + deletions + insertions)
              matches, insertions and deletions are read from the CIGAR string, the edit distance to the reference is
                read from the NM tag. Edit distance to the reference is defined as: Number of differences (mismatches
                plus inserted and deleted bases) between the sequence and reference
        list of str: read name from first field of each line
    """
    sam = np.zeros((len(lines), 3), dtype=int)
    reads = []
    c = 0
    for l in lines:
        fields = l.split("\t")
        # read, _, ref_acc, 1-based-leftmost-mapping-position, _. CIGAR-string, *_
        if fields[2] != "*":  # reference accesion
            edit_distance = int(NM_REGEX.search(l).group()[5:])
            sam[c, 0] = REFIDS_LOOKUP.get(
                REFID_REGEX.search(fields[2]).group(1), -1
            )  # refid code from lookup (gi or av) or -1
            sam[c, 1] = int(fields[3])  # POS
            # MATCH Score -> LEN - (NM + INDELS)
            sam[c, 2] = calculate_match_score(cigar_string=fields[5], edit_distance=edit_distance)
            reads.append(fields[0])
            c = c + 1
    return sam[:c, :], reads


def parse_lines_extended(lines):
    sam = np.zeros((len(lines), 3), dtype=int)
    reads = []
    c = 0
    for l in lines:
        fields = l.split("\t")
        if fields[2] != "*":
            cig = {"I": 0, "D": 0, "M": 0, "=": 0, "X": 0}
            for val, ci in re.findall(r"(\d+)([IDM=X])", fields[5]):
                cig[ci] += int(val)
            sam[c, 0] = REFIDS_LOOKUP.get(
                REFID_REGEX.search(fields[2]).group(1), -1
            )  # refid code from lookup (gi or av) or -1
            sam[c, 1] = int(fields[3])  # POS
            sam[c, 2] = (cig["X"] + cig["="] + cig["M"] + cig["I"]) - (
                cig["X"] + cig["D"] + cig["I"]
            )  # MATCHES -> LEN - (NM + INDELS)
            reads.append(fields[0])
            c = c + 1
    return sam[:c, :], reads


def parse_sam(sam_file, sam_format, rl, reference_mode, threads):
    """
    Read sam file into two arrays.

    :param sam_file: path to sam file
    :param sam_format: 'nm' or 'ex', SAM file format ['nm': sam file with standard
    cigar string plus NM flag (NM:i:[0-9]*) for mismatches count | 'ex': just the extended cigar string].
    :param rl: refids_lookup, dictionary, key: source reference IDs (as in sam), value: internal reference ID (int)
    :param reference_mode: ['gi', 'av']
    :param threads: number of threads
    :return: numpy.array, numpy.array:
        first array: 4 columns:
            - 0: int, 'RefID': matched reference ID in 'rl' dict
            - 1: int, 'MatchPosStart': start position of match in reference sequence
            - 2: int, 'MatchScore': score of the match
            - 3: int, 'ReadID': ID of the read
        second array: 2 columns
            - 0: int, dudes internal id for the reference accession
            - 1: int, LN attribute of reference
    """

    global REFIDS_LOOKUP
    global REFID_REGEX

    # Global refids lookup
    REFIDS_LOOKUP = rl

    if reference_mode == "gi":
        REFID_REGEX = re.compile(r"gi\|(\d*)")
    elif reference_mode == "up":
        REFID_REGEX = re.compile(r"\w{2}\|([^|]*)\|")
    else:
        REFID_REGEX = re.compile(r"([A-Z\d_.]*)")  # without ">" from fasta regex

    refs = []
    reads = defaultdict(lambda: len(reads))
    pool = mp.Pool(threads)

    file = open(sam_file, "r")
    while True:
        last_pos = file.tell()
        l = file.readline()
        fields = l.split("\t")
        if fields[0] == "@SQ":  # list of references on the database
            # example fields:
            # 0: @SQ
            # 1: SN:NC_010571.1
            # 2: LN:5957605
            # refs: list of lists, each sublist is a pair of
            # 1. dudes internal id for the reference accession, if the reference accession
            #     is not found in refids lookup: -1
            # 2. LN attribute: reference sequence length
            refs.append(
                [
                    REFIDS_LOOKUP.get(
                        REFID_REGEX.search(fields[1][3:]).group(1), -1
                    ),
                    int(fields[2][3:]),
                ]
            )
        elif l[0] == "@":
            continue  # other headers
        else:
            first_aln_line = l
            break
    file.seek(last_pos)  # return to the last iterated line (first sam result)

    n_lines = 1000
    res = []

    # Check for NM flag
    if sam_format == "nm" and not NM_REGEX.search(first_aln_line):
        print(
            "\n -- Warning: NM flag not found. Switching to extended CIGAR mode (-i 'ex')"
        )
        sam_format = "ex"

    if sam_format == "nm":
        while True:
            next_n_lines = list(islice(file, n_lines))  # select chunk of lines
            if not next_n_lines:
                break
            res.append(
                pool.apply_async(parse_lines, args=(next_n_lines,))
            )  # start the process in the pool
    elif sam_format == "ex":
        while True:
            next_n_lines = list(islice(file, n_lines))  # select chunk of lines
            if not next_n_lines:
                break
            res.append(
                pool.apply_async(parse_lines_extended, args=(next_n_lines,))
            )  # start the process in the pool

    sam_all = []
    for r in res:
        # assign unique id for each read (reads defaultdict)
        if r.get()[0].any():
            rd = np.array([reads[read] for read in r.get()[1]])  # READID
            sm = r.get()[0]  # sam data
            sam_all.append(
                np.insert(sm, sm[0, :].size, rd, axis=1)
            )  # append readid column in the data array and add to the results

    pool.close()
    pool.join()
    file.close()

    return np.concatenate(sam_all), np.array(refs)
