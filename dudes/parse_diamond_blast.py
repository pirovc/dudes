from collections import defaultdict
from pathlib import Path
from typing import Dict, Union

import numpy as np
import pandas as pd
from pandas import DataFrame


def parse_uniprot_accession(raw_accession):
    return raw_accession.split("|")[1]


def read_blast_tsv(f: Union[str, Path]) -> DataFrame:
    """Read blast output file into dataframe.

    The sseqid column is converted to uniprot accession format ("sp|wasdef|bar" -> "wasdef")

    The returned dataframe has the following columns:
    - qseqid: object, Query Seq - id
    - sseqid: object, Subject Seq - id
    - slen: int, Subject sequence length
    - sstart: int, Start of alignment in query
    - cigar: object, CIGAR string
    - pident: float, Percentage of identical matches
    - evalue: float, expected value of the hit quantifies the number of alignments of similar or better quality that you
     expect to find searching this query against a database of random sequences the same size as the actual target
     database.

    Args:
        f: path to blast output file, required columns:
            0: qseqid,
            1: sseqid,
            2: slen,
            3: sstart,
            7: evalue

    Returns:
        pandas DataFrame
    """
    dtypes = {
        "qseqid": object,
        "sseqid": object,
        "slen": int,
        "sstart": int,
        "evalue": float,
    }
    df = pd.read_table(
        f,
        usecols=list(range(0, len(dtypes))),
        names=list(dtypes.keys()),
        # work around pandas ParserWarning if "sseqid" is in dtype and converters
        dtype={k: v for k, v in dtypes.items() if k != "sseqid"},
        converters={"sseqid": parse_uniprot_accession}
    )
    return df


def parse_reference_lengths(blast_df: pd.DataFrame, refid_lookup: Dict[str, int]) -> np.array:
    """Build array of reference sequence lengths.

    Args:
        blast_df: dataframe with required columns "sseqid" and "slen".
        refid_lookup: dictionary with reference accessions as keys and their dudes ID as values

    Returns:
        np.array with two columns. First column is the dudes ID of the reference, second column is the length of the
        reference sequence
    """

    df = blast_df.drop_duplicates(subset=["sseqid"]).copy()
    df.loc[:, "refid"] = df["sseqid"].apply(lambda x: refid_lookup.get(x, -1))
    df = df[df["refid"] != -1]
    return df[["refid", "slen"]].to_numpy()


def transform_blast_df_into_sam_array(blast_df: pd.DataFrame, refid_lookup: Dict[str, int]) -> np.array:
    """Build query array.

    Args:
        blast_df: dataframe with required columns "sseqid", "sstart", "cigar", "mismatch", and "qseqid".
        refid_lookup: dictionary with reference accessions as keys and their dudes ID as values

    Returns:
        np.array with 4 columns.
              - 0: int, 'RefID': matched reference ID in 'refid_lookup'
              - 1: int, 'MatchPosStart': start position of match in reference sequence
              - 2: int, 'MatchScore': score of the match
              - 3: int, 'ReadID': ID of the read
    """
    read_id_lookup = defaultdict(lambda: len(read_id_lookup))
    sam_series = blast_df.apply(lambda row: [
        refid_lookup.get(row["sseqid"], -1),
        row["sstart"],
        int(-10 * np.log10(row["evalue"])),
        read_id_lookup[row["qseqid"]]
    ], axis=1)
    return np.stack(sam_series)


def parse_custom_blast(custom_blast_file: str, refid_lookup: Dict[str, int], threads: int):
    """Read custom BLAST file into two arrays.

    Args:
        custom_blast_file: path to custom blast file
        refid_lookup: dictionary, key: source reference IDs, value: internal reference ID (int)
        threads: number of threads

    Returns:
        numpy.array, numpy.array:
            first array: 4 columns:
                - 0: int, 'RefID': matched reference ID in 'refid_lookup' dict
                - 1: int, 'MatchPosStart': start position of match in reference sequence
                - 2: int, 'MatchScore': score of the match
                - 3: int, 'ReadID': ID of the read
            second array: 2 columns
                - 0: int, dudes internal id for the reference accession
                - 1: int, length of the reference
    """
    blast_df = read_blast_tsv(custom_blast_file)
    reference_lengths = parse_reference_lengths(blast_df, refid_lookup)
    sam_array = transform_blast_df_into_sam_array(blast_df, refid_lookup)
    return sam_array, reference_lengths
