from collections import defaultdict
from os import PathLike
from typing import Union

import numpy as np
import pandas as pd
from pandas import DataFrame


def parse_uniprot_accession(raw_accession):
    return raw_accession.split("|")[1]


def read_blast_tsv(f: Union[str, PathLike[str]]) -> DataFrame:
    """Read blast output file into dataframe.

    The sseqid column is converted to uniprot accession format ("sp|wasdef|bar" -> "wasdef")

    The dataframe has the following columns:
    - qseqid: object, Query Seq - id
    - sseqid: object, Subject Seq - id
    - slen: int, Subject sequence length
    - sstart: int, Start of alignment in query
    - cigar: object, CIGAR string
    - pident: float, Percentage of identical matches

    Args:
        f: path to blast output file

    Returns:
        pandas DataFrame
    """
    dtypes = {
        "qseqid": object,
        "sseqid": object,
        "slen": int,
        "sstart": int,
        "cigar": object,
        "pident": float
    }
    df = pd.read_table(
        f,
        names=dtypes.keys(),
        dtype=dtypes,
        converters={"sseqid": parse_uniprot_accession}
    )
    return df


def parse_reference_lengths(blast_df: pd.DataFrame, refid_lookup: dict[str, int]) -> np.array:
    """Build array of reference sequence lengths.

    Args:
        blast_df: dataframe with required columns "sseqid" and "slen".
        refid_lookup: dictionary with reference accessions as keys and their dudes ID as values

    Returns:
        np.array with two columns. First column is the dudes ID of the reference, second column is the length of the
        reference sequence
    """

    df = blast_df.drop_duplicates(subset=["sseqid"])
    df["refid"] = df["sseqid"].apply(lambda x: refid_lookup.get(x, -1))
    df = df[df["refid"] != -1]
    return df[["refid", "slen"]].to_numpy()


def transform_blast_df_into_sam_array(blast_df: pd.DataFrame, refid_lookup: dict[str, int]) -> np.array:
    """Build query array.

    Args:
        blast_df: dataframe with required columns "sseqid" and "slen".
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
        refid_lookup[row["sseqid"]],
        row["sstart"],
        # compute_score(row["cigar"]),
        read_id_lookup[row["qseqid"]]
    ], axis=1)
    return np.stack(sam_series)
