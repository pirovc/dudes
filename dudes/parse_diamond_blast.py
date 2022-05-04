import numpy as np
import pandas as pd


def parse_uniprot_accession(raw_accession):
    return raw_accession.split("|")[1]


def read_blast_tsv(f):
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


def parse_blast_df_into_sam_array(blast_df: pd.DataFrame, refid_lookup: dict[str, int]) -> np.array:
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
    pass
