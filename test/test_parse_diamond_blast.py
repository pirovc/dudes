import numpy as np
import pandas as pd

from dudes.parse_diamond_blast import read_blast_tsv, parse_uniprot_accession, parse_reference_lengths


def test_read_into_dataframe(resource_dir):
    returned = read_blast_tsv(resource_dir / "diamond_blast-qseqid-sseqid-slen-sstart-cigar-pident.tsv")
    assert isinstance(returned, pd.DataFrame)


def test_parse_uniprot_accession():
    assert parse_uniprot_accession("sp|the_accession|asdfg") == "the_accession"


def test_parse_reference_lengths():
    df = pd.DataFrame(
        {
            "sseqid": ["abc", "abc", "def", "not_present"],
            "slen": [20, 20, 10, 12]
        }
    )
    refids_lookup = {"abc": 0, "def": 1}
    expected = np.array([[0, 20], [1, 10]])
    returned = parse_reference_lengths(df, refids_lookup)
    np.testing.assert_array_equal(returned, expected)

#
# def test_parse_blast_df_into_sam_array():
#     # build mini blast_df
#     # build expected_df
#     assert False
