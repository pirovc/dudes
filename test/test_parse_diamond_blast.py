import numpy as np
import pandas as pd
import pytest

from dudes.parse_diamond_blast import read_blast_tsv, parse_uniprot_accession, parse_reference_lengths, \
    transform_blast_df_into_sam_array, parse_custom_blast


@pytest.fixture
def custom_blast_test_file(resource_dir):
    return resource_dir / "diamond_blast_minimal-qseqid-sseqid-slen-sstart-cigar-pident-mismatch.tsv"


@pytest.fixture
def blast_df(custom_blast_test_file):
    return read_blast_tsv(custom_blast_test_file)


def test_read_into_dataframe(custom_blast_test_file):
    returned = read_blast_tsv(custom_blast_test_file)
    assert isinstance(returned, pd.DataFrame)
    assert returned.columns.values.tolist() == ["qseqid", "sseqid", "slen", "sstart", "cigar", "pident", "mismatch"]


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


def test_transform_blast_df_into_sam_array(blast_df):
    refid_lookup = {refid: i for i, refid in enumerate(blast_df["sseqid"].unique())}
    sam_array = transform_blast_df_into_sam_array(blast_df, refid_lookup)
    assert isinstance(sam_array, np.ndarray)
    assert len(np.unique(sam_array[:, 3])) > 1, "Expected multiple Read ID values, got 1"
    assert sam_array.shape[1] == 4, "Missing Score column in sam array"


def test_parse_custom_blast(custom_blast_test_file, blast_df):
    refid_lookup = {refid: i for i, refid in enumerate(blast_df["sseqid"].unique())}
    sam_array, reference_lengths = parse_custom_blast(custom_blast_test_file, refid_lookup, threads=1)
    assert sam_array.shape[1] == 4
    assert reference_lengths.shape[1] == 2
