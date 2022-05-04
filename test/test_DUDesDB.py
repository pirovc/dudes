import sys
from unittest.mock import patch

import numpy as np

import dudes
from dudes.database.main import load_refid2taxid_files, get_reference_identifiers, main
from test.helper_funcs import RESOURCE_DIR, DUDES_DIR, md5sum, set_directory


def test_load_refid2taxid_files_up():
    refids = {'Q6GZX4', 'Q6GZX3', 'Q197F8', 'Q197F7', 'Q6GZX2', 'Q6GZX1', 'Q197F5', 'Q6GZX0', 'Q91G88', 'Q6GZW9', 'Q6GZW8', 'Q197F3', 'Q197F2', 'Q6GZW6', 'Q91G85', 'Q6GZW5'}
    refid2taxid_files = [str(RESOURCE_DIR / "idmapping_selected-test.tab")]
    refid_taxid, refids_lookup = load_refid2taxid_files(refid2taxid_files, "up", refids)
    expected_refid_taxid = np.array([[0, 654924], [1, 654924], [2, 345201], [3, 345201], [4, 654924], [5, 654924],
                                     [6, 345201], [7, 654924], [8, 176652], [9, 654924], [10, 654924], [11, 345201],
                                     [12, 345201], [13, 654924], [14, 176652], [15, 654924]])
    expected_refids_lookup = {
        'Q6GZX4': 0, 'Q6GZX3': 1, 'Q197F8': 2, 'Q197F7': 3, 'Q6GZX2': 4, 'Q6GZX1': 5, 'Q197F5': 6, 'Q6GZX0': 7,
        'Q91G88': 8, 'Q6GZW9': 9, 'Q6GZW8': 10, 'Q197F3': 11, 'Q197F2': 12, 'Q6GZW6': 13, 'Q91G85': 14, 'Q6GZW5': 15}
    assert np.array_equal(expected_refid_taxid, refid_taxid)
    assert expected_refids_lookup == refids_lookup


def test_get_reference_identifiers_uniprot():
    expected_accs = {'Q6GZX4', 'Q6GZX3', 'Q197F8', 'Q197F7', 'Q6GZX2', 'Q6GZX1', 'Q197F5', 'Q6GZX0', 'Q91G88', 'Q6GZW9', 'Q6GZW8', 'Q197F3', 'Q197F2', 'Q6GZW6', 'Q91G85', 'Q6GZW5'}
    assert expected_accs == get_reference_identifiers([str(RESOURCE_DIR / "uniprot_sprot-head.fasta")], "up")


def test_DUDesDB(tmp_path):
    produced_db = tmp_path / "dudesdb.npz"
    expected_db = RESOURCE_DIR / "dudesdb.npz"
    args = [DUDES_DIR / "DUDesDB.py",
            "-m", "up",
            "-f", RESOURCE_DIR / "uniprot_sprot-head.fasta",
            "-g", RESOURCE_DIR / "idmapping_selected-test.tab",
            "-n", RESOURCE_DIR / "nodes.dmp",
            "-a", RESOURCE_DIR / "names.dmp",
            "-o", tmp_path / "dudesdb"]
    args_str = list(map(str, args))
    with set_directory(DUDES_DIR):
        with patch.object(sys, "argv", args_str):
            with patch.object(dudes.database.main, "VERSION", "0.09"):
                main()
    assert md5sum(produced_db) == md5sum(expected_db), \
        'wrong dudes database file produced'
