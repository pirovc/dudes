from dudes.main import main
from test.helper_funcs import RESOURCE_DIR, DUDES_DIR, SAMPLEDATA_DIR, md5sum, set_directory
from unittest.mock import patch
import sys


def test_DUDes_sam(tmp_path):
    produced_output = tmp_path / "dudes_profile_output"
    expected_output = RESOURCE_DIR / "dudes_profile_expected_output_on_sampledata.out"
    args = ["dudes",
            "-s", str(SAMPLEDATA_DIR.relative_to(DUDES_DIR) / "hiseq_accuracy_k60.sam"),
            "-d", str(SAMPLEDATA_DIR.relative_to(DUDES_DIR) / "arc-bac_refseq-cg_201503.npz"),
            "-o", str(produced_output)
            ]
    with set_directory(DUDES_DIR):
        with patch.object(sys, "argv", args):
            main()
    assert md5sum(produced_output.with_suffix(".out")) == md5sum(expected_output), \
        f'wrong dudes output produced: {produced_output.with_suffix(".out")}\n'\
        f'expected output: {expected_output}'


def test_DUDes_blast(tmp_path):
    produced_output = tmp_path / "dudes_profile_output"
    args = ["dudes",
            "-c", str(RESOURCE_DIR.relative_to(DUDES_DIR) / "diamond_blast_minimal-qseqid-sseqid-slen-sstart-evalue.tsv"),
            "-d", str(SAMPLEDATA_DIR.relative_to(DUDES_DIR) / "arc-bac_refseq-cg_201503.npz"),
            "-o", str(produced_output)
            ]
    with set_directory(DUDES_DIR):
        with patch.object(sys, "argv", args):
            main()
