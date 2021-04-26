import subprocess
from test.helper_funcs import RESSOURCE_DIR, DUDES_DIR, SAMPLEDATA_DIR, md5sum


def test_DUDes(tmp_path):
    produced_output = tmp_path / "dudes_profile_output"
    expected_output = RESSOURCE_DIR / "dudes_profile_expected_output_on_sampledata.out"
    args = [DUDES_DIR / "DUDes.py",
            "-s", SAMPLEDATA_DIR.relative_to(DUDES_DIR) / "hiseq_accuracy_k60.sam",
            "-d", SAMPLEDATA_DIR.relative_to(DUDES_DIR) / "arc-bac_refseq-cg_201503.npz",
            "-o", produced_output
            ]
    subprocess.call(args, cwd=DUDES_DIR)
    assert md5sum(produced_output.with_suffix(".out")) == md5sum(expected_output), \
        f'wrong dudes output produced: {produced_output.with_suffix(".out")}\n'\
        f'expected output: {expected_output}'
