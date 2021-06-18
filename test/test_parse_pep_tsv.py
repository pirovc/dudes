import numpy as np
import pandas as pd
import pytest

from dudes.parse_pep_tsv import (
    FastaExtension,
    Peptide2ReferenceTable,
    get_peptide_start_pos_in_protein_sequence,
    get_uniparc_to_ref_id_map,
    build_dfs,
    get_uniparc_to_uniprot_acc_map
)
from test.helper_funcs import RESSOURCE_DIR


def test_read_peptide_tsv():
    pep2ref = Peptide2ReferenceTable(RESSOURCE_DIR / "peptide.tsv")
    # pep2ref.df.to_csv(RESSOURCE_DIR / "pep2acc.csv", index=False)
    expected_df = pd.read_csv(
        RESSOURCE_DIR / "pep2acc.csv", converters={"Proteins": eval}
    )
    pd.testing.assert_frame_equal(pep2ref.df, expected_df)


def test_get_set_of_accs_from_pep_df():
    from .ressource.parse_pep_data import (
        expected_acc_set_get_set_of_accs_from_pep_df as expected_set,
    )
    pep2ref = Peptide2ReferenceTable(
        pep2ref_df=pd.read_csv(
            RESSOURCE_DIR / "pep2acc.csv", converters={"Proteins": eval}
        )
    )
    produced_set = pep2ref.get_all_accs()
    assert expected_set == produced_set, (
        f"Expected but missing accs: {expected_set - produced_set}\n"
        f"Unexpected accs: {produced_set - expected_set}"
    )


def test_get_peptides_matching_acc():
    df = pd.read_csv(RESSOURCE_DIR / "pep2acc.csv", converters={"Proteins": eval})
    expected_peps = {"AAAAGFEK", "AAEYMTHAPLGSLNSVGGVATEINAVNFVSPR"}
    pep2ref = Peptide2ReferenceTable(pep2ref_df=df)
    res = pep2ref.get_peptides_matching_acc("UPI0000132624")
    assert res == expected_peps


def test_get_pep_id():
    df = pd.read_csv(RESSOURCE_DIR / "pep2acc.csv", converters={"Proteins": eval})
    pep2ref = Peptide2ReferenceTable(pep2ref_df=df)
    assert pep2ref.get_pep_id("AAAEDAGLPLYR") == 9


def test_get_sequence_for_accession():
    fasta_path = RESSOURCE_DIR / "uniprot_sprot-head.fasta"
    fasta_obj = FastaExtension(fasta_path, key_function=lambda x: x.split("|")[1])

    assert (
        fasta_obj.get_sequence_for_accession("Q6GZX4")
        == "MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPS"
        "EKGLIVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRRYQGIQTCKIPGKVLSDLD"
        "AKIKAYNLTVEGVEGFVRYSRVTKQHVAAFLKELRHSKQYENVNLIHYILTDKRVDIQHL"
        "EKDLVKDFKALVESAHRMRQGHMINVKYILYQLLKKHGHGPDGPDILTVKTGSKGVLYDD"
        "SFRKIYTDLGWKFTPL"
    )


def test_get_length_of_reference_sequence():
    fasta_path = RESSOURCE_DIR / "uniprot_sprot-head.fasta"
    fasta_obj = FastaExtension(fasta_path, key_function=lambda x: x.split("|")[1])
    assert fasta_obj.get_length_of_reference_sequence("Q6GZX4") == 256


def test_get_peptide_start_pos_in_protein_sequence():
    assert get_peptide_start_pos_in_protein_sequence("asd", "wasd") == 2


def test_get_swissprot_accs_matching_uniparc_accs():
    id_file = RESSOURCE_DIR / "idmapping_uniparc-test.tsv"
    up_accs = ["UPI000012529B", "UPI00017B9579"]
    expected = {
        "UPI000012529B": {
            "A0A0D1KAE8",
            "A0A4U9YN62",
            "A6QHP5",
            "A8Z2M7",
            "Q2FG27",
            "Q2FXL5",
            "Q2YTC9",
            "Q5HF63",
            "Q6G8L6",
            "Q8NW53",
            "X5E067",
        },
        "UPI00017B9579": {
            "A0A095G6Y0",
            "A0A104KHE9",
            "A0A142PQF5",
            "A0A2G1SKG6",
            "A0A364GMC2",
            "B4ECN0",
        },
    }
    result = get_uniparc_to_uniprot_acc_map(up_accs, id_file)
    assert expected == result


@pytest.mark.xfail
def test_get_uniparc_to_ref_id_map():
    from .ressource.parse_pep_data import (
        expected_acc_set_get_set_of_accs_from_pep_df as up_accs,
    )

    # print("\n", "|".join(up_accs))
    npzfile = np.load(RESSOURCE_DIR / "dudesdb.npz", allow_pickle=True)
    refids_lookup = npzfile["refids_lookup"].item()
    id_mapping_file = RESSOURCE_DIR / "idmapping_uniparc-test.tsv"
    produced_dict = get_uniparc_to_ref_id_map(up_accs, id_mapping_file, refids_lookup)
    print(produced_dict)
    assert False


def test_build_dfs():
    fasta = RESSOURCE_DIR / "uniparc_filtered.fasta"
    refids_lookup = ""
    idmapping_file = RESSOURCE_DIR / "idmapping_uniparc-test.tsv"
    pep2ref = Peptide2ReferenceTable(
        pep2ref_df=pd.read_csv(
            RESSOURCE_DIR / "pep2acc.csv", converters={"Proteins": eval}
        )
    )
    fasta_extension_obj = FastaExtension(
        str(fasta), sequence_always_upper=True
    )
    # get match match start pos, pep-length (=Score) and reference length
    ref_len, df = build_dfs(pep2ref, fasta_extension_obj)
    uprc2uprt = get_uniparc_to_uniprot_acc_map(pep2ref.get_all_accs(), idmapping_file)
    print(len(uprc2uprt))
    df["SwissProtAccs"] = df["RefID"].map(lambda x: uprc2uprt.get(x, set()))

    print("\n", df)
    # replace uniparc acc by swissprot acc
    # replace swissprot acc by reference ID, missing accs by -1
    # return two tables
    # 1: 'RefID','MatchPosStart','MatchScore','ReadID'
    # 2: reference id (str), reference length (int)
