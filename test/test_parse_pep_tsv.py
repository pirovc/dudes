from pyfaidx import Fasta

from dudes.parse_pep_tsv import read_peptide_tsv, get_sequence_for_accession, get_length_of_reference_sequence, \
    get_peptide_start_pos_in_protein_sequence, get_swissprot_accs_matching_uniparc_accs
from test.helper_funcs import RESSOURCE_DIR


def test_read_peptide_tsv():
    df = read_peptide_tsv(RESSOURCE_DIR / "peptide.tsv")


def test_get_sequence_for_accession():
    fasta_path = RESSOURCE_DIR / "uniprot_sprot-head.fasta"
    fasta_obj = Fasta(str(fasta_path), key_function=lambda x: x.split("|")[1])

    assert (
        get_sequence_for_accession("Q6GZX4", fasta_obj)
        == "MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPS"
        "EKGLIVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRRYQGIQTCKIPGKVLSDLD"
        "AKIKAYNLTVEGVEGFVRYSRVTKQHVAAFLKELRHSKQYENVNLIHYILTDKRVDIQHL"
        "EKDLVKDFKALVESAHRMRQGHMINVKYILYQLLKKHGHGPDGPDILTVKTGSKGVLYDD"
        "SFRKIYTDLGWKFTPL"
    )


def test_get_length_of_reference_sequence():
    fasta_path = RESSOURCE_DIR / "uniprot_sprot-head.fasta"
    fasta_obj = Fasta(str(fasta_path), key_function=lambda x: x.split("|")[1])
    assert get_length_of_reference_sequence("Q6GZX4", fasta_obj) == 256


def test_get_peptide_start_pos_in_protein_sequence():
    assert get_peptide_start_pos_in_protein_sequence("asd", "wasd") == 2


def test_get_swissprot_accs_matching_uniparc_accs():
    id_file = RESSOURCE_DIR / "idmapping_uniparc-head.tsv"
    up_accs = ["UPI0000D83464", "UPI0000D83465"]
    assert get_swissprot_accs_matching_uniparc_accs(up_accs, id_file) == {
        "UPI0000D83464": {"Q197F8", "Q197F9"},
        "UPI0000D83465": {"Q197F7"},
    }
