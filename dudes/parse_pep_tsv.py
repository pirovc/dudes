from collections import defaultdict

import pandas as pd
from numpy import nan


def read_peptide_tsv(peptide_tsv_file):
    """

    :param peptide_tsv_file: peptide tsv file as produced by msfragger
    :return: pandas data frame with two columns:
        Peptide: str, Peptide sequence
        Proteins: str, ', '-separated accessions that contain the peptide sequence
    """
    df = pd.read_csv(
        peptide_tsv_file, sep="\t", usecols=["Peptide", "Protein", "Mapped Proteins"]
    )
    df["Proteins"] = df[["Protein", "Mapped Proteins"]].apply(
        lambda x: x["Protein"] if x["Mapped Proteins"] is nan else ", ".join(x), axis=1
    )
    df.drop(columns=["Protein", "Mapped Proteins"], inplace=True)
    return df


def get_sequence_for_accession(acc, fasta_obj):
    return fasta_obj[acc][:].seq


def get_length_of_reference_sequence(acc, fasta_obj):
    return len(fasta_obj[acc])


def get_peptide_start_pos_in_protein_sequence(peptide, protein):
    return protein.index(peptide)


def get_swissprot_accs_matching_uniparc_accs(uniparc_accs, idmapping_file):
    up_accs = set(uniparc_accs)
    up_to_sp = defaultdict(set)
    with open(idmapping_file, "rt") as f:
        for line in f:
            sp_acc, up_acc = line.strip("\n").split("\t")
            if up_acc in up_accs:
                up_to_sp[up_acc].add(sp_acc)
    return up_to_sp


def main():
    # fasta = ?
    # peptide_tsv = ?
    # refids_lookup = ?
    peptides = defaultdict(lambda: len(peptides))
    # return two tables
    # 1: 'RefID','MatchPosStart','MatchScore','ReadID'
    # 2: reference id (str), reference length (int)

    pass
