from collections import defaultdict
from pyfaidx import Fasta
from tqdm import tqdm
import pandas as pd
from numpy import nan
import numpy as np


class Peptide2ReferenceTable:
    def __init__(self, msfragger_tsv_file=None, pep2ref_df=None):
        self.input_file = msfragger_tsv_file
        self.df = (
            pep2ref_df
            if (pep2ref_df is not None)
            else self.read_peptide_tsv(self.input_file)
        )

    def read_peptide_tsv(self, msfragger_tsv_file):
        """
        Read msfragger peptide.tsv file into pandas dataframe. Dataframe has index Peptide and columns MatchScore, and
        Proteins.
        Proteins column is a list of accessions from peptide.tsv columns "Protein" and "Mapped Proteins".
        Proteins Accession starting with "rev_" are removed.

        :param msfragger_tsv_file: peptide tsv file as produced by msfragger
        :return: pandas data frame with two columns and Peptide sequence as index:
            MatchScore: int, converted probability p: int(p.replace(".", ""))
            Proteins: list, accessions of proteins that contain the peptide sequence
        """
        df = pd.read_csv(
            msfragger_tsv_file,
            sep="\t",
            usecols=["Peptide", "Probability", "Protein", "Mapped Proteins"],
            converters={"Probability": lambda x: int(x.replace(".", ""))},
            index_col="Peptide",
        )
        df["Proteins"] = df[["Protein", "Mapped Proteins"]].apply(
            lambda x: x["Protein"] if x["Mapped Proteins"] is nan else ", ".join(x),
            axis=1,
        )
        df["Proteins"] = df["Proteins"].apply(
            lambda x: [acc for acc in x.split(", ") if not acc.startswith("rev_")]
        )
        df.drop(columns=["Protein", "Mapped Proteins"], inplace=True)
        df.rename({"Probability": "MatchScore"}, axis=1, inplace=True)
        return df

    def get_all_accs(self):
        return {acc for acc_lst in self.df["Proteins"] for acc in acc_lst}

    def get_peptides_matching_acc(self, acc):
        return set(self.df.index[self.df["Proteins"].map(lambda x: acc in x)])

    def get_pep_id(self, pep):
        """
        get index number of peptide (0-based)

        :param pep:
        :return: int
        """
        return self.df.index.get_loc(pep)

    def get_pep_score(self, pep):
        """
        get score of peptide sequence

        :param pep: str
        :return: int
        """
        return self.df.loc[pep, "MatchScore"]


class FastaExtension:
    def __init__(self, fasta_path, *args, **kwargs):
        self.fasta = Fasta(str(fasta_path), *args, **kwargs)

    def get_sequence_for_accession(self, acc):
        return self.fasta[acc][:].seq

    def get_length_of_reference_sequence(self, acc):
        return len(self.fasta[acc])


def get_peptide_start_pos_in_protein_sequence(
    peptide_seq, protein_seq, equate_i_and_l=True
):
    if equate_i_and_l:
        peptide_seq = peptide_seq.replace("I", "L")
        protein_seq = protein_seq.replace("I", "L")
    return protein_seq.index(peptide_seq) + 1


def get_uniparc_to_uniprot_acc_map(uniparc_accs, idmapping_file):
    lst_dfs = []
    chunksize = 10 ** 6
    with pd.read_csv(
        idmapping_file,
        sep="\t",
        header=0,
        usecols=[0, 10],
        names=["uniprot", "uniparc"],
        chunksize=chunksize,
        iterator=True,
    ) as reader:
        for chunk in tqdm(reader, total=int(214971037 / chunksize) + 1):
            lst_dfs.append(chunk[chunk["uniparc"].map(lambda x: x in uniparc_accs)])
    df = pd.concat(lst_dfs, ignore_index=True)
    df = df.groupby("uniparc").agg(set)
    uniparc_to_uniprot_map = df.squeeze().to_dict()
    return uniparc_to_uniprot_map


def get_uniparc_to_ref_id_map(uniparc_accs, idmapping_file, refid_lookup):
    """
    get a map of uniparc accessions to dudes internal reference ids

    :param uniparc_accs: iterable of strings, each string is a uniparc accession
    :param idmapping_file: path to uniprot idmapping-selected.tsv
    :param refid_lookup: dict, key: reference accession, value: dudes internal id of the reference
    :return: dictionary: key: uniparc accession, values: set of reference ids matching the uniparc accession
    """
    uprc_acc2sp_accs = get_uniparc_to_uniprot_acc_map(uniparc_accs, idmapping_file)
    uprc_acc2refids = {}
    for uprc_acc, sp_accs in uprc_acc2sp_accs.items():
        uprc_acc2refids[uprc_acc] = [
            refid_lookup[a] for a in sp_accs if a in refid_lookup
        ]
    return uprc_acc2refids


def build_dfs(pep2ref, fasta_extension_obj, refids_lookup, idmapping_file):
    """

    :param pep2ref: Peptide2ReferenceTable object
    :param fasta_extension_obj: a FastaExtension object
    :param refids_lookup: refid lookup table from dudesdb
    :return: read_table: rows: 'RefID','MatchPosStart','MatchScore','ReadID'
    """
    ref_lengths = []
    lst_read_table_rows = []
    set_of_matched_prot_accs = pep2ref.get_all_accs()
    uprc2uprt = get_uniparc_to_uniprot_acc_map(pep2ref.get_all_accs(), idmapping_file)
    # get pep position in each reference sequence (add ref sequence length to ref_lengths)
    for acc in set_of_matched_prot_accs:
        if acc in fasta_extension_obj.fasta:
            prot_seq = fasta_extension_obj.get_sequence_for_accession(acc)
            ref_lengths.append([acc, len(prot_seq)])
            for pep_seq in pep2ref.get_peptides_matching_acc(acc):
                match_pos_start = get_peptide_start_pos_in_protein_sequence(
                    pep_seq, prot_seq
                )
                lst_read_table_rows.append(
                    [
                        acc,
                        match_pos_start,
                        pep2ref.get_pep_score(pep_seq),
                        pep2ref.get_pep_id(pep_seq),
                    ]
                )
    df_reads = pd.DataFrame(
        lst_read_table_rows, columns=["RefID", "MatchPosStart", "MatchScore", "ReadID"]
    )
    # replace accs by target db accs in dataframe and ref_lengths
    # replace uniparc_accs by uniprot_accs in dataframe
    df_reads["RefID"] = df_reads["RefID"].map(lambda x: uprc2uprt.get(x, set()))
    # explode set of uniprot_accs, one acc per row
    df_reads = df_reads.explode("RefID", ignore_index=True)
    # replace uniprot_accs by refids
    df_reads["RefID"] = df_reads["RefID"].map(lambda acc: refids_lookup.get(acc, -1))

    # replace uniparc_accs by uniprot_refids in ref_lengths
    ref_lengths = [
        [refids_lookup.get(uprt_acc, -1), ref_len]
        for uprc_acc, ref_len in ref_lengths
        for uprt_acc in uprc2uprt.get(uprc_acc, set())
    ]
    return np.array(df_reads), np.array(ref_lengths)
