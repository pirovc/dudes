import argparse
from pyfaidx import Fasta
from tqdm import tqdm
import pandas as pd
from numpy import nan
import numpy as np
import sys
import time


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
        :return: pandas data frame with Peptide sequence as index and three columns:
            MatchScore: int, converted probability p: int(p.replace(".", ""))
            Proteins: list, accessions of proteins that contain the peptide sequence
            Spectral Count: int, spectral count of the peptide
        """
        df = pd.read_csv(
            msfragger_tsv_file,
            sep="\t",
            usecols=["Peptide", "Probability", "Spectral Count", "Protein", "Mapped Proteins"],
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

    def get_pep_count(self, pep):
        """
        get spectral count of peptide sequence

        :param pep: str
        :return: int
        """
        return self.df.loc[pep, "Spectral Count"]


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
    :return: 2 np.arrays:
        1) read_table: columns: 'RefID','MatchPosStart','MatchScore','ReadID'
        2) reference lengths: columns: 'RefID', 'ReferenceSequenceLength'
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
                        pep_seq,
                        pep2ref.get_pep_count(pep_seq)
                    ]
                )
    df_reads = pd.DataFrame(
        lst_read_table_rows, columns=["RefID", "MatchPosStart", "MatchScore", "PeptideSeq", "ReadCount"]
    )

    grouped = df_reads.groupby("PeptideSeq", sort=False)
    pep_id = 0
    dfs_with_pep_id = []
    for _, df in grouped:
        for i in range(df["ReadCount"].iloc[0]):
            df_c = df.copy()
            df_c["ReadID"] = pep_id
            dfs_with_pep_id.append(df_c)
            pep_id += 1
    df_reads = pd.concat(dfs_with_pep_id, axis=0, ignore_index=True)

    df_reads = df_reads[["RefID", "MatchPosStart", "MatchScore", "ReadID"]]

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


def parse_args():
    parser = argparse.ArgumentParser(
        description="Parse MsFragger peptide identification results into dudes-readable format. Input accessessions "
                    "are expected to be uniparc, mapped to uniprot (swissprot and trembl)"
    )
    parser.add_argument(
        "-m",
        required=True,
        metavar="<msfragger_file>",
        dest="msfragger_file",
        help="MsFragger peptide.tsv result file",
    )
    parser.add_argument(
        "-d",
        required=True,
        metavar="<database_file>",
        dest="database_file",
        help="Database file (output from DUDesDB [.npz])",
    )
    parser.add_argument(
        "-f",
        required=True,
        metavar="<fasta_file>",
        dest="fasta_file",
        help="Uniparc fasta file(s) for mapping the position of peptide sequences and retrieving reference sequence"
             "lengths. Each sequence header '>' should contain an "
        "identifier in the same format as in the MsFragger input file.",
    )
    parser.add_argument(
        "-i",
        required=True,
        metavar="<idmapping_file>",
        dest="idmapping_file",
        help="reference idmapping file, used to map uniparc accessions to uniprot accessions: "
        "'idmapping_selected.tab[.gz]'"
        "[from https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz",
    )
    parser.add_argument(
        "-o",
        metavar="<output_prefix>",
        dest="output_prefix",
        default="",
        help="Output prefix. Default: STDOUT",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    fasta = args.fasta_file
    idmapping_file = args.idmapping_file
    output_prefix = args.output_prefix

    refids_lookup = np.load(args.database_file, allow_pickle=True)["refids_lookup"].item()
    pep2ref = Peptide2ReferenceTable(msfragger_tsv_file=args.msfragger_file)
    fasta_extension_obj = FastaExtension(
        fasta, sequence_always_upper=True
    )
    reads, ref_len = build_dfs(pep2ref, fasta_extension_obj, refids_lookup, idmapping_file)

    # Save processed peptide identifications
    sys.stdout.write("Saving database %s ..." % (output_prefix + ".npz"))
    tx = time.time()
    np.savez(output_prefix, reads=reads, reference_lengths=ref_len)
    sys.stdout.write(" Done. Elapsed time: " + str(time.time() - tx) + " seconds\n")


if __name__ == "__main__":
    main()
