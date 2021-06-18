#!/usr/bin/env python3
import argparse
import sys
from tqdm import tqdm

import pandas as pd

from dudes.parse_pep_tsv import Peptide2ReferenceTable


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p",
        required=True,
        metavar="<pep_file>",
        dest="pep_file",
        help="MSFragger peptide.tsv result file",
    )
    parser.add_argument(
        "-i",
        required=True,
        metavar="<idmapping_file>",
        dest="idmapping_file",
        help="reference id to taxid file: "
        "'idmapping_selected.tab[.gz]' --> 'up' mode "
        "[from https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz",
    )
    parser.add_argument(
        "-o",
        metavar="<output>",
        dest="output",
        default=sys.stdout,
        help="Output file. Default: STDOUT",
    )
    return parser.parse_args()


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


def replace_uprc_w_uprt(acc_string, uprc_to_uprt_map):
    if pd.isna(acc_string):
        return acc_string
    uprt_accs = []
    for acc in acc_string.split(", "):
        if acc in uprc_to_uprt_map:
            uprt_accs.extend(uprc_to_uprt_map[acc])
    return ", ".join(set(uprt_accs))


def main():
    args = parse_args()
    df_pep_raw = pd.read_csv(args.pep_file, sep="\t")
    df_pep_filt = Peptide2ReferenceTable(pep2ref_df=df_pep_raw)
    uniparc_accs = df_pep_filt.get_all_accs()
    uprc_to_uprt_map = get_uniparc_to_uniprot_acc_map(uniparc_accs, args.idmapping_file)
    for col in ["Protein", 'Mapped Proteins']:
        df_pep_raw[col] = df_pep_raw[col].map(lambda x: replace_uprc_w_uprt(x, uprc_to_uprt_map))
    df_pep_raw.to_csv(args.output, index=False, sep="\t")


if __name__ == "__main__":
    main()
