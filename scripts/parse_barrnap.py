
"""Filters the barrnap rRNA gene prediction results table based on a user specified E-value. Furthermore all contig_ids of the input text file are filtered for contigs not containing rRNA based on the barrnap results and the these contig ids are output into a text file."""

import sys
import pandas as pd

barrnap_table = sys.argv[1]
e_value = sys.argv[2]
all_contig_ids = sys.argv[3]
non_rna_contigs_output = sys.argv[4]

def filter_barrnap_table(barrnap_table):
    """Filters the barrnap results for the user-specified E-value and returns a filtered pandas dataframe."""
    barrnap_output_df = pd.read_csv(barrnap_table, sep="\t", header=None)
    barrnap_output_df.columns = [
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ]

    barrnap_output_df["score"] = barrnap_output_df["score"].astype(float)
    filtered_df = barrnap_output_df[barrnap_output_df["score"] <= float(e_value)]
    return filtered_df

filtered_df = filter_barrnap_table(barrnap_table)
contigs_with_rna = set(filtered_df["seqid"].tolist())


with open(all_contig_ids, "r") as file:
    all_contig_ids = set(file.read().splitlines())

non_rna_contigs = all_contig_ids - contigs_with_rna

with open(non_rna_contigs_output, "w") as file:
    for contig in non_rna_contigs:
        file.write(contig + "\n")
