""" 
This script is used to filter the best GVOG hit per gene from the per-sequence-hits table of the hmmsearch against the Giant Virus Database. 
The hits will be first filtered for E-values lower than a user-defined value and bitscores need to be higher than the bias value. 
If one gene still has two significant hits against different GVOGs, the hit with the higher bitscore as well as the smaller E-value and bias is chosen.
The filtered hits will be written to a three new csv files (core3 hits, extcore7 hits, other hits) with the columns: contig_id, gene_id, GVOG_id, e-value, bit-score, bias.
"""

import Bio.SearchIO
import sys
import pandas as pd


def read_input_table(input_table):
    """Reads the input per-sequence-hits table from hmmsearch against the GVDB and filters for a user-defined E-value as well as hits with bitscores higher than the bias. 
    The function returns a pandas dataframe with the columns: GVOG_id, gene_id, e-value, bit-score, bias."""
    gene_hit_list_filtered = []
    with open(input_table, "r") as handle:
        for GVOG in Bio.SearchIO.parse(handle, "hmmer3-tab"):
            for gene in GVOG.hits:  
                if (
                    float(gene.evalue) <= float(evalue_cut_off)
                    and gene.bitscore > gene.bias
                ):  
                    gene_hit_list_filtered.append(
                        [
                            GVOG.id,
                            gene.id,
                            gene.evalue,
                            gene.bitscore,
                            gene.bias,
                        ]
                    )

        gene_hits_filtered_df = pd.DataFrame(
            gene_hit_list_filtered,
            columns=["GVOG_id", "gene_id", "e-value", "bit-score", "bias"],
        )
    return gene_hits_filtered_df



def filter_best_hit(gene_hits_filtered_df):
    """Filters the dataframe for the single best GVOG hit for each gene. 
    If a gene has multiple hits with the same bitscore, the hit with the lower bias and e-value will be chosen."""
    gene_hits_filtered_sorted_df = gene_hits_filtered_df.sort_values(
        by=["gene_id", "bit-score", "bias", "e-value"],
        ascending=[True, False, True, True],
    )
    best_hit_df = gene_hits_filtered_sorted_df.groupby(["gene_id"], sort=False).first().reset_index()
    for index, row in best_hit_df.iterrows():
        contig_id = row["gene_id"].split("_")[0:2]
        contig_id = "_".join(contig_id)
        best_hit_df.at[index, "contig_id"] = contig_id
    best_hit_df = best_hit_df[
        ["contig_id", "gene_id", "GVOG_id", "e-value", "bit-score", "bias"]
    ]
    return best_hit_df


core3_list = ["GVOGm0095", "GVOGm0890", "GVOGm0054"]
ext_core7_list = [
    "GVOGm0003",
    "GVOGm0013",
    "GVOGm0022",
    "GVOGm0023",
    "GVOGm0172",
    "GVOGm0461",
    "GVOGm0760",
]


def core_hits(best_hit_df, core_file, core_list):
    """Filter the single-best hits dataframe for either core3 or extcore7 gene hits and write the output to a csv file with the columns: contig_id, gene_id, GVOG_id, e-value, bit-score, bias"""
    core_df = best_hit_df[best_hit_df["GVOG_id"].isin(core_list)]
    core_df.to_csv(core_file, sep="\t", index=False)


def other_hits(best_hit_df, other_hits_file):
    """Filters the remaining hits from the single-best hits dataframe and writes the results to a new file with the columns: contig_id, gene_id, GVOG_id, e-value, bit-score, bias"""
    other_hits_df = best_hit_df[~best_hit_df["GVOG_id"].isin(core3_list)]
    other_hits_df = other_hits_df[~other_hits_df["GVOG_id"].isin(ext_core7_list)]
    other_hits_df.to_csv(other_hits_file, sep="\t", index=False)


if __name__ == "__main__":
    input_table = sys.argv[1]

    core3_file = sys.argv[2]
    ext_core7_file = sys.argv[3]
    other_hits_file = sys.argv[4]

    evalue_cut_off = sys.argv[5]

    filtered_df = read_input_table(input_table)

    best_hit_df = filter_best_hit(filtered_df)

    core3_df = core_hits(best_hit_df, core3_file, core3_list)
    ext_core7_df = core_hits(best_hit_df, ext_core7_file, ext_core7_list)
    other_hits_df = other_hits(best_hit_df, other_hits_file)
