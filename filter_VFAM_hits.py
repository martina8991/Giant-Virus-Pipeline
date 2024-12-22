""" 
This script is used to find contigs which possibly belong to the Nucleocytoviricota phylum. 
Input is the per-sequence-hits table of the hmmsearch against the VFAMs of the Virus Orthologous Groups Database.
First the hits will be filtered for E-values lower than defined by the user and hits with higher bitscores than bias values. 
If one gene still has two or more significant hits against different VFAMs, the hit with the higher bitscore as well as the smaller E-value and bias is chosen.
The contigs of the filtered hits will be evaluated for the number of significant hits against VFAMs of different taxa. A contig is considered a possible giant virus if the highest number of significant VFAMS is against VFMAs of the Nucleocytoviricota phylum. These contigs are filtered and the count against each taxon of all taxonomical ranks is output as a csv file.
"""

import Bio.SearchIO
import sys
import pandas as pd
import numpy as np


def read_input_table(input_table):
    """Takes the per-sequence-hits table from the hmmsearch against the VFAMS of the VOGDB as input, filters the hits based on a user-defined E-value and bitscore higher than the bias and returns a pandas dataframe with the columns: VFAM_id, gene_id, e-value, bit-score, bias."""
    gene_hit_list_filtered = []
    with open(input_table, "r") as handle:
        for VFAM in Bio.SearchIO.parse(handle, "hmmer3-tab"):
            for gene in VFAM.hits: 
                if (
                    float(gene.evalue) <= float(evalue_cut_off)
                    and gene.bitscore > gene.bias
                ):  
                    gene_hit_list_filtered.append(
                        [VFAM.id, gene.id, gene.evalue, gene.bitscore, gene.bias]
                    )
    gene_hits_filtered_df = pd.DataFrame(
        gene_hit_list_filtered,
        columns=[
            "VFAM_id",
            "gene_id",
            "e-value",
            "bit-score",
            "bias",
        ],
    )
    return gene_hits_filtered_df


def filter_best_hit(gene_hits_filtered_df):
    """Filters the dataframe for the single best VFAM hit for each gene. 
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
        ["contig_id", "gene_id", "VFAM_id", "e-value", "bit-score", "bias"]
    ]
    return best_hit_df


def get_vfam_list(vfam_ids):
    """Reads the VFAMs ids and the respective taxonomical identification from the metadata file "vfam.lca.tsv" and returns a pandas dataframe with the taxonomical rankings split into separate columns."""
    vfam_ids_df = pd.read_csv(vfam_ids, sep="\t")
    vfam_ids_df = vfam_ids_df[["#GroupName", "LastCommonAncestor_Name"]]
    vfam_ids_df.rename(
        columns={"#GroupName": "VFAM_id", "LastCommonAncestor_Name": "Taxa"},
        inplace=True,
    )
    taxon_list = [
        "tax_rank_1",
        "tax_rank_2",
        "tax_rank_3",
        "tax_rank_4",
        "tax_rank_5",
        "tax_rank_6",
        "tax_rank_7",
        "tax_rank_8",
        "tax_rank_9",
        "tax_rank_10",
        "tax_rank_11",
        "tax_rank_12",
    ]
    split_taxa = vfam_ids_df["Taxa"].str.split(";", expand=True)
    split_taxa.fillna(value=np.nan)
    for column in split_taxa.columns:
        if isinstance(column, int):
            split_taxa.rename(columns={column: taxon_list[column]}, inplace=True)
    vfam_ids_df = pd.concat([vfam_ids_df, split_taxa], axis=1)
    vfam_ids_df.drop("Taxa", axis=1, inplace=True)
    vfam_ids_df.fillna(value=np.nan)
    return vfam_ids_df


def merge_df(best_hit_df, vfam_ids_df):
    """Merges the best-hit dataframe with the VFAM list dataframe and returns the merged dataframe."""
    merged_df = pd.merge(best_hit_df, vfam_ids_df, on="VFAM_id", how="left")
    return merged_df


def get_nucleocytoviricota(merged_df):
    """Filter for rows with contig_ids from the best-hit dataframe which could belong to the phylum Nucleocytoviricota. The contig must have the most amount of significant gene hits against VFAMs belonging to Nucleocytoviricota. The filtered dataframe is returned"""
    possible_gv = []
    grouped_by_contig_id = merged_df.groupby("contig_id")
    for contig, group in grouped_by_contig_id:
        taxa_level_trans = {
            "tax_rank_2": "realm",
            "tax_rank_3": "kingdom",
            "tax_rank_4": "phylum",
        }
        required_taxa = [
            "Nucleocytoviricota",
            "Bamfordvirae",
            "Varidnaviria",
            "Pithoviruses",
            "Pandoravirus",
        ]
        counts_realm = group["tax_rank_2"].value_counts().nlargest(1, keep="all")
        list_realm = counts_realm.index.tolist()
        if len(list_realm) == 1 and "Varidnaviria" in list_realm:
            counts_kingdom = group["tax_rank_3"].value_counts().nlargest(1, keep="all")
            list_kingdom = counts_kingdom.index.tolist()
            counts_phylum = group["tax_rank_4"].value_counts().nlargest(1, keep="all")
            list_pyhlum = counts_phylum.index.tolist()
            if (
                len(list_kingdom) == 1
                and len(list_pyhlum) == 1
                and "Bamfordvirae" in list_kingdom
                and "Nucleocytoviricota" in list_pyhlum
            ):
                rows = group.values.tolist()
                for r in rows:
                    possible_gv.append(r)
        elif (len(list_realm) == 1 and "Pithoviruses" in list_realm) or (
            len(list_realm) == 1 and "Pandoravirus" in list_realm
        ):
            rows = group.values.tolist()
            for r in rows:
                possible_gv.append(r)
    possible_gv_df = pd.DataFrame(possible_gv, columns=merged_df.columns)
    return possible_gv_df


def count_taxa_per_contig(possible_gv_df):
    """Returns a dictionary which counts for each contig of the "possible Nucleocytoviricota dataframe" the significant gene hits for all taxa of the different taxonomical ranks. """
    grouped_by_contig_id = possible_gv_df.groupby("contig_id")
    taxon_tax_ranks = [
        "tax_rank_1",
        "tax_rank_2",
        "tax_rank_3",
        "tax_rank_4",
        "tax_rank_5",
        "tax_rank_6",
        "tax_rank_7",
        "tax_rank_8",
        "tax_rank_9",
        "tax_rank_10",
        "tax_rank_11",
        "tax_rank_12",
    ]
    contig_taxa_count = dict()
    for name, group in grouped_by_contig_id:
        contig_taxa_count[name] = []
        for column in group.columns:
            if column in taxon_tax_ranks:
                counts = group[column].value_counts()
                column_dict = dict()
                for idx, value in enumerate(counts.index.tolist()):
                    column_dict[value] = counts.iloc[idx]
                contig_taxa_count[name].append(column_dict)
    return contig_taxa_count


def write_output_table(contig_taxa_count, output_table):
    """Write the contig-VFAM-count dictionary as a CSV file into output."""
    with open(output_table, "w") as f:
        f.write(
            f"contig_id\ttax_rank_1\ttax_rank_2\ttax_rank_3\ttax_rank_4\ttax_rank_5\ttax_rank_6\ttax_rank_7\ttax_rank_8\ttax_rank_9\ttax_rank_10\ttax_rank_11\ttax_rank_12\n"
        )
        for contig, taxa in contig_taxa_count.items():
            names = ""
            for taxon in taxa:
                if taxon == {}:
                    continue
                taxon_string = ""
                for key, value in taxon.items():
                    taxon_string = taxon_string + f"{key}:{value},"
                taxon_list = taxon_string.split(",")
                taxon_list.pop()
                taxon_string = ",".join(taxon_list)
                names = names + taxon_string + "\t"
            f.write(f"{contig}\t{names}\n")

if __name__ == "__main__":
    input_table = sys.argv[1]
    vfam_ids = sys.argv[2]
    output_table = sys.argv[3]
    evalue_cut_off = sys.argv[4]
    gene_hits_filtered_df = read_input_table(input_table)
    best_hit_df = filter_best_hit(gene_hits_filtered_df)
    vfam_ids_df = get_vfam_list(vfam_ids)
    merged_df = merge_df(best_hit_df, vfam_ids_df)
    possible_gv_df = get_nucleocytoviricota(merged_df)
    contig_taxa_count = count_taxa_per_contig(possible_gv_df)
    write_output_table(contig_taxa_count, output_table)
