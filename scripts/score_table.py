import sys
import numpy as np
import pandas as pd
import math

"""
This script calculates a likelihood score for possible giant virus contigs and outputs the results in a tab-delimited table. 
The script takes in 8 arguments:
1. A text file with the number of predicted proteins per contig.
2. A text file with the contig length and GC-content per contig.
3. A text file with the average read coverage per contig.
4. A text file with the GVOG_id of core3 gene hits per contig.
5. A text file with the GVOG_id of extcore7 gene hits per contig.
6. A text file with the GVOG_id of GVDB hits per contig, except core3 and extcore7 hits.
7. A text file listing contigs considered as possible Nuclecocytoviricota based on the search against the Virus Orthologous Groups Database.
8. A text file listing contigs without rRNA detected by Barrnap.
9. The minimum contig length for the score table results defined by the user.
Based on the number of significant hits against core3, extcore7, other GVOG hits and a "hit" against the VOGDB, a score is calculated for each contig.
The output file will be tab-delimited csv file with the following columns: contig_name, contig_len, GC%, average_coverage, all_genes, perc_GVOG_hits/all_genes, score, core3_count,ext_core7_count, other_GVOG_count, VOGDB_hit.
"""


def read_count_genes_per_contig(count_genes_per_contig_file):
    """Read the number of predicted genes per contig from a text file and store it in a dictionary"""
    count_genes_per_contig = dict()
    with open(count_genes_per_contig_file, "r") as file:
        for line in file:
            line = line.strip()
            count, gene = line.split()
            count_genes_per_contig[gene] = count
    return count_genes_per_contig


def read_hit_contig_file(hit_per_contig_file):
    """Reads the GVOG_id of GVDB hits per contig from a text file and stores it in a dictionary"""
    hits_per_contig = dict()
    with open(hit_per_contig_file, "r") as file:
        next(file)
        for line in file:
            line = line.strip()
            contig, gene, GVOG, rest = line.split(maxsplit=3)
            if contig in hits_per_contig:
                hits_per_contig[contig] += [GVOG]
            else:
                hits_per_contig[contig] = [GVOG]
    return hits_per_contig


def get_seq_length_and_gc(seq_length_gc_file):
    """Stores the sequence length and gc content of each contig as a dictionary from the input text file."""
    seq_length_gc = dict()
    with open(seq_length_gc_file, "r") as file:
        for line in file:
            line = line.strip()
            name, length, gc = line.split()
            name = name.split("_")[0:2]
            name = "_".join(name)
            seq_length_gc[name] = [length, gc]
    return seq_length_gc


def get_coverage(coverage_file):
    """Gets the coverage of each contig from an input text file and stores it as a dictionary."""
    coverage = dict()
    with open(coverage_file, "r") as file:
        next(file)
        for line in file:
            line = line.strip()
            name, av_cov, rest = line.split(maxsplit=2)
            name = name.split("_")[0:2]
            name = "_".join(name)
            av_cov = round(float(av_cov), 2)
            coverage[name] = av_cov
    return coverage


def count_hits_per_contig(hits_per_contig):
    """Counts the number of GVOG hits per contig from the input dictionary and stores the counts again as a dictionary."""
    counter_hits = dict()
    for contig in hits_per_contig:
        values = hits_per_contig.get(contig)
        count = len(set(values))
        counter_hits[contig] = count
    return counter_hits


def get_vogdb_hits(vogdb_file):
    """Reads the contigs considered a possible Nucleocytoviricota from the search against VOGDB and stores them as a set."""
    vogdb_contigs = set()
    with open(vogdb_file, "r") as file:
        next(file)
        for line in file:
            contigs, rest = line.strip().split("\t", maxsplit=1)
            vogdb_contigs.add(contigs)
    return vogdb_contigs


def read_contig_ids_without_rRNA(contig_ids_without_rRNA_file):
    """Reads the contig ids of contigs without rRNA from a text file and stores them in a list."""
    contig_ids_without_rRNA = []
    with open(contig_ids_without_rRNA_file, "r") as file:
        for contig_id in file:
            contig_id = contig_id.strip()
            contig_id = contig_id.split("_")[0:2]
            contig_id = "_".join(contig_id)
            contig_ids_without_rRNA.append(contig_id)
    return contig_ids_without_rRNA


def create_matrix(
    all_contigs,
    counter_other_GVOG_hits,
    counter_core3_hits,
    counter_ext_core7_hits,
    vodb_contigs,
    count_genes_per_contig,
):
    """Creates a matrix for the contigs and their scores based on whether the contig has hits against the VOGDB, core3, extcore7 or other GVOGs."""
    dict_matrix = dict()
    for contig in all_contigs:
        dict_matrix[contig] = np.zeros(4)

        if contig in counter_core3_hits:
            value_hits = int(counter_core3_hits.get(contig))
            dict_matrix[contig][0] = value_hits * 3  # max value 9

        if contig in counter_ext_core7_hits:
            value_hits = int(counter_ext_core7_hits.get(contig))
            dict_matrix[contig][1] = value_hits * 2  # max value 14

        if contig in counter_other_GVOG_hits:
            count_hits = counter_other_GVOG_hits.get(contig)
            proportion_of_hits = round(
                float(count_hits) / float(count_genes_per_contig.get(contig)), 2
            )
            value_hits = math.sqrt(proportion_of_hits) * 5  # max value 5
            dict_matrix[contig][2] = value_hits

        if contig in vodb_contigs:
            dict_matrix[contig][3] = 2
    return dict_matrix


def get_matrix_sum(dict_matrix):
    """Adds together the scores of the matrix for each contig and stores the results as a dictionary."""
    possible_gv = dict()
    for contig, values in dict_matrix.items():
        sum_values = np.sum(values)
        sum_values = round(sum_values, 2)
        possible_gv[contig] = float(sum_values)
    return possible_gv


def create_table(
    possible_gv,
    contig_length_gc,
    contig_coverage,
    count_genes_per_contig,
    counter_other_GVOG_hits,
    counter_core3_hits,
    counter_ext_core7_hits,
    vogdb_contigs,
):
    """Creates the score table as a list of lists. Each row is a contig containing information about the contig_name, contig length, GC-content, average read coverage, number of predicted genes per contig, pergenctage of GVOG hits compared to all predicted genes, the score, count of core3, extcore7 and other GVOG hits as well as a VOGDB hit."""
    possible_gv_list = []
    for contig, score in possible_gv.items():
        length_gc = contig_length_gc.get(contig)
        length = int(length_gc[0])
        if length < min_length:
            continue
        gc = length_gc[1]
        all_genes = count_genes_per_contig.get(contig)
        if all_genes == None:
            all_genes = 0
        cov = contig_coverage.get(contig)

        if contig in counter_core3_hits:
            core_count = counter_core3_hits.get(contig)
        if contig not in counter_core3_hits:
            core_count = 0

        if contig in counter_ext_core7_hits:
            core_ext_count = counter_ext_core7_hits.get(contig)
        if contig not in counter_ext_core7_hits:
            core_ext_count = 0

        if contig in counter_other_GVOG_hits:
            all_hits_count = counter_other_GVOG_hits.get(contig)
        if contig not in counter_other_GVOG_hits:
            all_hits_count = 0

        if contig in vogdb_contigs:
            vogdb = "Hit"
        if contig not in vogdb_contigs:
            vogdb = "None"

        sum_of_gvogs = core_count + core_ext_count + all_hits_count
        if all_genes != 0:
            perc_gvog_all_genes = round(
                (float(sum_of_gvogs) / float(all_genes)) * 100, 2
            )
        else:
            perc_gvog_all_genes = 0
        possible_gv_list.append(
            [
                contig,
                length,
                gc,
                cov,
                all_genes,
                perc_gvog_all_genes,
                score,
                core_count,
                core_ext_count,
                all_hits_count,
                vogdb,
            ]
        )
    return possible_gv_list


def create_pd_df(possible_gv_list):
    """Converts the score table into a pandas DataFrame and sorts by score and contig length in descending order."""
    score_table = pd.DataFrame(
        possible_gv_list,
        columns=[
            "contig_name",
            "contig_len",
            "GC%",
            "average_coverage",
            "all_genes",
            "perc_GVOG_hits/all_genes",
            "score",
            "core3_count",
            "ext_core7_count",
            "other_GVOG_count",
            "VOGDB_hit",
        ],
    )
    score_table_sorted = score_table.sort_values(
        by=["score", "contig_len"], ascending=[False, False]
    )
    return score_table_sorted


if __name__ == "__main__":
    count_genes_per_contig_file = sys.argv[1]
    contig_length_gc_file = sys.argv[2]
    contig_coverage_file = sys.argv[3]
    core3_per_gene_file = sys.argv[4]
    ext_core7_per_gene_file = sys.argv[5]
    other_hits_per_gene_file = sys.argv[6]
    vogdb_file = sys.argv[7]
    contig_ids_without_rRNA_file = sys.argv[8]
    score_table = sys.argv[9]
    min_length = int(sys.argv[10])

    count_genes_per_contig = read_count_genes_per_contig(count_genes_per_contig_file)
    contig_length_gc = get_seq_length_and_gc(contig_length_gc_file)
    contig_coverage = get_coverage(contig_coverage_file)

    all_hits_per_contig = read_hit_contig_file(other_hits_per_gene_file)
    counter_other_GVOG_hits = count_hits_per_contig(all_hits_per_contig)

    core_hits_per_contig = read_hit_contig_file(core3_per_gene_file)
    counter_core3_hits = count_hits_per_contig(core_hits_per_contig)

    ext_core_hits_per_contig = read_hit_contig_file(ext_core7_per_gene_file)
    counter_ext_core7_hits = count_hits_per_contig(ext_core_hits_per_contig)

    vogdb_contigs = get_vogdb_hits(vogdb_file)

    all_contigs = read_contig_ids_without_rRNA(contig_ids_without_rRNA_file)
    dict_matrix = create_matrix(
        all_contigs,
        counter_other_GVOG_hits,
        counter_core3_hits,
        counter_ext_core7_hits,
        vogdb_contigs,
        count_genes_per_contig,
    )
    possible_gv = get_matrix_sum(dict_matrix)
    possible_gv_list = create_table(
        possible_gv,
        contig_length_gc,
        contig_coverage,
        count_genes_per_contig,
        counter_other_GVOG_hits,
        counter_core3_hits,
        counter_ext_core7_hits,
        vogdb_contigs,
    )
    possible_gv_df_sorted = create_pd_df(possible_gv_list)
    possible_gv_df_sorted.to_csv(score_table, sep="\t", header=True, index=False)
