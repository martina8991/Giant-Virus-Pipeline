import pandas as pd
import sys

"""This script takes the score table as input and creates basic statistics of the score table as a csv output file.
The statistics include:
average number of annotated genes per contig
average contig length
maximum and minimum contig length
number of contigs with GVDB core3, extcore7 and other hits
number of contigs with hits against the VOGDB
etc.
"""

input_score_table = sys.argv[1]
score_table_stats_output = sys.argv[2]


score_table_df = pd.read_csv(input_score_table, sep="\t", header=0)
score_table_df.fillna("None", inplace=True)

if score_table_df.shape[0] != 0:
    row_count = len(score_table_df)
    average_ann_prot = round(score_table_df["all_genes"].mean())

    average_contig_len = round(score_table_df["contig_len"].mean())
    max_contig_len = score_table_df["contig_len"].max()
    min_contig_len = score_table_df["contig_len"].min()

    core3_hits = score_table_df[score_table_df["core3_count"] != 0].shape[0]
    ext_core7_hits = score_table_df[score_table_df["ext_core7_count"] != 0].shape[0]
    other_GVOG_hits = score_table_df[score_table_df["other_GVOG_count"] != 0].shape[0]

    num_cont_with_GVOG_hits = score_table_df[
        (score_table_df["core3_count"] != 0)
        | (score_table_df["ext_core7_count"] != 0)
        | (score_table_df["other_GVOG_count"] != 0)
    ].shape[0]

    num_cont_with_VOGDB_hits = score_table_df[score_table_df["VOGDB_hit"] == "Hit"].shape[0]

    num_cont_GVOG_hit_VOGDB_none = score_table_df[
        (score_table_df["VOGDB_hit"] != "Hit")
        & (
            (
                (score_table_df["core3_count"] != 0)
                | (score_table_df["ext_core7_count"] != 0)
                | (score_table_df["other_GVOG_count"] != 0)
            )
        )
    ].shape[0]

    num_cont_VOGDB_hit_GVOG_hit = score_table_df[
        (score_table_df["VOGDB_hit"] == "Hit")
        & (
            (
                (score_table_df["core3_count"] != 0)
                | (score_table_df["ext_core7_count"] != 0)
                | (score_table_df["other_GVOG_count"] != 0)
            )
        )
    ].shape[0]

    num_cont_VOGDB_hit_GVOG_none = score_table_df[
        (score_table_df["VOGDB_hit"] == "Hit")
        & (score_table_df["core3_count"] == 0)
        & (score_table_df["ext_core7_count"] == 0)
        & (score_table_df["other_GVOG_count"] == 0)
    ].shape[0]


if score_table_df.shape[0] == 0:
    row_count = 0
    average_ann_prot = 0

    average_contig_len = 0
    max_contig_len = 0
    min_contig_len = 0

    core3_hits = 0
    ext_core7_hits = 0
    other_GVOG_hits = 0

    num_cont_with_GVOG_hits = 0
    num_cont_with_VOGDB_hits = 0
    num_cont_GVOG_hit_VOGDB_none = 0
    num_cont_VOGDB_hit_GVOG_hit = 0
    num_cont_VOGDB_hit_GVOG_none = 0


score_table_stats_df = pd.DataFrame(
    [
        ["number of contigs:", row_count],
        ["average number of annotated genes per contig:", average_ann_prot],
        ["average contig length:", average_contig_len],
        ["maximum contig length:", max_contig_len],
        ["minimum contig length:", min_contig_len],
        ["number of contigs with GVDB core3 hits:", core3_hits],
        ["number of contigs with GVDB core7 hits:", ext_core7_hits],
        ["number of contigs with other GVDB hits:", other_GVOG_hits],
        ["number of all contigs with GVDB hits:", num_cont_with_GVOG_hits],
        ["number of all contigs with VOGDB hits:", num_cont_with_VOGDB_hits],
        [
            "number of all contigs with GVDB hits but no VOGDB hit:",
            num_cont_GVOG_hit_VOGDB_none,
        ],
        [
            "number of all contigs with both GVDB hits and VOGDB hit:",
            num_cont_VOGDB_hit_GVOG_hit,
        ],
        [
            "number of all contigs with no GVDB hits but VOGDB hit:",
            num_cont_VOGDB_hit_GVOG_none,
        ],
    ],
    columns=["Metric", "Value"],
)


score_table_stats_df.to_csv(score_table_stats_output, sep="\t", index=False)
