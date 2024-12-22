import yaml
import os
import sys


snakefile_dir = workflow.basedir
# print("snakefile_dir: ", snakefile_dir)


configfile: "config/config.yaml"


with open("config/resources.yaml", "r") as stream:
    resources = yaml.safe_load(stream)


include: "rules/definitions_configurations.smk"
include: "rules/qc_short_reads.smk"
include: "rules/qc_long_reads.smk"
include: "rules/assemblies.smk"
include: "rules/assembly_qc.smk"
include: "rules/gene_prediction.smk"
include: "rules/make_gv_table.smk"


rule all:
    input:
        qc_short_reads,
        qc_long_reads,
        bam_and_cov_files=expand(
            analysis_dir + "assembly_qc/{samples}_cor_sorted.bam",
            samples=mapping_files,
            file_type=["_cor_sorted.bam", ".cov"],
        ),
        bbmap_stats=expand(
            analysis_dir + "assembly_qc/bbmap/{sample}_bbmap_stats.txt",
            sample=combined_bases,
        ),
        gaas=expand(
            analysis_dir + "assembly_qc/gaas/{sample}_gaas_stats",
            sample=combined_bases,
        ),
        dist_plot=expand(
            analysis_dir + "assembly_qc/dist_plots/{sample}_dist_plot.png",
            sample=combined_bases,
        ),
        coverm_files=expand(
            analysis_dir + "assembly_qc/coverm/{sample}_coverm.txt",
            sample=mapping_files,
        ),
        contig_len_gc=expand(
            analysis_dir + "assembly_qc/contig_len_gc/{sample}_contig_len_gc.txt",
            sample=combined_bases,
        ),
        gv_table=expand(
            analysis_dir + "gv_table/{sample}/{sample}_gv_table.txt",
            sample=combined_bases,
        ),
        gv_table_stats=expand(
            analysis_dir + "gv_table/{sample}/{sample}_gv_table_stats.txt",
            sample=combined_bases,
        ),
