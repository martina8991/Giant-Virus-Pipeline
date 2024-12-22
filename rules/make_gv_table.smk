
rule gene_count:
    log:
        err=analysis_dir + "logs/gene_count/{sample}_err.log",
    group:
        "gv_table_details"
    input:
        faa_file=analysis_dir + "prodigal/{sample}/{sample}.faa",
    output:
        gene_count=temp(
            analysis_dir + "gv_table_details/{sample}/{sample}_gene_count.txt"
        ),
    shell:
        """
        grep ">" {input.faa_file} | cut -f 1 -d ' ' | cut -f 2 -d '>' | sort | uniq | cut -f 1-2 -d '_' | sort | uniq -c | tr -s ' ' | sort > {output.gene_count} 2> {log.err}
        """


rule gv_table_details:
    log:
        out=analysis_dir + "logs/gv_table_details/{sample}_out.log",
        err=analysis_dir + "logs/gv_table_details/{sample}_err.log",
    conda:
        snakefile_dir + "/envs/assembly_qc.yaml"
    group:
        "gv_table_details"
    input:
        remerged_gvdb_sorted_filtered=analysis_dir
        + "hmmsearch/hmmsearch_GVDB/{sample}/per_sequence_hits.txt",
        remerged_vogdb_sorted_filtered=analysis_dir
        + "hmmsearch/hmmsearch_VOGDB/{sample}/per_sequence_hits.txt",
    output:
        core3=analysis_dir + "gv_table_details/{sample}/{sample}_core3_GVOG_hits.txt",
        ext_core7=analysis_dir
        + "gv_table_details/{sample}/{sample}_ext_core7_GVOG_hits.txt",
        other_GVOG_hits=analysis_dir
        + "gv_table_details/{sample}/{sample}_other_GVOG_hits.txt",
        vogdb_hits=analysis_dir + "gv_table_details/{sample}/{sample}_vogdb_hits.txt",
    params:
        GVDB_evalue=config["hmmsearch_GVDB_evalue"],
        VOGDB_evalue=config["hmmsearch_VOGDB_evalue"],
        filter_GVOG_script=snakefile_dir + "/scripts/filter_GVOG_hits.py",
        filter_VOG_script=snakefile_dir + "/scripts/filter_VFAM_hits.py",
        vfam_ids=snakefile_dir + "/databases/VOGDB/vfam.lca.tsv",
    shell:
        """
        python {params.filter_GVOG_script} {input.remerged_gvdb_sorted_filtered} {output.core3} \
        {output.ext_core7} {output.other_GVOG_hits} {params.GVDB_evalue} > {log.out} 2> {log.err}
        python {params.filter_VOG_script} {input.remerged_vogdb_sorted_filtered} {params.vfam_ids} {output.vogdb_hits} {params.VOGDB_evalue} >> {log.out} 2>> {log.err}
        """


def get_mapping_dir(wildcards):
    if wildcards.sample in base_single_paired_samples:
        dir = analysis_dir + "assembly_qc/coverm/bowtie2/{sample}_coverm.txt"
        return dir
    if wildcards.sample in sample_base_long_reads:
        dir = analysis_dir + "assembly_qc/coverm/minimap2/{sample}_coverm.txt"
        return dir
    if wildcards.sample in hybrid_metaspades_base:
        dir = analysis_dir + "assembly_qc/coverm/bowtie2/{sample}_coverm.txt"
        return dir


rule make_gv_table:
    default_target: True
    log:
        out=analysis_dir + "logs/gv_table/{sample}_out.log",
        err=analysis_dir + "logs/gv_table/{sample}_err.log",
    conda:
        snakefile_dir + "/envs/assembly_qc.yaml"
    input:
        gene_count=analysis_dir + "gv_table_details/{sample}/{sample}_gene_count.txt",
        len_gc=analysis_dir + "assembly_qc/contig_len_gc/{sample}_contig_len_gc.txt",
        cov=get_mapping_dir,
        core3=analysis_dir + "gv_table_details/{sample}/{sample}_core3_GVOG_hits.txt",
        ext_core7=analysis_dir
        + "gv_table_details/{sample}/{sample}_ext_core7_GVOG_hits.txt",
        other_GVOG_hits=analysis_dir
        + "gv_table_details/{sample}/{sample}_other_GVOG_hits.txt",
        vogdb_hits=analysis_dir + "gv_table_details/{sample}/{sample}_vogdb_hits.txt",
        contigs_without_rRNA=analysis_dir
        + "contigs_without_rRNA/{sample}_non_rRNA_contig_ids.txt",
    output:
        gv_table=analysis_dir + "gv_table/{sample}/{sample}_gv_table.txt",
        gv_table_stats=analysis_dir + "gv_table/{sample}/{sample}_gv_table_stats.txt",
    params:
        min_length=config["min_length_contigs_for_gv_table"],
        make_gv_table_script=snakefile_dir + "/scripts/score_table.py",
        gv_stats_script=snakefile_dir + "/scripts/score_table_stats.py",
    shell:
        """
        python {params.make_gv_table_script} {input.gene_count} {input.len_gc} {input.cov} {input.core3} {input.ext_core7} \
        {input.other_GVOG_hits} {input.vogdb_hits} {input.contigs_without_rRNA} {output.gv_table} {params.min_length} > {log.out} 2> {log.err}
        python {params.gv_stats_script} {output.gv_table} {output.gv_table_stats} >> {log.out} 2>> {log.err}
        """
