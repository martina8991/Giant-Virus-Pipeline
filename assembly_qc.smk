def get_assembly_file(wildcards):
    if wildcards.sample in base_single_paired_samples:
        return (
            analysis_dir
            + short_read_assembler
            + "/{sample}/"
            + short_read_assembler_file
        )
    if wildcards.sample in sample_base_long_reads:
        return analysis_dir + "metaflye/{sample}/assembly.fasta"
    if wildcards.sample in hybrid_metaspades_base:
        return analysis_dir + "hybrid_metaspades/{sample}/scaffolds.fasta"


def get_assembly_dir(wildcards):
    if wildcards.sample in base_single_paired_samples:
        return analysis_dir + short_read_assembler + "/{sample}"
    if wildcards.sample in sample_base_long_reads:
        return analysis_dir + "metaflye/{sample}"
    if wildcards.sample in hybrid_metaspades_base:
        return analysis_dir + "hybrid_metaspades/{sample}"


rule filter_contigs_post_assembly:
    log:
        err=analysis_dir + "logs/filter_cont/{sample}_err.log",
    conda:
        snakefile_dir + "/envs/assembly_qc.yaml"
    input:
        assembly=get_assembly_file,
        assembly_dir=get_assembly_dir,
    output:
        assembly_filtered=analysis_dir + "filtered_assemblies/{sample}.fa",
    shell:
        """
        seqtk seq -L 1000 {input.assembly} > {output.assembly_filtered} 2> {log.err}
        """


def is_trimmed(qc_step):
    if qc_step == True:
        return "P"
    else:
        return ""


rule bowtie2_paired_end_short_read:
    resources:
        mem_mb=resources["bowtie2"]["mem_mb"],
        runtime=resources["bowtie2"]["runtime"],
        cpus_per_task=resources["bowtie2"]["cpus_per_task"],
    threads: int(resources["bowtie2"]["threads"])
    log:
        out=analysis_dir + "logs/assembly_qc/bowtie2/{sample_short}_out.log",
        err=analysis_dir + "logs/assembly_qc/bowtie2/{sample_short}_err.log",
    conda:
        snakefile_dir + "/envs/mapping_tools.yaml"
    input:
        forward_paired=get_forward_paired,
        reverse_paired=get_reverse_paired,
        assembly=analysis_dir + "filtered_assemblies/{sample_short}.fa",
    output:
        index_dir=temp(
            directory(analysis_dir + "assembly_qc/bowtie2/{sample_short}_index")
        ),
        sam=temp(analysis_dir + "assembly_qc/bowtie2/{sample_short}.sam"),
        bam=temp(analysis_dir + "assembly_qc/bowtie2/{sample_short}.bam"),
        bam_cor_sorted=analysis_dir
        + "assembly_qc/bowtie2/{sample_short}_cor_sorted.bam",
        cov=analysis_dir + "assembly_qc/bowtie2/{sample_short}.cov",
    params:
        forward_unpaired=get_forward_unpaired(sample_name_short),
        reverse_unpaired=get_reverse_unpaired(),
    shell:
        """
        mkdir -p {output.index_dir} > {log.out} 2> {log.err}
        bowtie2-build {input.assembly} {output.index_dir}/index --large-index -p {threads} >> {log.out} 2>> {log.err}
        if [[ -n "{params.forward_unpaired}" ]]; then
            echo "Unpaired reads provided"
            unpaired_merged={analysis_dir}{short_read_assembler}/{wildcards.sample_short}/unpaired_merged.{short_read_extension}
            echo $unpaired_merged
            bowtie2 -x {output.index_dir}/index -1 '{input.forward_paired}' -2 '{input.reverse_paired}' -U ${{unpaired_merged}} \
            -S {output.sam} -p {threads} >> {log.out} 2>> {log.err}
            rm $unpaired_merged
        else
            echo "No unpaired reads provided"
            bowtie2 -x {output.index_dir}/index -1 '{input.forward_paired}' -2 '{input.reverse_paired}' -S {output.sam} -p {threads} >> {log.out} 2>> {log.err}
        fi
        samtools view -b {output.sam} -o {output.bam} --threads {threads} >> {log.out} 2>> {log.err}
        samtools sort {output.bam} -o {output.bam_cor_sorted} --threads {threads} >> {log.out} 2>> {log.err}
        samtools coverage {output.bam_cor_sorted} -D -o {output.cov} >> {log.out} 2>> {log.err}
        """


rule bowtie2_paired_end_hybrid_metaspades:
    resources:
        mem_mb=resources["bowtie2"]["mem_mb"],
        runtime=resources["bowtie2"]["runtime"],
        cpus_per_task=resources["bowtie2"]["cpus_per_task"],
    threads: int(resources["bowtie2"]["threads"])
    log:
        out=analysis_dir
        + "logs/assembly_qc/bowtie2/{sample_short}-{sample_long}_out.log",
        err=analysis_dir
        + "logs/assembly_qc/bowtie2/{sample_short}-{sample_long}_err.log",
    conda:
        snakefile_dir + "/envs/mapping_tools.yaml"
    input:
        forward_paired=get_forward_paired,
        reverse_paired=get_reverse_paired,
        assembly=analysis_dir + "filtered_assemblies/{sample_short}-{sample_long}.fa",
    output:
        index_dir=temp(
            directory(
                analysis_dir + "assembly_qc/bowtie2/{sample_short}-{sample_long}_index"
            )
        ),
        sam=temp(analysis_dir + "assembly_qc/bowtie2/{sample_short}-{sample_long}.sam"),
        bam=temp(analysis_dir + "assembly_qc/bowtie2/{sample_short}-{sample_long}.bam"),
        bam_cor_sorted=analysis_dir
        + "assembly_qc/bowtie2/{sample_short}-{sample_long}_cor_sorted.bam",
        cov=analysis_dir + "assembly_qc/bowtie2/{sample_short}-{sample_long}.cov",
    params:
        forward_unpaired=get_forward_unpaired(sample_name_short),
        reverse_unpaired=get_reverse_unpaired(),
    shell:
        """
        mkdir -p {output.index_dir} > {log.out} 2> {log.err}
        bowtie2-build {input.assembly} {output.index_dir}/index --large-index -p {threads} >> {log.out} 2>> {log.err}

        if [[ -n "{params.forward_unpaired}" ]]; then
            echo "Unpaired reads provided"
            unpaired_merged={analysis_dir}hybrid_metaspades/{wildcards.sample_short}-{wildcards.sample_long}/unpaired_merged.{short_read_extension}
            echo $unpaired_merged
            bowtie2 -x {output.index_dir}/index -1 '{input.forward_paired}' -2 '{input.reverse_paired}' -U ${{unpaired_merged}} \
            -S {output.sam} -p {threads} >> {log.out} 2>> {log.err}
            rm $unpaired_merged
        else
            echo "No unpaired reads provided"
            bowtie2 -x {output.index_dir}/index -1 '{input.forward_paired}' -2 '{input.reverse_paired}' -S {output.sam} -p {threads} >> {log.out} 2>> {log.err}
        fi
        samtools view -b {output.sam} -o {output.bam} --threads {threads} >> {log.out} 2>> {log.err}
        samtools sort {output.bam} -o {output.bam_cor_sorted} --threads {threads} >> {log.out} 2>> {log.err}
        samtools coverage {output.bam_cor_sorted} -D -o {output.cov} >> {log.out} 2>> {log.err}
        """


rule bowtie2_single_end:
    resources:
        mem_mb=resources["bowtie2"]["mem_mb"],
        runtime=resources["bowtie2"]["runtime"],
        cpus_per_task=resources["bowtie2"]["cpus_per_task"],
    threads: int(resources["bowtie2"]["threads"])
    log:
        out=analysis_dir + "logs/assembly_qc/bowtie2/{sample_short}_out.log",
        err=analysis_dir + "logs/assembly_qc/bowtie2/{sample_short}_err.log",
    conda:
        snakefile_dir + "/envs/mapping_tools.yaml"
    input:
        reads=short_read_input_dir + "/{sample_short}." + short_read_extension,
        assembly=analysis_dir + "filtered_assemblies/{sample_short}.fa",
    output:
        index_dir=temp(
            directory(analysis_dir + "assembly_qc/bowtie2/{sample_short}_index")
        ),
        sam=temp(analysis_dir + "assembly_qc/bowtie2/{sample_short}.sam"),
        bam=temp(analysis_dir + "assembly_qc/bowtie2/{sample_short}.bam"),
        bam_cor_sorted=analysis_dir
        + "assembly_qc/bowtie2/{sample_short}_cor_sorted.bam",
        cov=analysis_dir + "assembly_qc/bowtie2/{sample_short}.cov",
    shell:
        """
        mkdir -p {output.index_dir} > {log.out} 2> {log.err}
        bowtie2-build {input.assembly} {output.index_dir}/index --large-index -p {threads} >> {log.out} 2>> {log.err}
        bowtie2 -x {output.index_dir}/index -U {input.reads} -S {output.sam} -p {threads} >> {log.out} 2>> {log.err}
        samtools view -b {output.sam} -o {output.bam} --threads {threads} >> {log.out} 2>> {log.err}
        samtools sort {output.bam} -o {output.bam_cor_sorted} --threads {threads} >> {log.out} 2>> {log.err}
        samtools coverage {output.bam_cor_sorted} -D -o {output.cov} >> {log.out} 2>> {log.err}
        """


ruleorder: bowtie2_paired_end_short_read > bowtie2_single_end


def get_minimap_settings():
    if config["long_read_type"]["Nanopore"] == True:
        return "-ax map-ont"
    if config["long_read_type"]["PacBio CLR"] == True:
        return "-ax map-pb"
    if config["long_read_type"]["PacBio HiFi"] == True:
        return "-ax map-hifi"


rule minimap2_metaflye:
    resources:
        mem_mb=resources["minimap2"]["mem_mb"],
        runtime=resources["minimap2"]["runtime"],
        cpus_per_task=resources["minimap2"]["cpus_per_task"],
    threads: int(resources["minimap2"]["threads"])
    log:
        out=analysis_dir + "logs/assembly_qc/minimap2/{sample_long}_out.log",
        err=analysis_dir + "logs/assembly_qc/minimap2/{sample_long}_err.log",
    conda:
        snakefile_dir + "/envs/mapping_tools.yaml"
    input:
        assembly=analysis_dir + "filtered_assemblies/{sample_long}.fa",
        reads=metaflye_input_dir + "/{sample_long}." + metaflye_extension,
    output:
        index=temp(analysis_dir + "assembly_qc/minimap2/{sample_long}.mmi"),
        sam=temp(analysis_dir + "assembly_qc/minimap2/{sample_long}.sam"),
        bam=temp(analysis_dir + "assembly_qc/minimap2/{sample_long}.bam"),
        bam_cor_sorted=analysis_dir
        + "assembly_qc/minimap2/{sample_long}_cor_sorted.bam",
        cov=analysis_dir + "assembly_qc/minimap2/{sample_long}.cov",
    params:
        minimap_settings=get_minimap_settings(),
    shell:
        """
        minimap2 -d {output.index} {input.assembly} > {log.out} 2> {log.err}
        minimap2 {params.minimap_settings} {output.index} {input.reads} -o {output.sam} -t {threads} >> {log.out} 2>> {log.err}
        samtools view -b {output.sam} -o {output.bam} --threads {threads} >> {log.out} 2>> {log.err}
        samtools sort {output.bam} -o {output.bam_cor_sorted} --threads {threads} >> {log.out} 2>> {log.err}
        samtools coverage {output.bam_cor_sorted} -D -o {output.cov} >> {log.out} 2>> {log.err}
        """


rule minimap2_hybrid_metaspades:
    resources:
        mem_mb=resources["minimap2"]["mem_mb"],
        runtime=resources["minimap2"]["runtime"],
        cpus_per_task=resources["minimap2"]["cpus_per_task"],
    threads: int(resources["minimap2"]["threads"])
    log:
        out=analysis_dir
        + "logs/assembly_qc/minimap2/{sample_short}-{sample_long}_out.log",
        err=analysis_dir
        + "logs/assembly_qc/minimap2/{sample_short}-{sample_long}_err.log",
    conda:
        snakefile_dir + "/envs/mapping_tools.yaml"
    input:
        assembly=analysis_dir + "filtered_assemblies/{sample_short}-{sample_long}.fa",
        reads=metaflye_input_dir + "/{sample_long}." + metaflye_extension,
    output:
        index=temp(
            analysis_dir + "assembly_qc/minimap2/{sample_short}-{sample_long}.mmi"
        ),
        sam=temp(analysis_dir + "assembly_qc/minimap2/{sample_short}-{sample_long}.sam"),
        bam=temp(analysis_dir + "assembly_qc/minimap2/{sample_short}-{sample_long}.bam"),
        bam_cor_sorted=analysis_dir
        + "assembly_qc/minimap2/{sample_short}-{sample_long}_cor_sorted.bam",
        cov=analysis_dir + "assembly_qc/minimap2/{sample_short}-{sample_long}.cov",
    params:
        minimap_settings=get_minimap_settings(),
    shell:
        """
        minimap2 -d {output.index} {input.assembly} > {log.out} 2> {log.err}
        minimap2 {params.minimap_settings} {output.index} {input.reads} -o {output.sam} -t {threads} >> {log.out} 2>> {log.err}
        samtools view -b {output.sam} -o {output.bam} --threads {threads} >> {log.out} 2>> {log.err}
        samtools sort {output.bam} -o {output.bam_cor_sorted} --threads {threads} >> {log.out} 2>> {log.err}
        samtools coverage {output.bam_cor_sorted} -D -o {output.cov} >> {log.out} 2>> {log.err}
        """


rule coverm:
    log:
        err=analysis_dir + "logs/assembly_qc/coverm/{sample}_err.log",
    conda:
        snakefile_dir + "/envs/assembly_qc.yaml"
    group:
        "assembly_qc"
    input:
        map_file=analysis_dir + "assembly_qc/{sample}_cor_sorted.bam",
    output:
        cov=analysis_dir + "assembly_qc/coverm/{sample}_coverm.txt",
    shell:
        """
        coverm contig -b {input.map_file} -m mean covered_bases variance length count reads_per_base > {output.cov} 2> {log.err}
        """


rule contig_len_gc:
    log:
        err=analysis_dir + "logs/assembly_qc/contig_len_gc/{sample}_err.log",
    conda:
        snakefile_dir + "/envs/assembly_qc.yaml"
    group:
        "assembly_qc"
    input:
        assembly=analysis_dir + "filtered_assemblies/{sample}.fa",
    output:
        len_gc=analysis_dir + "assembly_qc/contig_len_gc/{sample}_contig_len_gc.txt",
        len=temp(analysis_dir + "assembly_qc/contig_len_gc/{sample}_contig_len.txt"),
    shell:
        """
        infoseq -auto -nocolumns -delimiter ' ' -only -noheading -name -length -pgc {input.assembly} > {output.len_gc} 2> {log.err}
        cut -f 1-2 -d ' ' {output.len_gc} > {output.len} 2>> {log.err}
        """


rule bbmap:
    log:
        err=analysis_dir + "logs/assembly_qc/bbmap/{sample}_err.log",
    conda:
        snakefile_dir + "/envs/assembly_qc.yaml"
    group:
        "assembly_qc"
    input:
        assembly=analysis_dir + "filtered_assemblies/{sample}.fa",
    output:
        stats=analysis_dir + "assembly_qc/bbmap/{sample}_bbmap_stats.txt",
    shell:
        """
        stats.sh in={input.assembly} > {output.stats} 2> {log.err}
        """


rule gaas:
    log:
        out=analysis_dir + "logs/assembly_qc/gaas/{sample}_out.log",
        err=analysis_dir + "logs/assembly_qc/gaas/{sample}_err.log",
    conda:
        snakefile_dir + "/envs/gaas.yaml"
    group:
        "assembly_qc"
    input:
        assembly=analysis_dir + "filtered_assemblies/{sample}.fa",
    output:
        stats=directory(analysis_dir + "assembly_qc/gaas/{sample}_gaas_stats"),
    shell:
        """
        gaas_fasta_statistics.pl --f {input.assembly} --output="{output.stats}" > {log.out} 2> {log.err}
        """


rule dist_plot:
    log:
        out=analysis_dir + "logs/assembly_qc/dist_plot/{sample}_out.log",
        err=analysis_dir + "logs/assembly_qc/dist_plot/{sample}_err.log",
    conda:
        snakefile_dir + "/envs/assembly_qc.yaml"
    input:
        len=analysis_dir + "assembly_qc/contig_len_gc/{sample}_contig_len.txt",
    output:
        plot=analysis_dir + "assembly_qc/dist_plots/{sample}_dist_plot.png",
    params:
        plot_script=snakefile_dir + "/scripts/cont_length_dist_plot.py",
    shell:
        """
        python {params.plot_script} {input.len} {output.plot} > {log.out} 2> {log.err}
        """
