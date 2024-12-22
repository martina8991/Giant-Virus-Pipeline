
rule fastqc_raw_short:
    resources:
        mem_mb=resources["fastqc"]["mem_mb"],
        runtime=resources["fastqc"]["runtime"],
        cpus_per_task=resources["fastqc"]["cpus_per_task"],
    threads: int(resources["fastqc"]["threads"])
    log:
        out=analysis_dir + "logs/read_qc/raw/{sample_short}_out.log",
        err=analysis_dir + "logs/read_qc/raw/{sample_short}_err.log",
    conda:
        snakefile_dir + "/envs/read_qc.yaml"
    input:
        reads=input_short + "/{sample_short}." + sample_extension_short,
    output:
        zip=analysis_dir + "read_qc/raw/{sample_short}_fastqc.zip",
        html=analysis_dir + "read_qc/raw/{sample_short}_fastqc.html",
    params:
        outdir=analysis_dir + "read_qc/raw",
    shell:
        """
        fastqc '{input.reads}' -t {threads} --outdir={params.outdir} > {log.out} 2> {log.err}
        """


rule fastp_paired_end:
    resources:
        mem_mb=resources["fastp"]["mem_mb"],
        runtime=resources["fastp"]["runtime"],
        cpus_per_task=resources["fastp"]["cpus_per_task"],
    threads: int(resources["fastp"]["threads"])
    log:
        out=analysis_dir + "logs/fastp/{sample_short}_out.log",
        err=analysis_dir + "logs/fastp/{sample_short}_err.log",
    conda:
        snakefile_dir + "/envs/read_qc.yaml"
    input:
        file_forward=input_short
        + "/{sample_short}"
        + sample_forward_id
        + "."
        + sample_extension_short,
        file_reverse=input_short
        + "/{sample_short}"
        + sample_reverse_id
        + "."
        + sample_extension_short,
    output:
        paired_forwardP=(
            analysis_dir + "fastp/{sample_short}" + sample_forward_id + "_P.fastq.gz"
            if config["keep_intermediate_files"] == True
            else temp(
                analysis_dir
                + "fastp/{sample_short}"
                + sample_forward_id
                + "_P.fastq.gz"
            )
        ),
        paired_reverseP=(
            analysis_dir + "fastp/{sample_short}" + sample_reverse_id + "_P.fastq.gz"
            if config["keep_intermediate_files"] == True
            else temp(
                analysis_dir
                + "fastp/{sample_short} "
                + sample_reverse_id
                + "_P.fastq.gz"
            )
        ),
        paired_forwardU=(
            analysis_dir + "fastp/{sample_short}" + sample_forward_id + "_U.fastq.gz"
            if config["keep_intermediate_files"] == True
            else temp(
                analysis_dir
                + "fastp/{sample_short}"
                + sample_forward_id
                + "_U.fastq.gz"
            )
        ),
        paired_reverseU=(
            analysis_dir + "fastp/{sample_short}" + sample_reverse_id + "_U.fastq.gz"
            if config["keep_intermediate_files"] == True
            else temp(
                analysis_dir
                + "fastp/{sample_short}"
                + sample_reverse_id
                + "_U.fastq.gz"
            )
        ),
    params:
        phred=phred_cut_off_short,
    shell:
        """
        fastp --in1 {input.file_forward} --in2 {input.file_reverse} --out1 {output.paired_forwardP} --out2 {output.paired_reverseP} \
        --unpaired1 {output.paired_forwardU} --unpaired2 {output.paired_reverseU} --detect_adapter_for_pe --trim_poly_g --compression 9 --overrepresentation_analysis \
        --dont_overwrite --average_qual {params.phred} --length_required 50 --cut_front --cut_front_window_size 1 --cut_front_mean_quality {params.phred} --cut_tail \
        --cut_tail_window_size 1 --cut_tail_mean_quality {params.phred} --cut_right --cut_right_mean_quality {params.phred} --trim_tail1 2 --trim_tail2 2 --thread {threads} \
        --json {analysis_dir}fastp/{wildcards.sample_short}.json --html {analysis_dir}fastp/{wildcards.sample_short}.html
         """


def get_fastqc_trimmed_file(wildcards):
    if wildcards.sample_short in paired_end_base:
        ending = expand(
            "{id}{end}.fastq.gz",
            id=[sample_forward_id, sample_reverse_id],
            end=["P", "U"],
        )
        file_list = [
            analysis_dir + "fastp/" + wildcards.sample_short + i for i in ending
        ]
        return file_list
    elif wildcards.sample_short in single_end_base:
        file = analysis_dir + "fastp/" + wildcards.sample_short + ".fastq.gz"
        return file


rule fastp_single_end:
    resources:
        mem_mb=resources["fastp"]["mem_mb"],
        runtime=resources["fastp"]["runtime"],
        cpus_per_task=resources["fastp"]["cpus_per_task"],
    threads: int(resources["fastp"]["threads"])
    log:
        out=analysis_dir + "logs/fastp/{sample_short}_out.log",
        err=analysis_dir + "logs/fastp/{sample_short}_err.log",
    conda:
        snakefile_dir + "/envs/read_qc.yaml"
    input:
        single_reads=input_short + "/{sample_short}." + sample_extension_short,
    output:
        trimmed_se=(
            analysis_dir + "fastp/{sample_short}.fastq.gz"
            if config["keep_intermediate_files"] == True
            else temp(analysis_dir + "fastp/{sample_short}.fastq.gz")
        ),
    params:
        phred=phred_cut_off_short,
    shell:
        """
        fastp -i {input.single_reads} -o {output.trimmed_se} --trim_poly_g --compression 9 --dont_overwrite --average_qual {params.phred} --overrepresentation_analysis \
        --length_required 50 --cut_front --cut_front_window_size 1 --cut_front_mean_quality {params.phred} --cut_tail --cut_tail_window_size 1 \
        --cut_tail_mean_quality {params.phred} --cut_right --cut_right_mean_quality {params.phred} --trim_tail1 2 --trim_tail2 2 --thread {threads} \
        --json {analysis_dir}fastp/{wildcards.sample_short}.json --html {analysis_dir}fastp/{wildcards.sample_short}.html
        """


ruleorder: fastp_paired_end > fastp_single_end


rule fastqc_trimmed_short:
    resources:
        mem_mb=resources["fastqc"]["mem_mb"],
        runtime=resources["fastqc"]["runtime"],
        cpus_per_task=resources["fastqc"]["cpus_per_task"],
    threads: int(resources["fastqc"]["threads"])
    log:
        out=analysis_dir + "logs/read_qc/trimmed/{sample_short}_out.log",
        err=analysis_dir + "logs/read_qc/trimmed/{sample_short}_err.log",
    conda:
        snakefile_dir + "/envs/read_qc.yaml"
    input:
        reads=analysis_dir + "fastp/{sample_short}.fastq.gz",
    output:
        zip=analysis_dir + "read_qc/trimmed/{sample_short}_fastqc.zip",
        html=analysis_dir + "read_qc/trimmed/{sample_short}_fastqc.html",
    params:
        outdir=analysis_dir + "read_qc/trimmed",
    shell:
        """
        fastqc '{input.reads}' -t {threads} --outdir={params.outdir} > {log.out} 2> {log.err}
        """
