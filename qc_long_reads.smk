rule fastqc_raw_long:
    resources:
        mem_mb=resources["fastqc"]["mem_mb"],
        runtime=resources["fastqc"]["runtime"],
        cpus_per_task=resources["fastqc"]["cpus_per_task"],
    threads: int(resources["fastqc"]["threads"])
    log:
        out=analysis_dir + "logs/read_qc/raw/{sample_long}_out.log",
        err=analysis_dir + "logs/read_qc/raw/{sample_long}_err.log",
    conda:
        snakefile_dir + "/envs/read_qc.yaml"
    input:
        reads=input_long + "/{sample_long}." + sample_extension_long,
    output:
        zip=analysis_dir + "read_qc/raw/{sample_long}_fastqc.zip",
        html=analysis_dir + "read_qc/raw/{sample_long}_fastqc.html",
    params:
        outdir=analysis_dir + "read_qc/raw",
    shell:
        """
        fastqc {input.reads} -t {threads} --outdir={params.outdir} > {log.out} 2> {log.err}
        """


rule nanoplot_raw:
    resources:
        mem_mb=resources["nanoplot"]["mem_mb"],
        runtime=resources["nanoplot"]["runtime"],
        cpus_per_task=resources["nanoplot"]["cpus_per_task"],
    threads: int(resources["nanoplot"]["threads"])
    log:
        out=analysis_dir + "logs/nanoplot/raw/{sample_long}_out.log",
        err=analysis_dir + "logs/nanoplot/raw/{sample_long}_err.log",
    conda:
        snakefile_dir + "/envs/nanoplot.yaml"
    input:
        reads=input_long + "/{sample_long}." + sample_extension_long,
    output:
        report=analysis_dir
        + "nanoplot/raw/{sample_long}/{sample_long}_NanoPlot-report.html",
    params:
        outdir=analysis_dir + "nanoplot/raw/{sample_long}",
    shell:
        """
        NanoPlot --fastq {input.reads} --outdir {params.outdir} --prefix {wildcards.sample_long}_ -t {threads} --verbose > {log.out} 2> {log.err}
        """


rule porechop:
    resources:
        mem_mb=resources["porechop"]["mem_mb"],
        runtime=resources["porechop"]["runtime"],
        cpus_per_task=resources["porechop"]["cpus_per_task"],
    threads: int(resources["porechop"]["threads"])
    log:
        out=analysis_dir + "logs/porechop/{sample_long}_out.log",
        err=analysis_dir + "logs/porechop/{sample_long}_err.log",
    conda:
        snakefile_dir + "/envs/read_qc.yaml"
    input:
        reads=input_long + "/{sample_long}." + sample_extension_long,
    output:
        reads=(
            analysis_dir + "porechop/{sample_long}.fastq.gz"
            if config["keep_intermediate_files"] == True
            else temp(analysis_dir + "porechop/{sample_long}.fastq.gz")
        ),
    shell:
        """
        porechop -i {input.reads} -o {output.reads} -t {threads} --verbosity 1 > {log.out} 2> {log.err}
        """


rule chopper:
    resources:
        mem_mb=resources["chopper"]["mem_mb"],
        runtime=resources["chopper"]["runtime"],
        cpus_per_task=resources["chopper"]["cpus_per_task"],
    threads: int(resources["chopper"]["threads"])
    log:
        err=analysis_dir + "logs/chopper/{sample_long}_err.log",
    conda:
        snakefile_dir + "/envs/read_qc.yaml"
    input:
        reads=(
            analysis_dir + "porechop/{sample_long}.fastq.gz"
            if config["long_read_type"]["Nanopore"] == True
            else input_long + "/{sample_long}.fastq.gz"
        ),
    output:
        reads=(
            analysis_dir + "chopper/{sample_long}.fastq.gz"
            if config["keep_intermediate_files"] == True
            else temp(analysis_dir + "chopper/{sample_long}.fastq.gz")
        ),
    params:
        phred=phred_cut_off_long,
    shell:
        """
        chopper -l 1000 -q {params.phred} -t {threads} -i {input.reads} | gzip -c -9 > {output.reads} 2> {log.err}
        """


rule fastqc_trimmed_long:
    resources:
        mem_mb=resources["fastqc"]["mem_mb"],
        runtime=resources["fastqc"]["runtime"],
        cpus_per_task=resources["fastqc"]["cpus_per_task"],
    threads: int(resources["fastqc"]["threads"])
    log:
        out=analysis_dir + "logs/read_qc/trimmed/{sample_long}_out.log",
        err=analysis_dir + "logs/read_qc/trimmed/{sample_long}_err.log",
    conda:
        snakefile_dir + "/envs/read_qc.yaml"
    input:
        reads=analysis_dir + "chopper/{sample_long}.fastq.gz",
    output:
        zip=analysis_dir + "read_qc/trimmed/{sample_long}_fastqc.zip",
        html=analysis_dir + "read_qc/trimmed/{sample_long}_fastqc.html",
    params:
        outdir=analysis_dir + "read_qc/trimmed",
    shell:
        """
        fastqc {input.reads} -t {threads} --outdir={params.outdir} > {log.out} 2> {log.err}
        """


rule nanoplot_trimmed:
    resources:
        mem_mb=resources["nanoplot"]["mem_mb"],
        runtime=resources["nanoplot"]["runtime"],
        cpus_per_task=resources["nanoplot"]["cpus_per_task"],
    threads: int(resources["nanoplot"]["threads"])
    log:
        out=analysis_dir + "logs/nanoplot/trimmed/{sample_long}_out.log",
        err=analysis_dir + "logs/nanoplot/trimmed/{sample_long}_err.log",
    conda:
        snakefile_dir + "/envs/nanoplot.yaml"
    input:
        reads=analysis_dir + "chopper/{sample_long}.fastq.gz",
    output:
        report=analysis_dir
        + "nanoplot/trimmed/{sample_long}/{sample_long}_NanoPlot-report.html",
    params:
        outdir=analysis_dir + "nanoplot/trimmed/{sample_long}",
    shell:
        """
        NanoPlot --fastq {input.reads} --outdir {params.outdir} --prefix {wildcards.sample_long}_ -t {threads} --verbose > {log.out} 2> {log.err}
        """
