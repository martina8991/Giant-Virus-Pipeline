

rule megahit_paired_end:
    resources:
        mem_mb=resources["megahit"]["mem_mb"],
        runtime=resources["megahit"]["runtime"],
        cpus_per_task=resources["megahit"]["cpus_per_task"],
    threads: int(resources["megahit"]["threads"])
    log:
        out=analysis_dir + "logs/megahit/{sample_short}_out.log",
        err=analysis_dir + "logs/megahit/{sample_short}_err.log",
    conda:
        snakefile_dir + "/envs/megahit.yaml"
    input:
        forward_paired=get_forward_paired,
        reverse_paired=get_reverse_paired,
    output:
        assembly=(
            analysis_dir + "megahit/{sample_short}/final.contigs.fa"
            if config["keep_intermediate_files"] == True
            else temp(analysis_dir + "megahit/{sample_short}/final.contigs.fa")
        ),
        assembly_tmp_dir=temp(directory(analysis_dir + "megahit/{sample_short}/tmp")),
        assembly_dir=(
            directory(analysis_dir + "megahit/{sample_short}")
            if config["keep_intermediate_files"] == True
            else temp(directory(analysis_dir + "megahit/{sample_short}"))
        ),
    params:
        forward_unpaired=get_forward_unpaired(sample_name_short),
        reverse_unpaired=get_reverse_unpaired(),
    shell:
        """
        if [[ -n "{params.forward_unpaired}" ]]; then
            echo "Unpaired reads are provided"
            forward_unpaired={short_read_input_dir}/{wildcards.sample_short}{params.forward_unpaired}
            echo $forward_unpaired
            reverse_unpaired={short_read_input_dir}/{wildcards.sample_short}{params.reverse_unpaired}
            if [[ -f $reverse_unpaired ]]; then
                echo $reverse_unpaired
                cat $forward_unpaired $reverse_unpaired > {analysis_dir}megahit/{wildcards.sample_short}/unpaired_merged.{short_read_extension}
            else
                cat $forward_unpaired > {analysis_dir}megahit/{wildcards.sample_short}/unpaired_merged.{short_read_extension} 
            fi
            unpaired_merged={analysis_dir}megahit/{wildcards.sample_short}/unpaired_merged.{short_read_extension}
            echo $unpaired_merged
            megahit -1 {input.forward_paired} -2 {input.reverse_paired} -r ${{unpaired_merged}} \
            --out-dir {output.assembly_tmp_dir} -t {threads} --verbose --min-contig-len 1000 > {log.out} 2> {log.err}
        else
            echo "No unpaired reads are provided"
            megahit -1 {input.forward_paired} -2 {input.reverse_paired} \
            --out-dir {output.assembly_tmp_dir} -t {threads} --verbose --min-contig-len 1000 > {log.out} 2> {log.err}
        fi
        mv -f {output.assembly_tmp_dir}/* {output.assembly_dir} >> {log.out} 2>> {log.err}
        """


rule megahit_single_end:
    resources:
        mem_mb=resources["megahit"]["mem_mb"],
        runtime=resources["megahit"]["runtime"],
        cpus_per_task=resources["megahit"]["cpus_per_task"],
    threads: int(resources["megahit"]["threads"])
    log:
        out=analysis_dir + "logs/megahit/{sample_short}_out.log",
        err=analysis_dir + "logs/megahit/{sample_short}_err.log",
    conda:
        snakefile_dir + "/envs/megahit.yaml"
    input:
        reads=short_read_input_dir + "/{sample_short}." + short_read_extension,
    output:
        assembly=(
            analysis_dir + "megahit/{sample_short}/final.contigs.fa"
            if config["keep_intermediate_files"] == True
            else temp(analysis_dir + "megahit/{sample_short}/final.contigs.fa")
        ),
        assembly_tmp_dir=temp(directory(analysis_dir + "megahit/{sample_short}/tmp")),
        assembly_dir=(
            directory(analysis_dir + "megahit/{sample_short}")
            if config["keep_intermediate_files"] == True
            else temp(directory(analysis_dir + "megahit/{sample_short}"))
        ),
    shell:
        """
        megahit -r {input.reads} --out-dir {output.assembly_tmp_dir} -t {threads} --verbose \
        --min-contig-len 1000 > {log.out} 2> {log.err}
        mv -f {output.assembly_tmp_dir}/* {output.assembly_dir} >> {log.out} 2>> {log.err}
        """


ruleorder: megahit_paired_end > megahit_single_end


rule metaspades_paired_end:
    resources:
        mem_mb=resources["metaspades"]["mem_mb"],
        runtime=resources["metaspades"]["runtime"],
        cpus_per_task=resources["metaspades"]["cpus_per_task"],
    threads: int(resources["metaspades"]["threads"])
    log:
        out=analysis_dir + "logs/metaspades/{sample_short}_out.log",
        err=analysis_dir + "logs/metaspades/{sample_short}_err.log",
    conda:
        snakefile_dir + "/envs/spades.yaml"
    input:
        forward_paired=get_forward_paired,
        reverse_paired=get_reverse_paired,
    output:
        assembly=(
            analysis_dir + "metaspades/{sample_short}/scaffolds.fasta"
            if config["keep_intermediate_files"] == True
            else temp(analysis_dir + "metaspades/{sample_short}/scaffolds.fasta")
        ),
        #assembly_tmp_dir=temp(directory(analysis_dir + "metaspades/{sample_short}/tmp")),
        assembly_dir=(
            directory(analysis_dir + "metaspades/{sample_short}")
            if config["keep_intermediate_files"] == True
            else temp(directory(analysis_dir + "metaspades/{sample_short}"))
        ),
    params:
        forward_unpaired=get_forward_unpaired(sample_name_short),
        reverse_unpaired=get_reverse_unpaired(),
    shell:
        """
        if [[ -n "{params.forward_unpaired}" ]]; then
            echo "Unpaired reads are provided"
            forward_unpaired={short_read_input_dir}/{wildcards.sample_short}{params.forward_unpaired}
            echo $forward_unpaired
            reverse_unpaired={short_read_input_dir}/{wildcards.sample_short}{params.reverse_unpaired}
            if [[ -f $reverse_unpaired ]]; then
                echo $reverse_unpaired
                cat $forward_unpaired $reverse_unpaired > {analysis_dir}metaspades/{wildcards.sample_short}/unpaired_merged.{short_read_extension}
            else
                cat $forward_unpaired > {analysis_dir}metaspades/{wildcards.sample_short}/unpaired_merged.{short_read_extension} 
            fi
            unpaired_merged={analysis_dir}metaspades/{wildcards.sample_short}/unpaired_merged.{short_read_extension}
            echo $unpaired_merged
            spades.py --meta -1 {input.forward_paired} -2 {input.reverse_paired} -s ${{unpaired_merged}} \
            -o {output.assembly_dir} -t {threads} -m {resources.mem_mb} > {log.out} 2> {log.err}
        else
            echo "No unpaired reads are provided"
            spades.py --meta -1 {input.forward_paired} -2 {input.reverse_paired} \
            -o {output.assembly_dir} -t {threads} -m {resources.mem_mb} > {log.out} 2> {log.err}
        fi
        """


rule metaspades_single_end:
    resources:
        mem_mb=resources["metaspades"]["mem_mb"],
        runtime=resources["metaspades"]["runtime"],
        cpus_per_task=resources["metaspades"]["cpus_per_task"],
    threads: int(resources["metaspades"]["threads"])
    log:
        out=analysis_dir + "logs/metaspades/{sample_short}_out.log",
        err=analysis_dir + "logs/metaspades/{sample_short}_err.log",
    conda:
        snakefile_dir + "/envs/spades.yaml"
    input:
        reads=short_read_input_dir + "/{sample_short}." + short_read_extension,
    output:
        assembly=(
            analysis_dir + "metaspades/{sample_short}/scaffolds.fasta"
            if config["keep_intermediate_files"] == True
            else temp(analysis_dir + "metaspades/{sample_short}/scaffolds.fasta")
        ),
        #assembly_tmp_dir=temp(directory(analysis_dir + "metaspades/{sample_short}/tmp")),
        assembly_dir=(
            directory(analysis_dir + "metaspades/{sample_short}")
            if config["keep_intermediate_files"] == True
            else temp(directory(analysis_dir + "metaspades/{sample_short}"))
        ),
    shell:
        """
        spades.py --s1 {input.reads} -o {output.assembly_dir} -t {threads} -m {resources.mem_mb} > {log.out} 2> {log.err}
        """


ruleorder: metaspades_paired_end > metaspades_single_end


is_nanopore = config["long_read_type"]["Nanopore"]
is_pacbio = (
    config["long_read_type"]["PacBio CLR"] or config["long_read_type"]["PacBio HiFi"]
)
is_metaflye = config["assembly_mode"]["metaflye"]
is_hybrid_metaspades = config["assembly_mode"]["hybrid_metaspades"]


def get_metaflye_read_type(is_nanopore, is_pacbio, is_metaflye):
    if is_nanopore == True and is_pacbio == False:
        setting = "--nano-raw"
        return setting
    elif is_nanopore == False and is_pacbio == True:
        setting = "--pacbio-raw"
        return setting


rule metaflye:
    resources:
        mem_mb=resources["metaflye"]["mem_mb"],
        runtime=resources["metaflye"]["runtime"],
        cpus_per_task=resources["metaflye"]["cpus_per_task"],
    threads: int(resources["metaflye"]["threads"])
    log:
        out=analysis_dir + "logs/metaflye/{sample_long}_out.log",
        err=analysis_dir + "logs/metaflye/{sample_long}_err.log",
    conda:
        snakefile_dir + "/envs/metaflye.yaml"
    input:
        reads=metaflye_input_dir + "/{sample_long}." + metaflye_extension,
    output:
        assembly=(
            analysis_dir + "metaflye/{sample_long}/assembly.fasta"
            if config["keep_intermediate_files"] == True
            else temp(analysis_dir + "metaflye/{sample_long}/assembly.fasta")
        ),
        assembly_dir=(
            directory(analysis_dir + "metaflye/{sample_long}")
            if config["keep_intermediate_files"] == True
            else temp(directory(analysis_dir + "metaflye/{sample_long}"))
        ),
    params:
        read_type=get_metaflye_read_type(is_nanopore, is_pacbio, is_metaflye),
    shell:
        """
        flye {params.read_type} {input.reads} -o {output.assembly_dir} --threads {threads} --meta > {log.out} 2> {log.err}
        """


def get_hybrid_metaspades_setting(is_nanopore, is_pacbio, is_hybrid_metaspades):
    if is_nanopore == True and is_pacbio == False:
        setting = "--nanopore"
        return setting
    elif is_nanopore == False and is_pacbio == True:
        setting = "--pacbio"
        return setting


rule hybrid_metaspades:
    resources:
        mem_mb=resources["hybrid_metaspades"]["mem_mb"],
        runtime=resources["hybrid_metaspades"]["runtime"],
        cpus_per_task=resources["hybrid_metaspades"]["cpus_per_task"],
    threads: int(resources["hybrid_metaspades"]["threads"])
    log:
        out=analysis_dir + "logs/hybrid_metaspades/{sample_short}-{sample_long}_out.log",
        err=analysis_dir + "logs/hybrid_metaspades/{sample_short}-{sample_long}_err.log",
    conda:
        snakefile_dir + "/envs/spades.yaml"
    input:
        long_reads=metaflye_input_dir + "/{sample_long}." + metaflye_extension,
        forward_paired=get_forward_paired,
        reverse_paired=get_reverse_paired,
    output:
        assembly=(
            analysis_dir
            + "hybrid_metaspades/{sample_short}-{sample_long}/scaffolds.fasta"
            if config["keep_intermediate_files"] == True
            else temp(
                analysis_dir
                + "hybrid_metaspades/{sample_short}-{sample_long}/scaffolds.fasta"
            )
        ),
        assembly_tmp_dir=temp(
            directory(
                analysis_dir + "hybrid_metaspades/{sample_short}-{sample_long}/temp"
            )
        ),
        assembly_dir=(
            directory(analysis_dir + "hybrid_metaspades/{sample_short}-{sample_long}")
            if config["keep_intermediate_files"] == True
            else temp(
                directory(
                    analysis_dir + "hybrid_metaspades/{sample_short}-{sample_long}"
                )
            )
        ),
    params:
        setting=get_hybrid_metaspades_setting(
            is_nanopore, is_pacbio, is_hybrid_metaspades
        ),
        forward_unpaired=get_forward_unpaired(sample_name_short),
        reverse_unpaired=get_reverse_unpaired(),
    shell:
        """
        if [[ -n "{params.forward_unpaired}" ]]; then
            echo "Unpaired reads are provided"
            forward_unpaired={short_read_input_dir}/{wildcards.sample_short}{params.forward_unpaired}
            echo $forward_unpaired
            reverse_unpaired={short_read_input_dir}/{wildcards.sample_short}{params.reverse_unpaired}
            if [[ -f $reverse_unpaired ]]; then
                echo $reverse_unpaired
                cat $forward_unpaired $reverse_unpaired > {analysis_dir}hybrid_metaspades/{wildcards.sample_short}-{wildcards.sample_long}/unpaired_merged.{short_read_extension}
            else
                cat $forward_unpaired > {analysis_dir}hybrid_metaspades/{wildcards.sample_short}-{wildcards.sample_long}/unpaired_merged.{short_read_extension} || true
            fi
            unpaired_merged={analysis_dir}hybrid_metaspades/{wildcards.sample_short}-{wildcards.sample_long}/unpaired_merged.{short_read_extension}
            echo $unpaired_merged
            spades.py --meta -1 {input.forward_paired} -2 {input.reverse_paired} -s ${{unpaired_merged}}  \
            {params.setting} {input.long_reads} -o {output.assembly_tmp_dir} -t {threads} -m {resources.mem_mb} > {log.out} 2> {log.err}
        else
            echo "No unpaired reads are provided"
            spades.py --meta -1 {input.forward_paired} -2 {input.reverse_paired} {params.setting} {input.long_reads} \
            -o {output.assembly_tmp_dir} -t {threads} -m {resources.mem_mb} > {log.out} 2> {log.err}
        fi
        mv -f {output.assembly_tmp_dir}/* {output.assembly_dir} >> {log.out} 2>> {log.err}
        """
