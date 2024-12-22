

rule barrnap:
    resources:
        mem_mb=resources["barrnap"]["mem_mb"],
        runtime=resources["barrnap"]["runtime"],
        cpus_per_task=resources["barrnap"]["cpus_per_task"],
    threads: int(resources["barrnap"]["threads"])
    log:
        err=analysis_dir + "logs/barrnap/{sample}_err.log",
    conda:
        snakefile_dir + "/envs/gene_prediction.yaml"
    input:
        assembly=analysis_dir + "filtered_assemblies/{sample}.fa",
    output:
        bac_gff=(
            analysis_dir + "barrnap/{sample}/{sample}_bac.gff"
            if config["keep_intermediate_files"] == True
            else temp(analysis_dir + "barrnap/{sample}/{sample}_bac.gff")
        ),
        euk_gff=(
            analysis_dir + "barrnap/{sample}/{sample}_euk.gff"
            if config["keep_intermediate_files"] == True
            else temp(analysis_dir + "barrnap/{sample}/{sample}_euk.gff")
        ),
        arc_gff=(
            analysis_dir + "barrnap/{sample}/{sample}_arc.gff"
            if config["keep_intermediate_files"] == True
            else temp(analysis_dir + "barrnap/{sample}/{sample}_arc.gff")
        ),
        mito_gff=(
            analysis_dir + "barrnap/{sample}/{sample}_mito.gff"
            if config["keep_intermediate_files"] == True
            else temp(analysis_dir + "barrnap/{sample}/{sample}_mito.gff")
        ),
    shell:
        """
        barrnap --kingdom bac --threads {threads} {input.assembly} > {output.bac_gff} 2>> {log.err}
        barrnap --kingdom euk --threads {threads} {input.assembly} > {output.euk_gff} 2>> {log.err}
        barrnap --kingdom arc --threads {threads} {input.assembly} > {output.arc_gff} 2>> {log.err}
        barrnap --kingdom mito --threads {threads} {input.assembly} > {output.mito_gff} 2>> {log.err}
        """


rule remove_rRNA_contigs:
    log:
        out=analysis_dir + "logs/contigs_without_rRNA/{sample}_out.log",
        err=analysis_dir + "logs/contigs_without_rRNA/{sample}_err.log",
    conda:
        snakefile_dir + "/envs/assembly_qc.yaml"
    input:
        assembly=analysis_dir + "filtered_assemblies/{sample}.fa",
        barrnap_gff=expand(
            analysis_dir + "barrnap/{{sample}}/{{sample}}_{species}.gff",
            species=["bac", "euk", "arc", "mito"],
        ),
    output:
        all_contig_ids=temp(
            analysis_dir + "contigs_without_rRNA/{sample}_all_contig_ids.txt"
        ),
        non_rna_contig_ids=temp(
            analysis_dir + "contigs_without_rRNA/{sample}_non_rRNA_contig_ids.txt"
        ),
        merged_gff=temp(analysis_dir + "contigs_without_rRNA/{sample}_merged.gff"),
        merged_sorted_gff=temp(
            analysis_dir + "contigs_without_rRNA/{sample}_merged_sorted.gff"
        ),
        non_rRNA_contigs=(
            analysis_dir + "contigs_without_rRNA/{sample}_non_rRNA.fa"
            if config["keep_intermediate_files"] == True
            else temp(analysis_dir + "contigs_without_rRNA/{sample}_non_rRNA.fa")
        ),
    params:
        parse_script=snakefile_dir + "/scripts/parse_barrnap.py",
        e_value=config["barrnap_evalue"],
    shell:
        """
        cat {input.barrnap_gff} > {output.merged_gff} 2> {log.err}
        grep -v "^#" {output.merged_gff} | sort > {output.merged_sorted_gff} 2>> {log.err}
        grep ">" {input.assembly} | cut -f 2 -d '>' | cut -f 1 -d ' ' > {output.all_contig_ids} 2>> {log.err}
        python {params.parse_script} {output.merged_sorted_gff} {params.e_value} {output.all_contig_ids} {output.non_rna_contig_ids} > {log.out} 2>> {log.err}
        seqtk subseq {input.assembly} {output.non_rna_contig_ids} > {output.non_rRNA_contigs} 2>> {log.err}
        """


def get_num_range(num_of_splits):
    digit_count = len(str(num_of_splits))
    num_range = range(num_of_splits)
    num_range_correct_digit = [str(num).zfill(digit_count) for num in num_range]
    return num_range_correct_digit


num_range_splits = get_num_range(num_of_splits)
# print("num_range_splits: ", num_range_splits)


rule split_contigs:
    log:
        out=analysis_dir + "logs/split_contigs/{sample}_out.log",
        err=analysis_dir + "logs/split_contigs/{sample}_err.log",
    conda:
        snakefile_dir + "/envs/ucsc-fasplit.yaml"
    input:
        assembly=analysis_dir + "contigs_without_rRNA/{sample}_non_rRNA.fa",
    output:
        split_files=temp(
            expand(
                analysis_dir
                + "split_contigs/{{sample}}_split/{{sample}}_split_{id}.fa",
                id=num_range_splits,
            )
        ),
    params:
        splits=num_of_splits,
        out_dir=analysis_dir + "split_contigs/{sample}_split",
    shell:
        """
        mkdir -p {params.out_dir} > {log.out} 2> {log.err}
        faSplit sequence {input.assembly} {params.splits} {params.out_dir}/{wildcards.sample}_split_ >> {log.out} 2>> {log.err}
        """


rule prodigal:
    resources:
        mem_mb=resources["prodigal"]["mem_mb"],
        runtime=resources["prodigal"]["runtime"],
    log:
        out=analysis_dir + "logs/prodigal/{sample}/{sample}_split_{id}_out.log",
        err=analysis_dir + "logs/prodigal/{sample}/{sample}_split_{id}_err.log",
    conda:
        snakefile_dir + "/envs/gene_prediction.yaml"
    input:
        file=analysis_dir + "split_contigs/{sample}_split/{sample}_split_{id}.fa",
    output:
        gff_file=temp(analysis_dir + "prodigal/{sample}/{sample}_split_{id}.gff"),
        faa_file=temp(analysis_dir + "prodigal/{sample}/{sample}_split_{id}.faa"),
        fna_file=temp(analysis_dir + "prodigal/{sample}/{sample}_split_{id}.fna"),
    shell:
        """
        prodigal -c -f gff -i {input.file} -o {output.gff_file} -a {output.faa_file} -d {output.fna_file} -p meta > {log.out} 2> {log.err}
        """


def check_if_db_pressed(name, path):
    file_names = [name + ".h3f", name + ".h3i", name + ".h3m", name + ".h3p"]
    for root, dirs, files in os.walk(path):
        for file in file_names:
            file_path = os.path.join(root, file)
            # print("file_path_name", file_path)
            if os.path.exists(file_path):
                continue
            else:
                return "hmmpress " + path + "/" + name + " > {log.out} 2> {log.err}"
    return ""


rule hmm_press:
    log:
        out=analysis_dir + "logs/hmmpress/hmmpress_out.log",
        err=analysis_dir + "logs/hmmpress/hmmpress_err.log",
    conda:
        snakefile_dir + "/envs/gene_prediction.yaml"
    input:
        GVDB=snakefile_dir + "/databases/GVDB/GVOGs/gvog.complete.hmm",
        VOGDB=snakefile_dir + "/databases/VOGDB/hmm_merged/merged.hmm",
    output:
        GVDB=expand(
            snakefile_dir + "/databases/GVDB/GVOGs/gvog.complete.hmm.{extension}",
            extension=["h3f", "h3i", "h3m", "h3p"],
        ),
        VOGDB=expand(
            snakefile_dir + "/databases/VOGDB/hmm_merged/merged.hmm.{extension}",
            extension=["h3f", "h3i", "h3m", "h3p"],
        ),
    shell:
        """
        hmmpress {input.GVDB} > {log.out} 2> {log.err}
        hmmpress {input.VOGDB} >> {log.out} 2>> {log.err}
        """


rule hmmersearch_GVDB:
    resources:
        mem_mb=resources["hmmsearch_GVDB"]["mem_mb"],
        runtime=resources["hmmsearch_GVDB"]["runtime"],
        cpus_per_task=resources["hmmsearch_GVDB"]["cpus_per_task"],
    threads: int(resources["hmmsearch_GVDB"]["threads"])
    retries: 3
    log:
        out=analysis_dir + "logs/hmmsearch_GVDB/{sample}/{sample}_split_{id}_out.log",
        err=analysis_dir + "logs/hmmsearch_GVDB/{sample}/{sample}_split_{id}_err.log",
    conda:
        snakefile_dir + "/envs/gene_prediction.yaml"
    input:
        prodigal=analysis_dir + "prodigal/{sample}/{sample}_split_{id}.faa",
        pressed_GVDB=expand(
            snakefile_dir + "/databases/GVDB/GVOGs/gvog.complete.hmm.{extension}",
            extension=["h3f", "h3i", "h3m", "h3p"],
        ),
    output:
        out=temp(analysis_dir + "hmmsearch/hmmsearch_GVDB/{sample}/out_split_{id}.txt"),
        per_sequence_hits=temp(
            analysis_dir
            + "hmmsearch/hmmsearch_GVDB/{sample}/per_sequence_hits_split_{id}.txt"
        ),
        per_dom_hits=temp(
            analysis_dir
            + "hmmsearch/hmmsearch_GVDB/{sample}/per_dom_hits_split_{id}.txt"
        ),
    params:
        database=snakefile_dir + "/databases/GVDB/GVOGs/gvog.complete.hmm",
    shell:
        """
        hmmsearch --cpu {threads} --tblout {output.per_sequence_hits} --domtblout {output.per_dom_hits} \
        -o {output.out} {params.database} {input.prodigal} >> {log.out} 2>> {log.err}
        """


rule hmmersearch_VOGDB:
    resources:
        mem_mb=resources["hmmsearch_VOGDB"]["mem_mb"],
        runtime=resources["hmmsearch_VOGDB"]["runtime"],
        cpus_per_task=resources["hmmsearch_VOGDB"]["cpus_per_task"],
    threads: int(resources["hmmsearch_VOGDB"]["threads"])
    retries: 3
    log:
        out=analysis_dir + "logs/hmmsearch_VOGDB/{sample}/{sample}_split_{id}_out.log",
        err=analysis_dir + "logs/hmmsearch_VOGDB/{sample}/{sample}_split_{id}_err.log",
    conda:
        snakefile_dir + "/envs/gene_prediction.yaml"
    input:
        prodigal=analysis_dir + "prodigal/{sample}/{sample}_split_{id}.faa",
        pressed_VOGDB=expand(
            snakefile_dir + "/databases/VOGDB/hmm_merged/merged.hmm.{extension}",
            extension=["h3f", "h3i", "h3m", "h3p"],
        ),
    output:
        out=temp(analysis_dir + "hmmsearch/hmmsearch_VOGDB/{sample}/out_split_{id}.txt"),
        per_sequence_hits=temp(
            analysis_dir
            + "hmmsearch/hmmsearch_VOGDB/{sample}/per_sequence_hits_split_{id}.txt"
        ),
        per_dom_hits=temp(
            analysis_dir
            + "hmmsearch/hmmsearch_VOGDB/{sample}/per_dom_hits_split_{id}.txt"
        ),
    params:
        database=snakefile_dir + "/databases/VOGDB/hmm_merged/merged.hmm",
    shell:
        """
        hmmsearch --cpu {threads} --tblout {output.per_sequence_hits} --domtblout {output.per_dom_hits} \
        -o {output.out} {params.database} {input.prodigal} >> {log.out} 2>> {log.err}
        """


header_string = """#                                                               --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------"""


rule remerge_files:
    log:
        err=analysis_dir + "logs/remerge_files/{sample}/{sample}_err.log",
    input:
        split_files_vogdb_tblout=expand(
            analysis_dir
            + "hmmsearch/hmmsearch_VOGDB/{{sample}}/per_sequence_hits_split_{id}.txt",
            id=num_range_splits,
        ),
        split_files_vogdb_domtblout=expand(
            analysis_dir
            + "hmmsearch/hmmsearch_VOGDB/{{sample}}/per_dom_hits_split_{id}.txt",
            id=num_range_splits,
        ),
        split_files_vogdb_out=expand(
            analysis_dir + "hmmsearch/hmmsearch_VOGDB/{{sample}}/out_split_{id}.txt",
            id=num_range_splits,
        ),
        split_files_gvdb_tblout=expand(
            analysis_dir
            + "hmmsearch/hmmsearch_GVDB/{{sample}}/per_sequence_hits_split_{id}.txt",
            id=num_range_splits,
        ),
        split_files_gvdb_domtblout=expand(
            analysis_dir
            + "hmmsearch/hmmsearch_GVDB/{{sample}}/per_dom_hits_split_{id}.txt",
            id=num_range_splits,
        ),
        split_files_gvdb_out=expand(
            analysis_dir + "hmmsearch/hmmsearch_GVDB/{{sample}}/out_split_{id}.txt",
            id=num_range_splits,
        ),
        split_files_prodigal_faa=expand(
            analysis_dir + "prodigal/{{sample}}/{{sample}}_split_{id}.faa",
            id=num_range_splits,
        ),
        split_files_prodigal_fna=expand(
            analysis_dir + "prodigal/{{sample}}/{{sample}}_split_{id}.fna",
            id=num_range_splits,
        ),
        split_files_prodigal_gff=expand(
            analysis_dir + "prodigal/{{sample}}/{{sample}}_split_{id}.gff",
            id=num_range_splits,
        ),
    output:
        remerged_vogdb_tblout=analysis_dir
        + "hmmsearch/hmmsearch_VOGDB/{sample}/per_sequence_hits.txt",
        remerged_vogdb_domtblout=analysis_dir
        + "hmmsearch/hmmsearch_VOGDB/{sample}/per_dom_hits.txt",
        remerged_vogdb_out=analysis_dir + "hmmsearch/hmmsearch_VOGDB/{sample}/out.txt",
        remerged_gvdb_tblout=analysis_dir
        + "hmmsearch/hmmsearch_GVDB/{sample}/per_sequence_hits.txt",
        remerged_gvdb_domtblout=analysis_dir
        + "hmmsearch/hmmsearch_GVDB/{sample}/per_dom_hits.txt",
        remerged_gvdb_out=analysis_dir + "hmmsearch/hmmsearch_GVDB/{sample}/out.txt",
        remerged_prodigal_faa=analysis_dir + "prodigal/{sample}/{sample}.faa",
        remerged_prodigal_fna=analysis_dir + "prodigal/{sample}/{sample}.fna",
        remerged_prodigal_gff=analysis_dir + "prodigal/{sample}/{sample}.gff",
        header_file=temp(analysis_dir + "hmmsearch/header_{sample}.txt"),
    params:
        header=header_string,
        analysis_dir=analysis_dir,
    shell:
        """
        echo "{params.header}" > {output.header_file} 2> {log.err}
        cat {input.split_files_vogdb_tblout} > {params.analysis_dir}/hmmsearch/hmmsearch_VOGDB/{wildcards.sample}/per_sequence_hits_merged.txt  2>> {log.err}

        grep -v "^#" {params.analysis_dir}/hmmsearch/hmmsearch_VOGDB/{wildcards.sample}/per_sequence_hits_merged.txt | sort > \
        {params.analysis_dir}/hmmsearch/hmmsearch_VOGDB/{wildcards.sample}/per_sequence_hits_merged_sorted_filtered.txt 2>> {log.err}

        cat {output.header_file} {params.analysis_dir}/hmmsearch/hmmsearch_VOGDB/{wildcards.sample}/per_sequence_hits_merged_sorted_filtered.txt > \
        {output.remerged_vogdb_tblout} 2>> {log.err}

        rm {params.analysis_dir}/hmmsearch/hmmsearch_VOGDB/{wildcards.sample}/per_sequence_hits_merged*.txt  2>> {log.err}

        cat {input.split_files_vogdb_domtblout} > {params.analysis_dir}/hmmsearch/hmmsearch_VOGDB/{wildcards.sample}/per_dom_hits_merged.txt 2>> {log.err}

        grep -v "^#" {params.analysis_dir}/hmmsearch/hmmsearch_VOGDB/{wildcards.sample}/per_dom_hits_merged.txt | sort > \
        {params.analysis_dir}/hmmsearch/hmmsearch_VOGDB/{wildcards.sample}/per_dom_hits_merged_sorted_filtered.txt 2>> {log.err}

        cat {output.header_file} {params.analysis_dir}/hmmsearch/hmmsearch_VOGDB/{wildcards.sample}/per_dom_hits_merged_sorted_filtered.txt > \
        {output.remerged_vogdb_domtblout} 2>> {log.err}
        rm {params.analysis_dir}/hmmsearch/hmmsearch_VOGDB/{wildcards.sample}/per_dom_hits_merged*.txt  2>> {log.err}

        cat {input.split_files_vogdb_out} > {output.remerged_vogdb_out} 2>> {log.err}

        cat {input.split_files_gvdb_tblout} > {params.analysis_dir}/hmmsearch/hmmsearch_GVDB/{wildcards.sample}/per_sequence_hits_merged.txt 2>> {log.err}

        grep -v "^#" {params.analysis_dir}/hmmsearch/hmmsearch_GVDB/{wildcards.sample}/per_sequence_hits_merged.txt | sort > \
        {params.analysis_dir}/hmmsearch/hmmsearch_GVDB/{wildcards.sample}/per_sequence_hits_merged_sorted_filtered.txt 2>> {log.err}

        cat {output.header_file} {params.analysis_dir}/hmmsearch/hmmsearch_GVDB/{wildcards.sample}/per_sequence_hits_merged_sorted_filtered.txt > \
        {output.remerged_gvdb_tblout} 2>> {log.err}

        rm {params.analysis_dir}/hmmsearch/hmmsearch_GVDB/{wildcards.sample}/per_sequence_hits_merged*.txt 2>> {log.err}

        cat {input.split_files_gvdb_domtblout} > {params.analysis_dir}/hmmsearch/hmmsearch_GVDB/{wildcards.sample}/per_dom_hits_merged.txt 2>> {log.err}

        grep -v "^#" {params.analysis_dir}/hmmsearch/hmmsearch_GVDB/{wildcards.sample}/per_dom_hits_merged.txt | sort > \
        {params.analysis_dir}/hmmsearch/hmmsearch_GVDB/{wildcards.sample}/per_dom_hits_merged_sorted_filtered.txt 2>> {log.err}

        cat {output.header_file} {params.analysis_dir}/hmmsearch/hmmsearch_GVDB/{wildcards.sample}/per_dom_hits_merged_sorted_filtered.txt > \
        {output.remerged_gvdb_domtblout} 2>> {log.err}
        rm {params.analysis_dir}/hmmsearch/hmmsearch_GVDB/{wildcards.sample}/per_dom_hits_merged*.txt  2>> {log.err}

        cat {input.split_files_gvdb_out} > {output.remerged_gvdb_out} 2>> {log.err}


        cat {input.split_files_prodigal_faa} > {output.remerged_prodigal_faa} 2>> {log.err}
        cat {input.split_files_prodigal_fna} > {output.remerged_prodigal_fna} 2>> {log.err}
        cat {input.split_files_prodigal_gff} > {output.remerged_prodigal_gff} 2>> {log.err}

        """
