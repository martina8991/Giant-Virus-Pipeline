input_short = ""
sample_name_short = ""

if (
    config["assembly_mode"]["megahit"] == True
    or config["assembly_mode"]["metaspades"] == True
    or config["assembly_mode"]["hybrid_metaspades"] == True
):
    input_short = config["input_dir_short_reads"]
    if input_short == None:
        raise ValueError("input_dir_short_reads is not set in the config file")
    # input_dir = Path(config["input_dir"])
    # print("input_short: ", input_short)

    (sample_name_short,) = glob_wildcards(str(input_short) + "/{sample}")
    # print("sample_name_short: ", sample_name_short)

    if sample_name_short == []:
        raise ValueError("No samples found in the input directory for short reads")

input_long = ""
sample_name_long = ""

if (
    config["assembly_mode"]["metaflye"] == True
    or config["assembly_mode"]["hybrid_metaspades"] == True
):
    input_long = config["input_dir_long_reads"]
    if input_long == None:
        raise ValueError("input_dir_long_reads is not set in the config file")
    # print("input_long: ", input_long)

    (sample_name_long,) = glob_wildcards(str(input_long) + "/{sample}")
    # print("sample_name_long: ", sample_name_long)
    if sample_name_long == []:
        raise ValueError("No samples found in the input directory for long reads")


if config["analysis_folder_name"] == None:
    raise ValueError("analysis_folder_name is not set in the config file")
analysis_dir = config["analysis_folder_name"] + "/"


phred_cut_off_short = config["phred_minimum_quality_short_reads"]
if phred_cut_off_short == None:
    phred_cut_off_short = 28
elif isinstance(phred_cut_off_short, int) == False:
    raise ValueError("Invalid phred cut-off for short reads.")

phred_cut_off_long = config["phred_minimum_quality_long_reads"]
if phred_cut_off_long == None:
    phred_cut_off_long = 8
elif isinstance(phred_cut_off_long, int) == False:
    raise ValueError("Invalid phred cut-off for long reads.")


if config["paired_read_ids"]["forward"] == None:
    sample_forward_id = "_dummy1"
if config["paired_read_ids"]["reverse"] == None:
    sample_reverse_id = "_dummy2"
else:
    sample_forward_id = "_" + str(config["paired_read_ids"]["forward"])
    sample_reverse_id = "_" + str(config["paired_read_ids"]["reverse"])

short_read_paired_ids = [sample_forward_id, sample_reverse_id]
# print("sample_forward_id: ", sample_forward_id)
# print("sample_reverse_id: ", sample_reverse_id)
# print("short_read_paired_ids: ", short_read_paired_ids)


paired_name = str(config["paired_reads"]["paired_name"])
unpaired_name = str(config["paired_reads"]["unpaired_name"])


# print("paired_name: ", paired_name)
# print("unpaired_name: ", unpaired_name)


def get_paired_end_reads(sample_name):
    sample_paired = set()
    id_list = short_read_paired_ids
    for name in sample_name:  # can loop through "dummy" name
        for id in id_list:
            if id in name:
                name = name.split(".")[0]  # remove extension
                name = name.split("_", 1)[0]  # remove id
                # sample_1_paired, sample_2_unpaired, sample_1, sample_2, sample_1_P, sample_2_U  etc.
                sample_paired.add(name)
    sample_paired = list(sample_paired)
    return sample_paired


def get_single_end_reads(sample_name):
    sample_single = set()
    id_list = short_read_paired_ids
    for name in sample_name:  # can loop through "dummy" name
        id_match = [id for id in id_list if id in name]
        if id_match != [] or unpaired_name in name:
            continue
        else:
            name = name.split(".")[0]
            sample_single.add(name)
    sample_single = list(sample_single)
    return sample_single


def get_short_reads(sample_name):
    short_reads = set()
    for name in sample_name:  # can loop through "dummy" name
        name = name.split(".")[0]
        short_reads.add(name)
    short_reads = list(short_reads)
    return short_reads


def get_extension(sample_name):
    sample_extension_set = set()
    extension_list = [
        "fastq",
        "fastq.gz",
        "fasta",
        "fasta.gz",
        "fq",
        "fq.gz",
        "fa",
        "fa.gz",
    ]
    for name in sample_name:
        sample_ext = name.split(".", 1)[1:]
        sample_ext = "".join(sample_ext)
        for ext in extension_list:
            if ext == sample_ext:
                sample_extension_set.add(ext)
                break
    if len(sample_extension_set) == 1:
        sample_extension_list = list(sample_extension_set)
        sample_extension = sample_extension_list[0]
    elif len(sample_extension_set) > 1:
        raise ValueError("Multiple extensions found")
    elif sample_extension_set == set():
        raise ValueError("No extension found")
    return sample_extension


def get_long_reads(sample_name):
    sample_long = set()
    for name in sample_name:  # can loop through "dummy" name
        name = name.split(".")[0]
        sample_long.add(name)
    sample_long = list(sample_long)
    return sample_long


paired_end_base = get_paired_end_reads(sample_name_short)

# print("paired_end_base: ", paired_end_base)

sample_extension_short = "dummy"
sample_extension_long = "dummy"
if (
    config["assembly_mode"]["megahit"] == True
    or config["assembly_mode"]["metaspades"] == True
    or config["assembly_mode"]["hybrid_metaspades"] == True
):
    sample_extension_short = get_extension(sample_name_short)
    # print("sample_extension_short: ", sample_extension_short)
    if config["qc_step"] == True and sample_extension_short in [
        "fasta",
        "fa",
        "fasta.gz",
        "fa.gz",
    ]:
        raise ValueError("QC and read trimming cannot be done with fasta format.")

if (
    config["assembly_mode"]["metaflye"] == True
    or config["assembly_mode"]["hybrid_metaspades"] == True
):
    sample_extension_long = get_extension(sample_name_long)
    # print("sample_extension_long: ", sample_extension_long)
    if config["qc_step"] == True and sample_extension_long in [
        "fasta",
        "fa",
        "fasta.gz",
        "fa.gz",
    ]:
        raise ValueError("QC and read trimming cannot be done with fasta format.")


sample_base_long_reads = get_long_reads(sample_name_long)
# print("sample_base_long_reads: ", sample_base_long_reads)


single_end_base = get_single_end_reads(sample_name_short)

# print("single_end_base: ", single_end_base)

base_single_paired_samples = paired_end_base + single_end_base
# print("base_single_paired_samples: ", base_single_paired_samples)

base_with_id_single_paired = get_short_reads(sample_name_short)
# print("base_with_id_single_paired: ", base_with_id_single_paired)


def check_forward_unpaired(sample_name):
    for name in sample_name:
        if unpaired_name not in name:
            continue
        if unpaired_name in name and (
            sample_forward_id in name or sample_reverse_id in name
        ):
            return sample_forward_id + "_" + unpaired_name + "." + short_read_extension
        if unpaired_name in name and (
            sample_forward_id not in name and sample_reverse_id not in name
        ):
            return "_" + unpaired_name + "." + short_read_extension


def get_forward_unpaired(sample_name_short):
    if config["qc_step"] == True:
        forward_unpaired = (sample_forward_id + "_U." + short_read_extension,)
        return forward_unpaired
    if config["qc_step"] == False:
        if unpaired_name == "None":
            return ""
        else:
            forward_unpaired = check_forward_unpaired(sample_name_short)
            return forward_unpaired


def get_reverse_unpaired():
    if config["qc_step"] == True:
        reverse_unpaired = (sample_reverse_id + "_U." + short_read_extension,)
        return reverse_unpaired
    if config["qc_step"] == False:
        reverse_unpaired = (
            sample_reverse_id + "_" + unpaired_name + "." + short_read_extension
        )
        return reverse_unpaired


def get_forward_paired(wildcards):
    if config["qc_step"] == True:
        forward_paired = (
            short_read_input_dir
            + "/{sample_short}"
            + sample_forward_id
            + "_P."
            + short_read_extension,
        )
        return forward_paired
    if config["qc_step"] == False:
        if paired_name == "None":
            forward_paired = (
                short_read_input_dir
                + "/{sample_short}"
                + sample_forward_id
                + "."
                + short_read_extension
            )
            return forward_paired
        else:
            forward_paired = (
                short_read_input_dir
                + "/{sample_short}"
                + sample_forward_id
                + "_"
                + paired_name
                + "."
                + short_read_extension
            )
            return forward_paired


def get_reverse_paired(wildcards):
    if config["qc_step"] == True:
        reverse_paired = (
            short_read_input_dir
            + "/{sample_short}"
            + sample_reverse_id
            + "_P."
            + short_read_extension,
        )
        return reverse_paired
    if config["qc_step"] == False:
        if paired_name == "None":
            reverse_paired = (
                short_read_input_dir
                + "/{sample_short}"
                + sample_reverse_id
                + "."
                + short_read_extension
            )
            return reverse_paired
        else:
            reverse_paired = (
                short_read_input_dir
                + "/{sample_short}"
                + sample_reverse_id
                + "_"
                + paired_name
                + "."
                + short_read_extension
            )
            return reverse_paired


short_read_input_dir = ""
short_read_extension = ""
metaflye_input_dir = ""
metaflye_extension = ""

if config["qc_step"] == True:
    if (
        config["assembly_mode"]["megahit"] == True
        or config["assembly_mode"]["metaspades"] == True
        or config["assembly_mode"]["hybrid_metaspades"] == True
    ):
        short_read_input_dir = analysis_dir + "fastp"
        short_read_extension = "fastq.gz"
        # print("short_read_input_dir", short_read_input_dir)
        # print("short_read_extension", short_read_extension)
    if (
        config["assembly_mode"]["metaflye"] == True
        or config["assembly_mode"]["hybrid_metaspades"] == True
    ):
        metaflye_input_dir = analysis_dir + "chopper"
        metaflye_extension = "fastq.gz"
        # print("metaflye_input_dir", metaflye_input_dir)
        # print("metaflye_extension", metaflye_extension)
if config["qc_step"] == False:
    if (
        config["assembly_mode"]["megahit"] == True
        or config["assembly_mode"]["metaspades"] == True
        or config["assembly_mode"]["hybrid_metaspades"] == True
    ):
        short_read_input_dir = input_short
        short_read_extension = sample_extension_short
        # print("short_read_input_dir", short_read_input_dir)
        # print("short_read_extension", short_read_extension)
    if (
        config["assembly_mode"]["metaflye"] == True
        or config["assembly_mode"]["hybrid_metaspades"] == True
    ):
        metaflye_input_dir = input_long
        metaflye_extension = sample_extension_long
        # print("metaflye_input_dir", metaflye_input_dir)
        # print("metaflye_extension", metaflye_extension)

qc_short_reads = []


if config["qc_step"] == True and (
    config["assembly_mode"]["megahit"] == True
    or config["assembly_mode"]["metaspades"] == True
    or config["assembly_mode"]["hybrid_metaspades"] == True
):
    # print("QC step is set to True")
    fastqc_raw = expand(
        analysis_dir + "read_qc/raw/{sample_short}_fastqc.{extension}",
        sample_short=base_with_id_single_paired,
        extension=["html", "zip"],
    )
    qc_short_reads.extend(fastqc_raw)
    fastqc_trimmed_SE = expand(
        analysis_dir + "read_qc/trimmed/{sample_short}_fastqc.{extension}",
        sample_short=single_end_base,
        extension=["html", "zip"],
    )
    qc_short_reads.extend(fastqc_trimmed_SE)
    fastqc_trimmed_PE = expand(
        analysis_dir + "read_qc/trimmed/{sample_short}{id}_{ftype}_fastqc.{extension}",
        sample_short=paired_end_base,
        id=[sample_forward_id, sample_reverse_id],
        ftype=["P", "U"],
        extension=["html", "zip"],
    )
    qc_short_reads.extend(fastqc_trimmed_PE)


# print("qc_short_reads: ", qc_short_reads)

if config["assembly_mode"]["metaspades"] == True:
    short_read_assembler = "metaspades"
    short_read_assembler_file = "scaffolds.fasta"
if config["assembly_mode"]["megahit"] == True:
    short_read_assembler = "megahit"
    short_read_assembler_file = "final.contigs.fa"


qc_long_reads = []

if config["qc_step"] == True and (
    config["assembly_mode"]["metaflye"] == True
    or config["assembly_mode"]["hybrid_metaspades"] == True
):
    # print("QC step is set to True")
    fastqc = expand(
        analysis_dir + "read_qc/{dir}/{sample_long}_fastqc.{extension}",
        dir=["raw", "trimmed"],
        sample_long=sample_base_long_reads,
        extension=["html", "zip"],
    )
    qc_long_reads.extend(fastqc)

    nanoplot_report = expand(
        analysis_dir
        + "nanoplot/{dir}/{sample_long}/{sample_long}_NanoPlot-report.html",
        dir=["raw", "trimmed"],
        sample_long=sample_base_long_reads,
    )
    qc_long_reads.extend(nanoplot_report)

# print("qc_long_reads: ", qc_long_reads)


short_reads_1_for_hybrid = config["hybrid_assembly_pairs"]["short_reads_forward"]
short_reads_2_for_hybrid = config["hybrid_assembly_pairs"]["short_reads_reverse"]
long_reads_for_hybrid = config["hybrid_assembly_pairs"]["long_reads"]

if config["assembly_mode"]["hybrid_metaspades"] == True and (
    short_reads_1_for_hybrid == []
    or short_reads_2_for_hybrid == []
    or long_reads_for_hybrid == []
):
    raise ValueError(
        "Hybrid assembly is set to True, however no input files are provided in the config file"
    )

if config["assembly_mode"]["hybrid_metaspades"] == True and (
    len(short_reads_1_for_hybrid)
    != len(short_reads_2_for_hybrid)
    != len(long_reads_for_hybrid)
):
    raise ValueError("Number of input files for hybrid assembly is not equal")

short_reads_hybrid_base = [name.split("_")[0] for name in short_reads_1_for_hybrid]
# print("Short reads 1 for hybrid: ", short_reads_1_for_hybrid)
# print("Short reads 2 for hybrid: ", short_reads_2_for_hybrid)
# print("Long reads for hybrid: ", long_reads_for_hybrid)
# print("Short reads hybrid base: ", short_reads_hybrid_base)


hybrid_metaspades_base = []

if config["assembly_mode"]["hybrid_metaspades"] == True:
    hybrid_metaspades = expand(
        "{sample_short}-{sample_long}",
        zip,
        sample_short=short_reads_hybrid_base,
        sample_long=long_reads_for_hybrid,
    )
    hybrid_metaspades_base.extend(hybrid_metaspades)


# print("hybrid_metaspades_base: ", hybrid_metaspades_base)


combined_bases = []
mapping_files = []

if (
    config["assembly_mode"]["megahit"] == True
    or config["assembly_mode"]["metaspades"] == True
):
    combined_bases.extend(base_single_paired_samples)
    mapping_files_short = ["bowtie2/" + file for file in base_single_paired_samples]
    mapping_files.extend(mapping_files_short)


if config["assembly_mode"]["metaflye"] == True:
    combined_bases.extend(sample_base_long_reads)
    mapping_files_long = ["minimap2/" + file for file in sample_base_long_reads]
    mapping_files.extend(mapping_files_long)


if config["assembly_mode"]["hybrid_metaspades"] == True:
    combined_bases.extend(hybrid_metaspades_base)
    mapping_files_short_long = expand(
        "{mapping_dir}/{file}",
        mapping_dir=["bowtie2", "minimap2"],
        file=hybrid_metaspades_base,
    )
    mapping_files.extend(mapping_files_short_long)


# print("combined_bases: ", combined_bases)
# print("mapping_files: ", mapping_files)

if (
    config["assembly_mode"]["metaflye"] == True
    or config["assembly_mode"]["hybrid_metaspades"] == True
):
    if (
        config["long_read_type"]["Nanopore"] == False
        and config["long_read_type"]["PacBio CLR"] == False
        and config["long_read_type"]["PacBio HiFi"] == False
    ):
        raise ValueError("No long read type is set to True in the config file")

if config["long_read_type"]["Nanopore"] == True and (
    config["long_read_type"]["PacBio CLR"] == True
    or config["long_read_type"]["PacBio HiFi"] == True
):
    raise ValueError("Both Nanopore and PacBio cannot be used at the same time")

if (
    config["long_read_type"]["PacBio CLR"] == True
    and config["long_read_type"]["PacBio HiFi"] == True
):
    raise ValueError("Both PacBio CLR and PacBio HiFi cannot be used at the same time")

if config["barrnap_evalue"] == None:
    raise ValueError("barrnap_evalue is not set in the config file")


num_of_splits = config["num_of_file_splits"]

if num_of_splits == None:
    num_of_splits = 1
if num_of_splits == 0:
    num_of_splits = 1
if num_of_splits > 40:
    raise ValueError("num_of_splits should be less than 40")


if config["hmmsearch_GVDB_evalue"] == None:
    raise ValueError("hmmsearch_GVDB_evalue is not set in the config file")
if config["hmmsearch_VOGDB_evalue"] == None:
    raise ValueError("hmmsearch_VOGDB_evalue is not set in the config file")

if config["min_length_contigs_for_gv_table"] == None:
    raise ValueError("min_length_contigs_for_gv_table is not set in the config file")
if config["min_length_contigs_for_gv_table"] < 1000:
    raise ValueError("min_length_contigs_for_gv_table should be greater than 1000")

if config["keep_intermediate_files"] == None:
    raise ValueError("keep_intermediate_files is not set in the config file")

if (
    config["assembly_mode"]["metaspades"] == True
    and config["assembly_mode"]["megahit"] == True
):
    raise ValueError(
        "cannot assemble short reads with both Megahit and MetaSPAdes assembler simultaneously"
    )
