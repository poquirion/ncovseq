from packaging import version
from snakemake.logging import logger
import sys
from shutil import which
import datetime
#
# Verify that required versions of dependencies are installed.
#
MIN_AUGUR_VERSION = "7.0.2"

try:
    from augur.__version__ import __version__ as augur_version
except ModuleNotFoundError:
    logger.error("ERROR: Could not find augur. Follow installation instructions at https://nextstrain.org/docs/ and try again.")
    sys.exit(1)

if version.parse(augur_version) < version.parse(MIN_AUGUR_VERSION):
    logger.error("ERROR: Found version '%s' of augur, but version '%s' or greater is required" % (augur_version, MIN_AUGUR_VERSION))
    sys.exit(1)

SHELL_COMMANDS_NEEDED = ["augur", "iqtree", "mafft"]
for sh_cmd in SHELL_COMMANDS_NEEDED:
    if not which(sh_cmd):
        logger.error(f"ERROR: `{sh_cmd}` is not available as a shell command. Please follow installation instructions at https://nextstrain.org/docs/ and try again.")
        sys.exit(1)


def get_todays_date():
    from datetime import datetime
    date = datetime.today().strftime('%Y-%m-%d')
    return date

# For information on how to run 'regions' runs, see Snakefile_Regions

# Add new regions here!
REGIONS = ["_quebec","_global"]

wildcard_constraints:
    region = "|".join(REGIONS) + "||",
    date = "[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]"

configfile: "config/Snakefile.yaml"


# simple rule to call snakemake for outsider users
rule all:
    input:
        auspice_json = "auspice/ncov.json",
        tip_frequencies_json = "auspice/ncov_tip-frequencies.json"



rule download_latest_excludes:
    message: "get metadata from nextstrain"
    output:
        local_world_exclude = 'results/nextstrain_exclude{}.tsv'.format(datetime.date.today())
    shell:
        """
        wget {config[exclude_latest]} -O {output.local_world_exclude}
        """

rule data_setup:
    message: "Organize and merge metadata, fasta, exclude"
    input:
        fasta_path = config["fasta_path"],
        qc_meta = config['qc_meta'],
        world_meta = config['world_meta'],
        world_exclude = rules.download_latest_excludes.output.local_world_exclude
    output:
        sequences = "results/sequences.fasta",
        metadata = "results/merged_metadata.tsv",
        exclude = "results/exclude.txt",
        selected_lspq = "results/lspq_only_metadata.tsv",
        ordering = "results/ordering.tsv",
        lat_long = "results/lat_long.tsv"
    params:
        extra_fasta = config['extra_fasta']
    shell:
        """
        scripts/concat_fasta.py --input_dir {input.fasta_path} --output {output.sequences}\
         {params.extra_fasta}
        scripts/concat_meta.py --inspq_meta {input.qc_meta} --nextstrain_metadata {input.world_meta}  --output \
        {output.metadata} --fasta_dir  {input.fasta_path} --output_lspq_only {output.selected_lspq} \
        --out_order {output.ordering} --out_lat_long {output.lat_long}
        cp {input.world_exclude} {output.exclude} 
        """



rule filter:
    message:
        """
        Filtering to
          - excluding strains in {input.exclude}
          - minimum genome length of {params.min_length}
        """
    input:
        sequences = rules.data_setup.output.sequences,
         metadata = rules.data_setup.output.metadata,
         include = config["include"],
         exclude = rules.data_setup.output.exclude
    output:
        sequences = "results/filtered.fasta"
    params:
        min_length = config["min_length"],
        exclude_where = config["exclude_where"],
        group_by = config["group_by"],
        sequences_per_group = config["sequences_per_group"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --include {input.include} \
            --exclude {input.exclude} \
            --exclude-where {params.exclude_where}\
            --min-length {params.min_length} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --output {output.sequences}

            echo contatenated fasta: `grep '>' {output.sequences} | wc -l` 
        """

checkpoint partition_sequences:
    input:
        sequences = rules.filter.output.sequences
    output:
        split_sequences = directory("results/split_sequences/pre/")
    params:
        sequences_per_group = 150
    shell:
        """
        python3 scripts/partition-sequences.py \
            --sequences {input.sequences} \
            --sequences-per-group {params.sequences_per_group} \
            --output-dir {output.split_sequences}
        """

rule partitions_intermediate:
    message:
        """
        partitions_intermediate: Copying sequence fastas
        {wildcards.cluster}
        """
    input:
        "results/split_sequences/pre/{cluster}.fasta"
    output:
        "results/split_sequences/post/{cluster}.fasta"
    shell:
        "cp {input} {output}"

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - gaps relative to reference are considered real
        Cluster:  {wildcards.cluster}
        """
    input:
        sequences = rules.partitions_intermediate.output,
        reference = config["reference"]
    output:
        alignment = "results/split_alignments/{cluster}.fasta"
    threads: 8
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --nthreads {threads} \
            --remove-reference \
            --fill-gaps
        """

def _get_alignments(wildcards):
    checkpoint_output = checkpoints.partition_sequences.get(**wildcards).output[0]
    return expand("results/split_alignments/{i}.fasta",
                  i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fasta")).i)

rule aggregate_alignments:
    message: "Collecting alignments"
    input:
        alignments = _get_alignments
    output:
        alignment = "results/aligned.fasta"
    shell:
        """
        cat {input.alignments} > {output.alignment}
        """

rule mask:
    message:
        """
        Mask bases in alignment
          - masking {params.mask_from_beginning} from beginning
          - masking {params.mask_from_end} from end
          - masking other sites: {params.mask_sites}
        """
    input:
        alignment = rules.aggregate_alignments.output.alignment
    output:
        alignment = "results/masked.fasta"
    params:
        mask_from_beginning = 130,
        mask_from_end = 50,
        mask_sites = "18529 29849 29851 29853"
    shell:
        """
        python3 scripts/mask-alignment.py \
            --alignment {input.alignment} \
            --mask-from-beginning {params.mask_from_beginning} \
            --mask-from-end {params.mask_from_end} \
            --mask-sites {params.mask_sites} \
            --output {output.alignment}
        """

rule subsample:
    input:
        rules.mask.output.alignment
    output:
        "results/subsampled_alignment{region}.fasta"
    shell:
        """
        cp {input} {output}
        """

rule adjust_metadata:
    input:
        rules.data_setup.output.metadata
    output:
        "results/metadata_adjusted{region}.tsv"
    shell:
        """
        cp {input} {output}
        """

rule tree:
    message: "Building tree"
    input:
        alignment = "results/subsampled_alignment{region}.fasta"
    output:
        tree = "results/tree_raw{region}.nwk"
    threads: 8
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads {threads} \
            --substitution-model {config[substitution-model]} \
	        --method {config[tree_method]}
        """

def _get_alignments_for_tree(wildcards):
    """Global builds use the complete alignment of sequences while regional builds
    use a subsampled set of sequences.
    """
    if wildcards.region == "":
        return rules.mask.output.alignment
    else:
        return rules.subsample_regions.output.alignment

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
        """
    input:
        tree = rules.tree.output.tree,
        alignment = _get_alignments_for_tree,
        metadata = rules.adjust_metadata.output
    output:
        tree = "results/tree{region}.nwk",
        node_data = "results/branch_lengths{region}.json"
    threads: 1
    params:
        root = "--root Wuhan-Hu-1/2019 Wuhan/WH01/2019",
        clock_rate = 0.0008,
        clock_std_dev = 0.0004,
        coalescent = "skyline",
        date_inference = "marginal",
        divergence_unit = "mutations",
        clock_filter_iqd = 4
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            {params.root} \
            --timetree \
            --clock-rate {params.clock_rate} \
            --clock-std-dev {params.clock_std_dev} \
            --coalescent {params.coalescent} \
            --date-inference {params.date_inference} \
            --divergence-unit {params.divergence_unit} \
            --date-confidence \
            --no-covariance \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    message:
        """
        Reconstructing ancestral sequences and mutations
          - inferring ambiguous mutations
        """
    input:
        tree = "results/tree{region}.nwk",
        alignment = _get_alignments_for_tree
    output:
        node_data = "results/nt_muts{region}.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference} \
            --infer-ambiguous
        """


rule haplotype_status:
    message: "Annotating haplotype status relative to {params.reference_node_name}"
    input:
        nt_muts = rules.ancestral.output.node_data
    output:
        node_data = "results/haplotype_status{region}.json"
    params:
        reference_node_name = "Wuhan/WH01/2019"
    shell:
        """
        python3 scripts/annotate-haplotype-status.py \
            --ancestral-sequences {input.nt_muts} \
            --reference-node-name {params.reference_node_name:q} \
            --output {output.node_data}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = "results/tree{region}.nwk",
        node_data = rules.ancestral.output.node_data,
        reference = config["reference"]
    output:
        node_data = "results/aa_muts{region}.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
        """


def _get_sampling_trait_for_wildcards(wildcards):
    mapping = {}
    return mapping[wildcards.region] if wildcards.region in mapping else "location"



def _get_exposure_trait_for_wildcards(wildcards):
    mapping = {}
    return mapping[wildcards.region] if wildcards.region in mapping else "location_exposure"

rule traits:
    message:
        """
        Inferring ancestral traits for {params.columns!s}
          - increase uncertainty of reconstruction by {params.sampling_bias_correction} to partially account for sampling bias
        """
    input:
        tree = "results/tree{region}.nwk",
        metadata = rules.adjust_metadata.output
    output:
        node_data = "results/traits{region}.json",
    params:
        columns = _get_exposure_trait_for_wildcards,
        sampling_bias_correction = 2.5
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence \
            --sampling-bias-correction {params.sampling_bias_correction} \
        """

rule clades:
    message: "Adding internal clade labels"
    input:
        tree = "results/tree{region}.nwk",
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = config["clades"]
    output:
        clade_data = "results/clades{region}.json"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.clade_data}
        """

rule colors:
    message: "Constructing colors file"
    input:
        ordering = rules.data_setup.output.ordering,
        color_schemes = config["color_schemes"],
        metadata = rules.adjust_metadata.output
    output:
        colors = "results/colors{region}.tsv"
    shell:
        """
        python3 scripts/assign-colors.py \
            --ordering {input.ordering} \
            --color-schemes {input.color_schemes} \
            --output {output.colors} \
            --metadata {input.metadata}
        """

rule recency:
    message: "Use metadata on submission date to construct submission recency field"
    input:
        metadata = rules.adjust_metadata.output
    output:
        "results/recency{region}.json"
    shell:
        """
        python3 scripts/construct-recency-from-submission-date.py \
            --metadata {input.metadata} \
            --output {output}
        """

rule tip_frequencies:
    message: "Estimating censored KDE frequencies for tips"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.adjust_metadata.output
    output:
        tip_frequencies_json = "auspice/ncov{region}_tip-frequencies.json"
    params:
        min_date = 2020.0,
        pivot_interval = 1,
        narrow_bandwidth = 0.05,
        proportion_wide = 0.0
    shell:
        """
        augur frequencies \
            --method kde \
            --metadata {input.metadata} \
            --tree {input.tree} \
            --min-date {params.min_date} \
            --pivot-interval {params.pivot_interval} \
            --narrow-bandwidth {params.narrow_bandwidth} \
            --proportion-wide {params.proportion_wide} \
            --output {output.tip_frequencies_json}
        """

def export_title(wildcards):
    region = wildcards.region

    if not region:
        return "Genomic epidemiology of novel coronavirus"
    elif region == "_global":
        return "Genomic epidemiology of novel coronavirus - Global subsampling"
    else:
        region_title = region.lstrip("_").replace("-", " ").title()
        return f"Genomic epidemiology of novel coronavirus - {region_title}-focused subsampling"

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.adjust_metadata.output,
        branch_lengths = rules.refine.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        traits = rules.traits.output.node_data,
        auspice_config = config["auspice_config"],
        colors = rules.colors.output.colors,
        lat_longs = rules.data_setup.output.lat_long,
        description = config["description"],
        clades = "results/clades{region}.json",
        recency = rules.recency.output
    output:
        auspice_json = "results/ncov_with_accessions{region}.json"
    params:
        title = export_title
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} {input.traits} {input.clades} {input.recency} \
            --auspice-config {input.auspice_config} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --title {params.title:q} \
            --description {input.description} \
            --output {output.auspice_json}
        """



rule incorporate_travel_history:
    message: "Adjusting main auspice JSON to take into account travel history"
    input:
        auspice_json = rules.export.output.auspice_json,
        colors = rules.colors.output.colors,
        lat_longs = rules.data_setup.output.lat_long
    params:
        sampling = _get_sampling_trait_for_wildcards,
        exposure = _get_exposure_trait_for_wildcards
    output:
        auspice_json = "results/ncov_with_accessions_and_travel_branches{region}.json"
    shell:
        """
        python3 ./scripts/modify-tree-according-to-exposure.py \
            --input {input.auspice_json} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --sampling {params.sampling} \
            --exposure {params.exposure} \
            --output {output.auspice_json}
        """



rule fix_colorings:
    message: "Remove extraneous colorings for main build"
    input:
        auspice_json = rules.incorporate_travel_history.output.auspice_json,
        auspice_json_no_travel = rules.export.output.auspice_json
    output:
        auspice_json = "auspice/ncov{region}.json",
        auspice_json_no_travel = "auspice/ncov{region}_no_travel_history.json"
    shell:
        """
        python scripts/fix-colorings.py \
            --input {input.auspice_json} \
            --output {output.auspice_json}
        python scripts/fix-colorings.py \
            --input {input.auspice_json_no_travel} \
            --output {output.auspice_json_no_travel}
        """



rule dated_json:
    message: "Copying dated Auspice JSON"
    input:
        auspice_json = rules.fix_colorings.output.auspice_json,
        tip_frequencies_json = rules.tip_frequencies.output.tip_frequencies_json
    output:
        dated_auspice_json = "auspice/ncov{region}_{date}.json",
        dated_tip_frequencies_json = "auspice/ncov{region}_{date}_tip-frequencies.json"
    shell:
        """
        cp {input.auspice_json} {output.dated_auspice_json}
        cp {input.tip_frequencies_json} {output.dated_tip_frequencies_json}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"


rule reprod_zip:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
