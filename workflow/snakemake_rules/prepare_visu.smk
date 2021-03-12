rule ancestral:
    message:
        """
        Reconstructing ancestral sequences and mutations
          - inferring ambiguous mutations
        """
    input:
        tree = rules.refine.output.tree,
        alignment = _get_alignments_for_tree
    output:
        node_data = "results/nt_muts.json"
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
        node_data = "results/haplotype_status.json"
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
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = config["reference"]
    output:
        node_data = "results/aa_muts.json"
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
          - increase uncertainty of reconstruction by {params.sampling_bias_correction} to 
          partially account for sampling bias
        """
    input:
        tree = rules.refine.output.tree,
        metadata = rules.adjust_metadata.output
    output:
        node_data = "results/traits_{trait_rez}.json"
    params:
        columns = "{trait_rez}",
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
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = config["clades"]
    output:
        clade_data = "results/clades.json"
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
        colors = "results/colors.tsv"
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
        "results/recency.json"
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
        tip_frequencies_json = "auspice/ncov_tip-frequencies.json"
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
        return "Quebec focused Genomic epidemiology of novel coronavirus"

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.adjust_metadata.output,
        branch_lengths = rules.refine.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        traits = expand(rules.traits.output.node_data, trait_rez=TRAIT_REZ, allow_missing=True),
        auspice_config = config["auspice_config"],
        colors = rules.colors.output.colors,
        lat_longs = rules.lat_long.output.lat_long,
        description = config["description"],
        clades = "results/clades.json",
        recency = rules.recency.output
    output:
        auspice_json = "results/ncov_with_accessions.json"
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


rule fix_colorings:
    message: "Remove extraneous colorings for main build"
    input:
        auspice_json = rules.export.output.auspice_json,
    output:
        auspice_json = "auspice/ncov.json",
    shell:
        """
        python scripts/fix-colorings.py \
            --input {input.auspice_json} \
            --output {output.auspice_json}
        """



rule dated_json:
    message: "Copying dated Auspice JSON"
    input:
        auspice_json = rules.fix_colorings.output.auspice_json,
        tip_frequencies_json = rules.tip_frequencies.output.tip_frequencies_json
    output:
        dated_auspice_json = "auspice/ncov_{}.json".format(get_todays_date()),
        dated_tip_frequencies_json = "auspice/ncov_{}_tip-frequencies.json".format(get_todays_date())
    shell:
        """
        cp {input.auspice_json} {output.dated_auspice_json}
        cp {input.tip_frequencies_json} {output.dated_tip_frequencies_json}
        """
