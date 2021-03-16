
rule adjust_metadata:
    input:
        rules.data_setup.output.metadata
    output:
        "results/metadata_adjusted.tsv"
    shell:
        """
        cp {input} {output}
        """

rule tree:
    message: "Building tree"
    input:
        alignment = "results/subsampled_alignment.fasta"
    output:
        tree = "results/tree_raw.nwk"
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
    return rules.mask.output.alignment

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
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json"
    threads: 1
    params:
        root = "--root Canada/Qc-L00218056/2020",
        clock_rate = 0.0008,
        clock_std_dev = 0.0004,
        coalescent = "const",
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
