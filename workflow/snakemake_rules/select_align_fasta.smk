
rule add_canada:
    message:
        """
        # add a minimum amount of seq from Specific country
        """
    input:
        sequences = rules.data_setup.output.sequences,
        metadata = rules.data_setup.output.metadata,
        exclude = rules.data_setup.output.exclude
    output:
        sequences = "results/filtered_canada.fasta"
    params:
        min_length = config["min_length"],
        max_date = config["max_date"]
    threads: 1
    shell:
        """augur filter \
         --sequences-per-group 100 \
         --group-by country \
         --min-length {params.min_length} \
         --max-date {params.max_date} \
         --exclude-where  "country!=Canada" \
         --metadata {input.metadata} \
         --sequences {input.sequences} \
         --priority data/priority.txt  \
         -o {output.sequences}
         echo contatenated fasta: `grep '>' {output.sequences} | wc -l` 
        """


rule add_neighbour:
    message:
        """
        # add a minimum amount of seq from Specific country
        # the neigbour list is  set in scripts/concat_meta.py
        """
    input:
        sequences = rules.data_setup.output.sequences,
        metadata = rules.data_setup.output.metadata,
        exclude = rules.data_setup.output.exclude
    output:
        sequences = "results/filtered_neigbour.fasta"
    params:
        min_length = config["min_length"],
        max_date = config["max_date"]
    threads: 1
    shell:
        """augur filter \
         --sequences-per-group 15 \
         --group-by division \
         --min-length {params.min_length} \
         --max-date {params.max_date} \
         --exclude-where  "neighbour=no" \
         --metadata {input.metadata} \
         --sequences {input.sequences} \
         --priority data/priority.txt  \
         -o {output.sequences}
         echo contatenated fasta: `grep '>' {output.sequences} | wc -l` 
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
        sequences = "results/filtered_global.fasta"
    threads: 1
    params:
        min_length = config["min_length"],
        exclude_where = config["exclude_where"],
        group_by = config["group_by"],
        sequences_per_group = config["sequences_per_group"],
        max_date = config["max_date"]
    shell:
        """
        augur filter \
            --max-date {params.max_date} \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --include {input.include} \
            --exclude {input.exclude} \
            --exclude-where {params.exclude_where}\
            --min-length {params.min_length} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --output {output.sequences} \
            --priority data/priority.txt
            echo contatenated fasta: `grep '>' {output.sequences} | wc -l` 
        """



rule merge_fasta_selection:
    message:
        """
        Combine and deduplicate FASTAs
        """
    input:
        rules.filter.output.sequences,
        # rules.add_canada.output.sequences,
        # rules.add_neighbour.output.sequences
    output:
        sequences = "results/filtered.fasta"
    shell:
        """
        python3 scripts/combine-and-dedup-fastas.py \
            --input {input} \
            --output {output.sequences}
        """




checkpoint partition_sequences:
    input:
        sequences = rules.merge_fasta_selection.output.sequences
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
    threads: 2
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
        "results/subsampled_alignment.fasta"
    shell:
        """
        cp {input} {output}
        """
