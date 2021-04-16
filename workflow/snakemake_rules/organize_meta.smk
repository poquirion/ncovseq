
rule download_latest_excludes:
    message: "get metadata from nextstrain"
    output:
        local_world_exclude = 'results/nextstrain_exclude{}.tsv'.format(datetime.date.today())
    shell:
        """
        # wget {config[exclude_latest]} -O {output.local_world_exclude}
        touch {output.local_world_exclude} # nothing comming from the outside on that run
        """

rule data_setup:
    message: "Organize and merge metadata, fasta, exclude"
    input:
        fasta_path = config["fasta_path"],
        qc_meta = config['qc_meta'],
        world_meta = config['world_meta'],
        world_exclude = rules.download_latest_excludes.output.local_world_exclude,
        local_exclude = config['local_exclude']
    output:
        sequences = "results/sequences.fasta",
        metadata = "results/merged_metadata.tsv",
        exclude = "results/exclude.txt",
        selected_lspq = "results/lspq_only_metadata.tsv",
        ordering = "results/ordering.tsv",
    params:
        extra_fasta = config['extra_fasta']
    shell:
        """
        scripts/concat_fasta.py --input_dir {input.fasta_path} --output {output.sequences}\
         {params.extra_fasta}
        scripts/concat_meta.py --inspq_meta {input.qc_meta} --nextstrain_metadata {input.world_meta}  --output \
        {output.metadata} --fasta_dir  {input.fasta_path} --output_lspq_only {output.selected_lspq} \
        --out_order {output.ordering} --keep_all_meta
        cp {input.world_exclude} {output.exclude} 
        cat {input.local_exclude} >> {output.exclude}
        """

rule lat_long:
    message: "Organize and merge metadata, fasta, exclude"
    output:
        lat_long = "results/lat_long.tsv"
    shell:
        """
        scripts/concat_location.py --out_lat_long {output.lat_long}
        """

