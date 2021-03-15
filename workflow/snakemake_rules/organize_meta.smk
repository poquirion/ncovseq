
rule download_latest_excludes:
    message: "get metadata from nextstrain"
    output:
        local_world_exclude = 'results/nextstrain_exclude{}.tsv'.format(datetime.date.today())
    shell:
        """
        # wget {config[exclude_latest]} -O {output.local_world_exclude}
        touch {output.local_world_exclude} # nothing comming from the outside on that run
        """

# rule aling_all_data:
#     message: "align all data up to june 1 and freeze data too"
#
#     input:
#         raw_fastq = config['gsaid_fasta'],
#         our_data = config['fasta_path'],
#         wuhan = config['wuhan']
#     output:
#         all_alined_fastq = 'intermediate/all_aligned.fasta'
#     shell:
#         """
#         """

rule data_setup:
    message: "Organize and merge metadata, fasta, exclude"
    input:
        fasta_path = config["fasta_path"],
        world_meta = config['world_meta'],
        world_exclude = rules.download_latest_excludes.output.local_world_exclude,
        local_exclude = config['local_exclude']
    output:
        sequences = "results/sequences.fasta",
        metadata = "results/merged_metadata.tsv",
        exclude = "results/exclude.txt",
        ordering = "results/ordering.tsv",
    shell:
        """
        cp {input.world_exclude} {output.exclude} 
        cat {input.local_exclude} >> {output.exclude}
        head -n 1 {input.world_meta}  > {output.metadata}
        grep -i "Canada/QC" {input.world_meta} | grep Moreira >> {output.metadata}
        grep "Wuhan/WH01/2019" {input.world_meta} >> {output.metadata}
        cp {input.fasta_path}  {output.sequences}
        cp config/ordering.tsv {output.ordering}
        """

rule lat_long:
    message: "Organize and merge metadata, fasta, exclude"
    output:
        lat_long = "results/lat_long.tsv"
    shell:
        """
        scripts/concat_location.py --out_lat_long {output.lat_long}
        """

