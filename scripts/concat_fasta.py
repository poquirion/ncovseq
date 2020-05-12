#!/usr/bin/env python3
import argparse
import glob
import os
import shutil

def main():
    parser = argparse.ArgumentParser(
        description="concat_fasta",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input_dir', type=dir_path, required=True, help="dir with fasta")
    parser.add_argument('--extra_fasta_files', nargs='+', help="extra fasta file list to include in the " \
                                                                          "sequence")
    parser.add_argument('--output', type=str, required=True, help="output fasta")
    parser.add_argument('--id_format', default='Canada/Qc-{}/2020', type=str, help="format for publish "
                                                                                                 "id")

    args = parser.parse_args()

    in_dir = args.input_dir
    out_file = args.output
    in_extra = args.extra_fasta_files

    in_fasta = glob.glob("{}/*fasta".format(in_dir))
    id_format = args.id_format


    with open(out_file, 'w') as wfp:
        for f in in_fasta:
            with open(f, 'r') as fp:
                shutil.copyfileobj(fp, wfp)
        for f in in_extra:
            with(open(f, 'r')) as fp:
                shutil.copyfileobj(fp, wfp)


def dir_path(path):
    if os.access(path, os.R_OK):
        return path
    else:
        raise argparse.ArgumentTypeError("{path} is not readable")


if __name__ == '__main__':
    main()
