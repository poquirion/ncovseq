#!/home/poq/covid/nextstrain/env_augur/bin/python
import argparse

from Bio import SeqIO


def main():
    """
    extact subsample a fasta file using a list of ids from a file or from stdin
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("--output", "-o", help='output_fata fasta', type=argparse.FileType('w'), default='-')
    parser.add_argument("list", help='list of ids to extract', type=argparse.FileType('r'))
    parser.add_argument("input", help='input fasta', type=argparse.FileType('r'), default='-')
    args = parser.parse_args()

    all_ids = args.list
    records = SeqIO.to_dict(SeqIO.parse(args.input, format='fasta'))

    SeqIO.write([records[k.strip()] for k in all_ids if records.get(k.strip(), None)], args.output, format='fasta')



if __name__ == '__main__':
    main()
