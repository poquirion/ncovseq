#!/home/poq/covid/nextstrain/env_augur/bin/python
import argparse

from Bio import SeqIO


def main():
    """
    extact subsample from fasta file
    :return:
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("--output", "-o", help='output_fata fasta', type=argparse.FileType('w'), default='-')
    parser.add_argument("-v", help='Exclude selection', action='store_true')
    parser.add_argument("-n", help='stop after n match', type=int, default=0)
    parser.add_argument("input", help='input fasta', type=argparse.FileType('r'), default='-')
    parser.add_argument("condition", help='grep on > line')
    args = parser.parse_args()

    condition = args.condition
    n_max = args.n
    n = 0 

    if not args.v:
        for record in SeqIO.parse(args.input, format='fasta'):
            if condition in record.id:
                SeqIO.write(record, args.output, format='fasta')
                if n_max and n + 1 >= n_max:
                    break
                else:
                    n = n + 1
    else:
        for record in SeqIO.parse(args.input, format='fasta'):
            if condition not in record.id:
                SeqIO.write(record, args.output, format='fasta')
                if n_max and n + 1 >= n_max:
                    break
                else:
                    n = n + 1

if __name__ == '__main__':
    main()
