#!/usr/bin/env python3
import argparse
import datetime
import csv
import glob
import os
import random

import numpy as np
import pandas as pd
import unidecode
from country_list import countries_for_language



PUBLISH_META_COLUMN = ['strain', 'virus', 'date','date_submitted', 'country', 'location', 'province', 'location_exposure',
                       'ct', 'age', 'sex', 'originating_lab', 'submitting_lab', 'url', 'neighbour']

DIVISIONS = ['country', 'location', 'province']


def count_fasta_len(fastas):
    id_len = {}
    for fasta in fastas:
        with open(fasta) as fp:
            fasta_id = fp.readline().lstrip('>').rstrip('\n').split()[0]
            id_len[fasta_id] = len(fp.read())

    return id_len


def genetate_lat_long(df, file_name='results/lat_long.tsv', misc_position='config/misc_lat_long.tsv',
                      lat_long_list=None):
    if lat_long_list is None:
        lat_long_list = ['config/country_lat_long.tsv',
                         'config/canada_lat_long.tsv',
                         'config/province_lat_long.tsv']

    misc_df = pd.read_csv(misc_position, sep='\t')
    misc_df = misc_df[misc_df['type'] == 'country']
    misc_df.rename(columns={'id': 'location'}, inplace=True)

    country = pd.DataFrame()
    for pos_file in lat_long_list:
        country = country.append(pd.read_csv(pos_file, sep='\t', header=0))

    country = country.append(misc_df[['location', 'lat', 'long']][~misc_df['location'].isin(country['location'])],
                             sort=False)

    # remove file if exist because we are in 'a' mode
    try:
        os.remove(file_name)
    except OSError:
        pass

    for d in DIVISIONS:
        country.insert(0, column='type', value=d)
        country.to_csv(file_name, sep='\t', index=False, header=False, mode='a')
        del country['type']



def create_ordering(df, file_name='results/ordering.tsv'):
    tops = ['Quebec']
    # tops = ['Quebec', 'Canada', 'USA', 'Mexico', 'Italia', 'France']
    with open(file_name, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for d in DIVISIONS:
            places = set(df[d].unique())
            try:
                places.update(df[d+'_exposure'].unique())
            except KeyError:
                pass
            # clean up nan and sort
            places = sorted([p for p in places if p == p])
            for t in reversed(tops):
                if t in places:
                    places.remove(t)
                    places.insert(0, t)
            for p in places:
                writer.writerow([d, p])


def main():
    parser = argparse.ArgumentParser(
        description="concat_meta_tsv",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser = argparse.ArgumentParser()
    parser.add_argument("--inspq_meta", default='data/sgil_extract.tsv', help="The LNSPQ .tsv")
    parser.add_argument("--nextstrain_metadata", help="The .tsv that comes by "
                                                                                       "with nextstrain")
    parser.add_argument("--fasta_dir", default=None, help="The modified .tsv")
    parser.add_argument("--output", default='results/merged_metadata.tsv', help="The modified .tsv")
    parser.add_argument("--output_lspq_only", default='results/lspq_only_metadata.tsv', help="Aproved lspq .tsv")
    parser.add_argument("--subsample", '-s', type=int,
                        default=None, help="A metadata file with N sample and all QC ones")
    parser.add_argument('--id_format', default='Canada/Qc-{}/2020', type=str, help="format for publish "
                                                                                                  "id")
    parser.add_argument('--out_order', type=str, help="ordering for auspice visualisation")
    parser.add_argument('--out_lat_long', type=str, help="lat_long output file")

    args = parser.parse_args()
    in_path = args.nextstrain_metadata
    out_path = args.output
    subsample = args.subsample
    lnspq_path = args.inspq_meta
    fasta_dir = args.fasta_dir
    id_format = args.id_format
    approved_lspq_tsv = args.output_lspq_only
    out_order = args.out_order
    out_lat_long = args.out_lat_long


    gsaid_df = pd.read_csv(in_path, sep='\t')
    lnspq_df = pd.read_csv(lnspq_path, sep='\t')

    # traduction from fr to en
    traduc = {'NO_LSPQ': 'strain',
              'AGE': 'age',
              'SEX': 'sex',
              'MAX_CT': 'ct',
              'RSS_PATIENT': 'rss',
              'VOYAGE_PAYS_1': 'country_exposure',
              'DATE_PRELEV': 'date',
              'DATE_RECU': 'date_submitted',
              'CH': 'originating_lab',
              'POSTAL_CODE': 'rta'}
    lnspq_df.rename(columns=traduc, inplace=True)

    lnspq_df['country_exposure'] = [unidecode.unidecode(p.title()) for p in lnspq_df['country_exposure']]
    # lnspq_df['date_submitted'] = "{}s".format(datetime.datetime.today())
    pays_qc = {k: unidecode.unidecode(v) for k, v in countries_for_language('fr_CA')}

    pays_anglo = dict(countries_for_language('en'))
    # Fix non stadard names
    pays_anglo['US'] ='USA'
    pays_anglo['HK'] = 'Hong Kong'
    pays_anglo['CZ'] = 'Czech Republic'
    pays_anglo['CD'] = 'Democratic Republic of the Congo'

    trans = {pays_qc[code]: pays_anglo[code] for code in pays_qc.keys()}
    trans['Rep. Dominicaine'] = 'Dominican Republic'
    trans['cuba'] = 'Cuba'
    trans['Aucun_Voyage'] = '?'
    lnspq_df['country_exposure'].replace(trans, inplace=True)
    lnspq_df['location_exposure'] = lnspq_df['country_exposure']
    lnspq_df['province_exposure'] = lnspq_df['country_exposure']
    lnspq_df['neighbour'] = 'no'

    def updateid(old_id):
        return id_format.format(old_id)
    lnspq_df['strain'] = lnspq_df['strain'].apply(updateid)

    fastas = glob.glob("{}/*fasta".format(fasta_dir))
    fasta_id_len = count_fasta_len(fastas)
    for fid, l in fasta_id_len.items():
        lnspq_df.loc[lnspq_df['strain'] == fid, 'lenth'] = int(l)


    # Keep in only sample present in the fasta
    lnspq_df.drop(lnspq_df[lnspq_df['lenth'].isna()].index, inplace=True)

    lnspq_df['province'] = 'Quebec'
    lnspq_df['virus'] = 'ncov'
    lnspq_df['title'] = 'CoVSeQ - Covid Sequencing Quebec'
    lnspq_df['country'] = 'Canada'
    lnspq_df['location'] = 'Quebec'
    lnspq_df['division'] = 'Quebec'
    lnspq_df['region'] = 'North America'
    lnspq_df['submitting_lab'] = 'LSPQ'
    lnspq_df['date_submitted'] = datetime.datetime.today().strftime('%Y-%m-%d')



    if fasta_dir:
        lnspq_df['url'] = 'http://www.covseq.ca/data/'
    else:
        lnspq_df['url'] = ''

    # add rta and rss entry to world

    # still need to fix 'Iles Turques-Caiques' and 'Iles Vierges (E-U)',

    neighbour = ['New York', 'Ontario', 'Vermont', 'New Hampshire',
                  "Massachusetts", 'Maine', 'New Brunswick', 'Grand Princess']

    gsaid_df['neighbour'] = 'no'
    gsaid_df.loc[gsaid_df['division'].isin(neighbour), 'neighbour'] = 'yes'

    # we will do only location and province
    # gsaid_df.loc[gsaid_df['region'] != 'North America', 'rss'] = gsaid_df['country']
    # gsaid_df.loc[gsaid_df['region'] != 'North America', 'rta'] = gsaid_df['country']
    # gsaid_df.loc[gsaid_df['region'] == 'North America', 'rss'] = gsaid_df['country']
    # gsaid_df.loc[gsaid_df['region'] == 'North America', 'rta'] = gsaid_df['country']
    # gsaid_df.loc[gsaid_df['division'].isin(neighbourg), 'rss'] = gsaid_df['division']

    gsaid_df['location_exposure'] = gsaid_df['country_exposure']
    gsaid_df['province_exposure'] = gsaid_df['country_exposure']
    gsaid_df.loc[gsaid_df['country'] == 'Canada', 'province'] = gsaid_df['division']
    gsaid_df.loc[gsaid_df['country'] != 'Canada', 'province'] = gsaid_df['country']
    gsaid_df['location'] = gsaid_df['country']

    def only_valid(in_df):
        return in_df[PUBLISH_META_COLUMN]

    approved_lspq = only_valid(lnspq_df)
    approved_lspq.to_csv(approved_lspq_tsv, sep='\t', index=False)

    final_df = pd.concat([approved_lspq, gsaid_df], sort=False)
    final_df.drop_duplicates(subset='strain', keep="first", inplace=True)
    final_df.to_csv(out_path, sep='\t', index=False)
    create_ordering(final_df, file_name=out_order)
    genetate_lat_long(final_df, file_name=out_lat_long)

    if subsample:
        name = os.path.basename(out_path)
        path = os.path.dirname(out_path)
        s_path = '{}/sampled_{}'.format(path, name)
        print('subsample with {} point in {}'.format(subsample, s_path))
        s_df = gsaid_df.iloc[random.sample(range(len(gsaid_df)), subsample)]

        # make sure root virus is in data
        extra = []
        for s in open('../config/include.txt').read().splitlines():
            if s not in s_df['strain']:
                extra = extra + [gsaid_df.loc[gsaid_df['strain'] == s]]

        pd.concat([lnspq_df, s_df] + extra).to_csv(s_path, sep='\t', index=False)



if __name__ == '__main__':
    main()
