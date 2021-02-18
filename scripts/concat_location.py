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

DIVISIONS = ['country', 'location', 'province', 'rss']


def genetate_lat_long(file_name='results/lat_long.tsv', misc_position='config/misc_lat_long.tsv',
                      lat_long_list=None):
    if lat_long_list is None:
        lat_long_list = ['config/country_lat_long.tsv',
                         'config/canada_lat_long.tsv',
                         'config/province_lat_long.tsv',
                         'config/rss_lat_long.tsv']

    misc_df = pd.read_csv(misc_position, sep='\t')
    # pick only contry on misc
    misc_df = misc_df[misc_df['type'] == 'country']
    misc_df.rename(columns={'id': 'location'}, inplace=True)

    result_df = pd.DataFrame()
    for pos_file in lat_long_list:
        result_df = result_df.append(pd.read_csv(pos_file, sep='\t', header=0))

    # add misc not in the list on top of that.
    result_df = result_df.append(misc_df[['location', 'lat', 'long']][~misc_df['location'].isin(result_df['location'])],
                             sort=False)

    # remove file if exist because we are in 'a' mode
    try:
        os.remove(file_name)
    except OSError:
        pass

    for d in DIVISIONS:
        result_df.insert(0, column='type', value=d)
        result_df.to_csv(file_name, sep='\t', index=False, header=False, mode='a')
        del result_df['type']



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
    parser.add_argument('--out_lat_long', type=str, help="lat_long output file", default='results/lat_long.tsv')

    args = parser.parse_args()
    out_lat_long = args.out_lat_long

    genetate_lat_long(file_name=out_lat_long)

    # create_ordering(final_df, file_name=out_order)


if __name__ == '__main__':
    main()
