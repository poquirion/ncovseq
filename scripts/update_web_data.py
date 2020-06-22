#!/usr/bin/env python3
import argparse
import csv
import datetime
import glob
import os
import shutil
import tarfile
import zipfile

import pandas as pd
import unidecode
from country_list import countries_for_language

PUBLISH_META_COLUMN = ['strain', 'virus', 'date', 'date_submitted', 'province', 'age', 'sex', 'originating_lab', 'submitting_lab', 'url']

DIVISIONS = ['country', 'location', 'province']
SCRIPT_DIR = os.path.dirname(__file__)


def count_fasta_len(fastas):
    id_len = {}
    for fasta in fastas:
        with open(fasta) as fp:
            fasta_id = fp.readline().lstrip('>').rstrip('\n')
            id_len[fasta_id.split()[0]] = len(fp.read())

    return id_len



def main():
    parser = argparse.ArgumentParser(
        description="concat_meta_tsv",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser = argparse.ArgumentParser()
    parser.add_argument("--inspq_meta", default='/genfs/projects/COVID_consensus/metadata/sgil_extract.tsv'
                        , help="The LNSPQ .tsv")
    parser.add_argument("--fasta_dir", default='/genfs/projects/COVID_consensus/FinalPublished', help="The modified .tsv")
    parser.add_argument("--output_lspq", default='lspq_metadata.tsv', help="Approved lspq .tsv")
    parser.add_argument('--id_format', default='Canada/Qc-{}/2020', type=str, help="format for publish "
                                                                                                  "id")
    parser.add_argument('--all_fasta', type=str, help="all_fasta.fasta.tgz")
    parser.add_argument('--per_hospital_dir', default='/tmp/hospitals', type=str, help="fasta_split_per_hospital")


    args = parser.parse_args()
    lnspq_path = args.inspq_meta
    fasta_dir = args.fasta_dir
    id_format = args.id_format
    approved_lspq_tsv = args.output_lspq
    all_fasta = args.all_fasta
    per_hospital_dir = args.per_hospital_dir

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

    def updateid(old_id):
        return id_format.format(old_id)
    lnspq_df['strain'] = lnspq_df['strain'].apply(updateid)

    fastas = glob.glob("{}/*.fasta".format(fasta_dir))
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
    # lnspq_df['authors'] = 'Moira et al'
    lnspq_df['region'] = 'North America'
    lnspq_df['submitting_lab'] = 'LSPQ'
    lnspq_df['date_submitted'] = datetime.datetime.today().strftime('%Y-%m-%d')



    if fasta_dir:
        lnspq_df['url'] = 'http://www.covseq.ca/data-info'
    else:
        lnspq_df['url'] = ''

    def only_valid(in_df):
        return in_df[PUBLISH_META_COLUMN]

    approved_lspq = only_valid(lnspq_df)
    os.makedirs('{}/all/'.format(per_hospital_dir), exist_ok=True)
    meta_data = '{}/{}'.format(per_hospital_dir, approved_lspq_tsv)
    approved_lspq.to_csv(meta_data, sep='\t', index=False)
    o_lab_group = lnspq_df.groupby('originating_lab')
    all_tar = tarfile.open('{}/all/all_fasta.tgz'.format(per_hospital_dir), "w:gz")
    all_zipMe = zipfile.ZipFile('{}/all/all_fasta.zip'.format(per_hospital_dir), 'w')

    seq_andfasta_tar = tarfile.open('{}/all_fasta_and_meta.tgz'.format(per_hospital_dir), "w:gz")
    seq_andfasta_zipMe = zipfile.ZipFile('{}/all_fasta_and_meta.zip'.format(per_hospital_dir), 'w')
    seq_andfasta_zipMe.write(meta_data, compress_type=zipfile.ZIP_DEFLATED
                             , arcname=os.path.basename(meta_data).replace('.consensus', ''))
    seq_andfasta_tar.add(meta_data, arcname=os.path.basename(meta_data).replace('.consensus', ''))
    for hospital, val in o_lab_group:
        strains = val.strain
        ids = [s.split('-')[1].split('/')[0] for s in strains]
        hospital = unidecode.unidecode(hospital).lower()\
                                                        .replace('(', '')\
                                                        .replace(')', '').replace("'", '_')\
                                                        .replace('/', '_').replace(',', '')\
                                                        .replace(' ', '_').split('_hopital', 1)[0].rstrip('-/')

        hospital_path = '{}/{}'.format(per_hospital_dir, hospital)
        os.makedirs(hospital_path, exist_ok=True)
        tar = tarfile.open('{}/{}_fasta.tgz'.format(hospital_path, hospital), "w:gz")
        zipMe = zipfile.ZipFile('{}/{}_fasta.zip'.format(hospital_path, hospital), 'w')

        for virus in ids:
            try:
                fasta = glob.glob("{}/{}*.fasta".format(fasta_dir, virus))[0]
                dest_fasta = os.path.basename(fasta.replace('.consensus', ''))
                os.chmod(fasta, int('664', base=8))
                shutil.copy2(fasta, '{}/{}'.format(hospital_path, dest_fasta))
                zipMe.write(fasta, compress_type=zipfile.ZIP_DEFLATED
                            , arcname=os.path.basename(fasta).replace('.consensus', ''))
                tar.add(fasta, arcname=os.path.basename(fasta).replace('.consensus', ''))
                all_zipMe.write(fasta, compress_type=zipfile.ZIP_DEFLATED
                            , arcname=os.path.basename(fasta).replace('.consensus', ''))
                all_tar.add(fasta, arcname=os.path.basename(fasta).replace('.consensus', ''))
                seq_andfasta_zipMe.write(fasta, compress_type=zipfile.ZIP_DEFLATED
                            , arcname=os.path.basename(fasta).replace('.consensus', ''))
                seq_andfasta_tar.add(fasta, arcname=os.path.basename(fasta).replace('.consensus', ''))
            except IndexError:
                pass
        tar.close()
        zipMe.close()
    all_tar.close()
    all_zipMe.close()
    seq_andfasta_zipMe.close()
    seq_andfasta_tar.close()
    print(per_hospital_dir)

    # then
    # rsync -lrv --progress  hospitals/ centos@covseq.ca:/data



if __name__ == '__main__':
    main()
