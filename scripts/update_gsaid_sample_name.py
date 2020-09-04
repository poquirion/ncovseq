#!/usr/bin/env python
import argparse
import copy
import csv

import Bio.SeqIO
from augur.align import read_sequences, write_seqs

# The gisaid id have changes over time, we getting them up to date here.


china_list = ["Anhui", "Beijing", "Foshan", "Fujian", "Fuyang", "Fuzhou", "Ganzhou", "Guangdong", "Guangzhou",
              "Hangzhou", "Harbin", "Hefei", "Henan", "Jian", "Jiangsu", "Jiangxi", "Jingzhou", "Jiujiang", "Liaoning",
              "Lishui", "NanChang", "Nanchang", "Pingxiang", "Shandong", "Shanghai", "Shangrao", "Shaoxing", "Shenzhen",
              "Sichuan", "Tianmen", "Weifang", "Wuhan", "Xinyu", "Yichun", "Yingtan", "Yunnan", "Zhejiang"]

manual_fix = {"Argentina/PAIS_A007/2020": "Argentina/PAIS-A0007/2020",
"Argentina/PAIS_A006/2020": "Argentina/PAIS-A0006/2020",
"Singapore/320nan/2020": "",
"Andalusia/COV001884/2020": "Spain/AN-001884/2020",
"Indonesia/JKT-EIJK02/2020": "Indonesia/JK-EIJK-02/2020",
"Brazil/AP161167-IEC/2020": "Brazil/AP-IEC-161167/2020",
"Mexico/EdoMex-InDRE_03/2020": "Mexico/CMX-InDRE_03/2020",
"Indonesia/JKT-EIJK01/2020": "Indonesia/JK-EIJK-01/2020",
"Argentina/PAIS_A001/2020": "Argentina/PAIS-A0001/2020",
"Brazil/AC162535-IEC/2020": "Brazil/AC-IEC-162535/2020",
"Brazil/PA161548-IEC/2020": "Brazil/PA-IEC-161548/2020",
"Germany/HH-1/2020": "Germany/HH-UMC-01/2020",
"Germany/BavPat1/2020": "",
"Italy/Siena-1/2020": "Italy/TUS-UniSI-1/2020",
"Italy/IZSPB_45/2020": "Italy/APU-UniMI-45/2020",
"Spain/Madrid_H12_1804/2021": "Spain/Madrid-H12-1804/2020",
"Italy/484/2020": "Italy/ABR-IZSGC-484/2020",
"Italy/INMI6/2020": "Italy/LAZ-INMI-6/2020",
"Italy/UnivPM1/2020": "Italy/MAR-UnivPM-1/2020",
"Italy/INMI8/2020": "Italy/LAZ-INMI-8/2020",
"Italy/206/2020": "Italy/EMR-AMVRC-206/2020",
"Italy/INMI9/2020": "Italy/LAZ-INMI-9/2020",
"USA/CA-QDX-03-2/2020": "USA/CA-QDX-03/2020",
"USA/NY-QDX-01-2/2020": "USA/NY-QDX-01/202000",
"USA/MI-MDHHS-SC20007/2020": "USA/MI-SC2-0007/2020",
"USA/MI-MDHHS-SC20008/2020": "USA/MI-SC2-0008/2020",
"Sweden/COV001953/2020": "Spain/AN-001953/2020"}




def return_match_names(confused_ids, right_ids):
    correction_dico = {}

    working_copy = set(confused_ids)

    tweak_dico = {}
    for wc in confused_ids:
        orig_wc = copy.copy(wc)

        south_korea_wc = None
        brazil_wc = None
        dusseldorf_wc = None
        italie_wc = None
        india_wc = None
        bavaria_wc = None

        if wc.startswith('Korea'):
            south_korea_wc = 'South{}'.format(wc)
        elif wc.startswith('Brazil'):
            brazil_wc = wc.replace('BR', '')
        elif wc.startswith('Germany/NRW-'):
            dusseldorf_wc = wc.replace('/NRW-', '/NW-HHU-')
        elif wc.startswith('Germany'):
            bavaria_wc = wc.replace('MVP', 'MVP-')
        elif wc.startswith('Italy') and 'Feb' in wc:
            c, i, y = wc.split('/')
            italie_wc ='{}/{}/{}'.format(c, i[0:-6], y)
        elif wc.startswith('Italy') and 'Mar' in wc:
            c , i, y = wc.split('/')
            italie_wc ='{}/{}/{}'.format(c, i[0:-5], y)
        elif wc.startswith('India'):
            india_wc = '{}/-{}/{}'.format(*wc.split('/'))

        def cleaup(the_wc):
            correction_dico[wc] = the_wc
            working_copy.remove(wc)

        if wc in right_ids:
            cleaup(wc)
        elif wc in manual_fix.keys():
            cleaup(manual_fix[wc])
        elif south_korea_wc in right_ids:
            cleaup(south_korea_wc)
        elif brazil_wc in right_ids:
            cleaup(brazil_wc)
        elif india_wc in right_ids:
            cleaup(india_wc)
        elif dusseldorf_wc in right_ids:
            cleaup(dusseldorf_wc)
        elif italie_wc in right_ids:
            cleaup(italie_wc)
        elif italie_wc:
            tweak_dico[italie_wc] = wc
            working_copy.remove(wc)
            working_copy.add(italie_wc)
        elif india_wc:
            tweak_dico[india_wc] = wc
            working_copy.remove(wc)
            working_copy.add(india_wc)
        elif bavaria_wc:
            tweak_dico[bavaria_wc] = wc
            working_copy.remove(wc)
            working_copy.add(bavaria_wc)


    confused_split = [i.split('/') for i in working_copy]

    print(len(confused_split))
    # do a {country_name : { year: id  }  } split for the right ids
    from collections import defaultdict
    finder_dico = defaultdict(lambda: defaultdict(list))
    for right_id in right_ids:
        country, uid, year = right_id.split('/')
        finder_dico[country][year].append(uid)


    def resolver(cs):
        """
        :param cs: the confused id split in 3 [country, uid, year]
        :return: the right name "country/uid/year"
        """
        # select the right country and year
        china_region_dico = None
        country= cs[0]
        year = cs[2]
        if country == 'China':
            china_region_dico= {elem: k for k in china_list for elem in finder_dico[k][cs[2]]}
            id_list = china_region_dico.keys()
        else:
            id_list = finder_dico[country][year]
        # reverse match
        to_match = cs[1].replace('_', '-')

        for i in range(len(to_match)-1, -1, -1):
            chars = to_match[i:]
            match = [j for j in id_list if j.endswith(chars)]
            if len(match) == 1:
                # print('Match! {} =~ {}'.format(match, to_match))
                if china_region_dico is not None:
                    ret_str = '{}/{}/{}'.format(china_region_dico[match[0]], match[0], year)
                else:
                    ret_str = '{}/{}/{}'.format(country, match[0], year)
                return ret_str
            elif len(match) == 0:
                # try stategy for italy
                break
        print('found no match for {} looking for other strategies'.format('/'.join(cs)))
        return 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'

    for cs in confused_split:
        good_name = resolver(cs)  # watch out cs is modified
        bad_name = '/'.join(cs)
        if tweak_dico.get(bad_name):
            bad_name = tweak_dico[bad_name]
        correction_dico[bad_name] = good_name
        # print(bad_name, good_name)

    return correction_dico




def main():
    parser = argparse.ArgumentParser(
        description="concat_meta_tsv",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser = argparse.ArgumentParser()
    parser.add_argument("gsaid_meta", help="gsaid .tsv metadata")
    parser.add_argument("input_fasta", help="fasta file with outdated sample names")

    parsed = parser.parse_args()

    meta_path = parsed.gsaid_meta
    fasta_path = parsed.input_fasta

    with open(meta_path) as fp :
        meta_reader = csv.reader(fp, delimiter='\t')
        headers = next(meta_reader)

        # We reject all animal and environmental samples, they are the one in four pieces
        meta_id = set([row[0] for row in meta_reader if len(row[0].split('/')) <= 3])

    fastas = read_sequences(fasta_path)

    translator = return_match_names(fastas.keys(), meta_id)
    for k, v in translator.items():
        if k != v:
            fastas[k].name = v

    Bio.SeqIO.write(fastas.values(), '/tmp/sequence_gisad.fasta', 'fasta')


if __name__ == '__main__':
    main()