"""
Add column to metadata to denote 'focal' samples based on supplied region
Rewrite location, division and country for non-focal samples to be region
Rewrite division_exposure and country_exposure for non-focal samples to be region_exposure
"""

import argparse
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Add column to metadata to denote 'focal' samples based on supplied region",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata", type = str, required=True, help="metadata")
    parser.add_argument("--region", type=str, required=True, help="focal region")
    parser.add_argument("--output", type=str, required=True, help="adjusted metadata")
    args = parser.parse_args()

    fix_casing = {
        "quebec": "Quebec"
    }
    # focal_region = fix_casing[args.region]
    focal_region = "North America"
    print("Adjusting metadata for focal region", focal_region)

    metadata = pd.read_csv(args.metadata, delimiter='\t')
    metadata.insert(12, 'focal', True)

    metadata.loc[metadata.region != focal_region, 'focal'] = False
    metadata.loc[metadata.region != focal_region, 'province'] = metadata.region
    metadata.loc[metadata.region != focal_region, 'location'] = metadata.region
    metadata.loc[metadata.region != focal_region, 'country'] = metadata.country
    # metadata.loc[(metadata.division == focal_region) & (metadata.division_exposure != focal_region),
    #              'location_exposure'] = metadata.country_exposure
    # metadata.loc[(metadata.division == focal_region) & (metadata.division_exposure != focal_region), 'rss_exposure'] = metadata.country_exposure
    # metadata.loc[(metadata.division == focal_region) & (metadata.division_exposure != focal_region), 'country_exposure'] = metadata.country_exposure

    metadata.to_csv(args.output, index=False, sep="\t")
'''
EricF comment
    fix_casing = {
        "asia": "Asia",
        "africa": "Africa",
        "europe": "Europe",
        "north america": "North America",
        "oceania": "Oceania",
        "south america": "South America"
    }
    focal_region = fix_casing[args.region]
'''


'''
EricF comment
    metadata.loc[metadata.region != focal_region, 'focal'] = False
    metadata.loc[metadata.region != focal_region, 'location'] = ""
    metadata.loc[metadata.region != focal_region, 'division'] = metadata.region
    metadata.loc[metadata.region != focal_region, 'country'] = metadata.region
    metadata.loc[metadata.region != focal_region, 'division_exposure'] = metadata.region_exposure
    metadata.loc[metadata.region != focal_region, 'country_exposure'] = metadata.region_exposure
    metadata.loc[(metadata.region == focal_region) & (metadata.region_exposure != focal_region), 'division_exposure'] = metadata.region_exposure
    metadata.loc[(metadata.region == focal_region) & (metadata.region_exposure != focal_region), 'country_exposure'] = metadata.region_exposure
'''

