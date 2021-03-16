import pandas as pd


data = pd.read_csv("/home/poq/meta/metadataQC_postprocessing2.csv", sep=',')
data['B.1.350'] = 'other'
data['B.1.3'] = 'other'
data['B.1.147'] = 'other'
data.loc[data['cladePang'] == 'B.1.350','B.1.350'] = 'B.1.350'
data.loc[data['cladePang'] == 'B.1.3','B.1.3'] = 'B.1.3'
data.loc[data['cladePang'] == 'B.1.147','B.1.147'] = 'B.1.147'
data['pangolin_large'] = data['cladePang']
data.loc[(data['cladePang'] != 'B.1.147') & (data['cladePang'] != 'B.1.350') & (data['cladePang'] != 'B.1.3') & (data['cladePang'] != 'B.1' ) , 'pangolin_large'] = 'other'
data.to_csv("/home/poq/meta/withpango_switch.tsv", sep='\t', index=False)
