# conbine brekspoints in 24 samples "

import pandas as pd
import sys,os

data = pd.DataFrame()

for file in sys.argv[1:]:
	df = pd.read_table(file)
	ID = os.path.basename(file).split('.')[0].replace('_css','')
	df = df.rename(columns = {'freq': ID })

	if file == sys.argv[1]:
		data = df
	else:
		data = pd.merge(data, df, on = ['chr', 'position'], how= 'outer')

out = '/haplox/users/wangzy/04.FineCNV/09.SFClip/out/'
	
data.to_csv( out + '24_Combined.gDNA.dedup.sfc_breakpoints.location.txt',index = False, sep = '\t',na_rep ='0', float_format = '%d')
