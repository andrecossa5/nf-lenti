#!/usr/bin/python

from mito_utils.preprocessing import *


path_ = '/Users/IEO5505/Desktop/MI_TO/mito_preprocessing/scratch'
pp_method = 'freebayes'

os.chdir(os.path.join(path_, pp_method))

for x in os.listdir():
    
    os.chdir(x)
    print(x)

    afm = sc.read('afm.h5ad')
    afm.uns['scLT_system'] = 'MAESTER'
    afm.write('afm.h5ad')

    cells = pd.read_csv(f'../../maegatk/{x}/cells_filter_2.csv').iloc[:,0].to_list()
    afm.uns['cell_filter'] = 'custom_filter2'
    afm_filtered = afm[afm.obs_names.isin(cells),:].copy()
    afm_filtered.write('afm_filtered.h5ad')

    os.chdir('..')
