#supress warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import os
_stderr = sys.stderr
null = open(os.devnull,'wb')

import os
work_dir = 'NFYA'
import pycisTopic

#make a directory for to store the processed scRNA-seq data.
if not os.path.exists(os.path.join(work_dir, 'scATACc')):
    os.makedirs(os.path.join(work_dir, 'scATACc'))
tmp_dir = '/temp_work/ch250798/Scenic/tmp3'

fragments_dict = {'10x_NFYA': os.path.join(work_dir, 'data/atac_fragments.tsv.gz')}

import scanpy as sc
adata = sc.read_h5ad(os.path.join(work_dir, 'scRNA/cdata.h5ad'))
cell_data = adata.obs
cell_data['sample_id'] = '10x_NFYA'
cell_data['celltype'] = cell_data['celltype'].astype(str) # set data type of the celltype column to str, otherwise the export_pseudobulk function will complain.
del(adata)

import pyranges as pr
import requests
import pandas as pd
target_url='http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes'
chromsizes=pd.read_csv(target_url, sep='\t', header=None)
chromsizes.columns=['Chromosome', 'End']
chromsizes['Start']=[0]*chromsizes.shape[0]
chromsizes=chromsizes.loc[:,['Chromosome', 'Start', 'End']]
# Exceptionally in this case, to agree with CellRangerARC annotations
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].replace('v', '.') for x in range(len(chromsizes['Chromosome']))]
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].split('_')[1] if len(chromsizes['Chromosome'][x].split('_')) > 1 else chromsizes['Chromosome'][x] for x in range(len(chromsizes['Chromosome']))]
chromsizes=pr.PyRanges(chromsizes)


import pickle
cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATACc/cistopic_obj.pkl'), 'rb'))
from pycisTopic.cistopic_class import *
models=run_cgs_models(cistopic_obj,
                    n_topics=[2,4,10,16,32,48],
                    n_cpu=20,
                    n_iter=500,
                    random_state=555,
                    alpha=50,
                    alpha_by_topic=True,
                    eta=0.1,
                    eta_by_topic=False,
                    save_path=None,
                    _temp_dir = os.path.join(tmp_dir + 'ray_spill'))


import pickle
if not os.path.exists(os.path.join(work_dir, 'scATACc/models')):
    os.makedirs(os.path.join(work_dir, 'scATACc/models'))

pickle.dump(models,
            open(os.path.join(work_dir, 'scATACc/models/10x_NFYA_models_500_iter_LDA.pkl'), 'wb'))







