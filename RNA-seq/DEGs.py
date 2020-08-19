import glob,os
import pandas as pd
import numpy as np
import scipy.stats as ss
import collections
from statsmodels.sandbox.stats.multicomp import multipletests


def get_degs():
    result = {}
    df = pd.read_table('../TPM_log_qn.txt',index_col=0)
    df = df+0.1
    group = collections.defaultdict(list)
    for sample in df.columns:
        g = sample.split('_rep')[0]
        group[g].append(sample)

    genes = df.index
    for g_8 in ['BWM_8','Coel_8','Int_8']:
        tmp = {}
        g_1 = g_8.replace('_8','_1')
        for gene in genes:
            s = df.loc[gene]
            v1s = s[group[g_8]].tolist()
            v2s = s[group[g_1]].tolist()
            s,p = ss.ttest_ind(v1s, v2s)
            fc = np.mean(v1s)/np.mean(v2s)
            tmp[gene] =v2s+v1s+ [np.mean(v2s), np.mean(v1s),fc, p]
        tmp = pd.DataFrame.from_dict(tmp, orient='index')
        cols = []
        for idx in ['_1','_2','_3']:
            cols.append(g_1 + idx)
        for idx in ['_1','_2','_3']:
            cols.append(g_8 + idx)
        tmp.columns = cols + ['mean_{}'.format(g_1),'mean_{}'.format(g_8),'foldchange','pvalue']
        ps = tmp['pvalue']
        ps_adj = multipletests(ps, method='bonferroni')[1]
        tmp['pvalue_adj'] = ps_adj
        result[g_8] = tmp

    for g_8 in result.keys():
        fo = '{}_ttest.txt'.format(g_8.replace('_8',''))
        result[g_8].to_csv(fo,sep='\t')


def seperate_degs():
    fs = glob.glob('*_ttest.txt')
    for f in fs:
        sample = os.path.basename(f).split('_')[0]
        df = pd.read_table(f, index_col=0)
        ups = df[(df['pvalue']<0.05)&(df['foldchange']>2)].index; ups = pd.Series([g.split('|')[0] for g in ups])###
        ups.to_csv('degs/ensembl_ids/{}_upregulation.txt'.format(sample), sep='\t', index=False)
        downs = df[(df['pvalue']<0.05)&(df['foldchange']<0.5)].index; downs = pd.Series([g.split('|')[0] for g in downs])###
        downs.to_csv('degs/ensembl_ids/{}_downregulation.txt'.format(sample), sep='\t', index=False)



get_degs()
#seperate_degs()


