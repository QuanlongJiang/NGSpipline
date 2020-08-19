import glob,os
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats.mstats import zscore

featureCount = '/home/jiangquanlong/jql_han/packages/subread-1.6.3-Linux-x86_64/bin/featureCounts'
gtf = '/home/jiangquanlong/jql_han/0.Genome/gencode.v31/gencode.v31.annotation.gtf'
anno_f = '/home/jiangquanlong/jql_han/0.Genome/gencode.v31/GenomeAnnotation.txt'

def run_featureCount():
    fs = glob.glob('mapping/*/*_Aligned.out.bam')
    print len(fs)
    fs = ' '.join(fs)
    fo = 'featureCount/reads.txt'
    cmd = '{tool} -T 30 -s 0 -p -M -O -C --fraction -t exon  -g gene_id -a {gtf} -o {fo} {bams}'.format(tool=featureCount, gtf=gtf, bams=fs, fo=fo)    #p for pair-end -C 
    print cmd
    os.system(cmd)

def getSymbolID():
    id_dict = {}
    with open(anno_f) as r:
        next(r)
        for line in r:
            line = line.strip().split('\t')
            id_dict[line[0]] = line[1]
    return id_dict

def getAllReads():
    f = 'Mapping_Ratio.txt'
    df = pd.read_table(f,index_col=0)
    #df = df.loc['Number of input reads']   #### 
    df = df.loc['Number of total mapped reads']  ###
    N_dict = dict(df)
    return N_dict

def generate_FPKM():
    id_dict = getSymbolID()
    N_dict = getAllReads()
    f = 'featureCount/reads.txt'
    df = pd.read_table(f,header=1)
    df.index = df.iloc[:,0]
    len_dict = dict(df.iloc[:,5])
    df = df.iloc[:,6:]
    samples = df.columns
    for sample in samples:
        N = N_dict[os.path.basename(sample).split('_Aligned')[0]]
        N = float(N)
        fo = 'featureCount/'+ os.path.basename(sample).split('_Aligned')[0] + '.txt'
        s_dict = dict(df[sample])
        data = {};rpk_dict = {}
        for gene in s_dict:
            reads = s_dict[gene] ; len_gene = len_dict[gene]
            fpkm = 1000000 * float( reads ) / N / (float( len_gene )/1000)
            rpm = 1000000 * float( reads ) / N 
            data[gene] = [int(reads), rpm, fpkm]
            rpk = float(reads) / (float(len_gene)/1000)
            rpk_dict[gene] = rpk
        rpk_m = sum(rpk_dict.values()) / 1000000
        for gene in data:
            rpk = rpk_dict[gene]
            tpm = rpk / rpk_m
            data[gene].append(tpm)
        with open(fo,'w')  as o:
            first = ['symbol','read','RPM','FPKM','TPM']
            o.write('\t'.join(first)+'\n')
            for gene in data:
                if gene in id_dict:
                    gene_id = gene + '|' + id_dict[gene]
                else:
                    gene_id = gene + '|' + gene
                o.write( gene_id+ '\t' + '\t'.join(map(str,data[gene]))+'\n')

def choose_quantitative_method(f, t):
    sample = os.path.basename(f).split('.')[0]
    lines = open(f).readlines(); lines = map(lambda l:l.strip().split('\t'),lines)[1:]
    d = {}
    idx_d = {'count':0, 'rpm':1,'fpkm':2,'tpm':3}
    idx = idx_d[t]
    for line in lines:
        fpkm = line[-2]
        tpm = line[-1]
        count = line[-3]
        tmp = [count, fpkm, tpm]
        d[line[0]] = tmp[idx]
    return sample,d

def summary(minReads, quantitative_method):
    fs = glob.glob('featureCount/*.txt')
    fs = [f for f in fs if 'read' not in f]
    eles = Parallel(n_jobs=1)(delayed(choose_quantitative_method)(f,'count') for f in fs)
    count_data = {}
    for ele in eles:
        sample = ele[0]
        count_data[sample] = ele[1]
    count_data = pd.DataFrame.from_dict(count_data,orient='columns',dtype=float)
    count_data.to_csv('raw_counts.txt')
    count_filter = filter_low_expression(count_data, minReads)
    count_filter.to_csv('counts.txt',sep='\t')
#####
    eles = Parallel(n_jobs=1)(delayed(choose_quantitative_method)(f, quantitative_method) for f in fs)
    fpkm_data = {}
    for ele in eles:
        sample = ele[0]
        fpkm_data[sample] = ele[1]
    fpkm_data = pd.DataFrame.from_dict(fpkm_data,orient='columns',dtype=float)
    fpkm_data.to_csv('raw_{}.txt'.format(quantitative_method))
    fpkm_filter = fpkm_data.loc[count_filter.index]
    fpkm_filter.to_csv('{}.txt'.format(quantitative_method,sep='\t')

def filter_low_expression(df,minReads):
    passed_idxs = []
    for idx, row in df.iterrows():
        n = len(row[row>minReads])
        if n >= 1:   ### at least one sample
            passed_idxs.append(idx)
    df = df.loc[passed_idxs]
    return df

def get_gene_type():
    lines = map(lambda l: l.strip().split('\t'), open(anno_f).readlines())[1:]
    gene_type_dict = {'PC':[], 'nonPC':[]}
    for line in lines:
        gene_id = '|'.join(line[:2])
        if line[2] == 'protein_coding':
            gene_type_dict['PC'].append(gene_id)
        else:
            gene_type_dict['nonPC'].append(gene_id)
    return gene_type_dict

def seperate_pc():
    gene_type_dict = get_gene_type()
    f = 'TPM.txt'
    df = pd.read_table(f,index_col=0)
    for gene_type in gene_type_dict:
        genes = df.index & set(gene_type_dict[gene_type])
        print gene_type, len(genes)
        tmp_df = df.loc[genes]
        tmp_fo = f.replace('.txt','.{}.txt'.format(gene_type))
        tmp_df.to_csv(tmp_fo, sep='\t')


run_featureCount()
generate_FPKM() ### calculate RPM, FPKM, TPM,
summary(minReads=5, quantitative_method = 'TPM') #min reads in at least one sample
#seperate_pc()
