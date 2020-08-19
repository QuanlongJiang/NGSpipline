import glob,os
import pandas as pd
import numpy as np
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import multipletests
from joblib import Parallel,delayed
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import seaborn as sns
sns.set_style( "white" )

plt.rcParams[ "font.size" ] = 4.0
plt.rcParams[ "figure.dpi" ] = 100
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['axes.labelsize'] = 6 
plt.rcParams[ "figure.figsize" ] = ( 4*0.7,2.75*0.8 )
plt.rcParams[ "font.serif" ] = 'Arial'

def getdbg(f_term):   #default bg
    fbg = f_term
    termname = fbg.split('/')[-1].split('.')[0]  ### term name
    bg = []
    dict_term = {}
    for line in open(fbg):
        line = line.strip().split('\t')
        gene = line[0];go_id = line[1]
        bg.append( gene )
        if line[-1] not in dict_term:
            dict_term[go_id] = []
            dict_term[go_id].append(gene)
        else:
            dict_term[go_id].append(gene)
    go_ids = dict_term.keys()
    for go_id in go_ids:
        if len(dict_term[go_id])<3:
            del dict_term[go_id]
    bg = set(bg)
    return dict_term,bg,termname

def getcbg(f_term):   #coustom bg   
    start = 0 #strat from which line
    cbg = []
    lines = map(lambda l:l.strip().split('\t'),open(bg_file).readlines())[start:]
    for line in lines:
        cbg.append(line[0])
    cbg = set(cbg)
    fbg = f_term
    termname = os.path.basename(fbg).split('.')[0]  ### term name
    termname = '.'.join(termname.split('_'))
    dict_term = {}
    for line in open(fbg):
        line = line.strip().split('\t')
        gene = line[0]; go_id = line[1]
        if gene in cbg: # the gene in goterm should in the coustom bg
            if go_id not in dict_term:
                dict_term[go_id] = []
                dict_term[go_id].append(gene)
            else:
                dict_term[go_id].append(gene)   
        else:
            pass 
    go_ids = dict_term.keys()
    for go_id in go_ids:
        if len(dict_term[go_id])<3:  ###if member of term < 3, delete
            del dict_term[go_id]
    return dict_term,cbg,termname

def getin(bg,pathin,suffix):
    start = 0 ###
    f_list = glob.glob(pathin+suffix) ###
    dic_lstin = {}
    for fin in f_list:
        lst = []
        tag = os.path.basename(fin).split('.')[0]
        lines = map(lambda l:l.strip().split('\t'),open(fin).readlines())[start:]
        for line in lines:
            if line[0] in bg:       ### the gene in should in the 
                lst.append(line[0])
            else:
                pass
        lst = list(set(lst))
        dic_lstin[tag] = lst
    return dic_lstin


def para(tag,dic_lstin,term_dict,bg,M,patho,termname):
    dic_result = {}
    fo = patho + tag + '_%s.txt'%termname
    lstin = dic_lstin[tag]
    N = len(lstin)  #list total
    n_forpadj = 0
    for term in term_dict:
        termlst = term_dict[term]
        n = len(termlst)    #term 
        intersect_genes=[i for i in lstin if i in termlst]
        k = len(intersect_genes)
        p = hypergeom.sf(k-1,M, n, N)
        if k!=0:
            if p < 0.05:
                dic_result[term] = [term,str(k),str((float(k)/N)*100),p,','.join(intersect_genes),str(N),str(n),str(M)]
            else:
                n_forpadj+=1
        else:
            pass
    if dic_result != {}:
        df = pd.DataFrame.from_dict(dic_result,orient='index')
        df.columns = ['Term','Count','%','PValue','Genes','List Total','Pop Hits','Pop Total']
        plst = df['PValue'].tolist()
        nump = len(plst)
        plst = plst + [0.1]*n_forpadj
        padjlst=multipletests(plst, method='bonferroni')[1][:nump]            
        df['Bonferroni'] = padjlst
        df = df.sort_values(by='Bonferroni')
        df.to_csv(fo,sep='\t',index=False)
        plot_func(fo)
    else:
        pass

def plot_func(f, numbercut=10, xlabel="-log10(p)", title="Enriched GO Terms",width=0.6):
#    if 'up' in f:
#        title  = '{} Up-regulated Genes {}'.format(os.path.basename(f).split('_')[0],title)
#    else:
#        title  = '{} Down-regulated Genes {}'.format(os.path.basename(f).split('_')[0],title)
    title = os.path.basename(f).split('.')[0] +'\n'+ title
    fo = f.replace('.txt','.pdf')
    df = pd.read_table(f,header=0,sep='\t')
    #df = df[df['Bonferroni']<0.05] ###
    df=df[df['PValue']<0.01]
    if df.empty:
        return 
    df = df.sort_values(by='PValue')
    df['-log10(pvalue)'] = df[['PValue']].apply(lambda x:-np.log10(x),1) ###
    for idx in df.index:
        if df['-log10(pvalue)'][idx] > 50:
            df['-log10(pvalue)'][idx] = 50
    df = df[df['Count']>=2]
    idxs = df['Term']; idxs = [idx.split('|')[1] for idx in idxs]
    df.index = idxs
    s = df['-log10(pvalue)']
    if s.shape[0] > numbercut:
        s = s[:numbercut]
    s = s[::-1]
    ns = s.shape[0]*2
    xa = np.arange(0,ns,2); xa = [x+0.3 for x in xa]
    xb = np.arange(1,ns,2)
    #xb = [x-0.3 for x in xb]
    fig,ax = plt.subplots( )
    ax.barh(xa,s.values,width,color= '#3498db')
    for i,t in enumerate(s.index):
        ax.annotate(t,(0+0.1,xb[i]),fontsize=5)
    plt.xlabel( xlabel )
    plt.title(title)
    ax.set_yticks([])
    ax.set_ylim([0,ns])
    sns.despine()
    plt.tight_layout()
    plt.savefig( fo )

def main(ifbg,pathin,suffix):   # ifbg = 'd' # default; coustom
    for f_term in term_list:
        if ifbg == 'd':
            term_dict,bg,termname = getdbg(f_term)
            patho = pathin + '%s_defaultbg/'%termname
        else:
            term_dict,bg,termname = getcbg(f_term)
            patho = pathin + '%s_custombg/'%termname
        if not os.path.exists(patho):
            os.system('mkdir %s'%patho)
        M = len(bg)
        dic_lstin = getin(bg,pathin,suffix)
        Parallel(n_jobs = 8)(delayed(para)(tag,dic_lstin,term_dict,bg,M,patho,termname) for tag in dic_lstin)

if __name__ == '__main__':
    bg_file = '../background.txt'
    term_list = glob.glob('../../database/Wb.GO') 
    pathin = ''
    suffix = '*.txt'
    #main('c',pathin,suffix) #c: custom background, d: default background
    fs = glob.glob('Wb_custombg/*txt')
    [plot_func(f) for f in fs]
