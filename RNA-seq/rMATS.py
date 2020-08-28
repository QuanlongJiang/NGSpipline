import glob,os
from joblib import Parallel, delayed
import collections
import itertools

def prepare():
    fs = glob.glob('.mapping/*/*out.bam')  ###
    fs = sorted(fs)
    d = collections.defaultdict(list)
    for f in fs:
        g = '_'.join(os.path.basename(f).split('_')[:2])
        d[g].append(os.path.abspath(f))
    for g in d:
        fo = '{}.txt'.format(g) ###
        with open(fo,'w') as o:
            print fo
            for line in d[g]:
                o.write(line+',')

def call_splicing(b1,b2):
    tag  = os.path.basename(b1).split('.')[0] + '_' + os.path.basename(b2).split('.')[0]
    out_dir = '{}'.format(tag)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    gtf = 'gencode.v31.annotation.gtf'
    rmast = "rMATS.4.0.1/rMATS-turbo-Linux-UCS4/rmats.py"
    cmd = 'python {rmast} --gtf {gtf} --b1 {b1} --b2 {b2} --od {out} -t paired --readLength 150 --nthread 5 --tstat 5'.format(rmast=rmast, gtf=gtf, b1=b1,b2=b2, out=out_dir)
    print cmd
    os.system(cmd)

### prepare 
#prepare()

### run rmats
#bs_comb = itertools.combinations(bs,2)
bs_comb = [['YA.txt','OA.txt'],['OA.txt','PD.txt'],['PD0.txt','PD16.txt']]

Parallel(n_jobs=3)(delayed(call_splicing)(b1,b2) for (b1,b2) in bs_comb)

