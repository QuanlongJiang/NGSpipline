import glob,os
#import gseapy
from joblib import Parallel, delayed

def func(f,gmt):
    print f
    gmt_type = os.path.basename(gmt).split('.')[0]
    patho = os.path.basename(f).split('.')[0] + '/'
    if not os.path.exists(patho):
        os.system('mkdir {}'.format(patho))
    #cmd = 'gseapy prerank -r {} -g {} -o {} --no-plot --min-size 3 --max-size=1000 -p 2'.format(f,gmt,patho)
    cmd = 'gseapy prerank -r {} -g {} -o {} --graph 20 --permu-num 1000 --min-size 3 --max-size=1000 -p 4'.format(f,gmt,patho)
    print cmd
    os.system(cmd)
    #change output file name
    os.system('mv {}/gseapy.prerank.gene_sets.report.csv {}/{}.report.csv'.format(patho,patho, gmt_type))

def main():
    gmts=['Wb_GO.gmt']
    print gmts
    fs = glob.glob('*.rnk')
    for gmt in gmts:
        Parallel(n_jobs=10)(delayed(func)(f,gmt) for f in fs)

main()
