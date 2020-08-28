import glob,os
from joblib import Parallel,delayed 

def func(f):
    out_dir=os.path.basename(f).split('.')[0]
    motif='homer/data/knownTFs/all.motifs'
    cmd='findMotifs.pl {} worm {} -start -1000 -end 1000 -p 8 -bg refseq_ids/background -mknown {} -mcheck {}'.format(f,out_dir,motif,motif)
    print cmd
    os.system(cmd)

fs=glob.glob('refseq_ids/*txt')

Parallel(n_jobs=len(fs))(delayed(func)(f) for f in fs)
