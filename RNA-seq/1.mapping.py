import glob,os
from joblib import Parallel, delayed
import pandas as pd

STAR_index = '/home/jiangquanlong/jql_han/0.Genome/gencode.v31/STAR'
chrom_size = '/home/jiangquanlong/jql_han/0.Genome/gencode.v31/GRCh38.p12.genome.fa.fai'

def fastqc(f):
    out_dir = 'fastqc/'
    #out_dir = 'fastqc_cutadapt/'
    if not os.path.exists(out_dir):
        os.system('mkdir {}'.format(out_dir))
    cmd='fastqc -t 4 --noextract -o {} -f fastq {}'.format(out_dir, f)
    print cmd
    os.system(cmd)

def mapping_func(prefix, fqs):
    mapping_dir = 'mapping/'
    out_dir = mapping_dir + prefix + '/'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    prefix = out_dir + prefix + '_'
    bam = prefix + 'Aligned.out.bam'
    cmd = "STAR --runMode alignReads --quantMode GeneCounts --runThreadN 6 --genomeDir {genomeDir} --genomeLoad LoadAndKeep --readFilesIn {fqs} --readMatesLengthsIn Equal --readFilesCommand 'zcat -1' --clip5pNbases 15 --clip3pNbases 0 --outFileNamePrefix {prefix} --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --outSAMreadID Number --outFilterMultimapNmax 10 --outFilterMismatchNmax 9 --outFilterMismatchNoverLmax 0.1 --outFilterMatchNminOverLread 0.3 --outFilterScoreMinOverLread 0.3 --alignEndsType EndToEnd".format(genomeDir=STAR_index, fqs=' '.join(fqs), prefix=prefix)    
    print cmd
    os.system(cmd)
    sort_cmd = 'samtools sort -@ 4 -o {} {}'.format(bam, bam)
    os.system(sort_cmd)
    index_cmd = 'samtools index {}'.format(bam)
    os.system(index_cmd)
    STAR_pileup(bam)
    
def mapping_main():
    f1s = glob.glob('cutadapt/*_1.fastq.gz')
    fq_dict = {}
    for f1 in f1s:
        prefix = os.path.basename(f1).split('_1')[0]
        f2 = f1.replace('_1.fastq.gz','_2.fastq.gz')
        fq_dict[prefix] = [f1, f2]
    Parallel(n_jobs=4)(delayed(mapping_func)(prefix, fq_dict[prefix]) for prefix in fq_dict)

def bg2bw(bg):
    bg_sort = bg.replace('.bg','.sort.bg')
    cmd_sort = 'LC_COLLATE=C sort -k 1,1 -k 2,2n {} > {}'.format(bg, bg_sort)
    os.system(cmd_sort)
    bw = bg_sort.replace('.bg','.bw')
    cmd1 = "bedGraphToBigWig {bedgraph} {chromsize} {bigwig}".format(bedgraph=bg_sort, chromsize=chrom_size, bigwig=bw)
    os.system(cmd1)
    cmd2 = 'rm {}'.format(' '.join([bg, bg_sort]))
    os.system(cmd2)
    
def STAR_pileup(bam):
    prefix = bam.replace('Aligned.out.bam', '')
    cmd = "STAR --runMode inputAlignmentsFromBAM --runThreadN 6 --genomeDir {genomeDir} --genomeLoad LoadAndKeep --outFileNamePrefix {prefix} --outWigType bedGraph --outWigStrand Stranded --outWigNorm RPM --inputBAMfile {bam}".format(genomeDir=STAR_index, prefix=prefix, bam=bam)
    os.system(cmd)
    bg = prefix + 'Signal.UniqueMultiple.str2.out.bg'
    bg2bw(bg)

def mapping_ratio():
    fs = glob.glob('mapping/*/*_Log.final.out')
    data = {}
    for f in fs:
        prefix = os.path.basename(f).split('_Log')[0]
        data[prefix] = {}
        for line in open(f):
            line = line.strip().split('\t')
            if line[0].startswith('Number of input reads'):
                data[prefix]['Number of input reads'] = float(line[-1])
            if line[0].startswith('Uniquely mapped reads number'):
                data[prefix]['Uniquely mapped reads number'] = float(line[-1])
            if line[0].startswith('Uniquely mapped reads %'):
                data[prefix]['Uniquely mapped reads %'] = line[-1]
            if line[0].startswith('Number of reads mapped to multiple loci'):
                data[prefix]['Number of reads mapped to multiple loci'] = float(line[-1])
            if line[0].startswith('% of reads mapped to multiple loci'):
                data[prefix]['% of reads mapped to multiple loci'] = line[-1]
    data = pd.DataFrame.from_dict(data, orient='columns')
    cols = sorted(data.columns)
    data = data[cols]
    data = data.loc[['Number of input reads', 'Uniquely mapped reads number', 'Uniquely mapped reads %', 'Number of reads mapped to multiple loci', '% of reads mapped to multiple loci']]
    data.loc['Number of total mapped reads'] = data.loc['Uniquely mapped reads number'] + data.loc['Number of reads mapped to multiple loci']
    data.loc['Total mapped reads %'] = data.loc['Number of total mapped reads']/data.loc['Number of input reads']
    data.loc['Total mapped reads %'] = data.loc['Total mapped reads %'].apply(lambda x: format(x,'.2%'))
    data.to_csv('Mapping_Ratio.txt',sep='\t')
    
### fastqc
#fs = glob.glob('fastq/*/*.fq.gz')
#Parallel(n_jobs=6)(delayed(fastqc)(f) for f in fs)
###

mapping_main()

mapping_ratio()

