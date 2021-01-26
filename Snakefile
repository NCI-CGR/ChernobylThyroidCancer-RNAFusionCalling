# vim: ft=python
import os
import glob
import itertools
import pandas
from collections import defaultdict

configfile: 'fusion_config.yaml'
workdir: os.environ['PWD']
shell.executable('bash')

## This workflow takes a consensus approach with four fusion detection programs.
# Once candidate fusions have been identified, users can look closer at the output statistics to better identify real and false positives 
# (e.g. no spanning frags + no large anchor support = suspicious)


## ---- The parser may have to be customized for each run ---- ##
def parse_sampleID(filename):
    return filename.split('/')[-1].split('_')[0]

fastqs = glob.glob(config['fastq_dir'] + '/*fastq.gz')

d = defaultdict(list) # Pair the R1 and R2 files by sample ID
for key, value in itertools.groupby(fastqs, parse_sampleID):
    d[key] += list(value)

sampleIDs = d.keys()

def input_files(wildcards):
    return sorted(d[wildcards.sampleID])

def sorted_fusion(gene1, gene2):
    return '|'.join(sorted([gene1, gene2]))


include: 'ericscript_Snakefile'
include: 'mapsplice_Snakefile'
include: 'starfusion_Snakefile'
include: 'chimerascan_Snakefile'


# These rules run on the host node and are not submitted to the cluster.
localrules: all


rule all:
    input:
        'tables/all_scores.txt',
 

# Some tools require unzipped fastqs 
rule unzip_fq:
    input: input_files
    output: 
        'unzipped/{sampleID}_R1.fastq',
        'unzipped/{sampleID}_R2.fastq'
    run:
        shell('gunzip -c {input[0]} > {output[0]}')
        shell('gunzip -c {input[1]} > {output[1]}')


rule combine_scores:
    input: 
        eric = 'tables/all_es.txt',
        map = 'tables/all_ms.txt',
        star = 'tables/all_sf.txt',
        chim = 'tables/all_cs.txt'

    output: 'tables/all_scores.txt'
    params: meta = config['metadata']
    run:
        df_es = pandas.read_table(input.eric, sep='\t')
        df_ms = pandas.read_table(input.map, sep='\t')
        df_sf = pandas.read_table(input.star, sep='\t')
        df_cs = pandas.read_table(input.chim, sep='\t')

        cols = ['sampleID', 'fusion']
        df = df_es[cols + ['ericscript', 'es_score']].merge(df_ms[cols + ['mapsplice', 'ms_score']], on=cols, how='outer')
        df = df.merge(df_sf[cols + ['starfusion', 'sf_score']], on=cols, how='outer')
        df = df.merge(df_cs[cols + ['chimerascan', 'cs_score']], on=cols, how='outer')
        df.fillna(0, inplace=True)

        df['sum'] = df.ericscript + df.mapsplice + df.starfusion + df.chimerascan
        df['total_score'] = df.es_score + df.ms_score + df.sf_score + df.cs_score

        df.to_csv(output[0], sep='\t', index=False)
        

