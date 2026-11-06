'''
Build the small example model set (example/example_models) from the full AGENT-P models.

The example models only cover caudate, and only the SNPs in the example genotypes
(example/example_data/genotypes), so that example/run_example.py runs in minutes. They are a
demo subset, not for real analyses.

  PrediXcan / UTMOST / JTI: caudate weights (.db) and covariances, restricted to example SNPs
  FUSION: 300 caudate genes evenly spaced along the genome (weights .RDat, gene scores, .pos)
  LD reference (FUSION) and SNP lists: restricted to example SNPs

usage: python make_example_models.py [FULL_MODELS_DIR]   (default: the installed agentp/models)
Needs plink2 on PATH.
'''

import os
import sys
import gzip
import shutil
import sqlite3
import subprocess
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   ## repository root
FULL = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'agentp', 'models')
OUT = f'{ROOT}/example/example_models'
GENO = f'{ROOT}/example/example_data/genotypes'
N_FUS_GENES = 300

PREDICTDB = {'PDX': ('PrediXcan', 'en_Brain'), 'UTM': ('UTMOST_dzhou', 'UTMOST_Brain'), 'JTI': ('JTI', 'JTI_Brain')}
TISSUE = 'Caudate_basal_ganglia'

if os.path.exists(OUT):
    shutil.rmtree(OUT)
os.makedirs(OUT)
tmp = f'{OUT}/tmp'
os.makedirs(tmp)

## SNPs in the example genotypes
for c in range(1, 23):
    subprocess.run(['plink2', '--bgen', f'{GENO}/test_c{c}.bgen', 'ref-first', '--sample', f'{GENO}/test_c{c}.sample',
                    '--make-just-pvar', '--out', f'{tmp}/c{c}'], capture_output=True, check=True)
snps = set()
for c in range(1, 23):
    snps |= set(pd.read_table(f'{tmp}/c{c}.pvar', comment=None, usecols=['ID'], header=0,
                              names=['#CHROM', 'POS', 'ID', 'REF', 'ALT'])['ID'])
with open(f'{tmp}/snps.txt', 'w') as f:
    f.write('\n'.join(sorted(snps)) + '\n')
print(f'{len(snps)} example SNPs')

## PrediXcan / UTMOST / JTI: weights restricted to example SNPs, then covariances of the kept SNPs
for method, (folder, prefix) in PREDICTDB.items():
    name = f'{prefix}_{TISSUE}'
    os.makedirs(f'{OUT}/{folder}/weights')
    os.makedirs(f'{OUT}/{folder}/covariance')
    db = f'{OUT}/{folder}/weights/{name}.db'
    shutil.copy(f'{FULL}/{folder}/weights/{name}.db', db)
    with sqlite3.connect(db) as con:
        con.execute('CREATE TEMP TABLE keep (rsid TEXT PRIMARY KEY)')
        con.executemany('INSERT INTO keep VALUES (?)', [(s,) for s in snps])
        con.execute('DELETE FROM weights WHERE rsid NOT IN (SELECT rsid FROM keep)')
        con.execute('DELETE FROM extra WHERE gene NOT IN (SELECT DISTINCT gene FROM weights)')
        n_genes = con.execute('SELECT COUNT(DISTINCT gene) FROM weights').fetchone()[0]
        n_rows = con.execute('SELECT COUNT(*) FROM weights').fetchone()[0]
        con.commit()
    with sqlite3.connect(db) as con:
        con.execute('VACUUM')

    ## keep covariance lines whose gene and both SNPs remain (original separator kept)
    kept = 0
    with gzip.open(f'{FULL}/{folder}/covariance/{name}.txt.gz', 'rt') as fin, \
         gzip.open(f'{OUT}/{folder}/covariance/{name}.txt.gz', 'wt') as fout:
        fout.write(fin.readline())
        for line in fin:
            gene, r1, r2 = line.split(maxsplit=3)[:3]
            if (r1 in snps) and (r2 in snps):
                fout.write(line)
                kept += 1
    print(f'{method}: {n_genes} genes, {n_rows} weights, {kept} covariance entries')

## SNP lists used to slice genotypes for individual TWAS
for folder in ['snps_PDX_UTM_JTI', 'FUSION/snps']:
    os.makedirs(f'{OUT}/{folder}')
    for c in range(1, 23):
        ids = pd.read_csv(f'{FULL}/{folder}/c{c}.csv', header=None, dtype=str)[0]
        ids[ids.isin(snps)].to_csv(f'{OUT}/{folder}/c{c}.csv', index=False, header=False)

## FUSION: N_FUS_GENES caudate genes evenly spaced along the genome
name = f'GTExv8.EUR.Brain_{TISSUE}'
pos = pd.read_table(f'{FULL}/FUSION/WEIGHTS/{name}.pos')
pos = pos.sort_values(['CHR', 'P0']).reset_index(drop=True)
step = len(pos) / N_FUS_GENES
pos = pos.iloc[[int(i * step) for i in range(N_FUS_GENES)]]
os.makedirs(f'{OUT}/FUSION/WEIGHTS/{name}')
os.makedirs(f'{OUT}/FUSION/scores/caudate')
for wgt in pos['WGT']:
    gene = wgt.split('/')[1].split('.')[0]
    shutil.copy(f'{FULL}/FUSION/WEIGHTS/{wgt}', f'{OUT}/FUSION/WEIGHTS/{wgt}')
    shutil.copy(f'{FULL}/FUSION/scores/caudate/{gene}.score', f'{OUT}/FUSION/scores/caudate/{gene}.score')
pos.to_csv(f'{OUT}/FUSION/WEIGHTS/{name}.pos', sep='\t', index=False)
print(f'FUS: {len(pos)} genes')

## FUSION LD reference restricted to example SNPs
os.makedirs(f'{OUT}/FUSION/LDREF')
for c in range(1, 23):
    subprocess.run(['plink2', '--bfile', f'{FULL}/FUSION/LDREF/1000G.EUR.{c}', '--extract', f'{tmp}/snps.txt',
                    '--make-bed', '--out', f'{OUT}/FUSION/LDREF/1000G.EUR.{c}'], capture_output=True, check=True)
    os.remove(f'{OUT}/FUSION/LDREF/1000G.EUR.{c}.log')

shutil.rmtree(tmp)
with open(f'{OUT}/README.txt', 'w') as f:
    f.write('Small example model set made by scripts/make_example_models.py: caudate models of PrediXcan, '
            'UTMOST, JTI (restricted to the SNPs in example_data/genotypes) and FUSION (300 genes). '
            'For example runs only; use the full models for real analyses.\n')
print(f'Saved example models to {OUT}')
