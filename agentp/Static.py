'''
File containing static paths and variables. 

- Nhung, Sept 2024
'''

import os
import sys
## importlib.resources instead of pkg_resources (removed in setuptools >= 81)
from importlib.resources import files, as_file
import multiprocessing as _mp

def Pool(processes=None):
    '''
    Process pool that uses "fork" where available (Linux, macOS), as on Linux by default.
    macOS defaults to "spawn", which re-runs the user's script in every worker process
    (unless it has an `if __name__ == '__main__':` guard) and does not carry over settings
    such as set_models_dir(); e.g. the workers then rewrite the genotype files that the
    main process is reading.
    '''
    method = 'fork' if 'fork' in _mp.get_all_start_methods() else None
    return _mp.get_context(method).Pool(processes=processes)

class Paths:

    with as_file(files('agentp') / 'external') as p: SOFTWARE = str(p)
    with as_file(files('agentp') / 'models') as p: MODELS = str(p)
    ## gene positions for plots (GTEx v8 genes, GRCh38; columns chr, start, end, gene_id, gene_name)
    with as_file(files('agentp') / 'data' / 'gene_positions.csv') as p: GENE_POSITIONS = str(p)

class Scripts:

    PLINK2 = 'plink2'
    BGENIX = 'bgenix'
    GWAS = 'regenie'

    with as_file(files('agentp') / 'external' / 'MetaXcan' / 'SPrediXcan.py') as p: SUMMARY_TWAS = str(p)
    with as_file(files('agentp') / 'external' / 'MetaXcan' / 'PrediXcanAssociation.py') as p: GENOTYPE_TWAS = str(p)
    with as_file(files('agentp') / 'external' / 'MetaXcan' / 'Predict.py') as p: GREX_PREDICTION = str(p)

    with as_file(files('agentp') / 'external' / 'ldsc' / 'munge_sumstats.py') as p: FUS_LDSC = str(p)
    with as_file(files('agentp') / 'external' / 'fusion_twas' / 'FUSION.assoc_test.R') as p: FUS_TWAS = str(p)

class Models: 

    rlabel = ['dlpfc', 'ant-cingulate', 'amygdala', 'hippocampus', \
              'caudate', 'putamen', 'nuc-accumbens', 'cerebellar-hemi']

    tissue = ['Frontal_Cortex_BA9', 'Anterior_cingulate_cortex_BA24', \
              'Amygdala', 'Hippocampus', 'Caudate_basal_ganglia', \
              'Putamen_basal_ganglia', 'Nucleus_accumbens_basal_ganglia', \
              'Cerebellar_Hemisphere']

    REGION = {r:t for r,t in zip(rlabel, tissue)}

    method = ['PDX', 'FUS', 'UTM', 'JTI']
    folder = ['PrediXcan', 'FUSION', 'UTMOST_dzhou', 'JTI']
    prefix = ['en_Brain', 'GTExv8.EUR.Brain', 'UTMOST_Brain', 'JTI_Brain']

    FOLD_MAP = {m:f for m,f in zip(method, folder)}
    PREF_MAP = {m:p for m,p in zip(method, prefix)}

    ## user-added PredictDB-format model sets (see Project.add_model)
    ## k: (method, region), v: {'weight': .db, 'covariance': .txt.gz, 'snps': SNP list}
    CUSTOM = {}

    def is_custom(method, region):
        return (method, region) in Models.CUSTOM

    '''
    All available (method, region) model sets: built-in, then user-added.
    '''
    def model_sets():
        sets = [(m, r) for m in Models.method for r in Models.rlabel]
        return sets + [k for k in Models.CUSTOM if k not in sets]

    def check(method, region):
        if (method, region) in Models.model_sets():
            return
        msg = f'{MSG}MODEL ERROR: No "{method}" models for region "{region}"' \
              f'\n{TAB}Built-in approaches: {", ".join(Models.method)}' \
              f'\n{TAB}Built-in regions: {", ".join(Models.rlabel)}'
        if len(Models.CUSTOM) > 0:
            msg += f'\n{TAB}Added model sets: {", ".join(f"{m} {r}" for m, r in Models.CUSTOM)}'
        msg += f'\n{TAB}Use Project.add_model() to add PredictDB-format models'
        sys.exit(msg)

    def WEIGHT_DIR(method):
        if method == 'FUS': 
            return f'{Paths.MODELS}/FUSION/WEIGHTS'
        else: 
            fold = Models.FOLD_MAP[method] 
            return f'{Paths.MODELS}/{fold}/weights'

    def WEIGHT(method, region):
        if Models.is_custom(method, region):
            return Models.CUSTOM[(method, region)]['weight']
        fold = Models.FOLD_MAP[method]
        pref = Models.PREF_MAP[method]
        rmod = Models.REGION[region]
        if method == 'FUS': 
            return f'{Paths.MODELS}/FUSION/WEIGHTS/{pref}_{rmod}.pos'
        else: 
            return f'{Paths.MODELS}/{fold}/weights/{pref}_{rmod}.db'

    def COVARIANCE(method, region):
        if Models.is_custom(method, region):
            return Models.CUSTOM[(method, region)]['covariance']
        fold = Models.FOLD_MAP[method]
        pref = Models.PREF_MAP[method]
        rmod = Models.REGION[region]
        return f'{Paths.MODELS}/{fold}/covariance/{pref}_{rmod}.txt.gz'

    def SNPS(method): 
        if method == 'FUS':
            return f'{Paths.MODELS}/FUSION/snps'
        else: 
            return f'{Paths.MODELS}/snps_PDX_UTM_JTI'

    '''
    SNP list used to slice genotypes for one chromosome. User-added models use
    one list of all their SNPs (plink2 --extract ignores IDs not on the chromosome).
    '''
    def SNP_FILE(method, region, chrom):
        if Models.is_custom(method, region):
            return Models.CUSTOM[(method, region)]['snps']
        return f'{Models.SNPS(method)}/c{chrom}.csv'

    def FUS_SCORES(region):
        return f'{Paths.MODELS}/FUSION/scores/{region}'

    FUS_LDREF = f'{Paths.MODELS}/FUSION/LDREF/1000G.EUR.'

'''
Use a different models folder with the same layout as agentp/models (e.g., the small
example models in example/example_models). Call before creating TWAS objects.
'''
def set_models_dir(models_dir):
    if not os.path.isdir(models_dir):
        msg = f'{MSG}MODEL ERROR: Models folder "{models_dir}" does not exist'
        sys.exit(msg)
    Paths.MODELS = os.path.abspath(models_dir)
    Models.FUS_LDREF = f'{Paths.MODELS}/FUSION/LDREF/1000G.EUR.'
    print(f'{MSG}Using models in {Paths.MODELS}')

## string elements for logging
TAB = ' ' * 4
MSG = '\n> '

