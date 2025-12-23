'''
File containing static paths and variables. 

- Nhung, Sept 2024
'''

import os
import pkg_resources as pkgr

class Paths: 

    SOFTWARE = pkgr.resource_filename('agentp', 'external')
    MODELS = pkgr.resource_filename('agentp', 'models') 

class Scripts: 

    PLINK2 = 'plink2'
    BGENIX = 'bgenix'
    GWAS = 'regenie'

    SUMMARY_TWAS = pkgr.resource_filename('agentp', 'external/MetaXcan/SPrediXcan.py') 
    GENOTYPE_TWAS = pkgr.resource_filename('agentp', 'external/MetaXcan/PrediXcanAssociation.py')
    GREX_PREDICTION = pkgr.resource_filename('agentp', 'external/MetaXcan/Predict.py')

    FUS_LDSC = pkgr.resource_filename('agentp', 'external/ldsc/munge_sumstats.py')
    FUS_TWAS = pkgr.resource_filename('agentp', 'external/fusion_twas/FUSION.assoc_test.R')

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
    
    def WEIGHT_DIR(method):
        if method == 'FUS': 
            return f'{Paths.MODELS}/FUSION/WEIGHTS'
        else: 
            fold = Models.FOLD_MAP[method] 
            return f'{Paths.MODELS}/{fold}/weights'

    def WEIGHT(method, region): 
        fold = Models.FOLD_MAP[method] 
        pref = Models.PREF_MAP[method]
        rmod = Models.REGION[region]
        if method == 'FUS': 
            return f'{Paths.MODELS}/FUSION/WEIGHTS/{pref}_{rmod}.pos'
        else: 
            return f'{Paths.MODELS}/{fold}/weights/{pref}_{rmod}.db'

    def COVARIANCE(method, region): 
        fold = Models.FOLD_MAP[method] 
        pref = Models.PREF_MAP[method]
        rmod = Models.REGION[region]
        return f'{Paths.MODELS}/{fold}/covariance/{pref}_{rmod}.txt.gz'

    def SNPS(method): 
        if method == 'FUS':
            return f'{Paths.MODELS}/FUSION/snps'
        else: 
            return f'{Paths.MODELS}/snps_PDX_UTM_JTI'

    def FUS_SCORES(region):
        return f'{Paths.MODELS}/FUSION/scores/{region}'

    FUS_LDREF = f'{Paths.MODELS}/FUSION/LDREF/1000G.EUR.'

## string elements for logging 
TAB = ' ' * 4
MSG = '\n> '

