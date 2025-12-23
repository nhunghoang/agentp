'''
TWAS class
'''

import pandas as pd 
import os 
import sys
import subprocess
from multiprocessing import Pool

from .Static import TAB, MSG
from .Static import Paths, Scripts, Models

## Note: imported in IndividualTWAS class if needed:
#from GREX import GREX

class TWAS:
    def __init__(self, project, twas_approach, grex_region): 
        self.project = project 
        self.approach = twas_approach
        self.grex_region = grex_region
        self.weight_file = Models.WEIGHT(twas_approach, grex_region) 

########################################################################################

'''
For parsing GWAS files, as needed for SummaryTWAS. 
'''
def get_gwas_params(sep, snp_col, ref_col, alt_col, \
                    beta_col, beta_sign_col, beta_se_col, \
                    or_col, zscore_col, pval_col): 

    params = [] 
    if sep: params.extend(['--separator', sep])
    if snp_col: params.extend(['--snp_column', snp_col])
    if ref_col: params.extend(['--effect_allele_column', ref_col])
    if alt_col: params.extend(['--non_effect_allele_column', alt_col])
    if beta_col: params.extend(['--beta_column', beta_col])
    if beta_sign_col: params.extend(['--beta_sign_column', beta_sign_col])
    if beta_se_col: params.extend(['--se_column', beta_se_col])
    if or_col: params.extend(['--or_column', or_col])
    if zscore_col: params.extend(['--zscore_column', zscore_col])
    if pval_col: params.extend(['--pvalue_column', pval_col])

    ## if no args, set to regenie args 
    if len(params) == 0: 
        params = ['--separator', ',', \
                  '--snp_column', 'SNP', \
                  '--effect_allele_column', 'ALLELE0', \
                  '--non_effect_allele_column', 'ALLELE1', \
                  '--beta_column', 'BETA', \
                  '--se_column', 'SE', \
                  '--pvalue_column', 'PVAL']
    return params
            
########################################################################################

## TWAS using GWAS summary statistics 
class SummaryTWAS(TWAS):
    def __init__(self, project, twas_approach, grex_region): 
        super().__init__(project, twas_approach, grex_region)
        self.framework = 'summary'

        paths = project.set_summary_twas_paths(twas_approach, grex_region)
        self.output_dir = paths['output']

        self.gwas_inputs = {} ## k: phen, v: path to gwas file 
        self.gwas_params = {} ## k: phen, v: list of params for parsing gwas file

        if twas_approach == 'FUS':
            self.ldref_file = Models.FUS_LDREF
            self.ldsc_dir = paths['ldsc']
        else:
            self.covariance_file = Models.COVARIANCE(twas_approach, grex_region)

        ## load existing TWAS 
        self.results = {} ## k: phen, v: output file path 
        for tfile in os.listdir(self.output_dir): 
            phen = tfile.split('.')[0]
            self.results[phen] = f'{self.output_dir}/{tfile}'

    '''
    Save input GWAS in LD-Score format (currently just for FUSION). 
    TODO: add option for signed-sumstats that aren't BETA
    '''
    def __save_gwas_ldscore(self, phen, gwas_file, gwas_params): 
        output_gwas = f'{self.ldsc_dir}/{phen}' ## .sumstats.gz 
        if os.path.exists(f'{output_gwas}.sumstats.gz'):
            msg = f'{MSG}Loaded existing LD-Score GWAS file for phenotype "{phen}"'
            print(msg)
            return f'{output_gwas}.sumstats.gz'

        print(f'{MSG}Saving "{phen}" GWAS file in LD-Score format...')
        params = {gwas_params[i][2:]: gwas_params[i+1] for \
                  i in range(0, len(gwas_params), 2)}

        cmd = ['python', '-u', Scripts.FUS_LDSC] 
        cmd.extend(['--sumstats', gwas_file, \
                    '--delimiter', params['separator'], \
                    '--snp', params['snp_column'], \
                    '--a1', params['effect_allele_column'], \
                    '--a2', params['non_effect_allele_column'], \
                    '--p', params['pvalue_column'], \
                    '--info-min', '0.1', \
                    '--signed-sumstats', f'{params["beta_column"]},0', \
                    '--out', output_gwas])

        sp_params = {'stdout': subprocess.PIPE, 'stderr': subprocess.PIPE, \
                     'text': True, 'bufsize': 1, 'universal_newlines': True}
        res = subprocess.Popen(cmd, **sp_params)

        flag = False
        for line in res.stdout: 
            first = line.split(' ')[0]
            if first == 'Read': flag = True
            elif first == 'Writing': flag = False

            if flag: print(f'{TAB}{line.strip()}')

        res.wait()
        if res.returncode != 0: 
            err = res.stderr.read().splitlines()[-1].strip()
            msg = f'{TAB}{err}\n'
            sys.exit(msg)

        print(f'{TAB}Done.')
        return f'{output_gwas}.sumstats.gz'

    '''
    FUSION approach for sTWAS. 
    Runs separately per chrom, then concats as one phen output file. 
    '''
    def __run_fusion(self, phen): 
        output_chrom = f'{self.output_dir}/tmp/{phen}_c*.dat'        
        for chrom in range(1,23):
            out_file = output_chrom.replace('*', str(chrom))
            cmd = ['Rscript', Scripts.FUS_TWAS]
            cmd.extend(['--sumstats', self.gwas_inputs[phen], \
                        '--weights', self.weight_file, \
                        '--weights_dir', Models.WEIGHT_DIR('FUS'), \
                        '--ref_ld_chr', self.ldref_file, \
                        '--chr', str(chrom), \
                        '--out', out_file])
            ## allow all models regardless of % missing SNPs; will get filtered in later steps 
            cmd.extend(['--max_impute', '1.0'])

            try:
                res = subprocess.run(cmd, capture_output=True, text=True, check=True)
                msg = f'{TAB}({phen}, {self.approach} {self.grex_region} models): chrom {chrom} done'
                print(msg)
            except subprocess.CalledProcessError as e:
                print(f'{e.stderr}')
                sys.exit()

        ## concat and reformat  
        cnames = {'ID': 'gene', 'TWAS.Z': 'zscore', 'TWAS.P': 'pvalue', 'MODEL': 'best_model'}
        dtypes = {'ID': str, 'TWAS.Z': float, 'TWAS.P': float, 'MODEL': str}

        cfiles = [] 
        cfile_ = [output_chrom.replace('*', str(c)) for c in range(1,23)]
        cfile_.append(output_chrom.replace('*', '6') + '.MHC')
        for cfile in cfile_:
            df = pd.read_table(cfile, usecols=cnames.keys())
            try:
                mask = df['TWAS.P'].str.contains('NA')
                df = df.loc[~mask].astype(dtypes)
            except: 
                pass

            df = df.rename(columns=cnames)
            cfiles.append(df)

        pdf = pd.concat(cfiles)
        pdf['gene'] = pdf['gene'].apply(lambda x: x.split('.')[0])
        output_phen = f'{self.output_dir}/{phen}.csv'
        pdf.to_csv(output_phen, sep='\t', index=False)
        print(f'\n{TAB}({phen}, {self.approach} {self.grex_region} models): finished')

    '''
    Run sTWAS on one phenotype. 
    Note: Making this a pseudo-private method (Pool can't call private functions)
    TODO: decide how to communicate TWAS logs
    '''
    def _run_phenotype_twas(self, phen, overwrite): 
        output_file = f'{self.output_dir}/{phen}.csv'        
        if os.path.exists(output_file) and (not overwrite): 
            msg = f'{MSG}Skipping sTWAS on "{phen}" using {self.approach} {self.grex_region} models - already exists'
            print(msg)
            return

        if phen not in self.gwas_inputs.keys(): 
            msg = f'{MSG}sTWAS ERROR: No input GWAS file was added for phenotype "{phen}"'
            print(msg)
            return

        print(f'{MSG}Running sTWAS on "{phen}" using {self.approach} {self.grex_region} models...')
        if self.approach == 'FUS':
            res = self.__run_fusion(phen)
        else: 
            cmd = ['python', '-u', Scripts.SUMMARY_TWAS]
            cmd.extend(['--model_db_path', self.weight_file, \
                        '--covariance', self.covariance_file, \
                        '--gwas_file', self.gwas_inputs[phen], \
                        '--output_file', output_file])
            cmd.extend(self.gwas_params[phen])

            sp_params = {'stdout': subprocess.PIPE, 'stderr': subprocess.PIPE, \
                         'text': True, 'bufsize': 1, 'universal_newlines': True}
            res = subprocess.Popen(cmd, **sp_params)
            for line in res.stderr:
                if line[:4] != 'INFO': 
                    errs = res.stderr.read()
                    if len(errs) == 0: err = line.split('-')[1].strip()
                    else: err = errs.splitlines()[-1].strip()
                    msg = f'{TAB}{err}\n'
                    sys.exit(msg)

                if '%' not in line: continue 
                info = line.split('-')[1].strip()
                print(f'{TAB}{info}')

            res.wait()

    '''
    Track path and params of every input gwas.
    '''
    def add_gwas(self, phenotype, gwas_file, sep=None, \
                 snp_col=None, ref_col=None, alt_col=None, \
                 beta_col=None, beta_sign_col=None, beta_se_col=None, \
                 or_col=None, zscore_col=None, pval_col=None): 

        if not os.path.exists(gwas_file): 
            msg = f'{MSG}sTWAS ERROR: Input gwas file "{gwas_file}" does not exist'
            sys.exit(msg)

        params = get_gwas_params(sep, snp_col, ref_col, alt_col, \
                                 beta_col, beta_sign_col, beta_se_col, \
                                 or_col, zscore_col, pval_col)
        self.gwas_params[phenotype] = params

        if self.approach == 'FUS': 
            ldsc_file = self.__save_gwas_ldscore(phenotype, gwas_file, params) 
            self.gwas_inputs[phenotype] = ldsc_file
        else:
            self.gwas_inputs[phenotype] = gwas_file 
                    
    '''
    Multi-process run on list of phenotypes. 
    '''
    def run_twas(self, phens, overwrite=False, num_threads=1): 
        if not type(phens) == list: 
            self._run_phenotype_twas(phens, overwrite)
        else:
            params = [(phen, overwrite) for phen in phens]
            pool = Pool(processes=num_threads)
            pool.starmap(self._run_phenotype_twas, params)
        print('')

########################################################################################

## TWAS using individual (genotype) data 
class IndividualTWAS(TWAS):
    def __init__(self, project, twas_approach, grex_region): 
        super().__init__(project, twas_approach, grex_region)
        self.framework = 'individual'

        paths = project.set_individual_twas_paths(twas_approach, grex_region)
        self.twas_genotypes_dir = paths['twas_gen']
        self.grex_file = paths['grex'] 
        self.output_dir = paths['output']

        self.subjects_file = project.get_subjects_file()
        self.covariates_file = project.get_covariates_file()
        self.genotypes_dir = project.get_genotypes_dir()

        try:
            covar_cols = pd.read_table(self.covariates_file, nrows=0)
            self.covar_cols = covar_cols.columns.values[2:]
        except: 
            self.covar_cols = None

        ## load existing TWAS 
        self.results = {} ## k: phen, v: output file path 
        for tfile in os.listdir(self.output_dir): 
            phen = tfile.split('.')[0]
            self.results[phen] = f'{self.output_dir}/{tfile}'

    '''
    Run iTWAS on one phenotype. 
    TODO: decide how to communicate TWAS logs
    '''
    def _run_phenotype_twas(self, phen, overwrite):
        output_twas_file = f'{self.output_dir}/{phen}.csv'        
        if os.path.exists(output_twas_file) and (not overwrite): 
            msg = f'{MSG}Skipping iTWAS on "{phen}" using {self.approach} {self.grex_region} models - already exists'
            print(msg) 
            return

        msg = f'{MSG}Running iTWAS on "{phen}" using {self.approach} {self.grex_region} models...'
        print(msg)

        input_phen_file = self.project.get_phenotype_file(phen) 
        cmd = ['python', '-u', Scripts.GENOTYPE_TWAS]  
        cmd.extend(['--hdf5_expression_file', self.grex_file, \
                    '--input_phenos_file', input_phen_file, \
                    '--input_phenos_column', phen, \
                    '--output', output_twas_file, \
                    '--covariates_file', self.covariates_file, \
                    '--covariates'])
        cmd.extend(self.covar_cols)

        try:
            res = subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f'{e.stderr}')
            sys.exit()

    '''
    Predict GREX (FUS: gene scores into grex | JTI/PDX/UTM: dosages into grex)
    TODO: decide how to communicate GREX logs
    '''
    def predict_grex(self, overwrite=False, cleanup=False): 
        if os.path.exists(self.grex_file) and (not overwrite):
            msg = f'{MSG}Skipping {self.approach} {self.grex_region} models - GREX already exists'
            print(msg)
            return

        ## init grex object 
        from .GREX import GREX
        grex = GREX(self.approach, self.grex_region, self.twas_genotypes_dir, self.grex_file)

        ## slice genotypes to TWAS variants (generate .bed/.dosage files) 
        if self.approach == 'FUS': 
            exists = [os.path.exists(f'{self.twas_genotypes_dir}/c{i}.bed') for i in range(1,23)]
        else: 
            exists = [os.path.exists(f'{self.twas_genotypes_dir}/c{i}.dosage.txt.gz') for i in range(1,23)]
        if sum(exists) != 22: 
            msg = f'{MSG}Filtering genotypes to TWAS SNPs...'
            print(msg)
            grex.slice_bgen(self.genotypes_dir, cleanup) 

        ## apply model weights to generate GREX 
        msg = f'{MSG}Predicting GREX based on {self.approach} {self.grex_region} models...'
        print(msg)
        grex.apply_weights(self.project.get_subjects(), self.project.get_subjects_file())

    '''
    Multi-process run on list of phenotypes. 
    '''
    def run_twas(self, phens, overwrite=False, num_threads=1): 
        if self.covar_cols is None: 
            try:
                covar_cols = pd.read_table(self.covariates_file, nrows=0)
                self.covar_cols = covar_cols.columns.values[2:]
            except: 
                msg = f'{MSG}iTWAS ERROR: Add subject-level covariates to project first'
                sys.exit(msg) 

        if not type(phens) == list: 
            self._run_phenotype_twas(phens, overwrite)
        else:
            params = [(phen, overwrite) for phen in phens]
            pool = Pool(processes=num_threads)
            pool.starmap(self._run_phenotype_twas, params)
        print('')



