'''
TWAS class
'''

import pandas as pd
import os
import sys
import json
import subprocess
from .Static import Pool   ## fork-based pool, see Static.Pool

from .Static import TAB, MSG
from .Static import Paths, Scripts, Models

## Note: imported in IndividualTWAS class if needed:
#from GREX import GREX

'''
Raised when a TWAS run fails. Used instead of sys.exit() in code that runs
inside Pool workers: sys.exit() kills the worker without returning a result,
which makes Pool.starmap wait forever.
'''
class TWASError(Exception):
    pass

class TWAS:
    def __init__(self, project, twas_approach, grex_region):
        self.project = project
        self.approach = twas_approach
        self.grex_region = grex_region
        Models.check(twas_approach, grex_region)
        self.weight_file = Models.WEIGHT(twas_approach, grex_region)

    '''
    Run TWAS on one phenotype, returning an error message instead of raising,
    so that one failed phenotype does not stop the others in the Pool.
    Note: Making this a pseudo-private method (Pool can't call private functions)
    '''
    def _try_phenotype_twas(self, phen, overwrite):
        try:
            self._run_phenotype_twas(phen, overwrite)
            return None
        except Exception as e:
            return f'{type(e).__name__}: {e}'

    '''
    Multi-process run on a list of phenotypes; failures are reported after all finish.
    '''
    def _run_pool(self, phens, overwrite, num_threads):
        params = [(phen, overwrite) for phen in phens]
        with Pool(processes=num_threads) as pool:
            errs = pool.starmap(self._try_phenotype_twas, params)

        failed = [(phen, err) for phen, err in zip(phens, errs) if err is not None]
        if len(failed) > 0:
            msg = f'{MSG}TWAS ERROR: {len(failed)} of {len(phens)} phenotypes failed ' \
                  f'using {self.approach} {self.grex_region} models'
            for phen, err in failed:
                msg += f'\n{TAB}{phen}: {err}'
            raise TWASError(msg)

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

        ## the GWAS file and parameters used to make the LD-Score file, saved next to it,
        ## so that a cached file is only reused for the same input
        source = {'gwas_file': os.path.abspath(gwas_file), \
                  'size': os.path.getsize(gwas_file), \
                  'mtime': int(os.path.getmtime(gwas_file)), \
                  'params': list(gwas_params)}
        source_file = f'{output_gwas}.source.json'

        if os.path.exists(f'{output_gwas}.sumstats.gz'):
            if not os.path.exists(source_file):
                msg = f'{MSG}WARNING: Loaded existing LD-Score GWAS file for phenotype "{phen}", ' \
                      f'but the GWAS file and parameters it was made from were not recorded' \
                      f'\n{TAB}Delete {output_gwas}.sumstats.gz to rebuild it from "{gwas_file}"'
                print(msg)
                return f'{output_gwas}.sumstats.gz'
            with open(source_file) as f:
                cached = json.load(f)
            if cached == source:
                msg = f'{MSG}Loaded existing LD-Score GWAS file for phenotype "{phen}"'
                print(msg)
                return f'{output_gwas}.sumstats.gz'
            msg = f'{MSG}The existing LD-Score GWAS file for phenotype "{phen}" was made from a ' \
                  f'different GWAS file or parameters ({cached["gwas_file"]}); rebuilding it'
            print(msg)

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

        with open(source_file, 'w') as f:
            json.dump(source, f, indent=2)
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
                msg = f'FUSION failed on chrom {chrom} ({phen}, {self.approach} {self.grex_region} models)\n{e.stderr}'
                raise TWASError(msg)

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

        ## the MetaXcan scripts refuse to overwrite an existing output (and exit 0), so remove it first
        if os.path.exists(output_file): os.remove(output_file)
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
            last = ''
            for line in res.stderr:
                last = line
                if line[:4] != 'INFO':
                    errs = res.stderr.read()
                    if len(errs) == 0: err = line.split('-')[1].strip()
                    else: err = errs.splitlines()[-1].strip()
                    res.kill()
                    msg = f'({phen}, {self.approach} {self.grex_region} models): {err}'
                    raise TWASError(msg)

                if '%' not in line: continue
                info = line.split('-')[1].strip()
                print(f'{TAB}{info}')

            ## SPrediXcan also reports errors as INFO lines (e.g. "INFO - Unexpected error: ...")
            res.wait()
            if (res.returncode != 0) or (not os.path.exists(output_file)):
                err = last.split(' - ', 1)[-1].strip()
                msg = f'({phen}, {self.approach} {self.grex_region} models): ' \
                      f'SPrediXcan exited with code {res.returncode}: {err}'
                raise TWASError(msg)

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
            self._run_pool(phens, overwrite, num_threads)
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

        ## the MetaXcan scripts refuse to overwrite an existing output (and exit 0), so remove it first
        if os.path.exists(output_twas_file): os.remove(output_twas_file)
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
            msg = f'({phen}, {self.approach} {self.grex_region} models): ' \
                  f'PrediXcanAssociation exited with code {e.returncode} (see its output above)'
            raise TWASError(msg)
        ## PrediXcanAssociation can also log an error and exit 0, so check for the output file
        if not os.path.exists(output_twas_file):
            msg = f'({phen}, {self.approach} {self.grex_region} models): ' \
                  f'PrediXcanAssociation did not create {output_twas_file} (see its output above)'
            raise TWASError(msg)

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
            self._run_pool(phens, overwrite, num_threads)
        print('')

########################################################################################

'''
Run TWAS for a list of phenotypes with every model set (all approaches and GREx
regions, including models added with Project.add_model), then combine the results
per phenotype with an ACAT meta-analysis (see Meta.collect_twas).

project:     an instance of the Project class
phens:       string or list of strings, phenotype names
twas_type:   'stwas' (summary TWAS) or 'itwas' (individual TWAS)
gwas_files:  for summary TWAS, dict of {phenotype: GWAS file}
gwas_params: for summary TWAS, dict of add_gwas() parsing parameters used for every
             GWAS file (e.g., {'sep': ',', 'snp_col': 'rsid', ...}); default is REGENIE format
approaches / regions: lists of approach and region abbreviations to include (default: all)
overwrite:   rerun TWAS that already exist
num_threads: number of phenotypes run in parallel per model set
combine:     run collect_twas() for each phenotype after all TWAS finish

A model set that fails (e.g., missing model files) is reported and skipped, and the
remaining model sets still run. Returns a dict of {phenotype: combined table} if combine.
'''
def run_all_twas(project, phens, twas_type='stwas', gwas_files=None, gwas_params=None, \
                 approaches=None, regions=None, overwrite=False, num_threads=1, combine=True):
    from .Meta import collect_twas

    phens = [phens] if not type(phens) == list else phens
    gwas_params = {} if gwas_params is None else gwas_params
    if twas_type not in ['stwas', 'itwas']:
        msg = f'{MSG}TWAS ERROR: twas_type must be "stwas" or "itwas"'
        sys.exit(msg)

    ## check inputs before starting any TWAS
    if twas_type == 'stwas':
        missing = [p for p in phens if (gwas_files is None) or (p not in gwas_files)]
        if len(missing) > 0:
            msg = f'{MSG}sTWAS ERROR: No GWAS file given for phenotype(s): {", ".join(missing)}'
            sys.exit(msg)
        for p in phens:
            if not os.path.exists(gwas_files[p]):
                msg = f'{MSG}sTWAS ERROR: Input gwas file "{gwas_files[p]}" does not exist'
                sys.exit(msg)

    model_sets = [(app, reg) for (app, reg) in Models.model_sets() \
                  if (approaches is None or app in approaches) and (regions is None or reg in regions)]
    if len(model_sets) == 0:
        msg = f'{MSG}TWAS ERROR: No model sets match the given approaches and regions'
        sys.exit(msg)

    print(f'{MSG}Running {"summary" if twas_type == "stwas" else "individual"} TWAS on ' \
          f'{len(phens)} phenotype(s) with {len(model_sets)} model sets...')
    failed = []
    for k, (app, reg) in enumerate(model_sets):
        print(f'{MSG}[{k+1}/{len(model_sets)}] {app} {reg} models')
        try:
            if twas_type == 'stwas':
                twas = SummaryTWAS(project, app, reg)
                for p in phens:
                    twas.add_gwas(p, gwas_files[p], **gwas_params)
            else:
                twas = IndividualTWAS(project, app, reg)
                twas.predict_grex()
            twas.run_twas(phens, overwrite=overwrite, num_threads=num_threads)
        except TWASError as e:
            failed.append((app, reg, str(e).strip()))
            print(f'{TAB}WARNING: {app} {reg} models failed (skipped); see the summary at the end')

    if len(failed) > 0:
        msg = f'{MSG}WARNING: {len(failed)} of {len(model_sets)} model sets failed; ' \
              f'their results are missing from the combined tables:'
        for app, reg, err in failed:
            lines = [l.strip() for l in err.splitlines() if l.strip() and not l.startswith('> TWAS ERROR')]
            msg += f'\n{TAB}{app} {reg}: ' + f'\n{TAB}{TAB}'.join(lines)
        print(msg)

    if not combine:
        return None
    return {p: collect_twas(project, p, twas_type, approaches=approaches, regions=regions) \
            for p in phens}


