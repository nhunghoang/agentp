'''
GWAS wrapper around REGENIE. 
'''

import pandas as pd
import sys
import os 
import subprocess
import shutil
import glob
from .Static import Pool   ## fork-based pool, see Static.Pool

from .Static import MSG, TAB
from .Static import Scripts

class GWAS: 
    def __init__(self, project): 
        self.project = project

        paths = project.set_gwas_paths()
        self.output_dir = paths['output']

        self.covariates_file = project.get_covariates_file() 
        self.genotypes_dir = project.get_genotypes_dir()

    '''
    Combine all input phenotypes into one temporary file. 
    Phenotypes with existing GWAS are excluded. 
    '''
    def __gather_phenotypes(self, phenotypes, overwrite): 

        phens = {} ## k: phen name, v: phen file 
        excl_phens = [] 
        none_phens = [] 
        for phen in phenotypes: 
            
            ## check if phen file exists 
            phen_file = self.project.get_phenotype_file(phen) 
            if not os.path.exists(phen_file):
                none_phens.append(phen)
                phen_file = None

            ## check if phen GWAS exists 
            if phen_file: 
                gwas_file = f'{self.output_dir}/{phen}.regenie'

                ## add to phen list if no GWAS yet (or overwrite=True)
                if not os.path.exists(gwas_file): 
                    phens[phen] = phen_file 
                elif overwrite: 
                    phens[phen] = phen_file 
                else: 
                    excl_phens.append(phen) 

        if len(none_phens) > 0:
            msg = f'{MSG}GWAS WARNING: The following phenotypes do not exist in this project:'
            msg += f'\n{TAB}{", ".join(none_phens)}'
            print(msg)
        if len(excl_phens) > 0: 
            msg = f'{MSG}GWAS WARNING: The following phenotypes already have GWAS and will be ignored:'
            msg += f'\n{TAB}{", ".join(excl_phens)}'
            print(msg) 

        ## combine valid phenotypes into one temp file 
        if len(phens) == 0:
            msg = f'{MSG}GWAS: There are no remaining phenotypes to compute GWAS for.\n'
            print(msg)
            return [], None
        else:
            phen_data = [] 
            for phen, pfile in phens.items():
                pdata = pd.read_table(pfile, sep='\t', index_col=['FID', 'IID']) 
                phen_data.append(pdata) 

            tmp_file = f'{self.output_dir}/tmp/gwas_phens.csv' 
            ptable = pd.concat(phen_data, axis=1) 
            ptable.to_csv(tmp_file, sep='\t')

            msg = f'{MSG}GWAS will be computed for the following phenotypes:'
            msg += f'\n{TAB}{", ".join(phens.keys())}'
            print(msg)
            return ptable.columns.values, tmp_file

    '''
    Call REGENIE for the list of phens for one chrom.
    Note: Making this a pseudo-private method (Pool can't call private functions)
    '''
    def _call_regenie(self, phen_list, phen_file, chrom, bsize=1000, threads=1):
        bgen_file = f'{self.genotypes_dir}/c{chrom}.bgen'
        samp_file = f'{self.genotypes_dir}/c{chrom}.sample'
        covs_file = self.covariates_file

        temp_path = f'{self.output_dir}/tmp/mem_c{chrom}'
        out1_path = f'{self.output_dir}/tmp/c{chrom}_step1'
         
        pred_path = f'{out1_path}_pred.list'
        out2_path = f'{self.output_dir}/tmp/c{chrom}_step2' 

        ## regenie step 1
        cmd1 = [Scripts.GWAS, '--step', '1', \
                            '--bsize', str(bsize), \
                            '--threads', str(threads), \
                            '--lowmem']

        phens = ','.join(phen_list)
        cmd1.extend(['--bgen', bgen_file, \
                     '--covarFile', covs_file, \
                     '--phenoFile', phen_file, \
                     '--phenoColList', phens, \
                     '--lowmem-prefix', temp_path, \
                     '--out', out1_path])
        try:
            res = subprocess.run(cmd1, capture_output=True, text=True, check=True)
            print(f'{TAB}[CHR {chrom}] finished REGENIE step 1/2')
        except subprocess.CalledProcessError as e:
            print(f'{e.stderr}')
            sys.exit()

        ## regenie step 2
        cmd2 = [Scripts.GWAS, '--step', '2', \
                            '--ref-first', \
                            '--bsize', str(bsize), \
                            '--threads', str(threads)]

        cmd2.extend(['--bgen', bgen_file, \
                     '--sample', samp_file, \
                     '--covarFile', covs_file, \
                     '--phenoFile', phen_file, \
                     '--phenoColList', phens, \
                     '--pred', pred_path, \
                     '--out', out2_path])

        try:
            res = subprocess.run(cmd2, capture_output=True, text=True, check=True)
            print(f'{TAB}[CHR {chrom}] finished REGENIE step 2/2')
        except subprocess.CalledProcessError as e:
            print(f'{e.stderr}')
            sys.exit()

    '''
    For one phen, combine all its REGENIE chrom files. 
    '''
    def __concat_chrom_results(self, phen):
        phen_out = f'{self.output_dir}/{phen}.regenie'

        ## done in Python, not with shell brace expansion ({2..22}), which /bin/sh
        ## does not support on Ubuntu/Debian (dash): there, only chromosome 1 was kept
        chrom_files = [f'{self.output_dir}/tmp/c{c}_step2_{phen}.regenie' for c in range(1, 23)]
        missing = [f for f in chrom_files if not os.path.exists(f)]
        if len(missing) > 0:
            print(f'{TAB}GWAS ERROR: REGENIE results missing for "{phen}": {", ".join(missing)}')
            return 1

        ## first chromosome with its header, the rest without
        with open(phen_out, 'w') as out:
            for k, cfile in enumerate(chrom_files):
                with open(cfile) as inp:
                    header = inp.readline()
                    if k == 0: out.write(header)
                    shutil.copyfileobj(inp, out)
        return 0

    '''
    For one phen, reformat its REGENIE file. 
    '''
    def __reformat_gwas_output(self, phen): 
        phen_file = f'{self.output_dir}/{phen}.regenie' 
        df = pd.read_table(phen_file, sep=' ')

        df['PVAL'] = 10 ** (-df['LOG10P']) 
        df = df.rename(columns={'ID': 'SNP'})
        df = df.drop(columns=['EXTRA'])
        
        df.to_csv(phen_file, index=False, sep=',')

    '''
    Apply REGENIE to a list of phenotypes. 
    '''
    def run(self, phens, overwrite=False, cleanup=True, num_threads=1): 
        if not type(phens) == list:
            phens = [phens]

        phen_list, phen_file = self.__gather_phenotypes(phens, overwrite)
        if len(phen_list) == 0:
            return

        args =[(phen_list, phen_file, chrom) for chrom in range(1,23)]

        print(f'{MSG}Running GWAS...')
        pool = Pool(processes=num_threads)
        pool.starmap(self._call_regenie, args)
        print('')

        failed = []
        for phen in phen_list:
            exit_code = self.__concat_chrom_results(phen)
            if exit_code == 0:
                self.__reformat_gwas_output(phen)
                print(f'{TAB}Saved GWAS results for "{phen}"')
            else:
                failed.append(phen)

        if cleanup:
            for c in range(1, 23):
                for pattern in [f'c{c}_step2_*.regenie', f'c{c}_step1_*.loco', f'c{c}_step1_pred.list']:
                    for f in glob.glob(f'{self.output_dir}/tmp/{pattern}'):
                        os.remove(f)
            if os.path.exists(phen_file): os.remove(phen_file)

        if len(failed) > 0:
            msg = f'{MSG}GWAS ERROR: no results were saved for: {", ".join(failed)} (see the errors above)'
            sys.exit(msg)
            
