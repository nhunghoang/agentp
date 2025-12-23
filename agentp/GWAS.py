'''
GWAS wrapper around REGENIE. 
'''

import pandas as pd
import sys
import os 
import subprocess
from multiprocessing import Pool

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
            msg = f'{MSG}GWAS ERROR: There are no remaining phenotypes to compute GWAS for.\n' 
            sys.exit(msg)
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

        ## copy first chrom with its headers 
        cmd = ['cp', f'{self.output_dir}/tmp/c1_step2_{phen}.regenie', phen_out]
        exit_code1 = os.system(' '.join(cmd))

        ## copy the rest without headers 
        chrom_in = f'{self.output_dir}/tmp/c' + '{2..22}' + f'_step2_{phen}.regenie' 
        cmd = ['tail', '-n', '+2', '-q', chrom_in, '>>', phen_out]
        exit_code2 = os.system(' '.join(cmd))

        return exit_code1 + exit_code2

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
        args = [(phen_list, phen_file, chrom) for chrom in range(1,23)]

        print(f'{MSG}Running GWAS...')
        pool = Pool(processes=num_threads)
        pool.starmap(self._call_regenie, args)
        print('')

        for phen in phen_list:
            exit_code = self.__concat_chrom_results(phen)
            if exit_code == 0:
                self.__reformat_gwas_output(phen)
            print(f'{TAB}Saved GWAS results for "{phen}"')

        if cleanup: 
            chroms = f'{self.output_dir}/tmp/c' + '{1..22}' + f'_step2_*.regenie' 
            os.system('rm ' + chroms)

            locos = f'{self.output_dir}/tmp/c' + '{1..22}' + f'_step1_*.loco'
            os.system('rm ' + locos)

            preds = f'{self.output_dir}/tmp/c' + '{1..22}' + f'_step1_pred.list'
            os.system('rm ' + preds)

            os.system(f'rm {phen_file}')
            
