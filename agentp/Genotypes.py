'''
Quality control on input genotype files. 

- Nhung, Sept 2024
'''

import os
import sys
import subprocess

from .Static import MSG, TAB
from .Static import Scripts

class Genotypes: 
    def __init__(self): 
        self.file_type = None

class BGEN(Genotypes): 
    def __init__(self, input_dir, input_pattern, subjects_file, genotype_dir): 
        self.file_type = 'bgen'
        self.genotype_dir = genotype_dir
        #self.__get_cohort(input_dir, input_pattern, subjects_file)

    '''
    NOTE that the sliced genotypes are not necessarily in the project subject order
    '''
    def __get_cohort(self, input_dir, input_pattern, subjects_file): 
        for chrom in range(1, 23):

            cfile = f'{self.genotype_dir}/tmp/c{chrom}_0.pgen'
            if os.path.exists(cfile): continue

            cname = input_pattern.replace('*', str(chrom))
            cmd0 = [Scripts.PLINK2, \
                    '--bgen', f'{input_dir}/{cname}.bgen', 'ref-first', \
                    '--sample', f'{input_dir}/{cname}.sample', \
                    '--snps-only', \
                    '--keep', subjects_file, \
                    '--make-pgen', \
                    '--out', f'{self.genotype_dir}/tmp/c{chrom}_0']
            try:
                os.system(' '.join(cmd0))
                #res0 = subprocess.run(cmd0, capture_output=True, text=True, check=True)
                #out0 = res0.stdout.splitlines()
                #for i in range(3,10): 
                #    print(f'{TAB}{out1[i]}')
                #print(f'\n{TAB}{out1[15]}')
                #print(f'{TAB}{out1[40]}\n')

            except subprocess.CalledProcessError as e:
                print(f'{e.stderr}')
                sys.exit()

    def apply_qc(self, geno_hardcall=0.1, snp_missing=0.05, maf=0.01, hwe=0.00001): 
        for chrom in range(1, 23):
            tag = f'{"#" * 40} CHROMOSOME {chrom} {"#" * 40}'
            print(f'\n{tag}')
        
            ## genotype hardcall 
            print(f'{MSG}Genotype Hard-call...')
            cmd1 = [Scripts.PLINK2, \
                    '--pfile', f'{self.genotype_dir}/tmp/c{chrom}_0', \
                    '--hard-call-threshold', str(geno_hardcall), \
                    '--make-pgen', \
                    '--out', f'{self.genotype_dir}/tmp/c{chrom}_1']
            try:
                res1 = subprocess.run(cmd1, capture_output=True, text=True, check=True)
                out1 = res1.stdout.splitlines()
                for i in range(3,8): 
                    print(f'{TAB}{out1[i]}')
                print(f'\n{TAB}{out1[13]}')
                print(f'{TAB}{out1[14]}')
                print(f'{TAB}{out1[15]}')
                print(f'{TAB}{out1[16]}\n')

            except subprocess.CalledProcessError as e:
                print(f'{e.stderr}')
                sys.exit()
            
            ## snp missingness filter, then MAF 
            print(f'{MSG}SNP Missingness and MAF Filters...')
            cmd2 = [Scripts.PLINK2, \
                    '--pfile', f'{self.genotype_dir}/tmp/c{chrom}_1', \
                    '--geno', str(snp_missing), \
                    '--maf', str(maf), \
                    '--make-pgen', \
                    '--out', f'{self.genotype_dir}/tmp/c{chrom}_2']
            try:
                res2 = subprocess.run(cmd2, capture_output=True, text=True, check=True)
                out2 = res2.stdout.splitlines()
                for i in range(3,9): 
                    print(f'{TAB}{out2[i]}')
                print(f'\n{TAB}{out2[20]}')
                print(f'{TAB}{out2[21]}')
                print(f'{TAB}{out2[23]}\n')
                del2 = os.system(f'rm {self.genotype_dir}/tmp/c{chrom}_1.p*')

            except subprocess.CalledProcessError as e:
                print(f'{e.stderr}')
                sys.exit()

            ## remove snp duplicates (by id, keep first), then HWE filter 
            print(f'{MSG}SNP Duplicate and HWE Filters...')
            cmd3 = [Scripts.PLINK2, \
                    '--pfile', f'{self.genotype_dir}/tmp/c{chrom}_2', \
                    '--rm-dup', 'exclude-mismatch', \
                    '--hwe', str(hwe), \
                    '--export', 'bgen-1.2', 'bits=8', \
                    '--out', f'{self.genotype_dir}/c{chrom}']
            try:
                res3 = subprocess.run(cmd3, capture_output=True, text=True, check=True)
                out3 = res3.stdout.splitlines()
                for i in range(3,9): 
                    print(f'{TAB}{out3[i]}')
                print(f'\n{TAB}{out3[19]}')
                print(f'{TAB}{out3[21]}')
                print(f'{TAB}{out3[22]}\n')

            except subprocess.CalledProcessError as e:
                print(f'{e.stderr}')
                sys.exit()

            ## bgen index 
            cmd4 = [Scripts.BGENIX, \
                    '-g', f'{self.genotype_dir}/c{chrom}.bgen', \
                    '-index', \
                    '-clobber']            
            try:
                res4 = subprocess.run(cmd4, capture_output=True, text=True, check=True)
                del4 = os.system(f'rm {self.genotype_dir}/tmp/c{chrom}_2.p*')
            except subprocess.CalledProcessError as e:
                print(f'{e.stderr}')
                sys.exit()

        print('#' * len(tag))

class VCF(Genotypes):
    pass 
