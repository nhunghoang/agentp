'''
Quality control on input genotype files.

- Nhung, Sept 2024
'''

import os
import re
import sys
import shutil
import subprocess

from .Static import MSG, TAB
from .Static import Scripts

'''
Print the sample/variant counts and filter results from a plink2 log.
Lines are matched by content rather than position, since line numbers shift
between plink2 versions and machines (e.g., wrapped RAM and file-path lines).
'''
def print_plink_summary(stdout):
    keys = ['loaded from', 'removed', 'remaining', 'Skipping']

    ## join lines that plink2 wrapped at 80 characters back into sentences
    sents, cur = [], ''
    for line in stdout.splitlines():
        if line.strip() == '':
            sents.append(cur); cur = ''
            continue
        cur = f'{cur} {line.strip()}'.strip()
        if cur.endswith(('.', ':')):
            sents.append(cur); cur = ''
    sents.append(cur)

    for sent in sents:
        if any(k in sent for k in keys):
            sent = re.sub(r' loaded from \S+', ' loaded.', sent)   ## drop temporary file paths
            print(f'{TAB}{sent}')
    print('')

class Genotypes:
    def __init__(self):
        self.file_type = None

class BGEN(Genotypes):
    def __init__(self, input_dir, input_pattern, subjects_file, genotype_dir):
        self.file_type = 'bgen'
        self.genotype_dir = genotype_dir

        ## check that the external tools are available before starting
        for tool in [Scripts.PLINK2, Scripts.BGENIX]:
            if shutil.which(tool) is None:
                msg = f'{MSG}QC ERROR: "{tool}" was not found on PATH (see Set Up in the README)'
                sys.exit(msg)

        ## slice input genotypes to the project cohort (tmp/c*_0), the input of apply_qc
        self.__get_cohort(input_dir, input_pattern, subjects_file)

    '''
    NOTE that the sliced genotypes are not necessarily in the project subject order
    '''
    def __get_cohort(self, input_dir, input_pattern, subjects_file):
        for chrom in range(1, 23):

            cfile = f'{self.genotype_dir}/tmp/c{chrom}_0.pgen'
            if os.path.exists(cfile): continue

            cname = input_pattern.replace('*', str(chrom))
            for ext in ['bgen', 'sample']:
                if not os.path.exists(f'{input_dir}/{cname}.{ext}'):
                    msg = f'{MSG}QC ERROR: Input genotype file "{input_dir}/{cname}.{ext}" does not exist'
                    sys.exit(msg)

            print(f'{MSG}Extracting project cohort from chromosome {chrom}...')
            cmd0 = [Scripts.PLINK2, \
                    '--bgen', f'{input_dir}/{cname}.bgen', 'ref-first', \
                    '--sample', f'{input_dir}/{cname}.sample', \
                    '--snps-only', \
                    '--keep', subjects_file, \
                    '--make-pgen', \
                    '--out', f'{self.genotype_dir}/tmp/c{chrom}_0']
            try:
                res0 = subprocess.run(cmd0, capture_output=True, text=True, check=True)
                print_plink_summary(res0.stdout)

            except subprocess.CalledProcessError as e:
                print(f'{e.stdout}\n{e.stderr}')
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
                print_plink_summary(res1.stdout)

            except subprocess.CalledProcessError as e:
                print(f'{e.stdout}\n{e.stderr}')
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
                print_plink_summary(res2.stdout)
                del2 = os.system(f'rm {self.genotype_dir}/tmp/c{chrom}_1.p*')

            except subprocess.CalledProcessError as e:
                print(f'{e.stdout}\n{e.stderr}')
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
                print_plink_summary(res3.stdout)

            except subprocess.CalledProcessError as e:
                print(f'{e.stdout}\n{e.stderr}')
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
