'''
Class for handling GREX input genotypes and output predictions 
for the different TWAS gene models. 
'''

import pandas as pd
import numpy as np 
import h5py
import bgen_reader as bgr

from multiprocessing import Pool
from time import time
import gzip
import glob
import os
import subprocess

from .Static import TAB
from .Static import Scripts, Models

class GREX: 

    def __init__(self, twas_approach, grex_region, twas_genotypes_dir, output_grex_file): 
        self.twas_approach = twas_approach 
        self.grex_region = grex_region
        self.twas_genotypes_dir = twas_genotypes_dir
        self.out_path = output_grex_file.replace('.hdf5', '')

    '''
    Note: Making this a pseudo-private method (Pool can't call private functions)
    '''
    def _compute_dosage_subset(self, version, chrom, snps, gens, output_file): 
        dos_lines = []
        for snp, gen in zip(snps.itertuples(), gens):
            rsid = snp.rsid
            post = str(snp.pos)
            [ref, alt] = snp.allele_ids.split(',')
            dose = np.argmax(gen, axis=1) ## (samples,)

            maf = '9' ## not used
            rsid_info = '\t'.join([str(chrom), rsid, post, ref, alt, maf])
            dose_vals = '\t'.join(dose.astype(str))
            dose_line = f'{rsid_info}\t{dose_vals}\n'
            dos_lines.append(dose_line)

        with open(output_file, 'w') as f:
            f.writelines(dos_lines)
        print(f'{TAB}done converting chr {chrom}.{version}')

    ''' 
    '''
    def __bgen_to_dosage(self, input_pattern, cleanup, num_threads=10, subset_size=5000): 
        print(f'\n{TAB}Converting genotypes to dosages:')
        pool = Pool(processes=num_threads)
        for chrom in range(1,23):
            ifile = f'{self.twas_genotypes_dir}/{input_pattern}' 
            ifile = ifile.replace('*', str(chrom))
            opath = f'{self.twas_genotypes_dir}/c{chrom}.dosage.txt' #.gz'
            #if os.path.exists(opath + '.gz'): continue

            bgen = bgr.read_bgen(ifile, verbose=False)
            snps = bgen['variants']
            gens = bgen['genotype']

            nvar = len(snps)
            nblocks = nvar // subset_size
            snps = snps.repartition(npartitions=nblocks)
            snp_blocks = snps.to_delayed()

            num_versions = np.ceil(nblocks / num_threads).astype(int)
            print(f'\n{TAB}Chr {chrom} contains {nvar} variants')
            print(f'{TAB}-- {nblocks} blocks across {num_versions} multi-threaded runs expected')

            if nblocks != len(snp_blocks):
                lb = len(snp_blocks)
                print(f'{TAB}-- NOTE: {nblocks} blocks expected, {lb} blocks created')

            params = []
            start = 0
            end = None
            t0 = time()
            for i, block in enumerate(snp_blocks):
                snp_set = block.compute()
                num_snp = snp_set.shape[0]
                end = start + num_snp

                gen_set = [g.compute()['probs'] for g in gens[start:end]]
                start = end
                print(f'{TAB}processing chr {chrom} block {i}...')

                outp_file = f'{self.twas_genotypes_dir}/tmp_c{chrom}_v{i}.dosage.txt'
                pp = (i, chrom, snp_set, gen_set, outp_file)
                params.append(pp)

                if (len(params) == num_threads) or (i == (nblocks-1)):
                    pool.starmap(self._compute_dosage_subset, params)
                    params = []

                    tpass = time() - t0
                    hr = int(tpass // 3600)
                    mn = int((tpass % 3600) // 60)
                    sc = int((tpass % 3600) % 60)

                    print(f'{TAB}> {hr} hr {mn} mn {sc} sc < for this run')
                    t0 = time()

            ## concat versions
            t0 = time()
            vpath = f'{self.twas_genotypes_dir}/tmp_c{chrom}_v*.dosage.txt'
            os.system(f'cat {vpath} >> {opath}')
            os.system(f'gzip {opath}')

            tpass = time() - t0
            hr = int(tpass // 3600)
            mn = int((tpass % 3600) // 60)
            sc = int((tpass % 3600) % 60)
            print(f'{TAB}combined all blocks into one file ({hr} hr {mn} mn {sc} sc)')
            
            ## delete temp dosage files 
            if cleanup:
                cmd = f'rm {self.twas_genotypes_dir}/tmp_c{chrom}_*.dosage.txt'
                os.system(cmd)

    '''
    FUSION version of predict_grex() using scores:
    Apply plink2-score to a chrom-specific list of genes.
    TODO: use .vars files to output num genes used vs in model
    TODO: decide how to communicate logs
    '''
    def _apply_fusion_weights(self, chrom, chr_genes): 
        accepted = 0
        for i, gene in enumerate(chr_genes):
            fus_file = f'{Models.FUS_SCORES(self.grex_region)}/{gene}.score'
            out_file = f'{self.out_path}/tmp_{gene}' ## .sscore

            cmd = [Scripts.PLINK2, '--silent', \
                   '--bfile', f'{self.twas_genotypes_dir}/c{chrom}', \
                   '--score', fus_file, '1', '2', '4', 'list-variants', \
                   '--out', out_file]
            try:
                res = subprocess.run(cmd, capture_output=True, text=True, check=True)
                accepted += 1
            except subprocess.CalledProcessError as e: 
                pass
        print(f'{TAB}{accepted} / {len(chr_genes)} genes predicted in chrom {chrom}')
                
    '''
    '''
    def __apply_metaxcan_weights(self, subj_file): 
        weight_file = Models.WEIGHT(self.twas_approach, self.grex_region)
        dos_pattern = f'{self.twas_genotypes_dir}/c*.dosage.txt.gz'
        cmd = ['python', '-u', Scripts.GREX_PREDICTION]
        cmd.extend(['--model_db_path', weight_file, \
                    '--text_genotypes', dos_pattern, \
                    '--text_sample_ids', subj_file, \
                    '--prediction_output', f'{self.out_path}.hdf5', 'HDF5', \
                    '--prediction_summary_output', f'{self.out_path}.log'])
        res = os.system(' '.join(cmd))

    '''
    Slice input bgen files to just the SNPs and 
    format needed for the TWAS models of interest. 

    FUS: .bed files 
    PDX/UTM/JTI: sliced .bgen to .dosage.txt.gz files
    '''
    def slice_bgen(self, input_dir, cleanup):  
        model_snps = Models.SNPS(self.twas_approach) 
        for c in range(1,23): 
            #if os.path.exists(f'{self.twas_genotypes_dir}/tmp_c{c}.bgen'): continue
            cmd = [Scripts.PLINK2, \
                   '--bgen', f'{input_dir}/c{c}.bgen', 'ref-first', \
                   '--sample', f'{input_dir}/c{c}.sample', \
                   '--extract', f'{model_snps}/c{c}.csv']

            if self.twas_approach == 'FUS': 
                cmd.append('--make-bed')
                cmd.extend(['--out', f'{self.twas_genotypes_dir}/c{c}'])
            else:
                cmd.extend(['--export', 'bgen-1.2', 'bits=8'])
                cmd.extend(['--out', f'{self.twas_genotypes_dir}/tmp_c{c}'])
            try:
                res = subprocess.run(cmd, capture_output=True, text=True, check=True)
                out = res.stdout.splitlines()
            except subprocess.CalledProcessError as e:
                print(f'{e.stderr}')
                sys.exit()

        ## for PDX/UTM/JTI files, generate dosage files 
        if self.twas_approach in ['PDX', 'UTM', 'JTI']:
            self.__bgen_to_dosage('tmp_c*.bgen', cleanup)
            if cleanup:
                os.system(f'rm {self.twas_genotypes_dir}/tmp_c*.bgen*') 
                os.system(f'rm {self.twas_genotypes_dir}/tmp_c*.log') 
                os.system(f'rm {self.twas_genotypes_dir}/tmp_c*.sample') 
        return

    '''
    '''
    def apply_weights(self, samples, samp_file):
        if self.twas_approach == 'FUS':
            ## get FUS-valid list of genes and chroms
            weight_file = Models.WEIGHT('FUS', self.grex_region)
            fus_genes = pd.read_table(weight_file, sep='\t', usecols=['WGT', 'CHR'])
            fus_genes['gene'] = fus_genes['WGT'].apply(lambda x: x.split('/')[1].split('.')[0])
        
            ## apply plink2 score to each gene (parallelize by chrom)
            params = [(ch, cdf['gene'].values) for ch, cdf in fus_genes.groupby('CHR')]
            pool = Pool(processes=8) ## TODO: make num_threads an option
            pool.starmap(self._apply_fusion_weights, params)

            ## gather all .sscore files into one hdf5 file
            gene_files = glob.glob(f'{self.out_path}/tmp_*.sscore')
            gene_arr = [i.split('_')[-1].split('.')[0] for i in gene_files]
            grex_mat = np.zeros((len(gene_arr), samples.size))

            for i, gfile in enumerate(gene_files): 
                df = pd.read_table(gfile, sep='\t', \
                                   usecols=['IID', 'SCORE1_AVG'], \
                                   index_col=['IID']) 
                df.index = df.index.astype(str)
                df = df.reindex(samples) 
                grex_mat[i] = df['SCORE1_AVG'].values
            
            with h5py.File(f'{self.out_path}.hdf5', 'w') as f:
                f['genes'] = np.array(gene_arr).astype(bytes)
                f['pred_expr'] = grex_mat

        else:
            self.__apply_metaxcan_weights(samp_file)
        
