'''
Project: keeps track of inputs and file organization for one cohort 

(Future) Create an "Input"/"Cohort" class and allow multiple cohorts under one Project

- Nhung, Sept 2024
'''

import pandas as pd
import os
import sys
import json
import sqlite3

from .Static import TAB, MSG
from .Static import Models

##########################################################################################

class Project: 

    def __init__(self, project_dir): 
        if os.path.exists(project_dir): 
            self.path = os.path.abspath(project_dir)
            msg = f'{MSG}Project already exists. Results will be saved to:\n{TAB}{self.path}'
            print(msg)
        else:
            os.mkdir(project_dir)
            os.mkdir(f'{project_dir}/input_data')
            self.path = os.path.abspath(project_dir)
            msg = f'{MSG}Created new project directory:\n{TAB}{self.path}\n'
            print(msg)

        self.input_paths = {'subj_file': f'{self.path}/input_data/subject_IDs.csv', \
                            'covs_file': f'{self.path}/input_data/covariates.csv', \
                            'phen_dir': f'{self.path}/input_data/phenotypes', \
                            'geno_dir': f'{self.path}/input_data/genotypes', \
                            'ldsc_dir': f'{self.path}/input_data/gwas_ldsc', \
                            'tgen_dir': f'{self.path}/input_data/twas_genotypes', \
                            'model_dir': f'{self.path}/input_data/custom_models'}

        self.output_paths = {'gwas': f'{self.path}/GWAS', \
                             'stwas': f'{self.path}/sTWAS', 
                             'itwas': f'{self.path}/iTWAS', \
                             'grex': f'{self.path}/GREX', \
                             'webg': f'{self.path}/enrichment'}

        ## get ordered array of subjects (for iTWAS) if exists
        try: 
            self.subjects = pd.read_table(self.input_paths['subj_file'], \
                                          sep=' ', \
                                          dtype=str, \
                                          header=None)[0].values
        except FileNotFoundError as e:
            self.subjects = None

        ## register model sets previously added to this project (see add_model)
        self.__load_custom_models()


    '''
    Add a set of PredictDB-format expression models (e.g., a new brain region
    for an existing approach, or a new approach) for summary and individual TWAS.
    weights: PredictDB .db file with "weights" (gene, rsid, weight, ref_allele, eff_allele)
             and "extra" (gene, ...) tables
    covariance: matching SNP covariance file (.txt.gz), used by summary TWAS
    The model set is saved in the project, so it only needs to be added once.
    '''
    def add_model(self, twas_approach, grex_region, weights, covariance):
        if twas_approach == 'FUS':
            msg = f'{MSG}MODEL ERROR: Only PredictDB-format models can be added; ' \
                  f'use an approach name other than "FUS"'
            sys.exit(msg)
        if (twas_approach in Models.method) and (grex_region in Models.rlabel):
            msg = f'{MSG}MODEL ERROR: "{twas_approach} {grex_region}" is a built-in model set; ' \
                  f'choose a different approach or region name'
            sys.exit(msg)
        for mfile in [weights, covariance]:
            if not os.path.exists(mfile):
                msg = f'{MSG}MODEL ERROR: Model file "{mfile}" does not exist'
                sys.exit(msg)

        ## check the PredictDB format and get the model SNPs
        try:
            with sqlite3.connect(f'file:{os.path.abspath(weights)}?mode=ro', uri=True) as con:
                tables = [t for (t,) in con.execute("SELECT name FROM sqlite_master WHERE type='table'")]
                if not {'weights', 'extra'}.issubset(tables):
                    raise ValueError(f'missing "weights" and/or "extra" tables (found: {", ".join(tables)})')
                wcols = [c[1] for c in con.execute('PRAGMA table_info(weights)')]
                need = ['gene', 'rsid', 'weight', 'ref_allele', 'eff_allele']
                if not set(need).issubset(wcols):
                    raise ValueError(f'"weights" table needs columns {", ".join(need)} (found: {", ".join(wcols)})')
                snps = pd.read_sql('SELECT DISTINCT rsid FROM weights', con)['rsid']
                n_genes = con.execute('SELECT COUNT(DISTINCT gene) FROM weights').fetchone()[0]
        except (sqlite3.DatabaseError, ValueError) as e:
            msg = f'{MSG}MODEL ERROR: "{weights}" is not a PredictDB model database:\n{TAB}{e}'
            sys.exit(msg)

        ## save the SNP list (for slicing genotypes in individual TWAS) and the model paths
        model_dir = self.input_paths['model_dir']
        if not os.path.exists(model_dir):
            os.mkdir(model_dir)
        snp_file = f'{model_dir}/{twas_approach}_{grex_region}.snps.txt'
        snps.sort_values().to_csv(snp_file, index=False, header=False)

        ## the SNP list is saved relative to the project, so the project folder can be moved
        entry = {'approach': twas_approach, 'region': grex_region, \
                 'weight': os.path.abspath(weights), \
                 'covariance': os.path.abspath(covariance), \
                 'snps': os.path.basename(snp_file)}
        reg_file = f'{model_dir}/models.json'
        registry = json.load(open(reg_file)) if os.path.exists(reg_file) else {}
        key = f'{twas_approach}_{grex_region}'
        if (key in registry) and (registry[key] != entry):
            print(f'{MSG}WARNING: Replacing the "{twas_approach} {grex_region}" models added earlier ' \
                  f'(existing TWAS results for this model set are not rerun unless overwrite=True)')
        registry[key] = entry
        with open(reg_file, 'w') as f:
            json.dump(registry, f, indent=2)

        Models.CUSTOM[(twas_approach, grex_region)] = dict(entry, snps=snp_file)
        msg = f'{MSG}Added {twas_approach} {grex_region} models ({n_genes} genes, {len(snps)} SNPs)\n'
        print(msg)

    '''
    Register the model sets saved in this project by add_model().
    '''
    def __load_custom_models(self):
        reg_file = f'{self.input_paths["model_dir"]}/models.json'
        if not os.path.exists(reg_file):
            return
        for entry in json.load(open(reg_file)).values():
            entry['snps'] = f'{self.input_paths["model_dir"]}/{entry["snps"]}'
            missing = [entry[k] for k in ['weight', 'covariance', 'snps'] if not os.path.exists(entry[k])]
            if len(missing) > 0:
                print(f'{MSG}WARNING: Skipping added {entry["approach"]} {entry["region"]} models - ' \
                      f'file(s) not found: {", ".join(missing)}')
                continue
            Models.CUSTOM[(entry['approach'], entry['region'])] = entry

    '''
    Set the subject order for the entire project.
    '''
    def set_subjects(self, subjectID_file): 
        if not (self.subjects is None):
            msg = f'{MSG}WARNING: Project cohort is already set - input file will be ignored\n'
            print(msg)
        else: 
            try:
                new_sids = pd.read_table(subjectID_file, header=None, dtype=str)[0].values
                self.subjects = new_sids
                new_vals = pd.DataFrame({0: new_sids, 1: new_sids}, dtype=str)
                new_vals.to_csv(self.input_paths['subj_file'], sep=' ', index=False, header=False)
                msg = f'{MSG}Saved {new_sids.size} IDs from {subjectID_file}\n'
                print(msg) 
            except Exception as e: 
                msg = f'{MSG}Error reading {subjectID_file}:\n{TAB}{e}\n'
                sys.exit(msg)

    '''
    Add covariates (in project cohort order).
    '''
    def add_covariates(self, covar_file, id_col='id', sep='\t', covar_cols=None):
        if self.subjects is None: 
            msg = f'{MSG}ERROR: Subject list has not been set for this project yet\n'
            sys.exit(msg)

        cur_file = self.input_paths['covs_file'] 
        if os.path.exists(cur_file): 
            msg = f'{MSG}WARNING: Project covariates are already set - input file will be ignored\n'
            print(msg)
            return
        try: 
            if covar_cols: covar_cols.append(id_col)
            cov_data = pd.read_table(covar_file, index_col=id_col, sep=sep, usecols=covar_cols)
        except Exception as e: 
            msg = f'{MSG}Error reading {covar_file}:\n{TAB}{e}\n'
            sys.exit(msg)

        cov_data.index = cov_data.index.astype(str)
        cov_data = cov_data.reindex(self.subjects)

        if not (cov_data.index == self.subjects).all():
            msg = f'{MSG}ERROR: Input covariates file does not include all subjects in this project\n'
            sys.exit(msg)

        ## add FID and IID columns (needed for regenie)
        cov_data.index.rename('FID', inplace=True)
        cov_data.insert(0, 'IID', cov_data.index)

        ## save
        cov_data.to_csv(cur_file, sep='\t')
        num_covs = cov_data.columns.size - 1
        msg = f'{MSG}Added {num_covs} covariates to this project\n'
        print(msg)

    '''
    Add genotypes and apply QC. 
    TODO: reorder subjects if needed / add VCF input feature
    TODO: let users define params
    '''
    def add_genotypes(self, input_geno_dir, file_pattern):
        if self.subjects is None: 
            msg = f'{MSG}ERROR: Subject list has not been set for this project yet\n'
            sys.exit(msg)

        proj_geno_dir = self.input_paths['geno_dir']
        exists = [os.path.exists(f'{proj_geno_dir}/c{i}.bgen') for i in range(1,23)]
        if sum(exists) == 22:
            msg = f'{MSG}WARNING: Project genotypes are already set - input files will be ignored\n'
            print(msg)
            return
        else:
            if not os.path.exists(proj_geno_dir):
                os.mkdir(proj_geno_dir) 
                os.mkdir(f'{proj_geno_dir}/tmp')

        ## apply QC and save 
        if file_pattern.split('.')[-1] == 'bgen':
            from .Genotypes import BGEN 
            bgen = BGEN(input_geno_dir, file_pattern[:-5], self.input_paths['subj_file'], proj_geno_dir)
            bgen.apply_qc()
            msg = f'{MSG}Applied QC to genotypes and added to project\n'
            print(msg)
        else: 
            msg = 'This genotype file format is not currently supported\n'
            sys.exit(msg)

    '''
    Add phenotypes (in project cohort order).
    Currently works by saving every column of an input file as a phen file in the project. 
    '''
    def add_phenotypes(self, phen_file, id_col='id', sep='\t', phen_cols=None):
        if self.subjects is None: 
            msg = f'{MSG}ERROR: Subject list has not been set for this project yet\n'
            sys.exit(msg)

        phen_dir = self.input_paths['phen_dir']
        if not os.path.exists(phen_dir): 
            os.mkdir(phen_dir) 

        try:
            if phen_cols: phen_cols.append(id_col) 
            phen_data = pd.read_table(phen_file, index_col=id_col, sep=sep, usecols=phen_cols)
        except Exception as e: 
            msg = f'{MSG}Error reading {phen_file}:\n{TAB}{e}\n'
            sys.exit(msg)

        phen_data.index = phen_data.index.astype(str)
        phen_data = phen_data.reindex(self.subjects)

        if not (phen_data.index == self.subjects).all():
            msg = f'{MSG}ERROR: Input phenotypes file does not include all subjects in this project\n'
            sys.exit(msg)

        ## add FID and IID columns (needed for regenie)
        phen_data.index.rename('FID', inplace=True)
        phen_data.insert(0, 'IID', phen_data.index)

        ## save
        print(f'{MSG}Parsing phenotype data...')
        for phen in phen_data.columns[1:]:
            new_phen_file = f'{phen_dir}/{phen}.csv'
            if os.path.exists(new_phen_file): 
                print(f'{TAB}WARNING: "{phen}" column ignored - phenotype with this name already exists in the project') 
                continue

            phen_data.to_csv(new_phen_file, columns=['IID', phen], sep='\t')
            print(f'{TAB}Added "{phen}" phenotype to this project')
        print('')

    '''
    Init GWAS analysis paths.
    '''
    def set_gwas_paths(self): 
        gwas_dir = self.output_paths['gwas'] 
        if not os.path.exists(gwas_dir): 
            os.mkdir(gwas_dir)
            os.mkdir(f'{gwas_dir}/tmp')
        return {'output': gwas_dir}

    '''
    Init TWAS analysis paths. 
    '''
    def set_twas_paths(self, twas_type, twas_approach, grex_region):
        twas_top = self.output_paths[twas_type] 
        twas_dir = f'{twas_top}/{twas_approach}_{grex_region}'
        if not os.path.exists(twas_top): 
            os.mkdir(twas_top) 
            os.mkdir(twas_dir)
        elif not os.path.exists(twas_dir): 
            os.mkdir(twas_dir)
        return twas_dir

    '''
    Init additional paths for Summary TWAS. 
    '''
    def set_summary_twas_paths(self, twas_approach, grex_region):
        twas_dir = self.set_twas_paths('stwas', twas_approach, grex_region)
        if twas_approach == 'FUS': 
            gwas_dir = self.input_paths['ldsc_dir']
            if not os.path.exists(gwas_dir): 
                os.mkdir(gwas_dir)
            if not os.path.exists(f'{twas_dir}/tmp'):
                os.mkdir(f'{twas_dir}/tmp')
            return {'ldsc': gwas_dir, 'output': twas_dir}
        else:
            return {'output': twas_dir}

    '''
    Init additional paths for Individual TWAS. 
    '''
    def set_individual_twas_paths(self, twas_approach, grex_region): 

        ## twas genotype dirs
        tgen_dir = self.input_paths['tgen_dir']
        fuse_dir = f'{tgen_dir}/bed_FUS'
        dose_dir = f'{tgen_dir}/dos_PDX_UTM_JTI'
        ## added models slice their own SNPs, so they get their own dosage files
        if Models.is_custom(twas_approach, grex_region):
            dose_dir = f'{tgen_dir}/dos_{twas_approach}_{grex_region}'
        if not os.path.exists(tgen_dir):
            os.mkdir(tgen_dir)
        twas_gen_dir = fuse_dir if (twas_approach == 'FUS') else dose_dir
        if not os.path.exists(twas_gen_dir):
            os.mkdir(twas_gen_dir)

        ## grex dirs and file
        grex_top = self.output_paths['grex']
        grex_dir = f'{grex_top}/{twas_approach}'
        if not os.path.exists(grex_top): 
            os.mkdir(grex_top)
        if not os.path.exists(grex_dir): 
            os.mkdir(grex_dir) 

        if twas_approach == 'FUS': 
            greg_dir = f'{grex_dir}/{grex_region}'
            if not os.path.exists(greg_dir): 
                os.mkdir(greg_dir)

        grex_file = f'{grex_dir}/{grex_region}.hdf5'
        twas_dir = self.set_twas_paths('itwas', twas_approach, grex_region)
        return {'twas_gen': twas_gen_dir, 'grex': grex_file, 'output': twas_dir}

    '''
    Init paths for Enrichment. 
    '''
    def set_enrichment_paths(self, twas_object, phen, pval_type, pval_alpha): 
        enr_dir = self.output_paths['webg']
        ref_dir = f'{enr_dir}/reference_sets'
        if not os.path.exists(enr_dir): 
            os.mkdir(enr_dir)
            os.mkdir(ref_dir)
        elif not os.path.exists(ref_dir): 
            os.mkdir(ref_dir)

        fram = twas_object.framework[0]
        appr = twas_object.approach
        grex = twas_object.grex_region
        out_dir = f'enr_dir/{phen}_{appr}_{grex}_{fram}TWAS'
        if not os.path.exists(out_dir): 
            os.mkdir(out_dir)

        out_dir2 = f'{out_dir}/{pval_type}_{pval_alpha}'
        if not os.path.exists(out_dir2): 
            os.mkdir(out_dir2)

        return {'ref_dir': ref_dir, 'out_dir': out_dir2}

    def set_custom_enrichment_paths(self): 
        enr_dir = self.output_paths['webg']
        ont_dir = f'{enr_dir}/custom_ontologies'
        if not os.path.exists(enr_dir): 
            os.mkdir(enr_dir)
            os.mkdir(ont_dir)
        elif not os.path.exists(ont_dir): 
            os.mkdir(ont_dir)
        return {'custom_ontology_dir': ont_dir}

    '''
    get() functions
    '''
    def get_subjects(self):
        return self.subjects

    def get_subjects_file(self):
        return self.input_paths['subj_file']

    def get_covariates(self):
        cfile = self.input_paths['covs_file']
        covs = pd.read_table(cfile, index_col='IID').drop(columns=['FID'])
        return covs

    def get_covariates_file(self): 
        return self.input_paths['covs_file']

    def get_genotypes_dir(self):
        return self.input_paths['geno_dir']

    def get_phenotypes(self, phenotype):
        phen_dir = self.input_paths['phen_dir']
        phen_file = f'{phen_dir}/{phenotype}.csv'
        phens = pd.read_table(phen_file, index_col='IID').drop(columns=['FID'])
        return phens

    def get_phenotype_file(self, phenotype):
        phen_dir = self.input_paths['phen_dir']
        return f'{phen_dir}/{phenotype}.csv'

    def get_phenotypes_dir(self): 
        return self.input_paths['phen_dir']

