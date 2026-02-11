## AGENT-P: The Association of Gene Expression to Neuroimaging Traits Pipeline

### Table of Contents

 - [Aims](#aims)
 - [Features](#features)
 - [Installation](#installation)
 - [Set Up](#set-up)
 - [Quick-start Example](#quick-start-example)
 - [Supported TWAS Models](#supported-twas-models)
 - [Developer Notes](#developer-notes)
 - [Usage Notes](#usage-notes)
     - [Starting a New Project](#starting-a-new-project)
     - [Running a Summary TWAS](#running-a-summary-twas)
     - [Running a Genotype TWAS](#running-a-genotype-twas)
     - [Running a GWAS](#running-a-gwas)
___
___
### Aims 
* develop a user-friendly python package for performing TWAS
* make it easier to test associations from different approaches and parameters
    * approaches (PrediXcan, FUSION, UTMOST, JTI)
    * TWAS based on summary stats vs genotype data
* standardize results from different TWAS approaches for comparative analyses
___
___
### Features 
* GWAS ([REGENIE](https://rgcgithub.github.io/regenie) wrapper)
* TWAS (pre-trained models from [MetaXcan](https://github.com/hakyimlab/MetaXcan), [FUSION](http://gusevlab.org/projects/fusion), [UTMOST](https://github.com/Joker-Jerome/UTMOST), and [JTI](https://github.com/gamazonlab/MR-JTI))
    * Models are available for the following brain regions: DLPFC, anterior cingulate, amygdala, hippocampus, caudate, putamen, nucleus accumbens, and cerebellar hemisphere 
* Enrichment Analyses ([WebGestalt](https://www.webgestalt.org/) wrapper)
* Visual Analytics
___
___

### Installation

This project supports multiple installation workflows depending on your environment and preferences:

**A)** Manual setup — Follow the instructions in `docs/installation.qmd` to install dependencies and run the program manually. This method provides the most transparency and flexibility for users who want to modify configurations or troubleshoot components.

**B)** Dockerfile setup — Follow the instructions in `docs/documentation.qmd` to build and run the provided Dockerfile. This approach ensures a reproducible and isolated runtime environment without requiring manual dependency installation.

**C)** Automated setup — Directly run `quick_run.sh` using following command to automatically download the required .sh files, install dependencies, and execute the example pipeline in one step. 
```
bash quick_run.sh
```
___
___

### Set Up

My machine runs the following specs: 
 - GCC/10.2.0 
 - OpenMPI/4.0.5 
 - Python/3.8.6
 - R/4.0.5

**A)** First, clone the AGENT-P repository (currently called `dev_agentp`). Note that for this beta version, you'll need to be able to import classes that exist in this repo. 
```
git clone --recurse-submodules git@github.com:nhunghoang/dev_agentp.git
```

**B)** Download and extract the TWAS `models.tar.gz` file into the main level of the repository. [[link to file]](https://drive.google.com/file/d/1llysnSFOut3E_g1omlgQfi-N5fptuuLa/view?usp=sharing)
```
cd dev_agentp
tar -xvzf models.tar.gz
```

**C)** I recommend running AGENT-P in a python virtual environment with the following python packages. If using `pip`, you can call `pip install -r dev_agentp/requirements.txt`. 

 - numpy >= 1.23.5
 - pandas >= 1.5.3
 - scipy >= 1.10.1
 - statsmodels >= 0.14.1
 - bgen_reader >= 4.0.8
 - pyliftover >= 0.4.1
 - h5py >= 3.7.0
 - patsy >= 0.5.6
 - bitarray >= 3.0.0

**D)** Download the appropriate alpha build of `plink2` from [here](https://www.cog-genomics.org/plink/2.0/) and unzip it into the `dev_agentp/software` folder.  

**E)** Next, you have to manually download `bgenix` and compile it. On your UNIX terminal:
```
cd dev_agentp/software
wget http://code.enkre.net/bgen/tarball/release/bgen.tgz
tar -xvzf bgen.tgz 
mv bgen.tgz bgen
cd bgen
./waf configure
./waf
```

**F)** Then, compile the provided REGENIE software: 
```
cd dev_agentp/software/regenie
make
```

**G)** Lastly, launch R from inside the main repo (`dev_agentp`) and run the following commands: 
```
install.packages(c('optparse','RColorBrewer'))
install.packages('software/plink2R/plink2R/',repos=NULL)
```

Hopefully that's it! 
___
___
### Quick-start Example
Here's example code that tests out all the features of AGENT-P. This code is available as a script (`example/run_example.py`) and will run from within the `dev_agentp` directory. This code uses the data provided in `example/example_data`. 
```
from Project import Project 
from TWAS import SummaryTWAS, GenotypeTWAS 
from GWAS import GWAS 

## init project 
main_project_folder = 'agentp_test'
project = Project(main_project_folder)

'''
init a summary twas that uses JTI models trained on hippocampus genes
'''
stwas = SummaryTWAS(project, 'JTI', 'hippocampus') 

## apply these models to a GWAS on amygdala volumes 
stwas.add_gwas('amygdala_volume', 'test_data/vol_mean_amygdala.regenie') 

## run the summary TWAS 
stwas.run_twas('amygdala_volume') 

## TWAS results are saved to
## dev_agentp/agentp_test/sTWAS/JTI_hippocampus/amygdala_volume.csv

'''
init a genotype twas that uses FUSION models trained on caudate genes
'''
gtwas = GenotypeTWAS(project, 'FUS', 'caudate')

## add subject data to the project 
project.set_subjects('test_data/subjects.txt') 
project.add_covariates('test_data/covariates.csv')
project.add_phenotypes('test_data/volumes.csv')
project.add_genotypes('test_data/genotypes', 'test_c*.bgen')

## predict FUSION-caudate grex for this cohort 
gtwas.predict_grex() 

## run the genotype TWAS on putamen volumes
gtwas.run_twas('putamen_volume')

## TWAS results are saved to
## dev_agentp/agentp_test/gTWAS/FUS_caudate/putamen_volume.csv

'''
init a gwas 
'''
gwas = GWAS(project) 

## no need to add the subject data again, since we already did for gTWAS
## run GWAS on nucleus hippocampus volumes 
gwas.run('hippocampus_volume')

## GWAS results are saved to
## dev_agentp/agentp_test/GWAS/hippocampus_volume.regenie
```
___
___
### Supported TWAS Models
When specifying TWAS approaches and brain regions, use the following abbreviations:  

| TWAS Approach / Brain Region | Abbreviation in AGENT-P |
|:--|:--|
| PrediXcan/MetaXcan | PDX |
| FUSION | FUS |
| UTMOST | UTM |
| Joint Tissue Imputation | JTI |
| dorsolateral prefrontal cortex | dlpfc |
| anterior cingulate | ant-cingulate |
| amygdala | amygdala |
| hippocampus | hippocampus |
| caudate | caudate |
| putamen | putamen |
| nucleus accumbens | nuc-accumbens |
| cerbellar hemisphere | cerebellar-hemi |

___
___
### Developer Notes 
Main Python Classes: 
* **Project** — manages input/output file organization 
* **Static** — list of static file paths and variables needed for this software
* **Genotypes** — applies quality control to input genotypes (currently supports BGEN, aiming for VCF next)
* **GWAS** — applies REGENIE to an input set of genotypes and traits 
* **TWAS** —  applies pre-trained GREx models to traits
    * **sTWAS** —  applies summary TWAS to a set of GWAS results
    * **gTWAS** — applies genotype TWAS to a set of genotypes and traits (see GREX class) 
* **GREX** — predicts GREx for an input set of genotypes, which then get passed to gTWAS 
___
___
### Usage Notes

### STARTING A NEW PROJECT 
Every AGENT-P analysis begins with a `Project` object, which creates a new directory for **one cohort** and will contain all GWAS and TWAS for this cohort. 

In a python script, call: 

    from Project import Project 
    project = Project('ukb_test')

This will create a directory called "ukb_test" in your current path. You can replace this string with any name you'd like to represent your cohort. 
___
### RUNNING A SUMMARY TWAS 
A summary TWAS requires three params: 1) a GWAS summary statistics file, 2) the TWAS approach you want to use, and 3) the brain region you want GREx models from. The  output file of associations will be saved using this file path pattern: `project/sTWAS/approach_region/phenotype.csv`

**Example Run:**  
```
# This will identify associations between 
# putamen volume and GREx from JTI caudate models. 

from TWAS import SummaryTWAS 
stwas = SummaryTWAS(project, 'JTI', 'caudate') 
stwas.add_gwas('putamen_volume', 'GWAS/putamen_volume.csv')
stwas.run_twas('putamen_volume')

# Output TWAS file will be saved as: 
# ukb_test/sTWAS/JTI_caudate/putamen_volume.csv
```

**Summary TWAS Methods:**
```
__init__(self, project:Project, twas_approach:str, grex_region:str)
```
 - **project**: an instance of the Project class
 - **twas_approach**: string abbreviation of a supported TWAS approach
 - **grex_region**: string abbreviation of a supported GREx brain region

```
add_gwas(self, phenotype:str, gwas_file:str)

# Optional parameters (all type str and default to None)
# sep, snp_col, ref_col, alt_col, beta_col, beta_sign_col, beta_se_col, or_col, zscore_col, pval_col
```
* **phenotype**: string name to identify the phenotype represented in the input GWAS 
    * Name does not have to match the input GWAS file name
* **gwas_file**: a tabulated file containing GWAS results for one phenotype 
    * Must contain at least these columns: *rsid*, *reference allele*, *alternative allele*, and a *test statistic* (e.g., p-value, zscore, effect size) 
* The optional parameters specify column names and are used to parse the GWAS file. If none of them are provided,  the file will be parsed like a REGENIE output file: 
    * **sep** = "," // **snp_col** = "SNP" // **ref_col** = "ALLELE0" // **alt_col** = "ALLELE1" // **beta_col** = "BETA" // **beta_se_col** = "SE" // **pval_col** = "PVAL" // None for all other params

```
run_twas(self, phens:str or list of str, overwrite:bool=False, num_threads:int=1)
```
* **phens**: string or list of strings specifying which phenotypes to run summary TWAS on 
    * String(s) must match the phenotype names used in `SummaryTWAS.add_gwas()`
* **overwrite**: boolean specifying whether or not to overwrite an existing TWAS file (default is False) 
* **num_threads**: number of threads to use when calling python's `multiprocessing.pool`  (default is 1) 
    * This is useful when calling `SummaryTWAS.run_twas()` on multiple phenotypes at once because it determines how many phenotype TWAS are ran in parallel 
___
### RUNNING A GENOTYPE TWAS 
A genotype TWAS requires 1) a list of subject IDs for the cohort of interest, 2) a covariate file, 3) a phenotype file, and 4) one genotype file each for chromosomes 1-22. Additionally, you will need to specify the TWAS approach you want to use and the brain region that you want GREx models from. The output file of associations will be saved using this file path pattern: `project/gTWAS/approach_region/phenotype.csv`

**Example Run:**
```
# This will identify associations between 
# putamen volume and GREx from JTI caudate models. 

# First, add the subject-level data to the project 
project.set_subjects('subject_ids.txt') 
project.add_covariates('covariates.csv', covar_cols=['age', 'sex', 'PC1'])
project.add_phenotypes('volumes.csv', phen_cols=['putamen_volume'])
project.add_genotypes('bgen_genotypes', 'chrom*.bgen')

# Then, run genotype TWAS
from TWAS import GenotypeTWAS 
gtwas = GenotypeTWAS(project, 'JTI', 'caudate') 
gtwas.predict_grex()
gtwas.run_twas('putamen_volume')

# Output TWAS file will be saved as: 
# ukb_test/gTWAS/JTI_caudate/putamen_volume.csv
```
**Project Methods needed for a Genotype TWAS:** 
```
set_subjects(self, subjectID_file:str)
```
Note: This method only needs to be called once per Project instance (subsequent calls will be ignored). 
* **subjectID_file**: string name of file that contains subject IDs for the cohort of interest 
    * Expected format is one subject ID per row 
```
add_covariates(self, covar_file:str, id_col:str='id', sep:str='\t', covar_cols:list_of_str=None)
```
Note: This method can only be called once per Project instance (subsequent calls will be ignored).
* **covar_file**: string name of tabulated file containing subject-level covariate data
    * Expected: rows represent subjects and columns represent covariates 
* **id_col**: string name of column that contains subject IDs (default is 'id')
    * Only rows with subject IDs that match the cohort of interest will be kept 
* **sep**: string separator/delimiter (default is tab) 
* **covar_cols**: string list of column names that represent your covariates of interest (default is None, meaning all covariates in the file will be included)
    * Only the specified covariates will be read and included in the TWAS
```
add_phenotypes(self, phen_file:str, id_col:str='id', sep:str='\t', phen_cols:list_of_str=None)
```
Note: This method can be called multiple times per Project instance, but phenotype column names that match existing phenotypes in this Project will be ignored. 
* **phen_file**: string name of tabulated file containing subject-level phenotype data
    * Multiple phenotypes in one file are allowed and will be parsed separately 
    * Expected: rows represent subjects and columns represent phenotypes  
* **id_col**: string name of column that contains subject IDs (default is 'id')
    * Only rows with subject IDs that match the cohort of interest will be kept 
* **sep**: string separator/delimiter (default is tab) 
* **covar_cols**: string list of column names that represent your phenotypes of interest (default is None, meaning all phenotypes in the file will be included)
    * Only the specified phenotypes will be read and included in the TWAS
```
add_genotypes(self, input_geno_dir:str, file_pattern:str)
```
Note: This method can only be called once per Project instance (subsequent calls will be ignored). 
* **input_geno_dir**: string path of the directory containing genotype files for the cohort of interest 
* **file_pattern**: string naming convention for files in *input_geno_dir* (with an asterisk (*) as the placeholder for chromosomes) 
    * Assumes the naming convention for chromosomes 1-22 are the same 
    * Example: if `file_pattern = 'chrom*.bgen'`, then genotype files are assumed to be at `input_geno_dir/chrom1.bgen`, `input_geno_dir/chrom2.bgen`, etc. 

**Genotype TWAS Methods:**
```
__init__(self, project:Project, twas_approach:str, grex_region:str)
```
 - **project**: an instance of the Project class
 - **twas_approach**: string abbreviation of a supported TWAS approach
 - **grex_region**: string abbreviation of a supported GREx brain region
```
predict_grex(self, overwrite:boolean=False, cleanup:boolean=False)
```
This method will predict GREx for all subjects in the cohort using the specified TWAS approach and models. You only need to call this method once for a given  approach-region pair. 
- **overwrite**: boolean specifying whether or not to overwrite the existing GREx file (default is False) 
- **cleanup**: boolean specifying whether or not to delete intermediary files that were created during GREx prediction (default is False) 
    - I suggest using `cleanup = True` when applying FUSION models, and `cleanup = False` when applying any other models (until you're done predicting GREx) because PrediXcan, UTMOST, and JTI models generate the same intermediary files. 
```
run_twas(self, phens:str or list of str, overwrite:bool=False, num_threads:int=1)
```
* **phens**: string or list of strings specifying which phenotypes to run genotype TWAS on 
    * String(s) must match the phenotype names used in `Project.add_phenotypes()`
* **overwrite**: boolean specifying whether or not to overwrite an existing TWAS file (default is False) 
* **num_threads**: number of threads to use when calling python's `multiprocessing.pool`  (default is 1) 
    * This is useful when calling `GenotypeTWAS.run_twas()` on multiple phenotypes at once because it determines how many phenotype TWAS are ran in parallel 
___
### RUNNING A GWAS 
Like a GenotypeTWAS, the required input files for a GWAS are 1) a list of subject IDs for the cohort of interest, 2) a covariate file, 3) a phenotype file, and 4) one genotype file each for chromosomes 1-22. The output file of associations will be saved using this file path pattern: `project/GWAS/phenotype.regenie`


 **Example Run:** 
```
# This will identify associations 
# between SNPs and putamen volumes and, 
# separately, SNPS and caudate volumes. 

# Add subject-level data to the project, if not already 
# (see GenotypeTWAS example above)

# Then run GWAS for putamen volume and caudate volume (in parallel)
from GWAS import GWAS
gwas = GWAS(project)
gwas.run(['putamen_volume', 'caudate_volume'])

# Output GWAS files will be saved as: 
# ukb_test/GWAS/putamen_volume.regenie
# ukb_test/GWAS/caudate_volume.regenie
```
**GWAS Methods:**
```
run(self, phens:str or list of str, overwrite:bool=False, cleanup:bool=True, num_threads:int=1)
```
Note: For multiple phenotypes of interest, it is more efficient to call this method on the list of phenotypes, rather than each phenotype separately. 

 - **phens**: string or list of strings specifying which phenotypes to run GWAS on 
     - String(s) must match the phenotype names used in `Project.add_phenotypes()`
 - **overwrite**: boolean specifying whether or not to overwrite existing GWAS file(s) for the input phenotypes (default is False)
 - **cleanup**: boolean specifying whether or not to delete intermediary files that were created during GWAS (default is True) 
 - **num_threads**: number of threads to use when calling python's `multiprocessing.pool`  (default is 1) 
     - The number of threads determines how many chromosomes are processed at once in parallel (regardless of the number of phenotypes). 

