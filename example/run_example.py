'''
Example run of AGENT-P on the example data, using the small example models in
example/example_models (caudate models; JTI is used here).

  - summary TWAS of amygdala volume, from a GWAS of the example cohort
    (example_data/amygdala_volume_gwas.regenie.gz, made with REGENIE from the example genotypes)
  - individual TWAS of amygdala volume, from the example genotypes, covariates, and phenotypes
  - a Manhattan plot of the individual TWAS, and a Miami plot of the summary (top) vs
    individual (bottom) TWAS

Runs in a few minutes. Needs plink2 and bgenix on PATH.
Output: example/example_run (TWAS results in sTWAS/ and iTWAS/, plots in figures/)

The example models also include PrediXcan ('PDX'), UTMOST ('UTM'), and FUSION ('FUS'; needs
Rscript) caudate models: change `approach` below to try them.

usage: python run_example.py
'''

import os

from agentp import Project, SummaryTWAS, IndividualTWAS, set_models_dir
from agentp import load_pvalues, manhattan_plot, miami_plot

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

## use the small example models (caudate only) instead of the full models
set_models_dir('example_models')

approach = 'JTI'
region = 'caudate'
phen = 'amygdala_volume'

########################################################################

## init project and add the subject-level data (genotypes go through QC)
project = Project('example_run')
project.set_subjects('example_data/subjects.txt')
project.add_covariates('example_data/covariates.csv')
project.add_phenotypes('example_data/volumes.csv')
project.add_genotypes('example_data/genotypes', 'test_c*.bgen')

## summary TWAS: apply the JTI caudate models to the GWAS of amygdala volume in the example cohort
stwas = SummaryTWAS(project, approach, region)
stwas.add_gwas(phen, 'example_data/amygdala_volume_gwas.regenie.gz')
stwas.run_twas(phen)

## individual TWAS: predict JTI caudate GREx for the cohort, then test amygdala volume
itwas = IndividualTWAS(project, approach, region)
itwas.predict_grex()
itwas.run_twas(phen)

## results are saved to example_run/sTWAS/JTI_caudate/amygdala_volume.csv
##                   and example_run/iTWAS/JTI_caudate/amygdala_volume.csv
s = load_pvalues(f'example_run/sTWAS/{approach}_{region}/{phen}.csv')
i = load_pvalues(f'example_run/iTWAS/{approach}_{region}/{phen}.csv')

## plots
os.makedirs('example_run/figures', exist_ok=True)
label = f'{approach} {region}'
manhattan_plot(i, label=f'Individual TWAS ({label})', title=phen, \
               out_file='example_run/figures/manhattan.png')
miami_plot(s, i, top_label=f'Summary TWAS ({label})', bottom_label=f'Individual TWAS ({label})', \
           title=phen, out_file='example_run/figures/miami.png')

print('\n> Done. Plots are in example_run/figures')
