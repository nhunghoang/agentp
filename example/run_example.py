from agentp import Project, SummaryTWAS, IndividualTWAS, GWAS

## init project 
main_project_folder = 'example_project'
project = Project(main_project_folder)

########################################################################

## init a summary twas that uses JTI models trained on hippocampus genes
stwas = SummaryTWAS(project, 'JTI', 'hippocampus') 

## apply these models to a GWAS on amygdala volumes 
stwas.add_gwas('amygdala_volume', 'example_data/vol_mean_amygdala.regenie') 

## run the summary TWAS 
stwas.run_twas('amygdala_volume') 

## TWAS results are saved to
## agentp_test/sTWAS/JTI_hippocampus/amygdala_volume.csv

########################################################################

## init a individual twas that uses FUSION models trained on caudate genes
gtwas = IndividualTWAS(project, 'FUS', 'caudate')

## add subject data to the project 
project.set_subjects('example_data/subjects.txt') 
project.add_covariates('example_data/covariates.csv')
project.add_phenotypes('example_data/volumes.csv')
project.add_genotypes('example_data/genotypes', 'test_c*.bgen')

## predict FUSION-caudate grex for this cohort 
gtwas.predict_grex() 

## run the individual TWAS on putamen volumes
gtwas.run_twas('putamen_volume')

## TWAS results are saved to
## agentp_test/gTWAS/FUS_caudate/putamen_volume.csv

########################################################################

## init a gwas 
gwas = GWAS(project) 

## no need to add the subject data again, since we already did for gTWAS
## run GWAS on nucleus hippocampus volumes 
gwas.run('hippocampus_volume')

## GWAS results are saved to
## agentp_test/GWAS/hippocampus_volume.regenie
