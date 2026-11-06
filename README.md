# AGENT-P

**AGENT-P** (Association of Gene Expression to Neuroimaging Traits Pipeline) is a Python package for **transcriptome-wide association studies (TWAS)** of brain traits.

It predicts how much each gene is expressed in the brain from a person's genotypes (*genetically regulated gene expression*, or GREx), and tests whether this predicted expression is associated with a trait, such as the volume of a brain region. AGENT-P includes ready-to-use expression models for **8 brain regions** from **4 TWAS methods**, and runs TWAS in two ways:

|                               | Input                                                          | AGENT-P class      |
| ----------------------------- | -------------------------------------------------------------- | ------------------ |
| **Summary-level TWAS**  | GWAS summary statistics of the trait                           | `SummaryTWAS`    |
| **Genotype-level TWAS** | genotypes, covariates, and trait values of individual subjects | `IndividualTWAS` |

AGENT-P can also combine the results of all models into one p-value per gene (ACAT), make Manhattan, Miami, and scatter plots, and run a GWAS with REGENIE.

**Try it online, without installing anything:** [AGENT-P example in Google Colab](https://colab.research.google.com/github/nhunghoang/agentp/blob/agentp-update-20261105/example/colab_example.ipynb)

### Contents

1. [Install](#1-install)
2. [Download the models](#2-download-the-models)
3. [Run the example](#3-run-the-example)
4. [Run your own analysis](#4-run-your-own-analysis)
5. [Results and plots](#5-results-and-plots)
6. [Run all models and combine them (ACAT)](#6-run-all-models-and-combine-them-acat)
7. [Available models](#7-available-models)
8. [More options](#8-more-options)
9. [Troubleshooting](#9-troubleshooting)

---

## 1. Install

AGENT-P runs on **Linux** and **macOS**. The easiest way to install it is with [conda](https://docs.conda.io/en/latest/miniconda.html) (if you do not have conda, install [Miniforge](https://github.com/conda-forge/miniforge) first). Open a terminal and run:

```bash
git clone https://github.com/nhunghoang/agentp.git
cd agentp
conda env create -f environment.yml      # takes 5-15 minutes
conda activate agentp
pip install --no-deps -e .
R CMD INSTALL agentp/external/plink2R/plink2R
```

This installs AGENT-P with everything it needs: Python packages, the genetics tools `plink2`, `bgenix`, and `regenie`, and R for the FUSION models. The commands in steps 2 and 3 are run from this `agentp` folder.

**On a Mac with an Apple chip (M1-M4)**, replace the `conda env create` line with the one below (the genetics tools are only available for Intel Macs, and run through Rosetta), and run `conda config --env --set subdir osx-64` after `conda activate agentp`:

```bash
CONDA_SUBDIR=osx-64 conda env create -f environment.yml
```

Each time you open a new terminal, run `conda activate agentp` before using AGENT-P.

---

## 2. Download the models

The expression models (about 2 GB) are not included in the GitHub repository. Download [`models.tar.gz`](https://drive.google.com/file/d/1llysnSFOut3E_g1omlgQfi-N5fptuuLa/view?usp=sharing) and extract it into the `agentp` folder inside the repository:

```bash
cd agentp                                # the "agentp" folder inside the downloaded repository
tar -xzf ~/Downloads/models.tar.gz       # change the path if you saved the file elsewhere
cd ..
```

The models are now in `agentp/models`. You can skip this step if you only want to run the example in step 3, which uses a small set of example models.

---

## 3. Run the example

The example uses a small dataset of 2,000 subjects that comes with AGENT-P. It runs both TWAS for one brain region (caudate) and one trait (amygdala volume):

```bash
cd example
python run_example.py
```

When it finishes, open the plots in `example/example_run/figures/`: a Manhattan plot (`manhattan.png`) and a Miami plot (`miami.png`). The [Colab notebook](https://colab.research.google.com/github/nhunghoang/agentp/blob/agentp-update-20261105/example/colab_example.ipynb) walks through the same example step by step.

---

## 4. Run your own analysis

Write the steps below in a Python script (for example `my_analysis.py`), and run it with `python my_analysis.py`.

### Step 1: Create a project

A project is a folder that keeps the inputs and results of one cohort together.

```python
from agentp import Project
project = Project('my_project')     # creates the folder "my_project"
```

### Step 2a: Summary-level TWAS (from a GWAS)

```python
from agentp import SummaryTWAS

stwas = SummaryTWAS(project, 'JTI', 'caudate')                  # model method and brain region
stwas.add_gwas('putamen_volume', 'putamen_volume_gwas.txt')     # trait name and GWAS file
stwas.run_twas('putamen_volume')
```

**The GWAS file** needs one row per SNP, with these columns: SNP ID (rsID), the two alleles, and the effect size (beta) with its standard error or a p-value. Output files of REGENIE (and of AGENT-P's own GWAS) are read without any extra settings. For other files, tell AGENT-P the column names:

```python
stwas.add_gwas('putamen_volume', 'putamen_volume_gwas.txt', sep='\t',
               snp_col='rsid', alt_col='a2', ref_col='a1',
               beta_col='beta', beta_se_col='se', pval_col='pval')
```

* `alt_col` is the **effect allele**: the allele that the beta refers to. `ref_col` is the other allele.
* If your file has a z-score or odds ratio instead of a beta, use `zscore_col` or `or_col`.

### Step 2b: Genotype-level TWAS (from individual-level data)

First add the subject data to the project (only once per project):

```python
project.set_subjects('subjects.txt')             # subject IDs
project.add_covariates('covariates.txt')         # e.g. age, sex, genetic principal components
project.add_phenotypes('phenotypes.txt')         # one or more traits
project.add_genotypes('genotypes', 'chr*.bgen')  # folder, and file names with * for the chromosome number
```

Then predict gene expression and test the traits:

```python
from agentp import IndividualTWAS

itwas = IndividualTWAS(project, 'JTI', 'caudate')
itwas.predict_grex()                         # predict gene expression (once per method and region)
itwas.run_twas('putamen_volume')             # or a list: ['putamen_volume', 'caudate_volume']
```

**Input files** (see `example/example_data` for examples):

| File        | Format                                                                                                               |
| ----------- | -------------------------------------------------------------------------------------------------------------------- |
| Subject IDs | one ID per line, no header                                                                                           |
| Covariates  | a table with an `id` column and one column per covariate (numbers only, e.g. sex as 0/1); tab-separated by default |
| Phenotypes  | a table with an `id` column and one column per trait; tab-separated by default                                     |
| Genotypes   | one BGEN file per chromosome (1-22), each with a matching `.sample` file, and SNP IDs as rsIDs                     |

* The subject IDs must be the same in all files (in the `.sample` files, both ID columns).
* To use only some columns, add `covar_cols=[...]` or `phen_cols=[...]`.
* AGENT-P cleans the genotypes before use: it keeps SNPs with less than 5% missing calls, a minor allele frequency of at least 1%, and a Hardy-Weinberg equilibrium p-value above 1e-5.
* Genotypes in VCF format can be converted to BGEN with `plink2` (see [Troubleshooting](#9-troubleshooting)).

---

## 5. Results and plots

Each TWAS saves one table per trait (tab-separated, one row per gene):

| TWAS           | Results file                                         |
| -------------- | ---------------------------------------------------- |
| Summary-level  |  `my_project/sTWAS/JTI_caudate/putamen_volume.csv` |
| Genotype-level | `my_project/iTWAS/JTI_caudate/putamen_volume.csv`  |

The most important columns are `gene` (Ensembl gene ID), `zscore` (direction and strength of the association: positive means that higher expression goes with a larger trait value), and `pvalue`. A common threshold for significance is a Bonferroni correction: p < 0.05 / (number of genes tested).

**Plots** (saved as image files):

```python
from agentp import load_pvalues, manhattan_plot, miami_plot, scatter_plot

s = load_pvalues('my_project/sTWAS/JTI_caudate/putamen_volume.csv')
g = load_pvalues('my_project/iTWAS/JTI_caudate/putamen_volume.csv')

manhattan_plot(g, label='Genotype-level TWAS', out_file='manhattan.png')        # one set of results
miami_plot(s, g, top_label='Summary-level TWAS', bottom_label='Genotype-level TWAS',
           out_file='miami.png')                                                # two sets, mirrored
scatter_plot(s, g, x_label='Summary-level TWAS', y_label='Genotype-level TWAS',
             out_file='scatter.png')                                            # two sets, gene by gene
```

Significant genes (Bonferroni) are colored, and the dashed lines show the significance thresholds.

---

## 6. Run all models and combine them (ACAT)

`run_all_twas` runs TWAS with every model (4 methods x 8 brain regions = 32 model sets), then combines the p-values of each gene into one p-value with ACAT (Aggregated Cauchy Association Test):

```python
from agentp import run_all_twas

# summary-level TWAS
run_all_twas(project, 'putamen_volume', twas_type='stwas',
             gwas_files={'putamen_volume': 'putamen_volume_gwas.txt'})

# genotype-level TWAS (after adding the subject data, step 2b)
run_all_twas(project, 'putamen_volume', twas_type='itwas')
```

The combined results are saved to `my_project/sTWAS/ACAT/putamen_volume.csv` (or `iTWAS/ACAT/`). This table has one p-value column per model set, `n_models` (how many models tested the gene), and the combined p-value (`sTWAS_ACAT` or `iTWAS_ACAT`). Genes tested by fewer than 4 models get no combined p-value.

If your GWAS file is not in REGENIE format, give the column names once for all files: `gwas_params={'sep': '\t', 'snp_col': 'rsid', 'alt_col': 'a2', ...}`. To combine results that you already ran one model at a time, use `collect_twas(project, 'putamen_volume')`.

---

## 7. Available models

| Method                        | Code    |  | Brain region                   | Code                |
| ----------------------------- | ------- | - | ------------------------------ | ------------------- |
| PrediXcan                     | `PDX` |  | Dorsolateral prefrontal cortex | `dlpfc`           |
| UTMOST                        | `UTM` |  | Anterior cingulate cortex      | `ant-cingulate`   |
| JTI (Joint-Tissue Imputation) | `JTI` |  | Amygdala                       | `amygdala`        |
| FUSION                        | `FUS` |  | Hippocampus                    | `hippocampus`     |
|                               |         |  | Caudate                        | `caudate`         |
|                               |         |  | Putamen                        | `putamen`         |
|                               |         |  | Nucleus accumbens              | `nuc-accumbens`   |
|                               |         |  | Cerebellar hemisphere          | `cerebellar-hemi` |

Any method can be combined with any region, e.g. `SummaryTWAS(project, 'FUS', 'hippocampus')`.

**Adding other models.** Models in PredictDB format (a `.db` weights file and a `.txt.gz` covariance file, e.g. other GTEx tissues) can be added and then used like the built-in ones:

```python
project.add_model('PDX', 'hypothalamus', weights='en_Brain_Hypothalamus.db',
                  covariance='en_Brain_Hypothalamus.txt.gz')
stwas = SummaryTWAS(project, 'PDX', 'hypothalamus')
```

---

## 8. More options

* **Run a GWAS** (needs the subject data from step 2b):
  ```python
  from agentp import GWAS
  GWAS(project).run(['putamen_volume', 'caudate_volume'])    # saved to my_project/GWAS/
  ```
* **Run again and replace old results:** add `overwrite=True`, e.g. `stwas.run_twas('putamen_volume', overwrite=True)`. Without it, AGENT-P skips results that already exist.
* **Several traits in parallel:** `run_twas([...], num_threads=4)` runs 4 traits at the same time.
* **Test the installation** with the full models: `python tests/test_pipeline.py` (add `--fusion` to include FUSION).

---

## 9. Troubleshooting

* **`git clone` fails with an error about "LFS" or "smudge filter".** Clone without the large files (the models come from step 2 anyway): `GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/nhunghoang/agentp.git`.
* **`conda env create` is very slow or runs out of memory.** Use [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) with the same file: `micromamba create -f environment.yml`.
* **`numpy` or `h5py` errors after installing.** AGENT-P was installed without `--no-deps`, and pip replaced conda's packages. Remove the environment (`conda env remove -n agentp`) and install again as in step 1.
* **On a Mac: "incompatible architecture".** The environment was created without `CONDA_SUBDIR=osx-64`. Remove it and create it again as in step 1. The first run on a Mac can take a few extra minutes while Rosetta translates the programs.
* **My genotypes are in VCF format.** Convert them to BGEN, one file per chromosome:
  ```bash
  for c in $(seq 1 22); do
      plink2 --vcf chr${c}.dose.vcf.gz dosage=DS --double-id --export bgen-1.2 bits=8 --out genotypes/chr${c}
  done
  ```

  SNP IDs must be rsIDs. If your VCF uses other IDs (e.g. `chr1:12345:A:G`), add rsIDs first (e.g. with `bcftools annotate` and a dbSNP file).
* **Very few genes are tested.** The SNP IDs in your GWAS or genotypes do not match the models: check that they are rsIDs.

---

### License

MIT (see `LICENSE`).
