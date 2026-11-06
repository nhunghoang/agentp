'''
End-to-end test of AGENT-P on the example data (example/example_data):

  1) project set up: subjects, covariates, phenotypes, genotype QC (plink2 + bgenix)
  2) GWAS with REGENIE (skipped if regenie is not on PATH; the example GWAS file is used)
  3) summary TWAS with every model set of the chosen approaches + ACAT meta-analysis
  4) individual TWAS with the same model sets + ACAT meta-analysis
  5) Miami and scatter plots of the summary vs individual ACAT p-values, and a Manhattan plot
     of the individual ACAT p-values

Each step is checked and reported as PASS / FAIL / SKIP; the script exits with code 1 if
any check fails. Results are written to tests/test_project (removed first, unless --keep).

usage: python test_pipeline.py [--fusion] [--keep]
  --fusion  also run FUSION models (needs Rscript; slower, FUSION runs per chromosome)
  --keep    keep an existing test_project and reuse finished steps
'''

import os
import sys
import shutil
import argparse
import pandas as pd

from agentp import Project, GWAS, run_all_twas, manhattan_plot, miami_plot, scatter_plot
from agentp.Static import Models, MSG, TAB

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), 'example', 'example_data')
PROJ = f'{HERE}/test_project'
PHEN = 'amygdala_volume'                              # in example_data/volumes.csv
GWAS_FILE = f'{DATA}/vol_mean_amygdala.regenie'       # example GWAS of amygdala volume

parser = argparse.ArgumentParser()
parser.add_argument('--fusion', action='store_true')
parser.add_argument('--keep', action='store_true')
args = parser.parse_args()
approaches = ['PDX', 'UTM', 'JTI'] + (['FUS'] if args.fusion else [])
n_sets = len([1 for (app, reg) in Models.model_sets() if app in approaches])

results = []
def check(name, ok, detail=''):
    status = 'PASS' if ok else 'FAIL'
    results.append((name, status, detail))
    print(f'{MSG}[{status}] {name}' + (f' -- {detail}' if detail else ''))

def skip(name, detail):
    results.append((name, 'SKIP', detail))
    print(f'{MSG}[SKIP] {name} -- {detail}')

## external tools
tools = ['plink2', 'bgenix'] + (['Rscript'] if args.fusion else [])
missing = [t for t in tools if shutil.which(t) is None]
check('external tools on PATH', len(missing) == 0, f'missing: {", ".join(missing)}' if missing else ', '.join(tools))
if missing:
    sys.exit(1)

if os.path.exists(PROJ) and not args.keep:
    shutil.rmtree(PROJ)
os.chdir(HERE)

## 1) project set up and genotype QC
project = Project(PROJ)
project.set_subjects(f'{DATA}/subjects.txt')
project.add_covariates(f'{DATA}/covariates.csv')
project.add_phenotypes(f'{DATA}/volumes.csv')
project.add_genotypes(f'{DATA}/genotypes', 'test_c*.bgen')
gdir = project.get_genotypes_dir()
n_bgen = sum(os.path.exists(f'{gdir}/c{c}.bgen') and os.path.exists(f'{gdir}/c{c}.bgen.bgi') for c in range(1, 23))
check('genotype QC', n_bgen == 22, f'{n_bgen}/22 chromosomes with .bgen + .bgi')

## 2) GWAS (its output is the input of summary TWAS; otherwise the example GWAS file is used)
if shutil.which('regenie') is None:
    skip('GWAS (REGENIE)', 'regenie is not on PATH; using the example GWAS file for summary TWAS')
else:
    GWAS(project).run(PHEN)
    gwas_out = f'{project.output_paths["gwas"]}/{PHEN}.regenie'
    n_snps = len(pd.read_csv(gwas_out, usecols=['SNP'])) if os.path.exists(gwas_out) else 0
    check('GWAS (REGENIE)', n_snps > 0, f'{n_snps} SNPs in {gwas_out}')
    if n_snps > 0:
        GWAS_FILE = gwas_out

## 3) summary TWAS + ACAT
print(f'{MSG}Summary TWAS input GWAS: {GWAS_FILE}')
stables = run_all_twas(project, PHEN, twas_type='stwas', gwas_files={PHEN: GWAS_FILE}, approaches=approaches)
st = stables[PHEN]
pcols = [c for c in st.columns if c.startswith('sTWAS_') and c != 'sTWAS_ACAT']
n_acat = int(st['sTWAS_ACAT'].notna().sum())
check('summary TWAS', len(pcols) == n_sets, f'{len(pcols)}/{n_sets} model sets, {len(st)} genes')
check('summary ACAT', n_acat > 0, f'{n_acat} genes with an ACAT p-value')

## 4) individual TWAS + ACAT
itables = run_all_twas(project, PHEN, twas_type='itwas', approaches=approaches)
it = itables[PHEN]
pcols = [c for c in it.columns if c.startswith('iTWAS_') and c != 'iTWAS_ACAT']
n_acat = int(it['iTWAS_ACAT'].notna().sum())
check('individual TWAS', len(pcols) == n_sets, f'{len(pcols)}/{n_sets} model sets, {len(it)} genes')
check('individual ACAT', n_acat > 0, f'{n_acat} genes with an ACAT p-value')

## 5) plots
fig_dir = f'{PROJ}/figures'
os.makedirs(fig_dir, exist_ok=True)
for name, func in [('miami', miami_plot), ('scatter', scatter_plot)]:
    out = f'{fig_dir}/{name}_{PHEN}.png'
    func(st['sTWAS_ACAT'], it['iTWAS_ACAT'], 'Summary TWAS', 'Individual TWAS', \
         title=f'{PHEN} (ACAT)' if name == 'miami' else None, out_file=out)
    check(f'{name} plot', os.path.exists(out) and os.path.getsize(out) > 0, out)
out = f'{fig_dir}/manhattan_{PHEN}.png'
manhattan_plot(it['iTWAS_ACAT'], 'Individual TWAS (ACAT)', title=PHEN, out_file=out)
check('manhattan plot', os.path.exists(out) and os.path.getsize(out) > 0, out)

## summary
print(f'\n{"=" * 60}\nAGENT-P end-to-end test ({", ".join(approaches)} models)')
for name, status, detail in results:
    print(f'{TAB}{status:4s}  {name}' + (f' -- {detail}' if detail else ''))
failed = [r for r in results if r[1] == 'FAIL']
print(f'{"=" * 60}\n{len(failed)} failed, {sum(r[1] == "PASS" for r in results)} passed, ' \
      f'{sum(r[1] == "SKIP" for r in results)} skipped')
sys.exit(1 if failed else 0)
