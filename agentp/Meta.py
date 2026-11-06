'''
Combine TWAS results across approaches and GREx regions.

collect_twas() merges the per-model TWAS outputs of one phenotype
(project/{sTWAS,iTWAS}/{approach}_{region}/{phen}.csv) into one table
with one p-value column per approach-region model set, and adds a
meta-analysis p-value per gene from the Aggregated Cauchy Association
Test (ACAT; Liu et al., 2019, AJHG).
'''

import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import cauchy

from .Static import TAB, MSG
from .Static import Models

'''
Aggregated Cauchy Association Test.
R source: https://github.com/yaowuliu/ACAT

Combines p-values with equal weights. NaNs are ignored; if fewer than
min_pvals p-values are available, the combined p-value is NaN.
Very small p-values (< 1e-15) use the asymptotic form 1/(p*pi) of
tan((0.5-p)*pi). The upper tail uses arctan(1/t)/pi for t > 0, which is
exact in double precision (1 - cauchy.cdf(t) rounds to 0 for strong
associations, and cauchy.sf(t) does too in scipy < 1.11 for t >= 1e16).
'''
def acat(pvals, min_pvals=4):
    pvals = np.asarray(pvals, dtype=float)
    pvals = pvals[~np.isnan(pvals)]
    if pvals.size < min_pvals:
        return np.nan

    if ((pvals < 0) | (pvals > 1)).any():
        msg = f'{MSG}ACAT ERROR: p-values must be between 0 and 1'
        sys.exit(msg)

    pset = np.zeros_like(pvals)
    small = (pvals < 1e-15)
    pset[~small] = np.tan((0.5 - pvals[~small]) * np.pi)
    pset[small] = 1 / pvals[small] / np.pi

    stat = np.mean(pset)
    if stat > 0:
        return np.arctan(1 / stat) / np.pi
    return cauchy.sf(stat)

'''
Gather TWAS p-values for one phenotype into a gene x model table.

project:    an instance of the Project class
phen:       phenotype name (as used in run_twas)
twas_type:  'stwas' (summary TWAS) or 'itwas' (individual TWAS)
approaches: TWAS approaches to include (default: all, including added models)
regions:    GREx regions to include (default: all, including added models)
min_pvals:  minimum number of available p-values for a gene's ACAT p-value
save:       write the table to project/{sTWAS,iTWAS}/ACAT/{phen}.csv

Columns: {sTWAS,iTWAS}_{approach}_{region} p-values, n_models (number of
available p-values), and {sTWAS,iTWAS}_ACAT. Gene IDs are Ensembl IDs
without version suffixes. Model sets with no output for this phenotype
are skipped with a warning.
'''
def collect_twas(project, phen, twas_type='stwas', approaches=None, \
                 regions=None, min_pvals=4, save=True):

    if twas_type not in ['stwas', 'itwas']:
        msg = f'{MSG}ACAT ERROR: twas_type must be "stwas" or "itwas"'
        sys.exit(msg)

    ## built-in and added model sets, optionally restricted to some approaches/regions
    model_sets = [(app, reg) for (app, reg) in Models.model_sets() \
                  if (approaches is None or app in approaches) and (regions is None or reg in regions)]
    prefix = {'stwas': 'sTWAS', 'itwas': 'iTWAS'}[twas_type]
    twas_top = project.output_paths[twas_type]

    print(f'{MSG}Combining {prefix} results for "{phen}"...')
    table = pd.DataFrame([], columns=['gene'])
    pcols = []
    for (app, reg) in model_sets:
        tfile = f'{twas_top}/{app}_{reg}/{phen}.csv'
        if not os.path.exists(tfile):
            print(f'{TAB}WARNING: no {prefix} results for {app} {reg} models (skipped)')
            continue

        data = pd.read_table(tfile, sep='\t', usecols=['gene', 'pvalue'], dtype={'gene': str})
        data['gene'] = data['gene'].str.split('.').str[0]
        data = data.dropna(subset=['pvalue']).drop_duplicates('gene')

        col = f'{prefix}_{app}_{reg}'
        data = data.rename(columns={'pvalue': col})
        table = table.merge(data, on='gene', how='outer')
        pcols.append(col)

    if len(pcols) == 0:
        msg = f'{MSG}ACAT ERROR: no {prefix} results found for "{phen}"'
        sys.exit(msg)

    table = table.set_index('gene')
    table['n_models'] = table[pcols].notna().sum(axis=1)
    table[f'{prefix}_ACAT'] = table[pcols].apply(lambda row: acat(row.values, min_pvals), axis=1)

    n_genes = table[f'{prefix}_ACAT'].notna().sum()
    print(f'{TAB}{len(pcols)} model sets, {len(table)} genes ({n_genes} with >= {min_pvals} p-values for ACAT)')

    if save:
        out_dir = f'{twas_top}/ACAT'
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        out_file = f'{out_dir}/{phen}.csv'
        table.to_csv(out_file)
        print(f'{TAB}Saved to {out_file}')

    return table
