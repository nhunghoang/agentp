'''
Plots for comparing gene-level association results (e.g., two TWAS, or TWAS vs a
gene-level GWAS such as MAGMA): a Miami plot and a -log10 p scatter plot.

Inputs are pandas Series of p-values indexed by Ensembl gene ID, e.g. columns of the
table from collect_twas(), or load_pvalues() for any other gene-level results file.
Each side is Bonferroni-thresholded on its own number of tested genes.

matplotlib is imported inside the plotting functions, so it is only needed for plots.
'''

import os
import sys
import numpy as np
import pandas as pd

from .Static import TAB, MSG
from .Static import Paths

## colors as in the AGENT-P paper figures
TOP_COLOR = '#984ea3'       ## significant in the top / x-axis results (purple)
BOTTOM_COLOR = '#ff7f00'    ## significant in the bottom / y-axis results (orange)
BOTH_COLOR = '#79502e'      ## significant in both (brown)
GRAYS = ['#c8c8c8', '#a7a7a7']

'''
Read a gene-level results file into a Series of p-values indexed by gene
(Ensembl IDs, version suffixes removed). Defaults read an AGENT-P TWAS output file.
'''
def load_pvalues(results_file, gene_col='gene', pval_col='pvalue', sep='\t'):
    if not os.path.exists(results_file):
        msg = f'{MSG}PLOT ERROR: Results file "{results_file}" does not exist'
        sys.exit(msg)
    df = pd.read_table(results_file, sep=sep, usecols=[gene_col, pval_col], dtype={gene_col: str})
    return pd.Series(df[pval_col].values, index=df[gene_col].values, name=pval_col)

'''
Clean a p-value Series: version suffixes removed from gene IDs, NaNs and duplicates dropped.
'''
def _clean(pvals, name):
    if not isinstance(pvals, pd.Series):
        msg = f'{MSG}PLOT ERROR: "{name}" must be a pandas Series of p-values indexed by gene'
        sys.exit(msg)
    pvals = pvals.dropna().astype(float)
    pvals.index = pd.Index(pvals.index.astype(str)).str.split('.').str[0]
    pvals = pvals[~pvals.index.duplicated()]
    if len(pvals) == 0:
        msg = f'{MSG}PLOT ERROR: "{name}" has no p-values'
        sys.exit(msg)
    return pvals

def _bonferroni(pvals, alpha):
    return -np.log10(alpha / len(pvals))

def _logp(pvals):
    ## p-values of exactly 0 are shown at the smallest positive double
    return -np.log10(pvals.clip(lower=np.finfo(float).tiny))

'''
Genome-wide x positions: chromosome offset + gene start (autosomes only).
'''
def _gene_positions():
    pos = pd.read_csv(Paths.GENE_POSITIONS, usecols=['chr', 'start', 'end', 'gene_id'], dtype={'chr': str})
    pos = pos[pos['chr'].isin([str(c) for c in range(1, 23)])].copy()
    pos['chr'] = pos['chr'].astype(int)
    size = pos.groupby('chr')['end'].max()
    offset = size.cumsum().shift(fill_value=0)
    pos['x'] = pos['start'] + pos['chr'].map(offset)
    centers = (offset + size / 2).to_dict()
    return pos.set_index('gene_id')[['chr', 'x']], centers

'''
Draw one set of gene p-values by genomic position (above zero if sign = 1, below if -1):
non-significant genes in alternating grays by chromosome, significant genes in color,
and a dashed Bonferroni line. Returns the largest -log10 p drawn (or the threshold).
'''
def _draw_genes(ax, pvals, sign, thr, color, label, pos):
    df = pd.DataFrame({'logp': _logp(pvals)}).join(pos, how='inner')
    n_missing = len(pvals) - len(df)
    if n_missing > 0:
        print(f'{TAB}{label}: {n_missing} of {len(pvals)} genes have no position and are not plotted')
    sig = df['logp'] >= thr
    gray = df['chr'].map(lambda c: GRAYS[c % 2])
    ax.scatter(df.loc[~sig, 'x'], sign * df.loc[~sig, 'logp'], c=gray[~sig], s=5, \
               alpha=0.5, edgecolor='none', rasterized=True)
    ax.scatter(df.loc[sig, 'x'], sign * df.loc[sig, 'logp'], s=6, facecolor=color, \
               edgecolor=(0, 0, 0, 0.45), linewidth=0.25, zorder=3)
    ax.axhline(sign * thr, color='k', lw=0.5, ls='dashed')
    return max(df['logp'].max() if len(df) else 0, thr)

'''
Chromosome ticks and labels on the x axis.
'''
def _chrom_axis(ax, pos, centers):
    ax.set_xlim(0, pos['x'].max())
    ax.set_xticks([centers[c] for c in range(1, 23)])
    ax.set_xticklabels([str(c) if (c <= 10) or (c % 2 == 0) else '' for c in range(1, 23)], fontsize=7)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel(r'$-\log_{10} p$')

'''
Miami plot: -log10 p of the top results above zero and of the bottom results below zero,
by genomic position. Dashed lines mark each side's Bonferroni threshold.

top, bottom:  Series of p-values indexed by gene (e.g., table['sTWAS_ACAT'])
top_label, bottom_label: labels shown in the plot corners
alpha:        family-wise error rate for the Bonferroni thresholds (default 0.05)
out_file:     save the figure here (e.g., .png or .pdf)
ax:           draw into an existing matplotlib axes instead of a new figure
Returns the matplotlib axes.
'''
def miami_plot(top, bottom, top_label='TWAS', bottom_label='GWAS', title=None, \
               alpha=0.05, out_file=None, ax=None):
    import matplotlib
    if (ax is None) and (out_file is not None):
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    top, bottom = _clean(top, top_label), _clean(bottom, bottom_label)
    thr_top, thr_bot = _bonferroni(top, alpha), _bonferroni(bottom, alpha)
    pos, centers = _gene_positions()

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 4))
    ymax = max(_draw_genes(ax, top, 1, thr_top, TOP_COLOR, top_label, pos), \
               _draw_genes(ax, bottom, -1, thr_bot, BOTTOM_COLOR, bottom_label, pos))

    ax.axhline(0, color='k', lw=0.3)
    lim = ymax * 1.08 + 1
    ax.set_ylim(-lim, lim)
    ticks = ax.get_yticks()
    ax.set_yticks(ticks)
    ax.set_yticklabels([f'{abs(t):g}' for t in ticks])
    ax.set_ylim(-lim, lim)
    _chrom_axis(ax, pos, centers)
    ax.text(0.01, 0.97, top_label, transform=ax.transAxes, va='top', ha='left', fontsize=8)
    ax.text(0.01, 0.03, bottom_label, transform=ax.transAxes, va='bottom', ha='left', fontsize=8)
    if title:
        ax.set_title(title, fontweight='bold')

    if out_file is not None:
        ax.figure.savefig(out_file, dpi=300, bbox_inches='tight')
        print(f'{MSG}Saved Miami plot to {out_file}')
    return ax

'''
Manhattan plot: -log10 p of one set of gene-level results by genomic position, with a
dashed line at the Bonferroni threshold (based on the number of genes tested).

pvals:     Series of p-values indexed by gene (e.g., table['iTWAS_ACAT'])
label:     label shown in the top-left corner
alpha, out_file, ax: as in miami_plot
Returns the matplotlib axes.
'''
def manhattan_plot(pvals, label='TWAS', title=None, alpha=0.05, out_file=None, ax=None):
    import matplotlib
    if (ax is None) and (out_file is not None):
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    pvals = _clean(pvals, label)
    thr = _bonferroni(pvals, alpha)
    pos, centers = _gene_positions()

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 3))
    ymax = _draw_genes(ax, pvals, 1, thr, TOP_COLOR, label, pos)
    ax.set_ylim(0, ymax * 1.08 + 1)
    _chrom_axis(ax, pos, centers)
    ax.text(0.01, 0.97, label, transform=ax.transAxes, va='top', ha='left', fontsize=8)
    if title:
        ax.set_title(title, fontweight='bold')

    if out_file is not None:
        ax.figure.savefig(out_file, dpi=300, bbox_inches='tight')
        print(f'{MSG}Saved Manhattan plot to {out_file}')
    return ax

'''
Scatter plot of -log10 p for genes tested in both results. Dashed lines mark each side's
Bonferroni threshold (based on all genes tested on that side). Colors: significant only on
x (purple), only on y (orange), in both (brown). The Pearson correlation of -log10 p is
shown in the plot.

x, y:      Series of p-values indexed by gene
x_label, y_label: axis labels
Returns the matplotlib axes.
'''
def scatter_plot(x, y, x_label='TWAS', y_label='GWAS', title=None, \
                 alpha=0.05, out_file=None, ax=None):
    import matplotlib
    if (ax is None) and (out_file is not None):
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    x, y = _clean(x, x_label), _clean(y, y_label)
    thr_x, thr_y = _bonferroni(x, alpha), _bonferroni(y, alpha)
    genes = x.index.intersection(y.index)
    if len(genes) < 3:
        msg = f'{MSG}PLOT ERROR: only {len(genes)} genes are tested in both results'
        sys.exit(msg)
    lx, ly = _logp(x[genes]), _logp(y[genes])
    sx, sy = lx >= thr_x, ly >= thr_y
    r = np.corrcoef(lx, ly)[0, 1]

    if ax is None:
        fig, ax = plt.subplots(figsize=(4, 4))
    ns = ~sx & ~sy
    ax.scatter(lx[ns], ly[ns], c=GRAYS[0], s=4, alpha=0.4, edgecolor='none', rasterized=True)
    for mask, color in [(sx & ~sy, TOP_COLOR), (~sx & sy, BOTTOM_COLOR), (sx & sy, BOTH_COLOR)]:
        ax.scatter(lx[mask], ly[mask], s=7, facecolor=color, edgecolor=(0, 0, 0, 0.45), \
                   linewidth=0.25, zorder=3)
    ax.axvline(thr_x, color='k', lw=0.5, ls='dashed')
    ax.axhline(thr_y, color='k', lw=0.5, ls='dashed')

    lim = max(lx.max(), ly.max(), thr_x, thr_y) * 1.05
    pad = 0.012 * lim
    ax.set_xlim(-pad, lim)
    ax.set_ylim(-pad, lim)
    ax.set_aspect('equal')
    ax.set_xlabel(rf'$-\log_{{10}} p$ ({x_label})')
    ax.set_ylabel(rf'$-\log_{{10}} p$ ({y_label})')
    ## white box so the text stays readable over points and threshold lines
    ax.text(0.04, 0.96, f'r = {r:.2f}', \
            transform=ax.transAxes, va='top', ha='left', fontsize=8, zorder=6, \
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='none', alpha=0.9))
    ax.set_title(title if title else f'{len(genes)} genes tested in both', fontsize=9)

    if out_file is not None:
        ax.figure.savefig(out_file, dpi=300, bbox_inches='tight')
        print(f'{MSG}Saved scatter plot to {out_file}')
    return ax
