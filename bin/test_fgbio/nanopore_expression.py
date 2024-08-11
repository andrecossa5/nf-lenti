"""
Correlation site-coverage and gene expression SCM-seq
"""

import os
import numpy as np
import pandas as pd
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


##


# Explore
sample = 'sAML1'
path_results = f'/Users/IEO5505/Desktop/example_mito/results/test_nanopore_{sample}'
path_counts = os.path.join(path_results, 'counts_table.csv') 
path_expr = os.path.join(f'/Users/IEO5505/Desktop/example_mito/scratch/nuclear_SNVs_expression_10x.csv')

allele_counts = pd.read_csv(path_counts)
expression = pd.read_csv(path_expr, index_col=0)
expression.index = expression.index.map(lambda x: x.split('-')[0])

df = allele_counts.query('MUT>0 or WT>0')
MUT = df.pivot(index='cell', columns='gene', values='MUT').fillna(0)
WT = df.pivot(index='cell', columns='gene', values='WT').fillna(0)
coverage = MUT+WT

cells = list(set(coverage.index) & set(expression.index))
df_ = coverage.loc[cells].mean(axis=0).to_frame('coverage_SCM').join(expression.loc[cells].mean(axis=0).to_frame('expression_10x'))

##

# All cells, aggregated with means
fig, ax = plt.subplots(figsize=(5,5))
sns.regplot(data=df_, x='expression_10x', y='coverage_SCM', scatter=False)
ax.plot(df_['expression_10x'], df_['coverage_SCM'], 'ko', markersize=10)
format_ax(ax=ax, title=f'Corr: {df_.corr().iloc[1,0]:.2f}', xlabel='n UMIs (expression, 10x)', ylabel='n UMIs (per target base, SCM-seq)')
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'corr_SCM_10x_all_cells_aggregated.png'), dpi=1000)

##

# All calls
df_ = (
    df.assign(nUMIs_SCM=lambda x: x['MUT']+x['WT']+x['MIS']).drop(columns=['MUT', 'WT', 'MIS'])
    .merge( 
        expression.reset_index().rename(columns={'index':'cell'})
        .melt(id_vars='cell', value_name='nUMIs_10x', var_name='gene'),
        on=['cell', 'gene']
    )
)

fig, ax = plt.subplots(figsize=(5,5))
sns.regplot(data=df_, x='nUMIs_10x', y='nUMIs_SCM', scatter=False)
ax.plot(df_['nUMIs_10x'], df_['nUMIs_SCM'], 'ko', markersize=2)
format_ax(ax=ax, title=f'Corr: {df_[["nUMIs_SCM", "nUMIs_10x"]].corr().iloc[1,0]:.2f}', 
          xlabel='n UMIs (expression, 10x)', ylabel='n UMIs (per target base, SCM-seq)')
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'corr_SCM_10x_all_cells.png'), dpi=1000)

plt.show()

##

# Cell type
path_meta = '/Users/IEO5505/Desktop/AML_clonal_reconstruction/data/meta/cells_meta.csv'
meta = pd.read_csv(path_meta, index_col=0)
df_ = df_.merge(meta[['aggregated_ct']].reset_index().rename(columns={'index':'cell'}), on='cell')


np.sum(coverage>0)
(WT['CSF3R']>0).sum()

L = []

for gene in df_['gene'].unique():
    df_gene = df_.query('gene==@gene')
    # tenx = (df_gene['nUMIs_SCM']>0).sum()
    # scm = (df_gene['nUMIs_10x']>0).sum()
    # print(f'{gene}: tenx {tenx}, scm {scm}')
    both = np.sum((df_gene['nUMIs_SCM']>0) & (df_gene['nUMIs_SCM']>0))
    print(f'{gene}: {both}')



df = pd.DataFrame(L, columns=['nUMIs_SCM', 'nUMIs_10x'], index=df_['gene'].unique())
df['nUMIs_SCM'] / df['nUMIs_10x']







    fig, axs = plt.subplots(1,2,figsize=(7,4.5))

    ax = axs[0]
    box(df_gene, 'aggregated_ct', 'nUMIs_10x', ax=ax)
    mean = df_gene['nUMIs_10x'].mean()
    std = df_gene['nUMIs_10x'].std()
    format_ax(ax=ax, rotx=90, ylabel='10x nUMIs', title=f'{mean:.2f}(+-{std:.2f})')
    ax.set_ylim((-2,20))

    ax = axs[1]
    box(df_gene, 'aggregated_ct', 'nUMIs_SCM', ax=ax)
    mean = df_gene['nUMIs_SCM'].mean()
    std = df_gene['nUMIs_SCM'].std()
    format_ax(ax=ax, rotx=90, ylabel='SCM-seq nUMIs', title=f'{mean:.2f}(+-{std:.2f})')
    ax.set_ylim((-4,50))

    fig.suptitle(f'{gene}: n cells {df_gene.shape[0]}')
    fig.tight_layout()
    fig.savefig(os.path.join(path_results, f'{gene}_nUMIs_by_celltype.png'), dpi=1000)


##