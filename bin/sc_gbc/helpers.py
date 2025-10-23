#!/usr/bin/python

"""
Utils for custom cell_assignement. Brings together elements from LARRY_2020, CR_2023 and
other (pre-)pubblication works.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import poisson
from sklearn.metrics.pairwise import pairwise_distances
from itertools import chain
from plotting_utils._plotting_base import *
from plotting_utils._utils import Timer


##


def read_data(path_sc, sample=None):
    """
    Create a table of CBC-UMI-GBC combinations from single-cell data, after correcting 
    GBCs with a bulk reference.
    """

    sc_df = pd.read_csv(path_sc, sep='\t', dtype='str')
    d_rev = {'A':'T', 'G':'C', 'T':'A', 'C':'G', 'N':'N'}
    sc_df['GBC'] = sc_df['GBC'].map(lambda x: ''.join([ d_rev[x] for x in reversed(x) ]))

    return sc_df


##


def to_numeric(X):
    return np.select([X=='A', X=='T', X=='C', X=='G', X=='N'], [1,2,3,4,5], default=0)


##


def hamming(bc1, bc2): 
    return np.sum([x1 != x2 for x1,x2 in zip(bc1,bc2)])


##

 
def map_GBCs(sc_df, bulk=None, bulk_correction_treshold=3, ncores=8):
    """
    Correct sc GBC sequences with bulk-DNA reference. Outputs a sequence map: sc:bulk.
    """

    # Correct with bulk reference
    sc = sc_df['GBC'].value_counts()
    sc_numeric = to_numeric(np.vstack(sc.index.map(lambda x: np.array(list(x)))))
    bulk_numeric = to_numeric(np.vstack(bulk.index.map(lambda x: np.array(list(x)))))
    D = pairwise_distances(
        sc_numeric, bulk_numeric, metric='hamming', n_jobs=int(ncores)
    ) * sc_numeric.shape[1]

    # Build a correction dict for sc GBCs at hamming distance <= bulk_sc_treshold from a bulk one.
    d_corr = (
        sc.to_frame('read_count')
        .assign(
            correct_GBC=[ bulk.index[i] for i in D.argmin(axis=1) ],
            hamming=D.min(axis=1),
        )
    )
    d_corr = d_corr.query('hamming<=@bulk_correction_treshold')['correct_GBC'].to_dict()

    return d_corr


##


def count_UMIs(sc_df, gbc_col='GBC'):
    counts = (
        sc_df.groupby(['CBC', gbc_col, 'UMI']).size() # .groupby(['CBC', gbc_col])['UMI'].nunique() 
        .reset_index(name='count')
    )
    counts['log'] = np.log(counts['count'])/np.log(10) # LARRY log transformation

    return counts


##


def get_combos(counts, gbc_col='GBC_reference-free'):
    """
    Compute CBC-GBC stats from filtered UMIs.
    """
    #filtered = counts.loc[lambda x: x['status']=='Retain']
    df_combos = (
        #filtered
        counts
        .rename(columns={gbc_col:'GBC'})
        .groupby(['CBC', 'GBC'])['UMI'].nunique()
        .to_frame('umi').reset_index()
        .assign(
            max_ratio=lambda x: \
            x.groupby('CBC')['umi'].transform(lambda x: x/x.max()),
            normalized_abundance=lambda x: \
            x.groupby('CBC')['umi'].transform(lambda x: x/x.sum())
        )
    )
    return df_combos


##


def filter_and_pivot(df_combos, umi_treshold, p_treshold, max_ratio_treshold, normalized_abundance_treshold):
    """
    Filter a CBC-GBC-UMI table and return it in large format.
    """

    # p of observing higher counts than the observed ones, assuming UMI counts are Poisson distributed
    avg_gbc = df_combos.groupby('GBC')['umi'].median()
    p_poisson = []
    for i in range(df_combos.shape[0]):
        x = df_combos.iloc[i,:]
        gbc = x['GBC']
        obs = x['umi']
        exp = avg_gbc.loc[gbc]
        p = 1-poisson.cdf(obs-1, mu=exp)
        p_poisson.append(p)
    df_combos['p'] = p_poisson
    df_combos['logp'] = -np.log10(p_poisson)
    df_combos['logumi'] = np.log10(df_combos['umi'])

    ##

    # Supported and unsupported CBC-GBC combinations
    print('-----------------------------------')
    print('df_combos[umi]_median=',df_combos['umi'].median(), 'th=',umi_treshold)
    print('df_combos[max_ratio]_median=',df_combos['max_ratio'].median(), 'th=',max_ratio_treshold)
    print('df_combos[normalized_abundance]_median=',df_combos['normalized_abundance'].median(), 'th=',normalized_abundance_treshold)
    print('df_combos[p]_median=',df_combos['p'].median(), 'th=',p_treshold)
    print('-----------------------------------')
    test = (df_combos['umi'] >= umi_treshold) & \
        (df_combos['p'] <= p_treshold) & \
        (df_combos['max_ratio'] >= max_ratio_treshold) & \
        (df_combos['normalized_abundance'] >= normalized_abundance_treshold)
    df_combos['status'] = np.where(test, 'supported', 'unsupported')

    # Filter and pivot
    M = (
        df_combos
        .query('status=="supported"')
        .pivot_table(index='CBC', columns='GBC', values='umi')
    )
    M[M.isna()] = 0

    # Print n 1-GBC cells
    print('Estimated MOI:')
    print((M>0).sum(axis=1).value_counts())

    return M, df_combos


##


def get_clones(M):
    """
    Get the clone dataframe for each sample.
    """
    cells_with_unique_combos = M.apply(lambda x: frozenset(M.columns[x>0]), axis=1)
    cells_with_unique_combos = cells_with_unique_combos.map(lambda x: ';'.join(x))
    clones = (
        cells_with_unique_combos
        .reset_index(drop=True)
        .value_counts(normalize=False)
        .to_frame('n cells').reset_index()
        .rename(columns={'index':'GBC_set'})
        .assign(clone=lambda x: [ f'C{i}' for i in range(x.shape[0]) ])
        .set_index('clone')
    )
    clones['prevalence'] = clones['n cells'] / clones['n cells'].sum()

    return clones


##


def cell_assignment_workflow(
    path_sc, bulk_correction_treshold,
    umi_treshold, p_treshold, max_ratio_treshold, normalized_abundance_treshold,sample=None, 
    path_bulk=None, path_sample_map=None, sample_params=None
    ):
    """
    Complete clone calling and cell assignment workflow.
    Read bulk and single-cell data, correct sc GBCs with the bulk reference, obtain a df with CBC-GBC combinations,
    filter them and assign cells to their clonal labels.
    """
    import numpy as np
    import plotly.graph_objs as go  
    import matplotlib.ticker as mticker

    T = Timer()

    T.start()
    f = open('clone_calling_summary.txt', 'w')

    if sample_params is not None:
        umi_treshold = sample_params['umi_treshold']
        p_treshold = sample_params['p_treshold']
        max_ratio_treshold = sample_params['max_ratio_treshold']
        normalized_abundance_treshold = sample_params['normalized_abundance_treshold']

    f.write(f'# Custom clone calling and cell assignment workflow, sample {sample}: \n')
    f.write('\n')
    f.write('Input params:\n')
    f.write(f'  * bulk_correction_treshold: {bulk_correction_treshold}\n')
    f.write(f'  * umi_treshold: {umi_treshold}\n')
    f.write(f'  * p_treshold: {p_treshold}\n')
    f.write(f'  * max_ratio_treshold: {max_ratio_treshold}\n')
    f.write(f'  * normalized_abundance_treshold: {normalized_abundance_treshold}\n')
    f.write(f'\n')

    ##

    # Read and count
    sc_df = read_data(path_sc, sample=sample)

    # Optional: correction with bulk reference
    bulk_correction = True if os.path.exists(path_sample_map) and os.path.exists(path_bulk) else False

    if bulk_correction:
        bulk = pd.read_csv(os.path.join(path_bulk, 'summary', 'bulk_GBC_reference.csv'), index_col=0)
        sample_map = pd.read_csv(path_sample_map, index_col=0)
        if sample in sample_map.index:
            ref = sample_map.loc[sample, 'reference']
            bulk_df = bulk.query('sample==@ref')
            assert bulk.shape[0]>0
            print(f'Found bulk GBC sequences for the {sample} sample, from ref {ref}.')
            # Create correction map and map sc GBCs to corrected bulk ones
            bulk_map = map_GBCs(sc_df, bulk=bulk_df, bulk_correction_treshold=bulk_correction_treshold)
            sc_df['GBC'] = sc_df['GBC'].map(bulk_map)
        else:
            raise KeyError(
                f'{sample} is not present in sample_map.csv index. Check errors.'
            )
    
    ##

    # Get CBC-GBC combos and their supporting nUMIs
    counts = count_UMIs(sc_df)
    df_combos = get_combos(counts)

    # Filter CBC-GBC combos
    M, _ = filter_and_pivot(
        df_combos, 
        umi_treshold=umi_treshold, 
        p_treshold=p_treshold, 
        max_ratio_treshold=max_ratio_treshold,
        normalized_abundance_treshold=normalized_abundance_treshold
    )

    ## 

    # General checks
    f.write(f'# General checks \n')
    if bulk_correction:
        f.write(f'- Unique, "good" GBCs, post correction (with bulk reference whitelist): {df_combos["GBC"].unique().size}\n')
    else:
        f.write(f'- Unique, "good" GBCs: {df_combos["GBC"].unique().size}\n')
    f.write(f'- n unsupported CBC-GBC combos (N.B. only good GBCs): {(df_combos["status"]=="unsupported").sum()}\n')
    f.write(f'- n supported CBC-GBC combos (N.B. only good GBCs): {(df_combos["status"]=="supported").sum()}\n')
    f.write(f'- Observed MOI (from supported CBC-GBC combos only): {(M>0).sum(axis=1).median():.2f} median, (+-{(M>0).sum(axis=1).std():.2f})\n\n')

    if bulk_correction:
        pseudobulk_sc = M.sum(axis=0) / M.sum(axis=0).sum()
        common = list(set(pseudobulk_sc.index) & set(bulk_df.index))
        if common:
            pseudobulk_sc = pseudobulk_sc.loc[common]
            bulk = bulk_df.loc[common]['read_count'] / bulk_df.loc[common]['read_count'].sum()
            corr = np.corrcoef(pseudobulk_sc, bulk)[0,1]
            f.write(f'# Individual "good" GBC sequences checks\n')
            f.write(f'- n "good" GBC sequences, bulk ({bulk_df.shape[0]}) vs sc ({n_good_gbcs})\n')
            f.write(f'- Fraction of bulk GBCs found in sc: {bulk_df.index.isin(df_combos["GBC"].unique()).sum()/bulk_df.shape[0]:.2f}\n')
            f.write(f'- Common GBCs abundance (normalized nUMIs from pseudobulk scRNA-seq vs normalized read counts from bulk DNA-seq) correlation: {corr:.2f}\n\n')

    sets = get_clones(M)
    GBC_set = list(chain.from_iterable(sets['GBC_set'].map(lambda x: x.split(';')).to_list()))
    redundancy = 1-np.unique(GBC_set).size/len(GBC_set)
    occurrences = pd.Series(GBC_set).value_counts().sort_values(ascending=False)
    f.write(f'# GBC sets (i.e. "good" GBCs combinations which are supported in a given cell subset) checks \n')
    f.write(f'- Unique GBCs sets: {sets.shape[0]}\n')
    f.write(f'- Unique GBCs in these sets: {np.unique(GBC_set).size}\n')
    f.write(f'- GBCs redundancy across sets: {redundancy:.2f}\n')
    f.write(f'- GBCs occurrences across sets: {occurrences.median():.2f} (+- {occurrences.std():.2f})\n\n')
    
    # Get clones (from cell with 1 single GBC only here) and cells tables
    unique_cells = (M>0).sum(axis=1).loc[lambda x: x==1].index
    filtered_M = M.loc[unique_cells]
    clones_df = get_clones(filtered_M)
    cells_df = (
        filtered_M
        .apply(lambda x: filtered_M.columns[x>0][0], axis=1)
        .to_frame('GBC')
    )
    cells_df.index = cells_df.index.map(lambda x: f'{x}_{sample}')

    # Final clones checks
    f.write(f'# Final clones (i.e., distinct populations of uniquely barcoded cells only) checks \n')
    f.write(f'- n starting CBC (STARSolo): {df_combos["CBC"].unique().size}\n')
    f.write(f'- n uniquely barcoded cells: {cells_df.shape[0]}\n')
    f.write(f'- n clones: {clones_df.shape[0]}\n')
    f.write(f'- n clones>=10 cells: {clones_df["n cells"].loc[lambda x:x>=10].size}\n')

    # Save
    df_combos.to_csv('CBC_GBC_combos.tsv.gz', sep='\t')
    clones_df.to_csv('clones_summary_table.csv')
    cells_df.to_csv('cells_summary_table.csv')

    # --- Plots ---
    # CBC-GBC combo support plot (static)
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.kdeplot(data=df_combos, x='max_ratio', y='normalized_abundance', fill=True, cmap="viridis", ax=ax)
    ax.axvline(x=max_ratio_treshold, linestyle='--', color='red', label=f"Max Ratio Threshold ({max_ratio_treshold})")
    ax.axhline(y=normalized_abundance_treshold, linestyle='--', color='blue', label=f"Normalized Abundance Threshold ({normalized_abundance_treshold})")
    ax.text(0.05, 0.93, f"Total CBC-GBC: {df_combos.shape[0]:,}", transform=ax.transAxes, fontsize=13)
    n_filtered = (df_combos['status'] == "supported").sum()
    ax.text(0.05, 0.89, f"Retained: {n_filtered:,} ({n_filtered/df_combos.shape[0]*100:.1f}%)", transform=ax.transAxes, fontsize=13)
    ax.set_xlabel("Max Ratio", fontsize=17)
    ax.set_ylabel("Normalized Abundance", fontsize=17)
    ax.set_title("CBC-GBC Combo Support", fontsize=18)
    ax.tick_params(axis='both', labelsize=15)
    ax.legend(loc="lower right", fontsize=13)
    fig.tight_layout()
    fig.savefig('CBC_GBC_combo_status.png', dpi=300)

    # CBC-GBC Pair Rank Plot (static)
    umi_vals = df_combos.loc[df_combos['status'] == "supported", 'umi'].sort_values(ascending=False).values
    ranks = np.arange(0, len(umi_vals) + 1)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(ranks, np.concatenate([[0], umi_vals]), color="steelblue", lw=2)
    median_umi = np.median(umi_vals)
    ax.axhline(median_umi, color='red', linestyle='--', label=f"Median = {median_umi:.0f}")
    ax.set_xscale('log')
    ax.set_yscale('log')
    import matplotlib.ticker as mticker
    y_min = umi_vals.min() if len(umi_vals) > 0 else 1
    y_max = umi_vals.max() if len(umi_vals) > 0 else 10
    yticks = np.geomspace(y_min, y_max, num=6)
    ax.set_yticks(yticks)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: '{:,.0f}'.format(x)))
    x_min = 0
    x_max = ranks.max() if len(ranks) > 0 else 10
    xticks = np.geomspace(1, x_max, num=6)
    xticks = np.insert(xticks, 0, 0)
    if len(xticks) > 1 and xticks[0] == xticks[1]:
        xticks = np.delete(xticks, 1)
    ax.set_xticks(xticks)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: '{:,.0f}'.format(x)))
    ax.set_title("UMI Count Distribution per CBC-GBC Pair", fontsize=14)
    ax.set_xlabel("CBC-GBC Pair Rank (log scale)")
    ax.set_ylabel("UMI Count (log scale)")
    ax.legend()
    ax.tick_params(axis='both', labelsize=13)
    # Ensure axes (spines) are visible and prominent
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    fig.tight_layout()
    fig.savefig('umi_distribution.png', dpi=300)

    # --- Interactive CBC-GBC Pair Rank Plot ---
    yticks_plotly = np.geomspace(y_min, y_max, num=6)
    xticks_plotly = np.geomspace(1, x_max, num=6)
    xticks_plotly = np.insert(xticks_plotly, 0, 0)
    if len(xticks_plotly) > 1 and xticks_plotly[0] == xticks_plotly[1]:
        xticks_plotly = np.delete(xticks_plotly, 1)
    fig_rank = go.Figure([
        go.Scatter(
            x=ranks, y=np.concatenate([[0], umi_vals]),
            mode='lines',
            line=dict(color="steelblue", width=2),
            name="UMI per CBC-GBC"
        ),
        go.Scatter(
            x=[ranks[0], ranks[-1]], y=[median_umi, median_umi],
            mode='lines',
            line=dict(color="red", dash="dash"),
            name=f"Median = {median_umi:.0f}"
        )
    ])
    fig_rank.update_layout(
        title="CBC-GBC Pair Rank Plot (UMI per CBC-GBC)",
        xaxis=dict(
            title="CBC-GBC Pair Rank (log scale)", 
            type="log", 
            ticks="outside",
            tickvals=xticks_plotly,
            ticktext=[f"{int(x):,.0f}" for x in xticks_plotly],
            tickfont=dict(size=16, family="Arial, sans-serif", color="black"),
            showline=True,
        ),
        yaxis=dict(
            title="UMI Count (log scale)",
            type="log",
            ticks="outside",
            tickvals=yticks_plotly,
            ticktext=[f"{int(y):,.0f}" for y in yticks_plotly],
            tickfont=dict(size=16, family="Arial, sans-serif", color="black"),
            showline=True,
        ),
        font=dict(family="Arial, sans-serif", size=15, color="black"),
        margin=dict(l=60, r=30, t=60, b=60),
        legend=dict(
            orientation="v",
            yanchor="top",
            y=0.98,
            xanchor="right",
            x=0.98,
            bgcolor="rgba(255,255,255,0.9)",
            bordercolor="black",
            borderwidth=1,
            font=dict(size=13, color="black")
        )
    )
    fig_rank.write_html('umi_distribution_interactive.html', include_plotlyjs="cdn", config={"displayModeBar": False, "displaylogo": False, "responsive": True})

    # MOI distribution (static)
    fig, ax = plt.subplots(figsize=(6, 4))
    moi_counts = (M > 0).sum(axis=1).value_counts().sort_index()
    sns.barplot(x=moi_counts.index, y=moi_counts.values, color="teal", ax=ax)
    for i, v in enumerate(moi_counts.values):
        ax.text(i, v + max(moi_counts.values)*0.01, f"{v:,}", ha="center", fontsize=9)
    ax.set_xlabel("Number of GBCs per Cell (MOI)")
    ax.set_yscale('log')
    ax.set_ylabel("Cell Count")
    ax.set_title("MOI Distribution", fontsize=14)
    y_min = moi_counts.values.min() if len(moi_counts.values) > 0 else 1
    y_max = moi_counts.values.max() if len(moi_counts.values) > 0 else 10
    yticks = np.geomspace(y_min, y_max, num=6)
    ax.set_yticks(yticks)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: '{:,.0f}'.format(x)))
    ax.set_xticks(moi_counts.index)
    ax.set_xticklabels([str(x) for x in moi_counts.index])
    ax.tick_params(axis='both', labelsize=13)
    fig.tight_layout()
    fig.savefig('MOI_distribution.png', dpi=300)

    # --- Interactive MOI Distribution ---
    yticks_plotly = np.geomspace(y_min, y_max, num=6)
    xticks_plotly = [str(x) for x in moi_counts.index]
    fig_moi = go.Figure([
        go.Bar(
            x=moi_counts.index.astype(str),
            y=moi_counts.values,
            marker_color="teal",
            text=[f"{v:,}" for v in moi_counts.values]
        )
    ])
    fig_moi.update_layout(
        barmode="group",
        title=dict(text="MOI Distribution", font=dict(size=16, color="black")),
        xaxis=dict(
            title="Number of GBCs per Cell (MOI)",
            ticks="outside",
            tickangle=15,
            tickfont=dict(size=16, family="Arial, sans-serif", color="black"),
            showgrid=False,
            zeroline=False,
            linecolor='black',
            linewidth=1,
            mirror=True,
            tickvals=xticks_plotly,
            ticktext=xticks_plotly,
            showline=True,
        ),
        yaxis=dict(
            title="Cell Count",
            type="log",
            ticks="outside",
            tickfont=dict(size=16, family="Arial, sans-serif", color="black"),
            showgrid=False,
            zeroline=False,
            linecolor='black',
            linewidth=1,
            mirror=True,
            tickvals=yticks_plotly,
            ticktext=[f"{int(y):,.0f}" for y in yticks_plotly],
            showline=True,
        ),
        font=dict(family="Arial, sans-serif", size=14, color="black"),
        margin=dict(l=60, r=30, t=60, b=60),
        plot_bgcolor="white",
        paper_bgcolor="white"
    )
    fig_moi.write_html('MOI_distribution_interactive.html', include_plotlyjs="cdn", config={"displayModeBar": False, "displaylogo": False, "responsive": True})

    # --- Cumulative Clone Size Distribution (static) ---
    sorted_sizes = np.sort(clones_df['n cells'].values)
    cum_clones = np.arange(1, len(sorted_sizes) + 1)[::-1]  
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(sorted_sizes, cum_clones, drawstyle='steps-post', color="darkorange")
    median_clone = np.median(sorted_sizes)
    ax.axvline(median_clone, color='red', linestyle='--', label=f"Median = {int(median_clone)}")

    if len(sorted_sizes) < 2 or np.all(sorted_sizes == sorted_sizes[0]) or np.all(cum_clones == cum_clones[0]):
        ax.set_xscale('linear')
        ax.set_yscale('linear')
        ax.set_xlim(0.5, 1.5 if len(sorted_sizes)==1 else sorted_sizes.max()*1.1)
        ax.set_ylim(0.5, 1.5 if len(cum_clones)==1 else cum_clones.max()*1.1)
        ax.set_xticks([sorted_sizes[0]] if len(sorted_sizes)==1 else sorted_sizes)
        ax.set_yticks([cum_clones[0]] if len(cum_clones)==1 else cum_clones)
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: '{:,.0f}'.format(x)))
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: '{:,.0f}'.format(x)))
    else:
        ax.set_xscale('log')
        ax.set_yscale('log')
        x_min = sorted_sizes.min()
        x_max = sorted_sizes.max()
        y_min = cum_clones.min()
        y_max = cum_clones.max()
        xticks = np.geomspace(x_min, x_max, num=6)
        yticks = np.geomspace(y_min, y_max, num=6)
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: '{:,.0f}'.format(x)))
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: '{:,.0f}'.format(x)))
    ax.set_xlabel("Clone Size (Number of Cells, log scale)")
    ax.set_ylabel("Number of Clones (log scale)")
    ax.set_title("Cumulative Clone Size Distribution", fontsize=14)
    ax.legend()
    ax.grid(True, which="both", ls="--", lw=0.5)
    ax.tick_params(axis='both', labelsize=13)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    fig.tight_layout()
    fig.savefig('clone_size_distribution.png', dpi=300)

    # --- Interactive Cumulative Clone Size Plot ---
    fig_cum = go.Figure([
        go.Scatter(
            x=sorted_sizes,
            y=cum_clones,
            mode='lines',
            line=dict(color="darkorange"),
            name="Cumulative"
        ),
        go.Scatter(
            x=[median_clone, median_clone],
            y=[cum_clones.min(), cum_clones.max()],
            mode='lines',
            line=dict(color="red", dash="dash"),
            name=f"Median = {int(median_clone)}"
        )
    ])
    x_min = sorted_sizes.min() if len(sorted_sizes) > 0 else 1
    x_max = sorted_sizes.max() if len(sorted_sizes) > 0 else 10
    y_min = cum_clones.min() if len(cum_clones) > 0 else 1
    y_max = cum_clones.max() if len(cum_clones) > 0 else 10
    if x_min <= 0: x_min = 1
    if x_max <= 0: x_max = 1
    if y_min <= 0: y_min = 1
    if y_max <= 0: y_max = 1
    xticks_plotly = np.geomspace(x_min, x_max, num=6)
    yticks_plotly = np.geomspace(y_min, y_max, num=6)
    fig_cum.update_layout(
        title="Cumulative Clone Size Distribution",
        xaxis=dict(
            title="Clone Size (Number of Cells, log scale)",
            type="log",
            showline=True,
            linecolor='black',
            linewidth=2,
            mirror=True,
            ticks="outside",
            tickvals=xticks_plotly,
            ticktext=[f"{int(x):,.0f}" for x in xticks_plotly],
            tickfont=dict(size=16, family="Arial, sans-serif", color="black"),
        ),
        yaxis=dict(
            title="Number of Clones (log scale)",
            type="log",
            showline=True,
            linecolor='black',
            linewidth=2,
            mirror=True,
            ticks="outside",
            tickvals=yticks_plotly,
            ticktext=[f"{int(y):,.0f}" for y in yticks_plotly],
            tickfont=dict(size=16, family="Arial, sans-serif", color="black"),
        ),
        font=dict(family="Arial, sans-serif", size=14, color="black"),
        margin=dict(l=60, r=30, t=60, b=60),
        plot_bgcolor="white",
        paper_bgcolor="white",
        legend=dict(
            orientation="v",
            yanchor="top",
            y=0.98,
            xanchor="right",
            x=0.98,
            bgcolor="rgba(255,255,255,0.9)",
            bordercolor="black",
            borderwidth=1,
            font=dict(size=13, color="black")
        )
    )
    fig_cum.write_html('clone_size_distribution_interactive.html', include_plotlyjs="cdn", config={"displayModeBar": False, "displaylogo": False, "responsive": True})

    f.write(f'- Execution time: {T.stop()}\n')
    f.close()



##
