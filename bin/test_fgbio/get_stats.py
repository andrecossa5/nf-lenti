#!/usr/bin/sh

import os
import pandas as pd


##


def main():

    d = {}
    d['median_filtered_base_umi_group_size'] = pd.read_csv('tables/median_filtered_base_umi_group_size.txt.gz', header=None, index_col=0).iloc[:,0].median()
    d['n_umis_filtered'] = pd.read_csv('tables/n_umis_filtered.txt.gz', header=None, index_col=0).iloc[:,0].median()
    d['median_filtered_base_consensus_error'] = pd.read_csv('tables/median_filtered_base_consensus_error.txt.gz', header=None, index_col=0).iloc[:,0].median()
    d['coverage'] = pd.read_csv('tables/coverage.txt.gz', header=None).iloc[:,2].median()
    d['n_umis_unfiltered'] = pd.read_csv('tables/n_umis_unfiltered.txt.gz', header=None, index_col=0).iloc[:,0].median()
    d['depth'] = pd.read_csv('tables/depth.txt.gz', header=None, index_col=0).iloc[:,0].median()
    d['median_filtered_read_quality'] = pd.read_csv('tables/median_filtered_read_quality.txt.gz', header=None, index_col=0).iloc[:,0].median()
    d['n_reads_unfiltered'] = pd.read_csv('tables/n_reads_unfiltered.txt.gz', header=None, index_col=0).iloc[:,0].median()
    d['n_reads_filtered'] = pd.read_csv('tables/n_reads_filtered.txt.gz', header=None, index_col=0).iloc[:,0].median()

    pd.Series(d).to_frame('value').reset_index().rename(columns={'index':'metric'}).to_csv('stats.csv')


##


if __name__ == '__main__':
    main()