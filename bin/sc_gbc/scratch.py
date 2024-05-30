import pysam
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches



##


#Parameter & data
bam_file = '/Users/ieo6943/Downloads/cell_folder/consensus.bam'
output_bam = '/Users/ieo6943/Downloads/cell_folder/filter_local.bam'
base_er = 0.1
quality_th = 30
er_th = 0.1
mask_th = 0.2
mask_th_GBC = 0
mean_qual_th = 30

bam_in = pysam.AlignmentFile(bam_file, "rb", check_sq=False)
bam_out = pysam.AlignmentFile(output_bam, "wb", template=bam_in)
df = pd.DataFrame(columns=['read', 'ce', 'quality','UMI','n_read_consensus', 'SUPPORT'])
GBC_bases = True


##


#Script
j=0
for read in bam_in:
    
    n_consensus_read = read.get_tag('cM:')
    mean_consensus_error = read.get_tag('cE:')
    base_consensus_error = np.array(read.get_tag('ce:'))/n_consensus_read
    high_consensus_err_index = np.where(base_consensus_error>base_er)[0]
    umi = read.get_tag('MI:')
    
    sequence = list(read.query_sequence)
    qualities = np.array(read.query_qualities)
    mean_qualities = qualities.sum()/qualities.shape[0]

    #MASKs for the bad bases #mettere np.where per la maschera
    low_quality_indices = np.where(qualities<quality_th)[0]
    # Replace low quality bases with 'N'
    for i in low_quality_indices:
        sequence[i] = 'N'
    for i in high_consensus_err_index:
        sequence[i] = 'N'
    bad_reads_density=sequence.count('N')/len(sequence)
    sequence = np.array(sequence)
    np.where(sequence=='N',1,0).sum()

    GBC_n_N = np.where(sequence[33:33+18]=='N',1,0).sum()
    #GBC_n_N = 0

    supported = 'not supported'

    if mean_consensus_error<er_th and mean_qualities>=mean_qual_th and bad_reads_density<mask_th :
        if GBC_bases:
            if GBC_n_N==mask_th_GBC:
                supported = 'supported'

        supported = 'supported'
    
    if supported=='supported':
        read.query_sequence = ''.join(sequence)
        bam_out.write(read)

    final_seq = ''.join(sequence)
    #save metadata for filtering analysis    
    df.loc[j] = [final_seq, base_consensus_error, qualities,umi,n_consensus_read, supported]
    j+=1
#if j==10: break
bam_in.close()
bam_out.close()


##


#VISUALIZATION


#SUPPORTED
#Histogram number of read for each read consensus
df_sup = df[df['SUPPORT']=='supported'].reset_index()
df_sup.rename(columns={'read': 'GBC'}, inplace=True)

plt.figure(figsize=(10, 6))
bar_plot = sns.barplot(x='index', y='n_read_consensus', data=df_sup)
plt.title('read count for each SUPPORTED consensus read')
plt.xlabel('# consensus read')
plt.ylabel('Number of read')
plt.xticks([])
#plt.show()
plt.savefig('/Users/ieo6943/Downloads/cell_folder/molecular_consensus_size.png')

#NOTSUPPORTED
##Histogram number of read for each read consensus
df_nosup = df[df['SUPPORT']=='not supported'].reset_index()
df_nosup.rename(columns={'read': 'GBC'}, inplace=True)

plt.figure(figsize=(10, 6))
bar_plot = sns.barplot(x='index', y='n_read_consensus', data=df_nosup)
plt.title('read count for each NOT SUPPORTED consensus read')
plt.xlabel('# consensus read')
plt.ylabel('Number of read')
plt.xticks([])
#plt.show()
plt.savefig('/Users/ieo6943/Downloads/cell_folder/molecular_consensus_size_NOSUP.png')







#SUPPORTED
#supported GBC data
GBC = df[df['SUPPORT']=='supported']['read']
GBC = GBC.apply(lambda x: x[33: 33+18])
df_sup['GBC'] = GBC.values

#NOT SUPPORTED
#not supported GBC data
GBC_nosup = df[df['SUPPORT']=='not supported']['read']
GBC_nosup = GBC_nosup.apply(lambda x: x[33: 33+18])
df_nosup['GBC'] = GBC_nosup.values



# GROUP SIZE DENSITY


#supported
# Calculate unique UMI counts per GBC
df_GBC = df_sup.groupby(['GBC'])['UMI'].nunique().to_frame('umi').reset_index()
#df_GBC = df_sup.groupby(['GBC'])['quality'].sum().to_frame('quality').reset_index()

# Calculate the number of 'N's in each GBC
n_N = df_GBC['GBC'].apply(lambda x: (pd.Series(list(x))=='N').sum())
df_GBC['#N'] = n_N

# Create the plot
plt.figure(figsize=(10, 6))
palette = sns.color_palette("viridis", n_colors=df_GBC['#N'].max()+1)  # Generate a color palette

# Map '#N' values to colors
color_map = {n: palette[n] for n in df_GBC['#N'].unique()}
bar_colors = [color_map[n] for n in df_GBC['#N']]

bar_plot = sns.barplot(x='GBC', y='umi', data=df_GBC, palette=bar_colors)
plt.title('UMI count for each GBC MASKon')
plt.xlabel('# consensus GBC')
plt.ylabel('Number of umi')
plt.xticks([])  # Hide x-tick labels
# Create a legend
handles = [mpatches.Patch(color=color_map[n], label=f'#N = {n}') for n in sorted(df_GBC['#N'].unique())]
plt.legend(handles=handles, title="Number of N's", title_fontsize='13', loc='upper right', fontsize='12', bbox_to_anchor=(1.05, 1), borderaxespad=0.)
#plt.show()  # Uncomment to show the plot
plt.savefig('/Users/ieo6943/Downloads/cell_folder/group_size_with_mask.png')


#not supported
df_GBC_no = df_nosup.groupby(['GBC'])['UMI'].nunique().to_frame('umi').reset_index()
n_N = df_GBC_no['GBC'].apply(lambda x: (pd.Series(list(x))=='N').sum())
df_GBC_no['#N'] = n_N
# Create the plot
plt.figure(figsize=(10, 6))
palette = sns.color_palette("viridis", n_colors=df_GBC_no['#N'].max()+1)  # Generate a color palette

# Map '#N' values to colors
color_map = {n: palette[n] for n in df_GBC_no['#N'].unique()}
bar_colors = [color_map[n] for n in df_GBC_no['#N']]

bar_plot = sns.barplot(x='GBC', y='umi', data=df_GBC_no, palette=bar_colors)
plt.title('UMI count for each GBC NOT SUPPORTED MASKon')
plt.xlabel('# consensus GBC')
plt.ylabel('Number of umi')
plt.xticks([])  # Hide x-tick labels
# Create a legend
handles = [mpatches.Patch(color=color_map[n], label=f'#N = {n}') for n in sorted(df_GBC_no['#N'].unique())]
plt.legend(handles=handles, title="Number of N's", title_fontsize='13', loc='upper right', fontsize='12', bbox_to_anchor=(1.05, 1), borderaxespad=0.)
#plt.show()  # Uncomment to show the plot
plt.savefig('/Users/ieo6943/Downloads/cell_folder/group_size_NOSUP_with_mask.png')



# Number of N histogram



gbc_position_df = df_GBC['GBC'].apply(lambda x: pd.Series(list(x)))
n_counts = (gbc_position_df == 'N').sum().values
n_counts_df = pd.DataFrame({'Position': range(1, 19), 'Count of N': n_counts})
n_counts_df['Count of N'] = n_counts_df['Count of N']
plt.figure(figsize=(10, 6))
sns.barplot(x='Position', y='Count of N', data=n_counts_df)
plt.title('Frequency of N at Each Position in GBC Strings')
plt.xlabel('Position')
plt.ylabel('Count of N')
#plt.show()

plt.savefig('/Users/ieo6943/Downloads/cell_folder/histogram_N.png')



# QUALITY OF EACH BASE
df_quality = pd.DataFrame(df_sup['quality'].to_list())
df_quality.columns = [f'Value_{i+1}' for i in range(df_quality.shape[1])]
data = pd.DataFrame({
    'Position': range(df_quality.mean().shape[0]),
    'quality': df_quality.mean(),
    'Error': df_quality.std()
})

# Create the scatter plot
plt.figure(figsize=(10, 6))
plot = sns.scatterplot(x='Position', y='quality', data=data)

# Add error bars
for idx, row in data.iterrows():
    plt.errorbar(row['Position'], row['quality'], yerr=row['Error'], fmt='o', color='blue')

plt.title('Scatter Plot with Error Bars')
plt.xlabel('Position')
plt.ylabel('quality')
plt.grid(True)
#plt.show()
plt.savefig('/Users/ieo6943/Downloads/cell_folder/scatter_quality_all.png')



#CONSENSUS ERROR DISTANCE

#SUPPORTED
df_ce = pd.DataFrame(df_sup['ce'].to_list())
df_ce.columns = [f'Value_{i+1}' for i in range(df_ce.shape[1])]
data = pd.DataFrame({
    'Position': range(df_ce.mean().shape[0]),
    'ce': df_ce.mean(),
    'Error': df_ce.std()
})

# Create the scatter plot
plt.figure(figsize=(10, 6))
plot = sns.scatterplot(x='Position', y='ce', data=data)

# Add error bars
for idx, row in data.iterrows():
    plt.errorbar(row['Position'], row['ce'], yerr=row['Error'], fmt='o', color='blue')

plt.title('Scatter Plot with Error Bars')
plt.xlabel('Position')
plt.ylabel('ce')
plt.grid(True)
#plt.show()
plt.savefig('/Users/ieo6943/Downloads/cell_folder/scatter_ce_all.png')



