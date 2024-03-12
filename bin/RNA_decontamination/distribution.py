
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plotting_utils._plotting_base import *
import umap
from scipy.stats import pearsonr
from mito_utils.kNN import *
from mito_utils.clustering import *

sys.path.append("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/sc_gbc")
from helpers import *

sample_params=None
correction_type='reference' 
sc_correction_treshold=3
ncores=8
filtering_method="fixed"
coverage_treshold=70 
umi_treshold=5
p_treshold=1
max_ratio_treshold=.8
normalized_abundance_treshold=.8
path = '/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/results/'
path_empty = '/Users/ieo6943/Documents/Guido/empty_bead/'
path_count='/Users/ieo6943/Documents/Guido/counts.pickle'

def majority_rule(row):
    max_val = row.max()
    return row.apply(lambda x: x if x == max_val else 0)

def Knn_rule(m_df, z):   
    for i in m_df.index:
        m_df.iloc[i][lambda x: x != z[i]] = 0

def set_to_one(x):
    return 1 if x != 0 else 0

def moi_counts(m_df):

    'It returns the counts sum of each cell'

    moi = list(m_df.sum(axis=1))
    return moi

def distance(m_df, m_df_dec):

    'It returns the percentage of the variation in the UMi counts and GBC counts after decontamination'

    m_discrete=m_df.copy()
    m_discrete[m_discrete>0]=1
    m_discrete_dec=m_df_dec.copy()
    m_discrete_dec[m_discrete_dec>0]=1
    m_discrete_dec[m_df_dec<0]=0
    diff_dis=m_discrete-m_discrete_dec
    dist_dis=[]
    diff=m_df-m_df_dec
    dist=[]
    assert (diff >= 0).all().all(), "The decontamination should only remove elements"
    assert (diff_dis >= 0).all().all(), "The decontamination should only remove elements"
    for i in range(m_df.shape[0]):
            dist_dis.append(sum(diff_dis.iloc[i]))
            dist.append(sum(diff.iloc[i]))
    UMI_variation = sum(dist)/(sum(sum(m_df.values[:])))
    GBC_variation = sum(dist_dis)/(sum(sum(m_discrete.values[:])))
    data=[GBC_variation, UMI_variation]
    return data

def complete_corr_relevance(m_df, m_df_dec):
    total_correction=m_df_dec[(m_df!=m_df_dec)].count().sum()
    X = m_df.apply(lambda x: x/x.sum(), axis=1)
    level_of_correction=X[(m_df!=m_df_dec) & (m_df_dec==0)].values.flatten()
    nan_mask = np.isnan(level_of_correction)
    level_of_correction = level_of_correction[~nan_mask]
    return total_correction, level_of_correction

def moi_frequency(moi):
    df = pd.DataFrame(moi, columns=['Valore'])

    # Calcolo della tabella di frequenza
    tabella_frequenza = df['Valore'].value_counts().reset_index()
    tabella_frequenza.columns = ['Valore', 'Frequenza']

    # Ordina la tabella di frequenza per valore
    tabella_frequenza = tabella_frequenza.sort_values('Valore')
    tabella_frequenza.to_csv('/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/results/MOI_'+coverage_treshold+'.csv', index=False) 
    return

def hist(list_, coverage_treshold, titolo):  
    
    plt.hist(list_, bins=10, color='skyblue', edgecolor='black')
    plt.savefig(f'/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/results/Istogramma dei valori{coverage_treshold}.jpg')

def double_hist(moi_m, moi_dec,coverage_treshold, title):
    fig, axs = plt.subplots(1,2)
    axs[0].hist(moi_m, bins=10, color='skyblue', edgecolor='black')
    # Aggiunta di titoli e label agli assi
    axs[0].set_title( title +f'{coverage_treshold}')
    axs[0].set_xlabel('Valore')
    axs[0].set_ylabel('Frequenza')

    axs[1].hist(moi_dec, bins=10, color='skyblue', edgecolor='black')
    # Aggiunta di titoli e label agli assi
    axs[1].set_title( title +f'{coverage_treshold}_dec')
    axs[1].set_xlabel('Valore')
    axs[1].set_ylabel('Frequenza')
    plt.savefig('/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/results/'+title+f'{coverage_treshold}.jpg')
    return
    
def calculateUmap(m_df, n_neighbors=15):
    seed=1234
    X = m_df.apply(lambda x: x/x.sum(), axis=1)
    X= np.log2(X+1)
    reducer = umap.UMAP(min_dist=0.01, spread=1, n_neighbors=n_neighbors)  # Initialize UMAP reducer
    embedding = reducer.fit_transform(X)  # Perform UMAP dimensionality reduction    
    return embedding

def calculateClonal_labels(m_df):
    clonal_labels=m_df.idxmax(axis=1).values
    return clonal_labels

def plot_double_UMAP(m_df,m_df_dec, coverage_treshold):
    
    # Plot the UMAP embeddings with cluster labels
    m_embedding=calculateUmap(m_df)
    m_clonal_labels=calculateClonal_labels(m_df)
    m_n_clones=np.unique(m_clonal_labels)

    m_dec_embedding=calculateUmap(m_df_dec)
    m_dec_clonal_labels=calculateClonal_labels(m_df_dec)
    m_dec_n_clones=np.unique(m_dec_clonal_labels)

    fig = plt.figure(figsize=(8,4))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    for i in m_n_clones:
       ax1.scatter(m_embedding[m_clonal_labels == i, 0], m_embedding[m_clonal_labels == i, 1], label=f'Cluster {i+1}',s=3)
    for i in m_dec_n_clones:
        ax2.scatter(m_dec_embedding[m_dec_clonal_labels == i, 0], m_dec_embedding[m_dec_clonal_labels  == i, 1], label=f'Cluster {i+1}',s=3)
    ax1.set_title(f'UMAP lentiviral colors thr={coverage_treshold}')
    ax2.set_title(f'UMAP dec. lentiviral colors thr={coverage_treshold}')

    fig.savefig(f'/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/results/UMAP_{coverage_treshold}.png', dpi=300) 

def one_hist(d,n_bins,coverage_treshold,title):
        plt.figure()
        plt.hist(d, bins=n_bins, color='skyblue', edgecolor='black')
        plt.xlabel('Valore')
        plt.ylabel('Frequenza')
        plt.title(title+str(coverage_treshold))
        #plt.show()
        plt.savefig('/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/results/'+title+f'_{coverage_treshold}.png')

def max_feature_percentage(m_df):
    feature_percentage=[]
    for i in range(len(m_df)):
        feature_percentage.append(max(m_df.iloc[i])/sum(m_df.iloc[i])*100)
    return feature_percentage

def scatter_correlation(other_feature, contamination,path):
    plt.figure()
    plt.scatter(other_feature, contamination, s=3)
    plt.gca().set_aspect('equal', adjustable='box')
    m, q = np.polyfit(other_feature, contamination, 1)
    plt.plot(np.array(other_feature), m*np.array(other_feature) + q, color='red')
    plt.xlim(0, 1)  
    plt.ylim(0, 1)
    plt.xlabel('Other barcodes ratio')
    plt.ylabel('DecontX contamination ratio')
    plt.title(f'Correlation with th= {coverage_treshold}')
    correlation_coef, p_value = pearsonr(other_feature, contamination)
    # Annotate plot with correlation coefficient and p-value
    plt.annotate(f'Correlation coefficient: {correlation_coef:.2f}\nP-value: {p_value:.2f}',
                xy=(np.mean(correlation_coef), np.mean(contamination)),  # Position of the annotation
                xytext=(0.4, 0.05),             # Offset of the text from the annotation position
                textcoords='axes fraction',  # Coordinate system of xytext
                fontsize=8,     )         # Font size of the annotation text
                #arrowprops=dict(facecolor='black', arrowstyle='->'))  # Arrow properties
        
    plt.savefig(path+f'correlation_contamination_{coverage_treshold}.png', dpi=300)


def filtering_cell(m):
    
    df_combos = m.reset_index().melt(id_vars='CBC', value_vars=m.columns, value_name='umi')



    df_combos = df_combos.assign(
                    max_ratio=lambda x: \
                    x.groupby('CBC')['umi'].transform(lambda x: x/x.max()),
                    normalized_abundance=lambda x: \
                    x.groupby('CBC')['umi'].transform(lambda x: x/x.sum()))


    (df_combos['normalized_abundance']>0.8).sum()

    M, _ = filter_and_pivot(
            df_combos, 
            umi_treshold=umi_treshold, 
            p_treshold=p_treshold, 
            max_ratio_treshold=max_ratio_treshold,
            normalized_abundance_treshold=normalized_abundance_treshold
        )

    # Get clones (from cell with 1 single GBC only here) and cells tables

    # Get 1-GBC CBCs
    unique_cells = (M>0).sum(axis=1).loc[lambda x: x==1].index
    filtered_M = M.loc[unique_cells]
    cells_df = (
        filtered_M
        .apply(lambda x: filtered_M.columns[x>0][0], axis=1)
        .to_frame('GBC_set')
    )
    return filtered_M, cells_df 

def get_filtered_labels(m):
        m.columns=list(range(m.shape[1]))
        cells_df = (
            m
            .apply(lambda x: m.columns[x>0][0], axis=1)
            .to_frame('GBC_set')
        )
        return cells_df

def label_similarity(M, M_dec):
    M_dis = M.applymap(set_to_one)  
    M_dec_dis = M_dec.applymap(set_to_one)
    min_l = min(M_dec_dis.shape[0],M_dis.shape[0])
    merged_pivot= M_dis.add(M_dec_dis, fill_value=0)
    check = (merged_pivot>1).sum().sum()/min_l
    return check

def merge_common_rows(df1, df2):
    """
    Returns two DataFrames containing only the rows that are common.

    Arguments:
    df1, df2: Input DataFrames.

    Returns:
    Tuple containing the two filtered DataFrames.
    """
    # Find common indices
    common_index = df1.index.intersection(df2.index)
    
    # Select only rows with common indices
    df1_filtered = df1.loc[common_index]
    df2_filtered = df2.loc[common_index]
    
    return df1_filtered, df2_filtered


def create_boxplots(data_list, path, title, list_of_value, y_ax ):
    """
    Create a grid of 2x3 boxplots from a list of data arrays.

    Parameters:
        data_list (list of arrays): List containing the data arrays for boxplots.

    Returns:
        None
    """
  # Check if the length of data_list is 6
    # Create a single plot
    plt.figure(figsize=(10, 6))

    # Create boxplot for each dataset
    positions = [0, 20, 40, 60, 80, 100]

    # Creazione del boxplot
    sns.boxplot(data=data_list, width=0.5)
    plt.xticks(range(len(positions)), positions)
    #plt.boxplot(data_list)

    # Set labels for x-axis
    #plt.xticks(range(1, 7), [f'{i}' for i in range(0,101,20)])

    # Set title and labels
    plt.title(title)
    plt.xlabel('Coverage threshold')
    plt.ylabel(y_ax)

    # Show the plot
    plt.grid(True)

    # Show the plot
    plt.savefig(path + f'{title}.png')

    
##


super_data = {'corrections_perc':[], 'complete_correction': []}

for coverage_treshold in list(range(0, 101, 20)):
    
    print(coverage_treshold)
    #path='/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/results/'
    path_count_pre='/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/results/M_'+str(coverage_treshold)+'.csv'
    path_count_dec='/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/results/M_'+str(coverage_treshold)+'_dec.csv'
    path_cluster_labels_dec='/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/results/clustering'+str(coverage_treshold)+'_dec.csv'
    path_umap_dec='/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/results/umap_'+str(coverage_treshold)+'_dec.csv'
    path_contamination=path+'contamination_'+str(coverage_treshold)+'.csv'
    m_df = pd.read_csv(path_count_pre, sep=',', index_col=0, dtype='str')
    m_df=m_df.apply(pd.to_numeric)
    m_df_old_columns=m_df.columns
    m_df.columns=list(range(m_df.shape[1]))
    m_df.index=list(range(len(m_df)))
    
    #umap_dec=pd.read_csv(path_umap_dec, sep=',', index_col=0).values
    #clustering_label_dec=pd.read_csv(path_cluster_labels_dec, sep=',', index_col=0).values[:,0]
    #umap_dec.index=clustering_label_dec.index
    m_df=m_df.apply(pd.to_numeric)
    
    m_df_dec = pd.read_csv(path_count_dec, sep=',', dtype='str')
    #m_df_dec.index = m_df.index
    m_df_dec.columns = list(range(m_df_dec.shape[1]))
    m_df_dec= m_df_dec.apply(pd.to_numeric)

    #combo_df_dec = m_df_dec.stack().reset_index(name='umi').rename(columns={'level_0': 'CBC', 'level_1': 'GBC'})
    #combo_df_dec = combo_df_dec[combo_df_dec['umi']!=0]


    #moi_df=moi_counts(m_df)
    #moi_df_dec=moi_counts(m_df_dec)
    #double_hist(moi_df, moi_df_dec, coverage_treshold)
    #plot_double_UMAP(m_df,m_df_dec,umap_dec, coverage_treshold)
    d=distance(m_df, m_df_dec)
    super_data['corrections_perc'].append(d)
    #one_hist(d, len(m_df.keys())+1,coverage_treshold,'number of correction for cell with coverage at')
    feature_percentage=max_feature_percentage(m_df)
    other_feature=list(1-np.array(max_feature_percentage(m_df))/100)
    contamination=list(pd.read_csv(path_contamination, sep=',').values.flatten())
    scatter_correlation(other_feature, contamination,path)
    total_n_correction, complete_correction_relevance=complete_corr_relevance(m_df, m_df_dec)
    percentage_complete_correction=complete_correction_relevance.shape[0]/total_n_correction
    super_data['complete_correction'].append([percentage_complete_correction, complete_correction_relevance])
    #one_hist(complete_correction_relevance,n_bins=10,coverage_treshold=coverage_treshold, title='complete_correction_relevance')

super_com_corr_rel=super_data['complete_correction']
data=[super_com_corr_rel[0][1],super_com_corr_rel[1][1],super_com_corr_rel[2][1],super_com_corr_rel[3][1],super_com_corr_rel[4][1],super_com_corr_rel[5][1]]
#create_boxplots(data)

with open(path_count, 'rb') as p:
    data= pickle.load(p)
median_data={}
for coverage_treshold in list(range(0, 101, 20)):
    print(coverage_treshold)
    # Mark noisy UMIs, from chosen correction method
    counts = mark_UMIs(data[correction_type], coverage_treshold=coverage_treshold)
    # Get CBC-GBC combos with selected UMIs
    df_combos = get_combos(counts, gbc_col=f'GBC_{correction_type}')
    median_data[coverage_treshold]=df_combos.groupby('CBC')['normalized_abundance'].median()
data=[median_data[0].values,median_data[20].values,median_data[40].values,median_data[60].values,median_data[80].values,median_data[100].values]
create_boxplots(data, path)








    
#super_sim_arr=np.array(list(super_similarity.values()))
corrections_perc=np.array(super_data['corrections_perc'])
#data=super_sim_arr
data= corrections_perc
plt.figure()
# Dati per le coppie di barre
labels = ['0', '20', '40', '60', '80', '100']  # Etichette delle barre
bar_width = 0.2  # Larghezza delle barre
index = np.arange(len(labels))  # Indici per posizionare le barre

# Dati per le prime 6 coppie di barre
data1 = data[:,0]
data2 = data[:,1]
#data3 = data[:,2]
#data4 = data[:,3]
#data5 = data[:,4]

# Creazione dell'istogramma
plt.bar(index, data1, bar_width, label='GBC')
plt.bar(index + bar_width, data2, bar_width, label='UMI')
#plt.bar(index + 2*bar_width, data3, bar_width, label='Cossa vs decontX and Knn')
# Aggiunta di etichette, titolo e legenda
plt.xlabel('coverage threshold')
plt.ylabel('loss ratio after correction')
#plt.title(' Similarity labels majority vs knn')
plt.xticks(index + bar_width, labels)
plt.legend()

# Mostra l'istogramma
plt.savefig(path +'Correction ratios.png')




############### Density plot decontamination

for coverage_treshold  in list(range(0, 101, 20)):
    path='/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/results/'
    path_count_pre='/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/results/M_'+str(coverage_treshold)+'.csv'
    path_count_dec='/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/results/M_'+str(coverage_treshold)+'_dec.csv'

    m_df = pd.read_csv(path_count_pre, sep=',', index_col=0, dtype='str')
    m_df=m_df.apply(pd.to_numeric)
    #m_df = m_df.apply(lambda x: x/x.sum(), axis=1)
    #m_df= np.log2(m_df+1)

    df_combos = m_df.stack().reset_index(name='umi').rename(columns={'level_0': 'CBC', 'level_1': 'GBC'})
    df_combos = df_combos[df_combos['umi']!=0]
    df_combos = df_combos.assign(
            max_ratio=lambda x: \
            x.groupby('CBC')['umi'].transform(lambda x: x/x.max()),
            normalized_abundance=lambda x: \
            x.groupby('CBC')['umi'].transform(lambda x: x/x.sum())
        )


    m_df_dec = pd.read_csv(path_count_dec, sep=',', dtype='str')
    m_df_dec.index = m_df.index
    m_df_dec= m_df_dec.apply(pd.to_numeric)

    #m_df_dec = m_df_dec.apply(lambda x: x/x.sum(), axis=1)
    #m_df_dec = np.log2(m_df_dec+1)

    df_combos_dec = m_df_dec.stack().reset_index(name='umi').rename(columns={'level_0': 'CBC', 'level_1': 'GBC'})
    df_combos_dec = df_combos_dec[df_combos_dec['umi']!=0]

    df_combos_dec = df_combos_dec.assign(
            max_ratio=lambda x: \
            x.groupby('CBC')['umi'].transform(lambda x: x/x.max()),
            normalized_abundance=lambda x: \
            x.groupby('CBC')['umi'].transform(lambda x: x/x.sum())
        )
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))

    sns.kdeplot(data=df_combos_dec, x = df_combos_dec['max_ratio'], y =df_combos_dec['normalized_abundance'], ax=axs[1])
    axs[1].set_title('After decontamination')
    sns.kdeplot(data=df_combos, x = df_combos['max_ratio'], y =df_combos['normalized_abundance'], ax=axs[0])
    axs[0].set_title('Before decontamination')

    fig.suptitle(f'density_plot_{coverage_treshold}', fontsize=16)

    fig.savefig(path+ f'density_plot_{coverage_treshold}')

############### Cluster distribution
path='/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/results/'
for coverage_treshold in list(range(0, 101, 20)):

    path_count_pre= path + '/M_'+str(coverage_treshold)+'.csv'
    path_count_dec= path + 'M_'+str(coverage_treshold)+'_dec.csv'

    m_df = pd.read_csv(path_count_pre, sep=',', index_col=0, dtype='str')
    m_df=m_df.apply(pd.to_numeric)
    #m_df = m_df.apply(lambda x: x/x.sum(), axis=1)
    #m_df= np.log2(m_df+1)

    df_combos = m_df.stack().reset_index(name='umi').rename(columns={'level_0': 'CBC', 'level_1': 'GBC'})
    df_combos = df_combos[df_combos['umi']!=0]
    df_combos = df_combos.assign(
            max_ratio=lambda x: \
            x.groupby('CBC')['umi'].transform(lambda x: x/x.max()),
            normalized_abundance=lambda x: \
            x.groupby('CBC')['umi'].transform(lambda x: x/x.sum())
        )


    m_df_dec = pd.read_csv(path_count_dec, sep=',', dtype='str')
    m_df_dec.index = m_df.index
    m_df_dec= m_df_dec.apply(pd.to_numeric)
    print('correzzioni cellule azzerate=',  (m_df_dec.sum(axis=1)==0).sum())
    #m_df_dec = m_df_dec.apply(lambda x: x/x.sum(), axis=1)
    #m_df_dec = np.log2(m_df_dec+1)

    df_combos_dec = m_df_dec.stack().reset_index(name='umi').rename(columns={'level_0': 'CBC', 'level_1': 'GBC'})
    df_combos_dec = df_combos_dec[df_combos_dec['umi']!=0]

    df_combos_dec = df_combos_dec.assign(
            max_ratio=lambda x: \
            x.groupby('CBC')['umi'].transform(lambda x: x/x.max()),
            normalized_abundance=lambda x: \
            x.groupby('CBC')['umi'].transform(lambda x: x/x.sum())
        )
 #  Take the information on the cluster from decontX
    z = pd.read_csv(path + f'z{coverage_treshold}.csv' , sep=',')
    df = pd.DataFrame({'CBC': list(m_df.index), 'cluster': list(z.values.flatten())})
    z_cluster = np.unique(z)
    df_combos_th= pd.merge(df_combos, df, on = 'CBC', how = 'inner')
    my_dict = dict(zip(m_df.columns, range(1,m_df.shape[1]+1)))
    df_combos_th['GBC'] = df_combos_th['GBC'].replace(my_dict)
    cluster={}
    cluster_data= {}
    palette = sns.color_palette("husl" ,len(z_cluster))
    dataset=[]
    plt.figure()
    for c in z_cluster:
        cluster[c]=df_combos_th[['normalized_abundance','GBC','CBC']][df_combos_th['cluster']== c]
        cluster[c]['normalized_abundance'] = round(cluster[c]['normalized_abundance']/len(np.unique(cluster[c]['CBC']))*1000)
        cluster_data[c] = cluster[c].groupby('GBC')['normalized_abundance'].sum().reset_index()
        data=[]
        for gbc in list(cluster_data[c]['GBC']):
            data += [gbc] * int(cluster_data[c][cluster_data[c]['GBC']==gbc]['normalized_abundance'].values)
        dataset.append(data)
        #sns.violinplot(data=data)
        #sns.histplot(data=cluster_data[c], x='normalized_abundance', hue='GBC', color='blue', label=f'Data{c}', alpha=0.5,multiple='stack', bins=5)
        #sns.histplot(data=cluster_data[c], color='blue', label=f'Data{c}', alpha=0.5, multiple='stack', bins=5)
        #sns.barplot(x=cluster_data[c]['GBC'], y=cluster_data[c]['normalized_abundance'],color=palette[c],label=f'cluster{c}', alpha=1-((c-1)/len(z_cluster)))
        #sns.lineplot(x=cluster_data[c]['GBC'], y=cluster_data[c]['normalized_abundance'],color=palette[c],label=f'cluster{c}', marker='o', linewidth=3)
    
    plt.figure(figsize=(6,4))
    positions = z_cluster
    sns.violinplot(data=dataset, positions=positions)
    plt.xlabel('clusters')
    plt.ylabel('Normalized abundance of GBC')
    plt.title(f'Coverage threshold ={coverage_treshold}')
    plt.savefig(path + f'istogrammi_sovrapposti{coverage_treshold}.png', dpi=300)
    

    # Create a violin plot for each dataset
    for i, data in enumerate(datasets, start=1):
        sns.violinplot(data=data, ax=plt.gca(), label=f'Dataset {i}')
    plt.savefig(path + 'istogrammi_sovrapposti.png')
        # Esempio di dati

    plt.figure()
    # Creazione del violin plot
    #ns.violinplot(data=cluster[2], x='GBC', y='normalized_abundance')
    sns.regplot(x=cluster[2]['GBC'], y=cluster[2]['normalized_abundance'], scatter=True)
    sns.lineplot(data=df, x='x', y='y', marker='o', linewidth=3)
    #sns.violinplot(x ='GBC', y= 'normalized_abundance', data=cluster[c])

    # Aggiunta di titolo e label agli assi
    plt.title('Violin Plot')
    plt.xlabel('Day')
    plt.ylabel('Total Bill')
    #sns.histplot(data=cluster_data[1],x='normalized_abundance', hue='GBC', color='blue', label=f'Data{1}', alpha=0.5,multiple='stack', bins=5)
    plt.savefig(path + 'istogrammi_sovrapposti.png')

plt.figure()

# Example datasets
data1 = [1, 2, 3, 4, 5]
data2 = [2, 3, 4, 5, 6]
data3 = [3, 4, 5, 6, 7]

# Set the positions on the x-axis for the violin plots
positions = [1, 2, 3]

# Create violin plots with specified positions
sns.violinplot(data=[data1, data2, data3], positions=positions)

# Show the plot
plt.show()

# Mostrare il plot
plt.savefig(path + 'istogrammi_sovrapposti.png')
  

### LEVEL of contamination
C_G = []
C_G_dec = []
G_G = []
C_C = []
path_count = '/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/data/counts.pickle'
with open(path_count, 'rb') as p:
    data= pickle.load(p)
for coverage_treshold in list(range(0, 101, 20)):

    path_count_pre= path + '/M_'+str(coverage_treshold)+'.csv'
    path_count_dec= path + 'M_'+str(coverage_treshold)+'_dec.csv'

    m_df = pd.read_csv(path_count_pre, sep=',', index_col=0, dtype='str')
    m_df=m_df.apply(pd.to_numeric)





    # Mark noisy UMIs, from chosen correction method
    counts = mark_UMIs(data[correction_type], coverage_treshold=coverage_treshold)

    # Get CBC-GBC combos with selected UMIs
    df_combos = get_combos(counts, gbc_col=f'GBC_{correction_type}')
    M, _ = filter_and_pivot(
    df_combos, 
    umi_treshold=umi_treshold, 
    p_treshold=p_treshold, 
    max_ratio_treshold=max_ratio_treshold,
    normalized_abundance_treshold=normalized_abundance_treshold
    )   

    G = m_df.max(axis=1)/m_df.sum(axis=1)
    C = 1 - m_df.max(axis=1)/m_df.sum(axis=1)
    C_G.append(list(C/G))
    m_df_dec = pd.read_csv(path_count_dec, sep=',', dtype='str')
    m_df_dec.index = m_df.index
    m_df_dec= m_df_dec.apply(pd.to_numeric)
    G_dec = m_df_dec.max(axis=1)/m_df_dec.sum(axis=1)
    C_dec = 1- m_df_dec.max(axis=1)/m_df_dec.sum(axis=1)
    C_G_dec.append(list(C_dec/G_dec))


    C_C.append(list(C_dec/C))
    G_G.append(G_dec/G)


create_boxplots(data_list=C_G_dec, path=path, title='Contamination vs major GBC after decontamination',list_of_value=list(range(0,101,20)), y_ax='C/G' )
create_boxplots(data_list=C_G, path=path, title='Contamination vs major GBC before decontamination',list_of_value=list(range(0,101,20)), y_ax='C/G' )
create_boxplots(data_list=G_G, path=path, title='Major GBC ratio after and before decontamination',list_of_value=list(range(0,101,20)), y_ax='G_dec/G' )
create_boxplots(data_list=C_C, path=path, title='Contamination ratio after and before decontamination',list_of_value=list(range(0,101,20)), y_ax='C_dec/C' )
 #  Take the information on the cluster from decontX
