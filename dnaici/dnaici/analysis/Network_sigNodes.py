import glob
import os
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, norm
import matplotlib as mlp
import matplotlib.pyplot as plt
import seaborn as sns
from dnaici.analysis.Enrichment_permutation import add_TvalSign2log10Pval


cluster_strs={0:'_cn1_',1:'_cn2_',2:'_cn3_',3:'_cn4_',4:'_cn5_',5:'_cn6_',6:'_cn7_',7:'_cn8_', 8:'_cn9_',9:'_cn10_'}

def find_cluster(x,cluster_strs):
    ''' find a cluster string (e.g., _cn1_, _cn2_, ..) for an input file name
    '''
    if x in cluster_strs.keys():
        x_out= cluster_strs[x]
    else:
        x_out='not_find'
    
    return x_out


def countNum(selected_df):
    
    chrom = ['chr'+str(i) for i in range(1,24)]
    count = []
    
    for i in range(1,24):
        c = selected_df.id.str.count('chr'+str(i)+'_').sum()
        count.append(c)
    
    num_df = pd.DataFrame({'chrom':chrom, 'count':count})
    
    return num_df


def find_tval2centrality(in_features_pval_file, in_nodes_pair_cluster_files,in_centrality_file,tval_str='tval'):
    ''' assign T-values/P-values of feautres from network clustering to nodes based on their cluster label,
    then convert expected pvalue of node centrality to abs(log10(pval)) for representing sig_centrality.
    Output, a dataframe by adding id and sig_centrality to loaded feature_pval_df and  centrality_df
    tval_str is used to select either tval or pval for the calculation
    '''
    #set a minimum pval for replacing 0 which is the same as min_p in plot_heatmap4tval
    minimum_pval=10**(-6)
    print('Read feature pval file:', '\n')
    print(in_features_pval_file, '\n')
    if os.path.exists(in_features_pval_file):
        #read feature pval
        feature_pval_df=pd.read_csv(in_features_pval_file,sep='\t',index_col=0)
        #read nodes pair
        print('Read nodes pair network cluster file:', '\n')
        print(in_nodes_pair_cluster_files, '\n')
        nodes_pair_clusters={}
        for fi in in_nodes_pair_cluster_files:
            nodes_pair_clusters[fi]=pd.read_csv(fi,sep='\t')
        #read nodes centrality
        #in_centrality_path=os.path.join(in_path0,  time_str , 'out_network_cluster')
        #in_centrality_file=glob.glob( os.path.join(in_centrality_path, chrom_str+ '_*_communities_info.tsv'))
        print('Read node centrality file:', '\n')
        print(in_centrality_file, '\n')
        centrality_df=pd.read_csv(in_centrality_file[0],sep='\t')
        
        #for each node find its features T-values/P-values
        centrality_df['cluster_str']=centrality_df.labels.apply(lambda x : find_cluster(x, cluster_strs))
        record_id_df=[]
        for idx, cs in centrality_df.iterrows():
            cs_id=cs.id
            cs_cluster=cs.cluster_str
            #find tval or pval for the corresponding network cluster
            out_df = add_TvalSign2log10Pval(feature_pval_df,tval_str,cluster_str=cs_cluster,is_transpose=False,min_p=minimum_pval) 
            if not out_df.empty:
                out_df.columns=[cs_id]
                if len(record_id_df)==0:
                    record_id_df=out_df.copy()
                else:
                    record_id_df=pd.concat([record_id_df,out_df],axis=1).copy()
        record_id_df=record_id_df.T
        record_id_df=(record_id_df-record_id_df.min())/(record_id_df.max()-record_id_df.min())
        #record_id_df = record_id_df.apply(zscore)
        record_id_df['id']=record_id_df.index
        out_df=pd.merge(centrality_df,record_id_df, on='id',how='left').copy()
        #replace zero pvalue by a minimum pval
        #out_df.loc[out_df.pval_1000sample10times<minimum_pval,'pval_1000sample10times']=minimum_pval
        out_df_scaled = out_df.copy()
        column = 'centrality'
        out_df_scaled['normalized_centrality'] = (out_df_scaled[column] - out_df_scaled[column].min()) / (out_df_scaled[column].max() - out_df_scaled[column].min())
        #out_df_scaled['normalized_centrality'] = (out_df_scaled[column] - out_df_scaled[column].mean()) / out_df_scaled[column].std(ddof=0)
        #out_df_min_max_scaled['sig_centrality']=np.abs(np.log10(out_df.pval_1000sample10times+minimum_pval))
    else:
        print('feature pval not find !', '\n')
        out_df_scaled=pd.DataFrame()
  
    return out_df_scaled.copy()


def compare_two_samples(out_df_t0, out_df_t1, selected_index_from=-9):
    ''' compare the features t-vals and sig_centrality between two samples
	by computing Ttest and Euclidean distance between the two vectors
	Output, nodes id, Euclidean distance of each node between two groups based on their features
    '''
    record_pval={}
    #loop in nodes of input dataframe where we assume both dataframe with the same number of nodes
    # and the last 9 columns are selected features and sig_centrality for nodes
    for i in out_df_t0.id.to_list():
        t0_df=out_df_t0[out_df_t0.id==i].copy()
        t1_df=out_df_t1[out_df_t1.id==i].copy()
        if not t0_df.empty and not t1_df.empty:
            t0_vec=t0_df.iloc[:,selected_index_from:].to_numpy().T
            t1_vec=t1_df.iloc[:,selected_index_from:].to_numpy().T
            tval,pval= ttest_ind(a=t0_vec,b=t1_vec,equal_var=False)
            Euclidean_dist=np.sqrt(np.sum((t0_vec-t1_vec)**2))
            if len(tval)>0:
                record_pval[i]=[tval[0],pval[0],Euclidean_dist]
            else:
                record_pval[i]=[np.nan,np.nan,Euclidean_dist]
        else:
            record_pval[i]=[np.nan,np.nan,np.nan] 
    
    new_df=pd.DataFrame(data=record_pval).T
    new_df.columns=['tval','pval','euclidean_dist']
    
    return new_df.copy()


def gaussian_distribution_of_distance(d,sigma=1):
    ''' fit Euclidean distance to a Gaussian probability function
    '''
    d=np.nan_to_num(d)
    prob=np.exp(-d**2)/sigma**2 / np.sum(np.exp(-d**2)/sigma**2)
    
    return prob


def compute_pval2dist(final_df):
    ''' fit Eucliean distance to a Gaussian dsitrution and get its p-value
    '''
    #remove nan from data
    d=final_df[~final_df.euclidean_dist.isna()].euclidean_dist.to_numpy()
    mean_d=np.mean(d)
    sigma_d=np.std(d)
    #fit distance to a gaussian distribution
    final_df['pval2dist']=final_df.euclidean_dist.apply(lambda x: norm.sf(x,mean_d,sigma_d))
    
    return final_df.copy()


def main(in_data_folder, 
         chromosome, 
         resolution, 
         out_data_folder,
         cohort1,
         cohort2,
         pval_cutoff,
         permutation,
         fig_dpi):
    
    bin_str = str(int(resolution/1000)) + 'kb'
    
    if chromosome == 'whole_genome':
        chrom_strs = ['chr'+str(i) for i in range(1,24)]
    else: 
        chrom_strs = chromosome
    chrom_len = len(chrom_strs)
    #tval_str='pval'
    tval_str = 'tval'

    percentage_str = '0'
    rm_global_str = '_rm_'
    
    out_folder = out_data_folder + '/differential_network_analysis/' + bin_str
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
    
    if tval_str=='pval':
        print('Use T-value signed log10(Pvalue) for comparison', '\n')
    else:
        print('Use Tvalue for comparison', '\n')
    
    record_all = []
    #compare two groups
    for chrom_str in chrom_strs:
        in_path = in_data_folder + '/hic_data/' + bin_str
        
        time_str = cohort1
        #t0 find files of enrichment of features in network clustesrs (e.g., histone, gene expression, DNas, et al)
        in_feature_path = os.path.join(in_path, 'hic_community_figures', time_str, chrom_str)
        in_features_pval_file = os.path.join(in_feature_path, chrom_str+'_'+str(permutation)+'_loopPgt_0.01_in_clusters_edges'+rm_global_str+percentage_str+'percent_noWeight_tval.tsv')
        in_nodes_pair_cluster_files = glob.glob(os.path.join(in_feature_path, chrom_str+'_cn*_'+percentage_str+'percent_'+'*info.tsv'))
        
        #t0 find files of significance of centrality of each node
        in_centrality_path = os.path.join(in_path, 'hic_community_data', time_str)
        in_centrality_file = glob.glob(os.path.join(in_centrality_path, chrom_str+'_*_'+percentage_str+'percent'+'_*_communities_info.tsv'))
        #find tval/pval2 and export a dataframe contains both feature enrichment and centritiy in each node 
        out_df_t0 = find_tval2centrality(in_features_pval_file, in_nodes_pair_cluster_files, in_centrality_file, tval_str=tval_str)
        print(out_df_t0, '\n')
        print('find Features for in nodes in %s, here node index strated from 0' %time_str, '\n')
        
        time_str = cohort2
        #t1 find feature files
        in_feature_path = os.path.join(in_path , 'hic_community_figures', time_str, chrom_str)
        in_features_pval_file = os.path.join(in_feature_path, chrom_str+'_'+str(permutation)+'_loopPgt_0.01_in_clusters_edges'+rm_global_str+percentage_str+'percent_noWeight_tval.tsv')
        in_nodes_pair_cluster_files = glob.glob(os.path.join(in_feature_path, chrom_str+'_cn*_'+percentage_str+'percent_'+'*info.tsv'))
        
        #t1 find files of significance of centrality of each node
        in_centrality_path = os.path.join(in_path, 'hic_community_data', time_str)
        in_centrality_file = glob.glob(os.path.join(in_centrality_path, chrom_str+'_*_'+percentage_str+'percent'+'_*_communities_info.tsv'))
        #ind tval/pval2 and export a dataframe contains both feature enrichment and centritiy in each node
        out_df_t1 = find_tval2centrality(in_features_pval_file, in_nodes_pair_cluster_files, in_centrality_file, tval_str=tval_str)
        print(out_df_t1, '\n')
        print('find Features for in nodes in %s, here node index strated from 0' %time_str, '\n')
        
        #only compare nodes of two samples with both have network clustering results
        #and return the Euclidean distance of a node between two samples based on its features
        if not out_df_t0.empty and not out_df_t1.empty:
            new_df = compare_two_samples(out_df_t0, out_df_t1)
            record_all.append(new_df.copy())
        else:
            print(chrom_str,'output df of one sample is empty, skip comparion!', '\n')
            #input('Click for continue')
        if out_df_t0.shape[0] != out_df_t1.shape[0]: 
            #here is only a warning if two dataframe with difference size!
            print(chrom_str, 'output df size of two samples are not the same, be carefule with the results!', '\n')
            print(out_df_t0.shape, '\n')
            print(out_df_t1.shape, '\n')
        
        print('Finish compare ', chrom_str, '\n') 
    
    #final data for all nodes with significant changes between two groups
    final_df = pd.concat(record_all).copy()
    #fit all nodes euclidean distance to a gaussian distribution for obtaining a P-value
    final_df = compute_pval2dist(final_df)
    #selected nodes for further study, here all nodes id start from 0, 1, ...
    
    #selected_df=final_df[final_df.pval2dist<pval_cutoff].copy()
    selected_df = final_df[final_df.pval2dist<pval_cutoff].copy()
    if not selected_df.empty:
        selected_df['id'] = selected_df.index.to_list()
        selected_chroms = selected_df.id.str.split('_',expand=True)[0].unique()
        print(selected_df, '\n')
        print(selected_df.euclidean_dist.min(), '\n')
        print(selected_chroms, '\n')
        
        #export figs of euclidean distance in all nodes
        final_df_nan=final_df.dropna()
        final_df_nan['P-value']='P>=%s'%pval_cutoff
        final_df_nan.loc[final_df_nan['pval2dist'] < pval_cutoff, 'P-value'] = 'P<%s'%pval_cutoff
        final_df_nan = final_df_nan.sort_values('pval2dist', ascending=False)
        min_val = min(final_df_nan[final_df_nan.pval2dist<pval_cutoff].euclidean_dist.values)
        
        mlp.rcParams['font.family'] = 'Arial'
        sns.set(font_scale=1.4)
        sns.set_style("whitegrid", {'axes.grid' : False})
        sns.displot(data=final_df_nan, x='euclidean_dist', hue='P-value')
        plt.axvline(x = min_val, color = 'gray', ls = '--')
        plt.text(min_val+0.1, 10*chrom_len, str(np.round(min_val,2)))
        plt.xlabel('Euclidean Distance')
        plt.ylabel('Number of Nodes')
        plt.tick_params(axis='both', labelsize=12)
        plt.savefig(out_folder+'/distribution_euclidean_distance_pval_'+str(pval_cutoff)+'.jpg', bbox_inches='tight', dpi=fig_dpi)
        
        out_fig_file = os.path.join(out_folder, '%s_vs_%s'%(cohort1, cohort2)+'.jpg')        
        #export selected nodes with significant changes between the two groups
        out_nodes_file = out_fig_file.replace('.jpg','_pval_ls'+str(pval_cutoff)+'_selected_sigNodes.tsv')
        out_df = selected_df[['id','euclidean_dist','pval2dist']].copy()
        out_df.to_csv(out_nodes_file, sep='\t', index=False)
        print(out_nodes_file, '\n')
        # count and and plot the distribution of significant nodes in each chromosome
        num_df = countNum(out_df)
        fig = plt.figure(figsize=(5,3))
        plt.grid(False)
        plt.bar(num_df['chrom'], num_df['count'], color='#EEC1A8')
        
        plt.tick_params(axis='x', labelsize=8)
        plt.tick_params(axis='y', labelsize=12)
        plt.xlabel('Chromosomes', fontsize=16)
        plt.ylabel('Number of Nodes', fontsize=16)
        plt.xticks(rotation = 60)
        plt.savefig(out_folder+'/num_nodes_in_chromosome.jpg', bbox_inches='tight', dpi=fig_dpi)
        print('Export :\n', out_folder+'/num_nodes_in_chromosome.jpg')
        fig.tight_layout()
        plt.close()
        
    else:
        print(chrom_strs, ' no significant nodes find! ', '\n')
        
    












