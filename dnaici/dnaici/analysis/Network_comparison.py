import pandas as pd
import glob
import os
import numpy as np
from scipy.stats import ttest_ind


def true2one(x):
    if x:
        out= 1
    else:
        out= 0
    return out


def find_nodes_pair_cluster_file(in_path0, time_str, si, percentage_str):
    in_feature_path = os.path.join(in_path0 , time_str, si)
    in_nodes_pair_cluster_files_t0 = glob.glob(os.path.join(in_feature_path, si + '_cn*_'+percentage_str+'percent_*info.tsv'))
    print(os.path.join(in_feature_path, si + '_cn*_'+percentage_str+'percent_*info.tsv'), '\n')
    
    return in_nodes_pair_cluster_files_t0.copy()


def find_pair_of_selected_nodes(in_nodes_pair_cluster_files_t0,tmp_selected_nodes):
    '''find a node pair with at least one node that fall in tmp_selected_nodes
    here, in_nodes_pair_cluster_files_t0 with id started from index 1 because they belong to bin/region index
    but in tmp_selected_nodes, the id started from index 0 beause they are nodes of network clusters  
    '''
    all_ids_t0=[]
    nodes_pair_clusters_t0={}
    #convert zero indexed network cluster to 1 indexed bin before searching
    tmp_selected_nodes4test=set([str(int(i)+1) for i in tmp_selected_nodes ])
    for fi in in_nodes_pair_cluster_files_t0:
        tmp_t0_df=pd.read_csv(fi,sep='\t')
        tmp_t0_df['nodes_pair']=tmp_t0_df.id.str.split(':').apply(lambda x: set(x))
        #find a node pair with at least one node that fall in tmp_selected_nodes
        nodes_pair_clusters_t0[fi]=tmp_t0_df[tmp_t0_df.nodes_pair.apply(lambda x: len(x&tmp_selected_nodes4test)>0)].copy()  
        all_ids_t0 += nodes_pair_clusters_t0[fi].id.to_list()
    
    return all_ids_t0.copy(), nodes_pair_clusters_t0.copy()


def add_common_unique4nodes(nodes_pair_clusters_t0,common_in_both_times, only_in_t0,col_name):
    '''assign common and unique lablel for each pair of nodes, which is common in both conditons and unique in current condition, respectively
    '''
    #for fi in in_nodes_pair_cluster_files_t0:
    if isinstance(nodes_pair_clusters_t0,dict):
        for fi in nodes_pair_clusters_t0.keys():
            if not nodes_pair_clusters_t0[fi].empty:
                nodes_pair_clusters_t0[fi]['common']= nodes_pair_clusters_t0[fi][col_name].apply(lambda x: x in common_in_both_times).apply(lambda x: true2one(x) )
                nodes_pair_clusters_t0[fi]['unique']= nodes_pair_clusters_t0[fi][col_name].apply(lambda x: x in only_in_t0).apply(lambda x: true2one(x))
    else:
        #assumes it is a dataframe
        nodes_pair_clusters_t0['common']= nodes_pair_clusters_t0[col_name].apply(lambda x: x in common_in_both_times).apply(lambda x: true2one(x) )
        nodes_pair_clusters_t0['unique']= nodes_pair_clusters_t0[col_name].apply(lambda x: x in only_in_t0).apply(lambda x: true2one(x))
    
    return nodes_pair_clusters_t0.copy()


def find_gene4selected_nodes(nodes_pair_clusters_t0):
    '''input dataframe of clusters to 
    obtain unique gene names in each cluster and all genes in all clusters
    '''
    out_uq_genes_in_clusters={}
    all_genes_t0=[]
    #loop in dictionary of clusters
    for ki in nodes_pair_clusters_t0.keys():
        print(ki, '\n')
        tmp_df=nodes_pair_clusters_t0[ki].copy()
        tmp_df['gene_expression'] = tmp_df['gene_expression'].apply(eval)
        tmp_genes=[]
        #loop in rows of dataframe of a cluster
        for idx,row in tmp_df.iterrows():
            tmp_genes += row['gene_expression']
        #unique name of genes in the cluster
        uniq_genes=list(set(tmp_genes))
        len_gene=len(uniq_genes)
        if len_gene>1:
            #split gene name and its expression level
            uq_df=pd.DataFrame(data=uniq_genes).squeeze().str.split(':',expand=True).squeeze()
            uq_df.columns=['gene','expression']
        elif len_gene==1:
            uq_df=pd.DataFrame(data= uniq_genes[0].split(':') )
            uq_df=uq_df.T
            uq_df.columns=['gene','expression']
        else:
            uq_df=pd.DataFrame()
        #record unique gene names of a cluster and append it to all gene names in the condition
        out_uq_genes_in_clusters[ki]=uq_df.copy()
        all_genes_t0.append(uq_df.copy())
    all_genes_t0_df=pd.concat(all_genes_t0).drop_duplicates(['gene','expression'])
    
    return out_uq_genes_in_clusters.copy(), all_genes_t0_df.copy()
    

def find_common_uniq_in2lists(all_ids_t0, all_ids_t1):
    #compare nodes pair between two samples
    common_in_both= set(all_ids_t0)& set(all_ids_t1)
    only_in_t0=set(all_ids_t0)-  set(all_ids_t1)
    only_in_t1=set(all_ids_t1)- set(all_ids_t0)
    
    return only_in_t0.copy(), only_in_t1.copy(),common_in_both.copy()


def calculate_rr4genes_in2samples(t0_exp_df,t1_exp_df, gonly_in_t0):
    ''' differential analysis of genes in selecte list such as calculating relative ratios of gene expressions between two samples
    '''
    if len(gonly_in_t0)>0:
        t0_exp4gonly_in_t0=t0_exp_df[t0_exp_df.gene.apply(lambda x: x in set(gonly_in_t0))].copy()
        t0_exp4gonly_in_t0.drop_duplicates('gene',inplace=True)
        
        t1_exp4gonly_in_t0=t1_exp_df[t1_exp_df.gene.apply(lambda x: x in set(gonly_in_t0))].copy()
        t1_exp4gonly_in_t0.drop_duplicates('gene',inplace=True)
            
        df2=pd.merge(t0_exp4gonly_in_t0,t1_exp4gonly_in_t0,on='gene',how='inner').copy()
        print('cohort 1 gene,','cohort 2 gene, ','merged gene', '\n') 
        print(t0_exp4gonly_in_t0.shape, t1_exp4gonly_in_t0.shape, df2.shape, '\n')
        tt= ttest_ind(a=df2.zscores_x.astype(float).to_numpy(),b=df2.zscores_y.astype(float).to_numpy())
        print(tt, '\n')
        rr=(df2.zscores_x.astype(float) -df2.zscores_y.astype(float))/(df2.zscores_x.astype(float) +df2.zscores_y.astype(float))*2
        df2['rr']=rr
        df2.sort_values(by='rr',inplace=True)
        print('percentages of genes abs(rr) >0.66',df2[df2.rr.apply(lambda x: np.abs(x)>0.66)].shape[0]/df2.shape[0], '\n')
    else:
        print('No genes find in selected ', gonly_in_t0, '\n')
        df2=pd.DataFrame()
    
    return df2.copy() 


def export_gene_and_nodes(out_folder2, out_gene_file,t0_df2,out_node_file, out_set,column_str):
    '''Export genes and nodes for selected sample
    '''
    out_file2=os.path.join(out_folder2, out_gene_file)
    print(out_file2, '\n')
    t0_df2.to_csv(out_file2,sep='\t',index=False)
    
    out_file3=os.path.join(out_folder2, out_node_file)
    print(out_file3, '\n')
    if len(out_set)>0:
        tmp_df=pd.DataFrame(data=list(out_set))
        tmp_df.columns=[column_str]
    else:
        tmp_df=pd.DataFrame()
    tmp_df.to_csv(out_file3,sep='\t',index=False)
    
    return out_file2, out_file3


def main(in_data_folder, 
         chromosome,
         resolution, 
         out_data_folder,
         cohort1,
         cohort2,
         pval_cutoff):
    
    bin_str = str(int(resolution/1000)) + 'kb'
        
    if chromosome == 'whole_genome':
        chrom_strs = set(['chr'+str(i) for i in range(1,24)])
    else: 
        chrom_strs = set(chromosome)
    
    in_node_folder = out_data_folder + '/differential_network_analysis/' + bin_str
    
    nodes_file = os.path.join(in_node_folder, '%s_vs_%s'%(cohort1,cohort2)+'_pval_ls'+str(pval_cutoff)+'_selected_sigNodes.tsv')
    if not os.path.exists(nodes_file):
        print('No DIEGs found in the input chromosome :( Please move to another chromosome')
    else:
        selected_df = pd.read_csv(nodes_file, sep='\t')
        selected_chroms = set(selected_df.id.str.split('_', expand=True)[0].unique())
        
        co_chroms = selected_chroms & chrom_strs
    
        in_path0 = out_data_folder + '/hic_data/' + bin_str + '/hic_community_figures/'
        percentage_str = '0'
        
        #input gene expression data
        in_exp_f1 = in_data_folder + '/expression_data/mcf7_%s_geneExp.bed'%cohort1
        in_exp_f2 = in_data_folder + '/expression_data/mcf7_%s_geneExp.bed'%cohort1
        t0_exp_df = pd.read_csv(in_exp_f1, sep='\t', header=None)
        t1_exp_df = pd.read_csv(in_exp_f2, sep='\t', header=None)
        column_names = ['chrom','starts','ends','zscores','gene','gene_id','strand','ensemble_id','count']
        t0_exp_df.columns = column_names
        t1_exp_df.columns = column_names
        
        #for each selected chroms find nodes pairs
        record_common_genes_nodes = {}
        record_t0_only_genes_nodes = {}
        record_t1_only_genes_nodes = {}
        record_selected_nodes_in_clusters = {}
        record_size = []
        
        if len(co_chroms) == 0:
            print('No DIEGs found in the input chromosome :( Please move to another chromosome')
        else:
            print('Start analysis...', '\n')
        
            for si in co_chroms:
                print(si, '\n')
                tmp_vec=[]
                #here nodes id starts from 0 index because it comes from network clustering
                tmp_df = selected_df[selected_df.id.str.find(si+'_')>=0].copy()
                tmp_df[['chrom','node']] = tmp_df.id.str.split('_', expand=True)
                
                #find network cluster files
                time_str1 = cohort1
                in_nodes_pair_cluster_files_t0 = find_nodes_pair_cluster_file(in_path0, time_str1, si, percentage_str)
                
                time_str2 = cohort2
                in_nodes_pair_cluster_files_t1 = find_nodes_pair_cluster_file(in_path0, time_str2, si, percentage_str)
                
                #find all nodes pairs involved in the selected nodes, we include all pairs in the clusters with at least one node fall in the selected nodes
                #here tmp selected nodes are zero indexed
                tmp_selected_nodes = set(tmp_df.node.to_list())
                #in_nodes_pair_cluster_files_t, id are 1 indexed because they are bin/region index
                all_ids_t0, nodes_pair_clusters_t0 = find_pair_of_selected_nodes(in_nodes_pair_cluster_files_t0, tmp_selected_nodes)
                all_ids_t1, nodes_pair_clusters_t1 = find_pair_of_selected_nodes(in_nodes_pair_cluster_files_t1, tmp_selected_nodes)
                
                #build a network based on all_ids_t0 or all_ids_t1??
                print('selected nodes, ','node-pair in %s, '%cohort1, 'node-pair in %s'%cohort2, '\n')
                print(len(tmp_selected_nodes), len(all_ids_t0), len(all_ids_t1), '\n')
                #input('Obtained all valid nodes pairs fall into selected nodes, click to continue!\n')
                tmp_vec += [len(tmp_selected_nodes), len(all_ids_t0), len(all_ids_t1)]
                
                #compare nodes pair between two samples
                only_in_t0, only_in_t1, common_in_both_times = find_common_uniq_in2lists(all_ids_t0, all_ids_t1)
                
                #add common or unique label for node pair in both two groups
                nodes_pair_clusters_t0 = add_common_unique4nodes(nodes_pair_clusters_t0, common_in_both_times, only_in_t0, 'id')
                nodes_pair_clusters_t1 = add_common_unique4nodes(nodes_pair_clusters_t1, common_in_both_times, only_in_t1, 'id')
                print('only %s, '%cohort1, 'only %s, '%cohort2, 'in common', '\n')
                print(len(only_in_t0), len(only_in_t1), len(common_in_both_times), '\n')
                #input('Added common or unique label for node pairs in network clusters, click to continue!\n')
                tmp_vec += [len(only_in_t0), len(only_in_t1), len(common_in_both_times)]
                
                #loop in each cluster to get an unique list of genes that associated with selected signifcant nodes pairs that common in both conditinos
                uq_genes_in_clusters_t0, all_genes_t0 = find_gene4selected_nodes(nodes_pair_clusters_t0)
                uq_genes_in_clusters_t1, all_genes_t1 = find_gene4selected_nodes(nodes_pair_clusters_t1)
                #find common and unique genes in the two condtioins
                gonly_in_t0, gonly_in_t1, gcommon_in_both_times = find_common_uniq_in2lists(all_genes_t0.gene.to_list(), all_genes_t1.gene.to_list())
                print('%s genes, '%cohort1, '%s genes, '%cohort2, 'common genes', '\n')
                print(len(gonly_in_t0),len(gonly_in_t1),len(gcommon_in_both_times), '\n')
                #input('Obtained gene names for unique in t0 or t1, and common in both t0 and t1, respectively, click for continue!\n')
                tmp_vec += [len(gonly_in_t0), len(gonly_in_t1), len(gcommon_in_both_times)]
                
                tmp_size_df = pd.DataFrame(data=tmp_vec,
                                           index=['selected nodes','node-pair %s'%cohort1,'node-pair %s'%cohort2,'node-pair %s only'%cohort1,'node-pair %s only'%cohort2,'node-pair in common','%s only genes'%cohort1,'%s only genes'%cohort2,'genes in common'],
                                           columns=[si])
                record_size.append(tmp_size_df.copy())
                
                t0_df2 = calculate_rr4genes_in2samples(t0_exp_df.copy(), t1_exp_df.copy(), gonly_in_t0)
                #input('t0 only, click \n')
                
                t1_df2 = calculate_rr4genes_in2samples(t0_exp_df.copy(), t1_exp_df.copy(), gonly_in_t1)
                #input('t1 only, click \n')
                
                common_df2 = calculate_rr4genes_in2samples(t0_exp_df.copy(), t1_exp_df.copy(), gcommon_in_both_times)
                #input('common genes, click \n')
                
                out_folder = os.path.join(in_node_folder, si)
                if not os.path.exists(out_folder):
                    os.mkdir(out_folder)
                
                #tmp_selected_nodes start from 0 index, which are selected significant nodes between two samples.
                #nodes_pair_clusters_ start from 1 index, which are dataframe of node pairs that meet one of selected significant node 
                #tmp_selected_nodes
                record_selected_nodes_in_clusters[si] = (tmp_selected_nodes.copy(), nodes_pair_clusters_t0.copy(), nodes_pair_clusters_t1.copy())
                out_file1 = os.path.join(out_folder, 'selected_nodes_in_clusters.tsv')
                print(out_file1, '\n')
                tmp_df = pd.DataFrame(data=list(tmp_selected_nodes))
                tmp_df.columns = ['selected_nodes_in_clusters']
                tmp_df.to_csv(out_file1,sep='\t', index=False)
                
                #t0_df2 is a dataframe for genes only in nodes pairsof t0,  only_in_t0 is a list of nodes pairs only in t0 (starts from index 1),  gonly_in_t0 is a list of genes only in t0
                #t0_df2, only_in_t0
                record_t0_only_genes_nodes[si] = (t0_df2.copy(), only_in_t0.copy(), gonly_in_t0.copy())
                export_gene_and_nodes(out_folder, '%s_only_genes.tsv'%cohort1, t0_df2, '%s_only_nodes_pairs.tsv'%cohort1, only_in_t0, '%s_only_nodes_pairs'%cohort1)
                
                record_t1_only_genes_nodes[si] = (t1_df2.copy(), only_in_t1.copy(), gonly_in_t1.copy())
                export_gene_and_nodes(out_folder, '%s_only_genes.tsv'%cohort2, t1_df2, '%s_only_nodes_pairs.tsv'%cohort2, only_in_t1, '%s_only_nodes_pairs'%cohort2)
                
                record_common_genes_nodes[si] = (common_df2.copy(), common_in_both_times.copy(), gcommon_in_both_times.copy())
                print(common_df2, common_in_both_times, '\n')
                export_gene_and_nodes(out_folder, 'common_genes.tsv', common_df2, 'common_nodes_pairs.tsv', common_in_both_times, 'common_nodes_pairs')
        
        record_size_df = pd.concat(record_size,axis=1).copy()
        out_file = os.path.join(in_node_folder, 'summary_of_nodes_count.tsv')
        print(out_file, '\n')
        record_size_df.to_csv(out_file, sep='\t')
    





