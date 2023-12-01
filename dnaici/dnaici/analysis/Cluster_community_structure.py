import pandas as pd
from dnaici.tools.data_functions import read_all_data_files, find_cluster_label4nodes, sort_clustereLabel_by_markers
from dnaici.tools.plot_tools import plot_clustered_heatmap, plot_network_clusteringLabel_heatmap
import numpy as np
import os
import warnings
warnings.filterwarnings("ignore")

def find_sort_plot_network_clusterings(chrom_str, feature_str, all_matrix, in_nodes_clusters_f, is_plotElbow, edge_df, 
                                       new_region_pair_df, edge_feature_str, out_path, selected_markers,
                                       minCluster_size, color_clip_value, record_edges_in_valid_clusters, fig_dpi):
    ''' This function is used to view cluster of each clustered network and export a heatmap plot of their hi-C adj-matrix,
    k-means clustering of histone data, where the order of clustered networks are sorted by represssion histone markers
    in descending order. Finally, a summary statistic of number of edges, percentage of edges in valied clusters,
    and modulirity score for clustered networks in each chromosome are listed.
    Input:
    in_nodes_clusters_f, 
        a file contains cluster label and the corresponding nodes .
    is_plotElbow, 
        True/False, whether to plot a elbow for each k-means clustering of network cluster. If true then it may get slow..
    e_df,
        a data frame of total number of edges before and after filtering and the modularity score of network clustering in each 
        omosome. It is inputed from 	in_edge_file, which is an exported file by make_inputFile4clustering_homer.py
    new_region_pair_df, 
        dataframe contains all histone and hi-C data in a chromosome, which is loaded by read_all_data_files from 
		multiple files with predefined winodow bin size such as 500kb and hi-c adj-matrix in hic_f. It will be updated in 
		this function by find_cluster_label4nodes() and add edge position in the adj-matrix.
	edge_feature_str, 
		a string of edge such as edges_rm_str + percent_str+ percent_noWeight 
    out_path,
		out file path for exported data
	selected_markers,		
		a list of selected histone names (e.g., ['h3k27me3', 'h3k9me3'] ) for sort the order of clusters of network clustering.
	minCluster_size,
		minimum size of edges in a valid network cluster such as minCluster_size=20
	color_clip_value,
		a clip value for maximum value in color coder for the heatmap such as equal to 3
	record_edges_in_valid_clusters,
		a dictionary of number of edges in valied clusters in each chromosome.
    Output:
	out_fig,
		a list of exported figure names
 	cluster_matrix, 
		adj-matrix for hic-interactions
	record_edges_in_valid_clusters,
		a dictionary of number of edges in valied clusters in each chromosome
    '''
    #read node cluster label file
    #in_nodes_clusters_f=chrom_str + edge_feature_str + '_4communities_class.tsv'
    #in_nodes_clusters_f=os.path.join(in_path,in_nodes_clusters_f)
    print('Read:\n', in_nodes_clusters_f, '\n')
    nodes_cluster_df = pd.read_csv(in_nodes_clusters_f, sep='\t')
        
    cutoff_value4edges = edge_df.loc[chrom_str, 'cutoff_value4edge']    
    #filering full edges data based on cutoff 
    if (not np.isnan(cutoff_value4edges)) and (cutoff_value4edges !=0) :
        new_region_pair_df=new_region_pair_df[new_region_pair_df[chrom_str+'_zscore']>cutoff_value4edges].copy()
        
    #find network clustering label for each nodes and subgraph based on network clustering,by using Hi-c data before the filtering
    num_of_cluster4network=len(nodes_cluster_df.labels.unique())
    #here edges/interaction across the clusters will not be shown
    node2clusters, tmp_record_clustered_sub_df = find_cluster_label4nodes(chrom_str, num_of_cluster4network, nodes_cluster_df, new_region_pair_df,
                                                                          out_path, edge_feature_str=edge_feature_str, plotElbow=is_plotElbow)
    
    label2mean_df, record_clustered_sub_df= \
          sort_clustereLabel_by_markers(tmp_record_clustered_sub_df,selected_markers,color_clip_value,clusterSize_clip_value=minCluster_size)
    
    total_clusters=set(nodes_cluster_df.labels.unique())
    valid_clusters=set(record_clustered_sub_df.keys())
    print('Total network clusters: ', len(total_clusters))
    print('Valid network clusters (clusters with more than 1 nodes and > %d pair of interactions ): %d' %(minCluster_size, len(valid_clusters)))
    total_edges_in_valid_clusters=label2mean_df.iloc[:,-1].sum().astype(int)
    print('Total number of edges in valid clusters: ', total_edges_in_valid_clusters, '\n')
    record_edges_in_valid_clusters[chrom_str]=[total_edges_in_valid_clusters, len(valid_clusters), len(total_clusters)]
    #select one of clustered networks for doing kmeans clustering then show its heatmap
    #here ci is from zero index
    #this shall be inputed manually
    num2clusters=[5,5]
    max_clusters=len(valid_clusters)
    if max_clusters >len(num2clusters):
        for i in range(len(num2clusters),max_clusters+1):
            num2clusters.append(3)
    
    out_fig_data_file=[]
    for ci in record_clustered_sub_df.keys():
        num_of_clusters=num2clusters[ci]
        color_type='GWR'
        #export figure of clustered heatmap
        out_fig= plot_clustered_heatmap(chrom_str+'_cn'+str(ci+1),feature_str,num_of_clusters, 
                                        record_clustered_sub_df[ci], color_clip_value,color_type,
                                        out_path,edge_feature_str=edge_feature_str,fig_dpi=fig_dpi)
        #export data of clustered heatmap, here index is a pair of nodes in an edge
        out_fig_data_path, out_fig_data_file0 =os.path.split(out_fig)
        out_fig_data_file0=out_fig_data_file0.replace('.jpg','_data.tsv')
        out_fig_data_file0=os.path.join(out_fig_data_path, out_fig_data_file0)
        record_clustered_sub_df[ci].to_csv(out_fig_data_file0,sep='\t')
        print('Output: \n', out_fig_data_file0, '\n')
        out_fig_data_file.append(out_fig_data_file0)
    #plot network clustering label in a heatmap based on Hi-C interaction matrix
    #assume the first one is Hi-C matrix, then build a cluster matrix based on network clustering labels
    #here the input matrix is the original hi-c matrix which will be filtered latter
    a_matrix=np.triu(all_matrix[0]).copy()
    
    #filtering hi-c interaction matrix based on cutoff
    if (not np.isnan(cutoff_value4edges)) and (cutoff_value4edges !=0) :
        a_matrix[a_matrix<cutoff_value4edges]=0
    
    #plot cluster label heatmap where only both nodes located in the same cluster are shown!
    out_fig, cluster_matrix= plot_network_clusteringLabel_heatmap(chrom_str,len(valid_clusters), a_matrix, new_region_pair_df , 
                                                                  record_clustered_sub_df,out_path,edge_feature_str=edge_feature_str,fig_dpi=fig_dpi)
    
    #export labels 2 new network node id after soring by features
    out_fig_node_file=out_fig.split('.jpg')[0] + '_ID2nodes.tsv'
    label2mean_df.columns=['sorted_id','original_network_cluster_node','mean_feature_score','number_of_edges']
    label2mean_df.number_of_edges=label2mean_df.number_of_edges.astype(int)
    label2mean_df.mean_feature_score=label2mean_df.mean_feature_score.map('{:5.3f}'.format)
    label2mean_df['sorted_network_cluster_node_cn']=label2mean_df.index.to_list()
    label2mean_df.to_csv(out_fig_node_file,sep='\t',index=False)
    print('Output: \n', out_fig_node_file, '\n')
    
    return out_fig, out_fig_data_file, cluster_matrix, record_edges_in_valid_clusters


def export_combined_df4valid_clusters(record_edges_in_valid_clusters, edge_df, out_path, minCluster_size, edges_rm_str, percent_str):
    '''Input:
	record_edges_in_valid_clusters, 
		a dictionary of chromosome specific number of edges and valid clusters
	edge_df,
		a data frame of total number of edges before and after filtering and the modularity score of network clustering in each 
        chromosome. It is inputed from in_edge_file, which is an exported file by make_inputFile4clustering_homer.py
	out_path,	
		output file path for exporting combined data from both record_edges_invalid_clusters and edge_df 
    '''
    valid_edges_df = pd.DataFrame.from_dict(record_edges_in_valid_clusters).T.copy()
    valid_edges_df.columns = ['edges_in_valid_clusters','num_of_valid_clusters','total_num_of_clustesr']
    combined_df = pd.merge(edge_df,valid_edges_df,left_index=True, right_index=True).copy()
    combined_df['percentage_edges_in_valid_clusters'] = combined_df.edges_in_valid_clusters/combined_df.edges_after_filtering
    out_file= os.path.join(out_path, 'min'+str(minCluster_size)+'edges4VaildCluster_num_of'+edges_rm_str+percent_str+ 'percentage.edges_in_network.tsv')
    print('Output: \n', out_file, '\n')
    combined_df.cutoff_value4edge=combined_df.cutoff_value4edge.map('{:5.3f}'.format)
    combined_df.percentage=combined_df.percentage.map('{:5.3f}'.format)
    combined_df.percentage_edges_in_valid_clusters=combined_df.percentage_edges_in_valid_clusters.map('{:5.3f}'.format)
    combined_df.to_csv(out_file,sep='\t')
    
    return out_file , combined_df.copy()


def main(in_data_folder, 
         cohort, 
         chromosome, 
         resolution, 
         out_data_folder,
         minCluster_size,
         fig_dpi):

    bin_str = str(int(resolution/1000)) + 'kb'
    
    if chromosome == 'whole_genome':
        chrom_strs = ['chr'+str(i) for i in range(1,24)]
    else: 
        chrom_strs = chromosome
    #dataframe format
    pd.options.display.float_format = '{:.2f}'.format    
    #plot parameters
    is_plotElbow = False
    color_clip_value = 3
    isFind_global_cutoff = False
    percent_str = '0' 
    feature_str = 'maxWeight'
    
    if isFind_global_cutoff:
        edges_rm_str='_edges_rmGlob_'
    else:
        edges_rm_str='_edges_rm_'

    in_path = in_data_folder + '/hic_data/' + bin_str + '/hic_community_data/' + cohort
    #this edge file is made by make_inputFile4clustering_homer.py
    in_edge_file = os.path.join(in_path, 'num_of'+edges_rm_str+percent_str+'percentage.edges')
    #hic files path
    hic_file_path = in_data_folder + '/hic_data/' + bin_str + '/hic_interaction_bed/' + cohort
    #hic file name postprefix string
    hic_f_postprefix_file_name = '_meanZscore.matrix'
    chrom_region_postprefix_file_name = '_' + bin_str + '_regions.tsv'
    #output path
    out_path = out_data_folder + '/hic_data/' + bin_str +  '/hic_community_figures/' + cohort
    if not os.path.exists(out_path):
        os.makedirs(out_path)
        print('Create output folder:\n', out_path, '\n')
    #repressors for sort
    selected_markers=['h3k27ac','h3k4me1']
    #load number of edges in each chrom and the corresponding cutoff values
    print('Read:\n', in_edge_file, '\n')
    edge_df = pd.read_csv(in_edge_file, sep='\t', index_col=0)

    record_edges_in_valid_clusters={}
    record_cluster_matrix={}
    record_out_fig=[]
    
    for chrom_str in chrom_strs:
        print('- - - - - - - - investigating %s - - - - - - - -' %chrom_str, '\n')
        edge_feature_str = edges_rm_str + percent_str + 'percent_noWeight'
        out_path1 = os.path.join(out_path, chrom_str)
        if not os.path.exists(out_path1):
            os.makedirs(out_path1)
            print('Creat output folder for %s:\n' %chrom_str, out_path1, '\n')
        
        hic_f = os.path.join(hic_file_path, chrom_str + hic_f_postprefix_file_name)
        chrom_region_file = os.path.join(hic_file_path, chrom_str + chrom_region_postprefix_file_name)
        #read all data files with full matrix, here we load origan hi-c matrix before the filtering and clustering
        all_matrix, new_region_pair_df, all_in_files = read_all_data_files(out_data_folder, bin_str, feature_str, chrom_str, 
                                                                           cohort, hic_f, chrom_region_file)
        #read node cluster information from network clustering
        in_nodes_clusters_f = chrom_str + edge_feature_str + '_4communities_class.tsv'
        in_nodes_clusters_f = os.path.join(in_path, in_nodes_clusters_f)
        
        out_fig, out_fig_data_file, cluster_matrix, record_edges_in_valid_clusters = find_sort_plot_network_clusterings(
            chrom_str, feature_str, all_matrix, in_nodes_clusters_f, is_plotElbow, edge_df, new_region_pair_df, edge_feature_str,
            out_path1, selected_markers, minCluster_size, color_clip_value , record_edges_in_valid_clusters, fig_dpi)
        record_cluster_matrix[chrom_str]= cluster_matrix.copy()
        record_out_fig.append(out_fig)
    
    #combine two dataframe and export the edge statistics after calculation
    if edge_df.shape[0] == len(record_edges_in_valid_clusters):
        out_file, combined_df = export_combined_df4valid_clusters(record_edges_in_valid_clusters, edge_df, out_path, 
                                                                  minCluster_size, edges_rm_str, percent_str)
    else:
        print('Output: \n', edge_df, '\n')
        print('Output: \n', record_edges_in_valid_clusters, '\n')



