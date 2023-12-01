#here contains funcations for data analysis
import pandas as pd
import numpy as np
from itertools import product
import os
from sklearn.cluster import KMeans
import scipy.stats as stats
import matplotlib.pyplot as plt


def sort_dict(deg_cen):
  temp_dict = {}
  for w in sorted(deg_cen, key=deg_cen.get, reverse=True):
     temp_dict[w] = deg_cen[w]
  return temp_dict.copy()


def read_all_data_files(out_data_folder, bin_str, feature_str, chrom_str, time_str, hic_f, chrom_region_f):
  '''
     read all data files for the data analysis 
     and return data matrix in dataframe, where the index of data matrix is correponding to the exported all_files list
     Note: here all data matrix are based on Hi-C interacton matrix such as adj-matrix, if there is not an interaction then the corresponding 
           fetuare will set to zero in the data matrix and dataframe. 
     output:
	all_matrix, contains all feaure elements in adj-matrix format where contain non-interaction elements
        new_region_pair_df, contains all feature elements for pair-wise interactions in dataframe formnat where non-interaction elements are removed         
  '''
  #input files
  exp_f = out_data_folder + '/expression_data/' + bin_str + '/out_plots/' + time_str + '/' + chrom_str +'_geneExp_'+ feature_str + '_zscore.matrix'
  dnas_f = out_data_folder + '/nucleosome_density_data/' + bin_str + '/out_plots/' + time_str + '/' + chrom_str +'_DNas_'+ feature_str + '_zscore.matrix'
  ctcf_f = out_data_folder + '/histone_data/' + bin_str + '/out_plots/' + time_str + '/' + chrom_str +'_ctcf_'+ feature_str + '_zscore.matrix'
  h3k4me1_f = out_data_folder + '/histone_data/' + bin_str + '/out_plots/' + time_str + '/' + chrom_str +'_h3k4me1_'+ feature_str + '_zscore.matrix'
  h3k4me3_f = out_data_folder + '/histone_data/' + bin_str + '/out_plots/' + time_str + '/' + chrom_str +'_h3k4me3_'+ feature_str + '_zscore.matrix'
  h3k27ac_f = out_data_folder + '/histone_data/' + bin_str + '/out_plots/' + time_str + '/' + chrom_str +'_h3k27ac_'+ feature_str + '_zscore.matrix'
  h3k27me3_f = out_data_folder + '/histone_data/' + bin_str + '/out_plots/' + time_str + '/' + chrom_str +'_h3k27me3_'+ feature_str + '_zscore.matrix'
  h3k9me3_f = out_data_folder + '/histone_data/' + bin_str + '/out_plots/' + time_str + '/' + chrom_str +'_h3k9me3_'+ feature_str + '_zscore.matrix'
  
  #file order
  all_files=[hic_f, exp_f,  dnas_f, ctcf_f, h3k4me1_f, h3k4me3_f, h3k27ac_f, h3k27me3_f, h3k9me3_f]

  #read feature files
  all_matrix=[]
  all_columns=['row_id','column_id']
  for ai in all_files:
    all_columns.append(os.path.basename(ai).replace('.matrix',''))
    all_matrix.append( np.loadtxt(ai) )
    #print(ai, all_matrix[-1].shape)

  #for hic-matrix only read up-diagnoal elements
  all_matrix[0]=np.triu(all_matrix[0])
 
  #read region files
  region_df=pd.read_csv(chrom_region_f,sep='\t',header=None)

  #make a matrix for pair-wise interaction positions
  list1=[str(i) for i in range(1,region_df.iloc[-1,3]+1)]
  list2=list1.copy()
  product_lists=list(product(list1,list2))
  region_pair_df=pd.DataFrame(data=product_lists).copy()

  #convert matrix to vector and add to pd
  for i in range(0,len(all_matrix)):
    #all matrix only read the updiagnola part elements
    tmp_matrix=np.triu(all_matrix[i]).copy()
    tmp_vector=tmp_matrix.reshape(np.product(tmp_matrix.shape),1).copy()
    region_pair_df[2+i]=tmp_vector.copy()

  #assign all feature matrix to a data frame before filtering rows with zero enteries
  region_pair_df.columns=all_columns

  #filtering rows with all zero entries
  new_region_pair_df= region_pair_df[~(np.abs(region_pair_df.iloc[:,2:]).sum(axis=1)==0)].copy()
  new_region_pair_df.index=new_region_pair_df['row_id']+':'+new_region_pair_df['column_id']
  new_region_pair_df.drop(['row_id','column_id'],axis=1,inplace=True)

  #filtering rows with zeros in hi-c interactions and assume the zero column is the hic-matrix
  new_region_pair_df=new_region_pair_df[new_region_pair_df.iloc[:,0]!=0]

  return all_matrix.copy(), new_region_pair_df.copy(), all_files


def read_data_files_without_interaction(in_data_folder, bin_str, chrom_str, time_str):
  ''' read full data of 8 features (expression, Dnas, histones) in a chromosome specific and window bin speicific format (e.g., _500kb_regions)
      return loaded dataframe in a dictoanry
  '''
  #input files
  exp_f = in_data_folder + '/expression_data/' + bin_str + '/out_data/' + time_str + '/' + chrom_str + '_' + bin_str + '_regions_geneExp_array.bed'
  dnas_f = in_data_folder + '/nucleosome_density_data/' + bin_str + '/out_data/' + time_str + '/' + chrom_str + '_' + bin_str + '_regions_DNas_array.bed'
  ctcf_f = in_data_folder + '/histone_data/' + bin_str + '/out_data/' + time_str + '/' + chrom_str + '_' + bin_str + '_regions_ctcf_array.bed'
  h3k4me1_f = in_data_folder + '/histone_data/' + bin_str + '/out_data/' + time_str + '/' + chrom_str + '_' + bin_str + '_regions_h3k4me1_array.bed'
  h3k4me3_f = in_data_folder + '/histone_data/' + bin_str + '/out_data/' + time_str + '/' + chrom_str + '_' + bin_str + '_regions_h3k4me3_array.bed'
  h3k27ac_f = in_data_folder + '/histone_data/' + bin_str + '/out_data/' + time_str + '/' + chrom_str + '_' + bin_str + '_regions_h3k27ac_array.bed'
  h3k27me3_f = in_data_folder + '/histone_data/' + bin_str + '/out_data/' + time_str + '/' + chrom_str + '_' + bin_str + '_regions_h3k27me3_array.bed'
  h3k9me3_f = in_data_folder + '/histone_data/' + bin_str + '/out_data/' + time_str + '/' + chrom_str + '_' + bin_str + '_regions_h3k9me3_array.bed'
  #file order
  all_files = [exp_f, dnas_f, ctcf_f, h3k4me1_f, h3k4me3_f, h3k27ac_f, h3k27me3_f, h3k9me3_f]
  #read feature files
  all_data = {}
  for ai in all_files:
    all_data[os.path.basename(ai).replace('.bed','')] = pd.read_csv(ai, sep='\t').copy() 
  
  return all_data.copy()


def plot_elbow(chrom_str, num_of_clusters,new_region_pair_df,out_path):
  #elbow curve for kmeans clustering
  print('Estimate the number of clusters ....')
  plt.clf()
  Nc = range(1, num_of_clusters)
  kmeans = [KMeans(n_clusters=i) for i in Nc]
  score = [kmeans[i].fit(new_region_pair_df).score(new_region_pair_df) for i in range(len(kmeans))]
  plt.plot(Nc,score)
  plt.vlines([5,10],plt.gca().get_yticks().min(), plt.gca().get_yticks().max(),colors='red')
  plt.xlabel('Number of Clusters')
  plt.ylabel('Score')
  plt.title('Elbow Curve')
  out_fig=os.path.join(out_path, chrom_str+'_Elbow_curve_kmeans.jpg')
  print(out_fig)
  plt.savefig(out_fig)
  return out_fig


def find_cluster_label4nodes(chrom_str,num_of_clusters4network,nodes_cluster_df,new_region_pair_df,out_path,edge_feature_str='',plotElbow=True):
  '''
    find network clustering label for all nodes in pair-wise interactions dataframe.
    In the dataframe, the index is a pair of nodes, the columns are histone or other features in a pair of interactions,
    plut ID, row_i,col_j position in adjancy matrix (position starts from 1, )
    For each selected clustering of networks, do kmeans elbow calculation
  '''

  #find a set of nodes for each cluster
  node2clusters={}
  for i in range(0,num_of_clusters4network):
     node2clusters[i]= nodes_cluster_df.nodes[nodes_cluster_df.labels==i].to_list()

  #renmove clusters with only one node
  new_node2clusters=node2clusters.copy()
  for i in new_node2clusters.copy().keys():
      if len(new_node2clusters[i])<=1:
         del new_node2clusters[i]

  #extract all interactions within the same clusters but ignore cross clusters interactions
  new_region_pair_df['ID']=new_region_pair_df.index.to_list()
  new_region_pair_df[['row_i','row_j']]=new_region_pair_df.ID.str.split(':',expand=True).astype(int)-1

  record_clustered_sub_df={}
  #here cli is from zero index
  print('Number of nodes in valid clusters:')
  for cli in new_node2clusters.keys():
    #assume both nodes in the same cluster
    sub_df=new_region_pair_df[( new_region_pair_df.row_i.apply(lambda x: x in new_node2clusters[cli]))
       & (new_region_pair_df.row_j.apply(lambda x: x in new_node2clusters[cli]))].copy()
    print(cli, sub_df.shape)
    sub_df.drop(['ID','row_i','row_j'], axis=1,inplace=True)
    num_of_sub_clusters=sub_df.shape[0]
    if num_of_sub_clusters>=30:
       num_of_sub_clusters=30
    if plotElbow:
        plot_elbow(chrom_str+'_'+str(cli+1)+edge_feature_str,num_of_sub_clusters,sub_df.copy(), out_path)
    record_clustered_sub_df[cli]=sub_df.copy()

  return node2clusters, record_clustered_sub_df.copy()


def sort_clustereLabel_by_markers(record_clustered_sub_df,selected_markers,color_clip_value,clusterSize_clip_value=20):
  ''' sort cluster labels based on mean profiles of selected markers such as activate markers.
    return sorted dict with new keys for new cluster labels
  '''
  #sort record_clustered_sub_df by the mean of h3k4me1, h3k3me3 and h3k27ac
  #set a maximum values in data to avoid outlier values may affect the mean of data
  maximum_clip_value=color_clip_value
  label2mean_h3k={}
  tmp_df_list=[]
  #loop in old cluster label
  for ki in record_clustered_sub_df.keys():
     tmp_df= record_clustered_sub_df[ki].copy()
     tmp_df[tmp_df>maximum_clip_value]=maximum_clip_value
     #loop in selected markers
     tmp_df_list=[]
     for si in selected_markers:
          tmp_df_list.append(tmp_df[tmp_df.columns[tmp_df.columns.str.find(si)>=0]].copy())
     joined_df=pd.concat(tmp_df_list,axis=1).copy()
     #record mean profiles of selected markers, and size of cluster
     label2mean_h3k[ki]=[joined_df.mean(axis=1).mean(), joined_df.shape[0]]

  #sort cluster label based on size of cluster and the mean of active markers
  label2mean_df=pd.DataFrame.from_dict(label2mean_h3k).T.sort_values(by=[0,1], ascending=True).copy()
  #reset the idex, the new index is the new sorted cluster label but the old index column is the old cluster label
  label2mean_df.reset_index(inplace=True)
  #print(label2mean_df) 
  #remove cluster with small size
  #print(label2mean_df)
  label2mean_df.drop(label2mean_df.index[label2mean_df[1]< clusterSize_clip_value],inplace=True)
  label2mean_df.reset_index(inplace=True)
   
  #move df to its new cluster label in dict
  sorted_record_clustered_sub_df={}
  for idx, row in label2mean_df.iterrows():
     #assigen new cluster label to the df
     sorted_record_clustered_sub_df[idx]=record_clustered_sub_df[row['index']].copy()
  #print(label2mean_df)
  return label2mean_df.copy(), sorted_record_clustered_sub_df.copy()


def find_cutoff4percentage(A,percentage):
  '''
    filtering matrix based on percentage of data
    first, sort data from minimum to maximum,
    then, based on selected lowest percentage to find the cutoff values for A matrix
  '''
  #find all non zeros in matrix A and sort them from minimum to maximum
  if (not np.isnan(percentage)) and (percentage !=0):
    non_zero_values=np.sort(A[A!=0])
    len_of_values=non_zero_values.shape[0]
    percentage_cutoff=non_zero_values[int(len_of_values*percentage):][0]
  else:
    percentage_cutoff=np.nan
  return percentage_cutoff


def call_permutation_test4a_feature(args):
    '''Here, we perform the permutation test for one of features in a cluster
    Input:
	ci, chrom string such as chr1
	fi, feature string such as _ctcf_
	tmp_df, dataframe of all rows in a cluster
	tmp_df_ids,	a set of all row IDs (e.g., 10:22) in tmp_df 
	out_record_features_dict, genome wide feaure dataframe where pair-wise nodes (e.g., 10:22) as ID
	loop4rtest,	number of loop for permutation tests
	pval_cutoff,	p-value cutoff for samples in a cluster is significantly higher/lower than randomly selected samples
	feature_str,	featture string for an edge such as maxWeight or mean
	ttest_column_str,	column label for ttest such as tval, pval
    '''
    ci, fi, tmp_df, tmp_df_ids, out_record_features_dict, loop4test, pval_cutoff, feature_str, ttest_column_str = args
    record_fi_tval_pval = {}
    #load data from a cluster and select column for feature fi
    tmp_df_fi_vect = tmp_df.loc[:,tmp_df.columns.str.contains(fi)].to_numpy()
    len_of_vect = len(tmp_df_fi_vect)
    #load full genome-wide data
    tmp_fi_data = out_record_features_dict[fi].copy()
    #remove rows already in tmp_df_fi_vect
    select_fi_data = tmp_fi_data[ ~( (tmp_fi_data.id.apply(lambda x: x in tmp_df_ids) ) & (tmp_fi_data.bin_chrom == int(ci.replace('chr','')) ))].copy()
    print(fi)
    record_ttest = []
    for i in range(0, loop4test):
        sampled_fi_data = select_fi_data.sample(len_of_vect).copy()
        if feature_str == 'maxWeight':
            sampled_fi_vect = sampled_fi_data.loc[:,'bin_max'].to_numpy()
        else:
            sampled_fi_vect = sampled_fi_data.loc[:,'bin_mean'].to_numpy()
        sampled_fi_vect = sampled_fi_vect.reshape(len_of_vect,1)
        #do ttest between selecte feature vector from a cluster and the randoḿ samples
        tval, pval = stats.ttest_ind(a=tmp_df_fi_vect, b=sampled_fi_vect, equal_var=False)
        record_ttest.append([tval[0], pval[0]])
    #
    ttest_df = pd.DataFrame(data=record_ttest)
    #ttest_df.columns=['tval','pval']
    ttest_df.columns = ttest_column_str
    expected_pval = 1 - ttest_df[ttest_df.pval<pval_cutoff].shape[0]/loop4test
    expected_tval = ttest_df.tval.mean()
    record_fi_tval_pval[fi] = [expected_tval, expected_pval]
    
    return record_fi_tval_pval.copy()

