import networkx as nx
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
import os
import json
from networkx.readwrite import json_graph


def read_json_file(filename):
    
    with open(filename) as f:
        js_graph = json.load(f)
    
    return json_graph.cytoscape_graph(js_graph)


def sort_dict(deg_cen):
    temp_dict = {}
    for w in sorted(deg_cen, key=deg_cen.get, reverse=True):
        temp_dict[w] = deg_cen[w]
    
    return temp_dict.copy()


def main(in_data_folder, 
         cohort, 
         chromosome, 
         resolution, 
         out_data_folder,
         permutation):
    
    bin_str = str(int(resolution/1000)) + 'kb'
    
    if chromosome == 'whole_genome':
        chrom_strs = ['chr'+str(i) for i in range(1,24)]
    else: 
        chrom_strs = chromosome

    percentage_str='0'
    edges_rm_str='_edges_rm_'
    
    file_path = in_data_folder + '/hic_data/' + bin_str + '/hic_community_data/' + cohort    
    
    print(file_path, '\n')
    
    cluster_label_g = {}
    chrom_g = {}
    degree_centrality = {}
    sorted_deg_cen = {}
    for ci in chrom_strs:
        #load network and network cluster label files
        f1 = ci + edges_rm_str + percentage_str + 'percent_noWeight_4communities_class.tsv'
        f2 = ci + edges_rm_str + percentage_str + 'percent_noWeight_network.json'
        #load cluster label file
        in_df = pd.read_csv(os.path.join(file_path, f1), sep='\t')
        cluster_label_g[ci] = in_df.copy()
        #load cytoscape network json file
        filename = os.path.join(file_path, f2)
        g = read_json_file(filename)
        chrom_g[ci] = g
        #compute degree of centrality of graph in each chromosome
        degree_centrality[g] = nx.degree_centrality(g)
        # let us now sort the degree centrality measure and identify the important nodes.
        sorted_deg_cen[ci] = sort_dict(degree_centrality[g])
    
    #collect all nodes degree centrity then assign a p-value to each of them
    all_df_list = []
    #remove duplicated nodes label
    for ki in sorted_deg_cen.keys():
        tmp_dict = sorted_deg_cen[ki]
        tmp_df = pd.DataFrame.from_dict(tmp_dict, orient='index')
        tmp_df.insert(0, 'id', tmp_df.index.to_list())
        tmp_df.id = ki + '_' + tmp_df.id.astype(str)
        tmp_df.columns = ['id', 'centrality']
        #remove duplicated id
        tmp_list = []
        for ii in tmp_df.id.unique():
            sel_df = tmp_df[tmp_df.id==ii].copy()
            if sel_df.shape[0] > 1:
                sel_df2 = sel_df[sel_df.centrality>0].copy()
                tmp_list.append(sel_df2.copy())
            else:
                #print(sel_df)
                tmp_list.append(sel_df.copy())
        #after removing duplicated nodes label then add this new df to list
        all_df_list.append(pd.concat(tmp_list).copy())
    
    #generate a df from list
    cen2genome_df=pd.concat(all_df_list)
    
    record_values = []
    loop = 10
    for idx,row in cen2genome_df.iterrows():
        #copy all values and remove the one is testing now  
        tmp_df = cen2genome_df.copy()
        tmp_df.drop(index=idx, inplace=True)
        tmp_pval = 0
        for li in range(0, loop):
            random_df = tmp_df.sample(n=permutation)
            tmp_pval += random_df[random_df.centrality>row.centrality].shape[0]/permutation
        tmp_pval = tmp_pval/loop
        tmp_values = [row.id, row.centrality, tmp_pval]
        record_values.append(tmp_values)
    
    #new dataframe conains p-values
    new_df = pd.DataFrame(data=record_values)
    new_df.columns = ['id','centrality','pval_'+ str(permutation)+'sample'+str(loop)+'times']
    
    #for each chrom, assigne pval to its node
    for ci in cluster_label_g.keys():
        tmp_df = cluster_label_g[ci]
        tmp_df['id'] = ci + '_' + tmp_df.nodes.astype(str) 
        new_tmp_df = pd.merge(tmp_df, new_df, how='left', on='id').copy()
        
        out_path = out_data_folder + '/hic_data/' + bin_str + '/hic_community_data/' + cohort    
        out_file = ci + edges_rm_str + percentage_str + 'percent_noWeight_network_communities_info.tsv'
        out_file_name = os.path.join(out_path, out_file)
        print(out_file_name)
        new_tmp_df.to_csv(out_file_name, sep='\t', index=False)






