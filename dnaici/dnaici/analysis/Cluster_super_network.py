import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import json
import os
import glob
from dnaici.tools.data_functions import sort_dict


def group_nodes_to_clusters(cluster_df, in_df, minCluster_size):
    '''
    compute super network based on clustered labels
    make a new adjancy matrix for super clusters
    cluster_df: dataframe with two columns, nodes and labels
    clusters_labels: numpy array of all cluster labels 
    in_df: dataframe with two columns, a pair of edges
    return
    new_edges_in_cluster: a dictionary of cluster->nodes
    new_matrix2: adj-matrix for clustered network
    '''
    clusters_labels = cluster_df.labels.unique()
    num_of_clusters = clusters_labels.shape[0]
    #create initial adj-matrix for clustered network
    adj_matrix4clusters = np.zeros((num_of_clusters,num_of_clusters))
    #initinal clustered nodes and edges
    nodes_in_clusters_dict = {}
    num_of_edges_in_cluster_dict = {}
    edges_in_cluster_dict = {}
    #loop in clusters labels and find all nodes for it.
    for i in clusters_labels:
        tmp = cluster_df.nodes[cluster_df.labels==i].to_list()
        nodes_in_clusters_dict[i] = set(tmp)
        num_of_edges_in_cluster_dict[i] = 0
        edges_in_cluster_dict[i] = []
    #for each pair of edges assign them to the new super cluster
    cluster_df.labels = cluster_df.labels.astype(int)
    #loop in each pair of edges to find its cluster label
    for idx,row in in_df.iterrows():
        label2node_row0 = cluster_df.labels[cluster_df.nodes==row[0]].to_list()[0]
        label2node_row1 = cluster_df.labels[cluster_df.nodes==row[1]].to_list()[0]
        if label2node_row0 == label2node_row1:
            #both nodes located in the same cluster
            num_of_edges_in_cluster_dict[label2node_row0] += 1
            edges_in_cluster_dict[label2node_row0].append([row[0],row[1]])
            #donot count self-loop
            #adj_matrix4clusters[label2node_row0,label2node_row1] +=1
        else:
            #add a weight to edge if two nodes located in two different clusters
            adj_matrix4clusters[label2node_row0,label2node_row1] += 1
            adj_matrix4clusters[label2node_row1,label2node_row0] += 1 
    #remove columns with all 0
    #new_matrix= adj_matrix4clusters[:,~(adj_matrix4clusters==0).all(axis=0)]
    #remove rows with all 0
    #new_matrix2= new_matrix[~(new_matrix==0).all(axis=1),:]
    #new_matrix2=adj_matrix4clusters.copy()
    
    #build a dictionary for cluster -> nodes 
    new_num_of_edges_in_cluster = {}
    new_edges_in_cluster = {}
    cluster_label_without_zero_edges = []
    for ki in num_of_edges_in_cluster_dict.keys():
        if num_of_edges_in_cluster_dict[ki] >= minCluster_size:
            cluster_label_without_zero_edges.append(int(ki))
            new_num_of_edges_in_cluster[int(ki)] = [num_of_edges_in_cluster_dict[ki]]
            new_edges_in_cluster[int(ki)] = edges_in_cluster_dict[ki]
    #remove cluster labels with zero edges from ajd-matrix
    cluster_label_without_zero_edges.sort(),
    new_matrix = adj_matrix4clusters[cluster_label_without_zero_edges,:].copy()
    new_matrix2 = new_matrix[:,cluster_label_without_zero_edges].copy()   
    #print(new_num_of_edges_in_cluster)
    #print(np.triu(new_matrix2))
    return new_num_of_edges_in_cluster.copy(), new_edges_in_cluster.copy(), new_matrix2.copy()


def plot_a_network_by_matrix(new_num_of_edges_in_cluster0, new_matrix2):
    '''
    build a network based on adj-matrix
    new_matrix2: numpy 2-D adj-matrix
    '''
    A = np.triu(new_matrix2).copy()
    G = nx.from_numpy_array(A)
    #make color for clusters
    color_code1 = {-1: np.array([0,0,0]), #black 
                   0: np.array([255,255,255]), #white
                   1: np.array([0, 255, 0]), # green
                   3: np.array([255, 0, 0]), # red
                   2: np.array([255, 255,0]), #yellow
                   4: np.array([153, 0, 153]), #purpe
                   5: np.array([51, 51,255]), #dark blue
                   6: np.array([255, 204, 204]), #pink
                   7: np.array([0, 0, 255]), # blue 
                   8: np.array([0, 255,255]), #light blue
                   9: np.array([153, 51, 255]),
                   10: np.array([160,160,160]), #gray
                   11: np.array([153,255,204]), #light green
                   12: np.array([153,0,76]),
                   13: np.array([0 ,128, 255]),
                   14: np.array([229,255,204]),
                   15: np.array([153,153,0]),
                   16: np.array([51, 0, 0 ]), # black	
                   17: np.array([255, 0, 255]), # magenta
                   18: np.array([192, 192, 192]), # silver
                   19: np.array([240, 230, 140]), # khaki
                   20: np.array([128, 0, 0]), # maroon
                   21: np.array([128, 128, 0]), # olive
                   22: np.array([165, 42, 42]), # brown
                   23: np.array([240, 128, 128]), # light coral
                   24: np.array([143, 188, 143]), # dark sea green
                   25: np.array([221, 160, 221]), # plum
                   26: np.array([210, 105, 30]), # chocolate
                   27: np.array([188, 143, 143]), # rosy brown
                   28: np.array([255, 255, 240]), # ivory
                   29: np.array([95, 158, 160]), # cadet blue
                   30: np.array([216, 191, 216]), # thistle
                   31: np.array([255, 250, 250]) # snow	     
                   }   
    color_code2 = {-1: np.array([0,0,0]), #black 
                   0: np.array([255,255,255]), #white
                   1: np.array([0, 255, 0]), # green
                   2: np.array([255, 0,0]) #red
                   }
    if len(new_num_of_edges_in_cluster0.keys()) > 2:
        color_code0 = color_code1
    else:
        color_code0 = color_code2
    #fill the rest of cluster as black color
    for i in range(len(color_code0), 100):
        color_code0[i] = np.array([0, 0 ,0])
    #normalize color map for figure legend purpose
    color_code = {}
    for ki in color_code0.keys():
        color_code[ki] = list(color_code0[ki]/255)+[1]
    #fig, ax=plt.subplots(1, figsize=(6,5))
    fig, ax = plt.subplots()
    plt.clf()
    #sort dict by its key , the cluster number, then
    #it has the same row/column order as new_matrix2 
    new_num_of_edges_in_cluster = dict(sorted(new_num_of_edges_in_cluster0.items()))
    colors = []
    node_size = []
    edge_weight = []
    for ki in new_num_of_edges_in_cluster.keys():
        t_ki = int(ki)
        #sorted new_num_of_edges_in_cluster has the same key order as the new_matrix2
        #then we can simply append the corresponding colors..
        colors.append(color_code[t_ki+1])
        node_size.append(new_num_of_edges_in_cluster[t_ki][0])
    #normalize node size based on maximum node
    node_size = list(np.array(node_size)/max(node_size)*1000)
    edges = G.edges()
    #normalize edge width baed on maximum width
    edge_weight = [G[u][v]['weight'] for u, v in edges]
    if len(edge_weight) == 0:
        edge_weight = [0]
    edge_weight = list( np.array(edge_weight)/max(edge_weight)*10)
    #set node label by original node id +1 because the original node id start from 0
    labeldict = {}
    num_of_row = A.shape[0]
    #sort dictionary by key if the key > G matrix then set the new key from num_of_row +1 ....
    new_num_of_edges_in_cluster = sort_dict(new_num_of_edges_in_cluster)
    loop = num_of_row-1
    for ki in new_num_of_edges_in_cluster.keys():
        if ki <= num_of_row-1:
            labeldict[ki] = ki+1
        else:
            labeldict[loop] = ki+1
            loop += 1
    #graph options
    options = {'node_color': colors,
               'node_size': node_size,
               'width': edge_weight,
               'labels': labeldict,
               'with_labels': True,
               'font_weight': 'bold',
               'edge_color': 'blue'
               }
    #nx.draw_spring(G, **options)
    nx.draw_circular(G, **options)
    #plt.show()
    return fig, ax, G

def export_supernetwork_data(G, chrom_str, edges_rm_str, percent_str, in_path, new_num_of_edges_in_cluster, new_edges_in_cluster):
    '''
    Export supernetwork data to cytoscape or json format
    input: G is a graph from the networkx 
    '''
    t = nx.cytoscape_data(G) 
    #change int to str for nodes names
    for ti in t['elements']['edges']:
        ti['data']['source'] = str(ti['data']['source'])
        ti['data']['target'] = str( ti['data']['target'])
    #export network and nodes size to cytoscape json format
    out_network_file = chrom_str + edges_rm_str + percent_str + 'percent_noWeight_supernetwork.json'
    out_node_f = chrom_str + edges_rm_str + percent_str + 'percent_noWeight_supernetwork_numOfedges_in_node.tsv'
    out_edge_f = chrom_str + edges_rm_str + percent_str + 'percent_noWeight_supernetwork_edges_in_node.json'
    
    out_network_file = os.path.join(in_path,out_network_file)
    #print(out_network_file)
    with open(out_network_file, 'w', encoding='utf-8') as f:
        json.dump(t, f, ensure_ascii=False, indent=4)
    
    node_df = pd.DataFrame().from_dict(new_num_of_edges_in_cluster).T
    node_df.insert(0,'node',node_df.index.to_list())
    node_df['node_size'] = node_df[0].copy()
    node_df.pop(0)
    out_node_file = os.path.join(in_path,out_node_f)
    #print(out_node_file)
    node_df.to_csv(out_node_file,sep='\t',index=False)
    
    out_edge_file = os.path.join(in_path,out_edge_f)
    #print(out_edge_file)
    
    class NpEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, np.integer):
                return int(obj)
            if isinstance(obj, np.floating):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return json.JSONEncoder.default(self, obj)
    
    with open(out_edge_file, "w") as outfile:
        json.dump(new_edges_in_cluster, outfile,cls=NpEncoder)
    
    return out_network_file, out_node_file, out_edge_file


def find_sorted_network_cluster_id(old2sorted_id, new_num_of_edges_in_cluster, new_edges_in_cluster, new_matrix2):
    '''
    find sorted network cluster id based on new order that is sorted by geomic features of clusters 
    '''
    #here used sorted id for network clustering if it is needed ??, currently, use all of nodes
    sorted_num_of_edges_in_cluster = {}
    sorted_edges_in_cluster = {}
    
    for ki in old2sorted_id.items():
        sorted_num_of_edges_in_cluster[ki[1]] = new_num_of_edges_in_cluster[ki[0]]
        sorted_edges_in_cluster[ki[1]] = new_edges_in_cluster[ki[0]]
    #new_num_of_edges_in_cluster.pop(ki[0])
    #new_edges_in_cluster.pop(ki[0])
    #reorder the matrix based on only sorted clustered cn ID
    m_shape = new_matrix2.shape
    sorted_matrix2 = np.zeros(m_shape)
    for ii in range(0, m_shape[0]):
        for jj in range(0, m_shape[1]):
            if ii in old2sorted_id.keys():
                new_ii = old2sorted_id[ii]
            else:
                new_ii = ii
        #
            if jj in old2sorted_id.keys():
                new_jj = old2sorted_id[jj]
            else:
                new_jj = jj
            #only record both nodes in the sorted df by features
            # print(new_ii, new_jj)
            # if (new_ii != -1) & (new_jj != -1) :
            sorted_matrix2[new_ii][new_jj] = new_matrix2[ii][jj]    
    #remove rows having all zeroes
    #sorted_matrix2  = sorted_matrix2[~np.all(sorted_matrix2 == 0, axis=1),:]
    #remove columns having all zeroes
    #sorted_matrix2  = sorted_matrix2[:,~np.all(sorted_matrix2 == 0, axis=0)]
    return sorted_num_of_edges_in_cluster, sorted_edges_in_cluster, sorted_matrix2


def main(in_data_folder, 
         cohort, 
         chromosome, 
         resolution, 
         out_data_folder,
         minCluster_size,
         fig_dpi):
    #Parameters of input files
    bin_str = str(int(resolution/1000)) + 'kb'
    
    if chromosome == 'whole_genome':
        chrom_strs = ['chr'+str(i) for i in range(1,24)]
    else: 
        chrom_strs = chromosome

    isFind_global_cutoff = False
    percent_str = '0'
    
    if isFind_global_cutoff:
        edges_rm_str = '_edges_rmGlob_'
    else:
        edges_rm_str = '_edges_rm_'
    
    in_path = in_data_folder + '/hic_data/' + bin_str + '/hic_community_data/' + cohort
    in_path_nodes = in_data_folder + '/hic_data/' + bin_str + '/hic_community_figures/' + cohort
    
    for chrom_str in chrom_strs:
        print('- - - - - - - - investigating %s - - - - - - - -' %chrom_str, '\n')
        #load network cluststering results for the nodes
        in_file = os.path.join(in_path, chrom_str + edges_rm_str + percent_str + 'percent_noWeight_4communities_class.tsv')
        cluster_df = pd.read_csv(in_file, sep='\t')
        print('Read: \n', in_file, '\n')
        #load sorted cluster label file
        in_path4nodes = os.path.join(in_path_nodes, chrom_str)
        in_id_file = glob.glob(os.path.join(in_path4nodes, chrom_str+'_*'+edges_rm_str+percent_str+'percent_*'+'_ID2nodes.tsv'))[0]
        cluster_id_df = pd.read_csv(in_id_file, sep='\t')
        print('Read: \n', in_id_file, '\n')
        #build a dictionary for original cluster id to sorted cluster id
        old2sorted_id={}
        for idx,row in cluster_id_df.iterrows():
            old2sorted_id[row.original_network_cluster_node] = int(row.sorted_network_cluster_node_cn)     
        #load hi-c interaction matrix which may already filtered in various analysis 
        f = os.path.join(in_path, chrom_str + edges_rm_str + percent_str +'percent_noWeight_4zscore.matrix')
        in_df = pd.read_csv(f, sep='\t',header=None)
        print('Read: \n', f, '\n')
        
        new_num_of_edges_in_cluster, new_edges_in_cluster, new_matrix2 = group_nodes_to_clusters(cluster_df, in_df, minCluster_size)
        #here used sorted id (_ID2nodes.tsv ) for network clustering if it is needed ??, currently, use all of clustered nodes
        sorted_num_of_edges_in_cluster, sorted_edges_in_cluster, sorted_matrix2 = find_sorted_network_cluster_id(old2sorted_id, new_num_of_edges_in_cluster, new_edges_in_cluster, new_matrix2)
        #figure plot based on sorted cluster id cn_
        #print('Valid communities: \n', sorted_num_of_edges_in_cluster, '\n')
        fig, ax, G = plot_a_network_by_matrix(sorted_num_of_edges_in_cluster, sorted_matrix2)
                #whether shall limite the plot nodes only for available in old2sorted_id ???
        #export network based on sorted cluster id cn_ in  s4_read_clustered_network_homer.py
        f1,f2,f3 = export_supernetwork_data(G, chrom_str, edges_rm_str, percent_str, in_path, sorted_num_of_edges_in_cluster, sorted_edges_in_cluster)
        print('Export network clusterirng based on sorted ID in features: \n', '1: ',f1,'\n', '2: ',f2,'\n', '3: ',f3,'\n')
        
        out_path = out_data_folder + '/hic_data/' + bin_str + '/hic_community_figures/' + cohort + '/' + chrom_str 
        out_fig = out_path + '/' + f1.split('/')[-1]
        out_fig = out_fig.replace('.json','.jpg')
        print('Export figure of super network: \n', out_fig, '\n')
        fig.savefig(out_fig, dpi=fig_dpi)
        plt.close('all')
        

 
